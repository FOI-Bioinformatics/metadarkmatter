"""
Base classifier implementation using ANI-weighted placement.

This module provides the ANIWeightedClassifier, the core classification
algorithm that combines BLAST alignment scores with ANI (Average Nucleotide
Identity) distances to identify novel bacterial taxa.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Iterator

import numpy as np
import polars as pl

from metadarkmatter.core.constants import (
    calculate_alignment_quality,
    calculate_confidence_score,
    calculate_discovery_score,
    calculate_identity_confidence,
    calculate_inferred_uncertainty,
    calculate_placement_confidence,
)
from metadarkmatter.core.io_utils import write_dataframe, write_dataframe_append
from metadarkmatter.core.parsers import BlastResultFast, StreamingBlastParser
from metadarkmatter.models.blast import BlastResult
from metadarkmatter.models.classification import (
    ReadClassification,
    TaxonomicCall,
    TAXONOMIC_TO_DIVERSITY,
)
from metadarkmatter.models.config import ScoringConfig

if TYPE_CHECKING:
    from metadarkmatter.core.aai_matrix_builder import AAIMatrix
    from metadarkmatter.core.classification.ani_matrix import ANIMatrix

class ANIWeightedClassifier:
    """
    Classifier for detecting novel microbial diversity using ANI-weighted placement.

    Combines BLAST competitive recruitment results with precomputed ANI matrices
    to calculate novelty index and placement uncertainty metrics, enabling
    detection of novel species and genera in environmental samples.

    Algorithm:
        1. For each read, identify top BLAST hit (highest bitscore)
        2. Calculate Novelty Index (N) = 100 - top_hit_identity
        3. Find secondary hits within 95% of top bitscore
        4. Calculate Placement Uncertainty (U) = 100 - max(ANI(top, secondary))
        5. Classify based on N and U thresholds
    """

    def __init__(
        self,
        ani_matrix: ANIMatrix,
        aai_matrix: AAIMatrix | None = None,
        config: ScoringConfig | None = None,
    ) -> None:
        """
        Initialize ANI-weighted classifier.

        Args:
            ani_matrix: Precomputed ANI matrix for genome comparisons
            aai_matrix: Optional AAI matrix for genus-level classification.
                AAI provides more reliable genus-level boundaries than ANI
                at high divergence (ANI becomes unreliable below ~80%).
            config: Scoring configuration (uses defaults if None).
                When alignment_mode is "protein", classification uses wider
                novelty thresholds appropriate for protein-level identity.
        """
        self.ani_matrix = ani_matrix
        self.aai_matrix = aai_matrix
        self.config = config or ScoringConfig()

        # Store effective thresholds based on alignment mode
        # This allows protein mode to use wider novelty thresholds
        self._effective_thresholds = self.config.get_effective_thresholds()

    def classify_read(
        self,
        blast_result: BlastResult,
        read_length: int | None = None,
    ) -> ReadClassification | None:
        """
        Classify a single read based on its BLAST hits.

        Args:
            blast_result: BLAST results for one read
            read_length: Optional explicit read length for coverage calculation.
                         If None, uses qlen from hits or qend as proxy.

        Returns:
            ReadClassification if classification successful, None if no hits
        """
        if not blast_result.hits:
            return None

        # Determine read_length for coverage calculation if not provided
        if read_length is None and self.config.coverage_weight_mode != "none":
            # Try to get qlen from first hit, otherwise use qend as proxy
            first_hit = blast_result.hits[0]
            if first_hit.qlen is not None:
                read_length = first_hit.qlen
            else:
                # Use max qend as conservative proxy for read length
                read_length = max(hit.qend for hit in blast_result.hits)

        # Get top hit using coverage-weighted scoring
        if self.config.coverage_weight_mode != "none" and read_length is not None:
            best_hit = blast_result.get_best_hit_weighted(
                read_length=read_length,
                mode=self.config.coverage_weight_mode,
                strength=self.config.coverage_weight_strength,
            )
        else:
            best_hit = blast_result.best_hit

        if best_hit is None:
            return None

        # pident is already clamped to [0, 100] by BlastHit validator
        top_hit_identity = best_hit.pident

        # Novelty Index measures divergence from reference
        # N = 0 means identical match, N = 100 means no similarity
        novelty_index = 100.0 - top_hit_identity
        best_genome = best_hit.genome_name

        # Find second-best hit (to a DIFFERENT genome)
        second_hit_identity: float | None = None
        identity_gap: float | None = None
        for hit in blast_result.hits:
            if hit.genome_name != best_genome:
                second_hit_identity = hit.pident
                # Clamp to 0 - secondary hit can have higher identity than top bitscore hit
                identity_gap = max(0.0, top_hit_identity - second_hit_identity)
                break

        # Calculate placement uncertainty and count ambiguous hits in single pass
        # Uses iterator with early termination for efficiency
        if self.config.coverage_weight_mode != "none" and read_length is not None:
            ambiguous_hits_iter = blast_result.iter_ambiguous_hits_weighted(
                read_length=read_length,
                mode=self.config.coverage_weight_mode,
                strength=self.config.coverage_weight_strength,
                threshold_pct=self.config.bitscore_threshold_pct,
            )
        else:
            ambiguous_hits_iter = blast_result.iter_ambiguous_hits(
                self.config.bitscore_threshold_pct
            )

        placement_uncertainty, num_ambiguous_hits = self._calculate_placement_uncertainty(
            best_genome,
            ambiguous_hits_iter,
        )

        # Calculate genus-level uncertainty using lower threshold
        # This captures hits from other species in the same genus
        genus_uncertainty, num_genus_hits = self._calculate_genus_uncertainty(
            best_genome,
            best_hit.bitscore,
            iter(blast_result.hits),  # Need fresh iterator
        )

        # Calculate AAI-based uncertainty for genus-level classification
        # AAI is more reliable than ANI at high divergence (>20% novelty)
        aai_uncertainty: float | None = None
        if (
            self.aai_matrix is not None
            and self.config.use_aai_for_genus
            and novelty_index >= self.config.novelty_novel_genus_min
        ):
            # Collect competing genomes from ambiguous hits
            competing_genomes = [
                hit.genome_name
                for hit in blast_result.hits
                if hit.genome_name != best_genome
            ]
            aai_uncertainty = self._calculate_aai_uncertainty(best_genome, competing_genomes)

        # Determine taxonomic classification
        taxonomic_call = self._classify_by_thresholds(
            novelty_index,
            placement_uncertainty,
            num_ambiguous_hits=num_ambiguous_hits,
            genus_uncertainty=genus_uncertainty,
            num_genus_hits=num_genus_hits,
            identity_gap=identity_gap,
            aai_uncertainty=aai_uncertainty,
        )

        # Calculate confidence score using margin-based calculation
        # Use effective thresholds to account for protein vs nucleotide mode
        # Scaling parameters are also mode-specific for proper confidence calculation
        eff = self._effective_thresholds
        confidence_score = calculate_confidence_score(
            novelty_index=novelty_index,
            placement_uncertainty=placement_uncertainty,
            num_ambiguous_hits=num_ambiguous_hits,
            identity_gap=identity_gap,
            top_hit_identity=top_hit_identity,
            taxonomic_call=taxonomic_call.value,
            novelty_known_max=eff["novelty_known_max"],
            novelty_novel_species_max=eff["novelty_novel_species_max"],
            novelty_novel_genus_max=eff["novelty_novel_genus_max"],
            uncertainty_confident_max=eff["uncertainty_known_max"],
            # Mode-specific scaling parameters
            margin_divisor_known=eff["margin_divisor_known"],
            margin_divisor_novel_species=eff["margin_divisor_novel_species"],
            margin_divisor_novel_genus=eff["margin_divisor_novel_genus"],
            identity_gap_thresholds=eff["identity_gap_thresholds"],
            identity_score_base=eff["identity_score_base"],
            identity_score_range=eff["identity_score_range"],
        )

        # Initialize enhanced scoring fields
        inferred_uncertainty: float | None = None
        uncertainty_type: str | None = None
        alignment_quality: float | None = None
        identity_confidence: float | None = None
        placement_confidence: float | None = None
        discovery_score: float | None = None

        # Calculate inferred uncertainty for single-hit reads if enabled
        if self.config.infer_single_hit_uncertainty or self.config.enhanced_scoring:
            if num_ambiguous_hits <= 1:
                inferred_uncertainty = calculate_inferred_uncertainty(novelty_index)
                uncertainty_type = "inferred"
            else:
                uncertainty_type = "measured"

        # Calculate enhanced scoring metrics if enabled
        if self.config.enhanced_scoring:
            # Calculate coverage for alignment quality
            if read_length is not None and read_length > 0:
                coverage = (best_hit.qend - best_hit.qstart + 1) / read_length
            else:
                coverage = min(1.0, (best_hit.qend - best_hit.qstart + 1) / best_hit.qend)
            coverage = max(0.0, min(1.0, coverage))

            # Alignment quality from BLAST statistics
            alignment_quality = calculate_alignment_quality(
                mismatch=best_hit.mismatch,
                gapopen=best_hit.gapopen,
                length=best_hit.length,
                coverage=coverage,
                evalue=best_hit.evalue,
            )

            # Identity confidence
            identity_confidence = calculate_identity_confidence(
                novelty_index=novelty_index,
                alignment_quality=alignment_quality,
            )

            # Placement confidence - use inferred uncertainty for single-hit reads
            effective_uncertainty = inferred_uncertainty if num_ambiguous_hits <= 1 else placement_uncertainty
            effective_type = uncertainty_type or ("inferred" if num_ambiguous_hits <= 1 else "measured")

            placement_confidence = calculate_placement_confidence(
                uncertainty=effective_uncertainty,
                uncertainty_type=effective_type,
                identity_gap=identity_gap,
                num_ambiguous_hits=num_ambiguous_hits,
            )

            # Discovery score for novel reads
            discovery_score = calculate_discovery_score(
                taxonomic_call=taxonomic_call.value,
                novelty_index=novelty_index,
                identity_confidence=identity_confidence,
                placement_confidence=placement_confidence,
                alignment_quality=alignment_quality,
            )

        return ReadClassification(
            read_id=blast_result.read_id,
            best_match_genome=best_genome,
            top_hit_identity=top_hit_identity,
            novelty_index=novelty_index,
            placement_uncertainty=placement_uncertainty,
            num_ambiguous_hits=num_ambiguous_hits,
            second_hit_identity=second_hit_identity,
            identity_gap=identity_gap,
            genus_uncertainty=genus_uncertainty,
            num_genus_hits=num_genus_hits,
            confidence_score=confidence_score,
            inferred_uncertainty=inferred_uncertainty,
            uncertainty_type=uncertainty_type,
            alignment_quality=alignment_quality,
            identity_confidence=identity_confidence,
            placement_confidence=placement_confidence,
            discovery_score=discovery_score,
            taxonomic_call=taxonomic_call,
        )

    def _calculate_placement_uncertainty(
        self,
        best_genome: str,
        ambiguous_hits_iter: Iterator,
    ) -> tuple[float, int]:
        """
        Calculate placement uncertainty from ANI values to secondary hits.

        Two modes are supported (controlled by config.uncertainty_mode):
        - 'max': U = 100 - max(ANI(top_hit, all_secondary_hits))
        - 'second': U = 100 - ANI(top_hit, second_best_hit)

        Performance optimization: Accepts an iterator with early termination
        and counts hits in a single pass, avoiding list materialization.

        Args:
            best_genome: Genome identifier for top BLAST hit
            ambiguous_hits_iter: Iterator of BlastHit objects within bitscore threshold

        Returns:
            Tuple of (placement_uncertainty, num_ambiguous_hits)
        """
        # Process iterator in single pass - collect ANI values and count
        max_ani = 0.0
        second_genome_ani = 0.0  # ANI to second-best genome (by bitscore)
        second_genome_found = False  # Flag for first secondary genome
        num_ambiguous_hits = 0

        for hit in ambiguous_hits_iter:
            num_ambiguous_hits += 1
            secondary_genome = hit.genome_name
            if secondary_genome != best_genome:
                ani = self.ani_matrix.get_ani(best_genome, secondary_genome)

                # Track first secondary genome (second-best by bitscore)
                if not second_genome_found:
                    second_genome_ani = ani
                    second_genome_found = True

                # Track max ANI (existing behavior)
                if ani > max_ani:
                    max_ani = ani

        if num_ambiguous_hits <= 1:
            # Only one genome hit, no uncertainty
            return 0.0, num_ambiguous_hits

        # Choose ANI based on uncertainty mode
        if self.config.uncertainty_mode == "second":
            effective_ani = second_genome_ani
        else:  # "max" mode (default)
            effective_ani = max_ani

        if effective_ani == 0.0:
            # Genome not in ANI matrix at all - maximum uncertainty
            # Note: Missing ANI values within the matrix return default_ani (70%)
            return 100.0, num_ambiguous_hits

        # Uncertainty is inverse of effective ANI
        return 100.0 - effective_ani, num_ambiguous_hits

    def _calculate_genus_uncertainty(
        self,
        best_genome: str,
        best_bitscore: float,
        hits_iter: Iterator,
    ) -> tuple[float | None, int]:
        """
        Calculate genus-level uncertainty using lower bitscore threshold.

        Uses genus_bitscore_threshold_pct (default 90%) to capture hits from
        other species in the same genus that would be missed by the stricter
        95% threshold. This helps identify reads hitting conserved genes
        shared across species within the genus.

        Args:
            best_genome: Genome identifier for top BLAST hit
            best_bitscore: Bitscore of top hit
            hits_iter: Iterator of BlastHit objects (all hits for this read)

        Returns:
            Tuple of (genus_uncertainty, num_genus_hits) where:
            - genus_uncertainty: 100 - max(ANI) for genus-level hits, or None
            - num_genus_hits: Number of hits within genus bitscore threshold
        """
        genus_bitscore_cutoff = best_bitscore * (
            self.config.genus_bitscore_threshold_pct / 100.0
        )

        max_ani = 0.0
        num_genus_hits = 0

        for hit in hits_iter:
            if hit.bitscore < genus_bitscore_cutoff:
                break  # Hits are sorted by bitscore
            num_genus_hits += 1
            secondary_genome = hit.genome_name
            if secondary_genome != best_genome:
                ani = self.ani_matrix.get_ani(best_genome, secondary_genome)
                if ani > max_ani:
                    max_ani = ani

        if num_genus_hits <= 1:
            # Only one genome hit at genus level
            return None, num_genus_hits

        if max_ani == 0.0:
            # Genome not in ANI matrix - maximum uncertainty
            return 100.0, num_genus_hits

        return 100.0 - max_ani, num_genus_hits

    def _calculate_aai_uncertainty(
        self,
        best_genome: str,
        competing_genomes: list[str],
    ) -> float | None:
        """
        Calculate genus-level uncertainty using AAI matrix.

        AAI provides more reliable genus-level classification than ANI
        at high divergence. ANI becomes unreliable below approximately 80%
        identity, while AAI maintains accuracy for genus-level boundaries.

        AAI genus boundaries (Riesco & Trujillo 2024):
        - AAI > 65%: Same genus (aai_genus_boundary_high)
        - AAI 58-65%: Genus boundary zone
        - AAI < 58%: Different genus (aai_genus_boundary_low)

        Args:
            best_genome: Genome identifier for top BLAST hit
            competing_genomes: List of genome identifiers from secondary hits

        Returns:
            AAI-based uncertainty (100 - max AAI), or None if AAI matrix
            is not available or genome not in matrix.
        """
        if self.aai_matrix is None:
            return None

        if not competing_genomes:
            return None

        max_aai = 0.0
        found_any = False

        for secondary_genome in competing_genomes:
            if secondary_genome != best_genome:
                aai = self.aai_matrix.get_aai(best_genome, secondary_genome)
                if aai > 0:
                    found_any = True
                    if aai > max_aai:
                        max_aai = aai

        if not found_any:
            # No AAI values available for these genomes
            return None

        # Uncertainty is inverse of maximum AAI
        return 100.0 - max_aai

    def _classify_by_thresholds(
        self,
        novelty_index: float,
        placement_uncertainty: float,
        num_ambiguous_hits: int = 1,
        genus_uncertainty: float | None = None,
        num_genus_hits: int | None = None,
        identity_gap: float | None = None,
        aai_uncertainty: float | None = None,
    ) -> TaxonomicCall:
        """
        Classify read based on novelty and uncertainty thresholds.

        Classification decision tree (evaluated in order):

        +------+-----------------------------------------------+------------------------+
        | Rule | Condition                                     | Classification         |
        +------+-----------------------------------------------+------------------------+
        |  0   | identity_gap < 2% AND multiple hits           | Ambiguous              |
        |  1   | U >= 5%                                       | Ambiguous              |
        |  2   | 2% <= U < 5%                                  | Species Boundary       |
        |  3   | U < 2% AND N < 5%                             | Known Species          |
        |  4   | U < 2% AND 5% <= N < 20%                      | Novel Species          |
        |  5   | U < 2% AND 20% <= N <= 25% AND AAI > 65%      | Novel Species (AAI)    |
        |  5a  | U < 2% AND 20% <= N <= 25% AND AAI 58-65%     | Ambiguous Within Genus |
        |  5b  | U < 2% AND 20% <= N <= 25% AND AAI < 58%      | Novel Genus            |
        |  5c  | U < 2% AND 20% <= N <= 25% AND genus_U >= 10% | Ambiguous Within Genus |
        |  5d  | U < 2% AND 20% <= N <= 25% (no AAI)           | Novel Genus            |
        |  6   | Otherwise                                     | Unclassified           |
        +------+-----------------------------------------------+------------------------+

        Threshold basis:
        - ANI thresholds (Jain et al. 2018): 95-96% ANI = species boundary
        - AAI thresholds (Riesco & Trujillo 2024): 58-65% AAI = genus boundary
        - N (novelty) = 100 - pident, so N < 5% means pident > 95%
        - U (uncertainty) = 100 - max(ANI/AAI between competing genomes)

        AAI provides more reliable genus-level classification than ANI because
        ANI becomes unreliable below approximately 80% identity.

        To customize thresholds, modify ScoringConfig values. Rules are
        applied in order; first matching rule determines classification.

        Args:
            novelty_index: Novelty index value (0-100)
            placement_uncertainty: Placement uncertainty value (0-100)
            num_ambiguous_hits: Number of hits within bitscore threshold
            genus_uncertainty: ANI-based uncertainty for genus-level hits (optional)
            num_genus_hits: Number of hits within genus bitscore threshold (optional)
            identity_gap: Gap between best and second-best hit identity (optional)
            aai_uncertainty: AAI-based uncertainty for genus-level classification (optional).
                When available and use_aai_for_genus is enabled, this is used instead
                of ANI-based genus_uncertainty for genus-level decisions.

        Returns:
            TaxonomicCall classification
        """
        cfg = self.config
        # Use effective thresholds which account for alignment mode (protein vs nucleotide)
        eff = self._effective_thresholds

        # Rule 0: Identity gap check (before ANI-based uncertainty)
        # When two hits have nearly identical identity (gap < 2%), the read
        # cannot be confidently placed regardless of ANI between those genomes.
        if (
            identity_gap is not None
            and num_ambiguous_hits > 1
            and identity_gap < cfg.identity_gap_ambiguous_max
        ):
            return TaxonomicCall.AMBIGUOUS

        # Rule 0b: Single-hit inferred uncertainty check
        # For single-hit reads, no ANI-based uncertainty is measurable.
        # Use inferred uncertainty based on novelty level to flag ambiguous cases.
        # This addresses the ~70% of environmental reads that have only one hit.
        if (
            cfg.use_inferred_for_single_hits
            and num_ambiguous_hits <= 1
            and novelty_index >= eff["novelty_novel_species_min"]  # Only for novel range
        ):
            inferred = calculate_inferred_uncertainty(novelty_index)
            if inferred >= cfg.single_hit_uncertainty_threshold:
                return TaxonomicCall.AMBIGUOUS

        # Uncertainty-based checks come first (ANI-based, more biologically meaningful)
        # High uncertainty indicates conserved region shared within genus
        # Note: Cross-genera conserved regions are handled in Polars classifier
        # with ambiguity_scope check
        if placement_uncertainty >= eff["uncertainty_conserved_min"]:
            return TaxonomicCall.AMBIGUOUS

        # Moderate uncertainty (2-5%): species boundary zone
        # Read matches multiple closely related species (95-98% ANI)
        if placement_uncertainty >= eff["uncertainty_novel_genus_max"]:
            return TaxonomicCall.SPECIES_BOUNDARY

        # Known species: low novelty, low uncertainty
        # N < 5% (nucleotide) or N < 10% (protein) means pident > species boundary
        if (
            novelty_index < eff["novelty_known_max"]
            and placement_uncertainty < eff["uncertainty_known_max"]
        ):
            return TaxonomicCall.KNOWN_SPECIES

        # Novel species: moderate novelty, low uncertainty
        # Nucleotide: 5% <= N < 20%; Protein: 10% <= N < 25%
        if (
            eff["novelty_novel_species_min"] <= novelty_index < eff["novelty_novel_species_max"]
            and placement_uncertainty < eff["uncertainty_novel_species_max"]
        ):
            return TaxonomicCall.NOVEL_SPECIES

        # Novel genus candidates: check genus-level uncertainty
        # Nucleotide: 20% <= N <= 25%; Protein: 25% <= N <= 40%
        if (
            eff["novelty_novel_genus_min"] <= novelty_index <= eff["novelty_novel_genus_max"]
            and placement_uncertainty < eff["uncertainty_novel_genus_max"]
        ):
            # Use AAI-based classification when available and enabled
            # AAI is more reliable than ANI at genus-level divergence
            if aai_uncertainty is not None and cfg.use_aai_for_genus:
                # AAI thresholds (Riesco & Trujillo 2024):
                # - aai_uncertainty < 35% (AAI > 65%): Same genus, likely novel species
                # - aai_uncertainty 35-42% (AAI 58-65%): Genus boundary zone, ambiguous
                # - aai_uncertainty > 42% (AAI < 58%): Different genus, novel genus
                aai_genus_uncertainty_low = 100.0 - cfg.aai_genus_boundary_high  # 35%
                aai_genus_uncertainty_high = 100.0 - cfg.aai_genus_boundary_low   # 42%

                if aai_uncertainty < aai_genus_uncertainty_low:
                    # High AAI (>65%) indicates same genus despite high novelty
                    # This read is likely hitting a conserved region within the genus
                    # Classify as novel species rather than novel genus
                    return TaxonomicCall.NOVEL_SPECIES
                elif aai_uncertainty <= aai_genus_uncertainty_high:
                    # AAI in boundary zone (58-65%), ambiguous genus classification
                    return TaxonomicCall.AMBIGUOUS_WITHIN_GENUS
                else:
                    # Low AAI (<58%) confirms different genus
                    return TaxonomicCall.NOVEL_GENUS

            # Fallback to ANI-based genus uncertainty when AAI not available
            # Check for genus-level ambiguity: multiple species hits with low ANI
            # This catches reads that hit conserved genes shared across species
            # within the same genus (e.g., hitting F. philomiragia and F. endociliophora
            # with 88% ANI between them)
            if (
                genus_uncertainty is not None
                and num_genus_hits is not None
                and num_genus_hits > 1
                and genus_uncertainty >= cfg.genus_uncertainty_ambiguous_min
            ):
                return TaxonomicCall.AMBIGUOUS_WITHIN_GENUS

            return TaxonomicCall.NOVEL_GENUS

        # Default: doesn't fit clear biological categories
        # Includes: strain variants (N=2-5), or very high divergence (N>25)
        return TaxonomicCall.UNCLASSIFIED

    def classify_blast_file(
        self,
        blast_path: Path,
    ) -> Iterator[ReadClassification]:
        """
        Classify all reads from a BLAST file using streaming.

        Memory-efficient processing for large BLAST result files.

        Args:
            blast_path: Path to BLAST tabular output

        Yields:
            ReadClassification objects for each successfully classified read
        """
        parser = StreamingBlastParser(blast_path)

        for blast_result in parser.iter_reads():
            classification = self.classify_read(blast_result)
            if classification is not None:
                yield classification

    def classify_to_dataframe(
        self,
        blast_path: Path,
    ) -> pl.DataFrame:
        """
        Classify reads and return results as Polars DataFrame.

        Args:
            blast_path: Path to BLAST tabular output

        Returns:
            DataFrame with classification results
        """
        classifications = list(self.classify_blast_file(blast_path))

        if not classifications:
            # Return empty DataFrame with correct schema
            return pl.DataFrame(
                schema={
                    "read_id": pl.Utf8,
                    "best_match_genome": pl.Utf8,
                    "top_hit_identity": pl.Float64,
                    "novelty_index": pl.Float64,
                    "placement_uncertainty": pl.Float64,
                    "num_ambiguous_hits": pl.Int64,
                    "second_hit_identity": pl.Float64,
                    "identity_gap": pl.Float64,
                    "confidence_score": pl.Float64,
                    "taxonomic_call": pl.Utf8,
                    "is_novel": pl.Boolean,
                }
            )

        # Convert to list of dicts
        data = [c.to_dict() for c in classifications]

        return pl.DataFrame(data)

    def write_classifications(
        self,
        blast_path: Path,
        output_path: Path,
        output_format: str = "csv",
    ) -> int:
        """
        Classify reads and write results to file.

        Args:
            blast_path: Path to BLAST tabular output
            output_path: Path for output file
            output_format: Output format - "csv" or "parquet" (default: "csv")
                          Parquet is 10x smaller and faster for large datasets.

        Returns:
            Number of reads classified
        """
        df = self.classify_to_dataframe(blast_path)
        write_dataframe(df, output_path, output_format)
        return len(df)

    # =========================================================================
    # High-Performance Classification Methods (Phase 2 Optimizations)
    # =========================================================================
    # These methods use lightweight NamedTuples and NumPy ANI lookups
    # for ~10x faster processing on large BLAST files.

    def classify_read_fast(
        self,
        result: BlastResultFast,
    ) -> dict | None:
        """
        Classify a read using lightweight data structures.

        Performance optimizations vs classify_read():
        - Uses BlastResultFast (NamedTuple) instead of Pydantic model
        - Pre-extracted genome_name avoids regex in hot path
        - Returns dict directly, bypassing Pydantic on output
        - ~10x faster per-read classification

        Args:
            result: BlastResultFast with pre-extracted genome names

        Returns:
            Dict with classification results, or None if no hits
        """
        if not result.hits:
            return None

        # Get top hit (hits are pre-sorted by bitscore descending)
        best_hit = result.hits[0]
        top_hit_identity = min(100.0, max(0.0, best_hit.pident))  # Clamp

        # Novelty Index: divergence from reference
        novelty_index = 100.0 - top_hit_identity
        best_genome = best_hit.genome_name

        # Find second-best hit (to a DIFFERENT genome) for identity gap
        second_hit_identity: float | None = None
        identity_gap: float | None = None
        for hit in result.hits:
            if hit.genome_name != best_genome:
                second_hit_identity = hit.pident
                # Clamp to 0 - secondary hit can have higher identity than top bitscore hit
                identity_gap = max(0.0, top_hit_identity - second_hit_identity)
                break

        # Calculate placement uncertainty using pre-extracted genome names
        bitscore_cutoff = best_hit.bitscore * (self.config.bitscore_threshold_pct / 100.0)

        max_ani = 0.0
        second_genome_ani = 0.0  # ANI to second-best genome (by bitscore)
        second_genome_found = False  # Flag for first secondary genome
        num_ambiguous_hits = 0

        for hit in result.hits:
            if hit.bitscore < bitscore_cutoff:
                break  # Early termination - hits are sorted

            num_ambiguous_hits += 1
            if hit.genome_name != best_genome:
                ani = self.ani_matrix.get_ani(best_genome, hit.genome_name)

                # Track first secondary genome (second-best by bitscore)
                if not second_genome_found:
                    second_genome_ani = ani
                    second_genome_found = True

                # Track max ANI
                if ani > max_ani:
                    max_ani = ani

        # Calculate placement uncertainty based on mode
        if num_ambiguous_hits <= 1:
            placement_uncertainty = 0.0
        else:
            # Choose ANI based on uncertainty mode
            if self.config.uncertainty_mode == "second":
                effective_ani = second_genome_ani
            else:  # "max" mode (default)
                effective_ani = max_ani

            if effective_ani == 0.0:
                # Genome not in ANI matrix - maximum uncertainty
                placement_uncertainty = 100.0
            else:
                placement_uncertainty = 100.0 - effective_ani

        # Calculate AAI-based uncertainty for genus-level classification
        # Only computed when AAI matrix is available and read is in genus-level novelty range
        aai_uncertainty: float | None = None
        if (
            self.aai_matrix is not None
            and self.config.use_aai_for_genus
            and novelty_index >= self.config.novelty_novel_genus_min
        ):
            competing_genomes = [
                hit.genome_name
                for hit in result.hits
                if hit.genome_name != best_genome
            ]
            aai_uncertainty = self._calculate_aai_uncertainty(best_genome, competing_genomes)

        # Determine classification
        taxonomic_call = self._classify_by_thresholds(
            novelty_index,
            placement_uncertainty,
            num_ambiguous_hits=num_ambiguous_hits,
            identity_gap=identity_gap,
            aai_uncertainty=aai_uncertainty,
        )

        # Calculate confidence score using effective thresholds
        # Scaling parameters are mode-specific for proper confidence calculation
        eff = self._effective_thresholds
        confidence_score = calculate_confidence_score(
            novelty_index=novelty_index,
            placement_uncertainty=placement_uncertainty,
            num_ambiguous_hits=num_ambiguous_hits,
            identity_gap=identity_gap,
            top_hit_identity=top_hit_identity,
            taxonomic_call=taxonomic_call.value,
            novelty_known_max=eff["novelty_known_max"],
            novelty_novel_species_max=eff["novelty_novel_species_max"],
            novelty_novel_genus_max=eff["novelty_novel_genus_max"],
            uncertainty_confident_max=eff["uncertainty_known_max"],
            # Mode-specific scaling parameters
            margin_divisor_known=eff["margin_divisor_known"],
            margin_divisor_novel_species=eff["margin_divisor_novel_species"],
            margin_divisor_novel_genus=eff["margin_divisor_novel_genus"],
            identity_gap_thresholds=eff["identity_gap_thresholds"],
            identity_score_base=eff["identity_score_base"],
            identity_score_range=eff["identity_score_range"],
        )

        # Build result dict
        result_dict = {
            "read_id": result.read_id,
            "best_match_genome": best_genome,
            "top_hit_identity": top_hit_identity,
            "novelty_index": novelty_index,
            "placement_uncertainty": placement_uncertainty,
            "num_ambiguous_hits": num_ambiguous_hits,
            "second_hit_identity": second_hit_identity,
            "identity_gap": identity_gap,
            "confidence_score": confidence_score,
            "taxonomic_call": taxonomic_call.value,
            "is_novel": taxonomic_call in (TaxonomicCall.NOVEL_SPECIES, TaxonomicCall.NOVEL_GENUS),
        }

        # Add enhanced scoring fields if enabled
        if self.config.infer_single_hit_uncertainty or self.config.enhanced_scoring:
            if num_ambiguous_hits <= 1:
                result_dict["inferred_uncertainty"] = calculate_inferred_uncertainty(novelty_index)
                result_dict["uncertainty_type"] = "inferred"
            else:
                result_dict["uncertainty_type"] = "measured"

        if self.config.enhanced_scoring:
            # Calculate coverage for alignment quality
            coverage = (best_hit.qend - best_hit.qstart + 1) / max(1, best_hit.qend)
            coverage = max(0.0, min(1.0, coverage))

            alignment_quality = calculate_alignment_quality(
                mismatch=best_hit.mismatch,
                gapopen=best_hit.gapopen,
                length=best_hit.length,
                coverage=coverage,
                evalue=best_hit.evalue,
            )
            result_dict["alignment_quality"] = alignment_quality

            identity_confidence = calculate_identity_confidence(
                novelty_index=novelty_index,
                alignment_quality=alignment_quality,
            )
            result_dict["identity_confidence"] = identity_confidence

            # Use inferred uncertainty for single-hit reads
            effective_uncertainty = result_dict.get("inferred_uncertainty", placement_uncertainty)
            effective_type = result_dict.get("uncertainty_type", "measured")

            placement_confidence = calculate_placement_confidence(
                uncertainty=effective_uncertainty,
                uncertainty_type=effective_type,
                identity_gap=identity_gap,
                num_ambiguous_hits=num_ambiguous_hits,
            )
            result_dict["placement_confidence"] = placement_confidence

            discovery_score = calculate_discovery_score(
                taxonomic_call=taxonomic_call.value,
                novelty_index=novelty_index,
                identity_confidence=identity_confidence,
                placement_confidence=placement_confidence,
                alignment_quality=alignment_quality,
            )
            if discovery_score is not None:
                result_dict["discovery_score"] = discovery_score

        return result_dict

    def classify_blast_file_fast(
        self,
        blast_path: Path,
    ) -> Iterator[dict]:
        """
        Classify reads using high-performance streaming.

        Performance optimizations vs classify_blast_file():
        - Uses iter_reads_fast() with NamedTuples (~50x faster object creation)
        - Vectorized genome extraction in Polars (~100x faster)
        - Returns dicts instead of Pydantic models

        Args:
            blast_path: Path to BLAST tabular output

        Yields:
            Dict with classification results for each read
        """
        parser = StreamingBlastParser(blast_path)

        for result in parser.iter_reads_fast():
            classification = self.classify_read_fast(result)
            if classification is not None:
                yield classification

    def classify_to_dataframe_fast(
        self,
        blast_path: Path,
    ) -> pl.DataFrame:
        """
        Classify reads and return Polars DataFrame directly.

        This is the fastest method for processing large BLAST files.
        Bypasses Pydantic on both input and output for maximum performance.

        Performance vs classify_to_dataframe():
        - ~10x faster for large files (100M+ reads)
        - ~5x less memory usage
        - Streaming-compatible

        Args:
            blast_path: Path to BLAST tabular output

        Returns:
            DataFrame with classification results
        """
        # Collect all classifications for DataFrame construction
        all_rows: list[dict] = list(self.classify_blast_file_fast(blast_path))

        if not all_rows:
            # Return empty DataFrame with correct schema
            return pl.DataFrame(
                schema={
                    "read_id": pl.Utf8,
                    "best_match_genome": pl.Utf8,
                    "top_hit_identity": pl.Float64,
                    "novelty_index": pl.Float64,
                    "placement_uncertainty": pl.Float64,
                    "num_ambiguous_hits": pl.Int64,
                    "second_hit_identity": pl.Float64,
                    "identity_gap": pl.Float64,
                    "confidence_score": pl.Float64,
                    "taxonomic_call": pl.Utf8,
                    "is_novel": pl.Boolean,
                }
            )

        return pl.DataFrame(all_rows)

    # =========================================================================
    # Phase 3: Scalability Optimizations
    # =========================================================================
    # These methods enable processing of very large files (100M+ reads)
    # with bounded memory and parallel execution.

    # Bounds for chunk_size parameters
    MIN_CHUNK_SIZE = 100
    MAX_CHUNK_SIZE = 10_000_000

    def stream_to_file_fast(
        self,
        blast_path: Path,
        output_path: Path,
        output_format: str = "parquet",
        chunk_size: int = 100_000,
        progress_callback: Callable[[int, int], None] | None = None,
    ) -> int:
        """
        Stream classification results directly to disk in chunks.

        Memory-efficient processing for very large files. Instead of
        accumulating all results in memory, writes batches to disk.

        For Parquet output, uses row group batching for optimal compression.
        For CSV output, appends chunks incrementally.

        Args:
            blast_path: Path to BLAST tabular output.
            output_path: Output file path.
            output_format: 'parquet' or 'csv' (default: parquet).
            chunk_size: Number of classifications per batch (default: 100K).
                        Valid range: 100 to 10,000,000.
            progress_callback: Optional callback(processed, total) for progress.

        Returns:
            Total number of reads classified.

        Raises:
            ValueError: If chunk_size is outside valid range.

        Memory usage:
            Each classification dict uses ~500 bytes. chunk_size=100K uses ~50MB.

        Example:
            def show_progress(done, total):
                print(f"Processed {done:,} reads...")

            classifier.stream_to_file_fast(
                blast_path,
                output_path,
                progress_callback=show_progress
            )
        """
        if not self.MIN_CHUNK_SIZE <= chunk_size <= self.MAX_CHUNK_SIZE:
            msg = (
                f"chunk_size must be between {self.MIN_CHUNK_SIZE:,} and "
                f"{self.MAX_CHUNK_SIZE:,}, got {chunk_size:,}"
            )
            raise ValueError(msg)
        total_classified = 0
        chunk_buffer: list[dict] = []
        first_chunk = True

        for classification in self.classify_blast_file_fast(blast_path):
            chunk_buffer.append(classification)
            total_classified += 1

            if len(chunk_buffer) >= chunk_size:
                self._write_chunk(
                    chunk_buffer,
                    output_path,
                    output_format,
                    first_chunk,
                )
                first_chunk = False
                chunk_buffer = []

                if progress_callback:
                    progress_callback(total_classified, 0)

        # Write remaining data
        if chunk_buffer:
            self._write_chunk(
                chunk_buffer,
                output_path,
                output_format,
                first_chunk,
            )

        if progress_callback:
            progress_callback(total_classified, total_classified)

        return total_classified

    def _write_chunk(
        self,
        chunk: list[dict],
        output_path: Path,
        output_format: str,
        is_first: bool,
    ) -> None:
        """Write a chunk of classifications to file."""
        df = pl.DataFrame(chunk)
        write_dataframe_append(df, output_path, output_format, is_first)


# =============================================================================
# Parallel Classification Worker Functions (Module Level for Pickling)
# =============================================================================


def _classify_chunk_worker(
    chunk_data: list[ReadChunk],
    ani_array: np.ndarray,
    genome_to_idx: dict[str, int],
    config_dict: dict[str, float],
) -> list[ClassificationDict]:
    """
    Worker function for parallel classification.

    Runs in a separate process. Receives raw data and config,
    performs classification, returns results.

    Args:
        chunk_data: List of (read_id, hits_data) tuples where hits_data is
            a list of (qseqid, sseqid, pident, bitscore, genome_name) tuples.
        ani_array: NumPy ANI matrix (shared memory).
        genome_to_idx: Genome name to index mapping.
        config_dict: Scoring configuration as dict.

    Returns:
        List of classification result dicts with keys:
            read_id, best_match_genome, top_hit_identity, top_bitscore,
            competing_genomes, novelty_index, placement_uncertainty,
            taxonomic_call.
    """
    results = []

    # Reconstruct config thresholds (classification boundaries)
    bitscore_pct = config_dict["bitscore_threshold_pct"]
    genus_bitscore_pct = config_dict["genus_bitscore_threshold_pct"]
    novelty_known_max = config_dict["novelty_known_max"]
    novelty_novel_species_min = config_dict["novelty_novel_species_min"]
    novelty_novel_species_max = config_dict["novelty_novel_species_max"]
    novelty_novel_genus_min = config_dict["novelty_novel_genus_min"]
    novelty_novel_genus_max = config_dict["novelty_novel_genus_max"]
    uncertainty_known_max = config_dict["uncertainty_known_max"]
    uncertainty_novel_species_max = config_dict["uncertainty_novel_species_max"]
    uncertainty_novel_genus_max = config_dict["uncertainty_novel_genus_max"]
    uncertainty_conserved_min = config_dict["uncertainty_conserved_min"]
    genus_uncertainty_ambiguous_min = config_dict["genus_uncertainty_ambiguous_min"]
    default_ani = config_dict["default_ani"]
    identity_gap_ambiguous_max = config_dict["identity_gap_ambiguous_max"]

    # Confidence score scaling parameters (mode-specific)
    margin_divisor_known = config_dict["margin_divisor_known"]
    margin_divisor_novel_species = config_dict["margin_divisor_novel_species"]
    margin_divisor_novel_genus = config_dict["margin_divisor_novel_genus"]
    identity_gap_thresholds = config_dict["identity_gap_thresholds"]
    identity_score_base = config_dict["identity_score_base"]
    identity_score_range = config_dict["identity_score_range"]

    # Enhanced scoring options
    enhanced_scoring = config_dict.get("enhanced_scoring", False)
    infer_single_hit_uncertainty = config_dict.get("infer_single_hit_uncertainty", False)

    # Single-hit classification options
    use_inferred_for_single_hits = config_dict.get("use_inferred_for_single_hits", False)
    single_hit_uncertainty_threshold = config_dict.get("single_hit_uncertainty_threshold", 10.0)

    for read_id, hits_data in chunk_data:
        if not hits_data:
            continue

        # Reconstruct hits (already sorted by bitscore)
        best_hit = hits_data[0]
        best_pident = best_hit[2]  # pident
        best_bitscore = best_hit[3]  # bitscore
        best_genome = best_hit[4]  # genome_name

        top_hit_identity = min(100.0, max(0.0, best_pident))
        novelty_index = 100.0 - top_hit_identity

        # Find second-best hit (to a DIFFERENT genome) for identity gap
        second_hit_identity: float | None = None
        identity_gap: float | None = None
        for hit in hits_data:
            if hit[4] != best_genome:  # genome_name
                second_hit_identity = hit[2]  # pident
                # Clamp to 0 - secondary hit can have higher identity than top bitscore hit
                identity_gap = max(0.0, top_hit_identity - second_hit_identity)
                break

        # Calculate placement uncertainty (95% threshold)
        bitscore_cutoff = best_bitscore * (bitscore_pct / 100.0)
        max_ani = 0.0
        second_genome_ani = 0.0  # ANI to second-best genome
        second_genome_found = False
        num_ambiguous_hits = 0

        best_idx = genome_to_idx.get(best_genome)
        uncertainty_mode = config_dict.get("uncertainty_mode", "max")

        for hit in hits_data:
            if hit[3] < bitscore_cutoff:  # bitscore
                break
            num_ambiguous_hits += 1
            secondary_genome = hit[4]  # genome_name
            if secondary_genome != best_genome:
                secondary_idx = genome_to_idx.get(secondary_genome)
                if best_idx is not None and secondary_idx is not None:
                    ani = float(ani_array[best_idx, secondary_idx])
                    # Use default_ani for missing values (0.0 means not computed)
                    if ani == 0.0:
                        ani = default_ani

                    # Track first secondary genome (second-best by bitscore)
                    if not second_genome_found:
                        second_genome_ani = ani
                        second_genome_found = True

                    if ani > max_ani:
                        max_ani = ani

        # Calculate genus-level uncertainty (90% threshold)
        genus_bitscore_cutoff = best_bitscore * (genus_bitscore_pct / 100.0)
        genus_max_ani = 0.0
        num_genus_hits = 0

        for hit in hits_data:
            if hit[3] < genus_bitscore_cutoff:  # bitscore
                break
            num_genus_hits += 1
            secondary_genome = hit[4]  # genome_name
            if secondary_genome != best_genome:
                secondary_idx = genome_to_idx.get(secondary_genome)
                if best_idx is not None and secondary_idx is not None:
                    ani = float(ani_array[best_idx, secondary_idx])
                    if ani == 0.0:
                        ani = default_ani
                    if ani > genus_max_ani:
                        genus_max_ani = ani

        # Calculate uncertainties based on mode
        if num_ambiguous_hits <= 1:
            placement_uncertainty = 0.0
        else:
            # Choose ANI based on uncertainty mode
            if uncertainty_mode == "second":
                effective_ani = second_genome_ani
            else:  # "max" mode (default)
                effective_ani = max_ani

            if effective_ani == 0.0:
                # Genome not in ANI matrix at all - maximum uncertainty
                placement_uncertainty = 100.0
            else:
                placement_uncertainty = 100.0 - effective_ani

        # Genus uncertainty
        genus_uncertainty: float | None = None
        if num_genus_hits > 1:
            if genus_max_ani == 0.0:
                genus_uncertainty = 100.0
            else:
                genus_uncertainty = 100.0 - genus_max_ani

        # Classify - Rule 0 (identity gap) first, then uncertainty-based checks
        # Rule 0: Identity gap check - when competing hits have nearly identical
        # BLAST identity (gap < 2%), the read cannot be confidently placed.
        if (
            identity_gap is not None
            and num_ambiguous_hits > 1
            and identity_gap < identity_gap_ambiguous_max
        ):
            call = "Ambiguous"
            is_novel = False
        # Rule 0b: Single-hit inferred uncertainty check
        # For single-hit reads, use inferred uncertainty to flag ambiguous cases.
        elif (
            use_inferred_for_single_hits
            and num_ambiguous_hits <= 1
            and novelty_index >= novelty_novel_species_min  # Only for novel range
        ):
            inferred = _calculate_inferred_uncertainty_inline(
                novelty_index, novelty_known_max, novelty_novel_species_max, novelty_novel_genus_max
            )
            if inferred >= single_hit_uncertainty_threshold:
                call = "Ambiguous"
                is_novel = False
            elif placement_uncertainty >= uncertainty_conserved_min:
                call = "Ambiguous"
                is_novel = False
            elif placement_uncertainty >= uncertainty_novel_genus_max:
                call = "Species Boundary"
                is_novel = False
            elif (
                novelty_novel_species_min <= novelty_index < novelty_novel_species_max
                and placement_uncertainty < uncertainty_novel_species_max
            ):
                call = "Novel Species"
                is_novel = True
            elif (
                novelty_novel_genus_min <= novelty_index <= novelty_novel_genus_max
                and placement_uncertainty < uncertainty_novel_genus_max
            ):
                if (
                    genus_uncertainty is not None
                    and num_genus_hits > 1
                    and genus_uncertainty >= genus_uncertainty_ambiguous_min
                ):
                    call = "Ambiguous Within Genus"
                    is_novel = False
                else:
                    call = "Novel Genus"
                    is_novel = True
            else:
                call = "Unclassified"
                is_novel = False
        elif placement_uncertainty >= uncertainty_conserved_min:
            # High uncertainty - conserved within genus (can't distinguish cross-genera here)
            call = "Ambiguous"
            is_novel = False
        elif placement_uncertainty >= uncertainty_novel_genus_max:
            # Moderate uncertainty (2-5%): species boundary zone
            call = "Species Boundary"
            is_novel = False
        elif (
            novelty_index < novelty_known_max
            and placement_uncertainty < uncertainty_known_max
        ):
            call = "Known Species"
            is_novel = False
        elif (
            novelty_novel_species_min <= novelty_index < novelty_novel_species_max
            and placement_uncertainty < uncertainty_novel_species_max
        ):
            call = "Novel Species"
            is_novel = True
        elif (
            novelty_novel_genus_min <= novelty_index <= novelty_novel_genus_max
            and placement_uncertainty < uncertainty_novel_genus_max
        ):
            # Check for genus-level ambiguity
            if (
                genus_uncertainty is not None
                and num_genus_hits > 1
                and genus_uncertainty >= genus_uncertainty_ambiguous_min
            ):
                call = "Ambiguous Within Genus"
                is_novel = False
            else:
                call = "Novel Genus"
                is_novel = True
        else:
            call = "Unclassified"
            is_novel = False

        # Calculate confidence score using margin-based approach
        # Pass mode-specific scaling parameters for protein vs nucleotide
        confidence_score = _calculate_confidence_score_inline(
            novelty_index=novelty_index,
            placement_uncertainty=placement_uncertainty,
            num_ambiguous_hits=num_ambiguous_hits,
            identity_gap=identity_gap,
            top_hit_identity=top_hit_identity,
            taxonomic_call=call,
            novelty_known_max=novelty_known_max,
            novelty_novel_species_max=novelty_novel_species_max,
            novelty_novel_genus_max=novelty_novel_genus_max,
            uncertainty_confident_max=uncertainty_known_max,
            # Mode-specific scaling parameters
            margin_divisor_known=margin_divisor_known,
            margin_divisor_novel_species=margin_divisor_novel_species,
            margin_divisor_novel_genus=margin_divisor_novel_genus,
            identity_gap_thresholds=identity_gap_thresholds,
            identity_score_base=identity_score_base,
            identity_score_range=identity_score_range,
        )

        result_dict = {
            "read_id": read_id,
            "best_match_genome": best_genome,
            "top_hit_identity": top_hit_identity,
            "novelty_index": novelty_index,
            "placement_uncertainty": placement_uncertainty,
            "num_ambiguous_hits": num_ambiguous_hits,
            "second_hit_identity": second_hit_identity,
            "identity_gap": identity_gap,
            "genus_uncertainty": genus_uncertainty,
            "num_genus_hits": num_genus_hits,
            "confidence_score": confidence_score,
            "taxonomic_call": call,
            "diversity_status": TAXONOMIC_TO_DIVERSITY[call],
            "is_novel": is_novel,
        }

        # Add enhanced scoring fields if enabled
        if infer_single_hit_uncertainty or enhanced_scoring:
            if num_ambiguous_hits <= 1:
                result_dict["inferred_uncertainty"] = _calculate_inferred_uncertainty_inline(
                    novelty_index, novelty_known_max, novelty_novel_species_max, novelty_novel_genus_max
                )
                result_dict["uncertainty_type"] = "inferred"
            else:
                result_dict["uncertainty_type"] = "measured"

        if enhanced_scoring:
            # Alignment quality (using basic coverage estimate from hit data)
            # Note: In parallel worker we don't have full BLAST fields, use defaults
            alignment_quality = 75.0  # Default quality for parallel processing
            result_dict["alignment_quality"] = alignment_quality

            identity_confidence = _calculate_identity_confidence_inline(
                novelty_index, alignment_quality
            )
            result_dict["identity_confidence"] = identity_confidence

            # Use inferred uncertainty for single-hit reads
            effective_uncertainty = result_dict.get("inferred_uncertainty", placement_uncertainty)
            effective_type = result_dict.get("uncertainty_type", "measured")

            placement_confidence = _calculate_placement_confidence_inline(
                effective_uncertainty, effective_type, identity_gap, num_ambiguous_hits
            )
            result_dict["placement_confidence"] = placement_confidence

            discovery_score = _calculate_discovery_score_inline(
                call, novelty_index, identity_confidence, placement_confidence, alignment_quality,
                novelty_novel_species_min, novelty_novel_genus_min
            )
            if discovery_score is not None:
                result_dict["discovery_score"] = discovery_score

        results.append(result_dict)

    return results


def _calculate_inferred_uncertainty_inline(
    novelty_index: float,
    novelty_known_max: float,
    novelty_novel_species_max: float,
    novelty_novel_genus_max: float,
) -> float:
    """
    Inline inferred uncertainty calculation for parallel workers.

    Duplicated from constants.calculate_inferred_uncertainty to avoid
    import issues with multiprocessing pickling.

    For single-hit reads, we cannot measure placement uncertainty from ANI
    between competing genomes. Instead, we infer uncertainty from novelty:
    higher novelty suggests either novel diversity OR database incompleteness.

    Args:
        novelty_index: Novelty index value (100 - top_hit_identity)
        novelty_known_max: Threshold for known species (default 5.0)
        novelty_novel_species_max: Threshold for novel species (default 20.0)
        novelty_novel_genus_max: Threshold for novel genus (default 25.0)

    Returns:
        Inferred uncertainty value (5-35%)
    """
    inferred_uncertainty_base = 5.0
    inferred_uncertainty_max = 35.0

    if novelty_index < novelty_known_max:
        # High identity: database likely complete for this species
        return inferred_uncertainty_base + novelty_index * 0.5  # 5-7.5%
    elif novelty_index < novelty_novel_species_max:
        # Novel species range: uncertain if truly novel or database gap
        return 7.5 + (novelty_index - novelty_known_max) * 1.0  # 7.5-22.5%
    elif novelty_index < novelty_novel_genus_max:
        # Novel genus range: high uncertainty
        return 22.5 + (novelty_index - novelty_novel_species_max) * 1.5  # 22.5-30%
    else:
        # Very high divergence: maximum uncertainty
        return inferred_uncertainty_max


def _calculate_confidence_score_inline(
    novelty_index: float,
    placement_uncertainty: float,
    num_ambiguous_hits: int,
    identity_gap: float | None,
    top_hit_identity: float,
    taxonomic_call: str,
    novelty_known_max: float,
    novelty_novel_species_max: float,
    novelty_novel_genus_max: float,
    uncertainty_confident_max: float,
    # Mode-specific scaling parameters (with nucleotide defaults)
    margin_divisor_known: float = 5.0,
    margin_divisor_novel_species: float = 7.5,
    margin_divisor_novel_genus: float = 5.0,
    identity_gap_thresholds: tuple[float, ...] = (5.0, 2.0, 1.0, 0.5),
    identity_score_base: float = 70.0,
    identity_score_range: float = 25.0,
) -> float:
    """
    Inline confidence score calculation for parallel workers.

    Duplicated from constants.calculate_confidence_score to avoid
    import issues with multiprocessing pickling.

    Supports both nucleotide and protein modes via parameterized scaling factors.
    """
    score = 0.0

    # Component 1: Margin from threshold boundaries (0-40 points)
    margin_score = 0.0

    if taxonomic_call == "Known Species":
        margin_from_boundary = novelty_known_max - novelty_index
        margin_score = min(40.0, (margin_from_boundary / margin_divisor_known) * 40.0)
    elif taxonomic_call == "Novel Species":
        margin_lower = novelty_index - novelty_known_max
        margin_upper = novelty_novel_species_max - novelty_index
        min_margin = min(margin_lower, margin_upper)
        margin_score = min(40.0, (min_margin / margin_divisor_novel_species) * 40.0)
    elif taxonomic_call == "Novel Genus":
        margin_lower = novelty_index - novelty_novel_species_max
        margin_upper = novelty_novel_genus_max - novelty_index
        min_margin = min(margin_lower, margin_upper)
        margin_score = min(40.0, (min_margin / margin_divisor_novel_genus) * 40.0)
    else:
        margin_score = 10.0

    if placement_uncertainty < uncertainty_confident_max:
        uncertainty_margin = uncertainty_confident_max - placement_uncertainty
        margin_score = min(
            40.0, margin_score + (uncertainty_margin / uncertainty_confident_max) * 10.0
        )

    score += max(0.0, margin_score)

    # Component 2: Placement certainty (0-40 points)
    placement_score = 0.0
    if num_ambiguous_hits <= 1:
        placement_score += 20.0
    elif num_ambiguous_hits <= 3:
        placement_score += 15.0
    elif num_ambiguous_hits <= 5:
        placement_score += 10.0
    elif num_ambiguous_hits <= 10:
        placement_score += 5.0

    # Use parameterized identity gap thresholds
    if identity_gap is not None:
        gap_high, gap_mid_high, gap_mid, gap_low = identity_gap_thresholds
        if identity_gap >= gap_high:
            placement_score += 20.0
        elif identity_gap >= gap_mid_high:
            placement_score += 15.0
        elif identity_gap >= gap_mid:
            placement_score += 10.0
        elif identity_gap >= gap_low:
            placement_score += 5.0
    else:
        if num_ambiguous_hits <= 1:
            placement_score += 20.0

    score += placement_score

    # Component 3: Alignment quality proxy (0-20 points)
    identity_score = min(
        20.0,
        max(0.0, (top_hit_identity - identity_score_base) / identity_score_range * 20.0)
    )
    score += identity_score

    return round(min(100.0, max(0.0, score)), 1)


def _calculate_identity_confidence_inline(
    novelty_index: float,
    alignment_quality: float,
) -> float:
    """
    Inline identity confidence calculation for parallel workers.

    Measures confidence in the identity measurement itself.
    """
    # Base score from novelty (inverted - higher identity = more confident)
    if novelty_index < 5:
        base = 80.0
    elif novelty_index < 15:
        base = 60.0 - (novelty_index - 5)
    elif novelty_index < 25:
        base = 50.0 - (novelty_index - 15) * 1.5
    else:
        base = 30.0

    # Alignment quality contribution (up to +20)
    quality_bonus = alignment_quality * 0.2

    return round(min(100.0, base + quality_bonus), 1)


def _calculate_placement_confidence_inline(
    uncertainty: float,
    uncertainty_type: str,
    identity_gap: float | None,
    num_ambiguous_hits: int,
) -> float:
    """
    Inline placement confidence calculation for parallel workers.

    Measures confidence in genome assignment.
    """
    # Base score from uncertainty level
    if uncertainty < 2:
        base = 80.0
    elif uncertainty < 5:
        base = 60.0
    elif uncertainty < 10:
        base = 40.0
    elif uncertainty < 20:
        base = 25.0
    else:
        base = 10.0

    # Penalty for inferred uncertainty
    if uncertainty_type == "inferred":
        base -= 15.0

    # Identity gap bonus (only meaningful for multi-hit)
    gap_bonus = 0.0
    if identity_gap is not None and num_ambiguous_hits > 1:
        if identity_gap >= 5:
            gap_bonus = 20.0
        elif identity_gap >= 2:
            gap_bonus = 10.0
        elif identity_gap >= 1:
            gap_bonus = 5.0

    return round(max(0.0, min(100.0, base + gap_bonus)), 1)


def _calculate_discovery_score_inline(
    taxonomic_call: str,
    novelty_index: float,
    identity_confidence: float,
    placement_confidence: float,
    alignment_quality: float,
    novelty_novel_species_min: float,
    novelty_novel_genus_min: float,
) -> float | None:
    """
    Inline discovery score calculation for parallel workers.

    Priority score for novel discoveries.
    """
    if taxonomic_call not in ("Novel Species", "Novel Genus"):
        return None

    # Novelty component (0-40)
    if taxonomic_call == "Novel Species":
        novelty_pts = 15.0 + (novelty_index - novelty_novel_species_min) * 1.5
    else:  # Novel Genus
        novelty_pts = 30.0 + (novelty_index - novelty_novel_genus_min) * 2.0
    novelty_pts = min(40.0, novelty_pts)

    # Quality component (0-30)
    quality_pts = alignment_quality * 0.3

    # Confidence component (0-30)
    conf_pts = (identity_confidence + placement_confidence) / 2.0 * 0.3

    return round(novelty_pts + quality_pts + conf_pts, 1)
