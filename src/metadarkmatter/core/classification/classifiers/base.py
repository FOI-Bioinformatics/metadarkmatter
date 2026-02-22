"""
Base classifier implementation using ANI-weighted placement.

This module provides the ANIWeightedClassifier, the core classification
algorithm that combines BLAST alignment scores with ANI (Average Nucleotide
Identity) distances to identify novel bacterial taxa.
"""

from __future__ import annotations

from collections.abc import Iterator
from typing import TYPE_CHECKING

from metadarkmatter.core.classification.bayesian import (
    BayesianClassifier,
    apply_stage2_refinement,
    entropy_to_confidence,
)
from metadarkmatter.core.constants import (
    calculate_alignment_quality,
    calculate_confidence_score,
    calculate_discovery_score,
    calculate_identity_confidence,
    calculate_inferred_uncertainty,
    calculate_placement_confidence,
)
from metadarkmatter.models.blast import BlastResult
from metadarkmatter.models.classification import (
    ReadClassification,
    TaxonomicCall,
)
from metadarkmatter.models.config import ScoringConfig

# Map from string category names to TaxonomicCall enum members
_CATEGORY_TO_ENUM: dict[str, TaxonomicCall] = {tc.value: tc for tc in TaxonomicCall}

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

        # Bayesian classifier for primary classification
        self._bayesian = BayesianClassifier(self.config)

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

        # Primary classification via Bayesian posteriors
        posterior = self._bayesian.compute_posteriors(
            novelty_index=novelty_index,
            placement_uncertainty=placement_uncertainty,
            num_hits=num_ambiguous_hits,
            identity_gap=identity_gap,
        )

        # Stage 2 discrete refinement (scalar version)
        import numpy as np
        refined_array = apply_stage2_refinement(
            [posterior.map_category],
            genus_uncertainty=np.array([genus_uncertainty or 0.0]),
            num_secondary=(
                np.array([num_genus_hits - 1]) if num_genus_hits is not None and num_genus_hits > 1
                else np.array([0])
            ),
            genus_uncertainty_threshold=self.config.genus_uncertainty_ambiguous_min,
        )
        bayesian_call = str(refined_array[0])

        # AAI-based override for genus-level classification
        # AAI is more reliable than ANI at high divergence (>20% novelty)
        eff = self._effective_thresholds
        if (
            aai_uncertainty is not None
            and self.config.use_aai_for_genus
            and eff["novelty_novel_genus_min"] <= novelty_index <= eff["novelty_novel_genus_max"]
            and placement_uncertainty < eff["uncertainty_novel_genus_max"]
        ):
            aai_genus_uncertainty_low = 100.0 - self.config.aai_genus_boundary_high
            aai_genus_uncertainty_high = 100.0 - self.config.aai_genus_boundary_low
            if aai_uncertainty < aai_genus_uncertainty_low:
                bayesian_call = "Novel Species"
            elif aai_uncertainty <= aai_genus_uncertainty_high:
                bayesian_call = "Ambiguous Within Genus"
            else:
                bayesian_call = "Novel Genus"

        taxonomic_call = _CATEGORY_TO_ENUM[bayesian_call]

        # Confidence score from posterior entropy
        confidence_score = float(entropy_to_confidence(posterior.entropy))

        # Initialize enhanced scoring fields
        inferred_uncertainty: float | None = None
        uncertainty_type: str | None = None
        alignment_quality: float | None = None
        identity_confidence: float | None = None
        placement_confidence: float | None = None
        discovery_score: float | None = None

        # Calculate inferred uncertainty for single-hit reads
        if num_ambiguous_hits <= 1:
            inferred_uncertainty = calculate_inferred_uncertainty(novelty_index)
            uncertainty_type = "inferred"
        else:
            uncertainty_type = "measured"

        # Calculate enhanced scoring metrics
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
            num_ambiguous_hits <= 1
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

