"""
Fully vectorized classifier using Polars for maximum performance.

This module provides VectorizedClassifier, which performs all classification
operations in Polars' native Rust backend for automatic parallelization
without multiprocessing overhead.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import polars as pl

from metadarkmatter.core.constants import calculate_confidence_score
from metadarkmatter.core.io_utils import write_dataframe, write_dataframe_append
from metadarkmatter.core.parsers import extract_genome_name_expr
from metadarkmatter.models.classification import TaxonomicCall, TAXONOMIC_TO_DIVERSITY
from metadarkmatter.models.config import ScoringConfig

if TYPE_CHECKING:
    from metadarkmatter.core.aai_matrix_builder import AAIMatrix
    from metadarkmatter.core.classification.ani_matrix import ANIMatrix
    from metadarkmatter.core.id_mapping import ContigIdMapping

class VectorizedClassifier:
    """
    Fully vectorized classifier using Polars for maximum performance.

    Unlike the iterator-based classifiers, this implementation performs
    all operations in Polars' native Rust backend, avoiding Python loops
    entirely. This provides automatic parallelization across CPU cores.

    Performance characteristics:
        - 5-10x faster than classify_to_dataframe_fast() for large files
        - Uses Polars' internal parallelism (no multiprocessing overhead)
        - Memory usage: ~1.5x input file size for in-memory operations
        - For 10M reads with 10 hits each: ~15GB memory
        - Use ParallelClassifier for memory-constrained systems

    Memory trade-offs:
        - Loads entire BLAST file into memory for vectorized operations
        - ANI lookup table: O(n^2) genome pairs, ~8MB for 1000 genomes
        - Symmetric lookup doubles memory but enables efficient joins

    Example:
        ani_matrix = ANIMatrix.from_file(ani_path)
        vectorized = VectorizedClassifier(ani_matrix, config)
        df = vectorized.classify_file(blast_path)
    """

    def __init__(
        self,
        ani_matrix: ANIMatrix,
        aai_matrix: AAIMatrix | None = None,
        config: ScoringConfig | None = None,
    ) -> None:
        """
        Initialize vectorized classifier.

        Args:
            ani_matrix: Precomputed ANI matrix
            aai_matrix: Optional AAI matrix for genus-level classification.
                Note: AAI integration in VectorizedClassifier is currently limited.
                Use ANIWeightedClassifier for full AAI-based genus classification.
            config: Scoring configuration (uses defaults if None).
                When alignment_mode is "protein", classification uses wider
                novelty thresholds appropriate for protein-level identity.
        """
        self.ani_matrix = ani_matrix
        self.aai_matrix = aai_matrix
        self.config = config or ScoringConfig()

        # Store effective thresholds based on alignment mode
        self._effective_thresholds = self.config.get_effective_thresholds()

        # Pre-build ANI lookup DataFrame for efficient joins
        self._build_ani_lookup()

        # Build AAI lookup if matrix provided
        if aai_matrix is not None:
            self._build_aai_lookup()
        else:
            self._aai_lookup = None
            self._aai_lookup_symmetric = None

    def _build_aai_lookup(self) -> None:
        """
        Build Polars DataFrame for vectorized AAI lookups.

        Similar to ANI lookup but for genus-level classification.
        """
        if self.aai_matrix is None:
            self._aai_lookup = None
            self._aai_lookup_symmetric = None
            return

        n = len(self.aai_matrix._genomes)
        if n == 0:
            self._aai_lookup = pl.DataFrame(
                schema={"genome1": pl.Utf8, "genome2": pl.Utf8, "aai": pl.Float64}
            )
            self._aai_lookup_symmetric = self._aai_lookup
            return

        # Get upper triangle indices using NumPy (vectorized)
        idx1, idx2 = np.triu_indices(n, k=1)

        # Extract genome names and AAI values using NumPy indexing
        genomes = np.array(self.aai_matrix._genomes)
        genome1_arr = genomes[idx1]
        genome2_arr = genomes[idx2]
        aai_arr = self.aai_matrix._aai_array[idx1, idx2].copy()

        # Replace missing AAI values (0.0) with default for distant organisms
        aai_arr[aai_arr == 0.0] = self.aai_matrix._default_aai

        # Build DataFrame directly from NumPy arrays (fast)
        self._aai_lookup = pl.DataFrame({
            "genome1": genome1_arr,
            "genome2": genome2_arr,
            "aai": aai_arr,
        })

        # Pre-compute symmetric version for faster joins
        self._aai_lookup_symmetric = pl.concat([
            self._aai_lookup,
            self._aai_lookup.select([
                pl.col("genome2").alias("genome1"),
                pl.col("genome1").alias("genome2"),
                pl.col("aai"),
            ]),
        ])

    def _build_ani_lookup(self) -> None:
        """
        Build Polars DataFrame for vectorized ANI lookups.

        Optimized using NumPy vectorization instead of Python loops.
        For 1000 genomes: 0.1s vs 10s with nested loops.
        """
        n = len(self.ani_matrix._genomes)
        if n == 0:
            self._ani_lookup = pl.DataFrame(
                schema={"genome1": pl.Utf8, "genome2": pl.Utf8, "ani": pl.Float64}
            )
            return

        # Get upper triangle indices using NumPy (vectorized)
        idx1, idx2 = np.triu_indices(n, k=1)

        # Extract genome names and ANI values using NumPy indexing
        genomes = np.array(self.ani_matrix._genomes)
        genome1_arr = genomes[idx1]
        genome2_arr = genomes[idx2]
        ani_arr = self.ani_matrix._ani_array[idx1, idx2].copy()

        # Replace missing ANI values (0.0) with default for distant organisms
        ani_arr[ani_arr == 0.0] = self.ani_matrix._default_ani

        # Build DataFrame directly from NumPy arrays (fast)
        self._ani_lookup = pl.DataFrame({
            "genome1": genome1_arr,
            "genome2": genome2_arr,
            "ani": ani_arr,
        })

        # Pre-compute symmetric version for faster joins
        self._ani_lookup_symmetric = pl.concat([
            self._ani_lookup,
            self._ani_lookup.select([
                pl.col("genome2").alias("genome1"),
                pl.col("genome1").alias("genome2"),
                pl.col("ani"),
            ]),
        ])

    def _coverage_weight_expr(
        self,
        coverage_col: pl.Expr,
        mode: str,
        strength: float,
    ) -> pl.Expr:
        """
        Generate Polars expression for coverage weighting.

        Args:
            coverage_col: Polars expression for coverage values (0.0-1.0)
            mode: Weighting mode ("linear", "log", "sigmoid")
            strength: Weight strength parameter (0.0-1.0)

        Returns:
            Polars expression that calculates coverage weight multiplier

        Raises:
            ValueError: If mode is unknown
        """
        min_weight = 1.0 - strength
        max_weight = 1.0 + strength
        weight_range = max_weight - min_weight

        if mode == "linear":
            normalized = coverage_col
        elif mode == "log":
            normalized = (1.0 + 9.0 * coverage_col).log() / np.log(10)
        elif mode == "sigmoid":
            normalized = 1.0 / (1.0 + (-10.0 * (coverage_col - 0.6)).exp())
        else:
            msg = f"Unknown coverage weight mode: {mode}"
            raise ValueError(msg)

        return min_weight + weight_range * normalized

    def classify_file(
        self,
        blast_path: Path,
        id_mapping: ContigIdMapping | None = None,
    ) -> pl.DataFrame:
        """
        Classify BLAST file using fully vectorized Polars operations.

        All operations run in Polars' Rust backend with automatic
        parallelization. No Python loops involved.

        Args:
            blast_path: Path to BLAST tabular output
            id_mapping: Optional ID mapping for external BLAST results.
                If provided, transforms sseqid values before genome extraction.

        Returns:
            DataFrame with classification results
        """
        from metadarkmatter.core.parsers import (
            StreamingBlastParser,
            extract_genome_name_expr,
        )

        parser = StreamingBlastParser(blast_path)

        # Read all data
        df = parser.parse_lazy().collect()

        # Apply ID transformation if mapping provided (for external BLAST results)
        if id_mapping is not None:
            df = id_mapping.transform_column(df, "sseqid")

        # Apply alignment quality filters (GTDB-compatible)
        if self.config.min_alignment_length > 0:
            df = df.filter(pl.col("length") >= self.config.min_alignment_length)

        # Note: alignment fraction filter requires read length which isn't
        # in standard BLAST output. Use qend - qstart + 1 as proxy.
        if self.config.min_alignment_fraction > 0:
            df = df.with_columns([
                ((pl.col("qend") - pl.col("qstart") + 1) / pl.col("qend")).alias("_approx_af")
            ]).filter(
                pl.col("_approx_af") >= self.config.min_alignment_fraction
            ).drop("_approx_af")

        # Add genome_name column (extracts accession from standardized sseqid)
        df = df.with_columns([extract_genome_name_expr()])

        if df.is_empty():
            return self._empty_dataframe()

        # Add read_length column for coverage calculation
        # Use qlen if available, otherwise fall back to qend as proxy
        if "qlen" in df.columns:
            df = df.with_columns([
                pl.when(pl.col("qlen").is_not_null())
                .then(pl.col("qlen"))
                .otherwise(pl.col("qend"))
                .alias("read_length")
            ])
        else:
            df = df.with_columns([
                pl.col("qend").alias("read_length")
            ])

        # Calculate coverage and weighted score if enabled
        if self.config.coverage_weight_mode != "none":
            df = df.with_columns([
                ((pl.col("qend") - pl.col("qstart") + 1) / pl.col("read_length"))
                .clip(0.0, 1.0)
                .alias("coverage"),
            ]).with_columns([
                self._coverage_weight_expr(
                    pl.col("coverage"),
                    self.config.coverage_weight_mode,
                    self.config.coverage_weight_strength,
                ).alias("coverage_weight"),
            ]).with_columns([
                (pl.col("bitscore") * pl.col("coverage_weight")).alias("weighted_bitscore"),
            ])
        else:
            # No weighting - use raw bitscore
            df = df.with_columns([
                pl.col("bitscore").alias("weighted_bitscore"),
            ])

        # Sort by weighted_bitscore DESC, pident DESC, genome_name ASC for deterministic tie-breaking
        # When multiple hits share max weighted_bitscore, this ensures consistent selection
        df = df.sort(["qseqid", "weighted_bitscore", "pident", "genome_name"], descending=[False, True, True, False])

        # Step 1: Find best hit per read (highest weighted_bitscore) and collect metrics for confidence
        best_hits = (
            df.group_by("qseqid", maintain_order=True)
            .agg([
                pl.col("weighted_bitscore").max().alias("max_weighted_bitscore"),
                # Keep original max_bitscore for compatibility
                pl.col("bitscore").max().alias("max_bitscore"),
                # Second-best weighted_bitscore for gap calculation
                pl.col("weighted_bitscore").get(1).alias("second_weighted_bitscore"),
                pl.col("pident")
                .filter(pl.col("weighted_bitscore") == pl.col("weighted_bitscore").max())
                .first()
                .alias("top_pident"),
                pl.col("genome_name")
                .filter(pl.col("weighted_bitscore") == pl.col("weighted_bitscore").max())
                .first()
                .alias("best_genome"),
                # Alignment length of best hit (for confidence calculation)
                pl.col("length")
                .filter(pl.col("weighted_bitscore") == pl.col("weighted_bitscore").max())
                .first()
                .alias("best_alignment_length"),
                # Approximate alignment fraction (qend proxy for read length)
                ((pl.col("qend") - pl.col("qstart") + 1) / pl.col("qend"))
                .filter(pl.col("weighted_bitscore") == pl.col("weighted_bitscore").max())
                .first()
                .alias("alignment_fraction"),
            ])
        )

        # Step 2: Calculate weighted_bitscore threshold per read
        threshold_pct = self.config.bitscore_threshold_pct / 100.0
        best_hits = best_hits.with_columns([
            (pl.col("max_weighted_bitscore") * threshold_pct).alias("weighted_bitscore_threshold"),
        ])

        # Step 3: Join back to find all ambiguous hits
        df_with_threshold = df.join(
            best_hits.select(["qseqid", "weighted_bitscore_threshold", "best_genome"]),
            on="qseqid",
            how="left",
        )

        # Filter to ambiguous hits only
        ambiguous = df_with_threshold.filter(
            pl.col("weighted_bitscore") >= pl.col("weighted_bitscore_threshold")
        )

        # Step 4: Count ambiguous hits per read
        hit_counts = ambiguous.group_by("qseqid").len().rename({"len": "num_ambiguous_hits"})

        # Step 4b: Calculate second_hit_identity (best pident to a DIFFERENT genome)
        # This is used for identity gap calculation to detect ambiguous placement
        secondary_hits = df_with_threshold.filter(
            pl.col("genome_name") != pl.col("best_genome")
        )
        if not secondary_hits.is_empty():
            # Get the best pident to a secondary genome per read
            second_hit_metrics = (
                secondary_hits
                .group_by("qseqid")
                .agg([
                    pl.col("pident").max().alias("second_hit_identity"),
                ])
            )
        else:
            second_hit_metrics = pl.DataFrame({
                "qseqid": [],
                "second_hit_identity": [],
            }).cast({
                "qseqid": pl.Utf8,
                "second_hit_identity": pl.Float64,
            })

        # Step 5: Find max ANI to secondary genomes and calculate phylogenetic context
        # Get secondary genomes (not the best genome)
        secondary = ambiguous.filter(pl.col("genome_name") != pl.col("best_genome"))

        # Join with ANI lookup (use pre-computed symmetric table)
        if not secondary.is_empty() and not self._ani_lookup.is_empty():
            secondary_with_ani = secondary.join(
                self._ani_lookup_symmetric,
                left_on=["best_genome", "genome_name"],
                right_on=["genome1", "genome2"],
                how="left",
            )

            # Max ANI per read (for genome-level placement uncertainty)
            max_ani_per_read = (
                secondary_with_ani
                .group_by("qseqid")
                .agg(pl.col("ani").max().fill_null(0.0).alias("max_ani"))
            )

            # Step 5b: Calculate phylogenetic context metrics
            # Count secondary genomes and their ANI relationship to best genome
            genus_threshold = 80.0  # Same genus if ANI >= 80%
            species_threshold = 95.0  # Same species if ANI >= 95%

            phylo_metrics = (
                secondary_with_ani
                .group_by("qseqid")
                .agg([
                    # Count distinct secondary genomes
                    pl.col("genome_name").n_unique().alias("num_secondary_genomes"),
                    # Count genomes in same genus as best hit (ANI >= 80%)
                    pl.col("ani").filter(pl.col("ani") >= genus_threshold).count().alias("same_genus_count"),
                    # Count genomes in same species as best hit (ANI >= 95%)
                    pl.col("ani").filter(pl.col("ani") >= species_threshold).count().alias("same_species_count"),
                    # Min ANI to secondary genomes (used for genus uncertainty)
                    pl.col("ani").min().fill_null(0.0).alias("min_ani_secondary"),
                ])
            )
        else:
            max_ani_per_read = pl.DataFrame({"qseqid": [], "max_ani": []}).cast({
                "qseqid": pl.Utf8, "max_ani": pl.Float64
            })
            phylo_metrics = pl.DataFrame({
                "qseqid": [],
                "num_secondary_genomes": [],
                "same_genus_count": [],
                "same_species_count": [],
                "min_ani_secondary": [],
            }).cast({
                "qseqid": pl.Utf8,
                "num_secondary_genomes": pl.UInt32,
                "same_genus_count": pl.UInt32,
                "same_species_count": pl.UInt32,
                "min_ani_secondary": pl.Float64,
            })

        # Step 6: Combine all metrics
        result = (
            best_hits
            .join(hit_counts, on="qseqid", how="left")
            .join(max_ani_per_read, on="qseqid", how="left")
            .join(phylo_metrics, on="qseqid", how="left")
            .join(second_hit_metrics, on="qseqid", how="left")
            .with_columns([
                pl.col("num_ambiguous_hits").fill_null(1),
                pl.col("max_ani").fill_null(0.0),
                pl.col("num_secondary_genomes").fill_null(0).cast(pl.Int64),
                pl.col("same_genus_count").fill_null(0).cast(pl.Int64),
                pl.col("same_species_count").fill_null(0).cast(pl.Int64),
                pl.col("min_ani_secondary").fill_null(0.0),
                # Fill nulls for confidence calculation inputs
                pl.col("second_weighted_bitscore").fill_null(0.0),
                pl.col("best_alignment_length").fill_null(0).cast(pl.Int64),
                pl.col("alignment_fraction").fill_null(0.0),
                # Second hit identity defaults to null (no secondary genome)
                pl.col("second_hit_identity"),
            ])
        )

        # Step 7: Calculate metrics
        # Use effective thresholds for protein vs nucleotide mode
        cfg = self.config
        eff = self._effective_thresholds
        result = result.with_columns([
            # Clamp pident
            pl.col("top_pident").clip(0.0, 100.0).alias("top_hit_identity"),
        ]).with_columns([
            # Novelty index
            (100.0 - pl.col("top_hit_identity")).alias("novelty_index"),
            # Identity gap: difference between best hit and second-best hit to DIFFERENT genome
            # Null if no secondary genome exists
            pl.when(pl.col("second_hit_identity").is_not_null())
            .then(pl.col("top_hit_identity") - pl.col("second_hit_identity"))
            .otherwise(pl.lit(None))
            .alias("identity_gap"),
            # Placement uncertainty (genome-level)
            pl.when(pl.col("num_ambiguous_hits") <= 1)
            .then(0.0)
            .when(pl.col("max_ani") == 0.0)
            .then(100.0)
            .otherwise(100.0 - pl.col("max_ani"))
            .alias("placement_uncertainty"),
            # Genus uncertainty: based on min ANI to secondary genomes
            # Low genus_uncertainty = all ambiguous hits in same genus
            # High genus_uncertainty = hits span multiple genera
            pl.when(pl.col("num_secondary_genomes") == 0)
            .then(0.0)  # No secondary genomes = no genus uncertainty
            .when(pl.col("min_ani_secondary") == 0.0)
            .then(100.0)  # No ANI data = max uncertainty
            .otherwise(100.0 - pl.col("min_ani_secondary"))
            .alias("genus_uncertainty"),
            # Number of competing genera: secondary genomes NOT in same genus as best hit
            # If all secondary genomes are in same genus (ANI >= 80%), competing_genera = 0
            (pl.col("num_secondary_genomes") - pl.col("same_genus_count"))
            .clip(0, None)
            .alias("num_competing_genera"),
        ]).with_columns([
            # Ambiguity scope: categorizes the phylogenetic scope of ambiguous hits
            # within_species: all secondary hits share >= 95% ANI with best genome
            # within_genus: all secondary hits share >= 80% ANI with best genome
            # across_genera: some secondary hits are from different genera
            pl.when(pl.col("num_secondary_genomes") == 0)
            .then(pl.lit("unambiguous"))  # Single genome hit
            .when(pl.col("same_species_count") == pl.col("num_secondary_genomes"))
            .then(pl.lit("within_species"))  # All hits same species (strain variants)
            .when(pl.col("same_genus_count") == pl.col("num_secondary_genomes"))
            .then(pl.lit("within_genus"))  # All hits same genus, confident genus call
            .otherwise(pl.lit("across_genera"))  # Hits span genera (conserved gene?)
            .alias("ambiguity_scope"),
        ])

        # Step 7b: Calculate confidence score (0-100)
        # Integrates multiple quality factors:
        # 1. Alignment quality (0-40 pts): length and coverage
        # 2. Placement certainty (0-40 pts): bitscore gap and hit count
        # 3. Phylogenetic context (0-20 pts): scope of ambiguity
        result = result.with_columns([
            # Bitscore gap: (best - second) / best * 100
            # Higher gap = more confident (single clear hit vs many similar hits)
            pl.when(pl.col("second_weighted_bitscore") == 0.0)
            .then(100.0)  # Only one hit = max gap
            .otherwise(
                ((pl.col("max_bitscore") - pl.col("second_weighted_bitscore")) / pl.col("max_bitscore") * 100.0)
                .clip(0.0, 100.0)
            )
            .alias("_bitscore_gap_pct"),
        ]).with_columns([
            # Component 1: Alignment quality score (0-40 points)
            # - Length: log-scaled, 100bp=20, 300bp=30, 500bp=35, 1000bp=40
            # - Fraction: linear, 0.5=20, 1.0=40
            (
                # Length component (0-40): log2(length/25) clamped to [0, 40]
                (pl.col("best_alignment_length").cast(pl.Float64).log(2.0) - 4.64)  # log2(25) ≈ 4.64
                .clip(0.0, 5.32)  # log2(1000/25) ≈ 5.32
                * (40.0 / 5.32)  # Scale to 0-40
            ).alias("_alignment_length_score"),
            # Fraction component (0-40): linear scaling
            (pl.col("alignment_fraction") * 40.0).clip(0.0, 40.0).alias("_alignment_frac_score"),
        ]).with_columns([
            # Combined alignment quality: average of length and fraction scores
            ((pl.col("_alignment_length_score") + pl.col("_alignment_frac_score")) / 2.0)
            .alias("_alignment_quality"),
        ]).with_columns([
            # Component 2: Placement certainty score (0-40 points)
            # - Bitscore gap: 0-20 pts (0%=0, 5%=10, 10%+=20)
            # - Hit count penalty: 1 hit=20, 2-5 hits=15, 6-10 hits=10, >10 hits=5
            (
                (pl.col("_bitscore_gap_pct") / 5.0).clip(0.0, 20.0)  # Gap score (0-20)
                + pl.when(pl.col("num_ambiguous_hits") <= 1).then(20.0)
                  .when(pl.col("num_ambiguous_hits") <= 5).then(15.0)
                  .when(pl.col("num_ambiguous_hits") <= 10).then(10.0)
                  .otherwise(5.0)  # Hit count score (5-20)
            ).alias("_placement_certainty"),
        ]).with_columns([
            # Component 3: Phylogenetic context score (0-20 points)
            # Rewards biologically expected ambiguity patterns
            pl.when(pl.col("ambiguity_scope") == "unambiguous").then(20.0)
            .when(pl.col("ambiguity_scope") == "within_species").then(18.0)  # Strain variation expected
            .when(pl.col("ambiguity_scope") == "within_genus").then(12.0)  # Genus confident
            .otherwise(5.0)  # across_genera = low confidence
            .alias("_phylo_context"),
        ]).with_columns([
            # Final confidence score: sum of all components (0-100)
            (
                pl.col("_alignment_quality")
                + pl.col("_placement_certainty")
                + pl.col("_phylo_context")
            ).clip(0.0, 100.0).round(1).alias("confidence_score"),
        ])

        # Step 8: Apply classification thresholds
        # Rule 0 (identity gap) comes first - when competing hits have nearly identical
        # BLAST identity (gap < 2%), the read cannot be confidently placed.
        # Uncertainty-based checks come next (ANI-based, more biologically meaningful)
        result = result.with_columns([
            # Rule 0: Identity gap check (before all other rules)
            # When two hits have nearly identical identity (gap < 2%), the read
            # cannot be confidently placed regardless of ANI between those genomes.
            pl.when(
                (pl.col("identity_gap").is_not_null())
                & (pl.col("num_ambiguous_hits") > 1)
                & (pl.col("identity_gap") < cfg.identity_gap_ambiguous_max)
            )
            .then(pl.lit("Ambiguous"))
            # Conserved Region: high uncertainty AND hits span multiple genera
            # This indicates a conserved gene present across genera (e.g., 16S, housekeeping)
            .when(
                (pl.col("placement_uncertainty") >= eff["uncertainty_conserved_min"])
                & (pl.col("ambiguity_scope") == "across_genera")
            )
            .then(pl.lit("Conserved Region"))
            # High uncertainty but NOT across genera = conserved within genus
            .when(pl.col("placement_uncertainty") >= eff["uncertainty_conserved_min"])
            .then(pl.lit("Ambiguous"))
            # Moderate uncertainty (2-5%): species boundary zone
            # Read matches multiple closely related species (95-98% ANI)
            .when(pl.col("placement_uncertainty") >= eff["uncertainty_novel_genus_max"])
            .then(pl.lit("Species Boundary"))
            .when(
                (pl.col("novelty_index") < eff["novelty_known_max"])
                & (pl.col("placement_uncertainty") < eff["uncertainty_known_max"])
            )
            .then(pl.lit("Known Species"))
            .when(
                (pl.col("novelty_index") >= eff["novelty_novel_species_min"])
                & (pl.col("novelty_index") < eff["novelty_novel_species_max"])
                & (pl.col("placement_uncertainty") < eff["uncertainty_novel_species_max"])
            )
            .then(pl.lit("Novel Species"))
            .when(
                (pl.col("novelty_index") >= eff["novelty_novel_genus_min"])
                & (pl.col("novelty_index") <= eff["novelty_novel_genus_max"])
                & (pl.col("placement_uncertainty") < eff["uncertainty_novel_genus_max"])
                # Check genus-level ambiguity: if genus_uncertainty >= threshold,
                # this indicates hits to multiple species with low ANI
                & (
                    (pl.col("genus_uncertainty") < cfg.genus_uncertainty_ambiguous_min)
                    | (pl.col("num_secondary_genomes") <= 0)
                )
            )
            .then(pl.lit("Novel Genus"))
            # Genus-level ambiguous: high novelty but competing species within genus
            .when(
                (pl.col("novelty_index") >= eff["novelty_novel_genus_min"])
                & (pl.col("novelty_index") <= eff["novelty_novel_genus_max"])
                & (pl.col("placement_uncertainty") < eff["uncertainty_novel_genus_max"])
                & (pl.col("genus_uncertainty") >= cfg.genus_uncertainty_ambiguous_min)
                & (pl.col("num_secondary_genomes") > 0)
            )
            .then(pl.lit("Ambiguous Within Genus"))
            .otherwise(pl.lit("Unclassified"))
            .alias("taxonomic_call")
        ])

        # Add diversity_status, is_novel, and low_confidence flags
        result = result.with_columns([
            pl.col("taxonomic_call").replace(TAXONOMIC_TO_DIVERSITY).alias("diversity_status"),
            pl.col("taxonomic_call").is_in(["Novel Species", "Novel Genus"]).alias("is_novel"),
            # Low confidence flag: indicates borderline classification quality
            # Useful for downstream filtering without changing the classification
            (pl.col("confidence_score") < cfg.confidence_threshold).alias("low_confidence"),
        ])

        # Select and rename final columns
        # Note: num_competing_genera removed (redundant with ambiguity_scope)
        return result.select([
            pl.col("qseqid").alias("read_id"),
            pl.col("best_genome").alias("best_match_genome"),
            "top_hit_identity",
            "novelty_index",
            "placement_uncertainty",
            "genus_uncertainty",
            "ambiguity_scope",
            "num_ambiguous_hits",
            "second_hit_identity",
            "identity_gap",
            "confidence_score",
            "taxonomic_call",
            "diversity_status",
            "is_novel",
            "low_confidence",
        ])

    def _empty_dataframe(self) -> pl.DataFrame:
        """Return empty DataFrame with correct schema."""
        return pl.DataFrame(
            schema={
                "read_id": pl.Utf8,
                "best_match_genome": pl.Utf8,
                "top_hit_identity": pl.Float64,
                "novelty_index": pl.Float64,
                "placement_uncertainty": pl.Float64,
                "genus_uncertainty": pl.Float64,
                "ambiguity_scope": pl.Utf8,
                "num_ambiguous_hits": pl.Int64,
                "second_hit_identity": pl.Float64,
                "identity_gap": pl.Float64,
                "confidence_score": pl.Float64,
                "taxonomic_call": pl.Utf8,
                "diversity_status": pl.Utf8,
                "is_novel": pl.Boolean,
                "low_confidence": pl.Boolean,
            }
        )

    def _classify_partition(self, df: pl.DataFrame) -> pl.DataFrame:
        """
        Classify a partition of BLAST hits.

        Internal method that processes a subset of data.
        Used by streaming methods to process in chunks.

        Args:
            df: DataFrame with BLAST hits (must have genome_name column)

        Returns:
            DataFrame with classification results
        """
        if df.is_empty():
            return self._empty_dataframe()

        # Add read_length column for coverage calculation
        if "qlen" in df.columns:
            df = df.with_columns([
                pl.when(pl.col("qlen").is_not_null())
                .then(pl.col("qlen"))
                .otherwise(pl.col("qend"))
                .alias("read_length")
            ])
        else:
            df = df.with_columns([
                pl.col("qend").alias("read_length")
            ])

        # Calculate coverage and weighted score if enabled
        if self.config.coverage_weight_mode != "none":
            df = df.with_columns([
                ((pl.col("qend") - pl.col("qstart") + 1) / pl.col("read_length"))
                .clip(0.0, 1.0)
                .alias("coverage"),
            ]).with_columns([
                self._coverage_weight_expr(
                    pl.col("coverage"),
                    self.config.coverage_weight_mode,
                    self.config.coverage_weight_strength,
                ).alias("coverage_weight"),
            ]).with_columns([
                (pl.col("bitscore") * pl.col("coverage_weight")).alias("weighted_bitscore"),
            ])
        else:
            df = df.with_columns([
                pl.col("bitscore").alias("weighted_bitscore"),
            ])

        # Sort for deterministic tie-breaking when multiple hits share max weighted_bitscore
        df = df.sort(["qseqid", "weighted_bitscore", "pident", "genome_name"], descending=[False, True, True, False])

        # Step 1: Find best hit per read
        best_hits = (
            df.group_by("qseqid", maintain_order=True)
            .agg([
                pl.col("weighted_bitscore").max().alias("max_weighted_bitscore"),
                pl.col("bitscore").max().alias("max_bitscore"),
                pl.col("pident").filter(
                    pl.col("weighted_bitscore") == pl.col("weighted_bitscore").max()
                ).first().alias("top_pident"),
                pl.col("genome_name").filter(
                    pl.col("weighted_bitscore") == pl.col("weighted_bitscore").max()
                ).first().alias("best_genome"),
            ])
        )

        # Step 2: Weighted bitscore threshold
        threshold_pct = self.config.bitscore_threshold_pct / 100.0
        best_hits = best_hits.with_columns([
            (pl.col("max_weighted_bitscore") * threshold_pct).alias("weighted_bitscore_threshold"),
        ])

        # Step 3: Find ambiguous hits
        df_with_threshold = df.join(
            best_hits.select(["qseqid", "weighted_bitscore_threshold", "best_genome"]),
            on="qseqid",
            how="left",
        )
        ambiguous = df_with_threshold.filter(
            pl.col("weighted_bitscore") >= pl.col("weighted_bitscore_threshold")
        )

        # Step 4: Count ambiguous hits
        hit_counts = ambiguous.group_by("qseqid").len().rename({"len": "num_ambiguous_hits"})

        # Step 4b: Calculate second_hit_identity (best pident to a DIFFERENT genome)
        secondary_all = df_with_threshold.filter(
            pl.col("genome_name") != pl.col("best_genome")
        )
        if not secondary_all.is_empty():
            second_hit_metrics = (
                secondary_all
                .group_by("qseqid")
                .agg([
                    pl.col("pident").max().alias("second_hit_identity"),
                ])
            )
        else:
            second_hit_metrics = pl.DataFrame({
                "qseqid": [],
                "second_hit_identity": [],
            }).cast({
                "qseqid": pl.Utf8,
                "second_hit_identity": pl.Float64,
            })

        # Step 5: Find max ANI to secondary genomes
        secondary = ambiguous.filter(pl.col("genome_name") != pl.col("best_genome"))

        has_ani_lookup = (
            hasattr(self, '_ani_lookup_symmetric')
            and not self._ani_lookup_symmetric.is_empty()
        )
        if not secondary.is_empty() and has_ani_lookup:
            secondary_with_ani = secondary.join(
                self._ani_lookup_symmetric,
                left_on=["best_genome", "genome_name"],
                right_on=["genome1", "genome2"],
                how="left",
            )
            max_ani_per_read = (
                secondary_with_ani
                .group_by("qseqid")
                .agg(pl.col("ani").max().fill_null(0.0).alias("max_ani"))
            )
        else:
            max_ani_per_read = pl.DataFrame({"qseqid": [], "max_ani": []}).cast({
                "qseqid": pl.Utf8, "max_ani": pl.Float64
            })

        # Step 6: Combine metrics
        result = (
            best_hits
            .join(hit_counts, on="qseqid", how="left")
            .join(max_ani_per_read, on="qseqid", how="left")
            .join(second_hit_metrics, on="qseqid", how="left")
            .with_columns([
                pl.col("num_ambiguous_hits").fill_null(1),
                pl.col("max_ani").fill_null(0.0),
                pl.col("second_hit_identity"),
            ])
        )

        # Step 7: Calculate metrics and classify
        # Use effective thresholds for protein vs nucleotide mode
        eff = self._effective_thresholds
        result = (
            result
            .with_columns([
                pl.col("top_pident").clip(0.0, 100.0).alias("top_hit_identity"),
            ])
            .with_columns([
                (100.0 - pl.col("top_hit_identity")).alias("novelty_index"),
                pl.when(pl.col("num_ambiguous_hits") <= 1)
                .then(0.0)
                .when(pl.col("max_ani") == 0.0)
                .then(100.0)
                .otherwise(100.0 - pl.col("max_ani"))
                .alias("placement_uncertainty"),
                # Identity gap: difference between best hit and second-best hit to DIFFERENT genome
                pl.when(pl.col("second_hit_identity").is_not_null())
                .then(pl.col("top_hit_identity") - pl.col("second_hit_identity"))
                .otherwise(pl.lit(None))
                .alias("identity_gap"),
            ])
            .with_columns([
                # Uncertainty-based checks first (ANI-based, more biologically meaningful)
                # High uncertainty - conserved within genus (no ambiguity_scope in fast mode)
                pl.when(pl.col("placement_uncertainty") >= eff["uncertainty_conserved_min"])
                .then(pl.lit("Ambiguous"))
                # Moderate uncertainty (2-5%): species boundary zone
                .when(pl.col("placement_uncertainty") >= eff["uncertainty_novel_genus_max"])
                .then(pl.lit("Species Boundary"))
                .when(
                    (pl.col("novelty_index") < eff["novelty_known_max"])
                    & (pl.col("placement_uncertainty") < eff["uncertainty_known_max"])
                )
                .then(pl.lit("Known Species"))
                .when(
                    (pl.col("novelty_index") >= eff["novelty_novel_species_min"])
                    & (pl.col("novelty_index") < eff["novelty_novel_species_max"])
                    & (pl.col("placement_uncertainty") < eff["uncertainty_novel_species_max"])
                )
                .then(pl.lit("Novel Species"))
                .when(
                    (pl.col("novelty_index") >= eff["novelty_novel_genus_min"])
                    & (pl.col("novelty_index") <= eff["novelty_novel_genus_max"])
                    & (pl.col("placement_uncertainty") < eff["uncertainty_novel_genus_max"])
                )
                .then(pl.lit("Novel Genus"))
                .otherwise(pl.lit("Unclassified"))
                .alias("taxonomic_call")
            ])
            .with_columns([
                pl.col("taxonomic_call").replace(TAXONOMIC_TO_DIVERSITY).alias("diversity_status"),
                pl.col("taxonomic_call").is_in(["Novel Species", "Novel Genus"]).alias("is_novel")
            ])
        )

        return result.select([
            pl.col("qseqid").alias("read_id"),
            pl.col("best_genome").alias("best_match_genome"),
            "top_hit_identity",
            "novelty_index",
            "placement_uncertainty",
            "num_ambiguous_hits",
            "second_hit_identity",
            "identity_gap",
            "taxonomic_call",
            "diversity_status",
            "is_novel",
        ])

    def stream_to_file(
        self,
        blast_path: Path,
        output_path: Path,
        output_format: str = "parquet",
        partition_size: int = 5_000_000,
        progress_callback: Callable[[int, int, float], None] | None = None,
    ) -> int:
        """
        Stream classification results to file for very large BLAST files.

        Processes the input file in partitions to maintain bounded memory
        usage. Suitable for files with 100M+ alignments.

        Memory usage: ~2GB per 5M alignments partition (default).

        Args:
            blast_path: Path to BLAST tabular output
            output_path: Output file path
            output_format: 'parquet' or 'csv'
            partition_size: Number of alignments per partition (default: 5M)
            progress_callback: Optional callback(rows_processed, total_reads, elapsed_secs)

        Returns:
            Total number of reads classified

        Example:
            vectorized = VectorizedClassifier(ani_matrix, config)

            def progress(rows, reads, elapsed):
                print(f"Processed {rows:,} rows -> {reads:,} reads in {elapsed:.1f}s")

            total = vectorized.stream_to_file(
                blast_path, output_path,
                progress_callback=progress
            )
        """
        import time

        from metadarkmatter.core.parsers import extract_genome_name_expr

        start_time = time.time()
        total_reads = 0
        total_rows = 0
        first_partition = True

        # Use scan_csv with collect_batches for streaming
        # Detect column count from file (12 for legacy, 13 for new format with qlen)
        from metadarkmatter.core.parsers import StreamingBlastParser
        parser = StreamingBlastParser(blast_path)

        batches = pl.scan_csv(
            blast_path,
            separator="\t",
            has_header=False,
            new_columns=parser.column_names,
        ).collect_batches(chunk_size=partition_size)

        # Track reads across partitions to handle boundary cases
        pending_reads: dict[str, list] = {}

        for partition_df in batches:
            if partition_df.is_empty():
                continue

            total_rows += len(partition_df)

            # Add genome_name column
            partition_df = partition_df.with_columns([extract_genome_name_expr()])

            # Get unique reads in this partition
            partition_reads = set(partition_df["qseqid"].unique().to_list())

            # Find reads that are complete (not the last read which may span partitions)
            if partition_df.height > 0:
                last_read = partition_df["qseqid"][-1]
                complete_reads = partition_reads - {last_read}

                # Filter to complete reads only
                if complete_reads:
                    complete_df = partition_df.filter(pl.col("qseqid").is_in(complete_reads))

                    # Add any pending data from previous partition
                    for read_id in list(pending_reads.keys()):
                        if read_id in complete_reads:
                            pending_data = pending_reads.pop(read_id)
                            pending_df = pl.DataFrame(pending_data)
                            complete_df = pl.concat([complete_df, pending_df])

                    # Classify this partition
                    if not complete_df.is_empty():
                        result_df = self._classify_partition(complete_df)
                        total_reads += len(result_df)

                        # Write results
                        self._write_partition(
                            result_df, output_path, output_format, first_partition
                        )
                        first_partition = False

                # Save incomplete reads for next partition
                incomplete_df = partition_df.filter(pl.col("qseqid") == last_read)
                if not incomplete_df.is_empty():
                    if last_read not in pending_reads:
                        pending_reads[last_read] = []
                    pending_reads[last_read].extend(incomplete_df.to_dicts())

            if progress_callback:
                elapsed = time.time() - start_time
                progress_callback(total_rows, total_reads, elapsed)

        # Process any remaining pending reads
        if pending_reads:
            remaining_data = []
            for rows in pending_reads.values():
                remaining_data.extend(rows)

            if remaining_data:
                remaining_df = pl.DataFrame(remaining_data)
                result_df = self._classify_partition(remaining_df)
                total_reads += len(result_df)
                self._write_partition(
                    result_df, output_path, output_format, first_partition
                )

        if progress_callback:
            elapsed = time.time() - start_time
            progress_callback(total_rows, total_reads, elapsed)

        return total_reads

    def _write_partition(
        self,
        df: pl.DataFrame,
        output_path: Path,
        output_format: str,
        is_first: bool,
    ) -> None:
        """Write a partition of results to file."""
        write_dataframe_append(df, output_path, output_format, is_first)


