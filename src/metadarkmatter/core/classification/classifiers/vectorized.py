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

from collections.abc import Callable

from metadarkmatter.core.classification.bayesian import (
    BayesianClassifier,
    apply_stage2_refinement,
    entropy_to_confidence,
)
from metadarkmatter.core.constants import (
    UNKNOWN_GENOME,
    calculate_confidence_score,
)
from metadarkmatter.core.io_utils import write_dataframe, write_dataframe_append
from metadarkmatter.core.parsers import extract_genome_name_expr
from metadarkmatter.models.classification import TaxonomicCall, TAXONOMIC_TO_DIVERSITY
from metadarkmatter.models.config import ScoringConfig

if TYPE_CHECKING:
    from metadarkmatter.core.aai_matrix_builder import AAIMatrix
    from metadarkmatter.core.classification.ani_matrix import ANIMatrix
    from metadarkmatter.core.classification.qc import QCMetrics
    from metadarkmatter.core.id_mapping import ContigIdMapping

class VectorizedClassifier:
    """
    Fully vectorized classifier using Polars for maximum performance.

    Unlike the iterator-based classifiers, this implementation performs
    all operations in Polars' native Rust backend, avoiding Python loops
    entirely. This provides automatic parallelization across CPU cores.

    Performance characteristics:
        - Optimized for large files with vectorized Polars operations
        - Uses Polars' internal parallelism (no multiprocessing overhead)
        - Memory usage: ~1.5x input file size for in-memory operations
        - For 10M reads with 10 hits each: ~15GB memory
        - Consider chunked processing for memory-constrained systems

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
        representative_mapping: dict[str, str] | None = None,
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
            representative_mapping: Optional mapping from genome accession to
                species representative accession. When provided, ANI lookups
                use representative genomes while output preserves actual hit
                genomes. This decouples alignment (all genomes) from
                classification (representative ANI matrix).
        """
        self.ani_matrix = ani_matrix
        self.aai_matrix = aai_matrix
        self.config = config or ScoringConfig()
        self.representative_mapping = representative_mapping

        # Store effective thresholds based on alignment mode
        self._effective_thresholds = self.config.get_effective_thresholds()

        # Pre-compute inferred uncertainty breakpoints for continuous piecewise formula
        eff = self._effective_thresholds
        self._inferred_known_break = 5.0 + eff["novelty_known_max"] * 0.5
        self._inferred_species_break = (
            self._inferred_known_break
            + (eff["novelty_novel_species_max"] - eff["novelty_known_max"]) * 1.0
        )

        # Initialize Bayesian classifier for primary classification
        self._bayesian = BayesianClassifier(self.config)

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
        blast_input: Path | pl.DataFrame,
        id_mapping: ContigIdMapping | None = None,
        compute_qc: bool = False,
    ) -> pl.DataFrame | tuple[pl.DataFrame, "QCMetrics"]:
        """
        Classify BLAST hits using fully vectorized Polars operations.

        All operations run in Polars' Rust backend with automatic
        parallelization. No Python loops involved.

        Args:
            blast_input: Path to BLAST tabular output, or a pre-loaded
                DataFrame with standard BLAST columns and genome_name.
            id_mapping: Optional ID mapping for external BLAST results.
                If provided, transforms sseqid values before genome extraction.
                Ignored when blast_input is a DataFrame.
            compute_qc: If True, return (DataFrame, QCMetrics) tuple.

        Returns:
            DataFrame with classification results, or (DataFrame, QCMetrics)
            tuple if compute_qc is True.
        """
        if isinstance(blast_input, pl.DataFrame):
            # DataFrame passed directly (e.g. from stream_to_file partitions).
            # Assumes genome_name column and alignment filters already applied.
            df = blast_input
            raw_df = df if compute_qc else None
        else:
            from metadarkmatter.core.parsers import (
                StreamingBlastParser,
                extract_genome_name_expr,
            )

            parser = StreamingBlastParser(blast_input)

            # Read all data
            df = parser.parse_lazy().collect()

            # Apply ID transformation if mapping provided (for external BLAST results)
            if id_mapping is not None:
                df = id_mapping.transform_column(df, "sseqid")

            # Capture raw count before filtering (for QC)
            raw_df = df if compute_qc else None

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

        # Family validation: partition hits by ANI matrix membership
        family_validation_active = self.config.target_family is not None
        family_metrics_df = None
        off_target_df = None

        if family_validation_active:
            if self.representative_mapping:
                in_family_genomes = set(self.representative_mapping.keys())
            else:
                in_family_genomes = set(self.ani_matrix._genomes)

            df = df.with_columns([
                pl.col("genome_name").is_in(in_family_genomes).alias("_is_in_family"),
            ])

            # Compute external best genome per read separately (avoids sort_by
            # length mismatch when filtering within group_by aggregation)
            external_hits = df.filter(~pl.col("_is_in_family"))
            if not external_hits.is_empty():
                external_best_per_read = (
                    external_hits
                    .sort(["qseqid", "bitscore"], descending=[False, True])
                    .group_by("qseqid")
                    .agg([
                        pl.col("genome_name").first().alias("external_best_genome"),
                        pl.col("sseqid").first().alias("_ext_sseqid"),
                    ])
                    # Use raw sseqid as fallback when genome_name is "unknown"
                    # (happens with external alignments where unmapped contigs
                    # have no ID mapping entry)
                    .with_columns(
                        pl.when(pl.col("external_best_genome") == UNKNOWN_GENOME)
                        .then(pl.col("_ext_sseqid"))
                        .otherwise(pl.col("external_best_genome"))
                        .alias("external_best_genome")
                    )
                    .drop("_ext_sseqid")
                )
            else:
                external_best_per_read = pl.DataFrame(
                    schema={"qseqid": pl.Utf8, "external_best_genome": pl.Utf8}
                )

            # Per-read family metrics
            read_metrics = (
                df.group_by("qseqid")
                .agg([
                    pl.col("bitscore").filter(pl.col("_is_in_family")).max().fill_null(0.0).alias("_best_if_bs"),
                    pl.col("bitscore").max().alias("_best_all_bs"),
                    pl.col("pident").filter(~pl.col("_is_in_family")).max().fill_null(0.0).alias("external_best_identity"),
                    pl.col("pident").filter(pl.col("_is_in_family")).max().fill_null(0.0).alias("_best_if_ident"),
                    pl.col("_is_in_family").sum().alias("_if_count"),
                    pl.len().alias("_total_count"),
                ])
                .join(external_best_per_read, on="qseqid", how="left")
                .with_columns([
                    (pl.col("_best_if_bs") / pl.col("_best_all_bs")).fill_nan(0.0).fill_null(0.0).alias("family_bitscore_ratio"),
                    (pl.col("_best_if_ident") - pl.col("external_best_identity")).alias("family_identity_gap"),
                    (pl.col("_if_count").cast(pl.Float64) / pl.col("_total_count").cast(pl.Float64)).alias("in_family_hit_fraction"),
                ])
            )

            family_metrics_df = read_metrics

            threshold = self.config.family_ratio_threshold
            off_target_reads = read_metrics.filter(
                (pl.col("family_bitscore_ratio") < threshold)
                | (pl.col("_best_if_bs") == 0.0)
            )

            if not off_target_reads.is_empty():
                off_target_read_ids = set(off_target_reads["qseqid"].to_list())

                off_target_cols = [
                    pl.col("qseqid").alias("read_id"),
                    pl.col("external_best_genome").fill_null("unknown").alias("best_match_genome"),
                    pl.col("external_best_identity").alias("top_hit_identity"),
                    (100.0 - pl.col("external_best_identity")).alias("novelty_index"),
                    pl.lit(0.0).alias("placement_uncertainty"),
                    pl.lit(0.0).alias("genus_uncertainty"),
                    pl.lit("unambiguous").alias("ambiguity_scope"),
                    pl.lit(0).cast(pl.Int64).alias("num_ambiguous_hits"),
                    pl.lit(None).cast(pl.Float64).alias("second_hit_identity"),
                    pl.lit(None).cast(pl.Float64).alias("identity_gap"),
                    pl.lit(0.0).alias("confidence_score"),
                    pl.lit("Off-target").alias("taxonomic_call"),
                    pl.lit("Uncertain").alias("diversity_status"),
                    pl.lit(False).alias("is_novel"),
                    pl.lit(False).alias("low_confidence"),
                    pl.lit(None).cast(pl.Float64).alias("inferred_uncertainty"),
                    pl.lit("none").alias("uncertainty_type"),
                    pl.col("family_bitscore_ratio"),
                    pl.col("family_identity_gap"),
                    pl.col("in_family_hit_fraction"),
                    pl.col("external_best_genome"),
                    pl.col("external_best_identity"),
                    pl.col("_best_if_ident").alias("in_family_identity"),
                    (100.0 - pl.col("_best_if_ident")).alias("in_family_novelty_index"),
                ]

                # Only include legacy score columns when they are
                # present in the classified reads (include_legacy_scores).
                if self.config.include_legacy_scores:
                    off_target_cols.extend([
                        pl.lit(0.0).alias("alignment_quality"),
                        pl.lit(0.0).alias("identity_confidence"),
                        pl.lit(0.0).alias("placement_confidence"),
                        pl.lit(None).cast(pl.Float64).alias("discovery_score"),
                    ])

                off_target_df = off_target_reads.select(off_target_cols)

                # Keep only in-family hits for non-off-target reads
                df = df.filter(
                    ~pl.col("qseqid").is_in(off_target_read_ids)
                    & pl.col("_is_in_family")
                )
            else:
                # No off-target reads - still filter to in-family hits only
                df = df.filter(pl.col("_is_in_family"))

            df = df.drop("_is_in_family")

        # Capture filtered state for QC before empty check
        filtered_df = df if compute_qc else None

        if df.is_empty():
            # If family validation produced off-target reads, return those
            if family_validation_active and off_target_df is not None and not off_target_df.is_empty():
                if compute_qc:
                    from metadarkmatter.core.classification.qc import (
                        QCMetrics,
                        compute_pre_qc,
                    )
                    qc = compute_pre_qc(raw_df, filtered_df, self.ani_matrix, self.config)
                    return off_target_df, qc
                return off_target_df
            empty = self._empty_dataframe()
            if compute_qc:
                from metadarkmatter.core.classification.qc import (
                    QCMetrics,
                    compute_pre_qc,
                )
                qc = compute_pre_qc(raw_df, filtered_df, self.ani_matrix, self.config)
                return empty, qc
            return empty

        # Add representative_genome column for ANI lookups
        # When a representative mapping is provided, ANI lookups use the
        # representative genome while output preserves the actual hit genome
        if self.representative_mapping:
            df = df.with_columns(
                pl.col("genome_name")
                .replace_strict(
                    self.representative_mapping,
                    default=pl.col("genome_name"),
                )
                .alias("representative_genome")
            )
        else:
            df = df.with_columns(
                pl.col("genome_name").alias("representative_genome")
            )

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
                # Representative of best genome (for ANI lookups)
                pl.col("representative_genome")
                .filter(pl.col("weighted_bitscore") == pl.col("weighted_bitscore").max())
                .first()
                .alias("best_representative"),
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
            best_hits.select(["qseqid", "weighted_bitscore_threshold", "best_genome", "best_representative"]),
            on="qseqid",
            how="left",
        )

        # Filter to ambiguous hits only
        ambiguous = df_with_threshold.filter(
            pl.col("weighted_bitscore") >= pl.col("weighted_bitscore_threshold")
        )

        # Step 4: Count ambiguous hits per read
        hit_counts = ambiguous.group_by("qseqid").len().rename({"len": "num_ambiguous_hits"})

        # Step 4b: Calculate second_hit_identity (best pident to a DIFFERENT representative)
        # This is used for identity gap calculation to detect ambiguous placement
        # Use representative_genome to collapse hits to same species
        secondary_hits = df_with_threshold.filter(
            pl.col("representative_genome") != pl.col("best_representative")
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
        # Get secondary genomes (not the best representative - collapses same-species hits)
        secondary = ambiguous.filter(pl.col("representative_genome") != pl.col("best_representative"))

        # Join with ANI lookup using representative genomes (use pre-computed symmetric table)
        if not secondary.is_empty() and not self._ani_lookup.is_empty():
            secondary_with_ani = secondary.join(
                self._ani_lookup_symmetric,
                left_on=["best_representative", "representative_genome"],
                right_on=["genome1", "genome2"],
                how="left",
            )

            # Max ANI per read (for genome-level placement uncertainty - "max" mode)
            max_ani_per_read = (
                secondary_with_ani
                .group_by("qseqid")
                .agg(pl.col("ani").max().fill_null(0.0).alias("max_ani"))
            )

            # Second genome ANI per read (for "second" mode)
            # Get ANI to the representative with highest bitscore among secondaries
            # Filter out same-representative hits since ANI lookup
            # doesn't have diagonal entries
            second_genome_ani_per_read = (
                secondary_with_ani
                .filter(pl.col("representative_genome") != pl.col("best_representative"))
                .sort(["qseqid", "bitscore"], descending=[False, True])
                .group_by("qseqid")
                .agg(pl.col("ani").first().fill_null(0.0).alias("second_genome_ani"))
            )

            # Join second_genome_ani to max_ani_per_read
            max_ani_per_read = max_ani_per_read.join(
                second_genome_ani_per_read, on="qseqid", how="left"
            ).with_columns([
                pl.col("second_genome_ani").fill_null(0.0),
            ])

            # Step 5b: Calculate phylogenetic context metrics
            # Count secondary representatives and their ANI relationship to best genome
            genus_threshold = 80.0  # Same genus if ANI >= 80%
            species_threshold = 95.0  # Same species if ANI >= 95%

            phylo_metrics = (
                secondary_with_ani
                .group_by("qseqid")
                .agg([
                    # Count distinct secondary representatives
                    pl.col("representative_genome").n_unique().alias("num_secondary_genomes"),
                    # Count representatives in same genus as best hit (ANI >= 80%)
                    pl.col("ani").filter(pl.col("ani") >= genus_threshold).count().alias("same_genus_count"),
                    # Count representatives in same species as best hit (ANI >= 95%)
                    pl.col("ani").filter(pl.col("ani") >= species_threshold).count().alias("same_species_count"),
                    # Min ANI to secondary representatives (used for genus uncertainty)
                    pl.col("ani").min().fill_null(0.0).alias("min_ani_secondary"),
                ])
            )
        else:
            max_ani_per_read = pl.DataFrame({
                "qseqid": [], "max_ani": [], "second_genome_ani": []
            }).cast({
                "qseqid": pl.Utf8, "max_ani": pl.Float64, "second_genome_ani": pl.Float64
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
                pl.col("second_genome_ani").fill_null(0.0),
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
            # Placement uncertainty (genome-level) - mode depends on config
            # "max" mode: uses maximum ANI to any competing genome
            # "second" mode: uses ANI to the second-best genome only
            pl.when(pl.col("num_ambiguous_hits") <= 1)
            .then(0.0)
            .when(pl.lit(cfg.uncertainty_mode == "second"))
            .then(
                pl.when(pl.col("second_genome_ani") == 0.0)
                .then(100.0)
                .otherwise(100.0 - pl.col("second_genome_ani"))
            )
            .otherwise(  # "max" mode (default)
                pl.when(pl.col("max_ani") == 0.0)
                .then(100.0)
                .otherwise(100.0 - pl.col("max_ani"))
            )
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

        # Step 7b: Bayesian-primary classification
        # Uses BayesianClassifier to compute 6-category posteriors and MAP
        # classification, then applies Stage 2 refinement for sub-categories.
        # The legacy threshold cascade is preserved as legacy_call for comparison.
        result = self._bayesian.classify_primary(result)

        # Compute legacy_call using the hard threshold cascade for side-by-side comparison
        from metadarkmatter.core.classification.thresholds import apply_legacy_thresholds
        result = apply_legacy_thresholds(result, self.config)

        # Add diversity_status, is_novel, and low_confidence flags
        bay_cfg = self.config.bayesian
        result = result.with_columns([
            pl.col("taxonomic_call").replace(TAXONOMIC_TO_DIVERSITY).alias("diversity_status"),
            pl.col("taxonomic_call").is_in(["Novel Species", "Novel Genus"]).alias("is_novel"),
            (pl.col("posterior_entropy") > bay_cfg.entropy_threshold).alias("low_confidence"),
        ])

        # Inferred uncertainty for single-hit reads
        result = result.with_columns([
            pl.when(pl.col("num_ambiguous_hits") <= 1)
            .then(
                pl.when(pl.col("novelty_index") < eff["novelty_known_max"])
                .then(5.0 + pl.col("novelty_index") * 0.5)
                .when(pl.col("novelty_index") < eff["novelty_novel_species_max"])
                .then(self._inferred_known_break + (pl.col("novelty_index") - eff["novelty_known_max"]) * 1.0)
                .when(pl.col("novelty_index") < eff["novelty_novel_genus_max"])
                .then(self._inferred_species_break + (pl.col("novelty_index") - eff["novelty_novel_species_max"]) * 1.5)
                .otherwise(35.0)
            )
            .otherwise(pl.lit(None))
            .alias("inferred_uncertainty"),
            pl.when(pl.col("num_ambiguous_hits") <= 1)
            .then(pl.lit("inferred"))
            .otherwise(pl.lit("measured"))
            .alias("uncertainty_type"),
        ])

        # Select and rename final columns
        base_columns = [
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
            "legacy_call",
            "diversity_status",
            "is_novel",
            "low_confidence",
        ]

        # Add posterior columns
        base_columns.extend([
            "p_known_species",
            "p_novel_species",
            "p_novel_genus",
            "p_species_boundary",
            "p_ambiguous",
            "p_unclassified",
            "posterior_entropy",
        ])

        # Add inferred uncertainty columns
        base_columns.extend(["inferred_uncertainty", "uncertainty_type"])

        # Optionally include legacy sub-scores
        if self.config.include_legacy_scores:
            result = self._compute_legacy_scores(result, eff)
            base_columns.extend([
                "alignment_quality",
                "identity_confidence",
                "placement_confidence",
                "discovery_score",
            ])

        classification_result = result.select(base_columns)

        # Add family validation columns if active
        if family_validation_active and family_metrics_df is not None:
            family_cols = family_metrics_df.select([
                pl.col("qseqid"),
                "family_bitscore_ratio",
                "family_identity_gap",
                "in_family_hit_fraction",
                "external_best_genome",
                "external_best_identity",
                pl.col("_best_if_ident").alias("in_family_identity"),
            ])
            classification_result = classification_result.join(
                family_cols,
                left_on="read_id",
                right_on="qseqid",
                how="left",
            )
            classification_result = classification_result.with_columns(
                (100.0 - pl.col("in_family_identity")).alias("in_family_novelty_index"),
            )

            if off_target_df is not None and not off_target_df.is_empty():
                classification_result = pl.concat(
                    [classification_result, off_target_df],
                    how="diagonal_relaxed",
                )

        if compute_qc:
            from metadarkmatter.core.classification.qc import (
                compute_post_qc,
                compute_pre_qc,
            )
            qc = compute_pre_qc(raw_df, filtered_df, self.ani_matrix, self.config)
            qc = compute_post_qc(qc, classification_result)
            return classification_result, qc

        return classification_result

    def _compute_legacy_scores(
        self,
        result: pl.DataFrame,
        eff: dict,
    ) -> pl.DataFrame:
        """Compute legacy sub-scores (opt-in via include_legacy_scores)."""
        result = result.with_columns([
            (pl.col("alignment_fraction") * 50.0 + 25.0).clip(0.0, 100.0).round(1)
            .alias("alignment_quality"),
        ]).with_columns([
            (
                pl.when(pl.col("novelty_index") < 5).then(80.0)
                .when(pl.col("novelty_index") < 15).then(60.0 - (pl.col("novelty_index") - 5))
                .when(pl.col("novelty_index") < 25).then(50.0 - (pl.col("novelty_index") - 15) * 1.5)
                .otherwise(30.0)
                + pl.col("alignment_quality") * 0.2
            ).clip(0.0, 100.0).round(1).alias("identity_confidence"),
        ]).with_columns([
            pl.when(pl.col("num_ambiguous_hits") <= 1)
            .then(
                pl.when(pl.col("inferred_uncertainty") < 2).then(65.0)
                .when(pl.col("inferred_uncertainty") < 5).then(45.0)
                .when(pl.col("inferred_uncertainty") < 10).then(25.0)
                .when(pl.col("inferred_uncertainty") < 20).then(10.0)
                .otherwise(0.0)
            )
            .otherwise(
                pl.when(pl.col("placement_uncertainty") < 2).then(80.0)
                .when(pl.col("placement_uncertainty") < 5).then(60.0)
                .when(pl.col("placement_uncertainty") < 10).then(40.0)
                .when(pl.col("placement_uncertainty") < 20).then(25.0)
                .otherwise(10.0)
                + pl.when((pl.col("identity_gap").is_not_null()) & (pl.col("identity_gap") >= 5)).then(20.0)
                .when((pl.col("identity_gap").is_not_null()) & (pl.col("identity_gap") >= 2)).then(10.0)
                .when((pl.col("identity_gap").is_not_null()) & (pl.col("identity_gap") >= 1)).then(5.0)
                .otherwise(0.0)
            )
            .clip(0.0, 100.0).round(1).alias("placement_confidence"),
        ]).with_columns([
            pl.when(pl.col("taxonomic_call") == "Novel Species")
            .then(
                (15.0 + (pl.col("novelty_index") - eff["novelty_novel_species_min"]) * 1.5).clip(15.0, 30.0)
                + pl.col("alignment_quality") * 0.3
                + (pl.col("identity_confidence") + pl.col("placement_confidence")) / 2.0 * 0.3
            )
            .when(pl.col("taxonomic_call") == "Novel Genus")
            .then(
                (30.0 + (pl.col("novelty_index") - eff["novelty_novel_genus_min"]) * 2.0).clip(30.0, 40.0)
                + pl.col("alignment_quality") * 0.3
                + (pl.col("identity_confidence") + pl.col("placement_confidence")) / 2.0 * 0.3
            )
            .otherwise(pl.lit(None))
            .round(1).alias("discovery_score"),
        ])
        return result

    def _empty_dataframe(self) -> pl.DataFrame:
        """Return empty DataFrame with correct schema."""
        schema = {
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
            "legacy_call": pl.Utf8,
            "diversity_status": pl.Utf8,
            "is_novel": pl.Boolean,
            "low_confidence": pl.Boolean,
            "p_known_species": pl.Float64,
            "p_novel_species": pl.Float64,
            "p_novel_genus": pl.Float64,
            "p_species_boundary": pl.Float64,
            "p_ambiguous": pl.Float64,
            "p_unclassified": pl.Float64,
            "posterior_entropy": pl.Float64,
            "inferred_uncertainty": pl.Float64,
            "uncertainty_type": pl.Utf8,
        }
        if self.config.include_legacy_scores:
            schema.update({
                "alignment_quality": pl.Float64,
                "identity_confidence": pl.Float64,
                "placement_confidence": pl.Float64,
                "discovery_score": pl.Float64,
            })
        return pl.DataFrame(schema=schema)

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
                        result_df = self.classify_file(complete_df)
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
                result_df = self.classify_file(remaining_df)
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


