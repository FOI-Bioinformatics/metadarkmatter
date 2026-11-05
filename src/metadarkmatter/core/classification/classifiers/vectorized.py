"""
Fully vectorized classifier using Polars for maximum performance.

This module provides VectorizedClassifier, which performs all classification
operations in Polars' native Rust backend for automatic parallelization
without multiprocessing overhead.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar

import numpy as np
import polars as pl

from metadarkmatter.core.classification.bayesian import (
    BayesianClassifier,
)
from metadarkmatter.core.constants import (
    UNKNOWN_GENOME,
)
from metadarkmatter.core.io_utils import OutputFormat, StreamingWriter
from metadarkmatter.core.parsers import extract_genome_name_expr
from metadarkmatter.models.classification import TAXONOMIC_TO_DIVERSITY
from metadarkmatter.models.config import EffectiveThresholds, ScoringConfig

if TYPE_CHECKING:
    from metadarkmatter.core.aai_matrix_builder import AAIMatrix
    from metadarkmatter.core.classification.ani_matrix import ANIMatrix
    from metadarkmatter.core.classification.qc import QCMetrics
    from metadarkmatter.core.id_mapping import ContigIdMapping

logger = logging.getLogger(__name__)


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
        vectorized = VectorizedClassifier(ani_matrix, config=config)
        df = vectorized.classify_file(blast_path)
    """

    # Lookup tables built in __init__. The ANI lookup is always a DataFrame;
    # the AAI lookups are None when no AAI matrix is provided.
    _ani_lookup: pl.DataFrame
    _ani_lookup_symmetric: pl.DataFrame
    _aai_lookup: pl.DataFrame | None
    _aai_lookup_symmetric: pl.DataFrame | None

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
        if (
            self.config.clade_novel_species
            and "clade_novel_species" in self.config.model_fields_set
            and self.config.alignment_mode != "nucleotide"
        ):
            logger.warning(
                "clade_novel_species is only validated for nucleotide mode and "
                "is not applied with alignment_mode=%r.",
                self.config.alignment_mode,
            )
        if representative_mapping:
            # Null/blank representatives fall back to the genome itself so
            # competitors do not vanish through a null representative.
            representative_mapping = {
                k: (v if isinstance(v, str) and v.strip() else k)
                for k, v in representative_mapping.items()
            }
        self.representative_mapping = representative_mapping
        # Diagnostics: hit genomes (or representatives) absent from the ANI
        # matrix and reads whose competitor ANI lookup failed (U defaults to 100).
        self.missing_ani_genomes: set[str] = set()
        self.missing_ani_reads: set[str] = set()

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
        # Float64 so uncertainty columns have a stable dtype (a float32
        # lookup made placement_uncertainty Float32 only in batches with
        # competing hits, which broke incremental Parquet writes).
        ani_arr = self.ani_matrix._ani_array[idx1, idx2].astype(np.float64)

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

    def _record_missing_ani(self, secondary_with_ani: pl.DataFrame) -> None:
        """Accumulate stats on hit genomes/reads with ANI pairs missing."""
        known = set(self.ani_matrix.genomes)
        hit_reps = set(secondary_with_ani["representative_genome"].unique().to_list())
        hit_reps |= set(secondary_with_ani["best_representative"].unique().to_list())
        self.missing_ani_genomes |= {g for g in hit_reps if g not in known}
        unresolved = secondary_with_ani.filter(pl.col("ani").is_null())
        self.missing_ani_reads |= set(unresolved["qseqid"].unique().to_list())

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

    def _check_memory(self, blast_path: Path) -> None:
        """Estimate memory usage and warn if it may exceed available memory."""
        try:
            file_size = blast_path.stat().st_size
        except OSError:
            return

        # Estimate in-memory size: gzipped files expand ~10x, then Polars adds ~1.5x overhead
        is_gzipped = blast_path.name.endswith((".gz", ".gzip"))
        estimated_bytes = file_size * 10 * 1.5 if is_gzipped else file_size * 1.5

        logger.debug(
            "Estimated memory usage for %s: %.1f GB",
            blast_path.name,
            estimated_bytes / (1024**3),
        )

        # Get total physical memory (stdlib only, Linux/macOS)
        try:
            page_size = os.sysconf("SC_PAGE_SIZE")
            page_count = os.sysconf("SC_PHYS_PAGES")
            total_memory = page_size * page_count
        except (ValueError, OSError, AttributeError):
            return

        available_memory = total_memory * 0.6  # conservative estimate

        if estimated_bytes > available_memory * 0.8:
            est_gb = estimated_bytes / (1024**3)
            avail_gb = available_memory / (1024**3)
            logger.warning(
                "Estimated memory usage (%.1f GB) may exceed available memory "
                "(%.1f GB). Consider --streaming (score classify) for bounded memory.",
                est_gb,
                avail_gb,
            )

    def _apply_alignment_filters(self, df: pl.DataFrame) -> pl.DataFrame:
        """Apply GTDB-compatible alignment-quality filters to raw BLAST rows.

        Shared by the in-memory (file) and streaming (DataFrame partition)
        paths so both classify an identical set of reads. Operates on raw
        BLAST columns (length, qstart/qend, evalue, pident, bitscore, qlen);
        does not require the genome_name column.
        """
        cfg = self.config

        # Normalize query coordinates: MMseqs2 reports minus-strand hits
        # with qstart > qend. Every downstream coverage expression uses
        # (qend - qstart + 1), which would otherwise be negative and
        # down-weight (or filter out) the true best hit of reverse reads.
        if "qstart" in df.columns and "qend" in df.columns:
            df = df.with_columns(
                pl.min_horizontal("qstart", "qend").alias("qstart"),
                pl.max_horizontal("qstart", "qend").alias("qend"),
            )

        n_before = df.height
        filter_counts: list[tuple[str, int]] = []

        def _record(name: str, frame: pl.DataFrame, before: int) -> None:
            if frame.height < before:
                filter_counts.append((name, before - frame.height))

        if cfg.min_alignment_length > 0:
            # BLASTX/DIAMOND report `length` in amino acids; in protein mode
            # compare against the aligned query span in nucleotides instead.
            if cfg.alignment_mode == "protein":
                len_expr = (pl.col("qend") - pl.col("qstart") + 1)
            else:
                len_expr = pl.col("length")
            h = df.height
            df = df.filter(len_expr >= cfg.min_alignment_length)
            _record(f"min_alignment_length={cfg.min_alignment_length}", df, h)

        # Alignment fraction = aligned query span / read length. Use qlen
        # when the 13/14-column format provides it; only fall back to the
        # qend proxy for legacy 12-column input (where the result depends on
        # where the alignment sits on the read).
        if cfg.min_alignment_fraction > 0:
            read_len = (
                pl.coalesce(pl.col("qlen"), pl.col("qend"))
                if "qlen" in df.columns
                else pl.col("qend")
            )
            h = df.height
            df = df.with_columns([
                ((pl.col("qend") - pl.col("qstart") + 1) / read_len).alias("_approx_af")
            ]).filter(
                pl.col("_approx_af") >= cfg.min_alignment_fraction
            ).drop("_approx_af")
            _record(f"min_alignment_fraction={cfg.min_alignment_fraction}", df, h)

        # E-value filter
        if cfg.max_evalue > 0:
            h = df.height
            df = df.filter(pl.col("evalue") <= cfg.max_evalue)
            _record(f"max_evalue={cfg.max_evalue}", df, h)

        # Percent identity filter
        if cfg.min_percent_identity > 0:
            h = df.height
            df = df.filter(pl.col("pident") >= cfg.min_percent_identity)
            _record(f"min_percent_identity={cfg.min_percent_identity}", df, h)

        # Bitscore filter
        if cfg.min_bitscore > 0:
            h = df.height
            df = df.filter(pl.col("bitscore") >= cfg.min_bitscore)
            _record(f"min_bitscore={cfg.min_bitscore}", df, h)

        # Read length and query coverage filters
        if cfg.min_read_length > 0 or cfg.min_query_coverage > 0:
            if "qlen" in df.columns:
                df = df.with_columns(
                    pl.when(pl.col("qlen").is_not_null())
                    .then(pl.col("qlen"))
                    .otherwise(pl.col("qend"))
                    .alias("_read_length")
                )
            else:
                df = df.with_columns(pl.col("qend").alias("_read_length"))

            if cfg.min_read_length > 0:
                h = df.height
                df = df.filter(pl.col("_read_length") >= cfg.min_read_length)
                _record(f"min_read_length={cfg.min_read_length}", df, h)

            if cfg.min_query_coverage > 0:
                h = df.height
                df = df.with_columns(
                    ((pl.col("qend") - pl.col("qstart") + 1)
                     / pl.col("_read_length") * 100.0)
                    .clip(0.0, 100.0)
                    .alias("_query_cov")
                ).filter(
                    pl.col("_query_cov") >= cfg.min_query_coverage
                ).drop("_query_cov")
                _record(f"min_query_coverage={cfg.min_query_coverage}", df, h)

            df = df.drop("_read_length")

        if n_before > 0 and df.height == 0:
            logger.warning(
                "All %d alignment rows were removed by the alignment filters "
                "(rows removed per filter: %s). Relax or disable the responsible "
                "filter(s); note that in protein mode min_alignment_length is "
                "compared against the query span in nucleotides.",
                n_before,
                ", ".join(f"{k}: {v}" for k, v in filter_counts) or "none",
            )

        return df

    def classify_file(
        self,
        blast_input: Path | pl.DataFrame,
        id_mapping: ContigIdMapping | None = None,
        compute_qc: bool = False,
    ) -> pl.DataFrame | tuple[pl.DataFrame, QCMetrics]:
        """Classify BLAST hits; see ``_classify_file_impl``.

        The result is conformed to a fixed output schema (column set, order
        and dtypes) so every code path -- including the early returns for
        empty or all-Off-target inputs -- yields identical columns. This is
        what keeps streaming batches aligned under one CSV header / Parquet
        schema.
        """
        out = self._classify_file_impl(blast_input, id_mapping, compute_qc)
        if isinstance(out, tuple):
            return self._conform_output(out[0]), out[1]
        return self._conform_output(out)

    _FAMILY_OUTPUT_SCHEMA: ClassVar[dict[str, Any]] = {
        "family_bitscore_ratio": pl.Float64,
        "family_identity_gap": pl.Float64,
        "in_family_hit_fraction": pl.Float64,
        "external_best_genome": pl.Utf8,
        "external_best_identity": pl.Float64,
        "in_family_identity": pl.Float64,
        "in_family_novelty_index": pl.Float64,
    }

    @property
    def _clade_rule_active(self) -> bool:
        """Clade-novelty rule enabled and applicable (nucleotide mode only)."""
        return (
            self.config.clade_novel_species
            and self.config.alignment_mode == "nucleotide"
        )

    def _output_schema(self) -> dict[str, Any]:
        schema = dict(self._empty_dataframe().schema)
        if self.config.target_family is not None:
            schema.update(self._FAMILY_OUTPUT_SCHEMA)
        if self._clade_rule_active:
            schema["clade_novel_species"] = pl.Boolean()
        return schema

    def _conform_output(self, df: pl.DataFrame) -> pl.DataFrame:
        schema = self._output_schema()
        missing = [
            pl.lit(None).cast(dtype).alias(name)
            for name, dtype in schema.items()
            if name not in df.columns
        ]
        if missing:
            df = df.with_columns(missing)
        return df.select(
            [pl.col(name).cast(dtype) for name, dtype in schema.items()]
        )

    def _classify_file_impl(
        self,
        blast_input: Path | pl.DataFrame,
        id_mapping: ContigIdMapping | None = None,
        compute_qc: bool = False,
    ) -> pl.DataFrame | tuple[pl.DataFrame, QCMetrics]:
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
            # DataFrame passed directly (e.g. from stream_to_file partitions),
            # with the genome_name column already added. Apply the same
            # alignment-quality filters as the file path so streaming and
            # in-memory modes classify an identical set of reads.
            df = blast_input
            raw_df = df if compute_qc else None
            df = self._apply_alignment_filters(df)
        else:
            self._check_memory(blast_input)

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
            df = self._apply_alignment_filters(df)

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
                    .sort(
                        ["qseqid", "bitscore", "pident", "sseqid"],
                        descending=[False, True, True, False],
                        maintain_order=True,
                    )
                    .group_by("qseqid", maintain_order=True)
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

                # For off-target reads we cannot compute a true placement
                # uncertainty (the matching genome is outside the ANI matrix),
                # so leave it null rather than reporting a misleading 0.0.
                # Confidence is intentionally low (10.0) rather than zero so
                # downstream filters can still rank off-target hits by their
                # external identity if needed.
                off_target_cols = [
                    pl.col("qseqid").alias("read_id"),
                    pl.col("external_best_genome").fill_null("unknown").alias("best_match_genome"),
                    pl.col("external_best_identity").alias("top_hit_identity"),
                    (100.0 - pl.col("external_best_identity")).alias("novelty_index"),
                    pl.lit(None).cast(pl.Float64).alias("placement_uncertainty"),
                    pl.lit(None).cast(pl.Float64).alias("genus_uncertainty"),
                    pl.lit("off_target").alias("ambiguity_scope"),
                    pl.lit(0).cast(pl.Int64).alias("num_ambiguous_hits"),
                    pl.lit(None).cast(pl.Float64).alias("second_hit_identity"),
                    pl.lit(None).cast(pl.Float64).alias("identity_gap"),
                    pl.lit(10.0).alias("confidence_score"),
                    pl.lit("Off-target").alias("taxonomic_call"),
                    pl.lit("Off-target").alias("legacy_call"),
                    pl.lit("Off-target").alias("diversity_status"),
                    pl.lit(False).alias("is_novel"),
                    pl.lit(True).alias("low_confidence"),
                    pl.lit(None).cast(pl.Float64).alias("inferred_uncertainty"),
                    pl.lit("off_target").alias("uncertainty_type"),
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

                off_target_df = off_target_reads.select(off_target_cols).sort("read_id")

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
                        compute_pre_qc,
                    )
                    # raw_df/filtered_df are set whenever compute_qc is True.
                    assert raw_df is not None and filtered_df is not None
                    qc = compute_pre_qc(raw_df, filtered_df, self.ani_matrix, self.config)
                    return off_target_df, qc
                return off_target_df
            empty = self._empty_dataframe()
            if compute_qc:
                from metadarkmatter.core.classification.qc import (
                    compute_pre_qc,
                )
                assert raw_df is not None and filtered_df is not None
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
                # Second-best weighted_bitscore for gap calculation. Use
                # slice(1, 1).first() rather than .get(1) so single-hit
                # groups yield null instead of raising on polars >= 1.37.
                pl.col("weighted_bitscore").slice(1, 1).first().alias("second_weighted_bitscore"),
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
                # Alignment fraction (qlen when available, else qend proxy)
                ((pl.col("qend") - pl.col("qstart") + 1) / pl.col("read_length"))
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
        if not secondary.is_empty():
            secondary_with_ani = secondary.join(
                self._ani_lookup_symmetric,
                left_on=["best_representative", "representative_genome"],
                right_on=["genome1", "genome2"],
                how="left",
            )

            # Track hit genomes / reads whose ANI pair is absent from the
            # matrix (they silently get U=100); reported by the runner.
            self._record_missing_ani(secondary_with_ani)

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
                # Same total order as best-hit selection (weighted bitscore,
                # then pident, then genome name). Sorting on raw bitscore
                # alone left ties to polars' unstable sort, so the chosen
                # "second" genome (and thus placement_uncertainty) depended
                # on batch size / thread scheduling, and streaming and
                # in-memory runs disagreed.
                .sort(
                    ["qseqid", "weighted_bitscore", "pident", "representative_genome"],
                    descending=[False, True, True, False],
                    maintain_order=True,
                )
                .group_by("qseqid", maintain_order=True)
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
            genus_threshold = self.config.same_genus_ani_threshold
            species_threshold = self.config.same_species_ani_threshold

            phylo_metrics = (
                secondary_with_ani
                .group_by("qseqid")
                .agg([
                    # Count distinct secondary representatives
                    pl.col("representative_genome").n_unique().alias("num_secondary_genomes"),
                    # Count DISTINCT representatives in the same genus / species
                    # as the best hit. These are compared with
                    # num_secondary_genomes (also distinct), so they must not
                    # count hit rows: several contigs of one genome would
                    # otherwise make the counts unequal and mislabel the read
                    # as across_genera (which can turn Ambiguous into
                    # Conserved Region).
                    pl.col("representative_genome")
                    .filter(pl.col("ani") >= genus_threshold)
                    .n_unique()
                    .alias("same_genus_count"),
                    pl.col("representative_genome")
                    .filter(pl.col("ani") >= species_threshold)
                    .n_unique()
                    .alias("same_species_count"),
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

        # Step 5c: AAI-based genus uncertainty (only when an AAI matrix is
        # supplied). Mirrors ANIWeightedClassifier: 100 - max AAI between the
        # best genome and any competing genome among the read's hits, using
        # only pairs with AAI > 0. Before this, --aai was loaded but never
        # used by the CLI classifier.
        aai_metrics: pl.DataFrame | None = None
        if (
            self._aai_lookup_symmetric is not None
            and not self._aai_lookup_symmetric.is_empty()
            and self.config.use_aai_for_genus
        ):
            competitors = df_with_threshold.filter(
                pl.col("representative_genome") != pl.col("best_representative")
            ).select(["qseqid", "best_representative", "representative_genome"])
            if not competitors.is_empty():
                aai_metrics = (
                    competitors.join(
                        self._aai_lookup_symmetric,
                        left_on=["best_representative", "representative_genome"],
                        right_on=["genome1", "genome2"],
                        how="inner",
                    )
                    .filter(pl.col("aai") > 0.0)
                    .group_by("qseqid")
                    .agg((100.0 - pl.col("aai").max()).cast(pl.Float64).alias("_aai_uncertainty"))
                )

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
            .when(pl.col("num_secondary_genomes") == 0)
            .then(0.0)  # All hits collapse to same representative = confident placement
            .when(pl.lit(cfg.uncertainty_mode == "second"))
            .then(
                pl.when(pl.col("second_genome_ani") == 0.0)
                .then(100.0)
                .otherwise(100.0 - pl.col("second_genome_ani"))
            )
            .otherwise(  # "max" mode (default is "second")
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

        # Optional clade-novelty refinement (nucleotide mode only): a read
        # divergent from every competing reference of one genus is a novel
        # species even when its placement among them is ambiguous.
        if self._clade_rule_active:
            from metadarkmatter.core.classification.bayesian import (
                apply_clade_novelty_rule,
            )
            result = apply_clade_novelty_rule(
                result, self.config.clade_novelty_min, self.config.clade_novelty_max
            )

        # Compute legacy_call using the hard threshold cascade for side-by-side comparison
        from metadarkmatter.core.classification.thresholds import apply_legacy_thresholds
        result = apply_legacy_thresholds(result, self.config)

        # AAI genus override (same rule as ANIWeightedClassifier): within the
        # Novel Genus novelty range and below the genus uncertainty cap,
        # AAI > aai_genus_boundary_high -> Novel Species, AAI in the boundary
        # zone -> Ambiguous Within Genus, AAI < aai_genus_boundary_low ->
        # Novel Genus. Applied to both calls before the derived flags.
        if aai_metrics is not None:
            result = result.join(aai_metrics, on="qseqid", how="left")
            in_genus_range = (
                pl.col("_aai_uncertainty").is_not_null()
                & (pl.col("novelty_index") >= eff["novelty_novel_genus_min"])
                & (pl.col("novelty_index") <= eff["novelty_novel_genus_max"])
                & (pl.col("placement_uncertainty") < eff["uncertainty_novel_genus_max"])
            )
            low = 100.0 - self.config.aai_genus_boundary_high
            high = 100.0 - self.config.aai_genus_boundary_low
            aai_call = (
                pl.when(pl.col("_aai_uncertainty") < low)
                .then(pl.lit("Novel Species"))
                .when(pl.col("_aai_uncertainty") <= high)
                .then(pl.lit("Ambiguous Within Genus"))
                .otherwise(pl.lit("Novel Genus"))
            )
            result = result.with_columns([
                pl.when(in_genus_range).then(aai_call).otherwise(pl.col(c)).alias(c)
                for c in ("taxonomic_call", "legacy_call")
            ]).drop("_aai_uncertainty")

        # Add diversity_status, is_novel, and low_confidence flags
        bay_cfg = self.config.bayesian
        result = result.with_columns([
            pl.col("taxonomic_call").replace(TAXONOMIC_TO_DIVERSITY).alias("diversity_status"),
            pl.col("taxonomic_call").is_in(["Novel Species", "Novel Genus"]).alias("is_novel"),
            (pl.col("posterior_entropy") > bay_cfg.entropy_threshold).alias("low_confidence"),
        ])

        # Inferred uncertainty for single-hit reads and representative-collapsed reads
        # (where all ambiguous hits map to the same representative genome)
        _needs_inferred = (
            (pl.col("num_ambiguous_hits") <= 1)
            | (pl.col("num_secondary_genomes") == 0)
        )
        result = result.with_columns([
            pl.when(_needs_inferred)
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
            pl.when(_needs_inferred)
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
        if self._clade_rule_active:
            base_columns.append("clade_novel_species")

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
            assert raw_df is not None and filtered_df is not None
            qc = compute_pre_qc(raw_df, filtered_df, self.ani_matrix, self.config)
            qc = compute_post_qc(qc, classification_result)
            return classification_result, qc

        return classification_result

    def _compute_legacy_scores(
        self,
        result: pl.DataFrame,
        eff: EffectiveThresholds,
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
        output_format: OutputFormat = "parquet",
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
            vectorized = VectorizedClassifier(ani_matrix, config=config)

            def progress(rows, reads, elapsed):
                print(f"Processed {rows:,} rows -> {reads:,} reads in {elapsed:.1f}s")

            total = vectorized.stream_to_file(
                blast_path, output_path,
                progress_callback=progress
            )
        """
        import time


        start_time = time.time()
        total_reads = 0
        total_rows = 0
        # Use scan_csv with collect_batches for streaming
        # Detect column count from file (12 for legacy, 13 for new format with qlen)
        from metadarkmatter.core.parsers import StreamingBlastParser
        parser = StreamingBlastParser(blast_path)

        # Reuse the parser's explicit schema and comment handling; letting
        # polars infer dtypes from the first rows crashed on integer-looking
        # bitscores followed by decimals, numeric read IDs, or '#' headers.
        batches = parser.parse_lazy().collect_batches(chunk_size=partition_size)

        # Track reads across partitions to handle boundary cases
        pending_reads: dict[str, list[Any]] = {}
        # 64-bit hashes of every read id already emitted (~8 bytes/read),
        # used to detect non-contiguous input without holding the ids.
        emitted_hashes = np.empty(0, dtype=np.uint64)

        # Incremental writer: each partition is appended exactly once (parquet
        # row group / CSV append), keeping total write cost O(rows) instead of
        # the O(partitions^2) read-concat-rewrite this previously used.
        with StreamingWriter(output_path, output_format) as writer:
            for partition_df in batches:
                if partition_df.is_empty():
                    continue

                total_rows += len(partition_df)

                # Add genome_name column
                partition_df = partition_df.with_columns([extract_genome_name_expr()])

                # Get unique reads in this partition
                partition_reads = set(partition_df["qseqid"].unique().to_list())

                # Streaming requires each read's hits to be contiguous (as
                # BLAST/MMseqs2 emit them). Detect violations instead of
                # silently classifying a read twice from partial hit sets:
                # (a) a read id appearing in more than one run within this
                # partition, (b) a read already emitted from an earlier
                # partition reappearing here.
                qids = partition_df["qseqid"]
                n_runs = int((qids != qids.shift(1)).fill_null(True).sum())
                uniq = qids.unique(maintain_order=True)
                uniq_hashes = uniq.hash(seed=0).to_numpy()
                seen_mask = np.isin(uniq_hashes, emitted_hashes)
                reappeared = uniq.filter(pl.Series(seen_mask)).to_list()
                if n_runs != len(partition_reads) or reappeared:
                    example = reappeared[0] if reappeared else None
                    msg = (
                        "Streaming mode requires alignment rows grouped by "
                        "read (query) id, but hits for the same read are not "
                        "contiguous"
                        + (f" (e.g. {example!r})" if example else "")
                        + ". Sort the file by the first column "
                        "(e.g. `sort -k1,1 -s`) or run without --streaming."
                    )
                    raise ValueError(msg)

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

                        # Classify this partition (no compute_qc -> DataFrame)
                        if not complete_df.is_empty():
                            result_df = self.classify_file(complete_df)
                            assert isinstance(result_df, pl.DataFrame)
                            total_reads += len(result_df)
                            writer.write(result_df)

                    if complete_reads:
                        done = pl.Series(sorted(complete_reads), dtype=pl.Utf8)
                        emitted_hashes = np.union1d(
                            emitted_hashes, done.hash(seed=0).to_numpy()
                        )

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
                    assert isinstance(result_df, pl.DataFrame)
                    total_reads += len(result_df)
                    writer.write(result_df)

        if progress_callback:
            elapsed = time.time() - start_time
            progress_callback(total_rows, total_reads, elapsed)

        return total_reads


