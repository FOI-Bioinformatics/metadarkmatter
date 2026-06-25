"""Execution runner for ANI-weighted placement classification.

This module holds the load -> classify -> finalize pipeline that the
``score classify`` CLI command (and the ``mdm run`` orchestrator) drive. It is
deliberately free of Rich/Typer dependencies: progress and status text are
delivered through a ``log`` callback and an optional streaming
``progress_callback`` supplied by the caller, so the same logic is reusable
in-process and is unit-testable without a CLI harness.

The caller is expected to have already resolved the :class:`ScoringConfig`
(see ``cli/score.build_classify_config``); ``run_classification`` performs no
flag/preset precedence resolution. It may, however, refine the config in two
documented ways recorded on the returned result: adaptive-threshold detection
and target-family inference from metadata.
"""

from __future__ import annotations

import logging
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import polars as pl

from metadarkmatter.core.aai_matrix_builder import AAIMatrix
from metadarkmatter.core.ani_placement import ANIMatrix, VectorizedClassifier
from metadarkmatter.core.exceptions import ConfigurationError
from metadarkmatter.core.id_mapping import ContigIdMapping
from metadarkmatter.core.io_utils import OutputFormat, write_dataframe
from metadarkmatter.core.metadata import GenomeMetadata
from metadarkmatter.core.parsers import StreamingBlastParser, extract_genome_name_expr
from metadarkmatter.models.config import ScoringConfig

logger = logging.getLogger(__name__)


def validate_ani_genome_coverage(
    blast_path: Path,
    ani_matrix: ANIMatrix,
    sample_rows: int = 10000,
    representative_mapping: dict[str, str] | None = None,
) -> tuple[int, int, float, set[str]]:
    """
    Early validation of ANI matrix coverage against BLAST file genomes.

    Samples the BLAST file and checks how many unique genomes are present
    in the ANI matrix. When representative_mapping is provided, maps BLAST
    genomes to their representatives before checking coverage.

    Args:
        blast_path: Path to BLAST file
        ani_matrix: Loaded ANI matrix
        sample_rows: Number of rows to sample (default 10000)
        representative_mapping: Optional mapping from genome accession to
            representative accession. When provided, checks representative
            coverage rather than raw genome coverage.

    Returns:
        Tuple of (matched_count, total_genomes, coverage_pct, missing_genomes)
    """
    # Read sample of BLAST file (auto-detect column count for backward compatibility)
    parser = StreamingBlastParser(blast_path)
    sample_df = pl.scan_csv(
        blast_path,
        separator="\t",
        has_header=False,
        new_columns=parser.column_names,
        n_rows=sample_rows,
    ).with_columns(
        extract_genome_name_expr()
    ).collect()

    # Get unique genome names from BLAST sample
    blast_genomes = set(sample_df["genome_name"].unique().to_list())
    blast_genomes.discard("unknown")

    # Get ANI matrix genomes
    ani_genomes = set(ani_matrix.genomes)

    # When representative mapping is active, only check family genomes
    # (those present in the mapping). Non-family genomes from broad databases
    # are excluded to avoid inflating the denominator.
    if representative_mapping:
        family_blast_genomes = {g for g in blast_genomes if g in representative_mapping}
        check_genomes = {representative_mapping[g] for g in family_blast_genomes}
    else:
        check_genomes = blast_genomes

    # Calculate coverage
    matched = check_genomes & ani_genomes
    missing = check_genomes - ani_genomes

    total_genomes = len(check_genomes)
    matched_count = len(matched)
    coverage_pct = (100.0 * matched_count / total_genomes) if total_genomes > 0 else 0.0

    return matched_count, total_genomes, coverage_pct, missing


def finalize_classification(
    classification_df: pl.DataFrame | None,
    genome_metadata: Any,
    output: Path,
    output_format: OutputFormat,
) -> int:
    """
    Finalize classification: join metadata and write output.

    Args:
        classification_df: Classification results DataFrame.
        genome_metadata: Optional GenomeMetadata for species/genus joining.
        output: Output file path.
        output_format: Output format ('csv' or 'parquet').

    Returns:
        Number of classified reads.
    """
    if classification_df is None:
        return 0

    # Join species/genus metadata if provided
    if genome_metadata is not None:
        classification_df = genome_metadata.join_classifications(classification_df)

    num_classified = len(classification_df)
    if num_classified > 0:
        write_dataframe(classification_df, output, output_format)

    return num_classified


@dataclass(frozen=True)
class ClassificationRequest:
    """Inputs for a single classification run.

    The ``config`` is assumed fully resolved (precedence already applied). All
    paths must exist; existence/permission checks are the caller's job so that
    typed errors can be surfaced consistently at the CLI boundary.
    """

    alignment: Path
    ani: Path
    output: Path
    config: ScoringConfig
    output_format: OutputFormat = "csv"
    aai: Path | None = None
    metadata: Path | None = None
    genomes: Path | None = None
    id_mapping: Path | None = None
    streaming: bool = False
    chunk_size: int = 5_000_000
    strict_ani: bool = False
    adaptive_thresholds: bool = False
    compute_qc: bool = False


@dataclass(frozen=True)
class ClassificationRunResult:
    """Outcome of :func:`run_classification`.

    ``config`` is the possibly-refined config actually used (adaptive
    thresholds and/or inferred target family applied). ``classification_df`` is
    ``None`` in streaming mode, where rows are written incrementally to disk.
    """

    num_classified: int
    config: ScoringConfig
    classification_df: pl.DataFrame | None = None
    qc_metrics: Any | None = None
    coverage: tuple[int, int, float, set[str]] | None = field(default=None)


def _noop_log(_message: str) -> None:
    """Default log sink: discard."""


def run_classification(
    request: ClassificationRequest,
    *,
    log: Callable[[str], None] = _noop_log,
    verbose: bool = False,
    progress_callback: Callable[[int, int, float], None] | None = None,
) -> ClassificationRunResult:
    """Execute the load -> classify -> finalize pipeline for one alignment file.

    Args:
        request: Fully-resolved inputs (see :class:`ClassificationRequest`).
        log: Receives human-facing status lines (may contain Rich markup; the
            CLI renders it, other callers may strip or ignore it).
        verbose: When true, emits additional coverage-detail lines via ``log``.
        progress_callback: Streaming-mode progress hook ``(rows, reads,
            elapsed_seconds)``; ignored when ``streaming`` is false. The caller
            owns any live progress display and keeps it alive across this call.

    Returns:
        A :class:`ClassificationRunResult` with the count, the config used, and
        (non-streaming only) the in-memory DataFrame and optional QC metrics.

    Raises:
        ConfigurationError: For mutually-exclusive inputs, a missing adaptive
            dependency, or ID-mapping construction failures. Other I/O and
            Polars errors propagate unchanged to the centralized CLI handler.
    """
    config = request.config

    # Summarize the active alignment filters so the user can see what is applied.
    active_filters = []
    if config.min_alignment_length > 0:
        active_filters.append(f"length >= {config.min_alignment_length}bp")
    if config.min_alignment_fraction > 0:
        active_filters.append(f"fraction >= {config.min_alignment_fraction:.0%}")
    if config.max_evalue > 0:
        active_filters.append(f"evalue <= {config.max_evalue:.0e}")
    if config.min_percent_identity > 0:
        active_filters.append(f"identity >= {config.min_percent_identity:.1f}%")
    if config.min_bitscore > 0:
        active_filters.append(f"bitscore >= {config.min_bitscore:.1f}")
    if config.min_read_length > 0:
        active_filters.append(f"read length >= {config.min_read_length}bp")
    if config.min_query_coverage > 0:
        active_filters.append(f"query coverage >= {config.min_query_coverage:.1f}%")
    if active_filters:
        log(f"[dim]Alignment filters: {', '.join(active_filters)}[/dim]")

    # Load ANI matrix. Load errors (incl. typed ANIMatrixError with its
    # suggestion) propagate to the centralized CLI error handler.
    ani_matrix = ANIMatrix.from_file(
        request.ani, symmetry_check="strict" if request.strict_ani else "warn"
    )
    log(f"[green]Loaded ANI matrix for {len(ani_matrix.genomes)} genomes[/green]\n")

    # Adaptive threshold detection (override config if enabled)
    if request.adaptive_thresholds:
        try:
            from metadarkmatter.core.classification.adaptive import (
                build_adaptive_config,
                detect_species_boundary,
            )
        except ImportError:
            raise ConfigurationError(
                "scikit-learn required for adaptive thresholds.",
                suggestion="Install with: pip install metadarkmatter[adaptive]",
            ) from None

        adaptive = detect_species_boundary(ani_matrix)
        if adaptive.method == "gmm":
            config = build_adaptive_config(config, adaptive)
            log(
                f"[green]Adaptive threshold:[/green] species boundary at "
                f"{adaptive.species_boundary:.1f}% ANI "
                f"(novelty_known_max={adaptive.novelty_known_max:.1f}%, "
                f"confidence={adaptive.confidence:.2f})\n"
            )
        else:
            log(
                f"[yellow]Adaptive threshold detection fell back to default "
                f"({adaptive.species_boundary:.1f}% ANI)[/yellow]\n"
            )

    # Load AAI matrix if provided (for genus-level classification)
    aai_matrix: AAIMatrix | None = None
    if request.aai:
        aai_matrix = AAIMatrix.from_file(request.aai)
        log(f"[green]Loaded AAI matrix for {len(aai_matrix.genomes)} genomes[/green]")
        log("[dim]AAI will be used for genus-level classification (20-25% novelty)[/dim]\n")

    # Load genome metadata if provided
    genome_metadata: GenomeMetadata | None = None
    representative_mapping: dict[str, str] | None = None
    if request.metadata:
        genome_metadata = GenomeMetadata.from_file(request.metadata)
        log(
            f"[green]Loaded metadata for {genome_metadata.genome_count} genomes "
            f"({genome_metadata.species_count} species)[/green]"
        )
        if genome_metadata.has_representatives:
            representative_mapping = genome_metadata.build_representative_mapping()
            log(
                f"[dim]Representative mapping: {genome_metadata.genome_count} genomes -> "
                f"{genome_metadata.representative_count} representatives[/dim]"
            )
        log("")

    # Infer target family from metadata if not explicitly provided
    if config.target_family is None and genome_metadata is not None:
        inferred = genome_metadata.infer_target_family()
        if inferred:
            log(f"[dim]Inferred target family from metadata: {inferred}[/dim]")
            config = config.model_copy(update={"target_family": inferred})

    if config.target_family:
        log(
            f"[cyan]Family validation: target={config.target_family}, "
            f"threshold={config.family_ratio_threshold}[/cyan]\n"
        )

    # Load or generate ID mapping for external BLAST results
    contig_mapping: ContigIdMapping | None = None
    if request.genomes and request.id_mapping:
        raise ConfigurationError(
            "Cannot specify both --genomes and --id-mapping. Use one or the other."
        )

    if request.genomes:
        try:
            contig_mapping = ContigIdMapping.from_genome_dir(request.genomes)
        except (FileNotFoundError, ValueError, OSError) as e:
            raise ConfigurationError(f"Error generating ID mapping: {e}") from e
        log(f"[green]Generated ID mapping for {len(contig_mapping):,} contigs[/green]\n")
    elif request.id_mapping:
        try:
            contig_mapping = ContigIdMapping.from_tsv(request.id_mapping)
        except (FileNotFoundError, ValueError, pl.exceptions.PolarsError, OSError) as e:
            raise ConfigurationError(f"Error loading ID mapping: {e}") from e
        log(f"[green]Loaded ID mapping with {len(contig_mapping):,} entries[/green]\n")

    # Early validation: check genome coverage between alignment file and ANI matrix.
    # When representative mapping is active, validate representative coverage.
    coverage: tuple[int, int, float, set[str]] | None = None
    try:
        matched, total, coverage_pct, missing = validate_ani_genome_coverage(
            request.alignment, ani_matrix, representative_mapping=representative_mapping,
        )
        coverage = (matched, total, coverage_pct, missing)

        if verbose:
            if representative_mapping:
                log(
                    f"[dim]Representative coverage: {matched}/{total} representatives in ANI matrix "
                    f"({coverage_pct:.1f}%)[/dim]"
                )
            else:
                log(
                    f"[dim]Genome coverage: {matched}/{total} genomes in ANI matrix "
                    f"({coverage_pct:.1f}%)[/dim]"
                )

        if coverage_pct < 50.0 and total > 0:
            log(
                f"[yellow]Warning: Low genome coverage ({coverage_pct:.1f}%)[/yellow]\n"
                f"  Only {matched} of {total} genomes in BLAST file are in ANI matrix.\n"
                f"  This may indicate mismatched input files.\n"
            )
            if missing and len(missing) <= 5:
                log(f"  Missing genomes: {', '.join(sorted(missing))}\n")
            elif missing:
                sample = sorted(missing)[:5]
                log(f"  Missing genomes (first 5): {', '.join(sample)}...\n")
    except (FileNotFoundError, ValueError, pl.exceptions.PolarsError, OSError) as e:
        # Always log warnings - don't silently swallow validation failures
        logger.warning("Could not validate genome coverage: %s", e)
        if verbose:
            log(f"[dim]Warning: Could not validate genome coverage: {e}[/dim]")

    # Run classification - keep DataFrame in memory to avoid redundant I/O.
    # Classification errors (Polars, file-system, memory) propagate to the
    # centralized CLI error handler, which shows a full traceback only under
    # --debug.
    classification_df: pl.DataFrame | None = None
    qc_metrics = None

    if request.streaming:
        # Streaming mode writes directly to file - no in-memory DataFrame.
        log("[cyan]Streaming mode: processing in 5M alignment partitions[/cyan]\n")

        if genome_metadata:
            log(
                "[yellow]Note: --metadata is not supported with --streaming mode.[/yellow]\n"
                "[yellow]Species columns will not be added to output.[/yellow]\n"
            )

        if contig_mapping:
            log(
                "[yellow]Note: --genomes/--id-mapping is not supported with --streaming mode.[/yellow]\n"
                "[yellow]ID transformation will not be applied. Remove --streaming to enable it.[/yellow]\n"
            )

        vectorized = VectorizedClassifier(
            ani_matrix=ani_matrix, aai_matrix=aai_matrix, config=config,
            representative_mapping=representative_mapping,
        )
        num_classified = vectorized.stream_to_file(
            blast_path=request.alignment,
            output_path=request.output,
            output_format=request.output_format,
            partition_size=request.chunk_size,
            progress_callback=progress_callback,
        )
    else:
        # All non-streaming classification uses VectorizedClassifier.
        vectorized = VectorizedClassifier(
            ani_matrix=ani_matrix, aai_matrix=aai_matrix, config=config,
            representative_mapping=representative_mapping,
        )
        result = vectorized.classify_file(
            request.alignment,
            id_mapping=contig_mapping,
            compute_qc=request.compute_qc,
        )
        if request.compute_qc:
            assert isinstance(result, tuple)
            classification_df, qc_metrics = result
        else:
            assert isinstance(result, pl.DataFrame)
            classification_df = result

        num_classified = finalize_classification(
            classification_df, genome_metadata, request.output, request.output_format
        )

    log(f"[green]Classified {num_classified:,} reads[/green]\n")

    return ClassificationRunResult(
        num_classified=num_classified,
        config=config,
        classification_df=classification_df,
        qc_metrics=qc_metrics,
        coverage=coverage,
    )
