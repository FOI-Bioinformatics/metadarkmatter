"""
Score command for ANI-weighted placement classification.

This is the core command that most users will interact with, as it implements
the ANI-weighted placement uncertainty algorithm on existing BLAST results.
"""

from __future__ import annotations

import gzip
import json
import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import polars as pl
import typer

logger = logging.getLogger(__name__)
from rich.console import Console
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)
from rich.table import Table

from metadarkmatter.cli.utils import QuietConsole
from metadarkmatter.core.ani_placement import (
    ANIMatrix,
)
from metadarkmatter.core.classification.runner import (
    ClassificationRequest,
    ClassificationRunResult,
    run_classification,
    validate_ani_genome_coverage,
)
from metadarkmatter.core.exceptions import ConfigurationError, MetadarkmatterError
from metadarkmatter.core.io_utils import OutputFormat, read_dataframe
from metadarkmatter.core.provenance import (
    build_provenance,
    provenance_path_for,
    scan_classifications,
    summarize_classifications,
    write_provenance,
)
from metadarkmatter.core.read_export import extract_reads, normalize_read_id
from metadarkmatter.core.runtime import is_dry_run
from metadarkmatter.models.classification import TaxonomicSummary
from metadarkmatter.models.config import ScoringConfig

# Threshold presets for common use-cases.
# Based on Parks et al. 2018, 2020 and GTDB methodology.
THRESHOLD_PRESETS: dict[str, ScoringConfig] = {
    "default": ScoringConfig(),  # built-in defaults (nucleotide mode)
    "gtdb-strict": ScoringConfig(
        # 95% species boundary (GTDB standard)
        novelty_known_max=5.0,
        novelty_novel_species_min=5.0,
        novelty_novel_species_max=15.0,
        novelty_novel_genus_min=15.0,
        novelty_novel_genus_max=22.0,  # ~78% (GTDB lower genus range)
        # Strict alignment quality (GTDB requires AF >= 50%)
        min_alignment_length=100,
        min_alignment_fraction=0.5,
    ),
    "gtdb-relaxed": ScoringConfig(
        # 97% upper bound (GTDB permits up to 97% for name retention)
        novelty_known_max=3.0,
        novelty_novel_species_min=3.0,
        novelty_novel_species_max=15.0,
        novelty_novel_genus_min=15.0,
        novelty_novel_genus_max=25.0,
        # Relaxed alignment quality
        min_alignment_length=50,
        min_alignment_fraction=0.3,
    ),
    "conservative": ScoringConfig(
        # Conservative genus boundary (more restrictive)
        novelty_novel_species_max=15.0,
        novelty_novel_genus_min=15.0,
        novelty_novel_genus_max=25.0,
        # Standard alignment quality
        min_alignment_length=100,
    ),
    # Literature-strict preset based on comprehensive literature review
    # See docs/CLASSIFICATION_STATISTICS.md for references:
    # - Jain et al. 2018 (Nature Communications): 95-96% ANI species boundary
    # - Riesco & Trujillo 2024: Genus delineation considerations
    # - GTDB standards (Parks et al. 2020)
    "literature-strict": ScoringConfig(
        # 85% identity = conservative novel species cutoff
        novelty_novel_species_max=15.0,
        novelty_novel_genus_min=15.0,
        # 78% identity = genus boundary (more conservative than default 75%)
        novelty_novel_genus_max=22.0,
        # Higher confidence threshold
        confidence_threshold=60.0,
        # GTDB-compatible alignment requirements
        min_alignment_length=100,
        min_alignment_fraction=0.5,
    ),
    "coverage-strict": ScoringConfig(
        coverage_weight_mode="sigmoid",
        coverage_weight_strength=0.7,
        novelty_novel_species_max=15.0,
        novelty_novel_genus_min=15.0,
        novelty_novel_genus_max=25.0,
    ),
    # Adaptive preset: intended for use with --adaptive-thresholds to auto-detect
    # species boundary via GMM. Starts from conservative defaults.
    "adaptive": ScoringConfig(
        min_alignment_length=100,
        min_alignment_fraction=0.3,
    ),
}


def validate_output_format_extension(
    output: Path,
    output_format: str,
    console: Console | QuietConsole,
) -> Path:
    """
    Validate that output file extension matches the specified format.

    If there's a mismatch, warns the user and returns the corrected path.

    Args:
        output: Original output path
        output_format: Specified format ('csv' or 'parquet')
        console: Console for printing warnings

    Returns:
        Corrected output path with appropriate extension
    """
    actual_ext = output.suffix.lower()

    # Check for common mismatches
    if output_format == "parquet" and actual_ext in (".csv", ".tsv"):
        corrected = output.with_suffix(".parquet")
        console.print(
            f"[yellow]Warning: Output extension '{actual_ext}' doesn't match "
            f"format 'parquet'[/yellow]\n"
            f"  Correcting to: {corrected.name}"
        )
        return corrected

    if output_format == "csv" and actual_ext == ".parquet":
        corrected = output.with_suffix(".csv")
        console.print(
            f"[yellow]Warning: Output extension '.parquet' doesn't match "
            f"format 'csv'[/yellow]\n"
            f"  Correcting to: {corrected.name}"
        )
        return corrected

    return output


app = typer.Typer(
    name="score",
    help="Classify BLAST results using ANI-weighted placement",
    no_args_is_help=True,
)

console = Console()


@dataclass(frozen=True)
class _ClassifyOptions:
    """CLI option values that feed scoring-config resolution (modes lowercased)."""

    alignment_mode: str
    bitscore_threshold: float
    min_alignment_length: int
    min_alignment_fraction: float
    max_evalue: float
    min_percent_identity: float
    min_bitscore: float
    min_read_length: int
    min_query_coverage: float
    coverage_weight_mode: str
    coverage_weight_strength: float
    uncertainty_mode: str
    single_hit_uncertainty_threshold: float
    target_family: str | None
    family_ratio_threshold: float
    include_legacy_scores: bool
    bayesian: bool
    config_file: Path | None
    preset: str | None
    clade_novel_species: bool | None = None
    single_hit_novel_species: bool | None = None


def _is_gzip(path: Path) -> bool:
    """True when the file starts with the gzip magic bytes."""
    try:
        with path.open("rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def _explicit_flag_overrides(opts: _ClassifyOptions) -> dict[str, Any]:
    """Return the ScoringConfig fields whose CLI flags differ from their defaults.

    Used to layer explicitly-set flags on top of a ``--config`` YAML or a
    ``--preset`` without discarding fields that only the YAML/preset sets.
    """
    overrides: dict[str, Any] = {}
    if opts.alignment_mode != "nucleotide":
        overrides["alignment_mode"] = opts.alignment_mode
    if opts.bitscore_threshold != 95.0:
        overrides["bitscore_threshold_pct"] = opts.bitscore_threshold
    if opts.min_alignment_length != 100:
        overrides["min_alignment_length"] = opts.min_alignment_length
    if opts.min_alignment_fraction != 0.3:
        overrides["min_alignment_fraction"] = opts.min_alignment_fraction
    if opts.coverage_weight_mode != "linear":
        overrides["coverage_weight_mode"] = opts.coverage_weight_mode
    if opts.coverage_weight_strength != 0.5:
        overrides["coverage_weight_strength"] = opts.coverage_weight_strength
    if opts.uncertainty_mode != "second":
        overrides["uncertainty_mode"] = opts.uncertainty_mode
    if opts.single_hit_uncertainty_threshold != 10.0:
        overrides["single_hit_uncertainty_threshold"] = opts.single_hit_uncertainty_threshold
    if opts.target_family is not None:
        overrides["target_family"] = opts.target_family
    if opts.family_ratio_threshold != 0.8:
        overrides["family_ratio_threshold"] = opts.family_ratio_threshold
    if opts.max_evalue != 0.0:
        overrides["max_evalue"] = opts.max_evalue
    if opts.min_percent_identity != 0.0:
        overrides["min_percent_identity"] = opts.min_percent_identity
    if opts.min_bitscore != 0.0:
        overrides["min_bitscore"] = opts.min_bitscore
    if opts.min_read_length != 0:
        overrides["min_read_length"] = opts.min_read_length
    if opts.min_query_coverage != 0.0:
        overrides["min_query_coverage"] = opts.min_query_coverage
    if opts.include_legacy_scores:
        overrides["include_legacy_scores"] = True
    if opts.clade_novel_species is not None:
        overrides["clade_novel_species"] = opts.clade_novel_species
    if opts.single_hit_novel_species is not None:
        overrides["single_hit_novel_species"] = opts.single_hit_novel_species
    return overrides


def build_classify_config(opts: _ClassifyOptions) -> tuple[ScoringConfig, list[str]]:
    """Resolve the scoring config with precedence --config > flags > preset > defaults.

    Pure aside from reading the optional ``--config`` YAML: returns the resolved
    ``ScoringConfig`` plus a list of human-readable info/deprecation messages
    (with Rich markup) for the caller to print. Raises
    :class:`~metadarkmatter.core.exceptions.ConfigurationError` for an unknown
    preset; a malformed/missing config file surfaces the underlying error.
    """
    messages: list[str] = []

    # Deprecation warnings for flags superseded by --config YAML.
    deprecated = {
        "--bitscore-threshold": (opts.bitscore_threshold, 95.0),
        "--min-alignment-length": (opts.min_alignment_length, 100),
        "--min-alignment-fraction": (opts.min_alignment_fraction, 0.3),
        "--coverage-weight-mode": (opts.coverage_weight_mode, "linear"),
        "--coverage-weight-strength": (opts.coverage_weight_strength, 0.5),
        "--uncertainty-mode": (opts.uncertainty_mode, "second"),
        "--single-hit-threshold": (opts.single_hit_uncertainty_threshold, 10.0),
        "--family-ratio-threshold": (opts.family_ratio_threshold, 0.8),
        "--max-evalue": (opts.max_evalue, 0.0),
        "--min-percent-identity": (opts.min_percent_identity, 0.0),
        "--min-bitscore": (opts.min_bitscore, 0.0),
        "--min-read-length": (opts.min_read_length, 0),
        "--min-query-coverage": (opts.min_query_coverage, 0.0),
    }
    if opts.config_file:
        for flag_name, (val, default) in deprecated.items():
            if val != default:
                messages.append(
                    f"[yellow]Warning: {flag_name} is deprecated when using --config. "
                    f"CLI flags override YAML values.[/yellow]"
                )

    if opts.bayesian:
        messages.append(
            "[yellow]Warning: --bayesian is deprecated. Bayesian posteriors are now "
            "computed by default. Use --include-legacy-scores for the old sub-scores.[/yellow]"
        )

    # Initialize configuration: --config > CLI flags > preset > defaults
    if opts.config_file:
        config = ScoringConfig.from_yaml(opts.config_file)
        messages.append(f"[dim]Loaded config from: {opts.config_file}[/dim]")

        # CLI flags override YAML values
        overrides = _explicit_flag_overrides(opts)
        if overrides:
            config = config.model_copy(update=overrides)
        return config, messages

    if opts.preset:
        preset_lower = opts.preset.lower()
        if preset_lower not in THRESHOLD_PRESETS:
            raise ConfigurationError(
                message=f"Unknown preset '{opts.preset}'",
                suggestion=f"Available presets: {', '.join(THRESHOLD_PRESETS)}.",
            )
        config = THRESHOLD_PRESETS[preset_lower]
        messages.append(f"[dim]Using preset: {preset_lower}[/dim]")
        # Layer explicitly-set flags on top of the preset; model_copy keeps
        # preset-only fields (Bayesian settings, gap thresholds, ...).
        preset_overrides = _explicit_flag_overrides(opts)
        if preset_overrides:
            config = config.model_copy(update=preset_overrides)
        return config, messages

    config = ScoringConfig(
        alignment_mode=opts.alignment_mode,
        bitscore_threshold_pct=opts.bitscore_threshold,
        min_alignment_length=opts.min_alignment_length,
        min_alignment_fraction=opts.min_alignment_fraction,
        max_evalue=opts.max_evalue,
        min_percent_identity=opts.min_percent_identity,
        min_bitscore=opts.min_bitscore,
        min_read_length=opts.min_read_length,
        min_query_coverage=opts.min_query_coverage,
        coverage_weight_mode=opts.coverage_weight_mode,
        coverage_weight_strength=opts.coverage_weight_strength,
        uncertainty_mode=opts.uncertainty_mode,
        single_hit_uncertainty_threshold=opts.single_hit_uncertainty_threshold,
        target_family=opts.target_family,
        family_ratio_threshold=opts.family_ratio_threshold,
        include_legacy_scores=opts.include_legacy_scores,
        **(
            {"clade_novel_species": opts.clade_novel_species}
            if opts.clade_novel_species is not None
            else {}
        ),
        **(
            {"single_hit_novel_species": opts.single_hit_novel_species}
            if opts.single_hit_novel_species is not None
            else {}
        ),
    )
    return config, messages


def _cli_command() -> list[str]:
    """Command line recorded in provenance: program basename plus arguments."""
    return [Path(sys.argv[0]).name, *sys.argv[1:]] if sys.argv else []


def _build_run_provenance(
    request: ClassificationRequest,
    result: ClassificationRunResult,
    *,
    summary_path: Path | None,
) -> dict[str, Any]:
    """Provenance dict for one finished classification run.

    ``result.config`` is used (not ``request.config``) so adaptive thresholds
    and the inferred target family are recorded as applied.
    """
    table: pl.DataFrame | pl.LazyFrame
    if result.classification_df is not None:
        table = result.classification_df
    else:
        # Streaming mode: summarise the on-disk output lazily.
        table = scan_classifications(request.output)
    return build_provenance(
        _cli_command(),
        {
            "alignment": request.alignment,
            "ani": request.ani,
            "metadata": request.metadata,
            "aai": request.aai,
            "id_mapping": request.id_mapping,
            "genomes": request.genomes,
        },
        result.config,
        {
            "outputs": {
                "classifications": str(request.output.resolve()),
                "summary": str(summary_path.resolve()) if summary_path else None,
            },
            "results": summarize_classifications(table),
        },
    )


def _write_summary_json(
    summary_obj: TaxonomicSummary, path: Path, provenance: dict[str, Any] | None
) -> None:
    """Write the summary JSON, embedding provenance under ``"provenance"``."""
    data = json.loads(summary_obj.model_dump_json())
    if provenance is not None:
        data["provenance"] = provenance
    path.write_text(json.dumps(data, indent=2) + "\n")



@app.command(name="classify")
def classify(
    alignment: Path = typer.Option(
        ...,
        "--alignment",
        help="Path to alignment file (BLAST or MMseqs2 tabular format, .tsv or .tsv.gz)",
        exists=True,
        dir_okay=False,
    ),
    ani: Path = typer.Option(
        ...,
        "--ani",
        "-a",
        help="Path to ANI matrix file (CSV or TSV)",
        exists=True,
        dir_okay=False,
    ),
    aai: Path | None = typer.Option(
        None,
        "--aai",
        help=(
            "Path to AAI matrix file (CSV) for genus-level classification. "
            "AAI provides more reliable genus boundaries than ANI at high divergence."
        ),
        exists=True,
        dir_okay=False,
    ),
    output: Path = typer.Option(
        ...,
        "--output",
        "-o",
        help="Output path for classification results (CSV)",
    ),
    summary: Path | None = typer.Option(
        None,
        "--summary",
        "-s",
        help="Output path for summary statistics (JSON)",
    ),
    metadata: Path | None = typer.Option(
        None,
        "--metadata",
        "-m",
        help="Path to genome_metadata.tsv for species-level aggregation",
        exists=True,
        dir_okay=False,
    ),
    genomes: Path | None = typer.Option(
        None,
        "--genomes",
        "-g",
        help="Genome FASTA directory for auto-generating ID mapping (external BLAST)",
        exists=True,
        file_okay=False,
    ),
    id_mapping: Path | None = typer.Option(
        None,
        "--id-mapping",
        help="Pre-computed ID mapping TSV file (alternative to --genomes)",
        exists=True,
        dir_okay=False,
    ),
    bitscore_threshold: float = typer.Option(
        95.0,
        "--bitscore-threshold",
        help="Percentage of top bitscore for ambiguous hits",
        min=0,
        max=100,
    ),
    preset: str | None = typer.Option(
        None,
        "--preset",
        help=(
            "Use predefined threshold preset: 'default', 'gtdb-strict' (95% ANI, AF>=50%), "
            "'gtdb-relaxed' (97% ANI), 'conservative' (narrower genus range), "
            "'literature-strict' (narrow genus, high confidence, AF>=50%), "
            "'coverage-strict' (sigmoid coverage weighting), "
            "'adaptive' (for use with --adaptive-thresholds)"
        ),
    ),
    alignment_mode: str = typer.Option(
        "nucleotide",
        "--alignment-mode",
        help=(
            "Alignment type: 'nucleotide' for BLASTN results (default), "
            "'protein' for BLASTX results. Protein mode uses wider novelty "
            "thresholds (10/25/40%) because amino acid sequences are more "
            "conserved than nucleotide sequences."
        ),
    ),
    min_alignment_length: int = typer.Option(
        100,
        "--min-alignment-length",
        help="Minimum alignment length in bp (0 = no filter)",
        min=0,
    ),
    min_alignment_fraction: float = typer.Option(
        0.3,
        "--min-alignment-fraction",
        help=(
            "Minimum fraction of read aligned (like GTDB's AF). "
            "Set to 0.5 for GTDB-compatible filtering. Default 0.3 filters "
            "short conserved domain hits."
        ),
        min=0.0,
        max=1.0,
    ),
    max_evalue: float = typer.Option(
        0.0,
        "--max-evalue",
        help="Maximum E-value for hit filtering. 0 disables (default).",
        min=0.0,
    ),
    min_percent_identity: float = typer.Option(
        0.0,
        "--min-percent-identity",
        help="Minimum percent identity for hit filtering. 0 disables (default).",
        min=0.0,
        max=100.0,
    ),
    min_bitscore: float = typer.Option(
        0.0,
        "--min-bitscore",
        help="Minimum bitscore for hit filtering. 0 disables (default).",
        min=0.0,
    ),
    min_read_length: int = typer.Option(
        0,
        "--min-read-length",
        help="Minimum read length in bp. 0 disables (default).",
        min=0,
    ),
    min_query_coverage: float = typer.Option(
        0.0,
        "--min-query-coverage",
        help="Minimum query coverage percentage (0-100). 0 disables (default).",
        min=0.0,
        max=100.0,
    ),
    coverage_weight_mode: str = typer.Option(
        "linear",
        "--coverage-weight-mode",
        help=(
            "Coverage weighting mode for hit selection: 'linear' (default, gradual weight increase), "
            "'none' (raw bitscore only), 'log' (diminishing returns), "
            "'sigmoid' (sharp 60%% threshold). Prioritizes longer alignments over short conserved domains."
        ),
    ),
    coverage_weight_strength: float = typer.Option(
        0.5,
        "--coverage-weight-strength",
        help=(
            "Coverage weight strength (0.0-1.0). Controls magnitude of coverage effect. "
            "Weight range = [1-strength, 1+strength]. At 0.5: weights range 0.5x to 1.5x."
        ),
        min=0.0,
        max=1.0,
    ),
    uncertainty_mode: str = typer.Option(
        "second",
        "--uncertainty-mode",
        help=(
            "Mode for calculating placement uncertainty: "
            "'max' uses maximum ANI to any competing genome, "
            "'second' uses ANI to the second-best genome only (default)."
        ),
    ),
    single_hit_uncertainty_threshold: float = typer.Option(
        10.0,
        "--single-hit-threshold",
        help=(
            "Inferred uncertainty threshold for single-hit Ambiguous classification "
            "in legacy_call (taxonomic_call does not use it). "
            "Default 10%% means reads with novelty ~7%% or higher become Ambiguous."
        ),
        min=0.0,
        max=100.0,
    ),
    output_format: str = typer.Option(
        "csv",
        "--format",
        "-f",
        help="Output format: 'csv' or 'parquet'. Parquet is 10x smaller/faster.",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Enable verbose output",
    ),
    streaming: bool = typer.Option(
        False,
        "--streaming",
        help="Use streaming mode for very large files (100M+ alignments)",
    ),
    chunk_size: int = typer.Option(
        1_000_000,
        "--chunk-size",
        help=(
            "Rows per streaming batch (default 1,000,000). Larger values "
            "trade memory for fewer batches; smaller values cap RAM use "
            "for memory-constrained machines. Range 100 to 100,000,000."
        ),
        min=100,
        max=100_000_000,
    ),
    strict_ani: bool = typer.Option(
        False,
        "--strict-ani",
        help=(
            "Refuse to load an ANI matrix whose forward and reverse entries "
            "disagree by more than 0.5 ANI units. Default is to warn and "
            "continue, which matches historical behaviour. Use this on "
            "untrusted matrices to catch silent overwrites of placement "
            "uncertainty inputs."
        ),
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Validate inputs and show what would be processed without running",
    ),
    qc_output: Path | None = typer.Option(
        None,
        "--qc-output",
        help="Output path for QC metrics (JSON). Computes alignment statistics, genome coverage, and classification quality checks.",
    ),
    adaptive_thresholds: bool = typer.Option(
        False,
        "--adaptive-thresholds",
        help=(
            "Detect species boundary from ANI matrix distribution using "
            "Gaussian Mixture Model. Overrides novelty_known_max with "
            "detected value. Requires scikit-learn."
        ),
    ),
    bayesian: bool = typer.Option(
        False,
        "--bayesian",
        help=(
            "Deprecated, no effect: Bayesian classification and the posterior "
            "columns (p_known_species ... p_unclassified) are always included."
        ),
    ),
    target_family: str | None = typer.Option(
        None,
        "--target-family",
        help=(
            "Target taxonomic family for family validation (e.g., 'f__Francisellaceae'). "
            "Enables off-target detection: reads with better hits outside the ANI matrix "
            "are classified as Off-target. If not set but --metadata is provided, the "
            "most common family is inferred."
        ),
    ),
    family_ratio_threshold: float = typer.Option(
        0.8,
        "--family-ratio-threshold",
        help=(
            "Bitscore ratio threshold for off-target detection. "
            "Reads with best_in_family / best_overall bitscore below this value "
            "are classified as Off-target. Default 0.8."
        ),
        min=0.0,
        max=1.0,
    ),
    config_file: Path | None = typer.Option(
        None,
        "--config",
        help=(
            "Path to YAML configuration file. Overrides default parameters. "
            "Generate an editable template with: score export-config"
        ),
        exists=True,
        dir_okay=False,
    ),
    include_legacy_scores: bool = typer.Option(
        False,
        "--include-legacy-scores",
        help=(
            "Include legacy sub-scores in output (alignment_quality, "
            "identity_confidence, placement_confidence, discovery_score). "
            "Off by default in the Bayesian-primary workflow."
        ),
    ),
    clade_novel_species: bool | None = typer.Option(
        None,
        "--clade-novel-species/--no-clade-novel-species",
        help=(
            "Stage 2 rule (nucleotide mode, on by default): call Ambiguous/"
            "Species Boundary reads Novel Species when all competitors are "
            "within one genus and the read is divergent from both the best and "
            "the second-best reference. Adds a clade_novel_species column. "
            "Default comes from --config/--preset (on). See docs/METHODS.md."
        ),
    ),
    single_hit_novel_species: bool | None = typer.Option(
        None,
        "--single-hit-novel-species/--no-single-hit-novel-species",
        help=(
            "Stage 2 rule (nucleotide mode, on by default): call Known Species "
            "reads with no competing reference Novel Species when their novelty "
            "is >= the Known Species threshold (4%%). Adds a "
            "single_hit_novel_species column. See docs/METHODS.md."
        ),
    ),
    provenance: bool = typer.Option(
        True,
        "--provenance/--no-provenance",
        help=(
            "Write <output>.provenance.json next to the classifications "
            "(tool version, command, input SHA-256 hashes, resolved config, "
            "result counts) and embed it in --summary. On by default."
        ),
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet",
        "-q",
        help="Suppress progress output (for scripting)",
    ),
) -> None:
    """
    Classify metagenomic reads using ANI-weighted placement uncertainty.

    This command takes BLAST competitive recruitment results and calculates
    novelty index and placement uncertainty metrics to detect novel bacterial
    diversity in environmental samples.

    For external BLAST results (run outside metadarkmatter), use --genomes
    or --id-mapping to transform contig IDs to the expected format.

    Example:

        metadarkmatter score classify \\
            --alignment results.blast.tsv.gz \\
            --ani ani_matrix.csv \\
            --output classifications.csv \\
            --summary summary.json

        # With species-level aggregation (requires genome_metadata.tsv):
        metadarkmatter score classify \\
            --alignment results.blast.tsv.gz \\
            --ani ani_matrix.csv \\
            --metadata genome_metadata.tsv \\
            --output classifications.csv \\
            --summary summary.json

        # External BLAST results - auto-generate ID mapping from genomes:
        metadarkmatter score classify \\
            --alignment external_results.blast.tsv.gz \\
            --ani ani_matrix.csv \\
            --genomes genomes/ \\
            --output classifications.csv

        # External BLAST results - use pre-computed ID mapping:
        metadarkmatter score classify \\
            --alignment external_results.blast.tsv.gz \\
            --ani ani_matrix.csv \\
            --id-mapping id_mapping.tsv \\
            --output classifications.csv

        # Use Parquet for 10x smaller files and faster I/O:
        metadarkmatter score classify \\
            --alignment results.blast.tsv.gz \\
            --ani ani_matrix.csv \\
            --output classifications.parquet \\
            --format parquet

        # Use GTDB-compatible settings (95% ANI, AF >= 50%):
        metadarkmatter score classify \\
            --alignment results.blast.tsv.gz \\
            --ani ani_matrix.csv \\
            --preset gtdb-strict \\
            --output classifications.csv

        # Custom alignment quality filters:
        metadarkmatter score classify \\
            --alignment results.blast.tsv.gz \\
            --ani ani_matrix.csv \\
            --min-alignment-length 150 \\
            --min-alignment-fraction 0.5 \\
            --output classifications.csv
    """
    # Create console wrapper (respects quiet mode)
    out = QuietConsole(console, quiet=quiet)

    out.print("\n[bold blue]Metadarkmatter ANI-Weighted Classification[/bold blue]\n")

    # Validate format option
    output_format = output_format.lower()
    if output_format not in ("csv", "parquet"):
        console.print(
            f"[red]Error: Invalid format '{output_format}'. "
            f"Use 'csv' or 'parquet'.[/red]"
        )
        raise typer.Exit(code=1) from None
    # Narrowed to the literal type after the membership check above.
    output_format = cast(OutputFormat, output_format)

    # Validate output extension matches format
    output = validate_output_format_extension(output, output_format, out)

    if verbose:
        mode_label = "streaming" if streaming else "vectorized"
        out.print(f"[dim]Processing mode: {mode_label}[/dim]")

    # Validate inputs
    if not alignment.exists():
        console.print(f"[red]Error: Alignment file not found: {alignment}[/red]")
        raise typer.Exit(code=1) from None

    if not ani.exists():
        console.print(f"[red]Error: ANI matrix file not found: {ani}[/red]")
        raise typer.Exit(code=1) from None

    # Dry-run mode: validate and show what would be processed
    if dry_run or is_dry_run():
        console.print("[bold cyan]DRY RUN MODE[/bold cyan] - Validating inputs only\n")

        # Show alignment file info
        alignment_size = alignment.stat().st_size
        console.print(f"[bold]Alignment File:[/bold] {alignment}")
        console.print(f"  Size: {alignment_size / 1024 / 1024:.1f} MB")

        # Count lines (estimate reads)
        try:
            opener = gzip.open if _is_gzip(alignment) else open
            with opener(alignment, "rb") as f:
                line_count = sum(1 for _ in f)
            console.print(f"  Alignments: ~{line_count:,}")
        except (OSError, EOFError):
            console.print("  Alignments: (unable to count)")

        # Show ANI matrix info
        console.print(f"\n[bold]ANI Matrix:[/bold] {ani}")
        try:
            ani_matrix = ANIMatrix.from_file(
                ani, symmetry_check="strict" if strict_ani else "warn"
            )
            console.print(f"  Genomes: {len(ani_matrix.genomes)}")
        except (FileNotFoundError, ValueError, pl.exceptions.PolarsError, OSError) as e:
            console.print(f"  [red]Error loading: {e}[/red]")
            raise typer.Exit(code=1) from None

        # Show genome coverage
        try:
            matched, total, coverage_pct, _missing = validate_ani_genome_coverage(
                alignment, ani_matrix
            )
            console.print("\n[bold]Genome Coverage:[/bold]")
            console.print(f"  Matched: {matched}/{total} ({coverage_pct:.1f}%)")
            if coverage_pct < 50.0:
                console.print("  [yellow]Warning: Low coverage[/yellow]")
        except (FileNotFoundError, pl.exceptions.PolarsError, OSError):
            logger.debug("Could not validate genome coverage", exc_info=True)

        # Show output configuration
        console.print("\n[bold]Output:[/bold]")
        console.print(f"  File: {output}")
        console.print(f"  Format: {output_format}")
        console.print(f"  Mode: {'streaming' if streaming else 'vectorized'}")
        if summary:
            console.print(f"  Summary: {summary}")

        console.print("\n[green]Validation complete. Ready to process.[/green]")
        raise typer.Exit(code=0)

    # Create output directory if needed
    output.parent.mkdir(parents=True, exist_ok=True)

    if summary:
        summary.parent.mkdir(parents=True, exist_ok=True)

    # Validate alignment mode
    alignment_mode_lower = alignment_mode.lower()
    if alignment_mode_lower not in ("nucleotide", "protein"):
        console.print(
            f"[red]Error: Unknown alignment mode '{alignment_mode}'. "
            f"Use 'nucleotide' or 'protein'.[/red]"
        )
        raise typer.Exit(code=1) from None

    # Log alignment mode if protein (non-default)
    if alignment_mode_lower == "protein":
        out.print("[bold]Protein mode enabled[/bold] - using wider novelty thresholds")
        out.print(
            "[dim]  Known Species: N < 10%, Novel Species: 10-25%, Novel Genus: 25-40%[/dim]"
        )
        out.print(
            "[dim]  Note: Ensure input is from BLASTX (DNA query vs protein DB)[/dim]"
        )

    # Validate uncertainty mode
    uncertainty_mode_lower = uncertainty_mode.lower()
    if uncertainty_mode_lower not in ("max", "second"):
        console.print(
            f"[red]Error: Unknown uncertainty mode '{uncertainty_mode}'. "
            f"Use 'max' or 'second'.[/red]"
        )
        raise typer.Exit(code=1) from None

    # Log uncertainty mode if non-default
    if uncertainty_mode_lower == "max":
        out.print(
            "[bold]Uncertainty mode: max[/bold] - using maximum ANI to any competing genome"
        )

    # Resolve scoring config (--config > flags > preset > defaults) and emit
    # any deprecation/info messages it produced.
    config, config_messages = build_classify_config(
        _ClassifyOptions(
            alignment_mode=alignment_mode_lower,
            bitscore_threshold=bitscore_threshold,
            min_alignment_length=min_alignment_length,
            min_alignment_fraction=min_alignment_fraction,
            max_evalue=max_evalue,
            min_percent_identity=min_percent_identity,
            min_bitscore=min_bitscore,
            min_read_length=min_read_length,
            min_query_coverage=min_query_coverage,
            coverage_weight_mode=coverage_weight_mode,
            coverage_weight_strength=coverage_weight_strength,
            uncertainty_mode=uncertainty_mode_lower,
            single_hit_uncertainty_threshold=single_hit_uncertainty_threshold,
            target_family=target_family,
            family_ratio_threshold=family_ratio_threshold,
            include_legacy_scores=include_legacy_scores,
            bayesian=bayesian,
            config_file=config_file,
            preset=preset,
            clade_novel_species=clade_novel_species,
            single_hit_novel_species=single_hit_novel_species,
        )
    )
    for _msg in config_messages:
        out.print(_msg)

    # Build the request and execute the load -> classify -> finalize pipeline.
    # Status text is emitted via out.print; the streaming progress bar is owned
    # here and kept alive across the call (run_classification drives it through
    # the progress_callback). Input-validation and classification errors
    # (Polars, file-system, memory) propagate to the centralized CLI error
    # handler (cli/errors.py), which renders a friendly message and shows a full
    # traceback only under --debug.
    request = ClassificationRequest(
        alignment=alignment,
        ani=ani,
        output=output,
        config=config,
        output_format=output_format,
        aai=aai,
        metadata=metadata,
        genomes=genomes,
        id_mapping=id_mapping,
        streaming=streaming,
        chunk_size=chunk_size,
        strict_ani=strict_ani,
        adaptive_thresholds=adaptive_thresholds,
        compute_qc=qc_output is not None,
    )

    if streaming:
        # Rich progress with ETA; estimate total rows from file size (~100 bytes
        # per alignment).
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            TimeElapsedColumn(),
            TextColumn("•"),
            TimeRemainingColumn(),
            console=console,
            refresh_per_second=2,
        ) as progress:
            estimated_rows = alignment.stat().st_size // 100
            task = progress.add_task(
                "Classifying reads (streaming)...",
                total=estimated_rows,
            )

            def streaming_progress(rows: int, reads: int, elapsed: float) -> None:
                """Update progress bar during streaming."""
                progress.update(task, completed=rows)
                rate = rows / elapsed if elapsed > 0 else 0
                progress.update(
                    task,
                    description=f"Streaming: {reads:,} reads @ {rate:,.0f} rows/s",
                )

            result = run_classification(
                request,
                log=out.print,
                verbose=verbose,
                progress_callback=streaming_progress,
            )
    else:
        result = run_classification(request, log=out.print, verbose=verbose)

    num_classified = result.num_classified
    classification_df = result.classification_df
    qc_metrics = result.qc_metrics
    summary_written = False
    # The output was just rewritten; drop companions from an earlier run that
    # this run will not rewrite, so they cannot describe a different result.
    if summary and num_classified == 0:
        summary.unlink(missing_ok=True)
    if not provenance:
        provenance_path_for(output).unlink(missing_ok=True)

    # Provenance is built before the summary so it can be embedded there.
    prov_data: dict[str, Any] | None = None
    if provenance:
        prov_data = _build_run_provenance(
            request,
            result,
            summary_path=summary if (summary and num_classified > 0) else None,
        )

    # Generate summary if requested
    if summary and num_classified > 0:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console if not quiet else None,
        ) as progress:
            progress.add_task(description="Generating summary...", total=None)

            try:
                # Use in-memory DataFrame if available, otherwise read from output file
                if classification_df is not None:
                    summary_df = classification_df
                else:
                    # Streaming mode: read back from output file
                    summary_df = read_dataframe(output)

                summary_obj = _generate_summary(summary_df)
                _write_summary_json(summary_obj, summary, prov_data)
                summary_written = True

                # Display summary table (skip in quiet mode)
                if not quiet:
                    _display_summary_table(summary_obj)

            except (pl.exceptions.PolarsError, ValueError, OSError) as e:
                if prov_data is not None:
                    prov_data["outputs"]["summary"] = None
                summary.unlink(missing_ok=True)
                console.print(f"\n[yellow]Warning: Failed to generate summary: {e}[/yellow]")
                if verbose:
                    console.print_exception()

    if prov_data is not None:
        prov_file = provenance_path_for(output)
        write_provenance(prov_file, prov_data)
        out.print(f"[bold green]Provenance written to:[/bold green] {prov_file}")

    # Write QC metrics if requested
    if qc_output is not None and qc_metrics is not None:
        qc_output.parent.mkdir(parents=True, exist_ok=True)
        qc_output.write_text(json.dumps(qc_metrics.to_dict(), indent=2))
        out.print(f"[bold green]QC metrics written to:[/bold green] {qc_output}")

        # Display warnings
        if qc_metrics.warnings and not quiet:
            out.print()
            for warning in qc_metrics.warnings:
                out.print(f"[yellow]QC Warning:[/yellow] {warning}")

    out.print(f"\n[bold green]Results written to:[/bold green] {output}")
    if summary_written:
        out.print(f"[bold green]Summary written to:[/bold green] {summary}")
    elif summary:
        out.print("[yellow]No summary written (no reads classified or summary failed).[/yellow]")

    # Point the user at the natural next step.
    out.print(
        f"\n[dim]Next: generate an interactive report with[/dim]\n"
        f"  mdm report generate --classifications {output} --ani {ani} --output report.html"
    )
    out.print()


@app.command(name="export-config")
def export_config(
    output: Path = typer.Option(
        ...,
        "--output",
        "-o",
        help="Output path for the YAML config file",
    ),
    preset: str | None = typer.Option(
        None,
        "--preset",
        "-p",
        help=(
            "Base preset to export: "
            + ", ".join(f"'{k}'" for k in THRESHOLD_PRESETS)
            + ". Defaults to built-in defaults."
        ),
    ),
) -> None:
    """
    Export an editable YAML configuration file.

    Generates a YAML file from the specified preset (or built-in defaults)
    that can be customized and passed to ``score classify --config``.

    Example:

        # Export default config
        metadarkmatter score export-config --output my_config.yaml

        # Export from a preset
        metadarkmatter score export-config --preset gtdb-strict --output gtdb.yaml

        # Edit and use
        metadarkmatter score classify --config my_config.yaml ...
    """
    if preset:
        preset_lower = preset.lower()
        if preset_lower not in THRESHOLD_PRESETS:
            console.print(
                f"[red]Error: Unknown preset '{preset}'. "
                f"Available: {', '.join(THRESHOLD_PRESETS.keys())}[/red]"
            )
            raise typer.Exit(code=1) from None
        config = THRESHOLD_PRESETS[preset_lower]
    else:
        config = ScoringConfig()

    output.parent.mkdir(parents=True, exist_ok=True)
    config.to_yaml(output)

    source = f"preset '{preset_lower}'" if preset else "built-in defaults"
    console.print(f"[green]Config exported from {source} to:[/green] {output}")
    console.print("[dim]Edit the file and use with: score classify --config <path>[/dim]")


# Alignment-file suffixes stripped (longest first) to derive a sample name.
_ALIGNMENT_SUFFIXES: tuple[str, ...] = (
    ".blast.tsv.gz",
    ".mmseqs2.tsv.gz",
    ".blast.tsv",
    ".mmseqs2.tsv",
    ".tsv.gz",
    ".tsv",
)

# Default --pattern set for ``score batch``: BLAST and MMseqs2 outputs, gzipped
# or plain.
DEFAULT_BATCH_PATTERNS: tuple[str, ...] = (
    "*.blast.tsv.gz",
    "*.mmseqs2.tsv.gz",
    "*.blast.tsv",
    "*.mmseqs2.tsv",
)


def sample_name_from_alignment(path: Path) -> str:
    """Derive a sample name by stripping a known alignment suffix.

    ``s1.blast.tsv.gz`` and ``s1.mmseqs2.tsv`` both give ``s1``. Names without
    a known suffix fall back to the file name minus its last extension.
    """
    name = path.name
    for suffix in _ALIGNMENT_SUFFIXES:
        if name.endswith(suffix) and len(name) > len(suffix):
            return name[: -len(suffix)]
    return path.stem


def _find_sample_collisions(files: list[Path]) -> dict[str, list[Path]]:
    """Return sample names claimed by more than one alignment file."""
    by_name: dict[str, list[Path]] = {}
    for f in files:
        by_name.setdefault(sample_name_from_alignment(f), []).append(f)
    return {k: v for k, v in by_name.items() if len(v) > 1}


@app.command(name="batch")
def batch(
    alignment_dir: Path = typer.Option(
        ...,
        "--alignment-dir",
        "-b",
        help="Directory containing alignment files (BLAST or MMseqs2 tabular format)",
        exists=True,
        file_okay=False,
    ),
    ani: Path = typer.Option(
        ...,
        "--ani",
        "-a",
        help="Path to ANI matrix file (CSV or TSV)",
        exists=True,
        dir_okay=False,
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output-dir",
        "-o",
        help="Output directory for classification results",
    ),
    pattern: list[str] | None = typer.Option(
        None,
        "--pattern",
        "-p",
        help=(
            "Glob pattern for alignment files; repeat to give several. Default: "
            + ", ".join(DEFAULT_BATCH_PATTERNS)
        ),
    ),
    aai: Path | None = typer.Option(
        None,
        "--aai",
        help="Path to AAI matrix file (CSV) for genus-level classification",
        exists=True,
        dir_okay=False,
    ),
    metadata: Path | None = typer.Option(
        None,
        "--metadata",
        "-m",
        help="Path to genome_metadata.tsv (species aggregation, representatives, target family)",
        exists=True,
        dir_okay=False,
    ),
    genomes: Path | None = typer.Option(
        None,
        "--genomes",
        "-g",
        help="Genome FASTA directory for auto-generating ID mapping (external alignments)",
        exists=True,
        file_okay=False,
    ),
    id_mapping: Path | None = typer.Option(
        None,
        "--id-mapping",
        help="Pre-computed ID mapping TSV file (alternative to --genomes)",
        exists=True,
        dir_okay=False,
    ),
    bitscore_threshold: float = typer.Option(
        95.0,
        "--bitscore-threshold",
        help="Percentage of top bitscore for ambiguous hits",
        min=0,
        max=100,
    ),
    preset: str | None = typer.Option(
        None,
        "--preset",
        help=(
            "Use predefined threshold preset: "
            + ", ".join(f"'{k}'" for k in THRESHOLD_PRESETS)
        ),
    ),
    config_file: Path | None = typer.Option(
        None,
        "--config",
        help="Path to YAML configuration file (see: score export-config)",
        exists=True,
        dir_okay=False,
    ),
    alignment_mode: str = typer.Option(
        "nucleotide",
        "--alignment-mode",
        help=(
            "Alignment type: 'nucleotide' for BLASTN results (default), "
            "'protein' for BLASTX results with wider thresholds."
        ),
    ),
    min_alignment_length: int = typer.Option(
        100,
        "--min-alignment-length",
        help="Minimum alignment length in bp (0 = no filter)",
        min=0,
    ),
    min_alignment_fraction: float = typer.Option(
        0.3,
        "--min-alignment-fraction",
        help="Minimum fraction of read aligned (GTDB uses 0.5). Default 0.3.",
        min=0.0,
        max=1.0,
    ),
    max_evalue: float = typer.Option(
        0.0,
        "--max-evalue",
        help="Maximum E-value for hit filtering. 0 disables (default).",
        min=0.0,
    ),
    min_percent_identity: float = typer.Option(
        0.0,
        "--min-percent-identity",
        help="Minimum percent identity for hit filtering. 0 disables (default).",
        min=0.0,
        max=100.0,
    ),
    min_bitscore: float = typer.Option(
        0.0,
        "--min-bitscore",
        help="Minimum bitscore for hit filtering. 0 disables (default).",
        min=0.0,
    ),
    min_read_length: int = typer.Option(
        0,
        "--min-read-length",
        help="Minimum read length in bp. 0 disables (default).",
        min=0,
    ),
    min_query_coverage: float = typer.Option(
        0.0,
        "--min-query-coverage",
        help="Minimum query coverage percentage (0-100). 0 disables (default).",
        min=0.0,
        max=100.0,
    ),
    coverage_weight_mode: str = typer.Option(
        "linear",
        "--coverage-weight-mode",
        help=(
            "Coverage weighting mode: 'linear' (default), 'none', 'log', or 'sigmoid'. "
            "Prioritizes longer alignments over short conserved domains."
        ),
    ),
    coverage_weight_strength: float = typer.Option(
        0.5,
        "--coverage-weight-strength",
        help="Coverage weight strength (0.0-1.0). Weight range = [1-strength, 1+strength].",
        min=0.0,
        max=1.0,
    ),
    uncertainty_mode: str = typer.Option(
        "second",
        "--uncertainty-mode",
        help=(
            "Mode for calculating placement uncertainty: "
            "'max' uses maximum ANI to any competing genome, "
            "'second' uses ANI to the second-best genome only (default)."
        ),
    ),
    single_hit_uncertainty_threshold: float = typer.Option(
        10.0,
        "--single-hit-threshold",
        help="Inferred uncertainty threshold for single-hit reads in legacy_call.",
        min=0.0,
        max=100.0,
    ),
    target_family: str | None = typer.Option(
        None,
        "--target-family",
        help=(
            "Target taxonomic family for off-target detection (e.g. "
            "'f__Francisellaceae'). Inferred from --metadata when omitted."
        ),
    ),
    family_ratio_threshold: float = typer.Option(
        0.8,
        "--family-ratio-threshold",
        help="Bitscore ratio threshold for off-target detection. Default 0.8.",
        min=0.0,
        max=1.0,
    ),
    include_legacy_scores: bool = typer.Option(
        False,
        "--include-legacy-scores",
        help="Include legacy sub-scores in the output.",
    ),
    clade_novel_species: bool | None = typer.Option(
        None,
        "--clade-novel-species/--no-clade-novel-species",
        help="Toggle the clade-novelty Stage 2 rule (default from --config/--preset: on).",
    ),
    single_hit_novel_species: bool | None = typer.Option(
        None,
        "--single-hit-novel-species/--no-single-hit-novel-species",
        help="Toggle the single-hit novelty Stage 2 rule (default from --config/--preset: on).",
    ),
    streaming: bool = typer.Option(
        False,
        "--streaming",
        help="Use streaming mode for very large files (no species columns; no --id-mapping)",
    ),
    chunk_size: int = typer.Option(
        1_000_000,
        "--chunk-size",
        help="Rows per streaming batch (default 1,000,000).",
        min=100,
        max=100_000_000,
    ),
    strict_ani: bool = typer.Option(
        False,
        "--strict-ani",
        help="Refuse to load an ANI matrix whose forward and reverse entries disagree.",
    ),
    adaptive_thresholds: bool = typer.Option(
        False,
        "--adaptive-thresholds",
        help="Detect the species boundary from the ANI matrix (requires scikit-learn).",
    ),
    provenance: bool = typer.Option(
        True,
        "--provenance/--no-provenance",
        help="Write <output>.provenance.json for every sample (default: on).",
    ),
    output_format: str = typer.Option(
        "csv",
        "--format",
        "-f",
        help="Output format: 'csv' or 'parquet'. Parquet is 10x smaller/faster.",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Enable verbose output",
    ),
    workers: int = typer.Option(
        0,
        "--workers",
        "-w",
        hidden=True,
        help="Deprecated and ignored; samples are processed sequentially.",
        min=0,
    ),
) -> None:
    """
    Batch classify multiple alignment files.

    Processes every alignment file (BLAST or MMseqs2) matching the pattern(s)
    in the input directory through the same code path as ``score classify``,
    so each sample's output has identical columns to a single classify run
    with the same options. Per sample it writes
    ``<sample>_classifications.csv``, ``<sample>_summary.json`` and (unless
    --no-provenance) ``<sample>_classifications.csv.provenance.json``.

    Sample names come from the file name with ``.blast.tsv[.gz]`` /
    ``.mmseqs2.tsv[.gz]`` / ``.tsv[.gz]`` removed; two files that map to the
    same name are an error (narrow --pattern or rename a file).

    Samples are processed one after another. The matrices and metadata are
    re-read for each sample (they are small next to the alignments), which
    keeps peak memory at one sample's worth and the output order deterministic.

    Example:

        metadarkmatter score batch \\
            --alignment-dir results/alignments/ \\
            --ani ani_matrix.csv \\
            --metadata genome_metadata.tsv \\
            --output-dir results/classifications/

        # Only MMseqs2 outputs, Parquet, GTDB-compatible settings:
        metadarkmatter score batch \\
            --alignment-dir results/alignments/ \\
            --ani ani_matrix.csv \\
            --output-dir results/classifications/ \\
            --pattern "*.mmseqs2.tsv.gz" --format parquet --preset gtdb-strict
    """
    console.print("\n[bold blue]Metadarkmatter Batch Classification[/bold blue]\n")

    if workers:
        console.print(
            "[yellow]Warning: --workers is deprecated and has no effect; samples "
            "are processed sequentially (each sample already uses all cores "
            "through Polars).[/yellow]"
        )

    # Validate format option
    output_format = output_format.lower()
    if output_format not in ("csv", "parquet"):
        console.print(
            f"[red]Error: Invalid format '{output_format}'. "
            f"Use 'csv' or 'parquet'.[/red]"
        )
        raise typer.Exit(code=1) from None
    output_format = cast(OutputFormat, output_format)

    file_ext = ".parquet" if output_format == "parquet" else ".csv"

    # Find all alignment files (union of the patterns, deterministic order)
    patterns = list(pattern) if pattern else list(DEFAULT_BATCH_PATTERNS)
    alignment_files = sorted(
        {f for pat in patterns for f in alignment_dir.glob(pat) if f.is_file()}
    )

    if not alignment_files:
        shown = ", ".join(f"'{p}'" for p in patterns)
        console.print(
            f"[red]Error: No alignment files found matching pattern: {shown}[/red]\n"
            f"  Directory: {alignment_dir}\n"
            f"  Pattern: {shown}\n\n"
            f"Suggestions:\n"
            f"  - Check if the directory contains alignment files\n"
            f"  - Try a different pattern (e.g., '*.tsv', '*.blast.tsv.gz')"
        )
        raise typer.Exit(code=1) from None

    collisions = _find_sample_collisions(alignment_files)
    if collisions:
        console.print("[red]Error: several alignment files map to the same sample name:[/red]")
        for name, files in sorted(collisions.items()):
            console.print(f"  {name}: {', '.join(f.name for f in files)}")
        console.print(
            "  Narrow --pattern (for example '*.mmseqs2.tsv.gz') or rename the files "
            "so that each sample name is unique."
        )
        raise typer.Exit(code=1) from None

    # Option combinations that would fail identically for every sample.
    if genomes is not None and id_mapping is not None:
        raise ConfigurationError(
            "Cannot specify both --genomes and --id-mapping. Use one or the other."
        )
    if streaming and (genomes is not None or id_mapping is not None):
        raise ConfigurationError(
            "--id-mapping/--genomes is not supported with --streaming.",
            suggestion="Remove --streaming, or pre-transform the alignments.",
        )

    # Validate alignment / uncertainty modes
    alignment_mode_lower = alignment_mode.lower()
    if alignment_mode_lower not in ("nucleotide", "protein"):
        console.print(
            f"[red]Error: Unknown alignment mode '{alignment_mode}'. "
            f"Use 'nucleotide' or 'protein'.[/red]"
        )
        raise typer.Exit(code=1) from None
    uncertainty_mode_lower = uncertainty_mode.lower()
    if uncertainty_mode_lower not in ("max", "second"):
        console.print(
            f"[red]Error: Unknown uncertainty mode '{uncertainty_mode}'. "
            f"Use 'max' or 'second'.[/red]"
        )
        raise typer.Exit(code=1) from None

    if alignment_mode_lower == "protein":
        console.print("[bold]Protein mode enabled[/bold] - using wider novelty thresholds")
        console.print(
            "[dim]  Known Species: N < 10%, Novel Species: 10-25%, Novel Genus: 25-40%[/dim]"
        )
        console.print("[dim]  Note: Ensure input is from BLASTX (DNA query vs protein DB)[/dim]")
    if uncertainty_mode_lower == "max":
        console.print(
            "[bold]Uncertainty mode: max[/bold] - using maximum ANI to any competing genome"
        )

    console.print(f"Found {len(alignment_files)} alignment files to process\n")

    # Same config resolution as 'score classify'
    config, config_messages = build_classify_config(
        _ClassifyOptions(
            alignment_mode=alignment_mode_lower,
            bitscore_threshold=bitscore_threshold,
            min_alignment_length=min_alignment_length,
            min_alignment_fraction=min_alignment_fraction,
            max_evalue=max_evalue,
            min_percent_identity=min_percent_identity,
            min_bitscore=min_bitscore,
            min_read_length=min_read_length,
            min_query_coverage=min_query_coverage,
            coverage_weight_mode=coverage_weight_mode,
            coverage_weight_strength=coverage_weight_strength,
            uncertainty_mode=uncertainty_mode_lower,
            single_hit_uncertainty_threshold=single_hit_uncertainty_threshold,
            target_family=target_family,
            family_ratio_threshold=family_ratio_threshold,
            include_legacy_scores=include_legacy_scores,
            bayesian=False,
            config_file=config_file,
            preset=preset,
            clade_novel_species=clade_novel_species,
            single_hit_novel_species=single_hit_novel_species,
        )
    )
    for _msg in config_messages:
        console.print(_msg)

    output_dir.mkdir(parents=True, exist_ok=True)

    def _log(message: str) -> None:
        # Per-sample runner chatter only with --verbose; warnings always.
        if verbose or "arning" in message:
            console.print(message)

    total_classified = 0
    failed_files: list[str] = []

    for i, alignment_file in enumerate(alignment_files, 1):
        sample_name = sample_name_from_alignment(alignment_file)
        output_file = output_dir / f"{sample_name}_classifications{file_ext}"
        summary_file = output_dir / f"{sample_name}_summary.json"

        console.print(f"[{i}/{len(alignment_files)}] Processing {sample_name}...")

        request = ClassificationRequest(
            alignment=alignment_file,
            ani=ani,
            output=output_file,
            config=config,
            output_format=output_format,
            aai=aai,
            metadata=metadata,
            genomes=genomes,
            id_mapping=id_mapping,
            streaming=streaming,
            chunk_size=chunk_size,
            strict_ani=strict_ani,
            adaptive_thresholds=adaptive_thresholds,
        )

        try:
            result = run_classification(request, log=_log, verbose=verbose)
            num_classified = result.num_classified

            prov_data: dict[str, Any] | None = None
            if provenance:
                prov_data = _build_run_provenance(
                    request,
                    result,
                    summary_path=summary_file if num_classified > 0 else None,
                )

            if num_classified > 0:
                summary_df = (
                    result.classification_df
                    if result.classification_df is not None
                    else read_dataframe(output_file)
                )
                _write_summary_json(_generate_summary(summary_df), summary_file, prov_data)
                console.print(f"  [green]Classified {num_classified:,} reads[/green]")
            else:
                summary_file.unlink(missing_ok=True)  # stale from an earlier run
                console.print(
                    "  [yellow]Warning: 0 reads classified (header-only output "
                    "written); check the alignment filters.[/yellow]"
                )
            if prov_data is not None:
                write_provenance(provenance_path_for(output_file), prov_data)
            else:
                provenance_path_for(output_file).unlink(missing_ok=True)
            total_classified += num_classified

        except (
            MetadarkmatterError,
            pl.exceptions.PolarsError,
            FileNotFoundError,
            ValueError,
            OSError,
        ) as e:
            console.print(f"  [red]Failed: {e}[/red]")
            failed_files.append(alignment_file.name)
            if verbose:
                console.print_exception()

    # Summary
    console.print("\n[bold]Batch Processing Complete[/bold]")
    console.print(f"Total reads classified: {total_classified:,}")
    console.print(f"Successful: {len(alignment_files) - len(failed_files)}/{len(alignment_files)}")

    if failed_files:
        console.print("\n[yellow]Failed files:[/yellow]")
        for fname in failed_files:
            console.print(f"  - {fname}")

    console.print(f"\n[bold green]Results written to:[/bold green] {output_dir}\n")

    if failed_files:
        raise typer.Exit(code=1)


def _generate_summary(df: pl.DataFrame) -> TaxonomicSummary:
    """Generate summary statistics from classification DataFrame."""
    total_reads = len(df)

    # Count each classification category
    counts = df.group_by("taxonomic_call").len().to_dict(as_series=False)
    count_dict = dict(zip(counts["taxonomic_call"], counts["len"], strict=False))

    # Get top 50 genome hits
    genome_counts = (
        df.group_by("best_match_genome")
        .len()
        .sort("len", descending=True)
        .head(50)
        .to_dict(as_series=False)
    )
    genome_hit_counts = dict(
        zip(genome_counts["best_match_genome"], genome_counts["len"], strict=False)
    )

    # Get species hit counts if species column is present (from metadata join)
    species_hit_counts: dict[str, int] = {}
    if "species" in df.columns:
        species_counts = (
            df.group_by("species")
            .len()
            .sort("len", descending=True)
            .head(50)
            .to_dict(as_series=False)
        )
        species_hit_counts = dict(
            zip(species_counts["species"], species_counts["len"], strict=False)
        )

    # High-level diversity grouping. These fields were never populated, so
    # the summary JSON always reported 0 known / novel / uncertain reads.
    diversity_counts: dict[str, int] = {}
    if "diversity_status" in df.columns:
        div = df.group_by("diversity_status").len().to_dict(as_series=False)
        diversity_counts = dict(zip(div["diversity_status"], div["len"], strict=False))

    return TaxonomicSummary(
        total_reads=total_reads,
        diversity_known=diversity_counts.get("Known", 0),
        diversity_novel=diversity_counts.get("Novel", 0),
        diversity_uncertain=diversity_counts.get("Uncertain", 0),
        known_species=count_dict.get("Known Species", 0),
        novel_species=count_dict.get("Novel Species", 0),
        novel_genus=count_dict.get("Novel Genus", 0),
        species_boundary=count_dict.get("Species Boundary", 0),
        ambiguous=count_dict.get("Ambiguous", 0),
        ambiguous_within_genus=count_dict.get("Ambiguous Within Genus", 0),
        conserved_regions=count_dict.get("Conserved Region", 0),
        unclassified=count_dict.get("Unclassified", 0),
        off_target=count_dict.get("Off-target", 0),
        mean_novelty_index=df["novelty_index"].mean(),
        mean_placement_uncertainty=df["placement_uncertainty"].mean(),
        genome_hit_counts=genome_hit_counts,
        species_hit_counts=species_hit_counts,
    )


def _display_summary_table(summary: TaxonomicSummary) -> None:
    """Display summary statistics as a Rich table."""
    table = Table(title="Classification Summary", show_header=True)

    table.add_column("Category", style="cyan", no_wrap=True)
    table.add_column("Count", justify="right", style="magenta")
    table.add_column("Percentage", justify="right", style="green")

    total = summary.total_reads or 1  # Avoid division by zero

    table.add_row(
        "Known Species",
        f"{summary.known_species:,}",
        f"{summary.known_species_pct:.2f}%",
    )
    table.add_row(
        "Novel Species",
        f"{summary.novel_species:,}",
        f"{summary.novel_species_pct:.2f}%",
        style="yellow",
    )
    table.add_row(
        "Novel Genus",
        f"{summary.novel_genus:,}",
        f"{summary.novel_genus_pct:.2f}%",
        style="yellow",
    )
    table.add_row(
        "Ambiguous",
        f"{summary.ambiguous:,}",
        f"{100.0 * summary.ambiguous / total:.2f}%",
        style="blue",
    )
    table.add_row(
        "Conserved Regions",
        f"{summary.conserved_regions:,}",
        f"{100.0 * summary.conserved_regions / total:.2f}%",
    )
    table.add_row(
        "Unclassified",
        f"{summary.unclassified:,}",
        f"{100.0 * summary.unclassified / total:.2f}%",
        style="dim",
    )
    if summary.off_target > 0:
        table.add_row(
            "Off-target",
            f"{summary.off_target:,}",
            f"{100.0 * summary.off_target / total:.2f}%",
            style="dim red",
        )
    table.add_section()
    table.add_row(
        "[bold]Total Novel Diversity[/bold]",
        f"[bold]{summary.novel_species + summary.novel_genus:,}[/bold]",
        f"[bold]{summary.novel_diversity_pct:.2f}%[/bold]",
        style="bold yellow",
    )

    console.print()
    console.print(table)
    console.print()

    # Additional metrics
    console.print(f"Mean Novelty Index: {summary.mean_novelty_index:.2f}")
    console.print(f"Mean Placement Uncertainty: {summary.mean_placement_uncertainty:.2f}")

    # Species information if available (from metadata)
    if summary.species_count > 0:
        console.print(f"Unique Species: {summary.species_count}")
        console.print()
        console.print("[bold]Top Species:[/bold]")
        for species, count in list(summary.species_hit_counts.items())[:5]:
            pct = 100.0 * count / summary.total_reads
            console.print(f"  {species}: {count:,} ({pct:.1f}%)")

    console.print()


@app.command(name="extract-novel")
def extract_novel(
    classifications: Path = typer.Option(
        ...,
        "--classifications", "-c",
        help="Path to classification results file (CSV or Parquet)",
        exists=True,
        dir_okay=False,
    ),
    output: Path = typer.Option(
        ...,
        "--output", "-o",
        help="Output path for novel species candidates (CSV)",
    ),
    read_ids: Path | None = typer.Option(
        None,
        "--read-ids", "-r",
        help="Output path for read IDs only (for seqtk extraction)",
    ),
    category: str = typer.Option(
        "all",
        "--category", "-t",
        help="Category to extract: 'species', 'genus', or 'all' (both)",
    ),
    min_novelty: float | None = typer.Option(
        None,
        "--min-novelty",
        help="Minimum novelty index threshold (overrides category defaults)",
        min=0.0,
        max=100.0,
    ),
    max_uncertainty: float = typer.Option(
        2.0,
        "--max-uncertainty",
        help=(
            "Maximum placement uncertainty (exclude ambiguous placements). "
            "Reads flagged clade_novel_species are exempt."
        ),
        min=0.0,
        max=100.0,
    ),
    group_by_genome: bool = typer.Option(
        True,
        "--group-by-genome/--no-group",
        help="Group results by closest reference genome (default: group)",
    ),
    min_reads: int = typer.Option(
        1,
        "--min-reads",
        help="Minimum reads per genome to include in output",
        min=1,
    ),
    reads_1: Path | None = typer.Option(
        None,
        "--reads-1",
        help=(
            "Original FASTQ/FASTA(.gz) (the R1 file, or single-end reads) to pull "
            "the selected reads from. Needs --fastq-out, --fasta-out or --per-lineage-dir."
        ),
        exists=True,
        dir_okay=False,
    ),
    reads_2: Path | None = typer.Option(
        None,
        "--reads-2",
        help="R2 FASTQ/FASTA(.gz). Outputs become {stem}_R1 / {stem}_R2 files.",
        exists=True,
        dir_okay=False,
    ),
    fastq_out: Path | None = typer.Option(
        None,
        "--fastq-out",
        help="Write the selected reads as FASTQ (gzip if the path ends in .gz).",
    ),
    fasta_out: Path | None = typer.Option(
        None,
        "--fasta-out",
        help="Write the selected reads as FASTA (gzip if the path ends in .gz).",
    ),
    per_lineage_dir: Path | None = typer.Option(
        None,
        "--per-lineage-dir",
        help=(
            "Also write one gzipped file per selected lineage (closest reference "
            "genome) into this directory, for per-lineage assembly."
        ),
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose", "-v",
        help="Enable verbose output",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet", "-q",
        help="Suppress progress output",
    ),
) -> None:
    """
    Extract candidate novel species or genera from classification results.

    Filters classification results to identify reads that represent potentially
    novel diversity, grouped by their closest reference genome. This facilitates
    targeted assembly or further validation of novel lineages.

    Categories:
      - species: Novel species candidates (5-20% novelty, <2% uncertainty)
      - genus: Novel genus candidates (20-25% novelty, <2% uncertainty)
      - all: Both novel species and genus candidates

    Example:

        # Extract all novel candidates grouped by genome
        metadarkmatter score extract-novel \\
            --classifications results.csv \\
            --output novel_candidates.csv \\
            --read-ids novel_read_ids.txt

        # Extract only novel genus candidates
        metadarkmatter score extract-novel \\
            --classifications results.csv \\
            --output novel_genera.csv \\
            --category genus

        # Custom novelty threshold
        metadarkmatter score extract-novel \\
            --classifications results.csv \\
            --output highly_novel.csv \\
            --min-novelty 20

    The output can be used with seqtk to extract reads for assembly:

        seqtk subseq reads.fastq.gz novel_read_ids.txt > novel_reads.fastq

    Or let this command stream the original reads itself (gzip-aware, constant
    memory). Read IDs are matched with the /1 and /2 mate suffix removed, so
    IDs from a paired MMseqs2 search match plain FASTQ headers; both mates of a
    selected pair are written:

        metadarkmatter score extract-novel \\
            --classifications results.csv --output novel.csv \\
            --reads-1 sample_R1.fastq.gz --reads-2 sample_R2.fastq.gz \\
            --fastq-out novel.fastq.gz --per-lineage-dir novel_by_lineage/

    This writes novel_R1.fastq.gz and novel_R2.fastq.gz (plus one file per
    lineage in novel_by_lineage/).
    """
    out = QuietConsole(console, quiet=quiet)

    export_requested = any(x is not None for x in (fastq_out, fasta_out, per_lineage_dir))
    if reads_2 is not None and reads_1 is None:
        console.print("[red]Error: --reads-2 requires --reads-1.[/red]")
        raise typer.Exit(code=1) from None
    if export_requested and reads_1 is None:
        console.print(
            "[red]Error: --fastq-out/--fasta-out/--per-lineage-dir require --reads-1.[/red]"
        )
        raise typer.Exit(code=1) from None
    if reads_1 is not None and not export_requested:
        console.print(
            "[red]Error: --reads-1 given but no output requested; add --fastq-out, "
            "--fasta-out or --per-lineage-dir.[/red]"
        )
        raise typer.Exit(code=1) from None

    out.print("\n[bold blue]Metadarkmatter Novel Species Extraction[/bold blue]\n")

    # Validate parameter combinations
    if min_novelty is not None and min_novelty > 50.0:
        console.print(
            f"[yellow]Warning: min_novelty={min_novelty}% is unusually high. "
            f"Typical novel species have 5-25% novelty.[/yellow]"
        )

    if max_uncertainty >= 5.0:
        console.print(
            f"[yellow]Warning: max_uncertainty={max_uncertainty}% is high. "
            f"Results may include conserved region hits.[/yellow]"
        )

    # Validate category
    category = category.lower()
    if category not in ("species", "genus", "all"):
        console.print(
            f"[red]Error: Invalid category '{category}'. "
            f"Use 'species', 'genus', or 'all'.[/red]"
        )
        raise typer.Exit(code=1) from None

    # Create output directory if needed
    output.parent.mkdir(parents=True, exist_ok=True)
    if read_ids:
        read_ids.parent.mkdir(parents=True, exist_ok=True)

    # Load classification results
    out.print("[bold]Step 1:[/bold] Loading classification results...")

    try:
        df = read_dataframe(classifications)
    except (FileNotFoundError, ValueError, pl.exceptions.PolarsError, OSError) as e:
        console.print(f"[red]Error reading classifications: {e}[/red]")
        raise typer.Exit(code=1) from None

    out.print(f"  [green]Loaded {len(df):,} classifications[/green]")

    if df.height == 0:
        # A header-only CSV (written when no read was classified) is read
        # with every column as String; give the filtered columns their types.
        typed = {
            "novelty_index": pl.Float64, "placement_uncertainty": pl.Float64,
            "top_hit_identity": pl.Float64, "clade_novel_species": pl.Boolean,
        }
        df = df.with_columns(
            pl.lit(None, dtype=dtype).alias(col)
            for col, dtype in typed.items() if col in df.columns
        )

    # Build filter expression based on category
    out.print(f"\n[bold]Step 2:[/bold] Filtering for {category} candidates...")

    if category == "species":
        filter_expr = pl.col("taxonomic_call") == "Novel Species"
        filter_desc = "Novel Species"
    elif category == "genus":
        filter_expr = pl.col("taxonomic_call") == "Novel Genus"
        filter_desc = "Novel Genus"
    else:  # all
        filter_expr = (
            (pl.col("taxonomic_call") == "Novel Species") |
            (pl.col("taxonomic_call") == "Novel Genus")
        )
        filter_desc = "Novel Species + Novel Genus"

    # --min-novelty narrows the selected category (intersection); the
    # uncertainty cap applies in every case.
    if min_novelty is not None:
        filter_expr = filter_expr & (pl.col("novelty_index") >= min_novelty)
        filter_desc += f", novelty >= {min_novelty}%"
    # Reads moved to Novel Species by the clade-novelty rule have high U by
    # construction (competitors within one genus, all distant), so U does not
    # measure doubt about their novelty; they are exempt from the cap.
    u_ok = pl.col("placement_uncertainty") < max_uncertainty
    if "clade_novel_species" in df.columns:
        u_ok = u_ok | pl.col("clade_novel_species").fill_null(False)
        filter_desc += f", uncertainty < {max_uncertainty}% (clade-rule reads exempt)"
    else:
        filter_desc += f", uncertainty < {max_uncertainty}%"
    filter_expr = filter_expr & u_ok

    # Apply filter
    novel_df = df.filter(filter_expr)

    if len(novel_df) == 0:
        console.print(f"[yellow]No reads found matching criteria: {filter_desc}[/yellow]")
        console.print("[yellow]Writing header-only output.[/yellow]")

    out.print(f"  [green]Found {len(novel_df):,} novel candidates[/green]")
    out.print(f"  [dim]Filter: {filter_desc}[/dim]")

    # --min-reads keeps reference genomes (lineages) with enough reads; it
    # applies to the grouped table, the per-read table, --read-ids and the
    # read export alike.
    selected_df = novel_df
    if min_reads > 1:
        kept = (
            novel_df.group_by("best_match_genome").len()
            .filter(pl.col("len") >= min_reads)["best_match_genome"]
        )
        selected_df = novel_df.filter(pl.col("best_match_genome").is_in(kept.to_list()))

    # Group by genome if requested
    if group_by_genome:
        out.print("\n[bold]Step 3:[/bold] Grouping by reference genome...")

        # Sort keys and the mode tie-break make the output order and the
        # dominant category reproducible across runs (polars group_by and
        # mode() do not preserve order).
        grouped = selected_df.group_by("best_match_genome").agg([
            pl.len().alias("read_count"),
            pl.col("novelty_index").mean().alias("mean_novelty"),
            pl.col("novelty_index").min().alias("min_novelty"),
            pl.col("novelty_index").max().alias("max_novelty"),
            pl.col("top_hit_identity").mean().alias("mean_identity"),
            pl.col("placement_uncertainty").mean().alias("mean_uncertainty"),
            pl.col("taxonomic_call").mode().sort().first().alias("dominant_category"),
        ]).sort(["read_count", "best_match_genome"], descending=[True, False])

        out.print(f"  [green]Found {len(grouped)} candidate lineages[/green]")

        # Display top candidates
        if not quiet:
            table = Table(title="Top Novel Lineages", show_header=True)
            table.add_column("Reference Genome", style="cyan")
            table.add_column("Reads", justify="right", style="magenta")
            table.add_column("Mean Novelty", justify="right", style="yellow")
            table.add_column("Mean Identity", justify="right", style="green")
            table.add_column("Category", style="blue")

            for row in grouped.head(10).iter_rows(named=True):
                table.add_row(
                    row["best_match_genome"],
                    f"{row['read_count']:,}",
                    f"{row['mean_novelty']:.1f}%",
                    f"{row['mean_identity']:.1f}%",
                    row["dominant_category"],
                )

            console.print()
            console.print(table)

            if len(grouped) > 10:
                console.print(f"  [dim]... and {len(grouped) - 10} more lineages[/dim]")
            console.print()

        # Write grouped output
        grouped.write_csv(output)
        out.print(f"[green]Grouped candidates written to:[/green] {output}")

    else:
        # Write individual reads
        selected_df.write_csv(output)
        out.print(f"[green]Novel reads written to:[/green] {output}")

    # Write read IDs for seqtk extraction
    if read_ids:
        out.print("\n[bold]Step 4:[/bold] Exporting read IDs...")

        selected_df.select("read_id").write_csv(read_ids, include_header=False)
        out.print(f"[green]Read IDs written to:[/green] {read_ids}")
        out.print(f"  [dim]Use with: seqtk subseq reads.fastq.gz {read_ids} > novel.fastq[/dim]")

    # Stream the original reads and write the selected records
    if reads_1 is not None:
        out.print("\n[bold]Step 5:[/bold] Extracting reads from FASTQ/FASTA...")
        # Sorted iteration makes the lineage of a pair whose mates hit
        # different genomes deterministic (the smallest genome name wins).
        lineage_of: dict[str, str] = {}
        for rid, genome in sorted(
            zip(
                selected_df["read_id"].to_list(),
                selected_df["best_match_genome"].to_list(),
                strict=True,
            ),
            key=lambda t: (normalize_read_id(t[0]), str(t[1])),
        ):
            lineage_of.setdefault(normalize_read_id(rid), str(genome))
        wanted = set(lineage_of)
        try:
            report = extract_reads(
                reads_1=reads_1,
                reads_2=reads_2,
                ids=wanted,
                fastq_out=fastq_out,
                fasta_out=fasta_out,
                per_lineage_dir=per_lineage_dir,
                lineage_of=lineage_of,
            )
        except ConfigurationError as e:
            console.print(f"[red]Error: {e}[/red]")
            raise typer.Exit(code=1) from None

        out.print(f"  Requested read IDs:  {report.requested:,}")
        out.print(f"  Found in reads-1:    {report.found_1:,} ({report.written_1:,} records)")
        if report.found_2 is not None:
            out.print(f"  Found in reads-2:    {report.found_2:,} ({report.written_2:,} records)")
        for path_out in [*report.outputs, *report.lineage_files]:
            if verbose:
                out.print(f"  [dim]wrote {path_out}[/dim]")
        for path_out in report.outputs:
            out.print(f"[green]Reads written to:[/green] {path_out}")
        if report.lineage_files:
            out.print(
                f"[green]Per-lineage files:[/green] {len(report.lineage_files)} "
                f"in {per_lineage_dir}"
            )
        for label, missing in (("reads-1", report.missing), ("reads-2", report.missing_2)):
            if missing and (label == "reads-1" or reads_2 is not None):
                console.print(
                    f"[yellow]Warning: {len(missing):,} of {report.requested:,} requested "
                    f"read IDs not found in {label} (first: "
                    f"{', '.join(missing[:5])}{'...' if len(missing) > 5 else ''}). "
                    "Check that the reads are the files that were aligned.[/yellow]"
                )

    # Summary
    out.print("\n[bold green]Extraction complete![/bold green]")
    out.print("\n[bold]Summary:[/bold]")
    out.print(f"  Total novel reads: {len(novel_df):,}")
    if group_by_genome:
        out.print(f"  Candidate lineages: {len(grouped)}")

    # Category breakdown
    category_counts = novel_df.group_by("taxonomic_call").len().to_dict(as_series=False)
    cat_dict = dict(zip(category_counts["taxonomic_call"], category_counts["len"], strict=False))
    if "Novel Species" in cat_dict:
        out.print(f"  Novel Species: {cat_dict['Novel Species']:,}")
    if "Novel Genus" in cat_dict:
        out.print(f"  Novel Genus: {cat_dict['Novel Genus']:,}")

    out.print()


@app.command(name="sensitivity")
def sensitivity(
    classifications: Path = typer.Option(
        ...,
        "--classifications",
        "-c",
        help=(
            "Classification CSV/Parquet produced by 'score classify'. Must "
            "contain the columns novelty_index, placement_uncertainty, "
            "num_ambiguous_hits, and identity_gap."
        ),
        exists=True,
        dir_okay=False,
    ),
    output: Path = typer.Option(
        ...,
        "--output",
        "-o",
        help="Output file path (TSV or JSON based on --format).",
    ),
    config_file: Path | None = typer.Option(
        None,
        "--config",
        help=(
            "Base ScoringConfig YAML to vary. Defaults to built-in defaults. "
            "Threshold ranges sweep around the novelty_known_max and "
            "uncertainty_known_max fields."
        ),
        exists=True,
        dir_okay=False,
    ),
    novelty_min: float = typer.Option(
        2.0,
        "--novelty-min",
        help="Lower bound of the novelty_known_max sweep (default 2.0).",
        min=0.0,
        max=100.0,
    ),
    novelty_max: float = typer.Option(
        8.0,
        "--novelty-max",
        help="Upper bound of the novelty_known_max sweep (default 8.0).",
        min=0.0,
        max=100.0,
    ),
    uncertainty_min: float = typer.Option(
        0.5,
        "--uncertainty-min",
        help="Lower bound of the uncertainty_known_max sweep (default 0.5).",
        min=0.0,
        max=100.0,
    ),
    uncertainty_max: float = typer.Option(
        3.0,
        "--uncertainty-max",
        help="Upper bound of the uncertainty_known_max sweep (default 3.0).",
        min=0.0,
        max=100.0,
    ),
    steps: int = typer.Option(
        9,
        "--steps",
        help="Number of threshold points to evaluate (default 9).",
        min=2,
        max=101,
    ),
    output_format: str = typer.Option(
        "tsv",
        "--format",
        "-f",
        help="Output format: 'tsv' (long form) or 'json' (case-insensitive).",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet",
        "-q",
        help="Suppress progress output.",
    ),
) -> None:
    """Sweep classification thresholds and report category counts.

    Re-classifies the rows of an existing classification file across a
    grid of (novelty_known_max, uncertainty_known_max) values, varying
    both together along a diagonal (point i uses the i-th novelty AND the
    i-th uncertainty value; this is not a 2D grid) in ``steps`` increments.
    Counts are for the legacy hard-threshold cascade (column ``legacy_call``),
    not the Bayesian ``taxonomic_call``; all categories present are listed so
    each point sums to the total read count.

    Example:

        metadarkmatter score sensitivity \\
            --classifications results.csv \\
            --output sensitivity.tsv \\
            --novelty-min 2 --novelty-max 10 --steps 17
    """
    import json

    from metadarkmatter.core.classification.sensitivity import (
        run_sensitivity_analysis,
    )

    if novelty_min >= novelty_max:
        console.print("[red]--novelty-min must be strictly less than --novelty-max[/red]")
        raise typer.Exit(code=2)
    if uncertainty_min >= uncertainty_max:
        console.print(
            "[red]--uncertainty-min must be strictly less than --uncertainty-max[/red]"
        )
        raise typer.Exit(code=2)
    output_format = output_format.lower()
    if output_format not in ("tsv", "json"):
        console.print(
            f"[red]Invalid --format '{output_format}'. Use 'tsv' or 'json'.[/red]"
        )
        raise typer.Exit(code=2)

    out = QuietConsole(console, quiet=quiet)
    out.print("\n[bold blue]Threshold sensitivity analysis[/bold blue]\n")

    df = read_dataframe(classifications)
    required = {"novelty_index", "placement_uncertainty", "num_ambiguous_hits", "identity_gap"}
    missing = required - set(df.columns)
    if missing:
        console.print(
            f"[red]Classification file is missing required columns: "
            f"{sorted(missing)}.[/red]\n"
            "[dim]Run 'score classify' to produce a compatible file.[/dim]"
        )
        raise typer.Exit(code=1)

    base_config = (
        ScoringConfig.from_yaml(config_file) if config_file is not None else ScoringConfig()
    )

    out.print(
        f"Sweeping novelty in [{novelty_min}, {novelty_max}] "
        f"and uncertainty in [{uncertainty_min}, {uncertainty_max}] "
        f"over {steps} points on {len(df):,} reads ..."
    )

    result = run_sensitivity_analysis(
        df=df,
        base_config=base_config,
        novelty_range=(novelty_min, novelty_max),
        uncertainty_range=(uncertainty_min, uncertainty_max),
        steps=steps,
    )

    output.parent.mkdir(parents=True, exist_ok=True)

    if output_format == "json":
        output.write_text(json.dumps(result.to_dict(), indent=2))
    else:
        # Long-form TSV: one row per (threshold_point, category).
        rows: list[dict[str, Any]] = []
        for i, (n_thr, u_thr) in enumerate(
            zip(result.novelty_thresholds, result.uncertainty_thresholds, strict=True)
        ):
            for category, counts in result.counts.items():
                rows.append(
                    {
                        "novelty_known_max": round(n_thr, 4),
                        "uncertainty_known_max": round(u_thr, 4),
                        "legacy_call": category,
                        "count": counts[i],
                    }
                )
        pl.DataFrame(rows).write_csv(output, separator="\t")

    out.print(f"\n[green]Sensitivity results written to:[/green] {output}")


@app.command(name="evaluate")
def evaluate(
    predictions: Path = typer.Option(
        ...,
        "--predictions",
        "-p",
        help=(
            "Classification CSV/Parquet produced by 'score classify'. "
            "Must contain a 'taxonomic_call' column; 'confidence_score' "
            "is required for calibration metrics."
        ),
        exists=True,
        dir_okay=False,
    ),
    truth_column: str = typer.Option(
        "true_category",
        "--truth-column",
        help="Name of the ground-truth column in the predictions file.",
    ),
    output: Path | None = typer.Option(
        None,
        "--output",
        "-o",
        help=(
            "Optional JSON path for the full evaluation result; without it "
            "only the console summary is printed."
        ),
    ),
    n_confidence_bins: int = typer.Option(
        10,
        "--bins",
        help="Number of equal-width confidence bins for ECE (default 10).",
        min=2,
        max=100,
    ),
    quiet: bool = typer.Option(
        False, "--quiet", "-q", help="Suppress progress and summary output."
    ),
) -> None:
    """Score classifications against a labelled ground-truth column.

    Reports top-1 accuracy, per-category precision and recall, a
    confusion matrix, and (when 'confidence_score' is present) the
    expected calibration error and a reliability histogram. Use this
    after 'score classify' to close the calibration loop produced by
    'scripts/calibrate_bayesian.py' and 'scripts/calibrate_entropy.py'.
    """
    import json

    from metadarkmatter.core.classification.evaluation import (
        evaluate_classifications,
    )

    out = QuietConsole(console, quiet=quiet)
    df = read_dataframe(predictions)
    if truth_column not in df.columns:
        console.print(
            f"[red]Truth column '{truth_column}' missing from "
            f"{predictions}.[/red]"
        )
        raise typer.Exit(code=1)
    if "taxonomic_call" not in df.columns:
        console.print(
            "[red]Predictions file must contain a 'taxonomic_call' column.[/red]"
        )
        raise typer.Exit(code=1)

    has_conf = "confidence_score" in df.columns
    result = evaluate_classifications(
        df,
        truth_column=truth_column,
        prediction_column="taxonomic_call",
        confidence_column="confidence_score" if has_conf else None,
        n_confidence_bins=n_confidence_bins,
    )

    out.print(f"\n[bold]Evaluation[/bold]  ({result.n_rows:,} rows)")
    out.print(f"  Overall accuracy: [bold]{result.overall_accuracy * 100:.2f}%[/bold]")
    if result.expected_calibration_error is not None:
        out.print(
            f"  Expected calibration error: "
            f"[bold]{result.expected_calibration_error * 100:.2f}%[/bold]"
        )

    table = Table(title="Per-category", show_header=True, header_style="bold")
    table.add_column("Category")
    table.add_column("Support", justify="right")
    table.add_column("Predicted", justify="right")
    table.add_column("Precision", justify="right")
    table.add_column("Recall", justify="right")
    for row in result.per_category:
        table.add_row(
            row.category,
            str(row.support),
            str(row.predicted),
            f"{row.precision * 100:.1f}%",
            f"{row.recall * 100:.1f}%",
        )
    out.print(table)

    if output is not None:
        payload = {
            "n_rows": result.n_rows,
            "overall_accuracy": result.overall_accuracy,
            "expected_calibration_error": result.expected_calibration_error,
            "per_category": [
                {
                    "category": r.category,
                    "support": r.support,
                    "predicted": r.predicted,
                    "correct": r.correct,
                    "precision": r.precision,
                    "recall": r.recall,
                }
                for r in result.per_category
            ],
            "confusion": result.confusion,
            "confidence_bins": {
                "centers": result.confidence_bin_centers,
                "accuracy": result.confidence_bin_accuracy,
                "counts": result.confidence_bin_counts,
            },
        }
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(payload, indent=2))
        out.print(f"[green]JSON written to {output}[/green]")
