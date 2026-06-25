"""
Score command for ANI-weighted placement classification.

This is the core command that most users will interact with, as it implements
the ANI-weighted placement uncertainty algorithm on existing BLAST results.
"""

from __future__ import annotations

import logging
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

from metadarkmatter.cli.utils import QuietConsole, extract_sample_name
from metadarkmatter.core.ani_placement import (
    ANIMatrix,
    VectorizedClassifier,
)
from metadarkmatter.core.classification.runner import (
    ClassificationRequest,
    run_classification,
    validate_ani_genome_coverage,
)
from metadarkmatter.core.exceptions import ConfigurationError
from metadarkmatter.core.io_utils import OutputFormat, read_dataframe, write_dataframe
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
        # Override bitscore/alignment-mode/weighting/uncertainty if explicitly provided
        if (
            opts.bitscore_threshold != 95.0
            or opts.alignment_mode != "nucleotide"
            or opts.coverage_weight_mode != "linear"
            or opts.uncertainty_mode != "second"
        ):
            config = ScoringConfig(
                alignment_mode=opts.alignment_mode,
                bitscore_threshold_pct=opts.bitscore_threshold,
                novelty_known_max=config.novelty_known_max,
                novelty_novel_species_min=config.novelty_novel_species_min,
                novelty_novel_species_max=config.novelty_novel_species_max,
                novelty_novel_genus_min=config.novelty_novel_genus_min,
                novelty_novel_genus_max=config.novelty_novel_genus_max,
                uncertainty_known_max=config.uncertainty_known_max,
                uncertainty_novel_species_max=config.uncertainty_novel_species_max,
                uncertainty_novel_genus_max=config.uncertainty_novel_genus_max,
                uncertainty_conserved_min=config.uncertainty_conserved_min,
                min_alignment_length=config.min_alignment_length,
                min_alignment_fraction=config.min_alignment_fraction,
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
            )
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
    )
    return config, messages


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
            "(only applies when --use-inferred-for-single-hits is enabled). "
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
            "Add Bayesian posterior probability columns to output. "
            "Provides continuous confidence (p_known_species, p_novel_species, "
            "p_novel_genus, p_ambiguous) rather than hard category labels."
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
    if dry_run:
        console.print("[bold cyan]DRY RUN MODE[/bold cyan] - Validating inputs only\n")

        # Show alignment file info
        alignment_size = alignment.stat().st_size
        console.print(f"[bold]Alignment File:[/bold] {alignment}")
        console.print(f"  Size: {alignment_size / 1024 / 1024:.1f} MB")

        # Count lines (estimate reads)
        try:
            with alignment.open("rb") as f:
                line_count = sum(1 for _ in f)
            console.print(f"  Alignments: ~{line_count:,}")
        except OSError:
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
                summary_obj.to_json(summary)

                # Display summary table (skip in quiet mode)
                if not quiet:
                    _display_summary_table(summary_obj)

            except (pl.exceptions.PolarsError, ValueError, OSError) as e:
                console.print(f"\n[yellow]Warning: Failed to generate summary: {e}[/yellow]")
                if verbose:
                    console.print_exception()

    # Write QC metrics if requested
    if qc_output is not None and qc_metrics is not None:
        import json

        qc_output.parent.mkdir(parents=True, exist_ok=True)
        qc_output.write_text(json.dumps(qc_metrics.to_dict(), indent=2))
        out.print(f"[bold green]QC metrics written to:[/bold green] {qc_output}")

        # Display warnings
        if qc_metrics.warnings and not quiet:
            out.print()
            for warning in qc_metrics.warnings:
                out.print(f"[yellow]QC Warning:[/yellow] {warning}")

    out.print(f"\n[bold green]Results written to:[/bold green] {output}")
    if summary:
        out.print(f"[bold green]Summary written to:[/bold green] {summary}")

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
    pattern: str = typer.Option(
        "*.blast.tsv.gz",
        "--pattern",
        "-p",
        help="Glob pattern for BLAST files",
    ),
    bitscore_threshold: float = typer.Option(
        95.0,
        "--bitscore-threshold",
        help="Percentage of top bitscore for ambiguous hits",
    ),
    preset: str | None = typer.Option(
        None,
        "--preset",
        help=(
            "Use predefined threshold preset: 'gtdb-strict', 'gtdb-relaxed', "
            "'conservative', 'literature-strict', or 'coverage-strict'"
        ),
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
        help="Number of worker processes (default: CPU count - 1)",
        min=0,
    ),
) -> None:
    """
    Batch classify multiple alignment files.

    Processes all alignment files (BLAST or MMseqs2) matching the pattern
    in the input directory and writes corresponding classification results.

    Example:

        metadarkmatter score batch \\
            --alignment-dir results/blast/ \\
            --ani ani_matrix.csv \\
            --output-dir results/classifications/ \\
            --pattern "*.tsv.gz"

        # Use Parquet for 10x smaller files:
        metadarkmatter score batch \\
            --alignment-dir results/blast/ \\
            --ani ani_matrix.csv \\
            --output-dir results/classifications/ \\
            --format parquet

        # Use GTDB-compatible settings:
        metadarkmatter score batch \\
            --alignment-dir results/blast/ \\
            --ani ani_matrix.csv \\
            --output-dir results/classifications/ \\
            --preset gtdb-strict
    """
    console.print("\n[bold blue]Metadarkmatter Batch Classification[/bold blue]\n")

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

    # Find all alignment files
    alignment_files = list(alignment_dir.glob(pattern))

    if not alignment_files:
        console.print(
            f"[red]Error: No alignment files found matching pattern: '{pattern}'[/red]\n"
            f"  Directory: {alignment_dir}\n"
            f"  Pattern: {pattern}\n\n"
            f"Suggestions:\n"
            f"  - Check if the directory contains alignment files\n"
            f"  - Try a different pattern (e.g., '*.tsv', '*.blast.tsv.gz')"
        )
        raise typer.Exit(code=1) from None

    console.print(f"Found {len(alignment_files)} alignment files to process\n")

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load ANI matrix once
    console.print("Loading ANI matrix...")
    try:
        ani_matrix = ANIMatrix.from_file(ani)
        num_genomes = len(ani_matrix.genomes)
    except (FileNotFoundError, ValueError, pl.exceptions.PolarsError, OSError) as e:
        console.print(f"[red]Error loading ANI matrix: {e}[/red]")
        raise typer.Exit(code=1) from None

    console.print(f"[green]Loaded ANI matrix for {num_genomes} genomes[/green]\n")

    # Early validation: check genome coverage using first BLAST file
    if alignment_files:
        try:
            matched, total, coverage_pct, missing = validate_ani_genome_coverage(
                alignment_files[0], ani_matrix
            )

            if verbose:
                console.print(
                    f"[dim]Genome coverage (sample): {matched}/{total} genomes in ANI matrix "
                    f"({coverage_pct:.1f}%)[/dim]"
                )

            if coverage_pct < 50.0 and total > 0:
                console.print(
                    f"[yellow]Warning: Low genome coverage ({coverage_pct:.1f}%)[/yellow]\n"
                    f"  Only {matched} of {total} genomes in BLAST files are in ANI matrix.\n"
                    f"  This may indicate mismatched input files.\n"
                )
                if missing and len(missing) <= 5:
                    console.print(f"  Missing genomes: {', '.join(sorted(missing))}\n")
                elif missing:
                    sample = sorted(missing)[:5]
                    console.print(
                        f"  Missing genomes (first 5): {', '.join(sample)}...\n"
                    )
        except (FileNotFoundError, ValueError, pl.exceptions.PolarsError, OSError) as e:
            # Always log warnings - don't silently swallow validation failures
            logger.warning("Could not validate genome coverage: %s", e)
            if verbose:
                console.print(f"[dim]Warning: Could not validate genome coverage: {e}[/dim]")

    # Validate alignment mode
    alignment_mode_lower = alignment_mode.lower()
    if alignment_mode_lower not in ("nucleotide", "protein"):
        console.print(
            f"[red]Error: Unknown alignment mode '{alignment_mode}'. "
            f"Use 'nucleotide' or 'protein'.[/red]"
        )
        raise typer.Exit(code=1) from None

    if alignment_mode_lower == "protein":
        console.print("[bold]Protein mode enabled[/bold] - using wider novelty thresholds")
        console.print(
            "[dim]  Known Species: N < 10%, Novel Species: 10-25%, Novel Genus: 25-40%[/dim]"
        )
        console.print(
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

    if uncertainty_mode_lower == "max":
        console.print(
            "[bold]Uncertainty mode: max[/bold] - using maximum ANI to any competing genome"
        )

    # Initialize configuration (preset or custom options)
    if preset:
        preset_lower = preset.lower()
        if preset_lower not in THRESHOLD_PRESETS:
            console.print(
                f"[red]Error: Unknown preset '{preset}'. "
                f"Available: {', '.join(THRESHOLD_PRESETS.keys())}[/red]"
            )
            raise typer.Exit(code=1) from None
        config = THRESHOLD_PRESETS[preset_lower]
        console.print(f"[dim]Using preset: {preset_lower}[/dim]")
        # Override with alignment mode if protein or non-default uncertainty mode
        if alignment_mode_lower != "nucleotide" or uncertainty_mode_lower != "second":
            config = ScoringConfig(
                alignment_mode=alignment_mode_lower,
                bitscore_threshold_pct=config.bitscore_threshold_pct,
                novelty_known_max=config.novelty_known_max,
                novelty_novel_species_min=config.novelty_novel_species_min,
                novelty_novel_species_max=config.novelty_novel_species_max,
                novelty_novel_genus_min=config.novelty_novel_genus_min,
                novelty_novel_genus_max=config.novelty_novel_genus_max,
                uncertainty_known_max=config.uncertainty_known_max,
                uncertainty_novel_species_max=config.uncertainty_novel_species_max,
                uncertainty_novel_genus_max=config.uncertainty_novel_genus_max,
                uncertainty_conserved_min=config.uncertainty_conserved_min,
                min_alignment_length=config.min_alignment_length,
                min_alignment_fraction=config.min_alignment_fraction,
                max_evalue=max_evalue,
                min_percent_identity=min_percent_identity,
                min_bitscore=min_bitscore,
                min_read_length=min_read_length,
                min_query_coverage=min_query_coverage,
                uncertainty_mode=uncertainty_mode_lower,
            )
    else:
        config = ScoringConfig(
            alignment_mode=alignment_mode_lower,
            bitscore_threshold_pct=bitscore_threshold,
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
        )

    # Log active filter settings
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
        console.print(f"[dim]Alignment filters: {', '.join(active_filters)}[/dim]")

    # Note: Batch mode does not currently support AAI matrix
    vectorized = VectorizedClassifier(ani_matrix=ani_matrix, aai_matrix=None, config=config)

    # Process each file
    total_classified = 0
    failed_files = []

    for i, alignment_file in enumerate(alignment_files, 1):
        sample_name = extract_sample_name(alignment_file)
        output_file = output_dir / f"{sample_name}_classifications{file_ext}"
        summary_file = output_dir / f"{sample_name}_summary.json"

        console.print(f"[{i}/{len(alignment_files)}] Processing {sample_name}...")

        try:
            df = vectorized.classify_file(alignment_file)
            assert isinstance(df, pl.DataFrame)
            num_classified = len(df)

            # Write results in specified format
            if num_classified > 0:
                write_dataframe(df, output_file, output_format)

                # Generate summary using in-memory DataFrame
                summary_obj = _generate_summary(df)
                summary_obj.to_json(summary_file)

            console.print(f"  [green]Classified {num_classified:,} reads[/green]")
            total_classified += num_classified

        except (pl.exceptions.PolarsError, FileNotFoundError, ValueError, OSError) as e:
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

    return TaxonomicSummary(
        total_reads=total_reads,
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
        help="Maximum placement uncertainty (exclude ambiguous placements)",
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
    """
    out = QuietConsole(console, quiet=quiet)

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

    # Build filter expression based on category
    out.print(f"\n[bold]Step 2:[/bold] Filtering for {category} candidates...")

    if min_novelty is not None:
        # Custom novelty threshold
        filter_expr = (
            (pl.col("novelty_index") >= min_novelty) &
            (pl.col("placement_uncertainty") < max_uncertainty)
        )
        filter_desc = f"novelty >= {min_novelty}%, uncertainty < {max_uncertainty}%"
    elif category == "species":
        filter_expr = pl.col("taxonomic_call") == "Novel Species"
        filter_desc = "Novel Species (5-20% novelty, <2% uncertainty)"
    elif category == "genus":
        filter_expr = pl.col("taxonomic_call") == "Novel Genus"
        filter_desc = "Novel Genus (20-25% novelty, <2% uncertainty)"
    else:  # all
        filter_expr = (
            (pl.col("taxonomic_call") == "Novel Species") |
            (pl.col("taxonomic_call") == "Novel Genus")
        )
        filter_desc = "Novel Species + Novel Genus"

    # Apply filter
    novel_df = df.filter(filter_expr)

    if len(novel_df) == 0:
        console.print(f"[yellow]No reads found matching criteria: {filter_desc}[/yellow]")
        raise typer.Exit(code=0)

    out.print(f"  [green]Found {len(novel_df):,} novel candidates[/green]")
    out.print(f"  [dim]Filter: {filter_desc}[/dim]")

    # Group by genome if requested
    if group_by_genome:
        out.print("\n[bold]Step 3:[/bold] Grouping by reference genome...")

        grouped = novel_df.group_by("best_match_genome").agg([
            pl.len().alias("read_count"),
            pl.col("novelty_index").mean().alias("mean_novelty"),
            pl.col("novelty_index").min().alias("min_novelty"),
            pl.col("novelty_index").max().alias("max_novelty"),
            pl.col("top_hit_identity").mean().alias("mean_identity"),
            pl.col("placement_uncertainty").mean().alias("mean_uncertainty"),
            pl.col("taxonomic_call").mode().first().alias("dominant_category"),
        ]).filter(
            pl.col("read_count") >= min_reads
        ).sort("read_count", descending=True)

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
        novel_df.write_csv(output)
        out.print(f"[green]Novel reads written to:[/green] {output}")

    # Write read IDs for seqtk extraction
    if read_ids:
        out.print("\n[bold]Step 4:[/bold] Exporting read IDs...")

        novel_df.select("read_id").write_csv(read_ids, include_header=False)
        out.print(f"[green]Read IDs written to:[/green] {read_ids}")
        out.print(f"  [dim]Use with: seqtk subseq reads.fastq.gz {read_ids} > novel.fastq[/dim]")

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
        help="Output format: 'tsv' (long form) or 'json'.",
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
    both proportionally in ``steps`` increments. The output reveals
    whether category counts are stable or sharply threshold-dependent.

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
                        "taxonomic_call": category,
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
