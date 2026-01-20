"""
Score command for ANI-weighted placement classification.

This is the core command that most users will interact with, as it implements
the ANI-weighted placement uncertainty algorithm on existing BLAST results.
"""

from __future__ import annotations

import logging
from enum import Enum
from pathlib import Path
from typing import Any

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
from metadarkmatter.core.aai_matrix_builder import AAIMatrix
from metadarkmatter.core.ani_placement import (
    ANIMatrix,
    ANIWeightedClassifier,
    VectorizedClassifier,
)
from metadarkmatter.core.id_mapping import ContigIdMapping
from metadarkmatter.core.io_utils import read_dataframe, write_dataframe
from metadarkmatter.core.metadata import GenomeMetadata
from metadarkmatter.core.parsers import StreamingBlastParser, extract_genome_name_expr
from metadarkmatter.models.classification import TaxonomicSummary
from metadarkmatter.models.config import ScoringConfig


# GTDB-compatible threshold presets
# Based on Parks et al. 2018, 2020 and GTDB methodology
THRESHOLD_PRESETS: dict[str, ScoringConfig] = {
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
        # 96% species boundary (middle ground)
        novelty_known_max=4.0,
        novelty_novel_species_min=4.0,
        novelty_novel_species_max=15.0,
        novelty_novel_genus_min=15.0,
        novelty_novel_genus_max=25.0,
        # Stricter placement uncertainty
        uncertainty_known_max=1.5,
        uncertainty_novel_species_max=1.5,
        uncertainty_novel_genus_max=1.5,
        # Standard alignment quality
        min_alignment_length=100,
        min_alignment_fraction=0.0,
    ),
    # Literature-strict preset based on comprehensive literature review
    # See docs/CLASSIFICATION_STATISTICS.md for references:
    # - Jain et al. 2018 (Nature Communications): 95-96% ANI species boundary
    # - Riesco & Trujillo 2024: Genus delineation considerations
    # - GTDB standards (Parks et al. 2020)
    "literature-strict": ScoringConfig(
        # 96% ANI species boundary (stricter than 95%)
        novelty_known_max=4.0,
        novelty_novel_species_min=4.0,
        # 85% identity = conservative novel species cutoff
        novelty_novel_species_max=15.0,
        novelty_novel_genus_min=15.0,
        # 78% identity = genus boundary (more conservative than default 75%)
        novelty_novel_genus_max=22.0,
        # Stricter placement confidence (98.5% ANI required)
        uncertainty_known_max=1.5,
        uncertainty_novel_species_max=1.5,
        uncertainty_novel_genus_max=1.5,
        # Higher confidence threshold
        confidence_threshold=60.0,
        # GTDB-compatible alignment requirements
        min_alignment_length=100,
        min_alignment_fraction=0.5,
    ),
    "default": ScoringConfig(),
    # Coverage-weighted presets
    # These prioritize longer alignments over short conserved domains
    "coverage-linear": ScoringConfig(
        coverage_weight_mode="linear",
        coverage_weight_strength=0.5,
    ),
    "coverage-strict": ScoringConfig(
        coverage_weight_mode="sigmoid",
        coverage_weight_strength=0.7,
        novelty_known_max=4.0,
        novelty_novel_species_min=4.0,
        novelty_novel_species_max=15.0,
        novelty_novel_genus_min=15.0,
        novelty_novel_genus_max=25.0,
        uncertainty_known_max=1.5,
        uncertainty_novel_species_max=1.5,
        uncertainty_novel_genus_max=1.5,
    ),
    "coverage-gentle": ScoringConfig(
        coverage_weight_mode="log",
        coverage_weight_strength=0.3,
    ),
    "gtdb-coverage": ScoringConfig(
        coverage_weight_mode="linear",
        coverage_weight_strength=0.5,
        min_alignment_length=100,
        min_alignment_fraction=0.5,
    ),
}


class ProcessingMode(str, Enum):
    """Processing mode for classification."""

    STANDARD = "standard"
    FAST = "fast"
    PARALLEL = "parallel"
    STREAMING = "streaming"


def validate_processing_modes(
    fast: bool,
    parallel: bool,
    streaming: bool,
) -> ProcessingMode:
    """
    Validate that only one processing mode is selected.

    Returns the selected mode or STANDARD if none specified.

    Raises:
        typer.BadParameter: If multiple modes are specified.
    """
    modes_selected = sum([fast, parallel, streaming])

    if modes_selected > 1:
        selected = []
        if fast:
            selected.append("--fast")
        if parallel:
            selected.append("--parallel")
        if streaming:
            selected.append("--streaming")

        raise typer.BadParameter(
            f"Processing modes are mutually exclusive. "
            f"You specified: {', '.join(selected)}. "
            f"Choose only one:\n\n"
            f"  --fast      : Single-threaded optimized (~3x faster)\n"
            f"  --parallel  : Polars vectorized (~16x faster, recommended)\n"
            f"  --streaming : Memory-efficient for 100M+ alignments"
        )

    if streaming:
        return ProcessingMode.STREAMING
    if parallel:
        return ProcessingMode.PARALLEL
    if fast:
        return ProcessingMode.FAST
    return ProcessingMode.STANDARD


def validate_ani_genome_coverage(
    blast_path: Path,
    ani_matrix: ANIMatrix,
    sample_rows: int = 10000,
) -> tuple[int, int, float, set[str]]:
    """
    Early validation of ANI matrix coverage against BLAST file genomes.

    Samples the BLAST file and checks how many unique genomes are present
    in the ANI matrix. This helps detect mismatched inputs early.

    Args:
        blast_path: Path to BLAST file
        ani_matrix: Loaded ANI matrix
        sample_rows: Number of rows to sample (default 10000)

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

    # Calculate coverage
    matched = blast_genomes & ani_genomes
    missing = blast_genomes - ani_genomes

    total_genomes = len(blast_genomes)
    matched_count = len(matched)
    coverage_pct = (100.0 * matched_count / total_genomes) if total_genomes > 0 else 0.0

    return matched_count, total_genomes, coverage_pct, missing


def validate_output_format_extension(
    output: Path,
    output_format: str,
    console: Console,
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


def _finalize_classification(
    classification_df: pl.DataFrame | None,
    genome_metadata: Any,
    output: Path,
    output_format: str,
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


app = typer.Typer(
    name="score",
    help="Classify BLAST results using ANI-weighted placement",
    no_args_is_help=True,
)

console = Console()


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
            "Use predefined threshold preset: 'gtdb-strict' (95% ANI, AF>=50%), "
            "'gtdb-relaxed' (97% ANI), 'conservative' (96% ANI), "
            "'literature-strict' (96% ANI, 1.5% uncertainty, high confidence), "
            "'coverage-linear' (linear coverage weighting), "
            "'coverage-strict' (sigmoid coverage weighting, stricter thresholds), "
            "'coverage-gentle' (log coverage weighting), "
            "'gtdb-coverage' (GTDB + linear coverage weighting), or 'default'"
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
        0.0,
        "--min-alignment-fraction",
        help=(
            "Minimum fraction of read aligned (like GTDB's AF). "
            "Set to 0.5 for GTDB-compatible filtering. Default 0.0 = no filter."
        ),
        min=0.0,
        max=1.0,
    ),
    coverage_weight_mode: str = typer.Option(
        "none",
        "--coverage-weight-mode",
        help=(
            "Coverage weighting mode for hit selection: 'none' (default, raw bitscore), "
            "'linear' (gradual weight increase), 'log' (diminishing returns), "
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
    fast: bool = typer.Option(
        False,
        "--fast",
        help="Use optimized fast path (~3x faster for large files)",
    ),
    parallel: bool = typer.Option(
        False,
        "--parallel",
        "-p",
        help="Use vectorized Polars processing (~16x faster, auto-parallel)",
    ),
    streaming: bool = typer.Option(
        False,
        "--streaming",
        help="Use streaming mode for very large files (100M+ alignments)",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Validate inputs and show what would be processed without running",
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

    # Validate output extension matches format
    output = validate_output_format_extension(output, output_format, out)

    # Validate processing mode exclusivity
    try:
        processing_mode = validate_processing_modes(fast, parallel, streaming)
    except typer.BadParameter as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(code=1) from None

    if verbose:
        out.print(f"[dim]Processing mode: {processing_mode.value}[/dim]")

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
        except Exception:
            console.print("  Alignments: (unable to count)")

        # Show ANI matrix info
        console.print(f"\n[bold]ANI Matrix:[/bold] {ani}")
        try:
            ani_matrix = ANIMatrix.from_file(ani)
            console.print(f"  Genomes: {len(ani_matrix.genomes)}")
        except Exception as e:
            console.print(f"  [red]Error loading: {e}[/red]")
            raise typer.Exit(code=1) from None

        # Show genome coverage
        try:
            matched, total, coverage_pct, missing = validate_ani_genome_coverage(
                alignment, ani_matrix
            )
            console.print("\n[bold]Genome Coverage:[/bold]")
            console.print(f"  Matched: {matched}/{total} ({coverage_pct:.1f}%)")
            if coverage_pct < 50.0:
                console.print("  [yellow]Warning: Low coverage[/yellow]")
        except Exception:
            logger.debug("Could not validate genome coverage", exc_info=True)

        # Show output configuration
        console.print("\n[bold]Output:[/bold]")
        console.print(f"  File: {output}")
        console.print(f"  Format: {output_format}")
        console.print(f"  Mode: {processing_mode.value}")
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
        out.print(f"[dim]Using preset: {preset_lower}[/dim]")
        # Override bitscore threshold and alignment mode if explicitly provided
        if bitscore_threshold != 95.0 or alignment_mode_lower != "nucleotide" or coverage_weight_mode != "none":
            config = ScoringConfig(
                alignment_mode=alignment_mode_lower,
                bitscore_threshold_pct=bitscore_threshold,
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
                coverage_weight_mode=coverage_weight_mode,
                coverage_weight_strength=coverage_weight_strength,
            )
    else:
        config = ScoringConfig(
            alignment_mode=alignment_mode_lower,
            bitscore_threshold_pct=bitscore_threshold,
            min_alignment_length=min_alignment_length,
            min_alignment_fraction=min_alignment_fraction,
            coverage_weight_mode=coverage_weight_mode,
            coverage_weight_strength=coverage_weight_strength,
        )

    # Log alignment filter settings if non-default
    if config.min_alignment_length > 0 or config.min_alignment_fraction > 0:
        out.print(
            f"[dim]Alignment filters: length >= {config.min_alignment_length}bp, "
            f"fraction >= {config.min_alignment_fraction:.0%}[/dim]"
        )

    # Load ANI matrix
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
    ) as progress:
        progress.add_task(description="Loading ANI matrix...", total=None)
        try:
            ani_matrix = ANIMatrix.from_file(ani)
            num_genomes = len(ani_matrix.genomes)
        except Exception as e:
            console.print(f"\n[red]Error loading ANI matrix: {e}[/red]")
            raise typer.Exit(code=1) from None

    out.print(f"[green]Loaded ANI matrix for {num_genomes} genomes[/green]\n")

    # Load AAI matrix if provided (for genus-level classification)
    aai_matrix: AAIMatrix | None = None
    if aai:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
        ) as progress:
            progress.add_task(description="Loading AAI matrix...", total=None)
            try:
                aai_matrix = AAIMatrix.from_file(aai)
                aai_genomes = len(aai_matrix.genomes)
            except Exception as e:
                console.print(f"\n[red]Error loading AAI matrix: {e}[/red]")
                raise typer.Exit(code=1) from None

        out.print(f"[green]Loaded AAI matrix for {aai_genomes} genomes[/green]")
        out.print(
            "[dim]AAI will be used for genus-level classification (20-25% novelty)[/dim]\n"
        )

    # Load genome metadata if provided
    genome_metadata: GenomeMetadata | None = None
    if metadata:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
        ) as progress:
            progress.add_task(description="Loading genome metadata...", total=None)
            try:
                genome_metadata = GenomeMetadata.from_file(metadata)
            except Exception as e:
                console.print(f"\n[red]Error loading metadata: {e}[/red]")
                raise typer.Exit(code=1) from None

        out.print(
            f"[green]Loaded metadata for {genome_metadata.genome_count} genomes "
            f"({genome_metadata.species_count} species)[/green]\n"
        )

    # Load or generate ID mapping for external BLAST results
    contig_mapping: ContigIdMapping | None = None
    if genomes and id_mapping:
        console.print(
            "[red]Error: Cannot specify both --genomes and --id-mapping. "
            "Use one or the other.[/red]"
        )
        raise typer.Exit(code=1) from None

    if genomes:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
        ) as progress:
            progress.add_task(description="Generating ID mapping from genomes...", total=None)
            try:
                contig_mapping = ContigIdMapping.from_genome_dir(genomes)
            except Exception as e:
                console.print(f"\n[red]Error generating ID mapping: {e}[/red]")
                raise typer.Exit(code=1) from None

        out.print(
            f"[green]Generated ID mapping for {len(contig_mapping):,} contigs[/green]\n"
        )

    elif id_mapping:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
        ) as progress:
            progress.add_task(description="Loading ID mapping...", total=None)
            try:
                contig_mapping = ContigIdMapping.from_tsv(id_mapping)
            except Exception as e:
                console.print(f"\n[red]Error loading ID mapping: {e}[/red]")
                raise typer.Exit(code=1) from None

        out.print(f"[green]Loaded ID mapping with {len(contig_mapping):,} entries[/green]\n")

    # Early validation: check genome coverage between alignment file and ANI matrix
    try:
        matched, total, coverage_pct, missing = validate_ani_genome_coverage(
            alignment, ani_matrix
        )

        if verbose:
            out.print(
                f"[dim]Genome coverage: {matched}/{total} genomes in ANI matrix "
                f"({coverage_pct:.1f}%)[/dim]"
            )

        if coverage_pct < 50.0 and total > 0:
            out.print(
                f"[yellow]Warning: Low genome coverage ({coverage_pct:.1f}%)[/yellow]\n"
                f"  Only {matched} of {total} genomes in BLAST file are in ANI matrix.\n"
                f"  This may indicate mismatched input files.\n"
            )
            if missing and len(missing) <= 5:
                out.print(f"  Missing genomes: {', '.join(sorted(missing))}\n")
            elif missing:
                sample = sorted(missing)[:5]
                out.print(
                    f"  Missing genomes (first 5): {', '.join(sample)}...\n"
                )
    except Exception as e:
        # Always log warnings - don't silently swallow validation failures
        logger.warning("Could not validate genome coverage: %s", e)
        if verbose:
            out.print(f"[dim]Warning: Could not validate genome coverage: {e}[/dim]")

    # Run classification - keep DataFrame in memory to avoid redundant I/O
    classification_df: pl.DataFrame | None = None
    num_classified = 0

    try:
        if streaming:
            # Use streaming mode for very large files (100M+ alignments)
            # Streaming writes directly to file - no in-memory DataFrame
            out.print("[cyan]Streaming mode: processing in 5M alignment partitions[/cyan]\n")

            if genome_metadata:
                out.print(
                    "[yellow]Note: --metadata is not supported with --streaming mode.[/yellow]\n"
                    "[yellow]Species columns will not be added to output.[/yellow]\n"
                )

            if contig_mapping:
                out.print(
                    "[yellow]Note: --genomes/--id-mapping is not supported with --streaming mode.[/yellow]\n"
                    "[yellow]ID transformation will not be applied. Use --parallel instead.[/yellow]\n"
                )

            vectorized = VectorizedClassifier(ani_matrix=ani_matrix, aai_matrix=aai_matrix, config=config)

            # Rich progress with ETA for streaming
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
                # Estimate total rows from file size (rough: ~100 bytes per alignment)
                file_size = alignment.stat().st_size
                estimated_rows = file_size // 100
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

                num_classified = vectorized.stream_to_file(
                    blast_path=alignment,
                    output_path=output,
                    output_format=output_format,
                    progress_callback=streaming_progress,
                )

        elif parallel:
            # Use vectorized classifier (Polars-native, auto-parallel)
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                console=console,
            ) as progress:
                progress.add_task(
                    description="Classifying reads (parallel, auto-parallel)...",
                    total=None,
                )
                vectorized = VectorizedClassifier(ani_matrix=ani_matrix, aai_matrix=aai_matrix, config=config)
                classification_df = vectorized.classify_file(alignment, id_mapping=contig_mapping)

            num_classified = _finalize_classification(
                classification_df, genome_metadata, output, output_format
            )

        elif fast:
            # Use fast single-threaded path
            if contig_mapping:
                out.print(
                    "[yellow]Note: --genomes/--id-mapping is not supported with --fast mode.[/yellow]\n"
                    "[yellow]ID transformation will not be applied. Use --parallel instead.[/yellow]\n"
                )

            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                console=console,
            ) as progress:
                progress.add_task(description="Classifying reads (fast mode)...", total=None)
                classifier = ANIWeightedClassifier(ani_matrix=ani_matrix, aai_matrix=aai_matrix, config=config)
                classification_df = classifier.classify_to_dataframe_fast(alignment)

            num_classified = _finalize_classification(
                classification_df, genome_metadata, output, output_format
            )

        else:
            # Use standard path
            if contig_mapping:
                out.print(
                    "[yellow]Note: --genomes/--id-mapping is not supported with standard mode.[/yellow]\n"
                    "[yellow]ID transformation will not be applied. Use --parallel instead.[/yellow]\n"
                )

            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                console=console,
            ) as progress:
                progress.add_task(description="Classifying reads...", total=None)
                classifier = ANIWeightedClassifier(ani_matrix=ani_matrix, aai_matrix=aai_matrix, config=config)
                classification_df = classifier.classify_to_dataframe(alignment)

            num_classified = _finalize_classification(
                classification_df, genome_metadata, output, output_format
            )

    except pl.exceptions.PolarsError as e:
        console.print(f"\n[red]Data processing error: {e}[/red]")
        console.print(
            "[dim]This may indicate malformed BLAST data or incompatible file format.[/dim]"
        )
        if verbose:
            console.print_exception()
        raise typer.Exit(code=1) from None
    except FileNotFoundError as e:
        console.print(f"\n[red]File not found: {e}[/red]")
        raise typer.Exit(code=1) from None
    except PermissionError as e:
        console.print(f"\n[red]Permission denied: {e}[/red]")
        console.print("[dim]Check file permissions and try again.[/dim]")
        raise typer.Exit(code=1) from None
    except MemoryError:
        console.print(
            "\n[red]Out of memory during classification.[/red]\n"
            "[dim]Try using --streaming mode for large files, or reduce chunk size.[/dim]"
        )
        raise typer.Exit(code=1) from None
    except Exception as e:
        console.print(f"\n[red]Unexpected error during classification: {e}[/red]")
        if verbose:
            console.print_exception()
        raise typer.Exit(code=1) from None

    out.print(f"[green]Classified {num_classified:,} reads[/green]\n")

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

            except Exception as e:
                console.print(f"\n[yellow]Warning: Failed to generate summary: {e}[/yellow]")
                if verbose:
                    console.print_exception()

    out.print(f"\n[bold green]Results written to:[/bold green] {output}")
    if summary:
        out.print(f"[bold green]Summary written to:[/bold green] {summary}")
    out.print()


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
            "'conservative', 'literature-strict', 'coverage-linear', "
            "'coverage-strict', 'coverage-gentle', 'gtdb-coverage', or 'default'"
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
        0.0,
        "--min-alignment-fraction",
        help="Minimum fraction of read aligned (GTDB uses 0.5)",
        min=0.0,
        max=1.0,
    ),
    coverage_weight_mode: str = typer.Option(
        "none",
        "--coverage-weight-mode",
        help=(
            "Coverage weighting mode: 'none' (default), 'linear', 'log', or 'sigmoid'. "
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
    fast: bool = typer.Option(
        False,
        "--fast",
        help="Use optimized fast path (~3x faster for large files)",
    ),
    parallel: bool = typer.Option(
        False,
        "--parallel",
        help="Use parallel processing per file (multi-core)",
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

    # Validate processing mode exclusivity (batch only supports fast and parallel)
    try:
        processing_mode = validate_processing_modes(fast, parallel, streaming=False)
    except typer.BadParameter as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(code=1) from None

    if verbose:
        console.print(f"[dim]Processing mode: {processing_mode.value}[/dim]")

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
    except Exception as e:
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
        except Exception as e:
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
        # Override with alignment mode if protein
        if alignment_mode_lower != "nucleotide":
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
            )
    else:
        config = ScoringConfig(
            alignment_mode=alignment_mode_lower,
            bitscore_threshold_pct=bitscore_threshold,
            min_alignment_length=min_alignment_length,
            min_alignment_fraction=min_alignment_fraction,
            coverage_weight_mode=coverage_weight_mode,
            coverage_weight_strength=coverage_weight_strength,
        )

    # Log alignment filter settings if non-default
    if config.min_alignment_length > 0 or config.min_alignment_fraction > 0:
        console.print(
            f"[dim]Alignment filters: length >= {config.min_alignment_length}bp, "
            f"fraction >= {config.min_alignment_fraction:.0%}[/dim]"
        )

    # Create appropriate classifier based on mode
    # Note: Batch mode does not currently support AAI matrix
    if parallel:
        vectorized = VectorizedClassifier(ani_matrix=ani_matrix, aai_matrix=None, config=config)
        mode_str = "vectorized (Polars-native parallel)"
    else:
        classifier = ANIWeightedClassifier(ani_matrix=ani_matrix, aai_matrix=None, config=config)
        mode_str = "fast" if fast else "standard"

    console.print(f"Processing mode: {mode_str}\n")

    # Process each file
    total_classified = 0
    failed_files = []

    for i, alignment_file in enumerate(alignment_files, 1):
        sample_name = extract_sample_name(alignment_file)
        output_file = output_dir / f"{sample_name}_classifications{file_ext}"
        summary_file = output_dir / f"{sample_name}_summary.json"

        console.print(f"[{i}/{len(alignment_files)}] Processing {sample_name}...")

        try:
            # Classify using the appropriate method
            if parallel:
                df = vectorized.classify_file(alignment_file)
            elif fast:
                df = classifier.classify_to_dataframe_fast(alignment_file)
            else:
                df = classifier.classify_to_dataframe(alignment_file)
            num_classified = len(df)

            # Write results in specified format
            if num_classified > 0:
                write_dataframe(df, output_file, output_format)

                # Generate summary using in-memory DataFrame
                summary_obj = _generate_summary(df)
                summary_obj.to_json(summary_file)

            console.print(f"  [green]Classified {num_classified:,} reads[/green]")
            total_classified += num_classified

        except Exception as e:
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
        ambiguous=count_dict.get("Ambiguous", 0),
        conserved_regions=count_dict.get("Conserved Region", 0),
        unclassified=count_dict.get("Unclassified", 0),
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
      - species: Novel species candidates (5-15% novelty, <2% uncertainty)
      - genus: Novel genus candidates (15-25% novelty, <2% uncertainty)
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
    except Exception as e:
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
        filter_desc = "Novel Species (5-15% novelty, <2% uncertainty)"
    elif category == "genus":
        filter_expr = pl.col("taxonomic_call") == "Novel Genus"
        filter_desc = "Novel Genus (15-25% novelty, <2% uncertainty)"
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
