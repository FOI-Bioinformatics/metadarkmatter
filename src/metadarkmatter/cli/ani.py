"""
ANI command for computing and validating ANI matrices.

Provides subcommands:
- compute: Build ANI matrix from genome directory using fastANI or skani
- validate: Check ANI matrix coverage against BLAST results
"""

from __future__ import annotations

from enum import Enum
from pathlib import Path

import typer
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, TimeElapsedColumn
from rich.table import Table

from metadarkmatter.cli.utils import QuietConsole
from metadarkmatter.core.ani_matrix_builder import (
    ani_dict_to_csv,
    parse_fastani_output,
    parse_skani_output,
    validate_ani_coverage,
)
from metadarkmatter.core.parsers import StreamingBlastParser
from metadarkmatter.external import ToolExecutionError


class ANIBackend(str, Enum):
    """ANI computation backend."""

    FASTANI = "fastani"
    SKANI = "skani"
    AUTO = "auto"


app = typer.Typer(
    name="ani",
    help="Compute and validate ANI matrices for placement uncertainty analysis",
    no_args_is_help=True,
)

console = Console()


def _detect_backend() -> ANIBackend:
    """Auto-detect available ANI backend.

    Preference order: skani (faster) > fastANI

    Returns:
        Detected backend.

    Raises:
        typer.Exit: If no backend is available.
    """
    from metadarkmatter.external import FastANI, Skani

    if Skani.check_available():
        return ANIBackend.SKANI

    if FastANI.check_available():
        return ANIBackend.FASTANI

    console.print("[red]Error: No ANI backend found[/red]")
    console.print("\nInstall one of:")
    console.print("  conda install -c bioconda skani")
    console.print("  conda install -c bioconda fastani")
    raise typer.Exit(code=1) from None


def _find_genome_files(
    genome_dir: Path,
    pattern: str,
) -> list[Path]:
    """Find genome files in directory, trying multiple patterns if needed."""
    genome_files = sorted(genome_dir.glob(pattern))

    if not genome_files:
        # Try alternative patterns
        for alt_pattern in ["*.fna", "*.fa", "*.fasta", "*.fna.gz", "*.fa.gz"]:
            if alt_pattern != pattern:
                genome_files = sorted(genome_dir.glob(alt_pattern))
                if genome_files:
                    break

    return genome_files


class ANISensitivity(str, Enum):
    """ANI computation sensitivity preset."""

    DEFAULT = "default"
    SENSITIVE = "sensitive"
    CUSTOM = "custom"


@app.command(name="compute")
def compute(
    genomes: Path = typer.Option(
        ...,
        "--genomes",
        "-g",
        help="Directory containing genome FASTA files",
        exists=True,
        file_okay=False,
    ),
    output: Path = typer.Option(
        ...,
        "--output",
        "-o",
        help="Output CSV file for ANI matrix",
    ),
    metadata_path: Path | None = typer.Option(
        None,
        "--metadata",
        "-m",
        help="Path to genome_metadata.tsv (required for --representatives-only)",
        exists=True,
        dir_okay=False,
    ),
    representatives_only: bool = typer.Option(
        False,
        "--representatives-only",
        help=(
            "Compute ANI only for species representative genomes (from metadata). "
            "Reduces O(N^2) comparisons for large families."
        ),
    ),
    backend: ANIBackend = typer.Option(
        ANIBackend.AUTO,
        "--backend",
        "-b",
        help="ANI computation backend: fastani, skani, or auto",
    ),
    sensitivity: ANISensitivity = typer.Option(
        ANISensitivity.SENSITIVE,
        "--sensitivity",
        "-s",
        help="Sensitivity preset: default (standard), sensitive (captures distant genomes)",
    ),
    min_af: float = typer.Option(
        None,
        "--min-af",
        help="Minimum aligned fraction (0-100). Lower = more distant comparisons. "
        "Default: 15 (default mode), 5 (sensitive mode)",
        min=0.0,
        max=100.0,
    ),
    screen_threshold: float = typer.Option(
        None,
        "--screen",
        help="Pre-filter threshold for ANI (0-100). Lower = attempt more comparisons. "
        "Default: 80 (default mode), 50 (sensitive mode). Only applies to skani.",
        min=0.0,
        max=100.0,
    ),
    threads: int = typer.Option(
        4,
        "--threads",
        "-t",
        help="Number of threads for parallel execution",
        min=1,
    ),
    genome_pattern: str = typer.Option(
        "*.fna",
        "--genome-pattern",
        "-p",
        help="Glob pattern for genome files",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Enable verbose output",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet",
        "-q",
        help="Suppress progress output",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Show commands without executing",
    ),
) -> None:
    """
    Compute ANI matrix from reference genomes.

    Builds an all-vs-all ANI matrix using either fastANI or skani,
    outputting in metadarkmatter CSV format for use with 'score classify'.

    Example:

        # Auto-detect backend
        metadarkmatter ani compute \\
            --genomes reference_genomes/ \\
            --output ani_matrix.csv \\
            --threads 16

        # Specify backend
        metadarkmatter ani compute \\
            --genomes reference_genomes/ \\
            --output ani_matrix.csv \\
            --backend skani

        # Sensitive mode for family-level comparisons
        metadarkmatter ani compute \\
            --genomes reference_genomes/ \\
            --output ani_matrix.csv \\
            --sensitivity sensitive

        # Custom settings for maximum sensitivity
        metadarkmatter ani compute \\
            --genomes reference_genomes/ \\
            --output ani_matrix.csv \\
            --min-af 5 --screen 50
    """
    import tempfile

    out = QuietConsole(console, quiet=quiet)

    # Determine sensitivity settings
    # Sensitive mode uses lower thresholds to capture more distant genome pairs
    if sensitivity == ANISensitivity.SENSITIVE:
        effective_min_af = min_af if min_af is not None else 5.0
        effective_screen = screen_threshold if screen_threshold is not None else 50.0
    else:  # DEFAULT or CUSTOM with explicit values
        effective_min_af = min_af if min_af is not None else 15.0
        effective_screen = screen_threshold if screen_threshold is not None else 80.0

    out.print("\n[bold blue]Metadarkmatter ANI Matrix Builder[/bold blue]\n")

    # Auto-detect backend if needed
    selected_backend = backend
    if backend == ANIBackend.AUTO:
        selected_backend = _detect_backend()
        out.print(f"[dim]Auto-detected backend: {selected_backend.value}[/dim]")

    # Check backend availability
    if selected_backend == ANIBackend.SKANI:
        from metadarkmatter.external import Skani

        if not Skani.check_available():
            console.print("[red]Error: skani not found in PATH[/red]")
            console.print("[dim]Install with: conda install -c bioconda skani[/dim]")
            raise typer.Exit(code=1) from None
    else:
        from metadarkmatter.external import FastANI

        if not FastANI.check_available():
            console.print("[red]Error: fastANI not found in PATH[/red]")
            console.print("[dim]Install with: conda install -c bioconda fastani[/dim]")
            raise typer.Exit(code=1) from None

    # Validate --representatives-only requires --metadata
    if representatives_only and metadata_path is None:
        console.print(
            "[red]Error: --representatives-only requires --metadata[/red]"
        )
        raise typer.Exit(code=1) from None

    # Find genome files
    genome_files = _find_genome_files(genomes, genome_pattern)

    if not genome_files:
        console.print(
            f"[red]Error: No genome files found matching '{genome_pattern}' "
            f"in {genomes}[/red]"
        )
        raise typer.Exit(code=1) from None

    # Filter to representative genomes only if requested
    total_genome_count = len(genome_files)
    if representatives_only and metadata_path is not None:
        from metadarkmatter.core.genome_utils import extract_accession_from_filename
        from metadarkmatter.core.metadata import GenomeMetadata

        gm = GenomeMetadata.from_file(metadata_path)
        rep_mapping = gm.build_representative_mapping()
        rep_accessions = set(rep_mapping.values())

        genome_files = [
            f for f in genome_files
            if extract_accession_from_filename(f.name) in rep_accessions
        ]

        if not genome_files:
            console.print(
                "[red]Error: No representative genome files found in directory[/red]"
            )
            console.print(
                f"[dim]Representatives: {len(rep_accessions)}, "
                f"genome files: {total_genome_count}[/dim]"
            )
            raise typer.Exit(code=1) from None

        out.print(
            f"[bold]Input:[/bold] {len(genome_files)} representative genomes "
            f"(of {total_genome_count} total) from {genomes}"
        )
    else:
        out.print(f"[bold]Input:[/bold] {len(genome_files)} genomes from {genomes}")
    out.print(f"[bold]Backend:[/bold] {selected_backend.value}")
    out.print(f"[bold]Threads:[/bold] {threads}")
    if selected_backend == ANIBackend.SKANI:
        out.print(f"[bold]Sensitivity:[/bold] {sensitivity.value} (min_af={effective_min_af}%, screen={effective_screen}%)")

    # Create output directory
    output.parent.mkdir(parents=True, exist_ok=True)

    if dry_run:
        out.print("\n[cyan]DRY RUN - Commands that would be executed:[/cyan]")
        if selected_backend == ANIBackend.SKANI:
            out.print(
                f"  skani triangle {genomes}/* -o <temp> --full-matrix -t {threads} "
                f"--min-af {effective_min_af} -s {effective_screen}"
            )
        else:
            out.print(
                f"  fastANI --ql <list> --rl <list> -o <temp> -t {threads} --matrix"
            )
        out.print(f"\nOutput would be written to: {output}")
        raise typer.Exit(code=0)

    # Run ANI computation
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            TimeElapsedColumn(),
            console=console if not quiet else None,
        ) as progress:
            if selected_backend == ANIBackend.SKANI:
                progress.add_task(
                    description=f"Computing ANI with skani ({len(genome_files)} genomes)...",
                    total=None,
                )

                from metadarkmatter.external import Skani

                skani = Skani()
                raw_output = tmp_path / "skani_output.txt"

                try:
                    result = skani.run_or_raise(
                        genomes=genome_files,
                        output=raw_output,
                        threads=threads,
                        full_matrix=True,
                        min_af=effective_min_af,
                        screen_threshold=effective_screen,
                    )
                except ToolExecutionError as e:
                    console.print(f"\n[red]skani failed:[/red]\n{e.message}")
                    raise typer.Exit(code=1) from None

                if verbose:
                    out.print(f"\n[dim]Command: {result.command_string}[/dim]")
                    if result.stderr.strip():
                        out.print(f"[dim]Tool output: {result.stderr.strip()[:300]}[/dim]")

                # Parse output
                ani_dict = parse_skani_output(raw_output, full_matrix=True)

            else:  # FastANI
                progress.add_task(
                    description=f"Computing ANI with fastANI ({len(genome_files)} genomes)...",
                    total=None,
                )

                from metadarkmatter.external import FastANI, create_genome_list_file

                # Create genome list file
                list_file = tmp_path / "genome_list.txt"
                create_genome_list_file(
                    genomes,
                    list_file,
                    patterns=(genome_pattern, "*.fna", "*.fa", "*.fasta"),
                )

                raw_output = tmp_path / "fastani_output.txt"

                fastani = FastANI()
                try:
                    result = fastani.run_or_raise(
                        query_list=list_file,
                        reference_list=list_file,
                        output=raw_output,
                        threads=threads,
                        matrix=True,
                    )
                except ToolExecutionError as e:
                    console.print(f"\n[red]fastANI failed:[/red]\n{e.message}")
                    raise typer.Exit(code=1) from None

                if verbose:
                    out.print(f"\n[dim]Command: {result.command_string}[/dim]")
                    if result.stderr.strip():
                        out.print(f"[dim]Tool output: {result.stderr.strip()[:300]}[/dim]")

                # Parse output
                ani_dict = parse_fastani_output(raw_output)

        # Convert to CSV (auto-detect compression from filename)
        out.print("\n[bold]Converting to metadarkmatter format...[/bold]")
        compress_output = str(output).endswith(".gz")
        n_genomes = ani_dict_to_csv(ani_dict, output, compress=compress_output)

    # Summary
    out.print("\n[bold green]ANI matrix computed successfully![/bold green]")
    out.print(f"\n[bold]Output:[/bold] {output}")
    out.print(f"[bold]Genomes:[/bold] {n_genomes}")

    file_size = output.stat().st_size / (1024 * 1024)
    out.print(f"[bold]Size:[/bold] {file_size:.2f} MB")

    if verbose:
        # Show sample values
        import polars as pl

        df = pl.read_csv(output)
        out.print(f"\n[dim]Matrix shape: {df.shape[0]} x {df.shape[1] - 1}[/dim]")

    out.print(f"\n[dim]Use with: metadarkmatter score classify --ani {output}[/dim]")
    out.print()


@app.command(name="validate")
def validate(
    ani: Path = typer.Option(
        ...,
        "--ani",
        "-a",
        help="ANI matrix CSV file to validate",
        exists=True,
        dir_okay=False,
    ),
    blast: Path = typer.Option(
        ...,
        "--blast",
        "-b",
        help="BLAST results file to check coverage against",
        exists=True,
        dir_okay=False,
    ),
    metadata_path: Path | None = typer.Option(
        None,
        "--metadata",
        "-m",
        help=(
            "Path to genome_metadata.tsv with representative column. "
            "When provided, validates that representative genomes are in ANI matrix "
            "rather than requiring all BLAST genomes."
        ),
        exists=True,
        dir_okay=False,
    ),
    sample_rows: int = typer.Option(
        100000,
        "--sample-rows",
        help="Number of BLAST rows to sample for coverage check",
        min=1000,
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="List missing genomes",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet",
        "-q",
        help="Suppress progress output",
    ),
) -> None:
    """
    Validate ANI matrix coverage against BLAST results.

    Checks that genomes appearing in BLAST results have corresponding
    entries in the ANI matrix. Reports missing genomes and coverage statistics.

    Example:

        metadarkmatter ani validate \\
            --ani ani_matrix.csv \\
            --blast sample.blast.tsv.gz
    """
    import polars as pl

    from metadarkmatter.core.ani_matrix_builder import extract_genome_name_from_path

    out = QuietConsole(console, quiet=quiet)

    out.print("\n[bold blue]Metadarkmatter ANI Matrix Validator[/bold blue]\n")

    # Load ANI matrix
    out.print(f"[bold]Loading ANI matrix:[/bold] {ani}")
    try:
        from metadarkmatter.core.ani_placement import ANIMatrix

        ani_matrix = ANIMatrix.from_file(ani)
        ani_genomes = ani_matrix.genomes
        out.print(f"  [green]Loaded {len(ani_genomes)} genomes[/green]")
    except Exception as e:
        console.print(f"[red]Error loading ANI matrix: {e}[/red]")
        raise typer.Exit(code=1) from None

    # Sample BLAST file for genomes
    out.print(f"\n[bold]Sampling BLAST file:[/bold] {blast}")
    out.print(f"  [dim]Sampling {sample_rows:,} rows...[/dim]")

    try:
        # Read BLAST output and extract genome names from subject IDs
        # Auto-detect column count for backward compatibility (12 or 13 columns)
        parser = StreamingBlastParser(blast)
        sample_df = pl.scan_csv(
            blast,
            separator="\t",
            has_header=False,
            new_columns=parser.column_names,
            n_rows=sample_rows,
        ).collect()

        # Extract genome names from subject sequence IDs
        blast_genomes: set[str] = set()
        for sseqid in sample_df["sseqid"].unique().to_list():
            genome_name = extract_genome_name_from_path(str(sseqid))
            if genome_name and genome_name != "unknown":
                blast_genomes.add(genome_name)

        out.print(f"  [green]Found {len(blast_genomes)} unique genomes in sample[/green]")

    except Exception as e:
        console.print(f"[red]Error reading BLAST file: {e}[/red]")
        raise typer.Exit(code=1) from None

    # Representative-aware validation when metadata is provided
    rep_mapping: dict[str, str] | None = None
    if metadata_path is not None:
        from metadarkmatter.core.metadata import GenomeMetadata as GM

        gm = GM.from_file(metadata_path)
        rep_mapping = gm.build_representative_mapping()
        if gm.has_representatives:
            rep_accessions = set(rep_mapping.values())
            out.print(
                f"\n[bold]Representative mode:[/bold] {len(rep_accessions)} "
                f"representative genomes for {gm.genome_count} total genomes"
            )

            # Check that all representatives are in ANI matrix
            matched_reps = rep_accessions & ani_genomes
            missing_reps = rep_accessions - ani_genomes

            out.print(f"  Representatives in ANI matrix: {len(matched_reps)}/{len(rep_accessions)}")
            if missing_reps:
                out.print(f"  [yellow]Missing representatives: {len(missing_reps)}[/yellow]")
                if verbose:
                    for rep in sorted(missing_reps)[:20]:
                        out.print(f"    - {rep}")

            # Map BLAST genomes to representatives for coverage check
            blast_reps = {rep_mapping.get(g, g) for g in blast_genomes}
            matched_count, total, coverage_pct, missing = validate_ani_coverage(
                ani_genomes, blast_reps
            )
            extra = ani_genomes - blast_reps

            out.print("\n[bold]Coverage Analysis (representative-mapped):[/bold]")

            table = Table(show_header=True, header_style="bold")
            table.add_column("Metric", style="cyan")
            table.add_column("Value", justify="right")

            table.add_row("BLAST genomes (sample)", f"{len(blast_genomes):,}")
            table.add_row("Mapped to representatives", f"{total:,}")
            table.add_row("ANI matrix genomes", f"{len(ani_genomes):,}")
            table.add_row("Matched", f"{matched_count:,}")
            table.add_row("Missing from ANI", f"{len(missing):,}")
            table.add_row("Extra in ANI", f"{len(extra):,}")
            table.add_row("Coverage", f"{coverage_pct:.1f}%")

            console.print(table)

            # Warnings/recommendations
            if coverage_pct < 50.0:
                console.print(f"\n[yellow]Warning: Low representative coverage ({coverage_pct:.1f}%)[/yellow]")
            elif coverage_pct < 90.0:
                console.print(f"\n[yellow]Note: Moderate coverage ({coverage_pct:.1f}%)[/yellow]")
            else:
                console.print(f"\n[green]Good representative coverage ({coverage_pct:.1f}%)[/green]")

            # Show missing
            if missing and verbose:
                console.print("\n[bold]Missing representative genomes (not in ANI matrix):[/bold]")
                for genome in sorted(missing)[:20]:
                    console.print(f"  - {genome}")
            elif missing:
                console.print(
                    f"\n[dim]Use --verbose to see list of {len(missing)} missing representatives[/dim]"
                )

            if coverage_pct < 50.0:
                raise typer.Exit(code=1) from None

            out.print()
            raise typer.Exit(code=0)

    # Standard (non-representative) coverage check
    # Calculate coverage
    matched_count, total, coverage_pct, missing = validate_ani_coverage(
        ani_genomes, blast_genomes
    )
    extra = ani_genomes - blast_genomes

    # Display results
    out.print("\n[bold]Coverage Analysis:[/bold]")

    table = Table(show_header=True, header_style="bold")
    table.add_column("Metric", style="cyan")
    table.add_column("Value", justify="right")

    table.add_row("BLAST genomes (sample)", f"{total:,}")
    table.add_row("ANI matrix genomes", f"{len(ani_genomes):,}")
    table.add_row("Matched", f"{matched_count:,}")
    table.add_row("Missing from ANI", f"{len(missing):,}")
    table.add_row("Extra in ANI", f"{len(extra):,}")
    table.add_row("Coverage", f"{coverage_pct:.1f}%")

    console.print(table)

    # Warnings/recommendations
    if coverage_pct < 50.0:
        console.print(f"\n[yellow]Warning: Low coverage ({coverage_pct:.1f}%)[/yellow]")
        console.print(
            "  This may indicate mismatched input files or incomplete ANI matrix."
        )
        console.print(
            "  Consider regenerating the ANI matrix with all reference genomes."
        )
    elif coverage_pct < 90.0:
        console.print(f"\n[yellow]Note: Moderate coverage ({coverage_pct:.1f}%)[/yellow]")
        console.print("  Some genomes in BLAST results are missing from ANI matrix.")
    else:
        console.print(f"\n[green]Good coverage ({coverage_pct:.1f}%)[/green]")

    # Show missing genomes
    if missing and verbose:
        console.print("\n[bold]Missing genomes (not in ANI matrix):[/bold]")
        for genome in sorted(missing)[:20]:
            console.print(f"  - {genome}")
        if len(missing) > 20:
            console.print(f"  ... and {len(missing) - 20} more")
    elif missing:
        console.print(
            f"\n[dim]Use --verbose to see list of {len(missing)} missing genomes[/dim]"
        )

    # Exit code based on coverage
    if coverage_pct < 50.0:
        raise typer.Exit(code=1) from None

    out.print()
