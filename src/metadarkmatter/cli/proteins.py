"""
Protein prediction commands using Prodigal.

Provides commands to predict protein-coding genes from genome FASTA files,
either for all genomes or only those missing protein files.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import typer
from rich.console import Console
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.table import Table

from metadarkmatter.cli.utils import QuietConsole

app = typer.Typer(
    name="proteins",
    help="Predict protein sequences from genome files using Prodigal",
    no_args_is_help=True,
)

console = Console()


def _find_genome_files(genome_dir: Path, pattern: str) -> list[Path]:
    """Find genome FASTA files in directory."""
    genome_files = sorted(genome_dir.glob(pattern))

    if not genome_files:
        # Try alternative patterns
        for alt_pattern in ["*.fna", "*.fa", "*.fasta", "*.fna.gz"]:
            if alt_pattern != pattern:
                genome_files = sorted(genome_dir.glob(alt_pattern))
                if genome_files:
                    break

    return genome_files


def _find_missing_proteins(
    genome_dir: Path,
    output_dir: Path,
    genome_pattern: str,
) -> tuple[list[Path], list[Path]]:
    """Find genomes that are missing protein files.

    Returns:
        Tuple of (genomes_missing_proteins, genomes_with_proteins)
    """
    genome_files = _find_genome_files(genome_dir, genome_pattern)

    missing = []
    existing = []

    for genome_file in genome_files:
        accession = genome_file.stem
        # Handle .fna.gz -> .faa
        if accession.endswith(".fna"):
            accession = accession[:-4]

        protein_file = output_dir / f"{accession}.faa"
        protein_file_gz = output_dir / f"{accession}.faa.gz"

        if protein_file.exists() or protein_file_gz.exists():
            existing.append(genome_file)
        else:
            missing.append(genome_file)

    return missing, existing


def _predict_single_genome(
    genome_file: Path,
    output_dir: Path,
    procedure: str,
) -> tuple[Path, int, str | None]:
    """Predict proteins for a single genome.

    Returns:
        Tuple of (genome_file, protein_count, error_message)
    """
    from metadarkmatter.external.prodigal import Prodigal

    accession = genome_file.stem
    if accession.endswith(".fna"):
        accession = accession[:-4]

    output_file = output_dir / f"{accession}.faa"

    try:
        prodigal = Prodigal()

        # Handle gzipped input
        if genome_file.suffix == ".gz":
            import gzip
            import tempfile

            with tempfile.NamedTemporaryFile(
                mode="wb", suffix=".fna", delete=False
            ) as tmp:
                tmp_path = Path(tmp.name)
                with gzip.open(genome_file, "rb") as gz_in:
                    tmp.write(gz_in.read())

            try:
                result = prodigal.predict_with_prefix(
                    input_file=tmp_path,
                    output_proteins=output_file,
                    genome_accession=accession,
                    procedure=procedure,
                    quiet=True,
                )
            finally:
                tmp_path.unlink()
        else:
            result = prodigal.predict_with_prefix(
                input_file=genome_file,
                output_proteins=output_file,
                genome_accession=accession,
                procedure=procedure,
                quiet=True,
            )

        if not result.success:
            return (genome_file, 0, result.stderr or "Unknown error")

        # Count proteins in output
        protein_count = 0
        if output_file.exists():
            with output_file.open("r") as f:
                protein_count = sum(1 for line in f if line.startswith(">"))

        return (genome_file, protein_count, None)

    except Exception as e:
        return (genome_file, 0, str(e))


@app.command(name="predict")
def predict(
    genomes: Path = typer.Option(
        ...,
        "--genomes",
        "-g",
        help="Directory containing genome FASTA files (.fna)",
        exists=True,
        file_okay=False,
    ),
    output_dir: Path = typer.Option(
        None,
        "--output-dir",
        "-o",
        help="Output directory for protein files (default: same as genomes)",
    ),
    missing_only: bool = typer.Option(
        False,
        "--missing-only",
        "-m",
        help="Only predict proteins for genomes without existing .faa files",
    ),
    threads: int = typer.Option(
        4,
        "--threads",
        "-t",
        help="Number of parallel Prodigal processes",
        min=1,
    ),
    procedure: str = typer.Option(
        "single",
        "--procedure",
        "-p",
        help="Prodigal procedure: 'single' for isolate genomes, 'meta' for metagenomes",
    ),
    genome_pattern: str = typer.Option(
        "*.fna",
        "--genome-pattern",
        help="Glob pattern for genome files",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Show detailed progress for each genome",
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
        help="Show what would be predicted without running",
    ),
) -> None:
    """
    Predict protein sequences from genome FASTA files using Prodigal.

    This command runs Prodigal gene prediction on genome files to generate
    protein FASTA files (.faa) suitable for AAI computation.

    Use --missing-only to fill in gaps from NCBI download (recommended):

        # After downloading genomes with --include-protein
        metadarkmatter proteins predict \\
            --genomes genomes/ \\
            --missing-only \\
            --threads 8

    Or predict proteins for all genomes:

        metadarkmatter proteins predict \\
            --genomes genomes/ \\
            --threads 8

    Prodigal procedures:
    - 'single': Best for isolate genomes (uses dynamic programming training)
    - 'meta': Best for metagenomic contigs (uses pre-trained models)

    Example workflow:

        # 1. Download genomes with protein files from NCBI
        metadarkmatter download genomes fetch \\
            --accessions genomes.tsv \\
            --output-dir genomes/ \\
            --include-protein

        # 2. Generate missing protein files with Prodigal
        metadarkmatter proteins predict \\
            --genomes genomes/ \\
            --missing-only \\
            --threads 8

        # 3. Now all genomes have protein files for AAI
        metadarkmatter aai compute \\
            --genomes genomes/ \\
            --output aai_matrix.csv
    """
    from metadarkmatter.external.prodigal import Prodigal

    out = QuietConsole(console, quiet=quiet)

    out.print("\n[bold blue]Metadarkmatter Protein Prediction[/bold blue]\n")

    # Check Prodigal availability
    if not Prodigal.check_available():
        console.print("[red]Error: Prodigal not found in PATH[/red]")
        console.print("[dim]Install with: conda install -c bioconda prodigal[/dim]")
        raise typer.Exit(code=1) from None

    # Set output directory
    if output_dir is None:
        output_dir = genomes
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find genomes to process
    if missing_only:
        genomes_to_process, existing = _find_missing_proteins(
            genomes, output_dir, genome_pattern
        )
        total_genomes = len(genomes_to_process) + len(existing)

        out.print(f"[bold]Mode:[/bold] Missing only")
        out.print(f"[bold]Total genomes:[/bold] {total_genomes}")
        out.print(f"[bold]Existing protein files:[/bold] {len(existing)}")
        out.print(f"[bold]Missing protein files:[/bold] {len(genomes_to_process)}")
    else:
        genomes_to_process = _find_genome_files(genomes, genome_pattern)
        out.print(f"[bold]Mode:[/bold] All genomes")
        out.print(f"[bold]Genomes found:[/bold] {len(genomes_to_process)}")

    out.print(f"[bold]Output directory:[/bold] {output_dir}")
    out.print(f"[bold]Procedure:[/bold] {procedure}")
    out.print(f"[bold]Threads:[/bold] {threads}")

    if not genomes_to_process:
        if missing_only:
            out.print("\n[green]All genomes already have protein files![/green]")
        else:
            console.print(
                f"\n[red]Error: No genome files found matching '{genome_pattern}' "
                f"in {genomes}[/red]"
            )
            raise typer.Exit(code=1) from None
        raise typer.Exit(code=0)

    # Load metadata for species information if available
    metadata_path = genomes.parent / "genome_metadata.tsv"
    if not metadata_path.exists():
        metadata_path = genomes / "genome_metadata.tsv"

    genome_species: dict[str, str] = {}
    if metadata_path.exists():
        try:
            import polars as pl
            metadata = pl.read_csv(metadata_path, separator="\t")
            for row in metadata.iter_rows(named=True):
                genome_species[row["accession"]] = row.get("species", "Unknown")
        except Exception:
            pass  # Metadata loading failed, continue without species info

    # Dry run
    if dry_run:
        out.print("\n[cyan]DRY RUN - Genomes that would be processed:[/cyan]")

        table = Table(show_header=True, header_style="bold")
        table.add_column("Accession", style="cyan")
        table.add_column("Species")
        table.add_column("Output File")

        for genome_file in genomes_to_process[:20]:
            accession = genome_file.stem
            if accession.endswith(".fna"):
                accession = accession[:-4]
            species = genome_species.get(accession, "")
            output_file = output_dir / f"{accession}.faa"
            table.add_row(accession, species, str(output_file.name))

        if len(genomes_to_process) > 20:
            table.add_row(
                f"... and {len(genomes_to_process) - 20} more",
                "",
                "",
                style="dim",
            )

        console.print(table)
        out.print("\n[green]Dry run complete. No files were created.[/green]")
        raise typer.Exit(code=0)

    # Run Prodigal in parallel
    out.print(f"\n[bold]Predicting proteins for {len(genomes_to_process)} genomes...[/bold]")

    results: list[tuple[Path, int, str | None]] = []
    errors: list[tuple[str, str]] = []

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        TimeElapsedColumn(),
        console=console if not quiet else None,
    ) as progress:
        task = progress.add_task(
            "Running Prodigal...",
            total=len(genomes_to_process),
        )

        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = {
                executor.submit(
                    _predict_single_genome,
                    genome_file,
                    output_dir,
                    procedure,
                ): genome_file
                for genome_file in genomes_to_process
            }

            for future in as_completed(futures):
                genome_file, protein_count, error = future.result()
                results.append((genome_file, protein_count, error))

                accession = genome_file.stem
                if accession.endswith(".fna"):
                    accession = accession[:-4]

                if error:
                    errors.append((accession, error))
                    if verbose:
                        out.print(f"  [red]Failed:[/red] {accession}: {error}")
                elif verbose:
                    species = genome_species.get(accession, "")
                    if species:
                        out.print(
                            f"  [green]Done:[/green] {accession} ({species}) - "
                            f"{protein_count:,} proteins"
                        )
                    else:
                        out.print(
                            f"  [green]Done:[/green] {accession} - {protein_count:,} proteins"
                        )

                progress.update(task, advance=1)

    # Summary
    successful = [r for r in results if r[2] is None]
    total_proteins = sum(r[1] for r in successful)

    out.print("\n[bold green]Prediction complete![/bold green]")
    out.print("\n[bold]Summary:[/bold]")
    out.print(f"  Genomes processed: {len(successful):,}")
    out.print(f"  Total proteins predicted: {total_proteins:,}")
    out.print(f"  Average proteins per genome: {total_proteins // max(len(successful), 1):,}")

    if errors:
        out.print(f"\n[yellow]Warnings: {len(errors)} genomes failed[/yellow]")
        if verbose:
            for accession, error in errors[:5]:
                out.print(f"  [dim]- {accession}: {error}[/dim]")
            if len(errors) > 5:
                out.print(f"  [dim]... and {len(errors) - 5} more[/dim]")
        else:
            out.print("[dim]Use --verbose to see error details[/dim]")

    # Verify all proteins now exist
    if missing_only:
        still_missing, now_existing = _find_missing_proteins(
            genomes, output_dir, genome_pattern
        )
        if still_missing:
            out.print(
                f"\n[yellow]Note: {len(still_missing)} genomes still missing "
                f"protein files (see errors above)[/yellow]"
            )
        else:
            out.print(
                f"\n[green]All {len(now_existing)} genomes now have protein files![/green]"
            )

    out.print(f"\n[dim]Output directory: {output_dir}[/dim]")
    out.print("[dim]Ready for: metadarkmatter aai compute --genomes " + str(output_dir) + "[/dim]")
    out.print()
