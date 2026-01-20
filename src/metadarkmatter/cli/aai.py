"""
AAI command for computing and validating AAI matrices.

Average Amino Acid Identity (AAI) provides a more reliable metric than ANI
for genus-level classification. While ANI becomes unreliable below approximately
80% identity, AAI maintains accuracy at the genus level.

Provides subcommands:
- compute: Build AAI matrix from genome protein files using Diamond
- validate: Check AAI matrix coverage against BLAST results
"""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, TimeElapsedColumn
from rich.table import Table

from metadarkmatter.cli.utils import QuietConsole
from metadarkmatter.external import ToolExecutionError

app = typer.Typer(
    name="aai",
    help="Compute and validate AAI matrices for genus-level classification",
    no_args_is_help=True,
)

console = Console()


def _find_protein_files(
    protein_dir: Path,
    pattern: str,
) -> list[Path]:
    """Find protein FASTA files in directory, trying multiple patterns if needed."""
    protein_files = sorted(protein_dir.glob(pattern))

    if not protein_files:
        # Try alternative patterns
        for alt_pattern in ["*.faa", "*.faa.gz", "*.fasta", "*.fasta.gz"]:
            if alt_pattern != pattern:
                protein_files = sorted(protein_dir.glob(alt_pattern))
                if protein_files:
                    break

    return protein_files


@app.command(name="compute")
def compute(
    genomes: Path = typer.Option(
        ...,
        "--genomes",
        "-g",
        help="Directory containing genome protein FASTA files (.faa)",
        exists=True,
        file_okay=False,
    ),
    output: Path = typer.Option(
        ...,
        "--output",
        "-o",
        help="Output CSV file for AAI matrix",
    ),
    threads: int = typer.Option(
        4,
        "--threads",
        "-t",
        help="Number of threads for Diamond",
        min=1,
    ),
    min_identity: float = typer.Option(
        30.0,
        "--min-identity",
        "-i",
        help="Minimum percent identity for Diamond hits",
        min=0.0,
        max=100.0,
    ),
    min_coverage: float = typer.Option(
        0.5,
        "--min-coverage",
        "-c",
        help="Minimum query coverage fraction for reciprocal best hits",
        min=0.0,
        max=1.0,
    ),
    evalue: float = typer.Option(
        1e-5,
        "--evalue",
        "-e",
        help="E-value threshold for Diamond",
        min=0.0,
    ),
    protein_pattern: str = typer.Option(
        "*.faa",
        "--protein-pattern",
        "-p",
        help="Glob pattern for protein FASTA files",
    ),
    keep_intermediates: bool = typer.Option(
        False,
        "--keep-intermediates",
        "-k",
        help="Keep intermediate files (Diamond output, RBH pairs)",
    ),
    sensitive: bool = typer.Option(
        False,
        "--sensitive",
        "-s",
        help="Use Diamond sensitive mode (slower but more accurate for divergent sequences)",
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
        help="Show what would be computed without executing",
    ),
) -> None:
    """
    Compute AAI matrix from genome protein files.

    AAI (Average Amino Acid Identity) is computed using Diamond blastp
    with reciprocal best hit (RBH) filtering. AAI provides more reliable
    genus-level classification than ANI, which becomes unreliable below
    approximately 80% identity.

    AAI genus boundaries (Riesco & Trujillo 2024):
    - Same genus: AAI > 65%
    - Genus boundary: AAI 58-65%
    - Different genus: AAI < 58%

    Requires protein FASTA files (.faa) for each genome. These can be
    downloaded from NCBI along with genome sequences, or generated using
    Prodigal gene prediction.

    Example:

        # Compute AAI matrix from protein files
        metadarkmatter aai compute \\
            --genomes proteins/ \\
            --output aai_matrix.csv \\
            --threads 16

        # Use with classification
        metadarkmatter score classify \\
            --blast results.blast.tsv.gz \\
            --ani ani_matrix.csv \\
            --aai aai_matrix.csv \\
            --output classifications.csv
    """
    out = QuietConsole(console, quiet=quiet)

    out.print("\n[bold blue]Metadarkmatter AAI Matrix Builder[/bold blue]\n")

    # Check Diamond availability
    from metadarkmatter.external.diamond import Diamond

    if not Diamond.check_available():
        console.print("[red]Error: Diamond not found in PATH[/red]")
        console.print("[dim]Install with: conda install -c bioconda diamond[/dim]")
        raise typer.Exit(code=1) from None

    # Find protein files
    protein_files = _find_protein_files(genomes, protein_pattern)

    if not protein_files:
        console.print(
            f"[red]Error: No protein files found matching '{protein_pattern}' "
            f"in {genomes}[/red]"
        )
        console.print(
            "\n[dim]Protein files can be obtained by:[/dim]\n"
            "  1. Download from NCBI with metadarkmatter:\n"
            "     metadarkmatter download genomes fetch --accessions genomes.tsv \\\n"
            "       --output-dir genomes/ --include-protein\n"
            "  2. Download from NCBI directly:\n"
            "     datasets download genome accession GCF_XXX --include protein\n"
            "  3. Generate with Prodigal:\n"
            "     prodigal -i genome.fna -a proteins.faa"
        )
        raise typer.Exit(code=1) from None

    out.print(f"[bold]Input:[/bold] {len(protein_files)} protein files from {genomes}")
    out.print(f"[bold]Threads:[/bold] {threads}")
    out.print(f"[bold]Parameters:[/bold] min_identity={min_identity}%, min_coverage={min_coverage}")

    # Create output directory
    output.parent.mkdir(parents=True, exist_ok=True)

    if dry_run:
        out.print("\n[cyan]DRY RUN - Steps that would be executed:[/cyan]")
        out.print("  1. Concatenate protein files with genome-prefixed headers")
        out.print("  2. Create Diamond database from concatenated proteins")
        out.print(f"  3. Run all-vs-all Diamond blastp (threads={threads})")
        out.print("  4. Compute reciprocal best hits (RBH)")
        out.print("  5. Calculate AAI from mean RBH identity per genome pair")
        out.print(f"\nOutput would be written to: {output}")
        raise typer.Exit(code=0)

    # Run AAI computation
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        TimeElapsedColumn(),
        console=console if not quiet else None,
    ) as progress:
        progress.add_task(
            description=f"Computing AAI with Diamond ({len(protein_files)} genomes)...",
            total=None,
        )

        try:
            from metadarkmatter.core.aai_matrix_builder import compute_aai_matrix

            aai_dict = compute_aai_matrix(
                protein_dir=genomes,
                output_path=output,
                threads=threads,
                min_identity=min_identity,
                min_coverage=min_coverage,
                evalue=evalue,
                keep_intermediates=keep_intermediates,
                sensitive=sensitive,
            )

        except FileNotFoundError as e:
            console.print(f"\n[red]Error: {e}[/red]")
            raise typer.Exit(code=1) from None
        except ToolExecutionError as e:
            console.print(f"\n[red]Diamond failed:[/red]\n{e.message}")
            raise typer.Exit(code=1) from None

    # Summary
    n_genomes = len(aai_dict)
    out.print("\n[bold green]AAI matrix computed successfully![/bold green]")
    out.print(f"\n[bold]Output:[/bold] {output}")
    out.print(f"[bold]Genomes:[/bold] {n_genomes}")

    file_size = output.stat().st_size / (1024 * 1024)
    out.print(f"[bold]Size:[/bold] {file_size:.2f} MB")

    if verbose:
        # Show sample AAI values
        import polars as pl

        df = pl.read_csv(output)
        out.print(f"\n[dim]Matrix shape: {df.shape[0]} x {df.shape[1] - 1}[/dim]")

        # Show AAI statistics
        all_values = []
        for col in df.columns[1:]:
            all_values.extend(df[col].to_list())
        # Exclude diagonal (100.0) and zeros
        non_diag = [v for v in all_values if v != 100.0 and v > 0]
        if non_diag:
            import statistics
            out.print(f"[dim]AAI range: {min(non_diag):.1f}% - {max(non_diag):.1f}%[/dim]")
            out.print(f"[dim]AAI mean: {statistics.mean(non_diag):.1f}%[/dim]")

    out.print(f"\n[dim]Use with: metadarkmatter score classify --aai {output}[/dim]")
    out.print()


@app.command(name="validate")
def validate(
    aai: Path = typer.Option(
        ...,
        "--aai",
        "-a",
        help="AAI matrix CSV file to validate",
        exists=True,
        dir_okay=False,
    ),
    ani: Path = typer.Option(
        ...,
        "--ani",
        help="ANI matrix CSV file to compare coverage against",
        exists=True,
        dir_okay=False,
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="List genomes with missing AAI/ANI pairs",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet",
        "-q",
        help="Suppress progress output",
    ),
) -> None:
    """
    Validate AAI matrix coverage against ANI matrix.

    Checks that genomes appearing in the ANI matrix have corresponding
    entries in the AAI matrix. Reports missing genomes and coverage statistics.

    For best results, both matrices should cover the same set of genomes.

    Example:

        metadarkmatter aai validate \\
            --aai aai_matrix.csv \\
            --ani ani_matrix.csv
    """
    out = QuietConsole(console, quiet=quiet)

    out.print("\n[bold blue]Metadarkmatter AAI Matrix Validator[/bold blue]\n")

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

    # Load AAI matrix
    out.print(f"\n[bold]Loading AAI matrix:[/bold] {aai}")
    try:
        from metadarkmatter.core.aai_matrix_builder import AAIMatrix

        aai_matrix = AAIMatrix.from_file(aai)
        aai_genomes = aai_matrix.genomes
        out.print(f"  [green]Loaded {len(aai_genomes)} genomes[/green]")
    except Exception as e:
        console.print(f"[red]Error loading AAI matrix: {e}[/red]")
        raise typer.Exit(code=1) from None

    # Calculate coverage
    matched = ani_genomes & aai_genomes
    missing_from_aai = ani_genomes - aai_genomes
    extra_in_aai = aai_genomes - ani_genomes

    coverage_pct = (100.0 * len(matched) / len(ani_genomes)) if len(ani_genomes) > 0 else 0.0

    # Display results
    out.print("\n[bold]Coverage Analysis:[/bold]")

    table = Table(show_header=True, header_style="bold")
    table.add_column("Metric", style="cyan")
    table.add_column("Value", justify="right")

    table.add_row("ANI matrix genomes", f"{len(ani_genomes):,}")
    table.add_row("AAI matrix genomes", f"{len(aai_genomes):,}")
    table.add_row("Matched", f"{len(matched):,}")
    table.add_row("Missing from AAI", f"{len(missing_from_aai):,}")
    table.add_row("Extra in AAI", f"{len(extra_in_aai):,}")
    table.add_row("Coverage", f"{coverage_pct:.1f}%")

    console.print(table)

    # Warnings/recommendations
    if coverage_pct < 50.0:
        console.print(f"\n[yellow]Warning: Low coverage ({coverage_pct:.1f}%)[/yellow]")
        console.print(
            "  Many genomes in ANI matrix are missing from AAI matrix."
        )
        console.print(
            "  Ensure protein files are available for all reference genomes."
        )
    elif coverage_pct < 90.0:
        console.print(f"\n[yellow]Note: Moderate coverage ({coverage_pct:.1f}%)[/yellow]")
        console.print("  Some genomes in ANI matrix are missing from AAI matrix.")
    else:
        console.print(f"\n[green]Good coverage ({coverage_pct:.1f}%)[/green]")

    # Show missing genomes
    if missing_from_aai and verbose:
        console.print("\n[bold]Missing from AAI matrix:[/bold]")
        for genome in sorted(missing_from_aai)[:20]:
            console.print(f"  - {genome}")
        if len(missing_from_aai) > 20:
            console.print(f"  ... and {len(missing_from_aai) - 20} more")
    elif missing_from_aai:
        console.print(
            f"\n[dim]Use --verbose to see list of {len(missing_from_aai)} missing genomes[/dim]"
        )

    # Show AAI statistics for matched genomes
    if matched and verbose:
        console.print("\n[bold]AAI Statistics (matched genomes):[/bold]")
        aai_values = []
        for g1 in list(matched)[:100]:  # Sample up to 100 genomes
            for g2 in list(matched)[:100]:
                if g1 != g2:
                    aai_val = aai_matrix.get_aai(g1, g2)
                    if aai_val > 0:
                        aai_values.append(aai_val)

        if aai_values:
            import statistics
            console.print(f"  Range: {min(aai_values):.1f}% - {max(aai_values):.1f}%")
            console.print(f"  Mean: {statistics.mean(aai_values):.1f}%")
            console.print(f"  Median: {statistics.median(aai_values):.1f}%")

    # Exit code based on coverage
    if coverage_pct < 50.0:
        raise typer.Exit(code=1) from None

    out.print()
