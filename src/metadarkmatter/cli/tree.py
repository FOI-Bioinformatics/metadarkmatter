"""
Tree command for building phylogenetic trees.

Provides subcommands:
- build: Build a phylogenetic tree from ANI matrix or genome assemblies
"""
from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console

from metadarkmatter.cli.utils import QuietConsole, spinner_progress
from metadarkmatter.core.phylogeny.tree_builder import TreeMethod

app = typer.Typer(
    name="tree",
    help="Build phylogenetic trees from genomes or ANI matrices",
    no_args_is_help=True,
)

console = Console()


@app.command(name="build")
def build(
    method: TreeMethod = typer.Option(
        ...,
        "--method",
        "-m",
        help="Tree building method: nj, upgma, or mashtree",
    ),
    ani: Path | None = typer.Option(
        None,
        "--ani",
        "-a",
        help="ANI matrix CSV file (required for nj/upgma methods)",
        exists=True,
        dir_okay=False,
    ),
    genomes: Path | None = typer.Option(
        None,
        "--genomes",
        "-g",
        help="Directory containing genome FASTA files (required for mashtree method)",
        exists=True,
        file_okay=False,
    ),
    output: Path = typer.Option(
        ...,
        "--output",
        "-o",
        help="Output Newick tree file",
    ),
    threads: int = typer.Option(
        4,
        "--threads",
        "-t",
        help="Number of threads (mashtree only)",
        min=1,
    ),
    genome_pattern: str = typer.Option(
        "*.fna",
        "--genome-pattern",
        "-p",
        help="Glob pattern for genome files (mashtree only)",
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
) -> None:
    """
    Build a phylogenetic tree.

    Supports three methods:

    - nj: Neighbor-joining from ANI distance matrix (BioPython)

    - upgma: UPGMA from ANI distance matrix (BioPython, ultrametric)

    - mashtree: Mash distance-based NJ from genome assemblies (requires mashtree)

    Examples:

        # NJ tree from ANI matrix
        metadarkmatter tree build --method nj --ani ani_matrix.csv --output tree.nwk

        # UPGMA tree
        metadarkmatter tree build --method upgma --ani ani_matrix.csv --output tree.nwk

        # Mashtree from genome files
        metadarkmatter tree build --method mashtree --genomes genomes/ --output tree.nwk -t 16
    """
    import polars as pl

    from metadarkmatter.core.phylogeny.tree_builder import build_tree

    out = QuietConsole(console, quiet=quiet)

    out.print("\n[bold blue]Metadarkmatter Tree Builder[/bold blue]\n")
    out.print(f"[bold]Method:[/bold] {method.value}")

    # Validate inputs
    if method in (TreeMethod.NJ, TreeMethod.UPGMA):
        if ani is None:
            console.print("[red]Error: --ani is required for nj/upgma methods[/red]")
            raise typer.Exit(code=1) from None

        # Load ANI matrix
        out.print(f"[bold]ANI matrix:[/bold] {ani}")
        try:
            ani_df = pl.read_csv(ani)
            # First column is the genome names (index)
            genome_names = ani_df.get_column(ani_df.columns[0]).to_list()
            ani_values = ani_df.select(ani_df.columns[1:])
            import pandas as pd

            ani_matrix = pd.DataFrame(
                ani_values.to_numpy(),
                index=genome_names,
                columns=genome_names,
            )
        except Exception as e:
            console.print(f"[red]Error loading ANI matrix: {e}[/red]")
            raise typer.Exit(code=1) from None

        n_genomes = len(ani_matrix)
        out.print(f"[bold]Genomes:[/bold] {n_genomes}")

        if n_genomes < 3:
            console.print(
                f"[red]Error: Need >= 3 genomes for tree building (got {n_genomes})[/red]"
            )
            raise typer.Exit(code=1) from None

        with spinner_progress(
            f"Building {method.value.upper()} tree from {n_genomes} genomes...",
            console,
            quiet,
        ):
            newick = build_tree(method, ani_matrix=ani_matrix)

    elif method == TreeMethod.MASHTREE:
        if genomes is None:
            console.print("[red]Error: --genomes is required for mashtree method[/red]")
            raise typer.Exit(code=1) from None

        out.print(f"[bold]Genomes:[/bold] {genomes}")
        out.print(f"[bold]Threads:[/bold] {threads}")

        with spinner_progress(
            f"Building Mashtree tree from {genomes}...",
            console,
            quiet,
        ):
            newick = build_tree(
                method,
                genome_dir=genomes,
                genome_pattern=genome_pattern,
                threads=threads,
            )
    else:
        console.print(f"[red]Error: Unknown method '{method}'[/red]")
        raise typer.Exit(code=1) from None

    if newick is None:
        console.print("[red]Error: Tree construction failed[/red]")
        raise typer.Exit(code=1) from None

    # Write output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(newick + "\n")

    out.print("\n[bold green]Tree built successfully![/bold green]")
    out.print(f"[bold]Output:[/bold] {output}")
    out.print(f"\n[dim]Use with: metadarkmatter report generate --tree {output}[/dim]")
    out.print()
