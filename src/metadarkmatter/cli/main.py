"""
Main CLI entry point for metadarkmatter.

Provides subcommands for each stage of the ANI-weighted placement pipeline:
- score: Classify BLAST results (core algorithm)
- kraken2: Taxonomic classification and read extraction
- blast: Run BLASTN competitive recruitment
- visualize: Generate recruitment plots
- run: Execute full pipeline
"""

from __future__ import annotations

import typer
from rich import print as rprint
from rich.console import Console

from metadarkmatter import __version__

app = typer.Typer(
    name="metadarkmatter",
    help="ANI-weighted placement for detecting novel microbial diversity in eDNA",
    add_completion=False,
    no_args_is_help=True,
)

console = Console()


def version_callback(value: bool) -> None:
    """Print version and exit."""
    if value:
        rprint(f"metadarkmatter version {__version__}")
        raise typer.Exit


@app.callback()
def main(
    version: bool = typer.Option(
        None,
        "--version",
        "-v",
        help="Show version and exit",
        callback=version_callback,
        is_eager=True,
    ),
) -> None:
    """
    Metadarkmatter: ANI-weighted placement for detecting novel microbial diversity.

    A tool for analyzing environmental DNA (eDNA) metagenomic data from air filters,
    water samples, and other environmental sources to characterize microbial diversity
    and detect novel bacterial taxa using whole-genome competitive read recruitment.
    """


# Import subcommands
from metadarkmatter.cli import aai, ani, blast, blastx, download, kraken2, mapping, mmseqs2, proteins, report, score, tree, visualize
from metadarkmatter.cli import map as map_cmd  # Alias to avoid shadowing builtin

# Register subcommands
app.add_typer(aai.app, name="aai")
app.add_typer(ani.app, name="ani")
app.add_typer(score.app, name="score")
app.add_typer(kraken2.app, name="kraken2")
app.add_typer(blast.app, name="blast")
app.add_typer(blastx.app, name="blastx")
app.add_typer(mmseqs2.app, name="mmseqs2")
app.add_typer(map_cmd.app, name="map")
app.add_typer(proteins.app, name="proteins")
app.add_typer(visualize.app, name="visualize")
app.add_typer(download.app, name="download")
app.add_typer(report.app, name="report")
app.add_typer(mapping.app, name="util")
app.add_typer(tree.app, name="tree")


if __name__ == "__main__":
    app()
