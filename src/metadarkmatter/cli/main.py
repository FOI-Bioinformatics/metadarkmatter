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
from metadarkmatter.cli.errors import wrap_app_commands
from metadarkmatter.core.logging_config import setup_logging
from metadarkmatter.core.runtime import set_debug, set_dry_run

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
    log_level: str = typer.Option(
        "WARNING",
        "--log-level",
        help="Logging level (DEBUG, INFO, WARNING, ERROR)",
    ),
    log_format: str = typer.Option(
        "text",
        "--log-format",
        help="Log format: text or json",
    ),
    log_file: str = typer.Option(
        None,
        "--log-file",
        help="Log file path",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help=(
            "Process-wide dry-run flag. Subcommands that support dry-run "
            "honour it; equivalent to setting MDM_DRY_RUN=1 in the "
            "environment."
        ),
    ),
    debug: bool = typer.Option(
        False,
        "--debug",
        help=(
            "Print full tracebacks on error instead of a friendly message. "
            "Equivalent to setting MDM_DEBUG=1 in the environment."
        ),
    ),
) -> None:
    """
    Metadarkmatter: ANI-weighted placement for detecting novel microbial diversity.

    A tool for analyzing environmental DNA (eDNA) metagenomic data from air filters,
    water samples, and other environmental sources to characterize microbial diversity
    and detect novel bacterial taxa using whole-genome competitive read recruitment.
    """
    setup_logging(log_level, log_format, log_file)
    if dry_run:
        set_dry_run(True)
    if debug:
        set_debug(True)


# Import subcommands
from metadarkmatter.cli import (
    aai,
    ani,
    blast,
    blastx,
    doctor,
    download,
    kraken2,
    mapping,
    mmseqs2,
    proteins,
    report,
    run,
    score,
    tree,
    validate,
    visualize,
)
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
app.add_typer(validate.app, name="validate")
app.add_typer(doctor.app, name="doctor")
app.add_typer(run.app, name="run")

# Wrap every registered command with the centralized error handler so that
# MetadarkmatterError subclasses and common library/runtime errors surface as
# friendly messages (full tracebacks only under --debug). Done once here to
# avoid per-command decoration churn across the subcommand modules.
for _sub_app in (
    aai.app, ani.app, score.app, kraken2.app, blast.app, blastx.app,
    mmseqs2.app, map_cmd.app, proteins.app, visualize.app, download.app,
    report.app, mapping.app, tree.app, validate.app, doctor.app, run.app,
):
    wrap_app_commands(_sub_app)


if __name__ == "__main__":
    app()
