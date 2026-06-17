"""
``mdm run`` - end-to-end pipeline orchestration.

Runs the full workflow (download -> kraken2 extract -> align -> ANI ->
classify -> report) with one command, reusing the existing ``mdm`` subcommands
as steps. Supports a structured output directory, per-step checkpointing,
``--from`` resume, ``--dry-run`` planning, and a JSON run manifest.
"""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from metadarkmatter import __version__
from metadarkmatter.cli.options import DryRun, Quiet, Threads
from metadarkmatter.cli.utils import QuietConsole
from metadarkmatter.core.pipeline import Pipeline, PipelineLayout, RunManifest, StepStatus
from metadarkmatter.core.runtime import is_dry_run
from metadarkmatter.models.pipeline_config import Aligner, PipelineConfig, PipelineStep

app = typer.Typer(
    name="run",
    help="Run the full metadarkmatter pipeline end-to-end",
    no_args_is_help=True,
)

console = Console()


@app.command(name="pipeline")
def run_pipeline(
    config: Path | None = typer.Option(
        None, "--config", "-c", help="Pipeline config YAML (CLI flags override its values).",
        exists=True, dir_okay=False,
    ),
    family: str | None = typer.Option(
        None, "--family", help="GTDB family to download references for (e.g. f__Francisellaceae)."
    ),
    family_taxid: int | None = typer.Option(
        None, "--family-taxid", help="NCBI taxid for the family (needed for Kraken2 extraction)."
    ),
    reads_1: Path | None = typer.Option(
        None, "--reads-1", "-1", help="R1 reads (FASTQ[.gz]).", exists=True, dir_okay=False,
    ),
    reads_2: Path | None = typer.Option(
        None, "--reads-2", "-2", help="R2 reads for paired-end input.", exists=True, dir_okay=False,
    ),
    kraken_db: Path | None = typer.Option(
        None, "--kraken-db", help="Kraken2 database directory.", file_okay=False,
    ),
    genomes: Path | None = typer.Option(
        None, "--genomes", "-g", help="Pre-downloaded reference genomes dir (skips download).",
        file_okay=False,
    ),
    kraken_output: Path | None = typer.Option(
        None, "--kraken-output", help="Pre-computed .kraken file (skips kraken classify).",
        dir_okay=False,
    ),
    extracted_reads_1: Path | None = typer.Option(
        None, "--extracted-reads-1", help="Pre-extracted family R1 reads (skips both kraken steps).",
        dir_okay=False,
    ),
    output_dir: Path = typer.Option(
        ..., "--output-dir", "-o", help="Root output directory for the run.",
    ),
    aligner: Aligner = typer.Option(
        Aligner.BLAST, "--aligner", help="Aligner for the alignment step.",
    ),
    no_report: bool = typer.Option(
        False, "--no-report", help="Skip the HTML report step.",
    ),
    from_step: PipelineStep | None = typer.Option(
        None, "--from", help="Resume from this step (earlier outputs must already exist).",
    ),
    force: bool = typer.Option(
        False, "--force", help="Re-run every step even if its outputs already exist.",
    ),
    threads: Threads = 8,
    quiet: Quiet = False,
    dry_run: DryRun = False,
) -> None:
    """Run the complete pipeline from reads (or extracted reads) to an HTML report.

    Examples:

        # Full run from raw reads (downloads references, extracts, aligns, classifies)
        mdm run pipeline --family f__Francisellaceae --family-taxid 34064 \\
            --reads-1 sample_R1.fastq.gz --reads-2 sample_R2.fastq.gz \\
            --kraken-db /db/kraken2 --output-dir results/run01

        # Re-use pre-downloaded genomes and pre-extracted reads
        mdm run pipeline --genomes genomes/ --extracted-reads-1 reads_R1.fastq.gz \\
            --output-dir results/run02

        # Plan only (no execution)
        mdm run pipeline --config run.yaml --dry-run
    """
    out = QuietConsole(console, quiet=quiet)
    effective_dry_run = dry_run or is_dry_run()

    overrides = {
        "family": family,
        "family_taxid": family_taxid,
        "reads_1": reads_1,
        "reads_2": reads_2,
        "kraken_db": kraken_db,
        "genomes_dir": genomes,
        "kraken_output": kraken_output,
        "extracted_reads_1": extracted_reads_1,
        "output_dir": output_dir,
        "aligner": aligner,
        "threads": threads,
        "report": False if no_report else None,
    }

    if config is not None:
        cfg = PipelineConfig.from_yaml(config, **overrides)
    else:
        cfg = PipelineConfig(**{k: v for k, v in overrides.items() if v is not None})

    layout = PipelineLayout(root=cfg.output_dir)
    pipeline = Pipeline(
        cfg,
        layout,
        dry_run=effective_dry_run,
        force=force,
        log=lambda msg: out.print(msg),
    )

    if effective_dry_run:
        out.print("\n[bold cyan]DRY RUN[/bold cyan] - planned steps (no execution)\n")
    else:
        out.print(f"\n[bold blue]metadarkmatter pipeline[/bold blue] -> {cfg.output_dir}\n")

    manifest = pipeline.run(from_step=from_step, mdm_version=__version__)

    _print_summary(out, manifest, layout, effective_dry_run, quiet)


def _print_summary(
    out: QuietConsole,
    manifest: RunManifest,
    layout: PipelineLayout,
    dry_run: bool,
    quiet: bool,
) -> None:
    if quiet:
        return

    table = Table(title="Pipeline steps", show_header=True, header_style="bold")
    table.add_column("Step")
    table.add_column("Status")
    table.add_column("Time (s)", justify="right")

    status_style = {
        StepStatus.COMPLETED.value: "green",
        StepStatus.CACHED.value: "cyan",
        StepStatus.SKIPPED.value: "dim",
        StepStatus.FAILED.value: "red",
    }
    for step in manifest.steps:
        style = status_style.get(step.status.value, "white")
        elapsed = f"{step.elapsed_s:.1f}" if step.elapsed_s else "-"
        table.add_row(step.name, f"[{style}]{step.status.value}[/{style}]", elapsed)

    console.print(table)

    if dry_run:
        out.print("\n[dim]Dry run complete. Re-run without --dry-run to execute.[/dim]")
        return

    out.print("\n[bold green]Pipeline complete.[/bold green]")
    if layout.classifications.exists():
        out.print(f"  Classifications: {layout.classifications}")
    if layout.report.exists():
        out.print(f"  Report:          {layout.report}")
    out.print(f"  Manifest:        {layout.manifest}")
