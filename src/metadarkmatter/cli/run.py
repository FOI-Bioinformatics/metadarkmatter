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
import yaml
from rich.console import Console
from rich.table import Table

from metadarkmatter import __version__
from metadarkmatter.cli.options import DryRun, Quiet
from metadarkmatter.cli.utils import QuietConsole
from metadarkmatter.core.pipeline import (
    Cohort,
    CohortLayout,
    CohortResult,
    Pipeline,
    PipelineLayout,
    RunManifest,
    StepStatus,
    read_samples_tsv,
)
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
    extracted_reads_2: Path | None = typer.Option(
        None, "--extracted-reads-2",
        help="Pre-extracted family R2 reads (paired-end; both aligners align both mates).",
        dir_okay=False,
    ),
    output_dir: Path | None = typer.Option(
        None, "--output-dir", "-o",
        help="Root output directory for the run (required unless set in --config).",
    ),
    aligner: Aligner | None = typer.Option(
        None, "--aligner", help="Aligner for the alignment step [default: blast].",
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
    threads: int | None = typer.Option(
        None, "--threads", "-t", min=1,
        help="Threads passed to each tool step [default: 8, or the --config value].",
    ),
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

    # CLI flags default to None so that values from --config are only
    # overridden when the user passes them explicitly; PipelineConfig supplies
    # the remaining defaults (aligner=blast, threads=8).
    if config is not None:
        cfg = PipelineConfig.from_yaml(config, **overrides)
        if extracted_reads_2 is None:
            extracted_reads_2 = _yaml_path(config, "extracted_reads_2")
    else:
        if output_dir is None:
            raise ValueError("Provide --output-dir (or output_dir in a --config YAML).")
        cfg = PipelineConfig(**{k: v for k, v in overrides.items() if v is not None})

    if extracted_reads_2 is not None and cfg.extracted_reads_1 is None:
        raise ValueError("--extracted-reads-2 requires --extracted-reads-1.")

    if cfg.genomes_dir is not None and not cfg.genomes_dir.is_dir():
        raise ValueError(f"--genomes path does not exist or is not a directory: {cfg.genomes_dir}")

    yaml_msgs = PipelineConfig.yaml_warnings(config) if config is not None else []

    layout = PipelineLayout(root=cfg.output_dir, family_taxid=cfg.family_taxid)
    pipeline = Pipeline(
        cfg,
        layout,
        dry_run=effective_dry_run,
        force=force,
        log=lambda msg: out.print(msg),
        extracted_reads_2=extracted_reads_2,
    )

    for msg in dict.fromkeys([*yaml_msgs, *getattr(pipeline, "warnings", [])]):
        console.print(f"[yellow]Warning:[/yellow] {msg}")

    if effective_dry_run:
        out.print("\n[bold cyan]DRY RUN[/bold cyan] - planned steps (no execution)\n")
    else:
        out.print(f"\n[bold blue]metadarkmatter pipeline[/bold blue] -> {cfg.output_dir}\n")

    manifest = pipeline.run(from_step=from_step, mdm_version=__version__)

    _print_summary(out, manifest, layout, effective_dry_run, quiet)


def _yaml_path(config: Path, key: str) -> Path | None:
    """Read an optional path key from the pipeline YAML (not a PipelineConfig field)."""
    raw = yaml.safe_load(config.read_text()) or {}
    value = raw.get(key) if isinstance(raw, dict) else None
    return Path(value) if value else None


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


@app.command(name="cohort")
def run_cohort(
    samples: Path = typer.Option(
        ..., "--samples", "-s", exists=True, dir_okay=False,
        help=(
            "Tab-separated sample sheet with columns sample, reads_1, reads_2 (optional) "
            "or extracted_reads_1, extracted_reads_2 (optional; skips Kraken2)."
        ),
    ),
    output_dir: Path = typer.Option(
        ..., "--output-dir", "-o",
        help="Cohort root: reference/, samples/<sample>/, cohort_summary.tsv, cohort_report.html.",
    ),
    family: str | None = typer.Option(
        None, "--family", help="GTDB family to download references for (e.g. f__Francisellaceae)."
    ),
    family_taxid: int | None = typer.Option(
        None, "--family-taxid", help="NCBI taxid for the family (needed for Kraken2 extraction)."
    ),
    kraken_db: Path | None = typer.Option(
        None, "--kraken-db", help="Kraken2 database directory.", file_okay=False,
    ),
    genomes: Path | None = typer.Option(
        None, "--genomes", "-g", help="Pre-downloaded reference genomes dir (skips download).",
        file_okay=False,
    ),
    aligner: Aligner = typer.Option(
        Aligner.BLAST, "--aligner", help="Aligner for the alignment step.",
    ),
    no_report: bool = typer.Option(
        False, "--no-report", help="Skip per-sample reports and the cohort report.",
    ),
    no_sample_reports: bool = typer.Option(
        False, "--no-sample-reports",
        help="Skip per-sample HTML reports but still build the cohort report.",
    ),
    from_step: PipelineStep | None = typer.Option(
        None, "--from", help="Resume from this step (earlier outputs must already exist).",
    ),
    force: bool = typer.Option(
        False, "--force", help="Re-run every step even if its outputs already exist.",
    ),
    threads: int = typer.Option(
        8, "--threads", "-t", min=1, help="Threads passed to each tool step.",
    ),
    quiet: Quiet = False,
    dry_run: DryRun = False,
) -> None:
    """Run the pipeline over many samples that share one reference.

    The reference (genomes, alignment database, ANI matrix) is built once under
    OUTPUT_DIR/reference; each sample then runs Kraken2, alignment, classification
    and an optional report under OUTPUT_DIR/samples/<sample>/. Samples run
    sequentially; a failing sample does not stop the rest, but the command then
    exits non-zero. Finishes with cohort_summary.tsv and cohort_report.html.

    Example:

        # samples.tsv: sample<TAB>reads_1<TAB>reads_2
        mdm run cohort --samples samples.tsv --family f__Francisellaceae \\
            --family-taxid 34064 --kraken-db /db/kraken2 --output-dir results/cohort01

        # Pre-extracted reads and pre-downloaded genomes (no Kraken2)
        mdm run cohort --samples extracted.tsv --genomes genomes/ --output-dir results/cohort02
    """
    out = QuietConsole(console, quiet=quiet)
    effective_dry_run = dry_run or is_dry_run()

    sheet = read_samples_tsv(samples)
    if genomes is not None and not genomes.is_dir():
        raise ValueError(f"--genomes path does not exist or is not a directory: {genomes}")

    shared: dict[str, object] = {
        "family": family,
        "family_taxid": family_taxid,
        "kraken_db": kraken_db,
        "genomes_dir": genomes,
        "threads": threads,
    }
    cohort = Cohort(
        sheet,
        {k: v for k, v in shared.items() if v is not None},
        CohortLayout(root=output_dir),
        aligner=aligner,
        dry_run=effective_dry_run,
        force=force,
        log=lambda msg: out.print(msg),
        sample_reports=not (no_report or no_sample_reports),
        cohort_report=not no_report,
    )

    if effective_dry_run:
        out.print("\n[bold cyan]DRY RUN[/bold cyan] - planned cohort steps (no execution)\n")
    else:
        out.print(
            f"\n[bold blue]metadarkmatter cohort[/bold blue] ({len(sheet)} samples) "
            f"-> {output_dir}\n"
        )

    result = cohort.run(from_step=from_step, mdm_version=__version__)
    _print_cohort_summary(out, result, cohort.layout, effective_dry_run, quiet)

    if not effective_dry_run and not result.ok:
        names = ", ".join(s.name for s in result.failed)
        if names:
            console.print(f"[red]Cohort finished with failed samples: {names}[/red]")
        if result.report_error:
            console.print(f"[red]Cohort report: {result.report_error}[/red]")
        raise typer.Exit(code=1)


def _print_cohort_summary(
    out: QuietConsole,
    result: CohortResult,
    layout: CohortLayout,
    dry_run: bool,
    quiet: bool,
) -> None:
    if quiet:
        return
    table = Table(title="Cohort samples", show_header=True, header_style="bold")
    table.add_column("Sample")
    table.add_column("Status")
    table.add_column("Detail")
    style = {"completed": "green", "failed": "red", "planned": "cyan"}
    for s in result.samples:
        table.add_row(s.name, f"[{style.get(s.status, 'white')}]{s.status}[/]", s.error or "")
    console.print(table)
    if dry_run:
        out.print("\n[dim]Dry run complete. Re-run without --dry-run to execute.[/dim]")
        return
    out.print(f"  Summary: {layout.summary_tsv}")
    if result.report_path is not None:
        out.print(f"  Report:  {result.report_path}")
