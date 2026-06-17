"""
Pipeline orchestrator for the ``mdm run`` command.

Productizes ``scripts/run_pipeline.sh`` as a first-class, resumable command.
Each step invokes an existing ``mdm`` subcommand (so there is no duplication of
per-command logic), wrapped with:

- a structured output directory layout (``PipelineLayout``),
- per-step checkpointing (skip when declared outputs already exist),
- ``--from`` resume and dry-run support,
- a JSON run manifest (``RunManifest``) recording provenance and per-step status.

The step runner is injectable so the orchestration logic can be unit-tested
without the external tools installed.
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from collections.abc import Callable, Sequence
from dataclasses import asdict, dataclass, field
from enum import Enum, StrEnum
from pathlib import Path

from metadarkmatter.core.exceptions import MetadarkmatterError
from metadarkmatter.models.pipeline_config import Aligner, PipelineConfig, PipelineStep

# Invoke the package CLI as a subprocess via the module entry point. This is
# robust to how the console script is installed (works under `uv run`, in the
# container, and from a venv) because cli/main.py has an `if __name__ ==
# "__main__"` guard.
MDM_INVOCATION: tuple[str, ...] = (sys.executable, "-m", "metadarkmatter.cli.main")

# Runner signature: (argv, dry_run) -> exit code.
StepRunner = Callable[[Sequence[str], bool], int]


class PipelineStepError(MetadarkmatterError):
    """Raised when a pipeline step fails (non-zero exit) or a precondition is unmet."""


class StepStatus(StrEnum):
    COMPLETED = "completed"
    CACHED = "cached"  # outputs already present, step skipped
    SKIPPED = "skipped"  # not applicable / disabled / dry-run
    FAILED = "failed"


@dataclass
class StepResult:
    name: str
    status: StepStatus
    command: str | None = None
    outputs: list[str] = field(default_factory=list)
    elapsed_s: float = 0.0
    detail: str | None = None

    def to_dict(self) -> dict[str, object]:
        d = asdict(self)
        d["status"] = self.status.value
        return d


@dataclass
class PipelineLayout:
    """Single source of truth for all pipeline file paths under one root."""

    root: Path

    @property
    def genomes_dir(self) -> Path:
        return self.root / "genomes"

    @property
    def genomes_tsv(self) -> Path:
        return self.root / "genomes.tsv"

    @property
    def metadata_tsv(self) -> Path:
        return self.root / "genome_metadata.tsv"

    @property
    def kraken_output(self) -> Path:
        return self.root / "sample.kraken"

    @property
    def extraction_dir(self) -> Path:
        return self.root / "extraction"

    @property
    def extracted_reads_1(self) -> Path:
        return self.extraction_dir / "reads_R1.fastq.gz"

    @property
    def extracted_reads_2(self) -> Path:
        return self.extraction_dir / "reads_R2.fastq.gz"

    def database(self, aligner: Aligner) -> Path:
        sub = "blastdb" if aligner == Aligner.BLAST else "mmseqs_db"
        return self.root / sub / "pangenome"

    def alignment(self, aligner: Aligner) -> Path:
        return self.root / f"sample.{aligner.value}.tsv.gz"

    @property
    def ani_matrix(self) -> Path:
        return self.root / "ani_matrix.csv"

    @property
    def aai_matrix(self) -> Path:
        return self.root / "aai_matrix.csv"

    @property
    def classifications(self) -> Path:
        return self.root / "classifications.csv"

    @property
    def summary_json(self) -> Path:
        return self.root / "summary.json"

    @property
    def report(self) -> Path:
        return self.root / "reports" / "report.html"

    @property
    def manifest(self) -> Path:
        return self.root / "run_manifest.json"


@dataclass
class RunManifest:
    mdm_version: str
    output_dir: str
    aligner: str
    config: dict[str, object]
    steps: list[StepResult] = field(default_factory=list)
    started_at: str | None = None
    finished_at: str | None = None

    def to_json(self) -> str:
        return json.dumps(
            {
                "mdm_version": self.mdm_version,
                "output_dir": self.output_dir,
                "aligner": self.aligner,
                "started_at": self.started_at,
                "finished_at": self.finished_at,
                "config": self.config,
                "steps": [s.to_dict() for s in self.steps],
            },
            indent=2,
            default=str,
        )


def _default_runner(argv: Sequence[str], dry_run: bool) -> int:
    """Run a subcommand, streaming its output to the caller's stdout/stderr."""
    if dry_run:
        return 0
    completed = subprocess.run(list(argv), check=False)
    return completed.returncode


def _is_present(path: Path) -> bool:
    """A path counts as a produced output if it exists and is non-empty.

    Directories must contain at least one file; this mirrors the ``[[ -s ]]``
    checkpoint test in run_pipeline.sh and the ``--skip-if-exists`` semantics.
    """
    if path.is_dir():
        return any(p.is_file() for p in path.rglob("*"))
    return path.exists() and path.stat().st_size > 0


def _db_present(db_prefix: Path) -> bool:
    """A makedb output is present if any file with the db prefix exists."""
    parent = db_prefix.parent
    if not parent.is_dir():
        return False
    return any(p.name.startswith(db_prefix.name) and p.stat().st_size > 0 for p in parent.iterdir())


@dataclass
class _Step:
    step: PipelineStep
    enabled: bool
    argv: list[str]
    outputs: list[Path]
    present: Callable[[], bool]


class Pipeline:
    """Build and execute the ordered set of pipeline steps for a config."""

    def __init__(
        self,
        config: PipelineConfig,
        layout: PipelineLayout,
        *,
        dry_run: bool = False,
        force: bool = False,
        runner: StepRunner = _default_runner,
        log: Callable[[str], None] | None = None,
    ) -> None:
        self.config = config
        self.layout = layout
        self.dry_run = dry_run
        self.force = force
        self._runner = runner
        self._log = log or (lambda _msg: None)

    # ------------------------------------------------------------------
    # Step planning
    # ------------------------------------------------------------------
    def plan(self) -> list[_Step]:
        """Build the ordered list of steps with their commands and outputs."""
        cfg = self.config
        lay = self.layout
        threads = str(cfg.threads)

        # Resolve effective inputs (provided vs produced by the pipeline).
        genomes_dir = cfg.genomes_dir or lay.genomes_dir
        kraken_output = cfg.kraken_output or lay.kraken_output
        reads_1 = cfg.extracted_reads_1 or lay.extracted_reads_1
        # When the user pre-supplies extracted reads, R2 is unknown to us; use
        # the pipeline's extracted R2 only when we run extraction ourselves.
        produced_reads_2 = lay.extracted_reads_2

        db = lay.database(cfg.aligner)
        alignment = lay.alignment(cfg.aligner)

        steps: list[_Step] = []

        # 1. Download references (only when no genomes_dir was supplied).
        download_argv = [
            *MDM_INVOCATION, "download", "genomes", "list", str(cfg.family or ""),
            "--output", str(lay.genomes_tsv), "--include-metadata",
        ]
        # Note: download is two subcommands; the orchestrator runs list then
        # fetch as a single logical step (see _run_download).
        steps.append(
            _Step(
                step=PipelineStep.DOWNLOAD,
                enabled=cfg.genomes_dir is None,
                argv=download_argv,
                outputs=[lay.genomes_dir],
                present=lambda: _is_present(lay.genomes_dir),
            )
        )

        # 2. Kraken2 classify (only when no .kraken provided and reads not pre-extracted).
        kc_argv = [
            *MDM_INVOCATION, "kraken2", "classify",
            "--reads-1", str(cfg.reads_1) if cfg.reads_1 else "",
            "--kraken-db", str(cfg.kraken_db) if cfg.kraken_db else "",
            "--output", str(lay.kraken_output), "--threads", threads,
        ]
        if cfg.reads_2:
            kc_argv += ["--reads-2", str(cfg.reads_2)]
        steps.append(
            _Step(
                step=PipelineStep.KRAKEN_CLASSIFY,
                enabled=cfg.extracted_reads_1 is None and cfg.kraken_output is None,
                argv=kc_argv,
                outputs=[lay.kraken_output],
                present=lambda: _is_present(lay.kraken_output),
            )
        )

        # 3. Kraken2 extract family reads.
        ke_argv = [
            *MDM_INVOCATION, "kraken2", "extract",
            "--kraken-output", str(kraken_output),
            "--reads-1", str(cfg.reads_1) if cfg.reads_1 else "",
            "--taxid", str(cfg.family_taxid) if cfg.family_taxid is not None else "",
            "--output", str(lay.extraction_dir),
        ]
        if cfg.reads_2:
            ke_argv += ["--reads-2", str(cfg.reads_2)]
        steps.append(
            _Step(
                step=PipelineStep.KRAKEN_EXTRACT,
                enabled=cfg.extracted_reads_1 is None,
                argv=ke_argv,
                outputs=[lay.extracted_reads_1],
                present=lambda: _is_present(lay.extracted_reads_1),
            )
        )

        # 4. Build alignment database.
        if cfg.aligner == Aligner.BLAST:
            makedb_argv = [
                *MDM_INVOCATION, "blast", "makedb",
                "--genomes", str(genomes_dir), "--output", str(db),
            ]
        else:
            makedb_argv = [
                *MDM_INVOCATION, "mmseqs2", "makedb",
                "--genomes", str(genomes_dir), "--output", str(db),
            ]
        steps.append(
            _Step(
                step=PipelineStep.MAKEDB,
                enabled=True,
                argv=makedb_argv,
                outputs=[db],
                present=lambda: _db_present(db),
            )
        )

        # 5. Align extracted reads.
        if cfg.aligner == Aligner.BLAST:
            align_argv = [
                *MDM_INVOCATION, "blast", "align",
                "--query", str(reads_1), "--database", str(db),
                "--output", str(alignment), "--threads", threads,
            ]
        else:
            align_argv = [
                *MDM_INVOCATION, "mmseqs2", "search",
                "--query-1", str(reads_1), "--database", str(db),
                "--output", str(alignment), "--threads", threads,
            ]
            if cfg.extracted_reads_1 is None and produced_reads_2 is not None:
                # paired-end available only when we extracted reads ourselves
                align_argv += ["--query-2", str(produced_reads_2)]
        steps.append(
            _Step(
                step=PipelineStep.ALIGN,
                enabled=True,
                argv=align_argv,
                outputs=[alignment],
                present=lambda: _is_present(alignment),
            )
        )

        # 6. Compute ANI matrix (use representative mapping when metadata present).
        ani_argv = [
            *MDM_INVOCATION, "ani", "compute",
            "--genomes", str(genomes_dir), "--output", str(lay.ani_matrix),
            "--threads", threads,
        ]
        if _is_present(lay.metadata_tsv):
            ani_argv += ["--metadata", str(lay.metadata_tsv), "--representatives-only"]
        steps.append(
            _Step(
                step=PipelineStep.ANI,
                enabled=True,
                argv=ani_argv,
                outputs=[lay.ani_matrix],
                present=lambda: _is_present(lay.ani_matrix),
            )
        )

        # 7. Classify.
        classify_argv = [
            *MDM_INVOCATION, "score", "classify",
            "--alignment", str(alignment), "--ani", str(lay.ani_matrix),
            "--output", str(lay.classifications), "--summary", str(lay.summary_json),
        ]
        if _is_present(lay.metadata_tsv):
            classify_argv += ["--metadata", str(lay.metadata_tsv)]
        steps.append(
            _Step(
                step=PipelineStep.CLASSIFY,
                enabled=True,
                argv=classify_argv,
                outputs=[lay.classifications],
                present=lambda: _is_present(lay.classifications),
            )
        )

        # 8. Report (optional).
        report_argv = [
            *MDM_INVOCATION, "report", "generate",
            "--classifications", str(lay.classifications), "--ani", str(lay.ani_matrix),
            "--output", str(lay.report), "--report-mode", "offline",
        ]
        if _is_present(lay.metadata_tsv):
            report_argv += ["--metadata", str(lay.metadata_tsv)]
        steps.append(
            _Step(
                step=PipelineStep.REPORT,
                enabled=cfg.report,
                argv=report_argv,
                outputs=[lay.report],
                present=lambda: _is_present(lay.report),
            )
        )

        return steps

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------
    def run(self, from_step: PipelineStep | None = None, *, mdm_version: str = "") -> RunManifest:
        """Execute the planned steps, honouring checkpoints and ``from_step``."""
        self.layout.root.mkdir(parents=True, exist_ok=True)
        plan = self.plan()

        order = PipelineStep.ordered()
        from_idx = order.index(from_step) if from_step is not None else 0
        from_label = from_step.value if from_step is not None else ""

        manifest = RunManifest(
            mdm_version=mdm_version,
            output_dir=str(self.layout.root),
            aligner=self.config.aligner.value,
            config={k: _jsonable(v) for k, v in self.config.model_dump().items()},
        )

        for step in plan:
            step_idx = order.index(step.step)

            if not step.enabled:
                manifest.steps.append(StepResult(step.step.value, StepStatus.SKIPPED,
                                                 detail="not applicable for provided inputs"))
                continue

            if step_idx < from_idx:
                # Resuming: earlier steps must already have their outputs.
                if not step.present():
                    raise PipelineStepError(
                        message=(
                            f"Cannot resume from '{from_label}': required output of "
                            f"earlier step '{step.step.value}' is missing."
                        ),
                        suggestion="Run without --from to produce the missing intermediates.",
                    )
                manifest.steps.append(StepResult(step.step.value, StepStatus.CACHED,
                                                 outputs=[str(p) for p in step.outputs]))
                continue

            command_str = " ".join(step.argv)

            if not self.force and step.present():
                self._log(f"[cached] {step.step.value}: outputs present, skipping")
                manifest.steps.append(StepResult(step.step.value, StepStatus.CACHED,
                                                 command=command_str,
                                                 outputs=[str(p) for p in step.outputs]))
                continue

            if self.dry_run:
                self._log(f"[dry-run] {step.step.value}: {command_str}")
                manifest.steps.append(StepResult(step.step.value, StepStatus.SKIPPED,
                                                 command=command_str, detail="dry-run"))
                continue

            self._log(f"[run] {step.step.value}")
            start = time.monotonic()
            code = self._run_step(step)
            elapsed = time.monotonic() - start

            if code != 0:
                manifest.steps.append(StepResult(step.step.value, StepStatus.FAILED,
                                                 command=command_str, elapsed_s=elapsed))
                self._write_manifest(manifest)
                raise PipelineStepError(
                    message=f"Pipeline step '{step.step.value}' failed (exit code {code}).",
                    suggestion=(
                        "Inspect the step output above. Re-run with --debug for tracebacks, "
                        "or re-run the step's command directly:\n  " + command_str
                    ),
                )

            manifest.steps.append(StepResult(step.step.value, StepStatus.COMPLETED,
                                             command=command_str,
                                             outputs=[str(p) for p in step.outputs],
                                             elapsed_s=round(elapsed, 2)))

        self._write_manifest(manifest)
        return manifest

    def _run_step(self, step: _Step) -> int:
        """Run a single step. DOWNLOAD is special (list + fetch)."""
        if step.step == PipelineStep.DOWNLOAD:
            code = self._runner(step.argv, self.dry_run)
            if code != 0:
                return code
            fetch_argv = [
                *MDM_INVOCATION, "download", "genomes", "fetch",
                "--accessions", str(self.layout.genomes_tsv),
                "--output-dir", str(self.layout.genomes_dir),
            ]
            return self._runner(fetch_argv, self.dry_run)
        return self._runner(step.argv, self.dry_run)

    def _write_manifest(self, manifest: RunManifest) -> None:
        if self.dry_run:
            return
        self.layout.manifest.parent.mkdir(parents=True, exist_ok=True)
        self.layout.manifest.write_text(manifest.to_json())


def _jsonable(value: object) -> object:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Enum):
        return value.value
    return value
