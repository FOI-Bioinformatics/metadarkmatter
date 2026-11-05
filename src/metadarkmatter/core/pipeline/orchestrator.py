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

import hashlib
import json
import subprocess
import sys
import time
from collections.abc import Callable, Sequence
from dataclasses import asdict, dataclass, field, replace
from datetime import UTC, datetime
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
    fingerprint: str | None = None

    def to_dict(self) -> dict[str, object]:
        d = asdict(self)
        d["status"] = self.status.value
        return d


@dataclass
class PipelineLayout:
    """Single source of truth for all pipeline file paths under one root.

    ``sample_name`` and ``family_taxid`` determine the Kraken2 output and
    extracted-read filenames, which must match what ``mdm kraken2 classify``
    (``{sample}.kraken`` / ``{sample}.kreport`` in its ``--output`` dir) and
    ``mdm kraken2 extract`` (``{sample}_taxid{taxid}_R1.fastq.gz``) write.
    """

    root: Path
    sample_name: str = "sample"
    family_taxid: int | None = None
    # Cohort runs keep the shared reference (genomes, alignment database, ANI
    # matrix) in one directory and per-sample files under ``root``. Unset, the
    # reference lives under ``root`` (single-sample layout).
    reference_root: Path | None = None

    @property
    def ref_root(self) -> Path:
        return self.reference_root if self.reference_root is not None else self.root

    @property
    def genomes_dir(self) -> Path:
        return self.ref_root / "genomes"

    @property
    def genomes_tsv(self) -> Path:
        return self.ref_root / "genomes.tsv"

    @property
    def metadata_tsv(self) -> Path:
        return self.ref_root / "genome_metadata.tsv"

    @property
    def kraken_dir(self) -> Path:
        """Directory passed as ``--output`` to ``kraken2 classify``."""
        return self.root

    @property
    def kraken_output(self) -> Path:
        return self.kraken_dir / f"{self.sample_name}.kraken"

    @property
    def kraken_report(self) -> Path:
        return self.kraken_dir / f"{self.sample_name}.kreport"

    @property
    def extraction_dir(self) -> Path:
        return self.root / "extraction"

    def _extracted_reads(self, mate: int) -> Path:
        taxid = f"_taxid{self.family_taxid}" if self.family_taxid is not None else ""
        return self.extraction_dir / f"{self.sample_name}{taxid}_R{mate}.fastq.gz"

    @property
    def extracted_reads_1(self) -> Path:
        return self._extracted_reads(1)

    @property
    def extracted_reads_2(self) -> Path:
        return self._extracted_reads(2)

    def database(self, aligner: Aligner) -> Path:
        sub = "blastdb" if aligner == Aligner.BLAST else "mmseqs_db"
        return self.ref_root / sub / "pangenome"

    def alignment(self, aligner: Aligner) -> Path:
        return self.root / f"{self.sample_name}.{aligner.value}.tsv.gz"

    @property
    def ani_matrix(self) -> Path:
        return self.ref_root / "ani_matrix.csv"

    @property
    def aai_matrix(self) -> Path:
        return self.ref_root / "aai_matrix.csv"

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

    @property
    def state_dir(self) -> Path:
        """Directory holding per-step completion markers."""
        return self.root / ".mdm_steps"

    def marker(self, step: PipelineStep) -> Path:
        return self.state_dir / f"{step.value}.done.json"

    def running_marker(self, step: PipelineStep) -> Path:
        """Present while a step is executing (or was killed mid-run)."""
        return self.state_dir / f"{step.value}.running"


@dataclass
class RunManifest:
    mdm_version: str
    output_dir: str
    aligner: str
    config: dict[str, object]
    steps: list[StepResult] = field(default_factory=list)
    started_at: str | None = None
    finished_at: str | None = None
    warnings: list[str] = field(default_factory=list)

    def to_json(self) -> str:
        return json.dumps(
            {
                "mdm_version": self.mdm_version,
                "output_dir": self.output_dir,
                "aligner": self.aligner,
                "started_at": self.started_at,
                "finished_at": self.finished_at,
                "config": self.config,
                "warnings": self.warnings,
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


def _db_present(db_prefix: Path, aligner: Aligner = Aligner.BLAST) -> bool:
    """A makedb output is present only when the finished database files exist.

    Intermediate files that share the prefix (e.g. the concatenated
    ``pangenome.fasta`` or ``pangenome_id_mapping.tsv`` written before the
    tool runs) do not count. BLAST: ``<prefix>.nhr`` / ``<prefix>.nal`` or a
    first volume ``<prefix>.00.nhr``. MMseqs2: ``<prefix>.dbtype``.
    """
    suffixes: tuple[str, ...] = (
        (".nhr", ".nal", ".00.nhr") if aligner == Aligner.BLAST else (".dbtype",)
    )
    for suffix in suffixes:
        candidate = db_prefix.parent / f"{db_prefix.name}{suffix}"
        if candidate.is_file() and candidate.stat().st_size > 0:
            return True
    return False


def _path_signature(path: Path) -> str:
    """Cheap change detector for an input: size + mtime (files) or a digest of both (dirs)."""
    try:
        if path.is_dir():
            h = hashlib.sha256()
            for p in sorted(path.rglob("*")):
                if p.is_file():
                    st = p.stat()
                    h.update(f"{p.relative_to(path)}:{st.st_size}:{st.st_mtime_ns}\n".encode())
            return "dir:" + h.hexdigest()
        st = path.stat()
        return f"{st.st_size}:{st.st_mtime_ns}"
    except OSError:
        return "missing"


@dataclass
class _Step:
    step: PipelineStep
    enabled: bool
    argv: list[str]
    outputs: list[Path]
    present: Callable[[], bool]
    # Input files/dirs whose size+mtime feed the step fingerprint.
    inputs: list[Path] = field(default_factory=list)
    # True when the step belongs to another pipeline of the same cohort run.
    out_of_scope: bool = False


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
        extracted_reads_2: Path | None = None,
        only_steps: frozenset[PipelineStep] | None = None,
    ) -> None:
        self.config = config
        # ``only_steps`` restricts the run to a subset of steps (cohort runs
        # split the shared-reference steps from the per-sample steps).
        self.only_steps = only_steps
        # Extracted-read filenames embed the taxid; derive it from the config
        # when the caller did not set it on the layout.
        if layout.family_taxid is None and config.family_taxid is not None:
            layout = replace(layout, family_taxid=config.family_taxid)
        self.layout = layout
        self.extracted_reads_2 = extracted_reads_2
        self.dry_run = dry_run
        self.force = force
        self._runner = runner
        self._log = log or (lambda _msg: None)
        self.warnings: list[str] = self._collect_warnings()

    def _collect_warnings(self) -> list[str]:
        cfg = self.config
        warnings: list[str] = []
        if cfg.compute_aai:
            warnings.append(
                "compute_aai is set but AAI computation is not implemented in "
                "'run pipeline'; it is ignored."
            )
        return warnings

    # ------------------------------------------------------------------
    # Completion markers / fingerprints
    # ------------------------------------------------------------------
    def _fingerprint(self, step: _Step) -> str:
        """Hash of the step's argv plus size+mtime of its input files."""
        payload = {
            "argv": step.argv,
            "inputs": {str(p): _path_signature(p) for p in step.inputs},
        }
        return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()

    def _read_marker(self, step: PipelineStep) -> dict[str, object] | None:
        try:
            data = json.loads(self.layout.marker(step).read_text())
        except (OSError, ValueError):
            return None
        return data if isinstance(data, dict) else None

    def _write_marker(self, step: _Step, fingerprint: str) -> None:
        marker = self.layout.marker(step.step)
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.write_text(
            json.dumps(
                {
                    "step": step.step.value,
                    "fingerprint": fingerprint,
                    "completed_at": _now_iso(),
                    "command": " ".join(step.argv),
                },
                indent=2,
            )
        )

    def _clear_marker(self, step: PipelineStep) -> None:
        self.layout.marker(step).unlink(missing_ok=True)

    def _mark_running(self, step: PipelineStep) -> None:
        running = self.layout.running_marker(step)
        running.parent.mkdir(parents=True, exist_ok=True)
        running.write_text(_now_iso())

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
        # Paired-end alignment input: the user-supplied extracted R2 when the
        # reads were pre-extracted, otherwise our own extracted R2 only when
        # the raw input was paired.
        reads_2: Path | None
        if cfg.extracted_reads_1 is not None:
            reads_2 = self.extracted_reads_2
        else:
            reads_2 = lay.extracted_reads_2 if cfg.reads_2 is not None else None

        # download writes genome_metadata.tsv next to genomes.tsv, so it will
        # exist by the time later steps run even though it is absent now.
        download_enabled = cfg.genomes_dir is None
        use_metadata = download_enabled or _is_present(lay.metadata_tsv)

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
                enabled=download_enabled,
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
            "--output", str(lay.kraken_dir), "--sample-name", lay.sample_name,
            "--threads", threads,
        ]
        if cfg.reads_2:
            kc_argv += ["--reads-2", str(cfg.reads_2)]
        steps.append(
            _Step(
                step=PipelineStep.KRAKEN_CLASSIFY,
                enabled=cfg.extracted_reads_1 is None and cfg.kraken_output is None,
                argv=kc_argv,
                outputs=[lay.kraken_output, lay.kraken_report],
                present=lambda: _is_present(lay.kraken_output) and _is_present(lay.kraken_report),
                inputs=[p for p in (cfg.reads_1, cfg.reads_2, cfg.kraken_db) if p is not None],
            )
        )

        # 3. Kraken2 extract family reads.
        ke_argv = [
            *MDM_INVOCATION, "kraken2", "extract",
            "--kraken-output", str(kraken_output),
            "--reads-1", str(cfg.reads_1) if cfg.reads_1 else "",
            "--taxid", str(cfg.family_taxid) if cfg.family_taxid is not None else "",
            "--output", str(lay.extraction_dir), "--sample-name", lay.sample_name,
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
                inputs=[
                    p for p in (kraken_output, cfg.reads_1, cfg.reads_2) if p is not None
                ],
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
                present=lambda: _db_present(db, cfg.aligner),
                inputs=[genomes_dir],
            )
        )

        # 5. Align extracted reads.
        if cfg.aligner == Aligner.BLAST:
            align_tool = "blast"
            align_sub = "align"
        else:
            align_tool = "mmseqs2"
            align_sub = "search"
        # Both aligners align both mates when R2 is available; mate IDs get
        # /1 and /2 suffixes so the classifier treats them as separate reads.
        query_args = (
            ["--query-1", str(reads_1), "--query-2", str(reads_2)]
            if reads_2 is not None
            else ["--query", str(reads_1)]
        )
        align_argv = [
            *MDM_INVOCATION, align_tool, align_sub,
            *query_args, "--database", str(db),
            "--output", str(alignment), "--threads", threads,
        ]
        steps.append(
            _Step(
                step=PipelineStep.ALIGN,
                enabled=True,
                argv=align_argv,
                outputs=[alignment],
                present=lambda: _is_present(alignment),
                inputs=[
                    reads_1,
                    *([reads_2] if reads_2 is not None else []),
                    db.parent,
                ],
            )
        )

        # 6. Compute ANI matrix (use representative mapping when metadata present).
        ani_argv = [
            *MDM_INVOCATION, "ani", "compute",
            "--genomes", str(genomes_dir), "--output", str(lay.ani_matrix),
            "--threads", threads,
        ]
        if use_metadata:
            ani_argv += ["--metadata", str(lay.metadata_tsv), "--representatives-only"]
        steps.append(
            _Step(
                step=PipelineStep.ANI,
                enabled=True,
                argv=ani_argv,
                outputs=[lay.ani_matrix],
                present=lambda: _is_present(lay.ani_matrix),
                inputs=[genomes_dir, *([lay.metadata_tsv] if use_metadata else [])],
            )
        )

        # 7. Classify.
        classify_argv = [
            *MDM_INVOCATION, "score", "classify",
            "--alignment", str(alignment), "--ani", str(lay.ani_matrix),
            "--output", str(lay.classifications), "--summary", str(lay.summary_json),
        ]
        if use_metadata:
            classify_argv += ["--metadata", str(lay.metadata_tsv)]
        steps.append(
            _Step(
                step=PipelineStep.CLASSIFY,
                enabled=True,
                argv=classify_argv,
                outputs=[lay.classifications],
                present=lambda: _is_present(lay.classifications),
                inputs=[alignment, lay.ani_matrix, *([lay.metadata_tsv] if use_metadata else [])],
            )
        )

        # 8. Report (optional).
        report_argv = [
            *MDM_INVOCATION, "report", "generate",
            "--classifications", str(lay.classifications), "--ani", str(lay.ani_matrix),
            "--output", str(lay.report), "--report-mode", "offline",
        ]
        if use_metadata:
            report_argv += ["--metadata", str(lay.metadata_tsv)]
        steps.append(
            _Step(
                step=PipelineStep.REPORT,
                enabled=cfg.report,
                argv=report_argv,
                outputs=[lay.report],
                present=lambda: _is_present(lay.report),
                inputs=[
                    lay.classifications,
                    lay.ani_matrix,
                    *([lay.metadata_tsv] if use_metadata else []),
                ],
            )
        )

        if self.only_steps is not None:
            for st in steps:
                if st.step not in self.only_steps:
                    st.enabled = False
                    st.out_of_scope = True

        return steps

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------
    def run(self, from_step: PipelineStep | None = None, *, mdm_version: str = "") -> RunManifest:
        """Execute the planned steps, honouring checkpoints and ``from_step``."""
        if not self.dry_run:
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
            started_at=_now_iso(),
            warnings=list(self.warnings),
        )

        # Once a step runs (or would run), every later step must re-run too.
        upstream_ran = False

        for step in plan:
            step_idx = order.index(step.step)

            if not step.enabled:
                detail = (
                    "handled by another stage of this run"
                    if step.out_of_scope
                    else "not applicable for provided inputs"
                )
                manifest.steps.append(StepResult(step.step.value, StepStatus.SKIPPED,
                                                 detail=detail))
                continue

            marker = self._read_marker(step.step)

            if step_idx < from_idx:
                # Resuming: earlier steps must have outputs and must not have
                # been interrupted mid-run. Outputs from runs that predate the
                # completion markers (no marker at all) are trusted here.
                if not step.present() or self.layout.running_marker(step.step).exists():
                    raise PipelineStepError(
                        message=(
                            f"Cannot resume from '{from_label}': required output of "
                            f"earlier step '{step.step.value}' is missing or was not "
                            "completed."
                        ),
                        suggestion="Run without --from to produce the missing intermediates.",
                    )
                fp_prev = marker.get("fingerprint") if marker else None
                manifest.steps.append(StepResult(
                    step.step.value, StepStatus.CACHED,
                    outputs=[str(p) for p in step.outputs],
                    fingerprint=fp_prev if isinstance(fp_prev, str) else None))
                continue

            command_str = " ".join(step.argv)
            fingerprint = self._fingerprint(step)

            # --from forces the selected step and every later one to re-run.
            forced = self.force or from_step is not None
            reason: str | None = None
            if forced:
                reason = "forced"
            elif upstream_ran:
                reason = "an upstream step was re-run"
            elif not step.present():
                reason = "outputs missing"
            elif marker is None:
                reason = "outputs have no completed-step marker (partial or pre-existing)"
            elif marker.get("fingerprint") != fingerprint:
                reason = "inputs or parameters changed since the last completed run"

            if reason is None:
                self._log(f"[cached] {step.step.value}: completed with identical inputs, skipping")
                manifest.steps.append(StepResult(step.step.value, StepStatus.CACHED,
                                                 command=command_str,
                                                 outputs=[str(p) for p in step.outputs],
                                                 fingerprint=fingerprint))
                continue

            upstream_ran = True

            if self.dry_run:
                self._log(f"[dry-run] {step.step.value}: {command_str}")
                manifest.steps.append(StepResult(step.step.value, StepStatus.SKIPPED,
                                                 command=command_str,
                                                 detail=f"dry-run ({reason})",
                                                 fingerprint=fingerprint))
                continue

            self._log(f"[run] {step.step.value} ({reason})")
            # Drop any stale marker first: a killed step must not look complete.
            self._clear_marker(step.step)
            self._mark_running(step.step)
            start = time.monotonic()
            code = self._run_step(step)
            elapsed = time.monotonic() - start

            if code != 0:
                manifest.steps.append(StepResult(step.step.value, StepStatus.FAILED,
                                                 command=command_str, elapsed_s=elapsed,
                                                 fingerprint=fingerprint))
                manifest.finished_at = _now_iso()
                self._write_manifest(manifest)
                raise PipelineStepError(
                    message=f"Pipeline step '{step.step.value}' failed (exit code {code}).",
                    suggestion=(
                        "Inspect the step output above. Re-run with --debug for tracebacks, "
                        "or re-run the step's command directly:\n  " + command_str
                    ),
                )

            # The step exited 0: only now is it marked complete. Re-fingerprint
            # so inputs produced during the step (e.g. downloads) are reflected.
            self._write_marker(step, self._fingerprint(step))
            self.layout.running_marker(step.step).unlink(missing_ok=True)
            manifest.steps.append(StepResult(step.step.value, StepStatus.COMPLETED,
                                             command=command_str,
                                             outputs=[str(p) for p in step.outputs],
                                             elapsed_s=round(elapsed, 2),
                                             detail=reason,
                                             fingerprint=fingerprint))

        manifest.finished_at = _now_iso()
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
        path = self.layout.manifest
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists():
            # Keep the previous run's manifest instead of overwriting it.
            stamp = datetime.fromtimestamp(path.stat().st_mtime, UTC).strftime("%Y%m%dT%H%M%SZ")
            archive = path.with_name(f"run_manifest.{stamp}.json")
            n = 1
            while archive.exists():
                archive = path.with_name(f"run_manifest.{stamp}.{n}.json")
                n += 1
            path.rename(archive)
        path.write_text(manifest.to_json())


def _now_iso() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds")


def _jsonable(value: object) -> object:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Enum):
        return value.value
    return value
