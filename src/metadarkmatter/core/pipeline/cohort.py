"""
Multi-sample cohort runs for ``mdm run cohort``.

A cohort shares one reference (downloaded genomes, alignment database, ANI
matrix) across many samples. The reference is built once under
``<output>/reference/``; each sample then runs its own per-sample steps
(Kraken2 classify/extract, align, classify, report) under
``<output>/samples/<sample>/``. Both stages reuse the single-sample
:class:`~metadarkmatter.core.pipeline.orchestrator.Pipeline` (and therefore its
fingerprint checkpointing), restricted to the relevant steps.

Samples run sequentially in input order so results are deterministic. A
failing sample does not stop the others; it is recorded as failed in the
cohort summary and the run as a whole reports failure.
"""

from __future__ import annotations

import csv
import json
import re
from collections.abc import Callable
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

import polars as pl

from metadarkmatter.core.exceptions import MetadarkmatterError
from metadarkmatter.core.pipeline.orchestrator import (
    MDM_INVOCATION,
    Pipeline,
    PipelineLayout,
    PipelineStepError,
    RunManifest,
    StepRunner,
    _default_runner,
)
from metadarkmatter.models.pipeline_config import Aligner, PipelineConfig, PipelineStep

# Steps that build the shared reference vs. the steps repeated per sample.
REFERENCE_STEPS: frozenset[PipelineStep] = frozenset(
    {PipelineStep.DOWNLOAD, PipelineStep.MAKEDB, PipelineStep.ANI}
)
SAMPLE_STEPS: frozenset[PipelineStep] = frozenset(
    set(PipelineStep.ordered()) - REFERENCE_STEPS
)

NOVEL_CATEGORIES = ("Novel Species", "Novel Genus")

_SAMPLE_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
# ``report multi`` derives sample names from file names by deleting these
# substrings anywhere in the stem, so a sample name containing one would be
# silently renamed in the multi-sample report.
_RESERVED_NAME_PARTS = ("_classifications", "classifications_", ".kraken", ".blast", ".tsv")
_SAMPLE_COLUMNS = ("sample", "reads_1", "reads_2", "extracted_reads_1", "extracted_reads_2")


class CohortError(MetadarkmatterError):
    """Raised for invalid cohort input or when the shared reference cannot be built."""


@dataclass(frozen=True)
class CohortSample:
    """One row of the cohort sample sheet."""

    name: str
    reads_1: Path | None = None
    reads_2: Path | None = None
    extracted_reads_1: Path | None = None
    extracted_reads_2: Path | None = None


def read_samples_tsv(path: Path) -> list[CohortSample]:
    """Parse and validate a cohort sample sheet.

    Columns: ``sample`` plus either ``reads_1`` (+ optional ``reads_2``) for
    raw reads or ``extracted_reads_1`` (+ optional ``extracted_reads_2``) for
    reads already extracted for the target family (skips Kraken2). Relative
    paths resolve against the sheet's directory. Blank lines and lines starting
    with ``#`` are ignored.
    """
    if not path.is_file():
        raise CohortError(message=f"Sample sheet not found: {path}")

    with path.open(newline="", encoding="utf-8-sig") as fh:
        lines = [ln for ln in fh if ln.strip() and not ln.lstrip().startswith("#")]
    reader = csv.DictReader(lines, delimiter="\t")
    header = [h.strip() for h in (reader.fieldnames or [])]
    if "sample" not in header:
        raise CohortError(
            message=f"Sample sheet {path} needs a tab-separated header with a 'sample' column.",
            suggestion=f"Expected columns: {', '.join(_SAMPLE_COLUMNS)}.",
        )
    unknown = [h for h in header if h not in _SAMPLE_COLUMNS]
    if unknown:
        raise CohortError(
            message=f"Unknown column(s) in sample sheet {path}: {', '.join(unknown)}.",
            suggestion=f"Expected columns: {', '.join(_SAMPLE_COLUMNS)}.",
        )

    base = path.parent
    samples: list[CohortSample] = []
    seen: set[str] = set()
    problems: list[str] = []

    def _path(row: dict[str, str | None], key: str) -> Path | None:
        value = (row.get(key) or "").strip()
        if not value:
            return None
        p = Path(value)
        return p if p.is_absolute() else base / p

    for lineno, raw in enumerate(reader, start=2):
        row = {(k or "").strip(): v for k, v in raw.items()}
        name = (row.get("sample") or "").strip()
        label = f"row {lineno} ({name or 'no sample name'})"
        if not name:
            problems.append(f"{label}: empty sample name")
            continue
        if not _SAMPLE_NAME_RE.match(name):
            problems.append(
                f"{label}: sample names may only contain letters, digits, '_', '.' and '-' "
                "and must start with a letter or digit"
            )
            continue
        reserved = [part for part in _RESERVED_NAME_PARTS if part in name]
        if reserved:
            problems.append(
                f"{label}: sample name must not contain {reserved[0]!r} "
                "(it would be altered in the multi-sample report)"
            )
            continue
        if name in seen:
            problems.append(f"{label}: duplicate sample name '{name}'")
            continue
        seen.add(name)

        sample = CohortSample(
            name=name,
            reads_1=_path(row, "reads_1"),
            reads_2=_path(row, "reads_2"),
            extracted_reads_1=_path(row, "extracted_reads_1"),
            extracted_reads_2=_path(row, "extracted_reads_2"),
        )
        if sample.reads_1 is None and sample.extracted_reads_1 is None:
            problems.append(f"{label}: provide reads_1 or extracted_reads_1")
        if sample.reads_1 is not None and sample.extracted_reads_1 is not None:
            problems.append(f"{label}: provide reads_1 or extracted_reads_1, not both")
        if sample.reads_2 is not None and sample.reads_1 is None:
            problems.append(f"{label}: reads_2 requires reads_1")
        if sample.extracted_reads_2 is not None and sample.extracted_reads_1 is None:
            problems.append(f"{label}: extracted_reads_2 requires extracted_reads_1")
        for column in ("reads_1", "reads_2", "extracted_reads_1", "extracted_reads_2"):
            file = getattr(sample, column)
            if file is not None and not file.is_file():
                problems.append(f"{label}: {column} file not found: {file}")
        samples.append(sample)

    if not samples and not problems:
        problems.append("no samples listed")
    if problems:
        raise CohortError(
            message=f"Invalid sample sheet {path}:\n  - " + "\n  - ".join(problems),
        )
    return samples


@dataclass
class CohortLayout:
    """File locations for a cohort run."""

    root: Path

    @property
    def reference_dir(self) -> Path:
        return self.root / "reference"

    @property
    def samples_dir(self) -> Path:
        return self.root / "samples"

    def sample_dir(self, name: str) -> Path:
        return self.samples_dir / name

    @property
    def multi_inputs_dir(self) -> Path:
        """Per-sample classification links consumed by ``report multi``."""
        return self.root / "multi_inputs"

    @property
    def summary_tsv(self) -> Path:
        return self.root / "cohort_summary.tsv"

    @property
    def report(self) -> Path:
        return self.root / "cohort_report.html"

    @property
    def manifest(self) -> Path:
        return self.root / "cohort_manifest.json"


@dataclass
class SampleOutcome:
    name: str
    status: str  # "completed", "failed", or "planned" (dry run)
    manifest: RunManifest | None = None
    error: str | None = None
    classifications: Path | None = None


@dataclass
class CohortResult:
    reference: RunManifest | None = None
    samples: list[SampleOutcome] = field(default_factory=list)
    summary_path: Path | None = None
    report_path: Path | None = None
    report_error: str | None = None

    @property
    def failed(self) -> list[SampleOutcome]:
        return [s for s in self.samples if s.status == "failed"]

    @property
    def ok(self) -> bool:
        return not self.failed and self.report_error is None


def summarize_classifications(path: Path) -> tuple[int, dict[str, int]]:
    """Return ``(total_reads, {taxonomic_call: count})`` for a classifications file."""
    df = pl.read_csv(path, columns=["taxonomic_call"])
    calls = df.get_column("taxonomic_call").fill_null("Unknown")
    counts: dict[str, int] = {}
    for call in calls:
        counts[str(call)] = counts.get(str(call), 0) + 1
    return df.height, counts


def write_cohort_summary(outcomes: list[SampleOutcome], path: Path) -> None:
    """Write ``cohort_summary.tsv``: one row per sample, categories sorted by name.

    A completed sample whose classifications cannot be read is marked failed.
    """
    rows: list[tuple[SampleOutcome, int, dict[str, int]]] = []
    for outcome in outcomes:
        total = 0
        counts: dict[str, int] = {}
        if outcome.status == "completed" and outcome.classifications is not None:
            try:
                total, counts = summarize_classifications(outcome.classifications)
            except (OSError, pl.exceptions.PolarsError) as exc:
                outcome.status = "failed"
                outcome.error = f"could not summarize classifications: {exc}"
        rows.append((outcome, total, counts))

    categories = sorted({c for _, _, counts in rows for c in counts})
    header = ["sample", "status", "total_reads", *categories, "novel_fraction", "error"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for outcome, total, counts in rows:
            if outcome.status != "completed":
                writer.writerow(
                    [outcome.name, outcome.status, "", *[""] * len(categories), "",
                     " ".join((outcome.error or "").split())]
                )
                continue
            novel = sum(counts.get(c, 0) for c in NOVEL_CATEGORIES)
            fraction = f"{novel / total:.6f}" if total else "0.000000"
            writer.writerow(
                [outcome.name, outcome.status, total, *[counts.get(c, 0) for c in categories],
                 fraction, ""]
            )


class Cohort:
    """Plan and execute a cohort run over a shared reference."""

    def __init__(
        self,
        samples: list[CohortSample],
        shared: dict[str, object],
        layout: CohortLayout,
        *,
        aligner: Aligner = Aligner.BLAST,
        dry_run: bool = False,
        force: bool = False,
        runner: StepRunner = _default_runner,
        log: Callable[[str], None] | None = None,
        sample_reports: bool = True,
        cohort_report: bool = True,
    ) -> None:
        if not samples:
            raise CohortError(message="A cohort needs at least one sample.")
        self.samples = samples
        self.layout = layout
        self.dry_run = dry_run
        self.force = force
        self._runner = runner
        self._log = log or (lambda _msg: None)
        self.cohort_report = cohort_report
        self.aligner = aligner

        base = {**shared, "aligner": aligner, "report": sample_reports}
        self._configs: dict[str, PipelineConfig] = {}
        for sample in samples:
            try:
                self._configs[sample.name] = PipelineConfig(
                    **base,
                    reads_1=sample.reads_1,
                    reads_2=sample.reads_2,
                    extracted_reads_1=sample.extracted_reads_1,
                    output_dir=layout.sample_dir(sample.name),
                )
            except ValueError as exc:
                raise CohortError(message=f"Sample '{sample.name}': {exc}") from exc

    def reference_pipeline(self) -> Pipeline:
        """Pipeline for the shared reference (download, makedb, ANI)."""
        cfg = self._configs[self.samples[0].name].model_copy(
            update={"output_dir": self.layout.reference_dir, "report": False}
        )
        layout = PipelineLayout(root=self.layout.reference_dir, family_taxid=cfg.family_taxid)
        return Pipeline(
            cfg, layout, dry_run=self.dry_run, force=self.force, runner=self._runner,
            log=lambda m: self._log(f"[reference] {m}"), only_steps=REFERENCE_STEPS,
        )

    def sample_pipeline(self, sample: CohortSample) -> Pipeline:
        """Pipeline for one sample's per-sample steps."""
        cfg = self._configs[sample.name]
        layout = PipelineLayout(
            root=self.layout.sample_dir(sample.name),
            sample_name=sample.name,
            family_taxid=cfg.family_taxid,
            reference_root=self.layout.reference_dir,
        )
        name = sample.name
        return Pipeline(
            cfg, layout, dry_run=self.dry_run, force=self.force, runner=self._runner,
            log=lambda m: self._log(f"[{name}] {m}"),
            extracted_reads_2=sample.extracted_reads_2, only_steps=SAMPLE_STEPS,
        )

    def run(self, from_step: PipelineStep | None = None, *, mdm_version: str = "") -> CohortResult:
        result = CohortResult()
        started = _now_iso()

        self._log("== shared reference ==")
        try:
            result.reference = self.reference_pipeline().run(
                from_step=from_step, mdm_version=mdm_version
            )
        except PipelineStepError as exc:
            raise CohortError(
                message=f"Shared reference could not be built: {exc.message}",
                suggestion=exc.suggestion,
            ) from exc

        for sample in self.samples:
            self._log(f"== sample {sample.name} ==")
            pipeline = self.sample_pipeline(sample)
            try:
                manifest = pipeline.run(from_step=from_step, mdm_version=mdm_version)
            except PipelineStepError as exc:
                result.samples.append(SampleOutcome(sample.name, "failed", error=exc.message))
                self._log(f"[{sample.name}] FAILED: {exc.message}")
                continue
            result.samples.append(
                SampleOutcome(
                    sample.name,
                    "planned" if self.dry_run else "completed",
                    manifest=manifest,
                    classifications=pipeline.layout.classifications,
                )
            )

        if self.dry_run:
            self._log("== cohort summary and report (planned) ==")
            self._log(f"[dry-run] write {self.layout.summary_tsv}")
            if self.cohort_report:
                self._log(f"[dry-run] {' '.join(self._multi_argv())}")
            return result

        write_cohort_summary(result.samples, self.layout.summary_tsv)
        result.summary_path = self.layout.summary_tsv
        if self.cohort_report:
            self._build_cohort_report(result)
        self._write_manifest(result, started, mdm_version)
        return result

    def _multi_argv(self) -> list[str]:
        return [
            *MDM_INVOCATION, "report", "multi",
            "--input-dir", str(self.layout.multi_inputs_dir),
            "--output", str(self.layout.report),
            "--title", "Cohort comparison",
        ]

    def _build_cohort_report(self, result: CohortResult) -> None:
        ok = [s for s in result.samples if s.status == "completed" and s.classifications]
        if not ok:
            result.report_error = "no sample completed; cohort report not generated"
            return
        inputs = self.layout.multi_inputs_dir
        inputs.mkdir(parents=True, exist_ok=True)
        # Rebuild the link set so samples that failed this time do not linger.
        for stale in inputs.glob("*_classifications.csv"):
            stale.unlink()
        for outcome in ok:
            assert outcome.classifications is not None
            link = inputs / f"{outcome.name}_classifications.csv"
            link.symlink_to(Path("..") / outcome.classifications.relative_to(self.layout.root))
        self._log("[cohort] report multi")
        code = self._runner(self._multi_argv(), False)
        if code != 0:
            result.report_error = f"'report multi' failed (exit code {code})"
        else:
            result.report_path = self.layout.report

    def _write_manifest(self, result: CohortResult, started: str, mdm_version: str) -> None:
        payload = {
            "mdm_version": mdm_version,
            "output_dir": str(self.layout.root),
            "aligner": self.aligner.value,
            "started_at": started,
            "finished_at": _now_iso(),
            "samples": [
                {
                    "name": s.name,
                    "status": s.status,
                    "error": s.error,
                    "directory": str(self.layout.sample_dir(s.name)),
                }
                for s in result.samples
            ],
            "summary": str(self.layout.summary_tsv),
            "report": str(result.report_path) if result.report_path else None,
            "report_error": result.report_error,
        }
        self.layout.root.mkdir(parents=True, exist_ok=True)
        self.layout.manifest.write_text(json.dumps(payload, indent=2))


def _now_iso() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds")
