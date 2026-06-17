"""Tests for the mdm run pipeline orchestrator (no external tools required)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from metadarkmatter.core.pipeline import (
    Pipeline,
    PipelineLayout,
    PipelineStepError,
    StepStatus,
)
from metadarkmatter.models.pipeline_config import Aligner, PipelineConfig, PipelineStep


def _compute_only_config(tmp_path: Path) -> PipelineConfig:
    """A config that skips download + kraken (genomes + extracted reads provided)."""
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    (genomes / "GCF_1.fna").write_text(">c1\nACGT\n")
    reads = tmp_path / "reads_R1.fastq.gz"
    reads.write_bytes(b"fake")
    return PipelineConfig(
        genomes_dir=genomes,
        extracted_reads_1=reads,
        output_dir=tmp_path / "out",
    )


class RecordingRunner:
    """Records argv of each step; always succeeds."""

    def __init__(self) -> None:
        self.calls: list[list[str]] = []

    def __call__(self, argv, dry_run: bool) -> int:
        self.calls.append(list(argv))
        return 0


def test_validation_requires_genomes_or_family(tmp_path: Path):
    with pytest.raises(ValueError, match="genomes.*family|family.*genomes"):
        PipelineConfig(output_dir=tmp_path, extracted_reads_1=tmp_path / "r.fq")


def test_validation_requires_taxid_when_extracting(tmp_path: Path):
    genomes = tmp_path / "g"
    genomes.mkdir()
    with pytest.raises(ValueError, match="family-taxid"):
        PipelineConfig(
            genomes_dir=genomes,
            reads_1=tmp_path / "r.fq",
            kraken_db=tmp_path / "db",
            output_dir=tmp_path / "out",
        )


def test_plan_skips_inapplicable_steps(tmp_path: Path):
    cfg = _compute_only_config(tmp_path)
    pipeline = Pipeline(cfg, PipelineLayout(root=cfg.output_dir))
    plan = {s.step: s for s in pipeline.plan()}

    # download + both kraken steps are disabled (inputs pre-supplied)
    assert not plan[PipelineStep.DOWNLOAD].enabled
    assert not plan[PipelineStep.KRAKEN_CLASSIFY].enabled
    assert not plan[PipelineStep.KRAKEN_EXTRACT].enabled
    # compute steps enabled
    assert plan[PipelineStep.MAKEDB].enabled
    assert plan[PipelineStep.ALIGN].enabled
    assert plan[PipelineStep.ANI].enabled
    assert plan[PipelineStep.CLASSIFY].enabled
    assert plan[PipelineStep.REPORT].enabled


def test_dry_run_executes_nothing_and_marks_skipped(tmp_path: Path):
    cfg = _compute_only_config(tmp_path)
    runner = RecordingRunner()
    pipeline = Pipeline(cfg, PipelineLayout(root=cfg.output_dir), dry_run=True, runner=runner)
    manifest = pipeline.run()

    assert runner.calls == []  # nothing executed
    statuses = {s.name: s.status for s in manifest.steps}
    assert statuses["makedb"] == StepStatus.SKIPPED
    assert statuses["classify"] == StepStatus.SKIPPED
    # manifest is not written in dry-run
    assert not (cfg.output_dir / "run_manifest.json").exists()


def test_full_run_executes_enabled_steps_and_writes_manifest(tmp_path: Path):
    cfg = _compute_only_config(tmp_path)
    layout = PipelineLayout(root=cfg.output_dir)
    runner = RecordingRunner()
    pipeline = Pipeline(cfg, layout, runner=runner)
    manifest = pipeline.run(mdm_version="9.9.9")

    executed = list(runner.calls)
    # 5 compute steps run (makedb, align, ani, classify, report)
    assert len(executed) == 5
    # all reference the module entry point
    assert all("metadarkmatter.cli.main" in " ".join(c) for c in executed)

    completed = [s for s in manifest.steps if s.status == StepStatus.COMPLETED]
    assert {s.name for s in completed} == {"makedb", "align", "ani", "classify", "report"}

    # manifest persisted with provenance
    assert layout.manifest.exists()
    data = json.loads(layout.manifest.read_text())
    assert data["mdm_version"] == "9.9.9"
    assert data["aligner"] == "blast"
    assert len(data["steps"]) == 8  # all steps recorded (3 skipped + 5 completed)


def test_checkpoint_skips_present_outputs(tmp_path: Path):
    cfg = _compute_only_config(tmp_path)
    layout = PipelineLayout(root=cfg.output_dir)
    # Pre-create the classification output so CLASSIFY is cached.
    layout.classifications.parent.mkdir(parents=True, exist_ok=True)
    layout.classifications.write_text("read_id,taxonomic_call\nr1,Known Species\n")

    runner = RecordingRunner()
    pipeline = Pipeline(cfg, layout, runner=runner)
    manifest = pipeline.run()

    statuses = {s.name: s.status for s in manifest.steps}
    assert statuses["classify"] == StepStatus.CACHED
    # classify command must NOT have been executed
    assert not any("classify" in " ".join(c) for c in runner.calls)


def test_force_reruns_even_when_present(tmp_path: Path):
    cfg = _compute_only_config(tmp_path)
    layout = PipelineLayout(root=cfg.output_dir)
    layout.classifications.parent.mkdir(parents=True, exist_ok=True)
    layout.classifications.write_text("data")

    runner = RecordingRunner()
    pipeline = Pipeline(cfg, layout, force=True, runner=runner)
    pipeline.run()

    assert any("classify" in " ".join(c) for c in runner.calls)


def test_from_step_errors_when_prior_outputs_missing(tmp_path: Path):
    cfg = _compute_only_config(tmp_path)
    layout = PipelineLayout(root=cfg.output_dir)
    runner = RecordingRunner()
    pipeline = Pipeline(cfg, layout, runner=runner)

    with pytest.raises(PipelineStepError, match="resume"):
        pipeline.run(from_step=PipelineStep.CLASSIFY)


def test_failing_step_raises_and_writes_manifest(tmp_path: Path):
    cfg = _compute_only_config(tmp_path)
    layout = PipelineLayout(root=cfg.output_dir)

    def failing_runner(argv, dry_run: bool) -> int:
        return 1

    pipeline = Pipeline(cfg, layout, runner=failing_runner)
    with pytest.raises(PipelineStepError, match="failed"):
        pipeline.run()

    # manifest records the failure
    data = json.loads(layout.manifest.read_text())
    assert any(s["status"] == "failed" for s in data["steps"])


def test_cli_dry_run_renders_plan(tmp_path: Path):
    """The `mdm run pipeline --dry-run` command runs without external tools."""
    from typer.testing import CliRunner

    from metadarkmatter.cli.main import app

    genomes = tmp_path / "genomes"
    genomes.mkdir()
    (genomes / "GCF_1.fna").write_text(">c1\nACGT\n")
    reads = tmp_path / "reads_R1.fastq.gz"
    reads.write_bytes(b"x")

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "run", "pipeline",
            "--genomes", str(genomes),
            "--extracted-reads-1", str(reads),
            "--output-dir", str(tmp_path / "out"),
            "--dry-run",
        ],
    )
    assert result.exit_code == 0, result.output
    assert "DRY RUN" in result.output
    assert "Pipeline steps" in result.output


def test_cli_invalid_inputs_exit_friendly(tmp_path: Path):
    """Missing required inputs surface as a friendly error (via the CLI wrapper)."""
    from typer.testing import CliRunner

    from metadarkmatter.cli.main import app

    runner = CliRunner()
    result = runner.invoke(app, ["run", "pipeline", "--output-dir", str(tmp_path / "out")])
    assert result.exit_code == 1
    assert "Invalid" in result.output or "Provide" in result.output


def test_mmseqs2_aligner_changes_paths(tmp_path: Path):
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    (genomes / "GCF_1.fna").write_text(">c1\nACGT\n")
    reads = tmp_path / "reads_R1.fastq.gz"
    reads.write_bytes(b"x")
    cfg = PipelineConfig(
        genomes_dir=genomes,
        extracted_reads_1=reads,
        output_dir=tmp_path / "out",
        aligner=Aligner.MMSEQS2,
    )
    layout = PipelineLayout(root=cfg.output_dir)
    plan = {s.step: s for s in Pipeline(cfg, layout).plan()}
    align_cmd = " ".join(plan[PipelineStep.ALIGN].argv)
    assert "mmseqs2" in align_cmd
    assert str(layout.alignment(Aligner.MMSEQS2)) in align_cmd
