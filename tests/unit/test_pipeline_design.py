"""Pipeline: --metadata with --genomes, dependency-based caching, tool versions."""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.main import app
from metadarkmatter.core.pipeline.orchestrator import Pipeline, PipelineLayout, PipelineStep

from .test_round2_cli_fixes import ProducingRunner, _cfg

runner = CliRunner()


@pytest.fixture(autouse=True)
def no_real_probes(monkeypatch):
    """Tool probing is tested separately; keep these tests independent of PATH."""
    from metadarkmatter.core import tool_versions

    monkeypatch.setattr(tool_versions, "probe_tool", lambda exe, args=None: ("", f"{exe} 1.0"))


def _step(plan, step):
    return next(s for s in plan if s.step == step)


class TestMetadata:
    def test_metadata_reaches_ani_classify_report(self, tmp_path):
        md = tmp_path / "genome_metadata.tsv"
        md.write_text("accession\tspecies\tgenus\tfamily\nGCF_1\tA b\tA\tF\n")
        cfg = _cfg(tmp_path, metadata=md)
        plan = Pipeline(cfg, PipelineLayout(root=cfg.output_dir)).plan()
        ani = _step(plan, PipelineStep.ANI)
        assert ani.argv[ani.argv.index("--metadata") + 1] == str(md)
        assert "--representatives-only" in ani.argv
        for st in (PipelineStep.CLASSIFY, PipelineStep.REPORT):
            argv = _step(plan, st).argv
            assert argv[argv.index("--metadata") + 1] == str(md)
            assert md in _step(plan, st).inputs

    def test_warning_without_metadata(self, tmp_path):
        cfg = _cfg(tmp_path)
        p = Pipeline(cfg, PipelineLayout(root=cfg.output_dir))
        assert any("without --metadata" in w for w in p.warnings)
        assert "--metadata" not in _step(p.plan(), PipelineStep.CLASSIFY).argv

    def test_no_warning_with_metadata(self, tmp_path):
        md = tmp_path / "md.tsv"
        md.write_text("accession\tspecies\nGCF_1\tA b\n")
        cfg = _cfg(tmp_path, metadata=md)
        assert not any("--metadata" in w for w in Pipeline(cfg, PipelineLayout(root=cfg.output_dir)).warnings)

    def test_changed_metadata_reruns_ani_but_not_align(self, tmp_path):
        md = tmp_path / "md.tsv"
        md.write_text("accession\tspecies\nGCF_1\tA b\n")
        cfg = _cfg(tmp_path, metadata=md)
        lay = PipelineLayout(root=cfg.output_dir)
        Pipeline(cfg, lay, runner=ProducingRunner()).run()
        md.write_text("accession\tspecies\nGCF_1\tA c changed\n")
        r = ProducingRunner()
        Pipeline(cfg, lay, runner=r).run()
        assert r.ran("ani") and r.ran("classify")
        assert not r.ran("align") and not r.ran("makedb")

    def test_cli_option(self, tmp_path):
        cfg = _cfg(tmp_path)
        md = tmp_path / "md.tsv"
        md.write_text("accession\tspecies\nGCF_1\tA b\n")
        res = runner.invoke(app, [
            "run", "pipeline", "--genomes", str(cfg.genomes_dir),
            "--extracted-reads-1", str(cfg.extracted_reads_1), "--metadata", str(md),
            "--output-dir", str(tmp_path / "o"), "--no-report", "--dry-run"],
            env={"COLUMNS": "2000"})
        assert res.exit_code == 0, res.output
        assert "without --metadata" not in res.output
        # Rich wraps long command lines; join them before searching.
        assert f"--metadata {md}" in res.output.replace("\n", "")


class TestDependencyCaching:
    def test_new_reads_keep_ani_and_makedb_cached(self, tmp_path):
        cfg = _cfg(tmp_path)
        lay = PipelineLayout(root=cfg.output_dir)
        Pipeline(cfg, lay, runner=ProducingRunner()).run()
        reads = Path(cfg.extracted_reads_1)
        reads.write_bytes(b"new reads, different size")
        r = ProducingRunner()
        m = Pipeline(cfg, lay, runner=r).run()
        status = {s.name: s.status.value for s in m.steps}
        assert status["ani"] == "cached" and status["makedb"] == "cached"
        assert status["align"] == "completed" and status["classify"] == "completed"
        assert status["report"] == "completed"

    def test_new_genome_reruns_everything_downstream(self, tmp_path):
        cfg = _cfg(tmp_path)
        lay = PipelineLayout(root=cfg.output_dir)
        Pipeline(cfg, lay, runner=ProducingRunner()).run()
        (Path(cfg.genomes_dir) / "GCF_2.fna").write_text(">c2\nACGT\n")
        m = Pipeline(cfg, lay, runner=ProducingRunner()).run()
        status = {s.name: s.status.value for s in m.steps}
        for step in ("makedb", "align", "ani", "classify", "report"):
            assert status[step] == "completed", (step, status)

    def test_from_still_forces_later_steps(self, tmp_path):
        cfg = _cfg(tmp_path)
        lay = PipelineLayout(root=cfg.output_dir)
        Pipeline(cfg, lay, runner=ProducingRunner()).run()
        m = Pipeline(cfg, lay, runner=ProducingRunner()).run(from_step=PipelineStep.ALIGN)
        status = {s.name: s.status.value for s in m.steps}
        assert status["align"] == status["ani"] == status["classify"] == "completed"
        assert status["makedb"] == "cached"


class TestToolVersions:
    def test_manifest_records_versions_of_used_tools(self, tmp_path):
        cfg = _cfg(tmp_path, ani_backend="fastani")
        lay = PipelineLayout(root=cfg.output_dir)
        Pipeline(cfg, lay, runner=ProducingRunner()).run()
        data = json.loads(lay.manifest.read_text())
        assert data["tool_versions"] == {
            "blastn": "blastn 1.0", "fastANI": "fastANI 1.0", "makeblastdb": "makeblastdb 1.0"}

    def test_mmseqs_and_auto_backend(self, tmp_path):
        from metadarkmatter.models.pipeline_config import Aligner

        cfg = _cfg(tmp_path, aligner=Aligner.MMSEQS2)
        lay = PipelineLayout(root=cfg.output_dir)
        Pipeline(cfg, lay, runner=ProducingRunner()).run()
        tools = json.loads(lay.manifest.read_text())["tool_versions"]
        assert set(tools) == {"mmseqs", "skani", "fastANI"}

    def test_dry_run_does_not_probe(self, tmp_path):
        cfg = _cfg(tmp_path)
        m = Pipeline(cfg, PipelineLayout(root=cfg.output_dir), dry_run=True).run()
        assert m.tool_versions == {}


def test_probe_tool_reads_first_version_line(tmp_path, monkeypatch):
    import importlib

    from metadarkmatter.core import tool_versions

    importlib.reload(tool_versions)  # undo the autouse stub for this test
    exe = tmp_path / "faketool"
    exe.write_text("#!/bin/sh\necho 'faketool 3.2.1'\necho other\n")
    exe.chmod(0o755)
    monkeypatch.setenv("PATH", f"{tmp_path}{os.pathsep}{os.environ['PATH']}")
    path, version = tool_versions.probe_tool("faketool", ["--version"])
    assert version == "faketool 3.2.1" and path.endswith("faketool")
    assert tool_versions.probe_tool("definitely-not-a-tool-xyz")[1] == "not found"
