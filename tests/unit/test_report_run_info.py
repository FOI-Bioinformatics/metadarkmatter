"""Tests for the Run information panel and the provenance sidecar loader."""

from __future__ import annotations

import json
import re
from pathlib import Path

import polars as pl
import pytest

pytest.importorskip("plotly")

from metadarkmatter.visualization.report.models import ReportConfig
from metadarkmatter.visualization.report.provenance import (
    default_provenance_path,
    human_size,
    load_provenance,
)
from tests.unit.report_helpers import make_ani, make_classifications, render

FIXTURE = Path(__file__).parent.parent / "fixtures" / "provenance_example.json"


@pytest.fixture(scope="module")
def prov() -> dict:
    data = load_provenance(FIXTURE)
    assert data is not None
    return data


def _panel(html_text: str) -> str:
    m = re.search(r'<div class="collapsible-panel[^"]*" id="panel-run-info">(.*?)</section>',
                  html_text, re.S)
    assert m is not None, "no run-info panel"
    return m.group(1)


class TestLoader:
    def test_fixture_matches_documented_schema(self, prov):
        assert prov["schema"] == "metadarkmatter.provenance/1"
        for key in ("metadarkmatter_version", "created_utc", "command", "python", "platform",
                    "seed", "inputs", "config", "outputs", "results"):
            assert key in prov
        for role in ("alignment", "ani", "metadata", "aai", "id_mapping"):
            assert role in prov["inputs"]
        assert {"total_reads", "taxonomic_call_counts"} <= set(prov["results"])
        # config holds ScoringConfig fields
        assert prov["config"]["uncertainty_mode"] in ("second", "max")

    def test_default_path_is_appended_suffix(self):
        assert default_provenance_path("/x/cls.csv") == Path("/x/cls.csv.provenance.json")

    def test_rejects_bad_files(self, tmp_path):
        bad = tmp_path / "bad.json"
        bad.write_text("{not json")
        assert load_provenance(bad) is None
        wrong = tmp_path / "wrong.json"
        wrong.write_text(json.dumps({"schema": "other/1"}))
        assert load_provenance(wrong) is None
        arr = tmp_path / "arr.json"
        arr.write_text("[1,2]")
        assert load_provenance(arr) is None
        assert load_provenance(tmp_path / "missing.json") is None

    def test_human_size(self):
        assert human_size(512) == "512 B"
        assert human_size(1536) == "1.5 KB"
        assert human_size(None) == "-"
        assert human_size(True) == "-"


@pytest.fixture(scope="module")
def panel(tmp_path_factory, prov):
    cfg = ReportConfig(input_names={"classifications": "cls.csv"})
    return _panel(render(tmp_path_factory.mktemp("ri"), provenance=prov, config=cfg))


class TestPanelWithProvenance:
    def test_shows_version_inputs_and_settings(self, panel):
        assert "0.2.0" in panel
        assert "sample.mmseqs2.tsv.gz" in panel and "997231" not in panel
        assert "973.9 KB" in panel  # size formatted
        assert "sha256" in panel and "1f0caaaaaaaa" in panel
        assert "Uncertainty mode" in panel and "second" in panel
        assert "Coverage weighting" in panel and "linear" in panel
        assert "Clade-novelty rule" in panel and ">on<" in panel
        assert "Min. percent identity = 70" in panel
        assert "Provenance was not available" not in panel

    def test_created_utc_only_inside_panel(self, tmp_path, prov):
        text = render(tmp_path, provenance=prov)
        assert text.count("2026-05-02T10:15:30Z") == 1
        assert "2026-05-02T10:15:30Z" in _panel(text)

    def test_matrix_completeness_and_uncertainty_counts(self, panel):
        assert "6 of 6 genome pairs computed (100.0%)" in panel
        assert "Uncertainty: measured / inferred" in panel
        assert "60 (66.7%) / 30 (33.3%)" in panel

    def test_read_count_mismatch_is_flagged(self, tmp_path, prov):
        prov2 = {**prov, "results": {**prov["results"], "total_reads": 123}}
        panel = _panel(render(tmp_path, provenance=prov2))
        assert "123" in panel and "differs from the table" in panel


class TestPanelWithoutProvenance:
    def test_says_unavailable_and_derives_values(self, tmp_path):
        ani = make_ani().with_columns(pl.when(pl.col("genome") == "GCF_000001.1")
                                      .then(None).otherwise(pl.col("GCF_000002.1"))
                                      .alias("GCF_000002.1"))
        cfg = ReportConfig(input_names={"classifications": "cls.csv", "ani": "ani.csv"})
        panel = _panel(render(tmp_path, ani_matrix=ani, config=cfg))
        assert "Provenance was not available" in panel
        assert "cls.csv" in panel and "ani.csv" in panel
        assert "Alignment mode (from --alignment-mode)" in panel
        assert "Reads in classification table" in panel and "<dd>90</dd>" in panel
        assert "genome pairs computed" in panel
        assert "created_utc" not in panel

    def test_html_escaping(self, tmp_path, prov):
        evil = {**prov, "command": ["mdm", "<script>alert(1)</script>"]}
        text = render(tmp_path, provenance=evil)
        assert "<script>alert(1)</script>" not in text

    def test_panel_is_collapsible_button(self, tmp_path):
        text = render(tmp_path)
        assert 'id="panel-header-run-info"' in text
        assert 'aria-controls="panel-body-run-info"' in text


def test_report_cli_autodetects_and_overrides_provenance(tmp_path, prov):
    from typer.testing import CliRunner

    from metadarkmatter.cli.main import app

    cls = tmp_path / "cls.csv"
    make_classifications(30).write_csv(cls)
    default_provenance_path(cls).write_text(json.dumps(prov))
    runner = CliRunner()

    out1 = tmp_path / "auto.html"
    res = runner.invoke(app, ["report", "generate", "--classifications", str(cls),
                              "--output", str(out1), "-q"])
    assert res.exit_code == 0, res.output
    assert "0.2.0" in _panel(out1.read_text(encoding="utf-8"))

    other = tmp_path / "other.json"
    other.write_text(json.dumps({**prov, "metadarkmatter_version": "9.9.9"}))
    out2 = tmp_path / "explicit.html"
    res = runner.invoke(app, ["report", "generate", "--classifications", str(cls),
                              "--output", str(out2), "--provenance", str(other), "-q"])
    assert res.exit_code == 0, res.output
    assert "9.9.9" in _panel(out2.read_text(encoding="utf-8"))

    default_provenance_path(cls).unlink()
    out3 = tmp_path / "none.html"
    res = runner.invoke(app, ["report", "generate", "--classifications", str(cls),
                              "--output", str(out3), "-q"])
    assert res.exit_code == 0, res.output
    assert "Provenance was not available" in _panel(out3.read_text(encoding="utf-8"))
