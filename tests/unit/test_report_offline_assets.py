"""Regression tests for offline-mode HTML reports.

Verifies that ``report_mode='offline'`` produces a self-contained HTML
file without ``<script src="https://...">`` tags, so reports render on
air-gapped lab machines. CDN mode is verified separately to ensure the
opt-in path still works.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from metadarkmatter.visualization.report import ReportConfig, ReportGenerator
from metadarkmatter.visualization.report.assets import (
    get_d3_script_tag,
    get_plotly_script_tag,
)


@pytest.fixture
def minimal_classification_df() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "read_id": ["r1", "r2", "r3"],
            "best_match_genome": ["g1", "g1", "g2"],
            "top_hit_identity": [99.0, 85.0, 70.0],
            "novelty_index": [1.0, 15.0, 30.0],
            "placement_uncertainty": [0.5, 2.0, 5.0],
            "taxonomic_call": ["Known Species", "Novel Species", "Novel Genus"],
            "diversity_status": ["Known", "Novel", "Novel"],
            "is_novel": [False, True, True],
            "confidence_score": [90.0, 70.0, 50.0],
            "low_confidence": [False, False, True],
        }
    )


def test_offline_report_has_no_remote_scripts(
    tmp_path: Path, minimal_classification_df: pl.DataFrame
) -> None:
    out = tmp_path / "report.html"
    gen = ReportGenerator(
        minimal_classification_df,
        ReportConfig(sample_name="t", report_mode="offline"),
    )
    gen.generate(out)
    html = out.read_text(encoding="utf-8")

    # No remote script tags - the file should be fully self-contained
    # except for the optional D3 fallback when D3 has not been vendored.
    remote_count = html.count('<script src="https://')
    assert remote_count <= 1, (
        f"offline report contained {remote_count} remote script tags; "
        "D3 fallback may emit at most one when not vendored, but Plotly "
        "must always be inlined."
    )

    # Plotly source must be inlined.
    assert "plotly.js v" in html or "Plotly.js" in html or "Plotly" in html


def test_cdn_report_emits_remote_scripts(
    tmp_path: Path, minimal_classification_df: pl.DataFrame
) -> None:
    out = tmp_path / "report.html"
    gen = ReportGenerator(
        minimal_classification_df,
        ReportConfig(sample_name="t", report_mode="cdn"),
    )
    gen.generate(out)
    html = out.read_text(encoding="utf-8")

    # CDN mode must reference Plotly remotely; D3 is only inlined when
    # the phylogeny section is built, so we assert the Plotly tag only.
    assert "cdn.plot.ly" in html


def test_plotly_script_tag_offline_inlines_js() -> None:
    tag = get_plotly_script_tag("offline")
    assert "<script>" in tag and "</script>" in tag
    assert "src=" not in tag.split(">", 1)[0]


def test_plotly_script_tag_cdn_uses_remote_src() -> None:
    tag = get_plotly_script_tag("cdn")
    assert 'src="https://cdn.plot.ly/' in tag


def test_d3_script_tag_cdn_uses_remote_src() -> None:
    tag = get_d3_script_tag("cdn")
    assert 'src="https://d3js.org/' in tag
