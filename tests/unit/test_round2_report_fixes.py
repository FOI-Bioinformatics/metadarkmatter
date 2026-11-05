"""Regression tests for round-2 HTML report fixes."""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

pytest.importorskip("plotly")

from metadarkmatter.visualization.plots.classification_charts import (
    DiversitySunburstChart,
)
from metadarkmatter.visualization.report import (
    ReportConfig,
    ReportGenerator,
)
from metadarkmatter.visualization.report.components.heatmap_builder import (
    ANI_HEATMAP_CONFIG,
    _build_similarity_heatmap,
    compute_missing_mask,
    compute_similarity_stats,
    count_missing_pairs,
)
from metadarkmatter.visualization.report.templates import METHODS_SECTION_TEMPLATE

G = ["GCF_1", "GCF_2", "GCF_3", "GCF_4"]


def _ani(with_zero: bool = True) -> pl.DataFrame:
    z = 0.0 if with_zero else 80.0
    return pl.DataFrame({
        "genome": G,
        "GCF_1": [100.0, 96.0, 85.0, z],
        "GCF_2": [96.0, 100.0, 84.0, z],
        "GCF_3": [85.0, 84.0, 100.0, 90.0],
        "GCF_4": [z, z, 90.0, 100.0],
    })


def _cls(n: int = 40, extra_cols: dict | None = None) -> pl.DataFrame:
    calls = ["Known Species"] * (n // 2) + ["Novel Species"] * (n - n // 2)
    data = {
        "read_id": [f"r{i}" for i in range(n)],
        "best_match_genome": [G[i % 3] for i in range(n)],
        "top_hit_identity": [99.0] * (n // 2) + [90.0] * (n - n // 2),
        "novelty_index": [1.0] * (n // 2) + [10.0] * (n - n // 2),
        "placement_uncertainty": [0.5] * n,
        "num_ambiguous_hits": [1, 2] * (n // 2),
        "taxonomic_call": calls,
        "legacy_call": calls,
        "confidence_score": [80.0] * n,
        "is_novel": [c == "Novel Species" for c in calls],
    }
    data.update(extra_cols or {})
    return pl.DataFrame(data)


def _html(df=None, ani=None, **cfg) -> str:
    gen = ReportGenerator(df if df is not None else _cls(), config=ReportConfig(**cfg), ani_matrix=ani)
    return gen._build_html()


# 1. sunburst
def test_sunburst_without_off_target_has_no_orphans():
    summary = {"total_reads": 10, "known_species": 6, "novel_species": 4, "off_target": 0}
    fig = DiversitySunburstChart(summary).create_figure()
    tr = fig.data[0]
    ids = set(tr.ids)
    assert all(p == "" or p in ids for p in tr.parents)
    assert "Off-target" not in ids
    assert not any(i.startswith("Off-target-") for i in ids)
    assert all(v > 0 for i, v in zip(tr.ids, tr.values, strict=True) if "-" in i)


# 2. methods
def test_methods_panel_describes_bayesian_call():
    m = METHODS_SECTION_TEMPLATE
    assert "taxonomic_call" in m and "legacy_call" in m
    assert "log<sub>2</sub>(6)" in m
    assert "6.2%" in m and "18.6%" in m and "27.6%" in m
    assert 'uncertainty_mode="second"' in m
    assert "margin_score" not in m


# 3/4. ANI stats and missing pairs
def test_similarity_stats_count_unordered_pairs_and_missing():
    z = np.array([[100, 96, 0], [96, 100, 85], [0, 85, 100]], dtype=float)
    miss = compute_missing_mask(z)
    stats = compute_similarity_stats(z, 70.0, 95.0, 93.0, missing_mask=miss)
    assert stats.total_pairs == 3
    assert stats.missing_count == 1
    assert stats.same_category_count == 1
    assert stats.different_category_count == 1
    assert stats.same_category_count + stats.boundary_count + stats.different_category_count == 2


def test_count_missing_pairs_dataframe():
    assert count_missing_pairs(_ani()) == (2, 6)
    assert count_missing_pairs(_ani(with_zero=False)) == (0, 6)


def test_heatmap_missing_cells_are_null_not_70():
    fig, stats, _ = _build_similarity_heatmap(_ani(), {g: g for g in G}, ANI_HEATMAP_CONFIG)
    z = fig.data[0].z
    flat = [v for row in z for v in row]
    assert flat.count(None) == 4  # 2 pairs x 2 triangles
    assert 70.0 not in flat
    assert stats.missing_count == 2 and stats.total_pairs == 6
    assert "not computed" in fig.data[1].hovertemplate


def test_report_notes_missing_pairs_on_reference_and_summary():
    out = _html(ani=_ani())
    assert out.count("have no computed ANI") == 2  # summary QC + reference tab
    assert "treated as 70% ANI" in out
    assert "2 of 6 genome pairs (33.3%)" in out
    assert "have no computed ANI" not in _html(ani=_ani(with_zero=False))


# 5. discovery card
def test_discovery_card_hidden_without_scores():
    out = _html()
    assert "Mean Discovery Score" not in out
    assert "previously reported" not in out
    df = _cls(extra_cols={"discovery_score": [50.0] * 20 + [None] * 20, "alignment_quality": [1.0] * 40})
    assert "Mean Discovery Score" in _html(df)


# 6. multi-sample
def test_multi_sample_includes_all_categories_and_offline_plotly():
    from metadarkmatter.visualization.report.multi_generator import (
        MultiSampleConfig,
        MultiSampleReportGenerator,
    )

    calls = ["Known Species"] * 5 + ["Ambiguous"] * 3 + ["Species Boundary"] * 2 + ["Off-target"] * 2 + ["Unclassified"]
    n = len(calls)
    df = pl.DataFrame({
        "taxonomic_call": calls,
        "novelty_index": [1.0] * n,
        "placement_uncertainty": [0.5] * n,
    })
    gen = MultiSampleReportGenerator({"a": df, "b": df}, MultiSampleConfig())
    s = gen.summaries["a"]
    assert s["ambiguous"] == 3 and s["off_target"] == 2 and s["other"] == 0
    from metadarkmatter.visualization.plots.multi_sample import MultiSampleBarChart

    fig = MultiSampleBarChart(gen.summaries, stacked=True).create_figure()
    assert {t.name for t in fig.data} >= {"Ambiguous", "Species Boundary", "Unclassified", "Off-target"}
    totals = sum(t.y[0] for t in fig.data)
    assert totals == pytest.approx(100.0)
    html_out = gen._build_html()
    assert '<script src="https://cdn.plot.ly' not in html_out
    assert "<th>Ambiguous</th>" in html_out and "<th>Off-target</th>" in html_out


# 7. D3 warning
def test_offline_d3_warning(monkeypatch, tmp_path):
    from typer.testing import CliRunner

    from metadarkmatter.cli.main import app
    from metadarkmatter.visualization.report import assets

    monkeypatch.setattr(assets, "is_d3_vendored", lambda: False)
    cls = tmp_path / "c.csv"
    _cls().write_csv(cls)
    ani = tmp_path / "a.csv"
    _ani(with_zero=False).write_csv(ani)
    res = CliRunner().invoke(
        app,
        ["report", "generate", "--classifications", str(cls), "--ani", str(ani),
         "--output", str(tmp_path / "r.html")],
    )
    assert res.exit_code == 0, res.output
    assert "vendor-d3" in res.output
    assert "vendor-d3" in (tmp_path / "r.html").read_text()


# 8. responsive plots
def test_plots_have_no_fixed_width():
    import json

    gen = ReportGenerator(_cls(), config=ReportConfig(), ani_matrix=_ani())
    gen._build_html()
    assert gen._plot_data
    for plot_id, js in gen._plot_data.items():
        layout = json.loads(js)["layout"]
        assert layout.get("width") is None, plot_id
        assert layout.get("autosize") is True, plot_id


# 9. data table
def test_data_table_optional_columns_and_clear():
    out = _html()
    assert 'onclick="clearAllFilters()"' in out and "clearFiltersTop" in out
    assert 'data-col="16"' in out and "Legacy Call" in out
    assert "var cols = [11, 12, 16, 17]" in out or "var cols = [12, 16, 17]" in out or "16, 17]" in out
    assert "13," not in out.split("var cols = ")[1].split("]")[0]


# 10. labels
def test_labels_titles_and_denominators():
    df = _cls(44)
    df = df.with_columns(
        pl.when(pl.int_range(pl.len()) >= 40).then(pl.lit("Off-target")).otherwise(pl.col("taxonomic_call")).alias("taxonomic_call"),
        pl.when(pl.int_range(pl.len()) >= 40).then(None).otherwise(pl.col("placement_uncertainty")).alias("placement_uncertainty"),
    )
    df = df.with_columns(pl.lit("Sp A").alias("species"))
    out = _html(df)
    assert 'title="Known: ' in out and "reads" in out
    assert "of on-target reads" in out
    assert "Top 20 Species" not in out
