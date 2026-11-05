"""Tests for the call-comparison table, species breakdown and CSV exports."""

from __future__ import annotations

import polars as pl
import pytest

pytest.importorskip("plotly")

from metadarkmatter.core.novel_diversity.models import NovelCluster
from metadarkmatter.visualization.report.breakdowns import (
    MAX_EXPORT_MATRIX_GENOMES,
    build_call_comparison_html,
    build_species_table_html,
    call_comparison_csv,
    clusters_csv,
    compute_call_comparison,
    compute_species_breakdown,
    matrix_csv,
    species_breakdown_csv,
)
from tests.unit.report_helpers import (
    embedded_json,
    make_ani,
    make_classifications,
    make_metadata,
    render,
)


def _calls_df() -> pl.DataFrame:
    tax = ["Known Species"] * 6 + ["Novel Species"] * 3 + ["Ambiguous"] * 1
    leg = (
        ["Known Species"] * 4 + ["Novel Species"] * 2   # 4 agree, 2 known->novel
        + ["Novel Species"] * 2 + ["Novel Genus"] * 1   # 2 agree, 1 species->genus
        + ["Ambiguous"]                                  # 1 agree
    )
    return pl.DataFrame({"taxonomic_call": tax, "legacy_call": leg})


class TestCallComparison:
    def test_counts_totals_agreement(self):
        cmp = compute_call_comparison(_calls_df())
        assert cmp is not None
        assert cmp.total == 10
        assert cmp.agree == 7
        assert cmp.agreement_pct == pytest.approx(70.0)
        assert cmp.counts[("Known Species", "Novel Species")] == 2
        assert cmp.row_total("Known Species") == 6
        assert cmp.col_total("Novel Species") == 4
        assert cmp.rows == ["Known Species", "Novel Species", "Ambiguous"]
        assert cmp.cols == ["Known Species", "Novel Species", "Novel Genus", "Ambiguous"]

    def test_missing_columns_or_empty(self):
        assert compute_call_comparison(pl.DataFrame({"taxonomic_call": ["a"]})) is None
        assert compute_call_comparison(_calls_df().head(0)) is None

    def test_null_legacy_call_is_labelled(self):
        df = pl.DataFrame({"taxonomic_call": ["Off-target", "Known Species"],
                           "legacy_call": [None, "Known Species"]})
        cmp = compute_call_comparison(df)
        assert cmp is not None and "(none)" in cmp.cols
        assert cmp.agree == 1

    def test_html_explains_both_calls_and_shows_totals(self):
        cmp = compute_call_comparison(_calls_df())
        out = build_call_comparison_html(cmp)
        assert "taxonomic_call" in out and "legacy_call" in out
        assert "do not coincide" in out
        assert "7 of 10 reads (70.0%)" in out
        assert 'class="diag"' in out
        assert '<td class="total">10</td>' in out

    def test_csv_has_totals_row_and_column(self):
        header, rows = call_comparison_csv(compute_call_comparison(_calls_df()))
        assert header[0] == "taxonomic_call" and header[-1] == "total"
        assert rows[-1][0] == "total" and rows[-1][-1] == 10
        assert all(sum(r[1:-1]) == r[-1] for r in rows)


class TestSpeciesBreakdown:
    def test_species_labels_and_counts(self):
        df = make_classifications(90)
        table, cats, uses_species = compute_species_breakdown(df, make_metadata())
        assert uses_species
        # Off-target reads are excluded; the two Allofrancisella genomes merge
        assert table["reads"].sum() == (df["taxonomic_call"] != "Off-target").sum()
        assert "Allofrancisella sp" in table["label"].to_list()
        assert "Off-target" not in cats
        for rec in table.iter_rows(named=True):
            assert sum(rec[c] for c in cats) == rec["reads"]
            assert rec["novel_reads"] == rec["Novel Species"] + rec["Novel Genus"]
        assert table["reads"].to_list() == sorted(table["reads"].to_list(), reverse=True)

    def test_accession_labels_without_metadata(self):
        table, _, uses_species = compute_species_breakdown(make_classifications(40), None)
        assert not uses_species
        assert set(table["label"].to_list()) <= {
            "GCF_000001.1", "GCF_000002.1", "GCF_000003.1", "GCF_000004.1"}

    def test_genome_listing_most_reads_first(self):
        table, _, _ = compute_species_breakdown(make_classifications(90), make_metadata())
        row = table.filter(pl.col("label") == "Allofrancisella sp").row(0, named=True)
        assert row["n_genomes"] == 2
        assert row["genomes"].count("GCF_") == 2

    def test_ties_sorted_by_label(self):
        df = pl.DataFrame({
            "best_match_genome": ["B", "A", "C"],
            "taxonomic_call": ["Known Species"] * 3,
            "top_hit_identity": [99.0, 98.0, 97.0],
        })
        table, _, _ = compute_species_breakdown(df, None)
        assert table["label"].to_list() == ["A", "B", "C"]

    def test_empty_inputs(self):
        table, cats, _ = compute_species_breakdown(make_classifications(9).head(0), None)
        assert len(table) == 0 and cats == []

    def test_html_table_sortable_and_escaped(self):
        df = pl.DataFrame({
            "best_match_genome": ["G1", "G2"],
            "taxonomic_call": ["Novel Species", "Known Species"],
            "top_hit_identity": [90.0, 99.0],
            "species": ["<b>evil</b>", "Fine sp"],
        })
        table, cats, uses = compute_species_breakdown(df, None)
        out = build_species_table_html(table, cats, uses_species=uses)
        assert "<b>evil</b>" not in out and "&lt;b&gt;evil" in out
        assert 'data-sortable="1"' in out and 'aria-sort="descending"' in out
        assert "button" in out

    def test_csv_rows_cover_all_species(self):
        table, cats, uses = compute_species_breakdown(make_classifications(90), make_metadata())
        header, rows = species_breakdown_csv(table, cats, uses_species=uses)
        assert header[:4] == ["species", "best_match_genomes", "n_genomes", "reads"]
        assert len(rows) == len(table)
        assert all(len(r) == len(header) for r in rows)


class TestCsvPayloads:
    def test_matrix_csv_roundtrip(self):
        header, rows = matrix_csv(make_ani())
        assert header[0] == "genome" and len(header) == 5
        assert rows[0] == ["GCF_000001.1", 100.0, 95.5, 85.2, 78.0]

    def test_matrix_csv_without_genome_column(self):
        m = make_ani().drop("genome")
        _, rows = matrix_csv(m)
        assert rows[1][0] == "GCF_000002.1"

    def test_matrix_csv_too_large_is_skipped(self):
        n = MAX_EXPORT_MATRIX_GENOMES + 1
        big = pl.DataFrame({"genome": [f"g{i}" for i in range(n)]}).with_columns(
            *[pl.lit(1.0).alias(f"g{i}") for i in range(3)]
        )
        assert matrix_csv(big) is None

    def test_matrix_csv_nan_becomes_none(self):
        m = pl.DataFrame({"genome": ["a", "b"], "a": [100.0, float("nan")], "b": [None, 100.0]})
        _, rows = matrix_csv(m)
        assert rows[1][1] is None and rows[0][2] is None

    def test_clusters_csv(self):
        c = NovelCluster(
            cluster_id="NSP_001", taxonomic_call="Novel Species", nearest_genome="G1",
            nearest_species="Sp A", nearest_genus="Gen", nearest_family="Fam", novelty_band=5,
            read_count=7, mean_novelty_index=8.0, novelty_min=6.0, novelty_max=9.0,
            mean_placement_uncertainty=0.5, suggested_name="Gen sp. nov.", confidence="High",
            phylogenetic_context="ctx", read_ids=["r1"],
        )
        header, rows = clusters_csv([c])
        assert rows[0][0] == "NSP_001" and rows[0][2] == 7
        assert len(rows[0]) == len(header)


@pytest.fixture(scope="module")
def html_text(tmp_path_factory):
    return render(tmp_path_factory.mktemp("exp"))


class TestEmbeddedExports:
    def test_expected_exports_present(self, html_text):
        for export_id in ("call-comparison", "species-breakdown", "ani-matrix"):
            payload = embedded_json(html_text, f"mdm-export-{export_id}")
            assert payload["filename"].endswith(".csv")
            assert payload["header"] and payload["rows"]
            assert f"mdmDownloadCsv('{export_id}')" in html_text

    def test_ani_export_matches_matrix(self, html_text):
        payload = embedded_json(html_text, "mdm-export-ani-matrix")
        assert payload["header"][0] == "genome"
        assert payload["rows"][0][0] == "GCF_000001.1"

    def test_call_comparison_in_confidence_panel(self, html_text):
        assert 'id="call-comparison"' in html_text
        assert "Bayesian call vs legacy threshold call" in html_text

    def test_species_stack_replaces_plain_bar(self, html_text):
        assert "Stacked by Classification" in html_text
        assert 'id="species-breakdown"' in html_text

    def test_no_legacy_column_means_no_comparison(self, tmp_path):
        out = render(tmp_path, make_classifications(30).drop("legacy_call"))
        assert 'id="call-comparison"' not in out
