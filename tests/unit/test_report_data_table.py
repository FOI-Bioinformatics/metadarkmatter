"""Tests for the client-side Data tab table payload and sampling."""

from __future__ import annotations

import re

import polars as pl
import pytest

pytest.importorskip("plotly")

from metadarkmatter.visualization.report.data_table import (
    active_columns,
    build_table_payload,
    select_table_rows,
)
from metadarkmatter.visualization.report.models import ReportConfig
from tests.unit.report_helpers import (
    embedded_json,
    make_classifications,
    make_metadata,
    render,
)


def _decode(payload: dict, i: int) -> dict:
    """Resolve row ``i`` of a payload back to plain values."""
    cols = payload["cols"]
    out = {}
    for c in payload["columns"]:
        k, kind = c["key"], c["kind"]
        if kind == "genome":
            out[k] = payload["genomes"][cols["genome"][i]]
        elif kind == "species":
            out[k] = payload["species"][cols["genome"][i]]
        elif kind == "genus":
            out[k] = payload["genus"][cols["genome"][i]]
        elif kind == "dict":
            out[k] = payload["dicts"][k][cols[k][i]]
        else:
            out[k] = cols[k][i]
    return out


class TestPayload:
    def test_roundtrip_matches_dataframe(self):
        df = make_classifications(40)
        payload = build_table_payload(df, make_metadata(), page_size=25)
        assert payload["n"] == 40
        assert payload["pageSize"] == 25
        for i, row in enumerate(df.iter_rows(named=True)):
            dec = _decode(payload, i)
            assert dec["read_id"] == row["read_id"]
            assert dec["best_match"] == row["best_match_genome"]
            assert dec["classification"] == row["taxonomic_call"]
            assert dec["legacy_call"] == row["legacy_call"]
            assert dec["identity"] == pytest.approx(row["top_hit_identity"], abs=0.005)
            if row["placement_uncertainty"] is None:
                assert dec["uncertainty"] is None
            else:
                assert dec["uncertainty"] == pytest.approx(row["placement_uncertainty"], abs=0.005)
            assert dec["is_novel"] == int(row["is_novel"])
        # species/genus come from the metadata, one entry per distinct genome
        assert len(payload["species"]) == len(payload["genomes"])
        assert "Francisella tularensis" in payload["species"]

    def test_floats_are_rounded_and_strings_dictionary_encoded(self):
        df = make_classifications(60).with_columns(
            (pl.col("top_hit_identity") + 0.123456789).alias("top_hit_identity")
        )
        payload = build_table_payload(df, None, page_size=100)
        for v in payload["cols"]["identity"]:
            assert round(v, 2) == v
        assert len(payload["dicts"]["classification"]) <= 9
        assert all(isinstance(i, int) for i in payload["cols"]["classification"])

    def test_optional_columns_only_when_data_present(self):
        df = make_classifications(20).drop("legacy_call", "uncertainty_type")
        keys = [c.key for c in active_columns(df)]
        assert "legacy_call" not in keys and "uncertainty_type" not in keys
        assert "confidence_score" in keys
        assert "discovery_score" not in keys

    def test_species_from_dataframe_columns_without_metadata(self):
        df = make_classifications(8).with_columns(
            pl.lit("Francisella sp. X").alias("species"), pl.lit("Francisella").alias("genus")
        )
        payload = build_table_payload(df, None, page_size=10)
        assert set(payload["species"]) == {"Francisella sp. X"}

    def test_empty_frame(self):
        df = make_classifications(5).head(0)
        payload = build_table_payload(df, None, page_size=10)
        assert payload["n"] == 0
        assert payload["cols"]["read_id"] == []

    def test_payload_is_compact(self):
        from metadarkmatter.visualization.report.embed import json_for_script

        df = make_classifications(90)
        text = json_for_script(build_table_payload(df, make_metadata(), page_size=100), compact=True)
        # The server-side <tr> approach cost ~840 bytes per row.
        assert len(text) / 90 < 200


@pytest.fixture(scope="module")
def html_text(tmp_path_factory):
    return render(tmp_path_factory.mktemp("tbl"))


class TestEmbedding:
    def test_no_server_side_rows(self, html_text):
        assert re.search(r'<tbody id="tableBody">\s*</tbody>', html_text)
        payload = embedded_json(html_text, "mdm-table-data")
        assert payload["n"] == 90

    def test_script_close_in_names_is_escaped(self, tmp_path):
        df = make_classifications(6).with_columns(pl.lit("</script><b>x").alias("read_id"))
        out = render(tmp_path, df)
        payload = embedded_json(out, "mdm-table-data")
        assert payload["cols"]["read_id"][0] == "</script><b>x"
        assert "</script><b>x" not in out

    def test_controls_present(self, html_text):
        for token in ("tableSearch", "exportTableTSV", "clearFiltersTop", "columnMenu",
                      "prevBtn", "nextBtn", 'data-col="read_id"'):
            assert token in html_text

    def test_no_truncation_notice_by_default(self, html_text):
        assert "Table truncated" not in html_text


class TestSampling:
    def _frame(self) -> pl.DataFrame:
        calls = ["Known Species"] * 500 + ["Novel Species"] * 30 + ["Novel Genus"] * 5 + ["Ambiguous"] * 200
        n = len(calls)
        return pl.DataFrame({
            "read_id": [f"r{i:05d}" for i in range(n)],
            "taxonomic_call": calls,
            "novelty_index": [float((i * 7) % 31) for i in range(n)],
        })

    def test_no_limit_or_fits_returns_everything(self):
        df = self._frame()
        assert select_table_rows(df, None).equals(df)
        assert select_table_rows(df, len(df)).equals(df)

    def test_exact_size_and_original_order(self):
        df = self._frame()
        out = select_table_rows(df, 100)
        assert len(out) == 100
        ids = out["read_id"].to_list()
        assert ids == sorted(ids)

    def test_rare_categories_kept_in_full(self):
        out = select_table_rows(self._frame(), 100)
        counts = dict(out.group_by("taxonomic_call").len().iter_rows())
        assert counts["Novel Genus"] == 5
        assert counts["Novel Species"] == 30
        # remaining budget is shared by the two large categories
        assert counts["Known Species"] == counts["Ambiguous"] == 32 or (
            counts["Known Species"] + counts["Ambiguous"] == 65
        )

    def test_highest_novelty_first_within_category(self):
        df = self._frame()
        out = select_table_rows(df, 100)
        kept = out.filter(pl.col("taxonomic_call") == "Known Species")["novelty_index"].min()
        dropped = df.filter(
            (pl.col("taxonomic_call") == "Known Species") & ~pl.col("read_id").is_in(out["read_id"].implode())
        )["novelty_index"].max()
        assert kept >= dropped

    def test_deterministic_and_row_order_independent(self):
        df = self._frame()
        a = select_table_rows(df, 80)["read_id"].to_list()
        b = select_table_rows(df, 80)["read_id"].to_list()
        shuffled = df.sample(fraction=1.0, shuffle=True, seed=3)
        c = sorted(select_table_rows(shuffled, 80)["read_id"].to_list())
        assert a == b
        # ties on novelty are broken by read_id, so the chosen set is stable
        assert c == sorted(a)

    def test_limit_smaller_than_category_count(self):
        out = select_table_rows(self._frame(), 3)
        assert len(out) == 3

    def test_report_states_truncation(self, tmp_path):
        out = render(
            tmp_path,
            config=ReportConfig(max_table_rows=30, classifications_name="cls.csv"),
        )
        assert "Table truncated" in out
        assert "30 of 90" in out
        assert "cls.csv" in out
        assert "--max-table-rows" in out
        assert embedded_json(out, "mdm-table-data")["n"] == 30
        # summary counts still describe all reads
        assert re.search(r'class="stat-value">90<', out)


class TestCli:
    def test_max_table_rows_option(self, cli_runner, temp_classification_file, temp_dir):
        from metadarkmatter.cli.main import app

        out_path = temp_dir / "r.html"
        result = cli_runner.invoke(
            app,
            ["report", "generate", "--classifications", str(temp_classification_file),
             "--output", str(out_path), "--max-table-rows", "2", "-q"],
        )
        assert result.exit_code == 0, result.output
        payload = embedded_json(out_path.read_text(encoding="utf-8"), "mdm-table-data")
        assert payload["n"] <= 2
