"""Tests for the subspecies summary and its report table."""

from __future__ import annotations

import polars as pl
import pytest

pytest.importorskip("plotly")

from metadarkmatter.core.subspecies import summarize_subspecies
from tests.unit.report_helpers import GENOMES, make_classifications, make_metadata, render

TITLE = "Subspecies (NCBI species within GTDB species)"
NOTE = (
    "Subspecies are assigned from the best-hit reference genomes of Known Species "
    "reads and are indicative for closely related taxa. Reads whose best hits span "
    "several subspecies are listed as unresolved."
)


def _frame(with_species: bool) -> pl.DataFrame:
    data = {
        "best_match_genome": [
            GENOMES[0], GENOMES[0], GENOMES[0], GENOMES[0], GENOMES[1], GENOMES[1], GENOMES[1],
        ],
        "taxonomic_call": ["Known Species"] * 6 + ["Novel Species"],
        "closest_subspecies": [
            "Francisella tularensis", "Francisella tularensis", "unresolved", None,
            "Francisella novicida", "Francisella novicida", "Francisella novicida",
        ],
    }
    if with_species:
        data["species"] = ["Francisella tularensis"] * 4 + ["Francisella novicida"] * 3
    return pl.DataFrame(data, schema_overrides={"closest_subspecies": pl.Utf8})


def test_summarize_subspecies_counts_and_pct():
    out = summarize_subspecies(_frame(True), None)
    assert out.columns == ["gtdb_species", "subspecies", "reads", "pct"]
    assert out.rows() == [
        ("Francisella novicida", "Francisella novicida", 2, 100.0),
        ("Francisella tularensis", "Francisella tularensis", 2, 50.0),
        ("Francisella tularensis", "unresolved", 1, 25.0),
    ]


def test_summarize_subspecies_uses_metadata_when_species_column_absent():
    out = summarize_subspecies(_frame(False), make_metadata())
    assert out["gtdb_species"].to_list() == [
        "Francisella novicida", "Francisella tularensis", "Francisella tularensis",
    ]
    assert out["pct"].to_list() == [100.0, 50.0, 25.0]


def test_summarize_subspecies_empty_without_column():
    out = summarize_subspecies(make_classifications(), None)
    assert out.is_empty()
    assert out.columns == ["gtdb_species", "subspecies", "reads", "pct"]


def _with_sub(df: pl.DataFrame) -> pl.DataFrame:
    return df.with_columns(
        pl.when(pl.col("taxonomic_call") == "Known Species")
        .then(pl.lit("Francisella <x>"))
        .otherwise(None)
        .alias("closest_subspecies")
    )


def test_report_has_subspecies_table(tmp_path):
    html_text = render(tmp_path, _with_sub(make_classifications()))
    assert TITLE in html_text
    assert NOTE in html_text
    assert "Francisella &lt;x&gt;" in html_text


def test_report_without_subspecies_unchanged(tmp_path):
    assert TITLE not in render(tmp_path)
