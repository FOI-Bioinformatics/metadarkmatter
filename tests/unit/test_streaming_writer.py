"""Tests for the incremental StreamingWriter (O(rows) partition append)."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from metadarkmatter.core.io_utils import StreamingWriter, read_dataframe


def _parts() -> list[pl.DataFrame]:
    return [
        pl.DataFrame({"read_id": ["a", "b"], "call": ["Known Species", "Novel Species"]}),
        pl.DataFrame({"read_id": ["c"], "call": ["Novel Genus"]}),
        pl.DataFrame({"read_id": ["d", "e"], "call": ["Known Species", "Ambiguous"]}),
    ]


def test_parquet_incremental_append_concatenates_all_rows(tmp_path: Path):
    out = tmp_path / "out.parquet"
    parts = _parts()
    with StreamingWriter(out, "parquet") as w:
        for p in parts:
            w.write(p)

    combined = pl.read_parquet(out)
    expected = pl.concat(parts)
    assert combined.shape == expected.shape
    assert combined["read_id"].to_list() == expected["read_id"].to_list()
    assert combined["call"].to_list() == expected["call"].to_list()


def test_csv_append_writes_header_once(tmp_path: Path):
    out = tmp_path / "out.csv"
    parts = _parts()
    with StreamingWriter(out, "csv") as w:
        for p in parts:
            w.write(p)

    text = out.read_text()
    # Header appears exactly once (first partition only).
    assert text.count("read_id,call") == 1
    df = read_dataframe(out)
    assert df["read_id"].to_list() == ["a", "b", "c", "d", "e"]


def test_single_partition_parquet(tmp_path: Path):
    out = tmp_path / "single.parquet"
    only = _parts()[0]
    with StreamingWriter(out, "parquet") as w:
        w.write(only)
    assert pl.read_parquet(out)["read_id"].to_list() == ["a", "b"]


def test_close_is_idempotent(tmp_path: Path):
    out = tmp_path / "x.parquet"
    w = StreamingWriter(out, "parquet")
    w.write(_parts()[0])
    w.close()
    w.close()  # second close must not raise
    assert out.exists()
