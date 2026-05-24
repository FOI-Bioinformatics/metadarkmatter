"""End-to-end determinism check.

Runs the score classify CLI twice on the same input with the same
METADARKMATTER_SEED and asserts the resulting CSVs are byte-identical.
This guards against future regressions where a non-deterministic
ordering or unseeded RNG sneaks into the pipeline.
"""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.main import app


@pytest.fixture
def blast_input(tmp_path: Path) -> Path:
    """A small BLAST tabular file covering a handful of reads/genomes."""
    rows = []
    # 5 reads, each with 2-3 hits to a 3-genome pangenome.
    for i in range(5):
        rows.extend(
            [
                {
                    "qseqid": f"read_{i}",
                    "sseqid": "GCF_000123456.1|c1",
                    "pident": 99.0 - i * 0.3,
                    "length": 150,
                    "mismatch": 1,
                    "gapopen": 0,
                    "qstart": 1,
                    "qend": 150,
                    "sstart": 1,
                    "send": 150,
                    "evalue": 1e-50,
                    "bitscore": 280.0 - i,
                },
                {
                    "qseqid": f"read_{i}",
                    "sseqid": "GCF_000789012.1|c1",
                    "pident": 95.0 - i * 0.3,
                    "length": 150,
                    "mismatch": 5,
                    "gapopen": 0,
                    "qstart": 1,
                    "qend": 150,
                    "sstart": 1,
                    "send": 150,
                    "evalue": 1e-40,
                    "bitscore": 270.0 - i,
                },
            ]
        )
    blast_path = tmp_path / "input.blast.tsv"
    pl.DataFrame(rows).write_csv(blast_path, separator="\t", include_header=False)
    return blast_path


@pytest.fixture
def ani_input(tmp_path: Path) -> Path:
    ani_path = tmp_path / "ani.csv"
    df = pl.DataFrame(
        {
            "genome": ["GCF_000123456.1", "GCF_000789012.1", "GCA_000111222.1"],
            "GCF_000123456.1": [100.0, 95.5, 80.0],
            "GCF_000789012.1": [95.5, 100.0, 82.0],
            "GCA_000111222.1": [80.0, 82.0, 100.0],
        }
    )
    df.write_csv(ani_path)
    return ani_path


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _run_classify(blast: Path, ani: Path, out: Path) -> None:
    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "score", "classify",
            "--alignment", str(blast),
            "--ani", str(ani),
            "--output", str(out),
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output


def test_same_seed_produces_identical_csv(
    blast_input: Path, ani_input: Path, tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    out_a = tmp_path / "run_a.csv"
    out_b = tmp_path / "run_b.csv"

    monkeypatch.setenv("METADARKMATTER_SEED", "42")
    _run_classify(blast_input, ani_input, out_a)
    _run_classify(blast_input, ani_input, out_b)

    assert _sha256(out_a) == _sha256(out_b), (
        "Two runs of score classify with the same input and seed "
        "produced different CSVs - a non-determinism regression has "
        "been introduced. Diff the two CSVs to find the offending field."
    )


def test_classification_columns_present(
    blast_input: Path, ani_input: Path, tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Sanity: the determinism test would be meaningless if the run produced
    no rows or missing columns. Keep this here so a silent classifier
    short-circuit cannot dress up as a passing determinism test.
    """
    out = tmp_path / "run.csv"
    monkeypatch.setenv("METADARKMATTER_SEED", "42")
    _run_classify(blast_input, ani_input, out)
    df = pl.read_csv(out)
    assert len(df) == 5  # one row per read
    for col in ("read_id", "taxonomic_call", "novelty_index", "placement_uncertainty"):
        assert col in df.columns
