"""Tests for the `score sensitivity` CLI command."""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.main import app


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


@pytest.fixture
def classification_file(tmp_path: Path) -> Path:
    """Minimal classification CSV with the four required metric columns."""
    df = pl.DataFrame(
        {
            "read_id": [f"r{i}" for i in range(12)],
            "novelty_index": [1.0, 2.5, 3.5, 5.0, 7.0, 10.0,
                              13.0, 16.0, 20.0, 23.0, 26.0, 30.0],
            "placement_uncertainty": [0.2, 0.7, 1.0, 1.5, 2.0, 3.0,
                                      4.0, 5.0, 6.5, 8.0, 10.0, 12.0],
            "num_ambiguous_hits": [1, 2, 3, 2, 1, 4, 2, 1, 3, 1, 5, 2],
            "identity_gap": [None, 1.0, 2.0, 3.0, None, 1.5,
                             4.0, None, 2.5, None, 1.0, 5.0],
            "taxonomic_call": ["Known Species"] * 12,
        }
    )
    out = tmp_path / "classifications.csv"
    df.write_csv(out)
    return out


def test_sensitivity_tsv_output_long_format(
    runner: CliRunner, classification_file: Path, tmp_path: Path
) -> None:
    out_tsv = tmp_path / "sens.tsv"
    result = runner.invoke(
        app,
        [
            "score", "sensitivity",
            "--classifications", str(classification_file),
            "--output", str(out_tsv),
            "--steps", "3",
        ],
    )
    assert result.exit_code == 0, result.output
    df = pl.read_csv(out_tsv, separator="\t")
    expected_cols = {"novelty_known_max", "uncertainty_known_max", "taxonomic_call", "count"}
    assert set(df.columns) == expected_cols
    # 3 threshold points * 7 categories = 21 rows
    assert len(df) == 21
    # Counts sum per threshold-point should equal total reads (12)
    for (n, u), group in df.group_by(["novelty_known_max", "uncertainty_known_max"]):
        assert group["count"].sum() == 12, f"row total wrong at ({n}, {u})"


def test_sensitivity_json_output_keys(
    runner: CliRunner, classification_file: Path, tmp_path: Path
) -> None:
    out_json = tmp_path / "sens.json"
    result = runner.invoke(
        app,
        [
            "score", "sensitivity",
            "--classifications", str(classification_file),
            "--output", str(out_json),
            "--steps", "4",
            "--format", "json",
        ],
    )
    assert result.exit_code == 0, result.output
    payload = json.loads(out_json.read_text())
    assert set(payload.keys()) == {"counts", "novelty_thresholds", "uncertainty_thresholds"}
    assert len(payload["novelty_thresholds"]) == 4
    assert len(payload["uncertainty_thresholds"]) == 4
    for counts in payload["counts"].values():
        assert len(counts) == 4


def test_sensitivity_rejects_missing_columns(
    runner: CliRunner, tmp_path: Path
) -> None:
    bad = tmp_path / "bad.csv"
    pl.DataFrame({"read_id": ["r1"], "novelty_index": [5.0]}).write_csv(bad)
    out = tmp_path / "out.tsv"
    result = runner.invoke(
        app,
        ["score", "sensitivity", "--classifications", str(bad), "--output", str(out)],
    )
    assert result.exit_code == 1
    assert "missing required columns" in result.output


def test_sensitivity_rejects_invalid_ranges(
    runner: CliRunner, classification_file: Path, tmp_path: Path
) -> None:
    out = tmp_path / "out.tsv"
    result = runner.invoke(
        app,
        [
            "score", "sensitivity",
            "--classifications", str(classification_file),
            "--output", str(out),
            "--novelty-min", "5.0",
            "--novelty-max", "3.0",
        ],
    )
    assert result.exit_code == 2
    assert "strictly less than" in result.output
