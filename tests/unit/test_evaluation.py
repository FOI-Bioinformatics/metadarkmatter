"""Tests for the evaluation module and the `score evaluate` CLI."""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.main import app
from metadarkmatter.core.classification.evaluation import (
    evaluate_classifications,
)


def _df(predictions: list[str], truths: list[str], confidences: list[float] | None = None) -> pl.DataFrame:
    data: dict[str, list] = {
        "taxonomic_call": predictions,
        "true_category": truths,
    }
    if confidences is not None:
        data["confidence_score"] = confidences
    return pl.DataFrame(data)


# ---------------------------------------------------------------------------
# Library: evaluate_classifications
# ---------------------------------------------------------------------------


def test_perfect_classifier_has_full_accuracy() -> None:
    df = _df(
        predictions=["A", "B", "A", "B"],
        truths=["A", "B", "A", "B"],
        confidences=[100.0, 100.0, 100.0, 100.0],
    )
    res = evaluate_classifications(df)
    assert res.overall_accuracy == 1.0
    assert res.expected_calibration_error == pytest.approx(0.0)


def test_per_category_precision_recall() -> None:
    # Truths: 3 A, 2 B. Predicted: 4 A, 1 B. One B mislabelled as A.
    df = _df(
        predictions=["A", "A", "A", "A", "B"],
        truths=["A", "A", "A", "B", "B"],
    )
    res = evaluate_classifications(df, confidence_column=None)
    by_cat = {c.category: c for c in res.per_category}
    # A: TP=3, predicted=4, support=3 -> precision=0.75, recall=1.0
    assert by_cat["A"].precision == pytest.approx(0.75)
    assert by_cat["A"].recall == pytest.approx(1.0)
    # B: TP=1, predicted=1, support=2 -> precision=1.0, recall=0.5
    assert by_cat["B"].precision == pytest.approx(1.0)
    assert by_cat["B"].recall == pytest.approx(0.5)


def test_confusion_matrix_is_keyed_by_truth_then_prediction() -> None:
    df = _df(
        predictions=["A", "B", "A"],
        truths=["A", "A", "B"],
    )
    res = evaluate_classifications(df, confidence_column=None)
    assert res.confusion["A"]["A"] == 1
    assert res.confusion["A"]["B"] == 1
    assert res.confusion["B"]["A"] == 1


def test_overconfident_classifier_has_high_ece() -> None:
    # All predictions wrong, all at 100% confidence -> ECE near 100%.
    df = _df(
        predictions=["A"] * 10,
        truths=["B"] * 10,
        confidences=[100.0] * 10,
    )
    res = evaluate_classifications(df, n_confidence_bins=5)
    # accuracy=0, mean confidence=1.0 -> gap=1.0 in the top bin.
    assert res.expected_calibration_error == pytest.approx(1.0)


def test_ece_skipped_when_confidence_column_missing() -> None:
    df = _df(
        predictions=["A", "B"],
        truths=["A", "B"],
    )
    res = evaluate_classifications(df, confidence_column=None)
    assert res.expected_calibration_error is None
    assert res.confidence_bin_centers == []


def test_missing_required_columns_raises() -> None:
    df = pl.DataFrame({"only_column": [1, 2]})
    with pytest.raises(ValueError, match="missing"):
        evaluate_classifications(df)


def test_empty_input_raises() -> None:
    df = pl.DataFrame({"taxonomic_call": [], "true_category": []})
    with pytest.raises(ValueError, match="empty"):
        evaluate_classifications(df)


# ---------------------------------------------------------------------------
# CLI: score evaluate
# ---------------------------------------------------------------------------


def test_cli_evaluate_prints_summary(tmp_path: Path) -> None:
    runner = CliRunner()
    df = _df(
        predictions=["A", "A", "B"],
        truths=["A", "B", "B"],
        confidences=[80.0, 60.0, 90.0],
    )
    pred = tmp_path / "p.csv"
    df.write_csv(pred)
    result = runner.invoke(app, ["score", "evaluate", "--predictions", str(pred)])
    assert result.exit_code == 0, result.output
    assert "Overall accuracy" in result.output
    assert "Expected calibration error" in result.output


def test_cli_evaluate_writes_json(tmp_path: Path) -> None:
    runner = CliRunner()
    df = _df(predictions=["A"] * 4, truths=["A"] * 3 + ["B"], confidences=[90.0] * 4)
    pred = tmp_path / "p.csv"
    df.write_csv(pred)
    out_json = tmp_path / "eval.json"

    result = runner.invoke(
        app,
        [
            "score", "evaluate",
            "--predictions", str(pred),
            "--output", str(out_json),
            "--bins", "5",
        ],
    )
    assert result.exit_code == 0, result.output
    payload = json.loads(out_json.read_text())
    assert set(payload.keys()) >= {
        "n_rows", "overall_accuracy", "expected_calibration_error",
        "per_category", "confusion", "confidence_bins",
    }
    assert payload["n_rows"] == 4


def test_cli_evaluate_exits_when_truth_column_missing(tmp_path: Path) -> None:
    runner = CliRunner()
    pl.DataFrame({"taxonomic_call": ["A"]}).write_csv(tmp_path / "p.csv")
    result = runner.invoke(
        app, ["score", "evaluate", "--predictions", str(tmp_path / "p.csv")]
    )
    assert result.exit_code == 1
    assert "Truth column" in result.output
