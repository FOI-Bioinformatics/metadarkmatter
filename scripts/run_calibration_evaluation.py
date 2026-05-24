#!/usr/bin/env python
"""
Compare baseline vs calibrated classifier on a labelled corpus.

Phase 4 helper. Given a single labelled metrics TSV (with both the
baseline and calibrated ``taxonomic_call`` / ``confidence_score``
predictions), produces a Markdown summary table covering top-1
accuracy, expected calibration error, and per-category precision and
recall, so the lab can decide whether the calibrated config meets the
shipment criteria documented in docs/CALIBRATION.md.

Inputs must be metrics TSVs as produced by ``build_metrics_tsv.py``
with ``--include-confidence`` (the default), one for the baseline run
and one for the calibrated run.

Usage:

    python scripts/run_calibration_evaluation.py \\
        --family Francisellaceae \\
        --baseline corpora/francisella/baseline_metrics.tsv \\
        --calibrated corpora/francisella/calibrated_metrics.tsv \\
        --output corpora/francisella/calibration_report.md
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
from pathlib import Path

import polars as pl

from metadarkmatter.core.classification.evaluation import (
    EvaluationResult,
    evaluate_classifications,
)


def _evaluate(path: Path) -> EvaluationResult:
    df = pl.read_csv(path, separator="\t")
    return evaluate_classifications(
        df,
        truth_column="true_category",
        prediction_column="taxonomic_call",
        confidence_column="confidence_score" if "confidence_score" in df.columns else None,
    )


def _format_pct(value: float | None) -> str:
    if value is None:
        return "-"
    return f"{value * 100:.2f}%"


def _delta(a: float | None, b: float | None) -> str:
    if a is None or b is None:
        return "-"
    return f"{(b - a) * 100:+.2f} pp"


def _accept(baseline: EvaluationResult, calibrated: EvaluationResult) -> tuple[bool, list[str]]:
    """Apply the shipment criteria from docs/CALIBRATION.md."""
    notes: list[str] = []
    ok = True

    # Top-1 accuracy must not regress more than 1pp.
    acc_delta = (calibrated.overall_accuracy - baseline.overall_accuracy) * 100
    if acc_delta < -1.0:
        ok = False
        notes.append(
            f"Top-1 accuracy regressed by {-acc_delta:.2f} pp (> 1 pp threshold)."
        )

    # ECE must drop by >= 5pp to justify shipping.
    if baseline.expected_calibration_error is not None and calibrated.expected_calibration_error is not None:
        ece_drop = (
            baseline.expected_calibration_error
            - calibrated.expected_calibration_error
        ) * 100
        if ece_drop < 5.0:
            ok = False
            notes.append(
                f"ECE drop is only {ece_drop:.2f} pp; need >= 5 pp to justify shipment."
            )

    # No per-category recall may regress by more than 5 pp.
    baseline_by_cat = {c.category: c for c in baseline.per_category}
    for cat in calibrated.per_category:
        b = baseline_by_cat.get(cat.category)
        if b is None:
            continue
        recall_delta = (cat.recall - b.recall) * 100
        if recall_delta < -5.0:
            ok = False
            notes.append(
                f"Category {cat.category} recall regressed by "
                f"{-recall_delta:.2f} pp (> 5 pp threshold)."
            )

    return ok, notes


def _render_markdown(
    family: str,
    baseline: EvaluationResult,
    calibrated: EvaluationResult,
) -> str:
    lines = [
        f"# Calibration evaluation: {family}",
        "",
        f"- Rows scored: {baseline.n_rows:,}",
        "",
        "## Overall",
        "",
        "| Metric | Baseline | Calibrated | Delta |",
        "|---|---|---|---|",
        "| Top-1 accuracy | "
        + _format_pct(baseline.overall_accuracy)
        + " | "
        + _format_pct(calibrated.overall_accuracy)
        + " | "
        + _delta(baseline.overall_accuracy, calibrated.overall_accuracy)
        + " |",
        "| Expected calibration error | "
        + _format_pct(baseline.expected_calibration_error)
        + " | "
        + _format_pct(calibrated.expected_calibration_error)
        + " | "
        + _delta(
            baseline.expected_calibration_error,
            calibrated.expected_calibration_error,
        )
        + " |",
        "",
        "## Per-category",
        "",
        "| Category | Support | Baseline prec | Calibrated prec | Δ prec | Baseline rec | Calibrated rec | Δ rec |",
        "|---|---|---|---|---|---|---|---|",
    ]
    by_cat_b = {c.category: c for c in baseline.per_category}
    by_cat_c = {c.category: c for c in calibrated.per_category}
    for cat in sorted(set(by_cat_b) | set(by_cat_c)):
        b = by_cat_b.get(cat)
        c = by_cat_c.get(cat)
        support = b.support if b else (c.support if c else 0)
        bp = b.precision if b else None
        cp = c.precision if c else None
        br = b.recall if b else None
        cr = c.recall if c else None
        lines.append(
            f"| {cat} | {support} | {_format_pct(bp)} | {_format_pct(cp)} | "
            f"{_delta(bp, cp)} | {_format_pct(br)} | {_format_pct(cr)} | "
            f"{_delta(br, cr)} |"
        )

    ok, notes = _accept(baseline, calibrated)
    lines.extend(["", "## Shipment decision", ""])
    if ok:
        lines.append("**PASS** - calibrated config meets all shipment criteria.")
    else:
        lines.append("**FAIL** - calibrated config does not meet shipment criteria:")
        lines.extend([f"- {note}" for note in notes])
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family", required=True)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--calibrated", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    baseline = _evaluate(args.baseline)
    calibrated = _evaluate(args.calibrated)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(_render_markdown(args.family, baseline, calibrated))

    # Echo a compact dict to stdout for downstream automation.
    import json

    print(
        json.dumps(
            {
                "family": args.family,
                "baseline": {
                    "accuracy": baseline.overall_accuracy,
                    "ece": baseline.expected_calibration_error,
                },
                "calibrated": {
                    "accuracy": calibrated.overall_accuracy,
                    "ece": calibrated.expected_calibration_error,
                },
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
