"""
Classification evaluation against a labelled corpus.

Given a DataFrame of predictions (with at minimum ``taxonomic_call``,
optionally ``confidence_score``) and a truth column or separate
labels, this module computes the diagnostics that make a calibration
loop closeable: top-1 accuracy, per-category precision/recall, a
confusion matrix, and an expected calibration error (ECE).

Used by ``metadarkmatter score evaluate`` and exposed as a small
programmatic API.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import polars as pl


@dataclass
class CategoryStats:
    category: str
    support: int  # number of true rows in this category
    predicted: int  # number of rows predicted as this category
    correct: int  # number of true positives
    precision: float  # TP / predicted (0 when predicted == 0)
    recall: float  # TP / support (0 when support == 0)


@dataclass
class EvaluationResult:
    """Summary of a classifier's behaviour on a labelled corpus."""

    n_rows: int
    overall_accuracy: float
    per_category: list[CategoryStats]
    confusion: dict[str, dict[str, int]]  # confusion[true_label][predicted_label] = count
    expected_calibration_error: float | None  # None when no confidence column
    confidence_bin_centers: list[float]
    confidence_bin_accuracy: list[float]
    confidence_bin_counts: list[int]


def evaluate_classifications(
    df: pl.DataFrame,
    *,
    truth_column: str = "true_category",
    prediction_column: str = "taxonomic_call",
    confidence_column: str | None = "confidence_score",
    n_confidence_bins: int = 10,
) -> EvaluationResult:
    """Score a classifier against a labelled corpus.

    Args:
        df: Polars DataFrame with at least ``prediction_column`` and
            ``truth_column``. ``confidence_column`` is optional; when
            absent or ``None``, the expected calibration error and
            confidence-vs-accuracy histogram are omitted.
        truth_column: Column holding the ground-truth category name.
        prediction_column: Column holding the predicted category name.
        confidence_column: Column holding a 0-100 confidence score for
            the prediction, or ``None`` to skip calibration metrics.
        n_confidence_bins: Number of equal-width bins (default 10) over
            the [0, 100] confidence range used for ECE and the
            reliability histogram.

    Returns:
        :class:`EvaluationResult`.
    """
    if prediction_column not in df.columns:
        msg = f"prediction column '{prediction_column}' missing"
        raise ValueError(msg)
    if truth_column not in df.columns:
        msg = f"truth column '{truth_column}' missing"
        raise ValueError(msg)

    n_rows = len(df)
    if n_rows == 0:
        msg = "evaluation DataFrame is empty"
        raise ValueError(msg)

    pred = df[prediction_column].to_numpy()
    truth = df[truth_column].to_numpy()
    correct_mask = pred == truth

    categories = sorted(set(truth.tolist()) | set(pred.tolist()))

    per_category: list[CategoryStats] = []
    confusion: dict[str, dict[str, int]] = {c: dict.fromkeys(categories, 0) for c in categories}
    for t, p in zip(truth, pred, strict=True):
        confusion[t][p] += 1

    for cat in categories:
        support = int((truth == cat).sum())
        predicted = int((pred == cat).sum())
        tp = int(((pred == cat) & (truth == cat)).sum())
        precision = tp / predicted if predicted > 0 else 0.0
        recall = tp / support if support > 0 else 0.0
        per_category.append(
            CategoryStats(
                category=cat,
                support=support,
                predicted=predicted,
                correct=tp,
                precision=precision,
                recall=recall,
            )
        )

    overall_accuracy = float(correct_mask.mean())

    ece: float | None = None
    bin_centers: list[float] = []
    bin_accs: list[float] = []
    bin_counts: list[int] = []
    if confidence_column is not None and confidence_column in df.columns:
        conf = df[confidence_column].to_numpy().astype(float)
        # Clip out-of-range confidences so binning is well-defined.
        conf = np.clip(conf, 0.0, 100.0)
        edges = np.linspace(0.0, 100.0, n_confidence_bins + 1)
        ece_sum = 0.0
        for i in range(n_confidence_bins):
            lo, hi = edges[i], edges[i + 1]
            # Right edge is inclusive only for the last bin so 100.0 lands somewhere.
            if i == n_confidence_bins - 1:
                mask = (conf >= lo) & (conf <= hi)
            else:
                mask = (conf >= lo) & (conf < hi)
            count = int(mask.sum())
            if count == 0:
                bin_acc = 0.0
                gap = 0.0
            else:
                bin_acc = float(correct_mask[mask].mean())
                # ECE uses the mean predicted confidence in the bin
                # (which may differ from the bin centre).
                bin_conf = float(conf[mask].mean()) / 100.0
                gap = abs(bin_acc - bin_conf)
            bin_centers.append(float((lo + hi) / 2))
            bin_accs.append(bin_acc)
            bin_counts.append(count)
            ece_sum += (count / n_rows) * gap
        ece = ece_sum

    return EvaluationResult(
        n_rows=n_rows,
        overall_accuracy=overall_accuracy,
        per_category=per_category,
        confusion=confusion,
        expected_calibration_error=ece,
        confidence_bin_centers=bin_centers,
        confidence_bin_accuracy=bin_accs,
        confidence_bin_counts=bin_counts,
    )
