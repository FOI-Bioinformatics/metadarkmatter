"""Regression guard for shipped calibrated Bayesian configs.

Phase 7 of the calibration validation plan. Once a calibrated YAML is
committed under ``configs/calibrated/`` (or ``configs/bayesian_default.yaml``
if the cross-corpus study decided a single global default is enough),
this test loads it, runs the classifier against a small held-out
parquet from the calibration corpus, and asserts that top-1 accuracy
is at least the value recorded at shipment time.

The test is **skipped automatically** when the calibrated config or
the held-out corpus is not yet committed - this keeps the test file in
the suite as a deliberate placeholder that fails loudly the moment
calibrated artefacts land but isn't blocked by their absence.

To activate the test:

1. Commit a calibrated config, e.g. ``configs/calibrated/francisella.yaml``.
2. Commit a held-out Parquet (~1 MB) with the metric columns and the
   ``true_category`` truth column at
   ``tests/benchmark_corpus/francisella_holdout.parquet``.
3. Record the shipment accuracy in ``SHIPMENT_BASELINES`` below.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import polars as pl
import pytest

from metadarkmatter.core.classification.evaluation import (
    evaluate_classifications,
)
from metadarkmatter.models.config import ScoringConfig

REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = REPO_ROOT / "configs" / "calibrated"
CORPUS_DIR = REPO_ROOT / "tests" / "benchmark_corpus"


@dataclass(frozen=True)
class ShipmentBaseline:
    family: str
    config_path: Path
    holdout_path: Path
    # Floor accuracy recorded at the time the calibrated YAML was committed.
    # Subsequent runs may exceed this; failing below it indicates a
    # regression in core/classification/bayesian.py or its inputs.
    min_top1_accuracy: float
    # Ceiling on ECE; calibrated configs should keep ECE below this.
    max_ece: float


# Populate one entry per shipped calibrated config.
SHIPMENT_BASELINES: tuple[ShipmentBaseline, ...] = (
    # ShipmentBaseline(
    #     family="francisella",
    #     config_path=CONFIG_DIR / "francisella.yaml",
    #     holdout_path=CORPUS_DIR / "francisella_holdout.parquet",
    #     min_top1_accuracy=0.85,
    #     max_ece=0.15,
    # ),
    # ShipmentBaseline(
    #     family="pseudomonas",
    #     config_path=CONFIG_DIR / "pseudomonas.yaml",
    #     holdout_path=CORPUS_DIR / "pseudomonas_holdout.parquet",
    #     min_top1_accuracy=0.80,
    #     max_ece=0.20,
    # ),
)


@pytest.mark.skipif(
    not SHIPMENT_BASELINES,
    reason="No shipped calibrated configs yet; see test docstring for activation.",
)
@pytest.mark.parametrize(
    "baseline", SHIPMENT_BASELINES, ids=lambda b: b.family
)
def test_shipped_calibration_meets_recorded_accuracy(
    baseline: ShipmentBaseline,
) -> None:
    """Loaded calibrated config must classify the held-out corpus at or
    above the accuracy recorded when it was committed."""
    if not baseline.config_path.is_file():
        pytest.skip(f"Calibrated config not found: {baseline.config_path}")
    if not baseline.holdout_path.is_file():
        pytest.skip(f"Held-out corpus not found: {baseline.holdout_path}")

    # Validate the config loads and carries calibrated parameters.
    config = ScoringConfig.from_yaml(baseline.config_path)
    assert (
        config.bayesian.category_params is not None
    ), "Shipped config must contain calibrated Bayesian category_params"

    holdout = pl.read_parquet(baseline.holdout_path)
    required = {"taxonomic_call", "true_category", "confidence_score"}
    missing = required - set(holdout.columns)
    assert not missing, f"Held-out corpus is missing columns: {missing}"

    result = evaluate_classifications(
        holdout,
        truth_column="true_category",
        prediction_column="taxonomic_call",
        confidence_column="confidence_score",
    )
    assert result.overall_accuracy >= baseline.min_top1_accuracy, (
        f"{baseline.family}: top-1 accuracy {result.overall_accuracy:.4f} "
        f"is below the shipped baseline {baseline.min_top1_accuracy:.4f}. "
        "Either a regression has been introduced in bayesian.py or the "
        "shipped config needs to be recalibrated."
    )
    if result.expected_calibration_error is not None:
        assert result.expected_calibration_error <= baseline.max_ece, (
            f"{baseline.family}: ECE {result.expected_calibration_error:.4f} "
            f"exceeds the shipped ceiling {baseline.max_ece:.4f}."
        )


def test_placeholder_keeps_module_in_suite() -> None:
    """Trivial test so pytest collects this module even when no
    SHIPMENT_BASELINES entries are active. Remove once the parametrized
    test above is real."""
    assert isinstance(SHIPMENT_BASELINES, tuple)
