"""Tests for the entropy -> confidence calibration pipeline (Phase 1.3)."""

from __future__ import annotations

import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from metadarkmatter.core.classification.bayesian import entropy_to_confidence
from metadarkmatter.models.config import BayesianConfig, ScoringConfig

SCRIPTS_DIR = Path(__file__).resolve().parents[2] / "scripts"
MAX_ENTROPY_6 = math.log2(6)


# ---------------------------------------------------------------------------
# entropy_to_confidence: linear default vs calibrated map
# ---------------------------------------------------------------------------


def test_default_mapping_is_linear() -> None:
    assert entropy_to_confidence(0.0) == pytest.approx(100.0)
    assert entropy_to_confidence(MAX_ENTROPY_6) == pytest.approx(0.0)
    # Midpoint of the entropy range maps to the midpoint of confidence.
    assert entropy_to_confidence(MAX_ENTROPY_6 / 2) == pytest.approx(50.0)


def test_calibration_overrides_default() -> None:
    cal = [[0.0, 100.0], [1.0, 80.0], [2.0, 30.0], [MAX_ENTROPY_6, 0.0]]
    # At a knot the calibrated value wins, even though the default differs.
    assert entropy_to_confidence(1.0, calibration=cal) == pytest.approx(80.0)
    # Between knots we expect linear interpolation.
    assert entropy_to_confidence(0.5, calibration=cal) == pytest.approx(90.0)
    assert entropy_to_confidence(1.5, calibration=cal) == pytest.approx(55.0)


def test_calibration_clips_to_unit_range() -> None:
    cal = [[0.0, 120.0], [MAX_ENTROPY_6, -10.0]]
    result = entropy_to_confidence(np.array([0.0, MAX_ENTROPY_6 / 2, MAX_ENTROPY_6]), calibration=cal)
    assert np.all(result >= 0.0) and np.all(result <= 100.0)


def test_calibration_handles_unsorted_input() -> None:
    cal = [[2.0, 30.0], [0.0, 100.0], [1.0, 80.0]]  # deliberately scrambled
    # The function must sort defensively, so this returns the same as the sorted version.
    assert entropy_to_confidence(0.5, calibration=cal) == pytest.approx(90.0)


def test_calibration_supports_array_input() -> None:
    cal = [[0.0, 100.0], [1.0, 80.0], [MAX_ENTROPY_6, 0.0]]
    out = entropy_to_confidence(np.array([0.0, 1.0, MAX_ENTROPY_6]), calibration=cal)
    assert out.shape == (3,)
    assert np.allclose(out, [100.0, 80.0, 0.0])


# ---------------------------------------------------------------------------
# Config round-trip
# ---------------------------------------------------------------------------


def test_yaml_round_trip_preserves_entropy_calibration(tmp_path: Path) -> None:
    cal = [[0.0, 100.0], [1.5, 50.0], [MAX_ENTROPY_6, 5.0]]
    cfg = ScoringConfig(bayesian=BayesianConfig(entropy_calibration=cal))
    out = tmp_path / "cal.yaml"
    cfg.to_yaml(out)
    reloaded = ScoringConfig.from_yaml(out)
    assert reloaded.bayesian.entropy_calibration == cal


# ---------------------------------------------------------------------------
# Script: calibrate_entropy.py
# ---------------------------------------------------------------------------


def _run_script(script: str, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / script), *args],
        capture_output=True,
        text=True,
        check=False,
    )


def test_calibrate_entropy_produces_monotone_decreasing_knots(tmp_path: Path) -> None:
    bench = tmp_path / "bench.tsv"
    yaml_out = tmp_path / "cal.yaml"

    r1 = _run_script(
        "build_synthetic_benchmark.py",
        "--output", str(bench), "--per-category", "200", "--seed", "5",
    )
    assert r1.returncode == 0, r1.stderr

    r2 = _run_script(
        "calibrate_entropy.py",
        "--benchmark", str(bench),
        "--output", str(yaml_out),
        "--bins", "15",
    )
    assert r2.returncode == 0, r2.stderr

    cfg = ScoringConfig.from_yaml(yaml_out)
    assert cfg.bayesian.entropy_calibration is not None
    assert len(cfg.bayesian.entropy_calibration) == 15

    # Knots must be (a) ordered in entropy and (b) non-increasing in confidence.
    confidences = [c for _, c in cfg.bayesian.entropy_calibration]
    entropies = [e for e, _ in cfg.bayesian.entropy_calibration]
    assert entropies == sorted(entropies)
    for prev, nxt in zip(confidences, confidences[1:], strict=False):
        assert nxt <= prev + 1e-9, (
            f"calibration is not monotone non-increasing: {confidences}"
        )
    # Confidences must stay within [0, 100].
    assert all(0.0 <= c <= 100.0 for c in confidences)
