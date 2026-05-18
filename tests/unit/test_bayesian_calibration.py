"""Tests for the Bayesian calibration pipeline (Phase 1.1-1.2)."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from metadarkmatter.core.classification.bayesian import build_category_params
from metadarkmatter.models.config import BayesianConfig, ScoringConfig

SCRIPTS_DIR = Path(__file__).resolve().parents[2] / "scripts"


# ---------------------------------------------------------------------------
# Config-layer: fitted parameters override defaults
# ---------------------------------------------------------------------------


def test_default_bayesian_config_has_no_category_params() -> None:
    cfg = ScoringConfig()
    assert cfg.bayesian.category_params is None


def test_fitted_category_params_override_defaults() -> None:
    fitted = {
        "Known Species": {
            "novelty_mean": 1.23,
            "uncertainty_mean": 0.45,
            "novelty_sigma": 2.0,
            "uncertainty_sigma": 1.0,
        }
    }
    cfg = ScoringConfig(bayesian=BayesianConfig(category_params=fitted))
    params = build_category_params(cfg)
    by_name = {p.name: p for p in params}
    assert by_name["Known Species"].novelty_mean == pytest.approx(1.23)
    assert by_name["Known Species"].uncertainty_mean == pytest.approx(0.45)
    assert by_name["Known Species"].novelty_sigma == pytest.approx(2.0)
    assert by_name["Known Species"].uncertainty_sigma == pytest.approx(1.0)
    # A category not in the fitted dict must keep its hand-tuned default.
    assert by_name["Novel Species"].novelty_mean != pytest.approx(1.23)


def test_yaml_round_trip_preserves_category_params(tmp_path: Path) -> None:
    fitted = {
        "Novel Genus": {
            "novelty_mean": 22.5,
            "uncertainty_mean": 1.0,
            "novelty_sigma": 1.5,
            "uncertainty_sigma": 0.7,
        }
    }
    cfg = ScoringConfig(bayesian=BayesianConfig(category_params=fitted))
    yaml_path = tmp_path / "out.yaml"
    cfg.to_yaml(yaml_path)
    reloaded = ScoringConfig.from_yaml(yaml_path)
    assert reloaded.bayesian.category_params == fitted


# ---------------------------------------------------------------------------
# Script-layer: build_synthetic_benchmark.py + calibrate_bayesian.py
# ---------------------------------------------------------------------------


def _run_script(script: str, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / script), *args],
        capture_output=True,
        text=True,
        check=False,
    )


def test_build_synthetic_benchmark_produces_labeled_tsv(tmp_path: Path) -> None:
    out = tmp_path / "bench.tsv"
    result = _run_script(
        "build_synthetic_benchmark.py",
        "--output", str(out),
        "--per-category", "50",
        "--seed", "3",
    )
    assert result.returncode == 0, result.stderr
    df = pl.read_csv(out, separator="\t")
    assert set(df.columns) == {"novelty_index", "placement_uncertainty", "true_category"}
    # 6 categories x 50 = 300 rows
    assert len(df) == 300
    assert set(df["true_category"].unique()) == {
        "Known Species", "Novel Species", "Novel Genus",
        "Species Boundary", "Ambiguous", "Unclassified",
    }


def test_calibrate_bayesian_writes_loadable_yaml(tmp_path: Path) -> None:
    bench = tmp_path / "bench.tsv"
    yaml_out = tmp_path / "cal.yaml"

    r1 = _run_script(
        "build_synthetic_benchmark.py",
        "--output", str(bench),
        "--per-category", "100",
        "--seed", "11",
    )
    assert r1.returncode == 0, r1.stderr

    r2 = _run_script(
        "calibrate_bayesian.py",
        "--benchmark", str(bench),
        "--output", str(yaml_out),
        "--min-samples", "20",
    )
    assert r2.returncode == 0, r2.stderr

    cfg = ScoringConfig.from_yaml(yaml_out)
    assert cfg.bayesian.category_params is not None
    assert len(cfg.bayesian.category_params) == 6
    # Sanity: every category dict has the four required keys
    for spec in cfg.bayesian.category_params.values():
        assert set(spec.keys()) >= {
            "novelty_mean", "uncertainty_mean",
            "novelty_sigma", "uncertainty_sigma",
        }


def test_calibration_skips_undersampled_categories(tmp_path: Path) -> None:
    # Build a benchmark with only Known Species rows; other categories
    # should be skipped, not crash.
    rng = np.random.default_rng(0)
    df = pl.DataFrame(
        {
            "novelty_index": rng.normal(2.0, 1.0, 50).tolist(),
            "placement_uncertainty": rng.normal(0.5, 0.3, 50).tolist(),
            "true_category": ["Known Species"] * 50,
        }
    )
    bench = tmp_path / "single_cat.tsv"
    df.write_csv(bench, separator="\t")

    yaml_out = tmp_path / "cal.yaml"
    result = _run_script(
        "calibrate_bayesian.py",
        "--benchmark", str(bench),
        "--output", str(yaml_out),
        "--min-samples", "20",
    )
    assert result.returncode == 0, result.stderr

    cfg = ScoringConfig.from_yaml(yaml_out)
    assert cfg.bayesian.category_params is not None
    assert "Known Species" in cfg.bayesian.category_params
    # Categories with no rows are simply absent from the output.
    assert "Novel Species" not in cfg.bayesian.category_params
