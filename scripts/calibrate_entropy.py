#!/usr/bin/env python
"""
Fit an empirical entropy -> confidence calibration.

Reads a labelled corpus, computes the Bayesian posterior entropy for
every row using the current (possibly previously calibrated) Gaussian
parameters, and then fits a monotone-decreasing isotonic regression
of accuracy on entropy. The fitted mapping is written into the
``bayesian.entropy_calibration`` field of a YAML scoring config, where
it replaces the default linear ``(1 - entropy / log2(6)) * 100`` map.

Requires scikit-learn for the isotonic fit (already in the
``[adaptive]`` optional dependency group). If scikit-learn is
unavailable, falls back to per-bin averaging.

Usage:

    python scripts/calibrate_entropy.py \
        --benchmark benchmarks/synthetic.tsv \
        --base-config configs/bayesian_calibrated.yaml \
        --output configs/bayesian_calibrated.yaml \
        --bins 20
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import polars as pl

from metadarkmatter.core.classification.bayesian import (
    _CATEGORY_NAMES_6,
    _MAX_ENTROPY_6,
    _shannon_entropy,
    build_category_params,
)
from metadarkmatter.models.config import ScoringConfig


def _posterior_entropy_and_pred(
    df: pl.DataFrame, config: ScoringConfig
) -> tuple[np.ndarray, np.ndarray]:
    """Return posterior entropy and predicted category per row.

    Re-implements just the likelihood + posterior steps so this script
    has no dependency on the full BayesianClassifier wiring.
    """
    params = build_category_params(config)
    priors = np.array(
        [config.bayesian.priors[k.lower().replace(" ", "_")] for k in _CATEGORY_NAMES_6]
    )

    n = len(df)
    log_lik = np.zeros((n, len(params)))
    novelty = df["novelty_index"].to_numpy()
    uncertainty = df["placement_uncertainty"].to_numpy()
    for j, p in enumerate(params):
        nz = (novelty - p.novelty_mean) / max(p.novelty_sigma, 1e-6)
        uz = (uncertainty - p.uncertainty_mean) / max(p.uncertainty_sigma, 1e-6)
        log_norm = -np.log(p.novelty_sigma * p.uncertainty_sigma)
        log_lik[:, j] = -0.5 * (nz * nz + uz * uz) + log_norm

    log_post = log_lik + np.log(priors + 1e-12)
    log_post -= log_post.max(axis=1, keepdims=True)
    post = np.exp(log_post)
    post /= post.sum(axis=1, keepdims=True)

    entropy = _shannon_entropy(post)
    pred = np.array([_CATEGORY_NAMES_6[i] for i in post.argmax(axis=1)])
    return entropy, pred


def _fit_isotonic(
    entropy: np.ndarray, correct: np.ndarray, bins: int
) -> list[list[float]]:
    """Return [[entropy, confidence], ...] knots from an isotonic fit.

    Uses sklearn's IsotonicRegression with increasing=False on accuracy
    vs entropy, then scales to a 0-100 confidence and downsamples to
    ``bins`` knots for compact serialization.
    """
    try:
        from sklearn.isotonic import IsotonicRegression
    except ImportError:
        return _fit_per_bin_average(entropy, correct, bins)

    iso = IsotonicRegression(increasing=False, y_min=0.0, y_max=1.0)
    iso.fit(entropy, correct.astype(float))

    grid = np.linspace(0.0, _MAX_ENTROPY_6, bins)
    confs = iso.predict(grid) * 100.0
    return [[float(e), float(c)] for e, c in zip(grid, confs, strict=True)]


def _fit_per_bin_average(
    entropy: np.ndarray, correct: np.ndarray, bins: int
) -> list[list[float]]:
    """Fallback fit: per-bin mean accuracy, enforced monotone by cummax-from-right."""
    edges = np.linspace(0.0, _MAX_ENTROPY_6, bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    confs = np.zeros(bins)
    for i in range(bins):
        mask = (entropy >= edges[i]) & (entropy < edges[i + 1])
        confs[i] = correct[mask].mean() * 100.0 if mask.any() else np.nan

    # Forward-fill NaNs, then enforce non-increasing in entropy.
    last = 100.0
    for i in range(bins):
        if np.isnan(confs[i]):
            confs[i] = last
        else:
            last = confs[i]
    # Non-increasing in entropy: as index grows, confidence cannot rise.
    confs = np.minimum.accumulate(confs)
    return [[float(e), float(c)] for e, c in zip(centers, confs, strict=True)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", required=True, type=Path)
    parser.add_argument(
        "--base-config",
        type=Path,
        default=None,
        help="Optional starting YAML (e.g. one from calibrate_bayesian.py).",
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--bins",
        type=int,
        default=20,
        help="Number of knots in the saved entropy -> confidence mapping (default 20).",
    )
    args = parser.parse_args()

    df = pl.read_csv(args.benchmark, separator="\t")
    required = {"novelty_index", "placement_uncertainty", "true_category"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Benchmark TSV missing columns: {sorted(missing)}")

    base_config = (
        ScoringConfig.from_yaml(args.base_config)
        if args.base_config is not None
        else ScoringConfig()
    )

    entropy, pred = _posterior_entropy_and_pred(df, base_config)
    correct = (pred == df["true_category"].to_numpy()).astype(int)

    print(f"Fitting on {len(df)} rows; raw accuracy = {correct.mean() * 100:.1f}%")
    knots = _fit_isotonic(entropy, correct, args.bins)

    # Show a compact preview
    print("Calibration knots:")
    for e, c in knots[:: max(1, len(knots) // 5)]:
        print(f"  entropy={e:5.3f}  -> confidence={c:6.2f}")

    new_bay = base_config.bayesian.model_copy(update={"entropy_calibration": knots})
    new_config = base_config.model_copy(update={"bayesian": new_bay})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    new_config.to_yaml(args.output)
    print(f"\nWrote calibrated config to {args.output}")


if __name__ == "__main__":
    main()
