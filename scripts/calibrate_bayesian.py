#!/usr/bin/env python
"""
Fit per-category Bayesian 2D Gaussian parameters to a labelled corpus.

Reads a TSV containing ``novelty_index``, ``placement_uncertainty``, and
``true_category`` columns and produces a YAML scoring config with
calibrated ``bayesian.category_params``. The output is consumable by
``metadarkmatter score classify --config <path>``.

The script also reports the empirical Pearson correlation between
novelty and uncertainty for each category, so users can judge whether
the diagonal-covariance Gaussian assumption is justified or whether a
follow-up full-covariance model would be warranted.

Usage:

    python scripts/calibrate_bayesian.py \
        --benchmark benchmarks/synthetic.tsv \
        --output configs/bayesian_calibrated.yaml \
        --min-samples 20
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import polars as pl

from metadarkmatter.models.config import BayesianConfig, ScoringConfig

EXPECTED_CATEGORIES = [
    "Known Species",
    "Novel Species",
    "Novel Genus",
    "Species Boundary",
    "Ambiguous",
    "Unclassified",
]


def _fit_category(
    df: pl.DataFrame,
    category: str,
    min_samples: int,
    sigma_floor: float,
) -> tuple[dict[str, float] | None, float, int]:
    """Fit (mean_n, mean_u, sigma_n, sigma_u) for one category.

    Returns the parameter dict (or None if too few samples), the empirical
    correlation between novelty and uncertainty, and the sample count.
    """
    sub = df.filter(pl.col("true_category") == category)
    n = len(sub)
    if n < min_samples:
        return None, float("nan"), n

    novelty = sub["novelty_index"].to_numpy()
    uncertainty = sub["placement_uncertainty"].to_numpy()

    params = {
        "novelty_mean": float(np.mean(novelty)),
        "uncertainty_mean": float(np.mean(uncertainty)),
        "novelty_sigma": float(max(np.std(novelty, ddof=1), sigma_floor)),
        "uncertainty_sigma": float(max(np.std(uncertainty, ddof=1), sigma_floor)),
    }
    if np.std(novelty) > 0 and np.std(uncertainty) > 0:
        corr = float(np.corrcoef(novelty, uncertainty)[0, 1])
    else:
        corr = 0.0
    return params, corr, n


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--benchmark",
        required=True,
        type=Path,
        help="Labelled TSV with novelty_index, placement_uncertainty, true_category.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output YAML config path.",
    )
    parser.add_argument(
        "--min-samples",
        type=int,
        default=20,
        help="Minimum per-category sample count to attempt a fit (default 20).",
    )
    parser.add_argument(
        "--sigma-floor",
        type=float,
        default=0.5,
        help="Lower bound on fitted sigmas (default 0.5).",
    )
    parser.add_argument(
        "--base-config",
        type=Path,
        default=None,
        help=(
            "Optional base YAML to inherit non-Bayesian settings from. "
            "Defaults to the package defaults."
        ),
    )
    parser.add_argument(
        "--corr-threshold",
        type=float,
        default=0.3,
        help=(
            "Absolute correlation above which to warn that the diagonal-"
            "covariance assumption may be too restrictive (default 0.3)."
        ),
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

    fitted: dict[str, dict[str, float]] = {}
    print(f"{'Category':<20} {'N':>6}  {'mu_N':>7} {'mu_U':>7} "
          f"{'sg_N':>6} {'sg_U':>6}  {'r(N,U)':>7}")
    print("-" * 72)
    for category in EXPECTED_CATEGORIES:
        params, corr, n = _fit_category(df, category, args.min_samples, args.sigma_floor)
        if params is None:
            print(f"{category:<20} {n:>6}  (skipped: below --min-samples)")
            continue
        fitted[category] = params
        flag = " !" if abs(corr) > args.corr_threshold else ""
        print(
            f"{category:<20} {n:>6}  "
            f"{params['novelty_mean']:>7.2f} {params['uncertainty_mean']:>7.2f} "
            f"{params['novelty_sigma']:>6.2f} {params['uncertainty_sigma']:>6.2f}  "
            f"{corr:>+7.3f}{flag}"
        )
    if not fitted:
        raise SystemExit(
            "No categories had enough samples to fit. Increase the corpus "
            "or lower --min-samples."
        )

    new_bayesian = base_config.bayesian.model_copy(
        update={"category_params": fitted}
    )
    new_config = base_config.model_copy(update={"bayesian": new_bayesian})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    new_config.to_yaml(args.output)

    print(
        f"\nWrote calibrated config to {args.output}\n"
        f"Categories fitted: {len(fitted)}/{len(EXPECTED_CATEGORIES)}\n"
        "Lines marked ' !' have |r(N,U)| above the threshold; consider a "
        "full-covariance model for those categories."
    )


if __name__ == "__main__":
    main()
