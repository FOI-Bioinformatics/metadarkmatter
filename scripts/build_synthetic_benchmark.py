#!/usr/bin/env python
"""
Generate a synthetic labelled benchmark for Bayesian calibration.

Produces a TSV with columns ``novelty_index``, ``placement_uncertainty``
and ``true_category`` by sampling from per-category 2D Gaussians around
the hand-tuned defaults baked into the codebase. This is useful as a
smoke-testing corpus for the calibration pipeline; it is **not** a
substitute for a real labelled dataset.

Usage:

    python scripts/build_synthetic_benchmark.py \
        --output benchmarks/synthetic.tsv \
        --per-category 500 \
        --seed 42
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import polars as pl

from metadarkmatter.core.classification.bayesian import build_category_params
from metadarkmatter.models.config import ScoringConfig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output TSV path.",
    )
    parser.add_argument(
        "--per-category",
        type=int,
        default=500,
        help="Number of samples per category (default 500).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="RNG seed (default 42).",
    )
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    config = ScoringConfig()
    params = build_category_params(config)

    rows: list[dict[str, object]] = []
    for cp in params:
        n_samples = args.per_category
        novelty = rng.normal(cp.novelty_mean, cp.novelty_sigma, n_samples)
        uncertainty = rng.normal(cp.uncertainty_mean, cp.uncertainty_sigma, n_samples)
        # Clip to physically meaningful ranges (0-100 for both axes).
        novelty = np.clip(novelty, 0.0, 100.0)
        uncertainty = np.clip(uncertainty, 0.0, 100.0)
        for n_val, u_val in zip(novelty, uncertainty, strict=True):
            rows.append(
                {
                    "novelty_index": float(n_val),
                    "placement_uncertainty": float(u_val),
                    "true_category": cp.name,
                }
            )

    df = pl.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(args.output, separator="\t")
    print(
        f"Wrote {len(df)} rows ({args.per_category} per category) to {args.output}"
    )


if __name__ == "__main__":
    main()
