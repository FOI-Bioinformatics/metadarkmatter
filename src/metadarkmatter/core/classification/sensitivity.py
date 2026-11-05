"""
Threshold sensitivity analysis.

Evaluates how classification counts change across a range of threshold
values, helping users assess whether results are robust or threshold-dependent.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import polars as pl

from metadarkmatter.core.classification.thresholds import apply_classification_thresholds
from metadarkmatter.models.config import ScoringConfig


@dataclass
class SensitivityResult:
    """Results of a threshold sensitivity sweep."""

    novelty_thresholds: list[float] = field(default_factory=list)
    uncertainty_thresholds: list[float] = field(default_factory=list)
    counts: dict[str, list[int]] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Convert to JSON-serializable dictionary."""
        return {
            "category_column": "legacy_call",
            "note": (
                "Counts are for the legacy hard-threshold cascade (legacy_call), "
                "not the Bayesian taxonomic_call. Novelty and uncertainty "
                "thresholds are swept together along a diagonal (point i of "
                "each list is used together)."
            ),
            "novelty_thresholds": [round(t, 2) for t in self.novelty_thresholds],
            "uncertainty_thresholds": [round(t, 2) for t in self.uncertainty_thresholds],
            "counts": self.counts,
        }


def run_sensitivity_analysis(
    df: pl.DataFrame,
    base_config: ScoringConfig,
    novelty_range: tuple[float, float] = (2.0, 8.0),
    uncertainty_range: tuple[float, float] = (0.5, 3.0),
    steps: int = 9,
) -> SensitivityResult:
    """
    Run threshold sensitivity analysis.

    Sweeps the novelty_known_max threshold across a range, re-classifying
    reads at each point with the legacy hard-threshold cascade (reported as
    ``legacy_call``). Novelty and uncertainty thresholds are swept together
    along a diagonal (not a 2D grid).
    This reveals whether classification counts are stable or threshold-sensitive.

    The DataFrame must contain pre-computed metrics from the classification
    pipeline: novelty_index, placement_uncertainty, num_ambiguous_hits,
    identity_gap.

    Args:
        df: Pre-computed metrics DataFrame (from classification pipeline).
        base_config: Base scoring configuration to vary.
        novelty_range: (min, max) for novelty_known_max sweep.
        uncertainty_range: (min, max) for uncertainty_known_max sweep.
        steps: Number of threshold points to evaluate.

    Returns:
        SensitivityResult with counts at each threshold.
    """
    import numpy as np

    result = SensitivityResult()

    novelty_values = np.linspace(novelty_range[0], novelty_range[1], steps).tolist()
    uncertainty_values = np.linspace(uncertainty_range[0], uncertainty_range[1], steps).tolist()

    result.novelty_thresholds = novelty_values
    result.uncertainty_thresholds = uncertainty_values

    # Initialize counts for each category
    categories = [
        "Known Species", "Novel Species", "Novel Genus",
        "Ambiguous", "Conserved Region", "Species Boundary",
        "Unclassified",
    ]
    # Any other category present (Ambiguous Within Genus, Off-target, ...) is
    # appended below so that per-point counts sum to the total read count.
    per_point: list[dict[str, int]] = []

    for novelty_thresh, uncertainty_thresh in zip(novelty_values, uncertainty_values, strict=False):
        # Create modified config with swept thresholds
        # Maintain continuous boundaries: species_min = known_max
        # Start from the full base config (keeps alignment_mode, genus and
        # identity-gap settings, ...); only the swept thresholds change.
        data = base_config.model_dump()
        data.update(
            novelty_known_max=novelty_thresh,
            novelty_novel_species_min=novelty_thresh,  # Continuous boundary
            uncertainty_known_max=uncertainty_thresh,
            uncertainty_novel_species_max=uncertainty_thresh,
            uncertainty_novel_genus_max=uncertainty_thresh,
        )
        config = ScoringConfig.model_validate(data)

        # Re-classify with modified thresholds
        classified = apply_classification_thresholds(df, config)

        # Count categories (the cascade writes its result to taxonomic_call
        # here; it is the legacy cascade, reported as legacy_call)
        counts = classified.group_by("taxonomic_call").len()
        count_dict = {row[0]: row[1] for row in counts.iter_rows()}

        per_point.append(count_dict)

    extras = sorted({c for d in per_point for c in d if c not in categories and c is not None})
    for cat in [*categories, *extras]:
        result.counts[cat] = [d.get(cat, 0) for d in per_point]

    return result
