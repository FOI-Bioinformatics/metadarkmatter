"""
Threshold sensitivity analysis.

Evaluates how classification counts change across a range of threshold
values, helping users assess whether results are robust or threshold-dependent.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import polars as pl

from metadarkmatter.core.classification.thresholds import apply_classification_thresholds
from metadarkmatter.models.config import ScoringConfig


@dataclass
class SensitivityResult:
    """Results of a threshold sensitivity sweep."""

    novelty_thresholds: list[float] = field(default_factory=list)
    uncertainty_thresholds: list[float] = field(default_factory=list)
    counts: dict[str, list[int]] = field(default_factory=dict)

    def to_dict(self) -> dict:
        """Convert to JSON-serializable dictionary."""
        return {
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
    reads at each point. The uncertainty thresholds are co-varied proportionally.
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
    for cat in categories:
        result.counts[cat] = []

    for novelty_thresh, uncertainty_thresh in zip(novelty_values, uncertainty_values):
        # Create modified config with swept thresholds
        # Maintain continuous boundaries: species_min = known_max
        config = ScoringConfig(
            bitscore_threshold_pct=base_config.bitscore_threshold_pct,
            novelty_known_max=novelty_thresh,
            novelty_novel_species_min=novelty_thresh,  # Continuous boundary
            novelty_novel_species_max=base_config.novelty_novel_species_max,
            novelty_novel_genus_min=base_config.novelty_novel_genus_min,
            novelty_novel_genus_max=base_config.novelty_novel_genus_max,
            uncertainty_known_max=uncertainty_thresh,
            uncertainty_novel_species_max=uncertainty_thresh,
            uncertainty_novel_genus_max=uncertainty_thresh,
            uncertainty_conserved_min=base_config.uncertainty_conserved_min,
            uncertainty_mode=base_config.uncertainty_mode,
            coverage_weight_mode=base_config.coverage_weight_mode,
            coverage_weight_strength=base_config.coverage_weight_strength,
            min_alignment_fraction=base_config.min_alignment_fraction,
            single_hit_uncertainty_threshold=base_config.single_hit_uncertainty_threshold,
        )

        # Re-classify with modified thresholds
        classified = apply_classification_thresholds(df, config)

        # Count categories
        counts = classified.group_by("taxonomic_call").len()
        count_dict = {row[0]: row[1] for row in counts.iter_rows()}

        for cat in categories:
            result.counts[cat].append(count_dict.get(cat, 0))

    return result
