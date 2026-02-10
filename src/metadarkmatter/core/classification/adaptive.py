"""
Adaptive threshold detection from ANI distribution.

Uses Gaussian Mixture Model to detect natural species boundary
from the ANI matrix distribution, rather than assuming a fixed
threshold (e.g., 95-96% ANI) for all taxa.

Requires scikit-learn (optional dependency: pip install metadarkmatter[adaptive]).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from metadarkmatter.core.classification.ani_matrix import ANIMatrix
    from metadarkmatter.models.config import ScoringConfig

logger = logging.getLogger(__name__)


@dataclass
class AdaptiveThresholds:
    """Result of adaptive threshold detection."""

    species_boundary: float  # Detected ANI boundary (e.g., 94.2%)
    novelty_known_max: float  # = 100 - species_boundary
    confidence: float  # GMM separation quality (0-1)
    method: str  # "gmm" or "fallback"
    ani_values_used: int  # Number of pairwise ANI values analyzed


def detect_species_boundary(
    ani_matrix: ANIMatrix,
    min_genomes: int = 5,
) -> AdaptiveThresholds:
    """
    Detect natural species boundary from ANI matrix distribution.

    Fits a 2-component Gaussian Mixture Model to upper-triangle ANI
    values. The two components typically correspond to:
    - Within-species comparisons (high ANI, narrow distribution)
    - Between-species comparisons (lower ANI, broader distribution)

    The species boundary is estimated as the midpoint between the
    two component means, weighted by their standard deviations.

    Falls back to the default threshold (96% ANI = 4% novelty) if:
    - Too few genomes (< min_genomes)
    - GMM doesn't converge
    - Components are not well separated

    Args:
        ani_matrix: ANI matrix with pairwise genome comparisons.
        min_genomes: Minimum number of genomes required for GMM.

    Returns:
        AdaptiveThresholds with detected or fallback boundary.
    """
    from metadarkmatter.core.constants import NOVELTY_KNOWN_MAX

    default_boundary = 100.0 - NOVELTY_KNOWN_MAX  # 96.0%
    genomes = sorted(ani_matrix.genomes)  # Convert set to sorted list for indexing

    if len(genomes) < min_genomes:
        logger.info(
            "Too few genomes (%d < %d) for adaptive threshold detection. "
            "Using default boundary %.1f%%.",
            len(genomes), min_genomes, default_boundary,
        )
        return AdaptiveThresholds(
            species_boundary=default_boundary,
            novelty_known_max=NOVELTY_KNOWN_MAX,
            confidence=0.0,
            method="fallback",
            ani_values_used=0,
        )

    # Extract upper-triangle ANI values (excluding diagonal)
    ani_values = []
    for i, g1 in enumerate(genomes):
        for j in range(i + 1, len(genomes)):
            g2 = genomes[j]
            ani = ani_matrix.get_ani(g1, g2)
            if ani > 0:
                ani_values.append(ani)

    if len(ani_values) < 3:
        logger.info(
            "Insufficient pairwise ANI values (%d). Using default boundary.",
            len(ani_values),
        )
        return AdaptiveThresholds(
            species_boundary=default_boundary,
            novelty_known_max=NOVELTY_KNOWN_MAX,
            confidence=0.0,
            method="fallback",
            ani_values_used=len(ani_values),
        )

    ani_array = np.array(ani_values).reshape(-1, 1)

    try:
        from sklearn.mixture import GaussianMixture

        gmm = GaussianMixture(
            n_components=2,
            covariance_type="full",
            max_iter=200,
            n_init=5,
            random_state=42,
        )
        gmm.fit(ani_array)

        means = gmm.means_.flatten()
        stds = np.sqrt(gmm.covariances_.flatten())

        # Order components: low ANI first, high ANI second
        order = np.argsort(means)
        low_mean, high_mean = means[order]
        low_std, high_std = stds[order]

        # Separation quality: how distinct are the two components?
        # Using the standardized gap between means
        separation = abs(high_mean - low_mean) / (low_std + high_std + 1e-10)

        if separation < 1.0:
            # Components not well separated - fallback
            logger.info(
                "GMM components not well separated (separation=%.2f). "
                "Using default boundary %.1f%%.",
                separation, default_boundary,
            )
            return AdaptiveThresholds(
                species_boundary=default_boundary,
                novelty_known_max=NOVELTY_KNOWN_MAX,
                confidence=separation / 3.0,  # Normalize to ~0-1
                method="fallback",
                ani_values_used=len(ani_values),
            )

        # Species boundary: weighted midpoint between components
        # Weight by inverse of standard deviations (tighter distribution = more influence)
        w_low = 1.0 / (low_std + 1e-10)
        w_high = 1.0 / (high_std + 1e-10)
        boundary = (low_mean * w_low + high_mean * w_high) / (w_low + w_high)

        # Clamp to reasonable range (80-99% ANI)
        boundary = max(80.0, min(99.0, boundary))
        novelty_max = 100.0 - boundary

        confidence = min(1.0, separation / 3.0)

        logger.info(
            "Detected species boundary at %.1f%% ANI (novelty_known_max=%.1f%%, "
            "confidence=%.2f, method=gmm, n=%d).",
            boundary, novelty_max, confidence, len(ani_values),
        )

        return AdaptiveThresholds(
            species_boundary=boundary,
            novelty_known_max=novelty_max,
            confidence=confidence,
            method="gmm",
            ani_values_used=len(ani_values),
        )

    except ImportError:
        logger.warning(
            "scikit-learn not available. Install with: pip install metadarkmatter[adaptive]"
        )
        return AdaptiveThresholds(
            species_boundary=default_boundary,
            novelty_known_max=NOVELTY_KNOWN_MAX,
            confidence=0.0,
            method="fallback",
            ani_values_used=len(ani_values),
        )
    except Exception as e:
        logger.warning(
            "GMM fitting failed: %s. Using default boundary.", e,
        )
        return AdaptiveThresholds(
            species_boundary=default_boundary,
            novelty_known_max=NOVELTY_KNOWN_MAX,
            confidence=0.0,
            method="fallback",
            ani_values_used=len(ani_values),
        )


def build_adaptive_config(
    base_config: ScoringConfig,
    adaptive: AdaptiveThresholds,
) -> ScoringConfig:
    """
    Build a ScoringConfig with adaptive thresholds applied.

    Replaces the novelty_known_max and novelty_novel_species_min
    with the detected boundary. Other thresholds remain unchanged.

    Args:
        base_config: Base scoring configuration.
        adaptive: Detected adaptive thresholds.

    Returns:
        New ScoringConfig with adaptive novelty boundary.
    """
    from metadarkmatter.models.config import ScoringConfig

    return ScoringConfig(
        bitscore_threshold_pct=base_config.bitscore_threshold_pct,
        novelty_known_max=adaptive.novelty_known_max,
        novelty_novel_species_min=adaptive.novelty_known_max,  # Continuous boundary
        novelty_novel_species_max=base_config.novelty_novel_species_max,
        novelty_novel_genus_min=base_config.novelty_novel_genus_min,
        novelty_novel_genus_max=base_config.novelty_novel_genus_max,
        uncertainty_known_max=base_config.uncertainty_known_max,
        uncertainty_novel_species_max=base_config.uncertainty_novel_species_max,
        uncertainty_novel_genus_max=base_config.uncertainty_novel_genus_max,
        uncertainty_conserved_min=base_config.uncertainty_conserved_min,
        uncertainty_mode=base_config.uncertainty_mode,
        coverage_weight_mode=base_config.coverage_weight_mode,
        coverage_weight_strength=base_config.coverage_weight_strength,
        min_alignment_fraction=base_config.min_alignment_fraction,
        single_hit_uncertainty_threshold=base_config.single_hit_uncertainty_threshold,
    )
