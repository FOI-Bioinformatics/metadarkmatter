"""
Adaptive threshold detection from ANI distribution.

Uses Gaussian Mixture Models to detect natural species and genus boundaries
from the ANI matrix distribution, rather than assuming fixed thresholds
(e.g., 95-96% ANI for species, 80% ANI for genera).

Species boundary: 2-component GMM separating within-species from between-species pairs.
Genus boundary: 3-component GMM separating within-species, between-species-within-genus,
and between-genus pairs.

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


@dataclass
class AdaptiveGenusThreshold:
    """Result of adaptive genus boundary detection."""

    genus_boundary: float  # Detected ANI boundary (e.g., 82%)
    novelty_genus_min: float  # = 100 - genus_boundary
    confidence: float  # Separation quality (0-1)
    method: str  # "gmm_3component" or "fallback"
    inter_genus_ani_range: tuple[float, float] | None  # Observed min-max from metadata
    ani_values_used: int  # Number of pairwise ANI values analyzed


def detect_genus_boundary(
    ani_matrix: ANIMatrix,
    min_genomes: int = 8,
    genus_map: dict[str, str] | None = None,
) -> AdaptiveGenusThreshold:
    """
    Detect natural genus boundary from ANI matrix distribution.

    Fits a 3-component Gaussian Mixture Model to upper-triangle ANI
    values. The three components typically correspond to:
    - Between-genus comparisons (lowest ANI)
    - Between-species-within-genus comparisons (intermediate ANI)
    - Within-species comparisons (highest ANI)

    The genus boundary is estimated as the weighted midpoint between the
    lowest and middle component means, weighted by their standard deviations.

    When a genus_map is provided (genome -> genus name), empirical
    inter-genus ANI values are computed for validation. A warning is
    logged if the GMM boundary diverges from the empirical range by
    more than 5 percentage points.

    Falls back to the default threshold (80% ANI = 20% novelty) if:
    - Too few genomes (< min_genomes, default 8)
    - Too few pairwise ANI values (< 10)
    - GMM does not converge
    - Poor separation between components
    - scikit-learn is not available

    Args:
        ani_matrix: ANI matrix with pairwise genome comparisons.
        min_genomes: Minimum number of genomes required for 3-component GMM.
        genus_map: Optional mapping of genome accession to genus name.
            Used to compute empirical inter-genus ANI range for validation.

    Returns:
        AdaptiveGenusThreshold with detected or fallback boundary.
    """
    from metadarkmatter.core.constants import NOVELTY_NOVEL_GENUS_MIN

    default_boundary = 100.0 - NOVELTY_NOVEL_GENUS_MIN  # 80.0%
    genomes = sorted(ani_matrix.genomes)

    if len(genomes) < min_genomes:
        logger.info(
            "Too few genomes (%d < %d) for adaptive genus threshold detection. "
            "Using default boundary %.1f%%.",
            len(genomes), min_genomes, default_boundary,
        )
        return AdaptiveGenusThreshold(
            genus_boundary=default_boundary,
            novelty_genus_min=NOVELTY_NOVEL_GENUS_MIN,
            confidence=0.0,
            method="fallback",
            inter_genus_ani_range=None,
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

    if len(ani_values) < 10:
        logger.info(
            "Insufficient pairwise ANI values (%d < 10). Using default genus boundary.",
            len(ani_values),
        )
        return AdaptiveGenusThreshold(
            genus_boundary=default_boundary,
            novelty_genus_min=NOVELTY_NOVEL_GENUS_MIN,
            confidence=0.0,
            method="fallback",
            inter_genus_ani_range=None,
            ani_values_used=len(ani_values),
        )

    # Compute empirical inter-genus ANI range if genus_map is provided
    inter_genus_range: tuple[float, float] | None = None
    if genus_map is not None:
        inter_genus_anis = []
        for i, g1 in enumerate(genomes):
            genus1 = genus_map.get(g1)
            if genus1 is None:
                continue
            for j in range(i + 1, len(genomes)):
                g2 = genomes[j]
                genus2 = genus_map.get(g2)
                if genus2 is None:
                    continue
                if genus1 != genus2:
                    ani = ani_matrix.get_ani(g1, g2)
                    if ani > 0:
                        inter_genus_anis.append(ani)
        if inter_genus_anis:
            inter_genus_range = (min(inter_genus_anis), max(inter_genus_anis))

    ani_array = np.array(ani_values).reshape(-1, 1)

    try:
        from sklearn.mixture import GaussianMixture

        gmm = GaussianMixture(
            n_components=3,
            covariance_type="full",
            max_iter=200,
            n_init=5,
            random_state=42,
        )
        gmm.fit(ani_array)

        means = gmm.means_.flatten()
        stds = np.sqrt(gmm.covariances_.flatten())

        # Order components by mean: low (between-genus), mid (between-species), high (within-species)
        order = np.argsort(means)
        low_mean = means[order[0]]
        mid_mean = means[order[1]]
        low_std = stds[order[0]]
        mid_std = stds[order[1]]

        # Separation quality: gap between lowest and middle components
        separation = abs(mid_mean - low_mean) / (low_std + mid_std + 1e-10)

        if separation < 0.8:
            logger.info(
                "GMM components not well separated for genus boundary "
                "(separation=%.2f < 0.8). Using default boundary %.1f%%.",
                separation, default_boundary,
            )
            return AdaptiveGenusThreshold(
                genus_boundary=default_boundary,
                novelty_genus_min=NOVELTY_NOVEL_GENUS_MIN,
                confidence=separation / 3.0,
                method="fallback",
                inter_genus_ani_range=inter_genus_range,
                ani_values_used=len(ani_values),
            )

        # Genus boundary: weighted midpoint between low and mid components
        w_low = 1.0 / (low_std + 1e-10)
        w_mid = 1.0 / (mid_std + 1e-10)
        boundary = (low_mean * w_low + mid_mean * w_mid) / (w_low + w_mid)

        # Clamp to reasonable range (70-90% ANI)
        boundary = max(70.0, min(90.0, boundary))
        novelty_min = 100.0 - boundary

        confidence = min(1.0, separation / 3.0)

        # Warn if GMM boundary diverges from empirical inter-genus range
        if inter_genus_range is not None:
            empirical_mid = (inter_genus_range[0] + inter_genus_range[1]) / 2.0
            if abs(boundary - empirical_mid) > 5.0:
                logger.warning(
                    "GMM genus boundary (%.1f%%) diverges from empirical "
                    "inter-genus ANI midpoint (%.1f%%) by >5%%.",
                    boundary, empirical_mid,
                )

        logger.info(
            "Detected genus boundary at %.1f%% ANI (novelty_genus_min=%.1f%%, "
            "confidence=%.2f, method=gmm_3component, n=%d).",
            boundary, novelty_min, confidence, len(ani_values),
        )

        return AdaptiveGenusThreshold(
            genus_boundary=boundary,
            novelty_genus_min=novelty_min,
            confidence=confidence,
            method="gmm_3component",
            inter_genus_ani_range=inter_genus_range,
            ani_values_used=len(ani_values),
        )

    except ImportError:
        logger.warning(
            "scikit-learn not available for genus boundary detection. "
            "Install with: pip install metadarkmatter[adaptive]"
        )
        return AdaptiveGenusThreshold(
            genus_boundary=default_boundary,
            novelty_genus_min=NOVELTY_NOVEL_GENUS_MIN,
            confidence=0.0,
            method="fallback",
            inter_genus_ani_range=inter_genus_range,
            ani_values_used=len(ani_values),
        )
    except Exception as e:
        logger.warning(
            "3-component GMM fitting failed: %s. Using default genus boundary.", e,
        )
        return AdaptiveGenusThreshold(
            genus_boundary=default_boundary,
            novelty_genus_min=NOVELTY_NOVEL_GENUS_MIN,
            confidence=0.0,
            method="fallback",
            inter_genus_ani_range=inter_genus_range,
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
