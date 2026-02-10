"""
Bayesian confidence framework for classification.

Replaces hard threshold categories with posterior probabilities
P(category | novelty, uncertainty), giving users continuous confidence
rather than binary labels. Near threshold boundaries, posteriors will be
split (e.g., 55% Novel Species, 40% Known Species), making uncertainty
explicit.

The posterior entropy provides a single scalar measure of classification
confidence: low entropy (near 0) indicates a confident assignment to one
category, while high entropy (near 2.0 for 4 categories) indicates the
read lies in a region where multiple categories are similarly plausible.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import polars as pl

if TYPE_CHECKING:
    from metadarkmatter.models.config import ScoringConfig

logger = logging.getLogger(__name__)

# Category names in fixed order for vectorized computation
_CATEGORY_NAMES = ["Known Species", "Novel Species", "Novel Genus", "Ambiguous"]
_N_CATEGORIES = len(_CATEGORY_NAMES)


@dataclass
class PosteriorResult:
    """Posterior probabilities for a single read."""

    p_known_species: float
    p_novel_species: float
    p_novel_genus: float
    p_ambiguous: float
    map_category: str  # Maximum A Posteriori classification
    entropy: float  # Shannon entropy of posterior (0 = certain, 2.0 = uniform)


def _shannon_entropy(posteriors: np.ndarray) -> np.ndarray:
    """
    Compute Shannon entropy of posterior distributions.

    H = -sum(p * log2(p)) for p > 0.

    For 4 categories, entropy ranges from 0 (all mass on one category)
    to 2.0 (uniform distribution).

    Args:
        posteriors: Array of shape (n_reads, n_categories) with probabilities.

    Returns:
        Array of shape (n_reads,) with entropy values.
    """
    # Avoid log(0) by masking zeros
    safe_p = np.where(posteriors > 0, posteriors, 1.0)
    return -np.sum(posteriors * np.log2(safe_p), axis=1)


class BayesianClassifier:
    """
    Bayesian classification using Gaussian likelihood functions.

    Each category is modeled as a 2D Gaussian centered on its expected
    (novelty, uncertainty) values. The likelihood of a read belonging
    to each category is computed from its distance to each center,
    and posteriors are obtained via Bayes' rule with uniform priors.

    This provides continuous confidence scores near threshold boundaries,
    where hard thresholds would produce abrupt classification changes.
    """

    def __init__(self, config: ScoringConfig) -> None:
        self.config = config
        self._setup_category_params()

    def _setup_category_params(self) -> None:
        """
        Build vectorized parameter arrays for likelihood computation.

        Creates arrays of shape (4,) for novelty/uncertainty means and
        sigmas, in the order: Known Species, Novel Species, Novel Genus,
        Ambiguous. This enables fully vectorized broadcasting over reads.
        """
        cfg = self.config

        # Category parameters as dicts (for compute_posteriors scalar path)
        self._categories = {
            "Known Species": {
                "novelty_mean": cfg.novelty_known_max / 2.0,
                "uncertainty_mean": cfg.uncertainty_known_max / 2.0,
                "novelty_sigma": cfg.novelty_known_max / 2.0 + 0.5,
                "uncertainty_sigma": cfg.uncertainty_known_max / 2.0 + 0.5,
            },
            "Novel Species": {
                "novelty_mean": (cfg.novelty_novel_species_min + cfg.novelty_novel_species_max) / 2.0,
                "uncertainty_mean": cfg.uncertainty_novel_species_max / 2.0,
                "novelty_sigma": (cfg.novelty_novel_species_max - cfg.novelty_novel_species_min) / 3.0 + 0.5,
                "uncertainty_sigma": cfg.uncertainty_novel_species_max / 2.0 + 0.5,
            },
            "Novel Genus": {
                "novelty_mean": (cfg.novelty_novel_genus_min + cfg.novelty_novel_genus_max) / 2.0,
                "uncertainty_mean": cfg.uncertainty_novel_genus_max / 2.0,
                "novelty_sigma": (cfg.novelty_novel_genus_max - cfg.novelty_novel_genus_min) / 3.0 + 0.5,
                "uncertainty_sigma": cfg.uncertainty_novel_genus_max / 2.0 + 0.5,
            },
            "Ambiguous": {
                "novelty_mean": 15.0,
                "uncertainty_mean": cfg.uncertainty_conserved_min,
                "novelty_sigma": 10.0,
                "uncertainty_sigma": 5.0,
            },
        }

        # Vectorized parameter arrays (shape: (4,)) in _CATEGORY_NAMES order
        self._novelty_means = np.array([
            self._categories[c]["novelty_mean"] for c in _CATEGORY_NAMES
        ])
        self._uncertainty_means = np.array([
            self._categories[c]["uncertainty_mean"] for c in _CATEGORY_NAMES
        ])
        self._novelty_sigmas = np.array([
            self._categories[c]["novelty_sigma"] for c in _CATEGORY_NAMES
        ])
        self._uncertainty_sigmas = np.array([
            self._categories[c]["uncertainty_sigma"] for c in _CATEGORY_NAMES
        ])

    def compute_posteriors(
        self,
        novelty_index: float,
        placement_uncertainty: float,
        num_hits: int = 1,
        identity_gap: float | None = None,
    ) -> PosteriorResult:
        """
        Compute posterior probabilities for a single read.

        Uses Gaussian likelihoods centered on each category's expected values
        with uniform priors. Additional factors from num_hits and identity_gap
        modulate the ambiguous category probability.

        Args:
            novelty_index: Divergence from reference (100 - pident).
            placement_uncertainty: ANI-based uncertainty.
            num_hits: Number of ambiguous alignment hits.
            identity_gap: Gap between best and second-best hit identity.

        Returns:
            PosteriorResult with posterior probabilities, MAP category, and entropy.
        """
        likelihoods = {}

        for cat, params in self._categories.items():
            # 2D Gaussian log-likelihood
            z_n = (novelty_index - params["novelty_mean"]) / params["novelty_sigma"]
            z_u = (placement_uncertainty - params["uncertainty_mean"]) / params["uncertainty_sigma"]
            log_like = -0.5 * (z_n ** 2 + z_u ** 2)
            likelihoods[cat] = np.exp(log_like)

        # Boost Ambiguous likelihood for high-uncertainty or low-gap reads
        if identity_gap is not None and identity_gap < 2.0 and num_hits > 1:
            likelihoods["Ambiguous"] *= 3.0
        if placement_uncertainty > self.config.uncertainty_conserved_min:
            likelihoods["Ambiguous"] *= 2.0

        # Normalize to posteriors (uniform prior)
        total = sum(likelihoods.values())
        if total == 0:
            total = 1e-10

        posteriors = {cat: like / total for cat, like in likelihoods.items()}

        # MAP classification
        map_cat = max(posteriors, key=posteriors.get)

        # Shannon entropy
        probs = np.array([posteriors[c] for c in _CATEGORY_NAMES])
        safe_p = np.where(probs > 0, probs, 1.0)
        entropy = float(-np.sum(probs * np.log2(safe_p)))

        return PosteriorResult(
            p_known_species=posteriors["Known Species"],
            p_novel_species=posteriors["Novel Species"],
            p_novel_genus=posteriors["Novel Genus"],
            p_ambiguous=posteriors["Ambiguous"],
            map_category=map_cat,
            entropy=entropy,
        )

    def classify_dataframe(
        self,
        df: pl.DataFrame,
    ) -> pl.DataFrame:
        """
        Add posterior probability columns to a classification DataFrame.

        Uses fully vectorized numpy operations (no Python loops) for
        efficient computation over large DataFrames.

        Computes posteriors for each read and adds columns:
        - p_known_species, p_novel_species, p_novel_genus, p_ambiguous
        - bayesian_category (MAP classification)
        - posterior_entropy (Shannon entropy; 0 = confident, 2.0 = uniform)

        Args:
            df: DataFrame with novelty_index and placement_uncertainty columns.

        Returns:
            DataFrame with posterior and entropy columns added.
        """
        if df.is_empty():
            return df.with_columns([
                pl.lit(None).cast(pl.Float64).alias("p_known_species"),
                pl.lit(None).cast(pl.Float64).alias("p_novel_species"),
                pl.lit(None).cast(pl.Float64).alias("p_novel_genus"),
                pl.lit(None).cast(pl.Float64).alias("p_ambiguous"),
                pl.lit(None).cast(pl.Utf8).alias("bayesian_category"),
                pl.lit(None).cast(pl.Float64).alias("posterior_entropy"),
            ])

        n_reads = df.height

        # Extract input arrays
        novelty = df["novelty_index"].to_numpy().astype(np.float64)
        uncertainty = df["placement_uncertainty"].to_numpy().astype(np.float64)

        has_hits = "num_ambiguous_hits" in df.columns
        num_hits = (
            df["num_ambiguous_hits"].to_numpy().astype(np.float64)
            if has_hits
            else np.ones(n_reads)
        )

        has_gap = "identity_gap" in df.columns
        if has_gap:
            identity_gap_raw = df["identity_gap"].to_numpy().astype(np.float64)
            gap_is_valid = ~np.isnan(identity_gap_raw)
            identity_gap = np.where(gap_is_valid, identity_gap_raw, np.inf)
        else:
            identity_gap = None
            gap_is_valid = None

        # --- Vectorized likelihood computation ---
        # Shapes: novelty (n,), _novelty_means (4,) -> z_n (n, 4) via broadcasting
        z_n = (novelty[:, np.newaxis] - self._novelty_means[np.newaxis, :]) / self._novelty_sigmas[np.newaxis, :]
        z_u = (uncertainty[:, np.newaxis] - self._uncertainty_means[np.newaxis, :]) / self._uncertainty_sigmas[np.newaxis, :]

        # 2D Gaussian log-likelihoods -> likelihoods, shape (n, 4)
        log_likes = -0.5 * (z_n ** 2 + z_u ** 2)
        likelihoods = np.exp(log_likes)

        # Ambiguous boosting (column index 3 = Ambiguous)
        # Boost for small identity gap with multiple hits
        if identity_gap is not None:
            gap_boost = gap_is_valid & (identity_gap < 2.0) & (num_hits > 1)
            likelihoods[gap_boost, 3] *= 3.0

        # Boost for high placement uncertainty
        high_u = uncertainty > self.config.uncertainty_conserved_min
        likelihoods[high_u, 3] *= 2.0

        # --- Normalize to posteriors ---
        totals = likelihoods.sum(axis=1, keepdims=True)
        totals = np.where(totals == 0, 1e-10, totals)
        posteriors = likelihoods / totals  # shape (n, 4)

        # --- MAP classification ---
        map_indices = np.argmax(posteriors, axis=1)
        category_array = np.array(_CATEGORY_NAMES)
        map_cats = category_array[map_indices].tolist()

        # --- Shannon entropy ---
        entropy = _shannon_entropy(posteriors)

        return df.with_columns([
            pl.Series("p_known_species", np.round(posteriors[:, 0], 4)),
            pl.Series("p_novel_species", np.round(posteriors[:, 1], 4)),
            pl.Series("p_novel_genus", np.round(posteriors[:, 2], 4)),
            pl.Series("p_ambiguous", np.round(posteriors[:, 3], 4)),
            pl.Series("bayesian_category", map_cats),
            pl.Series("posterior_entropy", np.round(entropy, 4)),
        ])
