"""
Bayesian confidence framework for classification.

Replaces hard threshold categories with posterior probabilities
P(category | novelty, uncertainty), giving users continuous confidence
rather than binary labels. Near threshold boundaries, posteriors will be
split (e.g., 55% Novel Species, 40% Known Species), making uncertainty
explicit.
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


@dataclass
class PosteriorResult:
    """Posterior probabilities for a single read."""

    p_known_species: float
    p_novel_species: float
    p_novel_genus: float
    p_ambiguous: float
    map_category: str  # Maximum A Posteriori classification


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
        """Define category centers and spreads for likelihood computation."""
        cfg = self.config

        # Category centers: (novelty_mean, uncertainty_mean)
        # Spread: (novelty_sigma, uncertainty_sigma)
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
            PosteriorResult with posterior probabilities and MAP category.
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

        return PosteriorResult(
            p_known_species=posteriors["Known Species"],
            p_novel_species=posteriors["Novel Species"],
            p_novel_genus=posteriors["Novel Genus"],
            p_ambiguous=posteriors["Ambiguous"],
            map_category=map_cat,
        )

    def classify_dataframe(
        self,
        df: pl.DataFrame,
    ) -> pl.DataFrame:
        """
        Add posterior probability columns to a classification DataFrame.

        Computes posteriors for each read and adds columns:
        - p_known_species, p_novel_species, p_novel_genus, p_ambiguous
        - bayesian_category (MAP classification)

        Args:
            df: DataFrame with novelty_index and placement_uncertainty columns.

        Returns:
            DataFrame with posterior columns added.
        """
        if df.is_empty():
            return df.with_columns([
                pl.lit(None).cast(pl.Float64).alias("p_known_species"),
                pl.lit(None).cast(pl.Float64).alias("p_novel_species"),
                pl.lit(None).cast(pl.Float64).alias("p_novel_genus"),
                pl.lit(None).cast(pl.Float64).alias("p_ambiguous"),
                pl.lit(None).cast(pl.Utf8).alias("bayesian_category"),
            ])

        # Extract arrays for vectorized computation
        novelty = df["novelty_index"].to_numpy()
        uncertainty = df["placement_uncertainty"].to_numpy()
        num_hits = df["num_ambiguous_hits"].to_numpy() if "num_ambiguous_hits" in df.columns else np.ones(len(novelty))
        identity_gap = df["identity_gap"].to_numpy() if "identity_gap" in df.columns else None

        n_reads = len(novelty)
        p_known = np.zeros(n_reads)
        p_novel_s = np.zeros(n_reads)
        p_novel_g = np.zeros(n_reads)
        p_ambig = np.zeros(n_reads)
        map_cats = []

        for i in range(n_reads):
            ig = float(identity_gap[i]) if identity_gap is not None and not np.isnan(identity_gap[i]) else None
            post = self.compute_posteriors(
                novelty_index=float(novelty[i]),
                placement_uncertainty=float(uncertainty[i]),
                num_hits=int(num_hits[i]),
                identity_gap=ig,
            )
            p_known[i] = post.p_known_species
            p_novel_s[i] = post.p_novel_species
            p_novel_g[i] = post.p_novel_genus
            p_ambig[i] = post.p_ambiguous
            map_cats.append(post.map_category)

        return df.with_columns([
            pl.Series("p_known_species", p_known).round(4),
            pl.Series("p_novel_species", p_novel_s).round(4),
            pl.Series("p_novel_genus", p_novel_g).round(4),
            pl.Series("p_ambiguous", p_ambig).round(4),
            pl.Series("bayesian_category", map_cats),
        ])
