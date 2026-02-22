"""
Bayesian classification engine for taxonomic assignment.

Computes posterior probabilities over six core taxonomic categories using
independent 2D Gaussian likelihoods in (novelty_index, placement_uncertainty)
space. A second-stage discrete refinement applies phylogenetic context
(ambiguity scope, genus-level hit patterns) to distinguish biologically
specific sub-categories (Conserved Region, Ambiguous Within Genus).

The posterior entropy provides a single scalar confidence measure:
low entropy (near 0) indicates a confident assignment, while high entropy
(near log2(6) ~ 2.585) indicates the read lies in a region where multiple
categories are similarly plausible.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import polars as pl

if TYPE_CHECKING:
    from metadarkmatter.models.config import ScoringConfig

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Category names in fixed order for vectorized computation
# ---------------------------------------------------------------------------
# Stage 1: six Gaussian categories
_CATEGORY_NAMES_6 = [
    "Known Species",
    "Novel Species",
    "Novel Genus",
    "Species Boundary",
    "Ambiguous",
    "Unclassified",
]
_N_CATEGORIES_6 = len(_CATEGORY_NAMES_6)

# Legacy 4-category names (for backward-compatible column names)
_CATEGORY_NAMES_4 = ["Known Species", "Novel Species", "Novel Genus", "Ambiguous"]
_N_CATEGORIES_4 = len(_CATEGORY_NAMES_4)

# Maximum entropy for 6 categories
_MAX_ENTROPY_6 = math.log2(_N_CATEGORIES_6)  # ~2.585


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class CategoryParams:
    """Parameters for a single Gaussian category in (N, U) space."""

    name: str
    novelty_mean: float
    uncertainty_mean: float
    novelty_sigma: float
    uncertainty_sigma: float


@dataclass
class PosteriorResult:
    """Posterior probabilities for a single read (6-category model)."""

    p_known_species: float
    p_novel_species: float
    p_novel_genus: float
    p_species_boundary: float
    p_ambiguous: float
    p_unclassified: float
    map_category: str  # Maximum A Posteriori classification (Stage 1)
    entropy: float  # Shannon entropy (0 = certain, ~2.585 = uniform over 6)

    # Backward-compatible 4-category aliases
    @property
    def bayesian_category(self) -> str:
        return self.map_category


# ---------------------------------------------------------------------------
# Parameter factory
# ---------------------------------------------------------------------------
def build_category_params(config: ScoringConfig) -> list[CategoryParams]:
    """
    Derive Gaussian parameters from ScoringConfig thresholds.

    Means are placed at the midpoints of threshold ranges; sigmas are
    range/3 + epsilon to provide adequate overlap near boundaries.

    Args:
        config: Scoring configuration with threshold values.

    Returns:
        List of CategoryParams in _CATEGORY_NAMES_6 order.
    """
    eff = config.get_effective_thresholds()

    n_known_max = eff["novelty_known_max"]
    n_sp_min = eff["novelty_novel_species_min"]
    n_sp_max = eff["novelty_novel_species_max"]
    n_gen_min = eff["novelty_novel_genus_min"]
    n_gen_max = eff["novelty_novel_genus_max"]
    u_known_max = eff["uncertainty_known_max"]
    u_sp_max = eff["uncertainty_novel_species_max"]
    u_gen_max = eff["uncertainty_novel_genus_max"]
    u_conserved = eff["uncertainty_conserved_min"]

    return [
        # Known Species: low divergence, confident placement
        CategoryParams(
            name="Known Species",
            novelty_mean=n_known_max / 2.0,
            uncertainty_mean=u_known_max / 2.0,
            novelty_sigma=n_known_max / 2.0 + 0.5,
            uncertainty_sigma=u_known_max / 2.0 + 0.5,
        ),
        # Novel Species: moderate divergence, confident placement
        CategoryParams(
            name="Novel Species",
            novelty_mean=(n_sp_min + n_sp_max) / 2.0,
            uncertainty_mean=u_sp_max / 2.0,
            novelty_sigma=(n_sp_max - n_sp_min) / 3.0 + 0.5,
            uncertainty_sigma=u_sp_max / 2.0 + 0.5,
        ),
        # Novel Genus: high divergence, confident placement
        CategoryParams(
            name="Novel Genus",
            novelty_mean=(n_gen_min + n_gen_max) / 2.0,
            uncertainty_mean=u_gen_max / 2.0,
            novelty_sigma=(n_gen_max - n_gen_min) / 3.0 + 0.5,
            uncertainty_sigma=u_gen_max / 2.0 + 0.5,
        ),
        # Species Boundary: any divergence, moderate uncertainty
        CategoryParams(
            name="Species Boundary",
            novelty_mean=(n_sp_min + n_sp_max) / 2.0,
            uncertainty_mean=(u_gen_max + u_conserved) / 2.0,
            novelty_sigma=(n_sp_max - n_sp_min) / 3.0 + 0.5,
            uncertainty_sigma=(u_conserved - u_gen_max) / 3.0 + 0.5,
        ),
        # Ambiguous: broad catch-all, high uncertainty
        CategoryParams(
            name="Ambiguous",
            novelty_mean=15.0,
            uncertainty_mean=u_conserved + 2.0,
            novelty_sigma=10.0,
            uncertainty_sigma=5.0,
        ),
        # Unclassified: very high divergence tail
        CategoryParams(
            name="Unclassified",
            novelty_mean=n_gen_max + 5.0,
            uncertainty_mean=5.0,
            novelty_sigma=5.0,
            uncertainty_sigma=5.0,
        ),
    ]


# ---------------------------------------------------------------------------
# Entropy helpers
# ---------------------------------------------------------------------------
def _shannon_entropy(posteriors: np.ndarray) -> np.ndarray:
    """
    Compute Shannon entropy of posterior distributions.

    H = -sum(p * log2(p)) for p > 0.

    Args:
        posteriors: Array of shape (n_reads, n_categories) with probabilities.

    Returns:
        Array of shape (n_reads,) with entropy values.
    """
    safe_p = np.where(posteriors > 0, posteriors, 1.0)
    return -np.sum(posteriors * np.log2(safe_p), axis=1)


def entropy_to_confidence(entropy: np.ndarray | float) -> np.ndarray | float:
    """
    Convert posterior entropy to a 0-100 confidence score.

    confidence = (1 - entropy / log2(6)) * 100

    Args:
        entropy: Shannon entropy value(s).

    Returns:
        Confidence score(s) on 0-100 scale.
    """
    return (1.0 - np.asarray(entropy) / _MAX_ENTROPY_6) * 100.0


# ---------------------------------------------------------------------------
# Stage 2: discrete refinement
# ---------------------------------------------------------------------------
def apply_stage2_refinement(
    map_categories: np.ndarray | list[str],
    ambiguity_scope: np.ndarray | list[str] | None = None,
    genus_uncertainty: np.ndarray | None = None,
    num_secondary: np.ndarray | None = None,
    genus_uncertainty_threshold: float = 10.0,
) -> np.ndarray:
    """
    Refine Stage 1 MAP categories using phylogenetic context.

    Rules:
        - Ambiguous + ambiguity_scope == "across_genera" -> Conserved Region
        - Novel Genus + genus_uncertainty >= threshold + num_secondary > 0
          -> Ambiguous Within Genus

    Args:
        map_categories: Array of Stage 1 MAP category strings.
        ambiguity_scope: Per-read ambiguity scope strings (optional).
        genus_uncertainty: Per-read genus uncertainty values (optional).
        num_secondary: Per-read count of secondary genomes (optional).
        genus_uncertainty_threshold: Threshold for Ambiguous Within Genus.

    Returns:
        Array of refined category strings (9 possible values).
    """
    cats = np.asarray(map_categories, dtype=object)
    refined = cats.copy()

    # Ambiguous + across_genera -> Conserved Region
    if ambiguity_scope is not None:
        scope = np.asarray(ambiguity_scope, dtype=object)
        mask = (cats == "Ambiguous") & (scope == "across_genera")
        refined[mask] = "Conserved Region"

    # Novel Genus + high genus uncertainty + secondary genomes -> Ambiguous Within Genus
    if genus_uncertainty is not None and num_secondary is not None:
        gu = np.asarray(genus_uncertainty, dtype=np.float64)
        ns = np.asarray(num_secondary, dtype=np.int64)
        mask = (cats == "Novel Genus") & (gu >= genus_uncertainty_threshold) & (ns > 0)
        refined[mask] = "Ambiguous Within Genus"

    return refined


# ---------------------------------------------------------------------------
# Main classifier
# ---------------------------------------------------------------------------
class BayesianClassifier:
    """
    Bayesian classification using 2D Gaussian likelihood functions.

    Each of six core categories is modeled as a 2D Gaussian centered on
    its expected (novelty, uncertainty) values. Posteriors are obtained
    via Bayes' rule with configurable priors and optional prior modulation
    based on identity gap and single-hit patterns.

    Stage 2 discrete refinement converts MAP categories to the full
    9-category taxonomy using phylogenetic context columns.
    """

    def __init__(self, config: ScoringConfig) -> None:
        self.config = config
        self._category_params = build_category_params(config)
        self._setup_vectorized_params()

    def _setup_vectorized_params(self) -> None:
        """Build numpy arrays for vectorized likelihood computation."""
        params = self._category_params
        self._novelty_means = np.array([p.novelty_mean for p in params])
        self._uncertainty_means = np.array([p.uncertainty_mean for p in params])
        self._novelty_sigmas = np.array([p.novelty_sigma for p in params])
        self._uncertainty_sigmas = np.array([p.uncertainty_sigma for p in params])

        # Prior array from config
        bay = self.config.bayesian
        self._priors = np.array([
            bay.priors.get("known_species", 1.0 / 6),
            bay.priors.get("novel_species", 1.0 / 6),
            bay.priors.get("novel_genus", 1.0 / 6),
            bay.priors.get("species_boundary", 1.0 / 6),
            bay.priors.get("ambiguous", 1.0 / 6),
            bay.priors.get("unclassified", 1.0 / 6),
        ])

        # Also keep legacy 4-category parameter dict for compute_posteriors()
        self._categories = {
            p.name: {
                "novelty_mean": p.novelty_mean,
                "uncertainty_mean": p.uncertainty_mean,
                "novelty_sigma": p.novelty_sigma,
                "uncertainty_sigma": p.uncertainty_sigma,
            }
            for p in params
        }

    # -------------------------------------------------------------------
    # Scalar API (single read)
    # -------------------------------------------------------------------
    def compute_posteriors(
        self,
        novelty_index: float,
        placement_uncertainty: float,
        num_hits: int = 1,
        identity_gap: float | None = None,
    ) -> PosteriorResult:
        """
        Compute posterior probabilities for a single read.

        Uses 6-category 2D Gaussian likelihoods with configurable priors
        and prior modulation for ambiguity signals.

        Args:
            novelty_index: Divergence from reference (100 - pident).
            placement_uncertainty: ANI-based uncertainty.
            num_hits: Number of ambiguous alignment hits.
            identity_gap: Gap between best and second-best hit identity.

        Returns:
            PosteriorResult with 6 posteriors, MAP category, and entropy.
        """
        bay = self.config.bayesian
        likelihoods = np.zeros(_N_CATEGORIES_6)

        for i, params in enumerate(self._category_params):
            z_n = (novelty_index - params.novelty_mean) / params.novelty_sigma
            z_u = (placement_uncertainty - params.uncertainty_mean) / params.uncertainty_sigma
            likelihoods[i] = np.exp(-0.5 * (z_n ** 2 + z_u ** 2))

        # Prior modulation (multiply into effective priors)
        priors = self._priors.copy()
        amb_idx = 4  # Ambiguous index in _CATEGORY_NAMES_6

        if identity_gap is not None and identity_gap < self.config.identity_gap_ambiguous_max and num_hits > 1:
            priors[amb_idx] *= bay.identity_gap_boost

        if num_hits == 1 and novelty_index >= self.config.novelty_novel_species_min:
            priors[amb_idx] *= bay.single_hit_boost

        # Bayes' rule: posterior = likelihood * prior / evidence
        joint = likelihoods * priors
        total = joint.sum()
        if total == 0:
            total = 1e-10

        posteriors = joint / total

        # MAP classification
        map_idx = int(np.argmax(posteriors))
        map_cat = _CATEGORY_NAMES_6[map_idx]

        # Shannon entropy
        safe_p = np.where(posteriors > 0, posteriors, 1.0)
        entropy_val = float(-np.sum(posteriors * np.log2(safe_p)))

        return PosteriorResult(
            p_known_species=float(posteriors[0]),
            p_novel_species=float(posteriors[1]),
            p_novel_genus=float(posteriors[2]),
            p_species_boundary=float(posteriors[3]),
            p_ambiguous=float(posteriors[4]),
            p_unclassified=float(posteriors[5]),
            map_category=map_cat,
            entropy=entropy_val,
        )

    # -------------------------------------------------------------------
    # Vectorized DataFrame API (batch of reads)
    # -------------------------------------------------------------------
    def classify_dataframe(
        self,
        df: pl.DataFrame,
    ) -> pl.DataFrame:
        """
        Add posterior probability columns to a classification DataFrame.

        Computes 6-category posteriors and adds columns:
        - p_known_species, p_novel_species, p_novel_genus,
          p_species_boundary, p_ambiguous, p_unclassified
        - bayesian_category (Stage 1 MAP classification)
        - posterior_entropy (Shannon entropy; 0 = confident, ~2.585 = uniform)

        This method is backward-compatible: when called on a DataFrame that
        already has a taxonomic_call column (e.g., from the legacy overlay path),
        it adds posterior columns without modifying the existing classification.

        Args:
            df: DataFrame with novelty_index and placement_uncertainty columns.

        Returns:
            DataFrame with posterior and entropy columns added.
        """
        empty_cols = [
            pl.lit(None).cast(pl.Float64).alias("p_known_species"),
            pl.lit(None).cast(pl.Float64).alias("p_novel_species"),
            pl.lit(None).cast(pl.Float64).alias("p_novel_genus"),
            pl.lit(None).cast(pl.Float64).alias("p_species_boundary"),
            pl.lit(None).cast(pl.Float64).alias("p_ambiguous"),
            pl.lit(None).cast(pl.Float64).alias("p_unclassified"),
            pl.lit(None).cast(pl.Utf8).alias("bayesian_category"),
            pl.lit(None).cast(pl.Float64).alias("posterior_entropy"),
        ]

        if df.is_empty():
            return df.with_columns(empty_cols)

        posteriors, map_cats, entropy = self._compute_vectorized_posteriors(df)

        return df.with_columns([
            pl.Series("p_known_species", np.round(posteriors[:, 0], 4)),
            pl.Series("p_novel_species", np.round(posteriors[:, 1], 4)),
            pl.Series("p_novel_genus", np.round(posteriors[:, 2], 4)),
            pl.Series("p_species_boundary", np.round(posteriors[:, 3], 4)),
            pl.Series("p_ambiguous", np.round(posteriors[:, 4], 4)),
            pl.Series("p_unclassified", np.round(posteriors[:, 5], 4)),
            pl.Series("bayesian_category", map_cats),
            pl.Series("posterior_entropy", np.round(entropy, 4)),
        ])

    def classify_primary(
        self,
        df: pl.DataFrame,
    ) -> pl.DataFrame:
        """
        Bayesian-primary classification: posteriors determine taxonomic_call.

        This is the main entry point for the new Bayesian-primary workflow.
        It computes 6-category posteriors, applies Stage 2 refinement, and
        derives confidence_score from posterior entropy.

        Output columns:
        - taxonomic_call: Stage 2 refined MAP category (9 possible values)
        - p_known_species ... p_unclassified: 6 posterior columns
        - posterior_entropy: Shannon entropy
        - confidence_score: (1 - entropy/log2(6)) * 100

        Args:
            df: DataFrame with novelty_index, placement_uncertainty, and
                optionally ambiguity_scope, genus_uncertainty, num_secondary_genomes.

        Returns:
            DataFrame with Bayesian classification columns.
        """
        if df.is_empty():
            return df.with_columns([
                pl.lit(None).cast(pl.Utf8).alias("taxonomic_call"),
                pl.lit(None).cast(pl.Float64).alias("p_known_species"),
                pl.lit(None).cast(pl.Float64).alias("p_novel_species"),
                pl.lit(None).cast(pl.Float64).alias("p_novel_genus"),
                pl.lit(None).cast(pl.Float64).alias("p_species_boundary"),
                pl.lit(None).cast(pl.Float64).alias("p_ambiguous"),
                pl.lit(None).cast(pl.Float64).alias("p_unclassified"),
                pl.lit(None).cast(pl.Float64).alias("posterior_entropy"),
                pl.lit(None).cast(pl.Float64).alias("confidence_score"),
            ])

        posteriors, map_cats_list, entropy = self._compute_vectorized_posteriors(df)

        # Stage 2 refinement
        ambiguity_scope = (
            df["ambiguity_scope"].to_numpy() if "ambiguity_scope" in df.columns else None
        )
        genus_uncertainty = (
            df["genus_uncertainty"].to_numpy().astype(np.float64)
            if "genus_uncertainty" in df.columns else None
        )
        num_secondary = (
            df["num_secondary_genomes"].to_numpy().astype(np.int64)
            if "num_secondary_genomes" in df.columns else None
        )

        refined = apply_stage2_refinement(
            map_cats_list,
            ambiguity_scope=ambiguity_scope,
            genus_uncertainty=genus_uncertainty,
            num_secondary=num_secondary,
            genus_uncertainty_threshold=self.config.genus_uncertainty_ambiguous_min,
        )

        # Confidence from entropy
        confidence = entropy_to_confidence(entropy)

        return df.with_columns([
            pl.Series("taxonomic_call", refined.tolist()),
            pl.Series("p_known_species", np.round(posteriors[:, 0], 4)),
            pl.Series("p_novel_species", np.round(posteriors[:, 1], 4)),
            pl.Series("p_novel_genus", np.round(posteriors[:, 2], 4)),
            pl.Series("p_species_boundary", np.round(posteriors[:, 3], 4)),
            pl.Series("p_ambiguous", np.round(posteriors[:, 4], 4)),
            pl.Series("p_unclassified", np.round(posteriors[:, 5], 4)),
            pl.Series("posterior_entropy", np.round(entropy, 4)),
            pl.Series("confidence_score", np.round(np.clip(confidence, 0, 100), 1)),
        ])

    # -------------------------------------------------------------------
    # Internal vectorized computation
    # -------------------------------------------------------------------
    def _compute_vectorized_posteriors(
        self,
        df: pl.DataFrame,
    ) -> tuple[np.ndarray, list[str], np.ndarray]:
        """
        Compute posteriors for all reads in a DataFrame.

        Returns:
            Tuple of (posteriors [n, 6], map_category_list, entropy [n]).
        """
        n_reads = df.height
        bay = self.config.bayesian

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
        # Shapes: novelty (n,), means (6,) -> z (n, 6) via broadcasting
        z_n = (novelty[:, np.newaxis] - self._novelty_means[np.newaxis, :]) / self._novelty_sigmas[np.newaxis, :]
        z_u = (uncertainty[:, np.newaxis] - self._uncertainty_means[np.newaxis, :]) / self._uncertainty_sigmas[np.newaxis, :]

        log_likes = -0.5 * (z_n ** 2 + z_u ** 2)
        likelihoods = np.exp(log_likes)

        # --- Prior modulation ---
        # Start with broadcast priors for each read
        priors = np.tile(self._priors, (n_reads, 1))  # shape (n, 6)
        amb_idx = 4  # Ambiguous index

        # Boost Ambiguous prior for small identity gap with multiple hits
        if identity_gap is not None:
            gap_boost_mask = gap_is_valid & (identity_gap < self.config.identity_gap_ambiguous_max) & (num_hits > 1)
            priors[gap_boost_mask, amb_idx] *= bay.identity_gap_boost

        # Boost Ambiguous prior for single-hit reads in novel range
        eff = self.config.get_effective_thresholds()
        single_novel_mask = (num_hits == 1) & (novelty >= eff["novelty_novel_species_min"])
        priors[single_novel_mask, amb_idx] *= bay.single_hit_boost

        # --- Bayes' rule ---
        joint = likelihoods * priors
        totals = joint.sum(axis=1, keepdims=True)
        totals = np.where(totals == 0, 1e-10, totals)
        posteriors = joint / totals  # shape (n, 6)

        # --- MAP classification ---
        map_indices = np.argmax(posteriors, axis=1)
        category_array = np.array(_CATEGORY_NAMES_6)
        map_cats = category_array[map_indices].tolist()

        # --- Shannon entropy ---
        entropy = _shannon_entropy(posteriors)

        return posteriors, map_cats, entropy
