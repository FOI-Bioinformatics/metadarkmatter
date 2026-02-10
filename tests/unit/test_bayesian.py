"""
Unit tests for Bayesian confidence framework.

Tests posterior computation, vectorized DataFrame classification,
entropy calculation, and edge case handling.
"""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from metadarkmatter.core.classification.bayesian import (
    BayesianClassifier,
    PosteriorResult,
    _shannon_entropy,
)
from metadarkmatter.models.config import ScoringConfig


class TestBayesianClassifier:
    """Tests for BayesianClassifier scalar computation."""

    @pytest.fixture
    def classifier(self):
        config = ScoringConfig()
        return BayesianClassifier(config)

    def test_posteriors_sum_to_one(self, classifier):
        """Posterior probabilities should sum to 1.0."""
        result = classifier.compute_posteriors(
            novelty_index=5.0,
            placement_uncertainty=1.0,
        )
        total = (
            result.p_known_species
            + result.p_novel_species
            + result.p_novel_genus
            + result.p_ambiguous
        )
        assert total == pytest.approx(1.0, abs=1e-6)

    def test_known_species_highest_for_low_novelty(self, classifier):
        """Low novelty should give highest posterior for Known Species."""
        result = classifier.compute_posteriors(
            novelty_index=1.0,
            placement_uncertainty=0.5,
        )
        assert result.p_known_species > result.p_novel_species
        assert result.p_known_species > result.p_novel_genus
        assert result.map_category == "Known Species"

    def test_novel_species_highest_for_moderate_novelty(self, classifier):
        """Moderate novelty should give highest posterior for Novel Species."""
        result = classifier.compute_posteriors(
            novelty_index=12.0,
            placement_uncertainty=0.5,
            num_hits=3,
        )
        assert result.p_novel_species > result.p_known_species
        assert result.map_category == "Novel Species"

    def test_novel_genus_highest_for_high_novelty(self, classifier):
        """High novelty should give highest posterior for Novel Genus."""
        result = classifier.compute_posteriors(
            novelty_index=22.0,
            placement_uncertainty=0.5,
            num_hits=3,
        )
        assert result.p_novel_genus > result.p_known_species
        assert result.map_category == "Novel Genus"

    def test_ambiguous_boosted_by_small_gap(self, classifier):
        """Small identity gap should boost Ambiguous posterior."""
        result_gap = classifier.compute_posteriors(
            novelty_index=5.0,
            placement_uncertainty=1.0,
            num_hits=3,
            identity_gap=0.5,
        )
        result_no_gap = classifier.compute_posteriors(
            novelty_index=5.0,
            placement_uncertainty=1.0,
            num_hits=3,
            identity_gap=10.0,
        )
        assert result_gap.p_ambiguous > result_no_gap.p_ambiguous

    def test_ambiguous_boosted_by_high_uncertainty(self, classifier):
        """High placement uncertainty should boost Ambiguous posterior."""
        result_high_u = classifier.compute_posteriors(
            novelty_index=10.0,
            placement_uncertainty=8.0,
        )
        result_low_u = classifier.compute_posteriors(
            novelty_index=10.0,
            placement_uncertainty=0.5,
        )
        assert result_high_u.p_ambiguous > result_low_u.p_ambiguous

    def test_boundary_region_split_posteriors(self, classifier):
        """Near boundary, posteriors should be split (not dominated by one)."""
        # At the species boundary (novelty = 4.0 for default config)
        result = classifier.compute_posteriors(
            novelty_index=4.0,
            placement_uncertainty=0.5,
            num_hits=3,
        )
        # Neither category should dominate completely
        max_posterior = max(
            result.p_known_species,
            result.p_novel_species,
            result.p_novel_genus,
            result.p_ambiguous,
        )
        assert max_posterior < 0.95  # No single category should be >95%

    def test_posterior_result_dataclass(self, classifier):
        """PosteriorResult should have all required fields including entropy."""
        result = classifier.compute_posteriors(
            novelty_index=5.0,
            placement_uncertainty=1.0,
        )
        assert isinstance(result, PosteriorResult)
        assert hasattr(result, "p_known_species")
        assert hasattr(result, "p_novel_species")
        assert hasattr(result, "p_novel_genus")
        assert hasattr(result, "p_ambiguous")
        assert hasattr(result, "map_category")
        assert hasattr(result, "entropy")

    def test_entropy_low_for_confident_classification(self, classifier):
        """Entropy should be lower when one category dominates."""
        result = classifier.compute_posteriors(
            novelty_index=1.0,
            placement_uncertainty=0.3,
        )
        # Very low novelty, low uncertainty -> Known Species is dominant
        # With 4 categories (max entropy 2.0), entropy below 1.5 indicates
        # meaningful concentration on one category
        assert result.entropy < 1.5
        assert result.p_known_species > 0.5

    def test_entropy_higher_at_boundary(self, classifier):
        """Entropy should be higher at classification boundaries."""
        boundary = classifier.compute_posteriors(
            novelty_index=4.0,
            placement_uncertainty=0.5,
            num_hits=3,
        )
        confident = classifier.compute_posteriors(
            novelty_index=1.0,
            placement_uncertainty=0.3,
        )
        assert boundary.entropy > confident.entropy

    def test_entropy_range(self, classifier):
        """Entropy should be between 0 and 2.0 (log2(4))."""
        for ni in [0.5, 2.0, 5.0, 12.0, 22.0, 30.0]:
            for pu in [0.1, 0.5, 1.0, 3.0, 8.0]:
                result = classifier.compute_posteriors(
                    novelty_index=ni,
                    placement_uncertainty=pu,
                )
                assert 0.0 <= result.entropy <= 2.0 + 1e-6

    def test_extreme_high_novelty(self, classifier):
        """Extreme novelty values should not cause errors."""
        result = classifier.compute_posteriors(
            novelty_index=50.0,
            placement_uncertainty=0.5,
        )
        total = (
            result.p_known_species + result.p_novel_species
            + result.p_novel_genus + result.p_ambiguous
        )
        assert total == pytest.approx(1.0, abs=1e-6)

    def test_zero_novelty_zero_uncertainty(self, classifier):
        """Zero novelty and uncertainty should favor Known Species."""
        result = classifier.compute_posteriors(
            novelty_index=0.0,
            placement_uncertainty=0.0,
        )
        assert result.map_category == "Known Species"
        total = (
            result.p_known_species + result.p_novel_species
            + result.p_novel_genus + result.p_ambiguous
        )
        assert total == pytest.approx(1.0, abs=1e-6)


class TestBayesianAdaptiveIntegration:
    """Tests for Bayesian classification with adaptive thresholds."""

    def test_adaptive_thresholds_shift_category_centers(self):
        """Adaptive thresholds should shift the Bayesian category centers."""
        default_clf = BayesianClassifier(ScoringConfig())
        adaptive_clf = BayesianClassifier(ScoringConfig(
            novelty_known_max=6.0,
            novelty_novel_species_min=6.0,
        ))

        # Known Species center should be at novelty_known_max / 2
        assert default_clf._categories["Known Species"]["novelty_mean"] == 2.0  # 4.0 / 2
        assert adaptive_clf._categories["Known Species"]["novelty_mean"] == 3.0  # 6.0 / 2

    def test_wider_known_boundary_reclassifies_reads(self):
        """Widening Known Species boundary should shift reads from Novel to Known."""
        strict_clf = BayesianClassifier(ScoringConfig(
            novelty_known_max=3.0,
            novelty_novel_species_min=3.0,
        ))
        relaxed_clf = BayesianClassifier(ScoringConfig(
            novelty_known_max=6.0,
            novelty_novel_species_min=6.0,
        ))

        # A read at novelty=4.5 should be more likely Known with relaxed boundary
        strict_result = strict_clf.compute_posteriors(
            novelty_index=4.5, placement_uncertainty=0.5,
        )
        relaxed_result = relaxed_clf.compute_posteriors(
            novelty_index=4.5, placement_uncertainty=0.5,
        )
        assert relaxed_result.p_known_species > strict_result.p_known_species


class TestShannonEntropy:
    """Tests for the _shannon_entropy function."""

    def test_uniform_distribution(self):
        """Uniform distribution over 4 categories should give entropy = 2.0."""
        posteriors = np.array([[0.25, 0.25, 0.25, 0.25]])
        entropy = _shannon_entropy(posteriors)
        assert entropy[0] == pytest.approx(2.0, abs=1e-6)

    def test_certain_distribution(self):
        """All mass on one category should give entropy = 0."""
        posteriors = np.array([[1.0, 0.0, 0.0, 0.0]])
        entropy = _shannon_entropy(posteriors)
        assert entropy[0] == pytest.approx(0.0, abs=1e-6)

    def test_binary_split(self):
        """50/50 split over two categories should give entropy = 1.0."""
        posteriors = np.array([[0.5, 0.5, 0.0, 0.0]])
        entropy = _shannon_entropy(posteriors)
        assert entropy[0] == pytest.approx(1.0, abs=1e-6)

    def test_vectorized_over_multiple_reads(self):
        """Should compute entropy for multiple reads at once."""
        posteriors = np.array([
            [1.0, 0.0, 0.0, 0.0],   # certain
            [0.25, 0.25, 0.25, 0.25],  # uniform
            [0.5, 0.5, 0.0, 0.0],   # binary
        ])
        entropy = _shannon_entropy(posteriors)
        assert entropy[0] == pytest.approx(0.0, abs=1e-6)
        assert entropy[1] == pytest.approx(2.0, abs=1e-6)
        assert entropy[2] == pytest.approx(1.0, abs=1e-6)


class TestBayesianDataFrame:
    """Tests for vectorized DataFrame integration."""

    @pytest.fixture
    def classifier(self):
        config = ScoringConfig()
        return BayesianClassifier(config)

    def test_classify_dataframe_adds_columns(self, classifier):
        """Should add posterior and entropy columns to DataFrame."""
        df = pl.DataFrame({
            "novelty_index": [1.0, 10.0, 22.0],
            "placement_uncertainty": [0.5, 0.5, 0.5],
            "num_ambiguous_hits": [3, 3, 3],
            "identity_gap": [5.0, 5.0, 5.0],
        })

        result = classifier.classify_dataframe(df)

        assert "p_known_species" in result.columns
        assert "p_novel_species" in result.columns
        assert "p_novel_genus" in result.columns
        assert "p_ambiguous" in result.columns
        assert "bayesian_category" in result.columns
        assert "posterior_entropy" in result.columns

    def test_classify_dataframe_posteriors_sum(self, classifier):
        """Posteriors in DataFrame should sum to ~1.0 for each row."""
        df = pl.DataFrame({
            "novelty_index": [1.0, 10.0, 22.0],
            "placement_uncertainty": [0.5, 0.5, 0.5],
            "num_ambiguous_hits": [3, 3, 3],
            "identity_gap": [5.0, 5.0, 5.0],
        })

        result = classifier.classify_dataframe(df)

        sums = (
            result["p_known_species"]
            + result["p_novel_species"]
            + result["p_novel_genus"]
            + result["p_ambiguous"]
        )
        for s in sums.to_list():
            assert s == pytest.approx(1.0, abs=0.01)

    def test_classify_empty_dataframe(self, classifier):
        """Should handle empty DataFrame with entropy column."""
        df = pl.DataFrame({
            "novelty_index": [],
            "placement_uncertainty": [],
        }).cast({
            "novelty_index": pl.Float64,
            "placement_uncertainty": pl.Float64,
        })

        result = classifier.classify_dataframe(df)

        assert "p_known_species" in result.columns
        assert "posterior_entropy" in result.columns
        assert result.height == 0

    def test_classify_dataframe_preserves_existing_columns(self, classifier):
        """Should preserve all existing columns."""
        df = pl.DataFrame({
            "read_id": ["r1", "r2"],
            "novelty_index": [1.0, 10.0],
            "placement_uncertainty": [0.5, 0.5],
            "num_ambiguous_hits": [3, 3],
            "identity_gap": [5.0, 5.0],
        })

        result = classifier.classify_dataframe(df)

        assert "read_id" in result.columns
        assert result["read_id"].to_list() == ["r1", "r2"]

    def test_handles_null_identity_gap(self, classifier):
        """Should handle null identity_gap values."""
        df = pl.DataFrame({
            "novelty_index": [1.0, 10.0],
            "placement_uncertainty": [0.5, 0.5],
            "num_ambiguous_hits": [1, 3],
            "identity_gap": [None, 5.0],
        })

        result = classifier.classify_dataframe(df)

        assert result.height == 2
        assert result["p_known_species"][0] is not None
        assert result["posterior_entropy"][0] is not None

    def test_vectorized_matches_scalar(self, classifier):
        """Vectorized DataFrame results should match scalar compute_posteriors."""
        test_cases = [
            {"novelty_index": 1.0, "placement_uncertainty": 0.3, "num_ambiguous_hits": 1, "identity_gap": None},
            {"novelty_index": 12.0, "placement_uncertainty": 0.5, "num_ambiguous_hits": 3, "identity_gap": 5.0},
            {"novelty_index": 22.0, "placement_uncertainty": 0.5, "num_ambiguous_hits": 3, "identity_gap": 0.5},
            {"novelty_index": 5.0, "placement_uncertainty": 8.0, "num_ambiguous_hits": 2, "identity_gap": 1.0},
        ]

        df = pl.DataFrame({
            "novelty_index": [c["novelty_index"] for c in test_cases],
            "placement_uncertainty": [c["placement_uncertainty"] for c in test_cases],
            "num_ambiguous_hits": [c["num_ambiguous_hits"] for c in test_cases],
            "identity_gap": [c["identity_gap"] for c in test_cases],
        })

        vec_result = classifier.classify_dataframe(df)

        for i, case in enumerate(test_cases):
            scalar = classifier.compute_posteriors(
                novelty_index=case["novelty_index"],
                placement_uncertainty=case["placement_uncertainty"],
                num_hits=case["num_ambiguous_hits"],
                identity_gap=case["identity_gap"],
            )
            assert vec_result["p_known_species"][i] == pytest.approx(
                round(scalar.p_known_species, 4), abs=1e-4
            ), f"p_known_species mismatch at row {i}"
            assert vec_result["p_novel_species"][i] == pytest.approx(
                round(scalar.p_novel_species, 4), abs=1e-4
            ), f"p_novel_species mismatch at row {i}"
            assert vec_result["p_novel_genus"][i] == pytest.approx(
                round(scalar.p_novel_genus, 4), abs=1e-4
            ), f"p_novel_genus mismatch at row {i}"
            assert vec_result["p_ambiguous"][i] == pytest.approx(
                round(scalar.p_ambiguous, 4), abs=1e-4
            ), f"p_ambiguous mismatch at row {i}"
            assert vec_result["bayesian_category"][i] == scalar.map_category, (
                f"MAP category mismatch at row {i}"
            )
            assert vec_result["posterior_entropy"][i] == pytest.approx(
                round(scalar.entropy, 4), abs=1e-4
            ), f"entropy mismatch at row {i}"

    def test_dataframe_without_optional_columns(self, classifier):
        """Should work without num_ambiguous_hits and identity_gap columns."""
        df = pl.DataFrame({
            "novelty_index": [1.0, 12.0, 22.0],
            "placement_uncertainty": [0.5, 0.5, 0.5],
        })

        result = classifier.classify_dataframe(df)

        assert result.height == 3
        assert "bayesian_category" in result.columns
        assert "posterior_entropy" in result.columns
        sums = (
            result["p_known_species"]
            + result["p_novel_species"]
            + result["p_novel_genus"]
            + result["p_ambiguous"]
        )
        for s in sums.to_list():
            assert s == pytest.approx(1.0, abs=0.01)

    def test_entropy_column_values(self, classifier):
        """Entropy column should contain valid values in [0, 2.0]."""
        df = pl.DataFrame({
            "novelty_index": [0.5, 4.0, 12.0, 22.0, 50.0],
            "placement_uncertainty": [0.1, 1.0, 0.5, 0.5, 10.0],
        })

        result = classifier.classify_dataframe(df)

        for e in result["posterior_entropy"].to_list():
            assert 0.0 <= e <= 2.0 + 0.001

    def test_large_dataframe_performance(self, classifier):
        """Should handle large DataFrames without error (vectorized path)."""
        n = 10_000
        rng = np.random.default_rng(42)
        df = pl.DataFrame({
            "novelty_index": rng.uniform(0, 30, n).tolist(),
            "placement_uncertainty": rng.uniform(0, 10, n).tolist(),
            "num_ambiguous_hits": rng.integers(1, 10, n).tolist(),
            "identity_gap": rng.uniform(0, 15, n).tolist(),
        })

        result = classifier.classify_dataframe(df)

        assert result.height == n
        assert "posterior_entropy" in result.columns
        # Posteriors should still sum to 1
        sums = (
            result["p_known_species"]
            + result["p_novel_species"]
            + result["p_novel_genus"]
            + result["p_ambiguous"]
        )
        for s in sums.to_list():
            assert s == pytest.approx(1.0, abs=0.01)
