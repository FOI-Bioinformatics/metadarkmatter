"""
Unit tests for Bayesian confidence framework.

Tests posterior computation, MAP classification, and DataFrame integration.
"""

from __future__ import annotations

import polars as pl
import pytest

from metadarkmatter.core.classification.bayesian import (
    BayesianClassifier,
    PosteriorResult,
)
from metadarkmatter.models.config import ScoringConfig


class TestBayesianClassifier:
    """Tests for BayesianClassifier."""

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
        """PosteriorResult should have all required fields."""
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


class TestBayesianDataFrame:
    """Tests for DataFrame integration."""

    @pytest.fixture
    def classifier(self):
        config = ScoringConfig()
        return BayesianClassifier(config)

    def test_classify_dataframe_adds_columns(self, classifier):
        """Should add posterior columns to DataFrame."""
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
        """Should handle empty DataFrame."""
        df = pl.DataFrame({
            "novelty_index": [],
            "placement_uncertainty": [],
        }).cast({
            "novelty_index": pl.Float64,
            "placement_uncertainty": pl.Float64,
        })

        result = classifier.classify_dataframe(df)

        assert "p_known_species" in result.columns
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
