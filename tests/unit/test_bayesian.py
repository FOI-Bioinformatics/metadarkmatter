"""
Unit tests for Bayesian confidence framework.

Tests posterior computation, vectorized DataFrame classification,
entropy calculation, and edge case handling for the 6-category model.
"""

from __future__ import annotations

import math

import numpy as np
import polars as pl
import pytest

from metadarkmatter.core.classification.bayesian import (
    BayesianClassifier,
    PosteriorResult,
    _shannon_entropy,
    _MAX_ENTROPY_6,
    apply_stage2_refinement,
    entropy_to_confidence,
    build_category_params,
)
from metadarkmatter.models.config import ScoringConfig


def _sum_all_6(result: PosteriorResult) -> float:
    """Sum all 6 posterior probabilities."""
    return (
        result.p_known_species
        + result.p_novel_species
        + result.p_novel_genus
        + result.p_species_boundary
        + result.p_ambiguous
        + result.p_unclassified
    )


def _sum_all_6_df(result: pl.DataFrame) -> pl.Series:
    """Sum all 6 posterior columns in a DataFrame."""
    return (
        result["p_known_species"]
        + result["p_novel_species"]
        + result["p_novel_genus"]
        + result["p_species_boundary"]
        + result["p_ambiguous"]
        + result["p_unclassified"]
    )


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
        assert _sum_all_6(result) == pytest.approx(1.0, abs=1e-6)

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
            result.p_species_boundary,
            result.p_ambiguous,
            result.p_unclassified,
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
        assert hasattr(result, "p_species_boundary")
        assert hasattr(result, "p_ambiguous")
        assert hasattr(result, "p_unclassified")
        assert hasattr(result, "map_category")
        assert hasattr(result, "entropy")

    def test_entropy_low_for_confident_classification(self, classifier):
        """Entropy should be lower when one category dominates."""
        result = classifier.compute_posteriors(
            novelty_index=1.0,
            placement_uncertainty=0.3,
        )
        # Very low novelty, low uncertainty -> Known Species is dominant
        # With 6 categories (max entropy ~2.585), entropy below 2.0 indicates
        # meaningful concentration on one category
        assert result.entropy < 2.0
        assert result.p_known_species > 0.3

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
        """Entropy should be between 0 and log2(6) ~ 2.585."""
        for ni in [0.5, 2.0, 5.0, 12.0, 22.0, 30.0]:
            for pu in [0.1, 0.5, 1.0, 3.0, 8.0]:
                result = classifier.compute_posteriors(
                    novelty_index=ni,
                    placement_uncertainty=pu,
                )
                assert 0.0 <= result.entropy <= _MAX_ENTROPY_6 + 1e-6

    def test_extreme_high_novelty(self, classifier):
        """Extreme novelty values should not cause errors."""
        result = classifier.compute_posteriors(
            novelty_index=50.0,
            placement_uncertainty=0.5,
        )
        assert _sum_all_6(result) == pytest.approx(1.0, abs=1e-6)

    def test_zero_novelty_zero_uncertainty(self, classifier):
        """Zero novelty and uncertainty should favor Known Species."""
        result = classifier.compute_posteriors(
            novelty_index=0.0,
            placement_uncertainty=0.0,
        )
        assert result.map_category == "Known Species"
        assert _sum_all_6(result) == pytest.approx(1.0, abs=1e-6)


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

    def test_uniform_distribution_6(self):
        """Uniform distribution over 6 categories should give entropy = log2(6)."""
        p = 1.0 / 6
        posteriors = np.array([[p, p, p, p, p, p]])
        entropy = _shannon_entropy(posteriors)
        assert entropy[0] == pytest.approx(_MAX_ENTROPY_6, abs=1e-6)

    def test_uniform_distribution_4(self):
        """Uniform distribution over 4 categories should give entropy = 2.0."""
        posteriors = np.array([[0.25, 0.25, 0.25, 0.25]])
        entropy = _shannon_entropy(posteriors)
        assert entropy[0] == pytest.approx(2.0, abs=1e-6)

    def test_certain_distribution(self):
        """All mass on one category should give entropy = 0."""
        posteriors = np.array([[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]])
        entropy = _shannon_entropy(posteriors)
        assert entropy[0] == pytest.approx(0.0, abs=1e-6)

    def test_binary_split(self):
        """50/50 split over two categories should give entropy = 1.0."""
        posteriors = np.array([[0.5, 0.5, 0.0, 0.0, 0.0, 0.0]])
        entropy = _shannon_entropy(posteriors)
        assert entropy[0] == pytest.approx(1.0, abs=1e-6)

    def test_vectorized_over_multiple_reads(self):
        """Should compute entropy for multiple reads at once."""
        p = 1.0 / 6
        posteriors = np.array([
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],     # certain
            [p, p, p, p, p, p],                     # uniform over 6
            [0.5, 0.5, 0.0, 0.0, 0.0, 0.0],       # binary
        ])
        entropy = _shannon_entropy(posteriors)
        assert entropy[0] == pytest.approx(0.0, abs=1e-6)
        assert entropy[1] == pytest.approx(_MAX_ENTROPY_6, abs=1e-6)
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
        assert "p_species_boundary" in result.columns
        assert "p_ambiguous" in result.columns
        assert "p_unclassified" in result.columns
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

        sums = _sum_all_6_df(result)
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
            assert vec_result["p_species_boundary"][i] == pytest.approx(
                round(scalar.p_species_boundary, 4), abs=1e-4
            ), f"p_species_boundary mismatch at row {i}"
            assert vec_result["p_ambiguous"][i] == pytest.approx(
                round(scalar.p_ambiguous, 4), abs=1e-4
            ), f"p_ambiguous mismatch at row {i}"
            assert vec_result["p_unclassified"][i] == pytest.approx(
                round(scalar.p_unclassified, 4), abs=1e-4
            ), f"p_unclassified mismatch at row {i}"
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
        sums = _sum_all_6_df(result)
        for s in sums.to_list():
            assert s == pytest.approx(1.0, abs=0.01)

    def test_entropy_column_values(self, classifier):
        """Entropy column should contain valid values in [0, log2(6)]."""
        df = pl.DataFrame({
            "novelty_index": [0.5, 4.0, 12.0, 22.0, 50.0],
            "placement_uncertainty": [0.1, 1.0, 0.5, 0.5, 10.0],
        })

        result = classifier.classify_dataframe(df)

        for e in result["posterior_entropy"].to_list():
            assert 0.0 <= e <= _MAX_ENTROPY_6 + 0.001

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
        sums = _sum_all_6_df(result)
        for s in sums.to_list():
            assert s == pytest.approx(1.0, abs=0.01)


class TestStage2Refinement:
    """Tests for apply_stage2_refinement() discrete post-processing."""

    def test_ambiguous_across_genera_becomes_conserved_region(self):
        """Ambiguous + ambiguity_scope 'across_genera' should refine to Conserved Region."""
        categories = ["Ambiguous"]
        genus_uncertainty = np.array([0.0])
        num_secondary = np.array([0])
        ambiguity_scope = np.array(["across_genera"], dtype=object)

        refined = apply_stage2_refinement(
            categories,
            ambiguity_scope=ambiguity_scope,
            genus_uncertainty=genus_uncertainty,
            num_secondary=num_secondary,
        )

        assert refined[0] == "Conserved Region"

    def test_novel_genus_high_uncertainty_becomes_ambiguous_within_genus(self):
        """Novel Genus + genus_uncertainty >= threshold + secondary > 0 -> Ambiguous Within Genus."""
        categories = ["Novel Genus"]
        genus_uncertainty = np.array([12.0])
        num_secondary = np.array([2])

        refined = apply_stage2_refinement(
            categories,
            genus_uncertainty=genus_uncertainty,
            num_secondary=num_secondary,
            genus_uncertainty_threshold=10.0,
        )

        assert refined[0] == "Ambiguous Within Genus"

    def test_other_categories_pass_through_unchanged(self):
        """Categories not matching refinement rules should pass through."""
        for cat in ["Known Species", "Novel Species", "Species Boundary", "Unclassified"]:
            refined = apply_stage2_refinement(
                [cat],
                ambiguity_scope=np.array(["across_genera"], dtype=object),
                genus_uncertainty=np.array([15.0]),
                num_secondary=np.array([3]),
            )
            assert refined[0] == cat, f"Expected '{cat}' to pass through, got '{refined[0]}'"

    def test_batch_multiple_categories(self):
        """Should correctly refine a mixed batch of categories."""
        categories = [
            "Ambiguous",       # across_genera -> Conserved Region
            "Novel Genus",     # high genus_uncertainty + secondary -> Ambiguous Within Genus
            "Known Species",   # unchanged
            "Novel Species",   # unchanged
            "Ambiguous",       # within_genus scope -> unchanged (not across_genera)
        ]
        ambiguity_scope = np.array(
            ["across_genera", "within_genus", "within_genus", "within_genus", "within_genus"],
            dtype=object,
        )
        genus_uncertainty = np.array([0.0, 15.0, 0.0, 0.0, 0.0])
        num_secondary = np.array([0, 3, 0, 0, 0])

        refined = apply_stage2_refinement(
            categories,
            ambiguity_scope=ambiguity_scope,
            genus_uncertainty=genus_uncertainty,
            num_secondary=num_secondary,
            genus_uncertainty_threshold=10.0,
        )

        assert refined[0] == "Conserved Region"
        assert refined[1] == "Ambiguous Within Genus"
        assert refined[2] == "Known Species"
        assert refined[3] == "Novel Species"
        assert refined[4] == "Ambiguous"

    def test_without_ambiguity_scope(self):
        """Should work without ambiguity_scope (no Conserved Region refinement)."""
        categories = ["Ambiguous", "Novel Genus"]
        genus_uncertainty = np.array([0.0, 15.0])
        num_secondary = np.array([0, 2])

        refined = apply_stage2_refinement(
            categories,
            genus_uncertainty=genus_uncertainty,
            num_secondary=num_secondary,
            genus_uncertainty_threshold=10.0,
        )

        # Without ambiguity_scope, Ambiguous stays Ambiguous (no Conserved Region rule)
        assert refined[0] == "Ambiguous"
        # Novel Genus refinement still applies
        assert refined[1] == "Ambiguous Within Genus"


class TestEntropyToConfidence:
    """Tests for entropy_to_confidence() scalar conversion."""

    def test_zero_entropy_gives_full_confidence(self):
        """Entropy of 0 should yield confidence of 100."""
        conf = entropy_to_confidence(0.0)
        assert float(conf) == pytest.approx(100.0, abs=1e-6)

    def test_max_entropy_gives_zero_confidence(self):
        """Entropy equal to log2(6) should yield confidence near 0."""
        conf = entropy_to_confidence(_MAX_ENTROPY_6)
        assert float(conf) == pytest.approx(0.0, abs=1e-6)

    def test_intermediate_entropy(self):
        """Intermediate entropy should produce proportional confidence."""
        half_max = _MAX_ENTROPY_6 / 2.0
        conf = entropy_to_confidence(half_max)
        assert float(conf) == pytest.approx(50.0, abs=1e-6)

    def test_numpy_array_input(self):
        """Should accept and return numpy arrays."""
        entropies = np.array([0.0, _MAX_ENTROPY_6 / 2.0, _MAX_ENTROPY_6])
        confidences = entropy_to_confidence(entropies)

        assert isinstance(confidences, np.ndarray)
        assert confidences[0] == pytest.approx(100.0, abs=1e-6)
        assert confidences[1] == pytest.approx(50.0, abs=1e-6)
        assert confidences[2] == pytest.approx(0.0, abs=1e-6)


class TestBuildCategoryParams:
    """Tests for build_category_params() parameter factory."""

    def test_default_config_produces_six_categories(self):
        """Default ScoringConfig should produce exactly 6 categories."""
        config = ScoringConfig()
        params = build_category_params(config)

        assert len(params) == 6
        names = [p.name for p in params]
        assert names == [
            "Known Species",
            "Novel Species",
            "Novel Genus",
            "Species Boundary",
            "Ambiguous",
            "Unclassified",
        ]

    def test_known_species_mean_is_half_known_max(self):
        """Known Species novelty mean should equal novelty_known_max / 2."""
        config = ScoringConfig()
        params = build_category_params(config)
        known = params[0]

        expected_mean = config.novelty_known_max / 2.0
        assert known.novelty_mean == pytest.approx(expected_mean, abs=1e-6)

    def test_all_sigmas_positive(self):
        """All category sigma values should be strictly positive."""
        config = ScoringConfig()
        params = build_category_params(config)

        for p in params:
            assert p.novelty_sigma > 0, f"{p.name} novelty_sigma should be positive"
            assert p.uncertainty_sigma > 0, f"{p.name} uncertainty_sigma should be positive"

    def test_custom_config_shifts_centers(self):
        """Custom thresholds should shift category center positions."""
        default_params = build_category_params(ScoringConfig())
        custom_params = build_category_params(ScoringConfig(
            novelty_known_max=8.0,
            novelty_novel_species_min=8.0,
            novelty_novel_species_max=20.0,
            novelty_novel_genus_min=20.0,
        ))

        # Known Species mean should shift
        assert custom_params[0].novelty_mean == pytest.approx(4.0, abs=1e-6)  # 8.0 / 2
        assert custom_params[0].novelty_mean > default_params[0].novelty_mean  # 4.0 > 2.0

        # Novel Species mean should shift (8+20)/2 = 14 vs (4+20)/2 = 12
        assert custom_params[1].novelty_mean == pytest.approx(14.0, abs=1e-6)
        assert custom_params[1].novelty_mean > default_params[1].novelty_mean

    def test_category_params_are_frozen(self):
        """CategoryParams dataclass instances should be frozen (immutable)."""
        config = ScoringConfig()
        params = build_category_params(config)
        with pytest.raises(AttributeError):
            params[0].novelty_mean = 999.0


class TestPriorModulation:
    """Tests for Bayesian prior modulation mechanics in BayesianClassifier."""

    @pytest.fixture
    def default_classifier(self):
        return BayesianClassifier(ScoringConfig())

    def test_identity_gap_boost_increases_ambiguous(self, default_classifier):
        """Small identity gap + multiple hits should boost Ambiguous posterior."""
        # identity_gap < 2.0 (default threshold) AND num_hits > 1
        result_boosted = default_classifier.compute_posteriors(
            novelty_index=5.0,
            placement_uncertainty=1.0,
            num_hits=3,
            identity_gap=0.5,  # below threshold of 2.0
        )
        # identity_gap above threshold -> no boost
        result_no_boost = default_classifier.compute_posteriors(
            novelty_index=5.0,
            placement_uncertainty=1.0,
            num_hits=3,
            identity_gap=5.0,  # above threshold of 2.0
        )

        assert result_boosted.p_ambiguous > result_no_boost.p_ambiguous

    def test_single_hit_boost_for_novel_reads(self, default_classifier):
        """Single hit + novelty >= novel_species_min should boost Ambiguous."""
        # num_hits=1, novelty above novel_species_min (4.0 by default)
        result_single_novel = default_classifier.compute_posteriors(
            novelty_index=10.0,
            placement_uncertainty=1.0,
            num_hits=1,
        )
        # num_hits > 1 with same novelty (no single-hit boost)
        result_multi = default_classifier.compute_posteriors(
            novelty_index=10.0,
            placement_uncertainty=1.0,
            num_hits=3,
            identity_gap=10.0,  # large gap to avoid identity_gap boost
        )

        assert result_single_novel.p_ambiguous > result_multi.p_ambiguous

    def test_custom_higher_boost_values_increase_effect(self):
        """Higher boost values in BayesianConfig should amplify the effect."""
        from metadarkmatter.models.config import BayesianConfig

        default_clf = BayesianClassifier(ScoringConfig())
        boosted_clf = BayesianClassifier(ScoringConfig(
            bayesian=BayesianConfig(
                identity_gap_boost=5.0,
                single_hit_boost=4.0,
            ),
        ))

        # Identity gap boost comparison
        default_result = default_clf.compute_posteriors(
            novelty_index=5.0,
            placement_uncertainty=1.0,
            num_hits=3,
            identity_gap=0.5,
        )
        boosted_result = boosted_clf.compute_posteriors(
            novelty_index=5.0,
            placement_uncertainty=1.0,
            num_hits=3,
            identity_gap=0.5,
        )
        assert boosted_result.p_ambiguous > default_result.p_ambiguous

        # Single hit boost comparison
        default_single = default_clf.compute_posteriors(
            novelty_index=10.0,
            placement_uncertainty=1.0,
            num_hits=1,
        )
        boosted_single = boosted_clf.compute_posteriors(
            novelty_index=10.0,
            placement_uncertainty=1.0,
            num_hits=1,
        )
        assert boosted_single.p_ambiguous > default_single.p_ambiguous

    def test_no_boost_when_conditions_not_met(self, default_classifier):
        """No prior modulation should occur when boost conditions are not met."""
        # Large identity gap + multiple hits -> no identity_gap boost
        # num_hits > 1 + low novelty -> no single_hit boost
        result_large_gap = default_classifier.compute_posteriors(
            novelty_index=2.0,
            placement_uncertainty=0.5,
            num_hits=3,
            identity_gap=10.0,  # well above 2.0 threshold
        )
        # No identity_gap at all (None) -> no identity_gap boost
        result_no_gap = default_classifier.compute_posteriors(
            novelty_index=2.0,
            placement_uncertainty=0.5,
            num_hits=3,
        )

        # Both should produce the same posteriors (no boost applied in either case)
        assert result_large_gap.p_ambiguous == pytest.approx(
            result_no_gap.p_ambiguous, abs=1e-6
        )

    def test_single_hit_no_boost_below_novel_threshold(self, default_classifier):
        """Single hit with low novelty should not receive single_hit_boost."""
        # novelty below novel_species_min (4.0 default) -> no single_hit_boost
        result_single_low = default_classifier.compute_posteriors(
            novelty_index=2.0,
            placement_uncertainty=0.5,
            num_hits=1,
        )
        # Multiple hits with same parameters (no boosts expected for either)
        result_multi_low = default_classifier.compute_posteriors(
            novelty_index=2.0,
            placement_uncertainty=0.5,
            num_hits=3,
            identity_gap=10.0,  # large gap to avoid identity_gap boost
        )

        # With low novelty, single-hit should not receive a boost,
        # so ambiguous posterior should be comparable or lower
        # (no single_hit_boost because novelty < novel_species_min)
        assert result_single_low.p_ambiguous <= result_multi_low.p_ambiguous + 0.01


class TestClassifyPrimary:
    """Tests for BayesianClassifier.classify_primary().

    The classify_primary method is the Bayesian-primary entry point that
    computes 6-category posteriors, applies Stage 2 refinement, and derives
    a confidence_score from posterior entropy.
    """

    _EXPECTED_COLUMNS = [
        "taxonomic_call",
        "p_known_species",
        "p_novel_species",
        "p_novel_genus",
        "p_species_boundary",
        "p_ambiguous",
        "p_unclassified",
        "posterior_entropy",
        "confidence_score",
    ]

    _POSTERIOR_COLUMNS = [
        "p_known_species",
        "p_novel_species",
        "p_novel_genus",
        "p_species_boundary",
        "p_ambiguous",
        "p_unclassified",
    ]

    @pytest.fixture
    def classifier(self):
        config = ScoringConfig()
        return BayesianClassifier(config)

    def test_classify_primary_adds_all_columns(self, classifier):
        """All expected output columns should be present after classify_primary."""
        df = pl.DataFrame({
            "novelty_index": [2.0, 12.0, 22.0],
            "placement_uncertainty": [0.5, 0.5, 0.5],
        })

        result = classifier.classify_primary(df)

        for col in self._EXPECTED_COLUMNS:
            assert col in result.columns, f"Missing expected column: {col}"

    def test_classify_primary_posteriors_sum_to_one(self, classifier):
        """Posterior probabilities should sum to 1.0 for each row."""
        df = pl.DataFrame({
            "novelty_index": [1.0, 5.0, 12.0, 22.0, 30.0],
            "placement_uncertainty": [0.3, 1.0, 0.5, 0.5, 8.0],
            "num_ambiguous_hits": [1, 3, 3, 3, 2],
            "identity_gap": [None, 5.0, 5.0, 5.0, 1.0],
        })

        result = classifier.classify_primary(df)

        sums = _sum_all_6_df(result)
        for i, s in enumerate(sums.to_list()):
            assert s == pytest.approx(1.0, abs=0.01), (
                f"Row {i}: posteriors sum to {s}, expected ~1.0"
            )

    def test_classify_primary_known_species_classification(self, classifier):
        """Low novelty with low uncertainty should classify as Known Species."""
        df = pl.DataFrame({
            "novelty_index": [1.0],
            "placement_uncertainty": [0.3],
            "num_ambiguous_hits": [3],
            "identity_gap": [8.0],
        })

        result = classifier.classify_primary(df)

        assert result["taxonomic_call"][0] == "Known Species"

    def test_classify_primary_novel_species_classification(self, classifier):
        """Moderate novelty with low uncertainty should classify as Novel Species."""
        df = pl.DataFrame({
            "novelty_index": [12.0],
            "placement_uncertainty": [0.5],
            "num_ambiguous_hits": [3],
            "identity_gap": [8.0],
        })

        result = classifier.classify_primary(df)

        assert result["taxonomic_call"][0] == "Novel Species"

    def test_classify_primary_novel_genus_classification(self, classifier):
        """High novelty with low uncertainty should classify as Novel Genus."""
        df = pl.DataFrame({
            "novelty_index": [22.0],
            "placement_uncertainty": [0.5],
            "num_ambiguous_hits": [3],
            "identity_gap": [8.0],
        })

        result = classifier.classify_primary(df)

        assert result["taxonomic_call"][0] == "Novel Genus"

    def test_classify_primary_conserved_region_refinement(self, classifier):
        """Ambiguous + ambiguity_scope='across_genera' should refine to Conserved Region."""
        # High uncertainty drives Stage 1 MAP to Ambiguous; Stage 2 refines
        # to Conserved Region when ambiguity_scope is across_genera.
        df = pl.DataFrame({
            "novelty_index": [15.0],
            "placement_uncertainty": [8.0],
            "num_ambiguous_hits": [5],
            "identity_gap": [0.5],
            "ambiguity_scope": ["across_genera"],
        })

        result = classifier.classify_primary(df)

        assert result["taxonomic_call"][0] == "Conserved Region"

    def test_classify_primary_confidence_score_range(self, classifier):
        """confidence_score should be between 0 and 100 for all rows."""
        rng = np.random.default_rng(99)
        n = 50
        df = pl.DataFrame({
            "novelty_index": rng.uniform(0, 35, n).tolist(),
            "placement_uncertainty": rng.uniform(0, 12, n).tolist(),
        })

        result = classifier.classify_primary(df)

        scores = result["confidence_score"].to_list()
        for i, score in enumerate(scores):
            assert 0.0 <= score <= 100.0, (
                f"Row {i}: confidence_score={score} outside [0, 100]"
            )

    def test_classify_primary_empty_dataframe(self, classifier):
        """Empty input should return empty output with all expected columns."""
        df = pl.DataFrame({
            "novelty_index": [],
            "placement_uncertainty": [],
        }).cast({
            "novelty_index": pl.Float64,
            "placement_uncertainty": pl.Float64,
        })

        result = classifier.classify_primary(df)

        assert result.height == 0
        for col in self._EXPECTED_COLUMNS:
            assert col in result.columns, f"Missing column in empty result: {col}"

    def test_classify_primary_side_by_side_with_legacy(self, classifier):
        """Regression test: both taxonomic_call and legacy_call should coexist."""
        from metadarkmatter.core.classification.thresholds import apply_legacy_thresholds

        config = ScoringConfig()

        # Build a DataFrame with columns required by both classify_primary
        # and apply_legacy_thresholds.
        df = pl.DataFrame({
            "novelty_index": [1.0, 12.0, 22.0, 5.0, 15.0],
            "placement_uncertainty": [0.3, 0.5, 0.5, 8.0, 0.5],
            "num_ambiguous_hits": [3, 3, 3, 5, 1],
            "identity_gap": [8.0, 8.0, 8.0, 0.5, 8.0],
            "ambiguity_scope": [
                "within_species",
                "within_species",
                "within_species",
                "across_genera",
                "within_species",
            ],
            "genus_uncertainty": [0.0, 0.0, 0.0, 0.0, 0.0],
            "num_secondary_genomes": [0, 0, 0, 0, 0],
        })

        # Run Bayesian-primary classification
        result = classifier.classify_primary(df)

        # Run legacy thresholds on the Bayesian-classified result
        result = apply_legacy_thresholds(result, config)

        # Both columns should be present
        assert "taxonomic_call" in result.columns
        assert "legacy_call" in result.columns

        # All values should be valid category strings (not null)
        valid_taxonomic_calls = {
            "Known Species",
            "Novel Species",
            "Novel Genus",
            "Species Boundary",
            "Ambiguous",
            "Unclassified",
            "Conserved Region",
            "Ambiguous Within Genus",
            "Off-target",
        }
        valid_legacy_calls = {
            "Known Species",
            "Novel Species",
            "Novel Genus",
            "Species Boundary",
            "Ambiguous",
            "Unclassified",
            "Conserved Region",
        }

        for val in result["taxonomic_call"].to_list():
            assert val in valid_taxonomic_calls, (
                f"Unexpected taxonomic_call value: {val}"
            )
        for val in result["legacy_call"].to_list():
            assert val in valid_legacy_calls, (
                f"Unexpected legacy_call value: {val}"
            )
