"""
Unit tests for threshold sensitivity analysis.

Tests apply_classification_thresholds and run_sensitivity_analysis
for correctness and expected behavior across threshold sweeps.
"""

from __future__ import annotations

import polars as pl
import pytest

from metadarkmatter.core.classification.sensitivity import (
    SensitivityResult,
    run_sensitivity_analysis,
)
from metadarkmatter.core.classification.thresholds import (
    apply_classification_thresholds,
)
from metadarkmatter.models.config import ScoringConfig


class TestApplyClassificationThresholds:
    """Tests for the apply_classification_thresholds function."""

    @pytest.fixture
    def metrics_df(self):
        """Pre-computed metrics DataFrame for threshold testing."""
        return pl.DataFrame({
            "novelty_index": [1.0, 2.0, 10.0, 15.0, 22.0, 30.0],
            "placement_uncertainty": [0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
            "num_ambiguous_hits": [3, 3, 3, 3, 3, 3],
            "identity_gap": [5.0, 5.0, 5.0, 5.0, 5.0, 5.0],
        })

    def test_known_species(self, metrics_df):
        """Low novelty with low uncertainty should be Known Species."""
        config = ScoringConfig()
        result = apply_classification_thresholds(metrics_df, config)

        calls = result["taxonomic_call"].to_list()
        # novelty=1.0, uncertainty=0.5 -> Known Species
        assert calls[0] == "Known Species"

    def test_novel_species(self, metrics_df):
        """Moderate novelty should be Novel Species."""
        config = ScoringConfig()
        result = apply_classification_thresholds(metrics_df, config)

        calls = result["taxonomic_call"].to_list()
        # novelty=10.0, uncertainty=0.5 -> Novel Species
        assert calls[2] == "Novel Species"

    def test_novel_genus(self, metrics_df):
        """High novelty should be Novel Genus."""
        config = ScoringConfig()
        result = apply_classification_thresholds(metrics_df, config)

        calls = result["taxonomic_call"].to_list()
        # novelty=22.0, uncertainty=0.5 -> Novel Genus
        assert calls[4] == "Novel Genus"

    def test_threshold_shift_changes_classification(self, metrics_df):
        """Changing novelty_known_max should shift classification boundaries."""
        # With default (4.0): novelty=2.0 -> Known Species
        config_default = ScoringConfig()
        result_default = apply_classification_thresholds(metrics_df, config_default)
        assert result_default["taxonomic_call"][1] == "Known Species"

        # With strict (1.5): novelty=2.0 -> Novel Species
        config_strict = ScoringConfig(
            novelty_known_max=1.5,
            novelty_novel_species_min=1.5,
        )
        result_strict = apply_classification_thresholds(metrics_df, config_strict)
        assert result_strict["taxonomic_call"][1] == "Novel Species"

    def test_diversity_status_column_added(self, metrics_df):
        """Should add diversity_status column."""
        config = ScoringConfig()
        result = apply_classification_thresholds(metrics_df, config)

        assert "diversity_status" in result.columns
        assert "is_novel" in result.columns

    def test_is_novel_flag(self, metrics_df):
        """is_novel should be True for Novel Species and Novel Genus."""
        config = ScoringConfig()
        result = apply_classification_thresholds(metrics_df, config)

        novels = result["is_novel"].to_list()
        calls = result["taxonomic_call"].to_list()

        for novel, call in zip(novels, calls):
            if call in ("Novel Species", "Novel Genus"):
                assert novel is True
            else:
                assert novel is False

    def test_high_uncertainty_gives_ambiguous(self):
        """High placement uncertainty should give Ambiguous."""
        df = pl.DataFrame({
            "novelty_index": [10.0],
            "placement_uncertainty": [8.0],
            "num_ambiguous_hits": [3],
            "identity_gap": [5.0],
        })
        config = ScoringConfig()
        result = apply_classification_thresholds(df, config)

        assert result["taxonomic_call"][0] in ("Ambiguous", "Conserved Region", "Species Boundary")

    def test_identity_gap_ambiguous(self):
        """Small identity gap should give Ambiguous."""
        df = pl.DataFrame({
            "novelty_index": [2.0],
            "placement_uncertainty": [0.5],
            "num_ambiguous_hits": [3],
            "identity_gap": [0.5],  # Very small gap
        })
        config = ScoringConfig()
        result = apply_classification_thresholds(df, config)

        assert result["taxonomic_call"][0] == "Ambiguous"


class TestSensitivityAnalysis:
    """Tests for run_sensitivity_analysis function."""

    @pytest.fixture
    def metrics_df(self):
        """DataFrame with reads at various novelty levels."""
        return pl.DataFrame({
            "novelty_index": [1.0, 3.0, 5.0, 10.0, 15.0, 22.0],
            "placement_uncertainty": [0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
            "num_ambiguous_hits": [3, 3, 3, 3, 3, 3],
            "identity_gap": [5.0, 5.0, 5.0, 5.0, 5.0, 5.0],
        })

    def test_returns_sensitivity_result(self, metrics_df):
        """Should return SensitivityResult dataclass."""
        config = ScoringConfig()
        result = run_sensitivity_analysis(metrics_df, config, steps=5)

        assert isinstance(result, SensitivityResult)
        assert len(result.novelty_thresholds) == 5
        assert len(result.uncertainty_thresholds) == 5

    def test_counts_sum_to_total(self, metrics_df):
        """At each threshold, total counts should equal number of reads."""
        config = ScoringConfig()
        result = run_sensitivity_analysis(metrics_df, config, steps=5)

        total_reads = metrics_df.height
        for i in range(5):
            step_total = sum(
                result.counts[cat][i] for cat in result.counts
            )
            assert step_total == total_reads

    def test_known_species_increases_with_relaxed_threshold(self, metrics_df):
        """More reads should be Known Species when novelty threshold is higher."""
        config = ScoringConfig()
        result = run_sensitivity_analysis(
            metrics_df, config,
            novelty_range=(2.0, 8.0),
            steps=5,
        )

        known_counts = result.counts["Known Species"]
        # Known Species should generally increase as threshold increases
        assert known_counts[-1] >= known_counts[0]

    def test_to_dict_serializable(self, metrics_df):
        """to_dict should produce JSON-serializable output."""
        import json

        config = ScoringConfig()
        result = run_sensitivity_analysis(metrics_df, config, steps=3)

        d = result.to_dict()
        # Should not raise
        json_str = json.dumps(d)
        assert "novelty_thresholds" in json_str

    def test_custom_ranges(self, metrics_df):
        """Should respect custom novelty and uncertainty ranges."""
        config = ScoringConfig()
        result = run_sensitivity_analysis(
            metrics_df, config,
            novelty_range=(1.0, 10.0),
            uncertainty_range=(0.5, 5.0),
            steps=3,
        )

        assert result.novelty_thresholds[0] == pytest.approx(1.0)
        assert result.novelty_thresholds[-1] == pytest.approx(10.0)
        assert result.uncertainty_thresholds[0] == pytest.approx(0.5)
        assert result.uncertainty_thresholds[-1] == pytest.approx(5.0)
