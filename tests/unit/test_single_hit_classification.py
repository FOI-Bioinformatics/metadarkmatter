"""
Unit tests for single-hit classification behavior.

Tests the Rule 0b implementation that uses inferred uncertainty
for single-hit reads to reduce novel species overestimation.
"""

from __future__ import annotations

import pytest

from metadarkmatter.core.ani_placement import ANIMatrix, ANIWeightedClassifier
from metadarkmatter.core.constants import calculate_inferred_uncertainty
from metadarkmatter.models.classification import TaxonomicCall
from metadarkmatter.models.config import ScoringConfig


class TestInferredUncertaintyCalculation:
    """Tests for the calculate_inferred_uncertainty function."""

    def test_low_novelty_returns_low_uncertainty(self):
        """Low novelty (known species range) should return low inferred uncertainty."""
        # Novelty < 5% (known species range)
        uncertainty = calculate_inferred_uncertainty(3.0)
        # Should be around 5-7.5% (base + small increment)
        assert 5.0 <= uncertainty <= 7.5

    def test_novel_species_range_returns_moderate_uncertainty(self):
        """Novel species range should return moderate inferred uncertainty."""
        # Novelty 10% (novel species range)
        uncertainty = calculate_inferred_uncertainty(10.0)
        # Should be around 12-15% (increasing with novelty)
        assert 10.0 <= uncertainty <= 18.0

    def test_novel_genus_range_returns_high_uncertainty(self):
        """Novel genus range should return high inferred uncertainty."""
        # Novelty 20% (novel genus range)
        uncertainty = calculate_inferred_uncertainty(20.0)
        # Should be around 17-25%
        assert 15.0 <= uncertainty <= 25.0

    def test_very_high_novelty_returns_max_uncertainty(self):
        """Very high novelty should return maximum inferred uncertainty."""
        # Novelty > 25%
        uncertainty = calculate_inferred_uncertainty(30.0)
        # Should be maximum (35%)
        assert uncertainty == 35.0

    def test_uncertainty_monotonically_increases(self):
        """Inferred uncertainty should increase with novelty."""
        uncertainties = [
            calculate_inferred_uncertainty(n) for n in [2.0, 5.0, 10.0, 15.0, 20.0, 25.0]
        ]
        for i in range(1, len(uncertainties)):
            assert uncertainties[i] >= uncertainties[i - 1], (
                f"Uncertainty should increase: {uncertainties[i-1]} -> {uncertainties[i]}"
            )


class TestSingleHitConfigParameters:
    """Tests for single-hit classification config parameters."""

    def test_default_use_inferred_for_single_hits_is_false(self):
        """Default config should have use_inferred_for_single_hits=False."""
        config = ScoringConfig()
        assert config.use_inferred_for_single_hits is False

    def test_default_single_hit_uncertainty_threshold(self):
        """Default threshold should be 10%."""
        config = ScoringConfig()
        assert config.single_hit_uncertainty_threshold == 10.0

    def test_can_enable_single_hit_gating(self):
        """Should be able to enable single-hit gating."""
        config = ScoringConfig(use_inferred_for_single_hits=True)
        assert config.use_inferred_for_single_hits is True

    def test_can_set_custom_threshold(self):
        """Should be able to set custom threshold."""
        config = ScoringConfig(
            use_inferred_for_single_hits=True,
            single_hit_uncertainty_threshold=15.0,
        )
        assert config.single_hit_uncertainty_threshold == 15.0

    def test_threshold_validation_min(self):
        """Threshold should not be negative."""
        from pydantic import ValidationError

        with pytest.raises(ValidationError):
            ScoringConfig(single_hit_uncertainty_threshold=-5.0)

    def test_threshold_validation_max(self):
        """Threshold should not exceed 100."""
        from pydantic import ValidationError

        with pytest.raises(ValidationError):
            ScoringConfig(single_hit_uncertainty_threshold=150.0)


class TestSingleHitClassification:
    """Tests for single-hit classification behavior with Rule 0b."""

    @pytest.fixture
    def ani_matrix(self, small_ani_dict) -> ANIMatrix:
        """Create ANI matrix from small_ani_dict fixture."""
        return ANIMatrix(small_ani_dict)

    def test_default_behavior_unchanged(self, ani_matrix):
        """Default config should produce Novel Species for single-hit reads in novel range."""
        # Default config: use_inferred_for_single_hits=False
        config = ScoringConfig()
        classifier = ANIWeightedClassifier(ani_matrix, config=config)

        # Single-hit read with novelty_index=12% (novel species range)
        # placement_uncertainty=0% (single hit, no competing genomes)
        call = classifier._classify_by_thresholds(
            novelty_index=12.0,
            placement_uncertainty=0.0,
            num_ambiguous_hits=1,
        )

        # Should be Novel Species (default behavior unchanged)
        assert call == TaxonomicCall.NOVEL_SPECIES

    def test_single_hit_becomes_ambiguous_when_enabled(self, ani_matrix):
        """Single-hit Novel Species with high inferred uncertainty should become Ambiguous."""
        config = ScoringConfig(
            use_inferred_for_single_hits=True,
            single_hit_uncertainty_threshold=10.0,
        )
        classifier = ANIWeightedClassifier(ani_matrix, config=config)

        # Single-hit read with novelty_index=12%
        # Inferred uncertainty at 12% novelty is ~12.5% (7.5 + 5*1.0)
        # 12.5% >= 10% threshold -> should be Ambiguous
        call = classifier._classify_by_thresholds(
            novelty_index=12.0,
            placement_uncertainty=0.0,
            num_ambiguous_hits=1,
        )

        assert call == TaxonomicCall.AMBIGUOUS

    def test_single_hit_known_species_unaffected(self, ani_matrix):
        """Single-hit Known Species reads should not be affected by Rule 0b."""
        config = ScoringConfig(
            use_inferred_for_single_hits=True,
            single_hit_uncertainty_threshold=10.0,
        )
        classifier = ANIWeightedClassifier(ani_matrix, config=config)

        # Single-hit read with novelty_index=3% (known species range)
        # Rule 0b only applies to novel range (novelty >= 5%)
        call = classifier._classify_by_thresholds(
            novelty_index=3.0,
            placement_uncertainty=0.0,
            num_ambiguous_hits=1,
        )

        # Should remain Known Species (Rule 0b only applies to novel range)
        assert call == TaxonomicCall.KNOWN_SPECIES

    def test_multi_hit_reads_unaffected(self, ani_matrix):
        """Multi-hit reads should not be affected by Rule 0b."""
        config = ScoringConfig(
            use_inferred_for_single_hits=True,
            single_hit_uncertainty_threshold=10.0,
        )
        classifier = ANIWeightedClassifier(ani_matrix, config=config)

        # Multi-hit read with novelty_index=12%
        # Rule 0b only applies to single-hit reads (num_ambiguous_hits <= 1)
        call = classifier._classify_by_thresholds(
            novelty_index=12.0,
            placement_uncertainty=0.5,  # Low uncertainty from ANI
            num_ambiguous_hits=3,
        )

        # Should be Novel Species (Rule 0b doesn't apply to multi-hit)
        assert call == TaxonomicCall.NOVEL_SPECIES

    def test_threshold_boundary_exact(self, ani_matrix):
        """Test behavior at exact threshold boundary."""
        config = ScoringConfig(
            use_inferred_for_single_hits=True,
            single_hit_uncertainty_threshold=12.5,  # Set threshold exactly
        )
        classifier = ANIWeightedClassifier(ani_matrix, config=config)

        # Novelty=12% gives inferred uncertainty ~12.5%
        # At exact boundary (>=), should be Ambiguous
        call = classifier._classify_by_thresholds(
            novelty_index=12.0,
            placement_uncertainty=0.0,
            num_ambiguous_hits=1,
        )

        assert call == TaxonomicCall.AMBIGUOUS

    def test_below_threshold_remains_novel(self, ani_matrix):
        """Single-hit below threshold should remain Novel Species."""
        config = ScoringConfig(
            use_inferred_for_single_hits=True,
            single_hit_uncertainty_threshold=20.0,  # Higher threshold
        )
        classifier = ANIWeightedClassifier(ani_matrix, config=config)

        # Novelty=8% gives inferred uncertainty ~10.5%
        # 10.5% < 20% threshold -> should remain Novel Species
        call = classifier._classify_by_thresholds(
            novelty_index=8.0,
            placement_uncertainty=0.0,
            num_ambiguous_hits=1,
        )

        assert call == TaxonomicCall.NOVEL_SPECIES

    def test_novel_genus_also_affected(self, ani_matrix):
        """Single-hit Novel Genus range should also be affected by Rule 0b."""
        config = ScoringConfig(
            use_inferred_for_single_hits=True,
            single_hit_uncertainty_threshold=10.0,
        )
        classifier = ANIWeightedClassifier(ani_matrix, config=config)

        # Single-hit read with novelty_index=22% (novel genus range)
        # Inferred uncertainty at 22% is high (~20%)
        call = classifier._classify_by_thresholds(
            novelty_index=22.0,
            placement_uncertainty=0.0,
            num_ambiguous_hits=1,
        )

        # Should be Ambiguous (high inferred uncertainty)
        assert call == TaxonomicCall.AMBIGUOUS


class TestBalancedConservativePreset:
    """Tests for the balanced-conservative preset."""

    def test_preset_exists(self):
        """balanced-conservative preset should exist."""
        from metadarkmatter.cli.score import THRESHOLD_PRESETS

        assert "balanced-conservative" in THRESHOLD_PRESETS

    def test_preset_enables_single_hit_gating(self):
        """balanced-conservative preset should enable single-hit gating."""
        from metadarkmatter.cli.score import THRESHOLD_PRESETS

        config = THRESHOLD_PRESETS["balanced-conservative"]
        assert config.use_inferred_for_single_hits is True

    def test_preset_has_stricter_thresholds(self):
        """balanced-conservative preset should have stricter thresholds than default."""
        from metadarkmatter.cli.score import THRESHOLD_PRESETS

        default = THRESHOLD_PRESETS["default"]
        conservative = THRESHOLD_PRESETS["balanced-conservative"]

        # Stricter novelty thresholds (96% species boundary)
        assert conservative.novelty_known_max <= default.novelty_known_max

        # Stricter uncertainty thresholds
        assert conservative.uncertainty_known_max <= default.uncertainty_known_max

        # Lower bitscore threshold to capture more competing hits
        assert conservative.bitscore_threshold_pct < default.bitscore_threshold_pct

    def test_preset_single_hit_threshold(self):
        """balanced-conservative preset should have specific single-hit threshold."""
        from metadarkmatter.cli.score import THRESHOLD_PRESETS

        config = THRESHOLD_PRESETS["balanced-conservative"]
        assert config.single_hit_uncertainty_threshold == 12.0

    def test_preset_alignment_fraction(self):
        """balanced-conservative preset should have moderate alignment fraction."""
        from metadarkmatter.cli.score import THRESHOLD_PRESETS

        config = THRESHOLD_PRESETS["balanced-conservative"]
        assert config.min_alignment_fraction == 0.3
