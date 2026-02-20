"""Unit tests for configuration model features.

Tests for family validation fields in ScoringConfig.
"""

from __future__ import annotations

import pytest

from pydantic import ValidationError

from metadarkmatter.models.config import ScoringConfig


class TestFamilyValidationConfig:
    """Tests for family validation config fields."""

    def test_default_target_family_is_none(self):
        """target_family should default to None (disabled)."""
        config = ScoringConfig()
        assert config.target_family is None

    def test_default_family_ratio_threshold(self):
        """family_ratio_threshold should default to 0.8."""
        config = ScoringConfig()
        assert config.family_ratio_threshold == 0.8

    def test_target_family_accepts_string(self):
        """target_family should accept a family name string."""
        config = ScoringConfig(target_family="f__Francisellaceae")
        assert config.target_family == "f__Francisellaceae"

    def test_family_ratio_threshold_custom(self):
        """family_ratio_threshold should accept custom values."""
        config = ScoringConfig(family_ratio_threshold=0.9)
        assert config.family_ratio_threshold == 0.9

    def test_family_ratio_threshold_too_high(self):
        """family_ratio_threshold > 1.0 should fail validation."""
        with pytest.raises(ValidationError):
            ScoringConfig(family_ratio_threshold=1.5)

    def test_family_ratio_threshold_negative(self):
        """family_ratio_threshold < 0 should fail validation."""
        with pytest.raises(ValidationError):
            ScoringConfig(family_ratio_threshold=-0.1)
