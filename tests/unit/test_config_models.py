"""Unit tests for new configuration models.

Tests for ExtractConfig, MapConfig, and VisualizeConfig.
"""

from __future__ import annotations

import pytest

from pydantic import ValidationError

from metadarkmatter.models.config import (
    ExtractConfig,
    MapConfig,
    ScoringConfig,
    VisualizeConfig,
)


class TestExtractConfig:
    """Tests for ExtractConfig model."""

    def test_default_values(self):
        """Should use sensible defaults."""
        config = ExtractConfig()
        assert config.include_children is True
        assert config.confidence == 0.0

    def test_custom_values(self):
        """Should accept custom values."""
        config = ExtractConfig(
            include_children=False,
            confidence=0.5,
        )
        assert config.include_children is False
        assert config.confidence == 0.5

    def test_confidence_range(self):
        """Should validate confidence is 0-1."""
        with pytest.raises(ValueError):
            ExtractConfig(confidence=1.5)

        with pytest.raises(ValueError):
            ExtractConfig(confidence=-0.1)

    def test_frozen(self):
        """Should be immutable."""
        config = ExtractConfig()
        with pytest.raises(Exception):  # ValidationError or AttributeError
            config.include_children = False


class TestMapConfig:
    """Tests for MapConfig model."""

    def test_default_values(self):
        """Should use sensible defaults."""
        config = MapConfig()
        assert config.mode == "local"
        assert config.very_sensitive is True
        assert config.max_alignments == 500
        assert config.no_unal is True
        assert config.no_discordant is True

    def test_end_to_end_mode(self):
        """Should accept end-to-end mode."""
        config = MapConfig(mode="end-to-end")
        assert config.mode == "end-to-end"

    def test_invalid_mode_rejected(self):
        """Should reject invalid alignment modes."""
        with pytest.raises(ValueError):
            MapConfig(mode="invalid")

    def test_max_alignments_positive(self):
        """Should require positive max_alignments."""
        with pytest.raises(ValueError):
            MapConfig(max_alignments=0)


class TestVisualizeConfig:
    """Tests for VisualizeConfig model."""

    def test_default_values(self):
        """Should use sensible defaults."""
        config = VisualizeConfig()
        assert config.min_identity == 70.0
        assert config.max_points == 100000
        assert config.known_species_threshold == 98.0
        assert config.novel_species_threshold == 85.0
        assert config.novel_genus_threshold == 75.0

    def test_custom_thresholds(self):
        """Should accept custom identity thresholds."""
        config = VisualizeConfig(
            known_species_threshold=99.0,
            novel_species_threshold=90.0,
            novel_genus_threshold=80.0,
        )
        assert config.known_species_threshold == 99.0
        assert config.novel_species_threshold == 90.0
        assert config.novel_genus_threshold == 80.0

    def test_identity_range(self):
        """Should validate identity thresholds are 0-100."""
        with pytest.raises(ValueError):
            VisualizeConfig(min_identity=101.0)

        with pytest.raises(ValueError):
            VisualizeConfig(known_species_threshold=-1.0)

    def test_max_points_minimum(self):
        """Should require at least 1000 points."""
        with pytest.raises(ValueError):
            VisualizeConfig(max_points=500)

    def test_frozen(self):
        """Should be immutable."""
        config = VisualizeConfig()
        with pytest.raises(Exception):
            config.min_identity = 80.0


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
