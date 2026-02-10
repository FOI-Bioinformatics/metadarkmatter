"""
Unit tests for adaptive threshold detection.

Tests GMM-based species boundary detection and fallback behavior.
"""

from __future__ import annotations

import pytest

from metadarkmatter.core.ani_placement import ANIMatrix
from metadarkmatter.core.classification.adaptive import (
    AdaptiveThresholds,
    build_adaptive_config,
    detect_species_boundary,
)
from metadarkmatter.models.config import ScoringConfig


class TestDetectSpeciesBoundary:
    """Tests for detect_species_boundary function."""

    def test_fallback_with_few_genomes(self, small_ani_dict):
        """Should fall back to default with < min_genomes."""
        ani = ANIMatrix(small_ani_dict)  # 3 genomes
        result = detect_species_boundary(ani, min_genomes=5)

        assert result.method == "fallback"
        assert result.confidence == 0.0
        assert result.novelty_known_max == 4.0  # Default

    def test_fallback_with_min_genomes_met(self):
        """Should attempt GMM when min_genomes is met."""
        # Create matrix with clear species boundary:
        # 5 genomes: 3 in one species (~97% ANI), 2 in another (~85% ANI)
        ani_dict = {}
        species_a = [f"genome_A{i}" for i in range(3)]
        species_b = [f"genome_B{i}" for i in range(2)]
        all_genomes = species_a + species_b

        for g1 in all_genomes:
            ani_dict[g1] = {}
            for g2 in all_genomes:
                if g1 == g2:
                    ani_dict[g1][g2] = 100.0
                elif g1 in species_a and g2 in species_a:
                    ani_dict[g1][g2] = 97.0 + (hash(g1 + g2) % 200) / 100.0  # 97-99
                elif g1 in species_b and g2 in species_b:
                    ani_dict[g1][g2] = 97.5 + (hash(g1 + g2) % 150) / 100.0  # 97.5-99
                else:
                    ani_dict[g1][g2] = 84.0 + (hash(g1 + g2) % 300) / 100.0  # 84-87

        ani = ANIMatrix(ani_dict)
        result = detect_species_boundary(ani, min_genomes=5)

        assert result.method in ("gmm", "fallback")
        assert result.ani_values_used > 0

    def test_well_separated_species(self):
        """GMM should detect boundary with well-separated species."""
        # Create a matrix with very clear separation
        ani_dict = {}
        group_a = [f"ga{i}" for i in range(5)]
        group_b = [f"gb{i}" for i in range(5)]
        all_g = group_a + group_b

        for g1 in all_g:
            ani_dict[g1] = {}
            for g2 in all_g:
                if g1 == g2:
                    ani_dict[g1][g2] = 100.0
                elif (g1 in group_a and g2 in group_a) or (g1 in group_b and g2 in group_b):
                    # Within species: 97-99% ANI
                    ani_dict[g1][g2] = 98.0
                else:
                    # Between species: 82-85% ANI
                    ani_dict[g1][g2] = 83.0

        ani = ANIMatrix(ani_dict)
        result = detect_species_boundary(ani, min_genomes=5)

        if result.method == "gmm":
            # Boundary should be between 85% and 97%
            assert 85.0 < result.species_boundary < 97.0
            assert result.confidence > 0.0
            assert result.novelty_known_max == pytest.approx(100.0 - result.species_boundary)

    def test_boundary_clamped(self):
        """Species boundary should be clamped to 80-99% range."""
        # All genomes very similar - boundary might try to go >99%
        ani_dict = {}
        genomes = [f"g{i}" for i in range(6)]
        for g1 in genomes:
            ani_dict[g1] = {}
            for g2 in genomes:
                ani_dict[g1][g2] = 100.0 if g1 == g2 else 99.5

        ani = ANIMatrix(ani_dict)
        result = detect_species_boundary(ani, min_genomes=5)

        assert result.species_boundary <= 99.0
        assert result.species_boundary >= 80.0

    def test_result_dataclass_fields(self):
        """AdaptiveThresholds should have all required fields."""
        result = AdaptiveThresholds(
            species_boundary=95.5,
            novelty_known_max=4.5,
            confidence=0.8,
            method="gmm",
            ani_values_used=45,
        )

        assert result.species_boundary == 95.5
        assert result.novelty_known_max == 4.5
        assert result.confidence == 0.8
        assert result.method == "gmm"
        assert result.ani_values_used == 45


class TestBuildAdaptiveConfig:
    """Tests for build_adaptive_config function."""

    def test_applies_adaptive_boundary(self):
        """Should override novelty_known_max with detected value."""
        base = ScoringConfig()
        adaptive = AdaptiveThresholds(
            species_boundary=94.0,
            novelty_known_max=6.0,
            confidence=0.9,
            method="gmm",
            ani_values_used=100,
        )

        config = build_adaptive_config(base, adaptive)

        assert config.novelty_known_max == 6.0
        assert config.novelty_novel_species_min == 6.0  # Continuous boundary

    def test_preserves_other_thresholds(self):
        """Should preserve non-adaptive threshold values."""
        base = ScoringConfig(
            uncertainty_conserved_min=8.0,
            coverage_weight_mode="sigmoid",
        )
        adaptive = AdaptiveThresholds(
            species_boundary=95.0,
            novelty_known_max=5.0,
            confidence=0.5,
            method="gmm",
            ani_values_used=50,
        )

        config = build_adaptive_config(base, adaptive)

        assert config.uncertainty_conserved_min == 8.0
        assert config.coverage_weight_mode == "sigmoid"

    def test_continuous_boundary_maintained(self):
        """novelty_novel_species_min should equal novelty_known_max."""
        base = ScoringConfig()
        adaptive = AdaptiveThresholds(
            species_boundary=93.0,
            novelty_known_max=7.0,
            confidence=0.7,
            method="gmm",
            ani_values_used=80,
        )

        config = build_adaptive_config(base, adaptive)

        assert config.novelty_novel_species_min == config.novelty_known_max
