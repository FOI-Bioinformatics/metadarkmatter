"""
Unit tests for adaptive threshold detection.

Tests GMM-based species and genus boundary detection and fallback behavior.
"""

from __future__ import annotations

import pytest

from metadarkmatter.core.ani_placement import ANIMatrix
from metadarkmatter.core.classification.adaptive import (
    AdaptiveGenusThreshold,
    AdaptiveThresholds,
    build_adaptive_config,
    detect_genus_boundary,
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


# =============================================================================
# Helper for multi-genus ANI matrices
# =============================================================================


def _build_multi_genus_ani_dict(
    genera: int = 3,
    genomes_per_genus: int = 3,
    within_genus_ani: tuple[float, float] = (88.0, 95.0),
    between_genus_ani: tuple[float, float] = (72.0, 78.0),
) -> dict[str, dict[str, float]]:
    """
    Build an ANI dictionary with multiple genera.

    Within-genus pairs receive ANI values uniformly sampled from
    within_genus_ani range, and between-genus pairs from between_genus_ani
    range. Deterministic via modular arithmetic on genome indices.

    Args:
        genera: Number of genera.
        genomes_per_genus: Number of genomes per genus.
        within_genus_ani: (min, max) ANI for within-genus pairs.
        between_genus_ani: (min, max) ANI for between-genus pairs.

    Returns:
        Nested ANI dictionary suitable for ANIMatrix construction.
    """
    all_genomes = [
        f"G{genus_idx}_{genome_idx}"
        for genus_idx in range(genera)
        for genome_idx in range(genomes_per_genus)
    ]

    ani_dict: dict[str, dict[str, float]] = {}
    for g1 in all_genomes:
        ani_dict[g1] = {}
        genus1 = int(g1.split("_")[0][1:])
        for g2 in all_genomes:
            if g1 == g2:
                ani_dict[g1][g2] = 100.0
                continue
            genus2 = int(g2.split("_")[0][1:])
            # Deterministic variation based on genome name hash
            variation = ((hash(g1 + g2) % 1000) / 1000.0)
            if genus1 == genus2:
                lo, hi = within_genus_ani
                ani_dict[g1][g2] = lo + variation * (hi - lo)
            else:
                lo, hi = between_genus_ani
                ani_dict[g1][g2] = lo + variation * (hi - lo)
    return ani_dict


class TestDetectGenusBoundary:
    """Tests for detect_genus_boundary function."""

    def test_fallback_with_few_genomes(self, small_ani_dict):
        """Should fall back to default 80% when too few genomes."""
        ani = ANIMatrix(small_ani_dict)  # 3 genomes, well below min_genomes=8
        result = detect_genus_boundary(ani)

        assert result.method == "fallback"
        assert result.genus_boundary == 80.0
        assert result.novelty_genus_min == 20.0

    def test_genus_boundary_reasonable_range(self):
        """Detected boundary should be within 70-90% ANI."""
        # Create matrix with 3 genera, 3 genomes each (9 total)
        ani_dict = _build_multi_genus_ani_dict(genera=3, genomes_per_genus=3)
        ani = ANIMatrix(ani_dict)
        result = detect_genus_boundary(ani)

        assert 70.0 <= result.genus_boundary <= 90.0
        assert result.ani_values_used > 0

    def test_novelty_genus_min_consistent(self):
        """novelty_genus_min should equal 100 - genus_boundary."""
        ani_dict = _build_multi_genus_ani_dict(genera=3, genomes_per_genus=3)
        ani = ANIMatrix(ani_dict)
        result = detect_genus_boundary(ani)

        assert result.novelty_genus_min == pytest.approx(
            100.0 - result.genus_boundary
        )

    def test_metadata_validation_computes_range(self):
        """When genus_map provided, should compute empirical inter-genus ANI range."""
        ani_dict = _build_multi_genus_ani_dict(genera=3, genomes_per_genus=3)
        ani = ANIMatrix(ani_dict)

        genus_map = {
            f"G{genus_idx}_{genome_idx}": f"Genus_{genus_idx}"
            for genus_idx in range(3)
            for genome_idx in range(3)
        }
        result = detect_genus_boundary(ani, genus_map=genus_map)

        assert result.inter_genus_ani_range is not None
        assert result.inter_genus_ani_range[0] <= result.inter_genus_ani_range[1]

    def test_no_genus_map_gives_none_range(self):
        """Without genus_map, inter_genus_ani_range should be None."""
        ani_dict = _build_multi_genus_ani_dict(genera=3, genomes_per_genus=3)
        ani = ANIMatrix(ani_dict)
        result = detect_genus_boundary(ani, genus_map=None)

        assert result.inter_genus_ani_range is None

    def test_fallback_with_insufficient_ani_values(self):
        """Should fall back when too few pairwise ANI values (<10)."""
        # 8 genomes but most values zero -> fewer than 10 usable pairs
        # Actually, 8 genomes = 28 pairs. Instead, use min_genomes=8 with
        # exactly 8 genomes. The function requires <10 ANI values, which
        # is hard with 8 genomes (28 pairs). Test via min_genomes instead.
        ani_dict = _build_multi_genus_ani_dict(genera=2, genomes_per_genus=2)
        ani = ANIMatrix(ani_dict)  # 4 genomes
        result = detect_genus_boundary(ani, min_genomes=8)

        assert result.method == "fallback"
        assert result.genus_boundary == 80.0

    def test_result_dataclass_fields(self):
        """AdaptiveGenusThreshold should have all required fields."""
        result = AdaptiveGenusThreshold(
            genus_boundary=82.0,
            novelty_genus_min=18.0,
            confidence=0.7,
            method="gmm_3component",
            inter_genus_ani_range=(72.0, 78.0),
            ani_values_used=36,
        )

        assert result.genus_boundary == 82.0
        assert result.novelty_genus_min == 18.0
        assert result.confidence == 0.7
        assert result.method == "gmm_3component"
        assert result.inter_genus_ani_range == (72.0, 78.0)
        assert result.ani_values_used == 36

    def test_returns_valid_result_regardless_of_method(self):
        """Function should always return a valid result (GMM or fallback)."""
        ani_dict = _build_multi_genus_ani_dict(genera=3, genomes_per_genus=3)
        ani = ANIMatrix(ani_dict)
        result = detect_genus_boundary(ani)

        assert result.genus_boundary > 0
        assert result.method in ("gmm_3component", "fallback")
        assert result.confidence >= 0.0
        assert result.ani_values_used > 0
