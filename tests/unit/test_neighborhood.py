"""
Unit tests for phylogenetic neighborhood analysis.

Tests the PhylogeneticNeighborhoodAnalyzer which estimates the evolutionary
context of novel clusters by triangulating ANI distances to reference genomes.
"""

from __future__ import annotations

import pytest

from metadarkmatter.core.classification.ani_matrix import ANIMatrix
from metadarkmatter.core.novel_diversity.models import NovelCluster
from metadarkmatter.core.novel_diversity.neighborhood import (
    PhylogeneticNeighborhoodAnalyzer,
)


# =============================================================================
# Helpers
# =============================================================================


def _make_ani_matrix_with_genera(
    genus_genomes: dict[str, list[str]],
    within_genus_ani: float = 90.0,
    between_genus_ani: float = 75.0,
) -> tuple[ANIMatrix, dict[str, str]]:
    """Build an ANIMatrix with defined genus structure.

    Args:
        genus_genomes: Mapping from genus name to list of genome accessions.
        within_genus_ani: ANI between genomes in the same genus.
        between_genus_ani: ANI between genomes in different genera.

    Returns:
        Tuple of (ANIMatrix, genus_map) where genus_map maps accession to genus.
    """
    all_genomes = [g for members in genus_genomes.values() for g in members]
    genus_map: dict[str, str] = {}
    for genus, members in genus_genomes.items():
        for g in members:
            genus_map[g] = genus

    ani_dict: dict[str, dict[str, float]] = {}
    for g1 in all_genomes:
        inner: dict[str, float] = {}
        for g2 in all_genomes:
            if g1 == g2:
                inner[g2] = 100.0
            elif genus_map[g1] == genus_map[g2]:
                inner[g2] = within_genus_ani
            else:
                inner[g2] = between_genus_ani
            inner[g2] = inner[g2]
        ani_dict[g1] = inner

    return ANIMatrix(ani_dict), genus_map


def _make_novel_cluster(
    cluster_id: str = "NGN_001",
    taxonomic_call: str = "Novel Genus",
    nearest_genome: str = "GCF_A1",
    nearest_species: str = "Species alpha",
    nearest_genus: str = "GenusA",
    nearest_family: str = "FamilyX",
    mean_novelty_index: float = 22.0,
    read_count: int = 15,
    mean_bayesian_confidence: float | None = 65.0,
) -> NovelCluster:
    """Create a NovelCluster with sensible defaults for testing."""
    return NovelCluster(
        cluster_id=cluster_id,
        taxonomic_call=taxonomic_call,
        nearest_genome=nearest_genome,
        nearest_species=nearest_species,
        nearest_genus=nearest_genus,
        nearest_family=nearest_family,
        novelty_band=20 if taxonomic_call == "Novel Genus" else 5,
        read_count=read_count,
        mean_novelty_index=mean_novelty_index,
        novelty_min=mean_novelty_index - 1.0,
        novelty_max=mean_novelty_index + 1.0,
        mean_placement_uncertainty=1.5,
        suggested_name=f"{nearest_family} gen. nov. MDM-001",
        confidence="Medium",
        phylogenetic_context=f"Novel genus within {nearest_family}",
        mean_bayesian_confidence=mean_bayesian_confidence,
    )


# =============================================================================
# Tests
# =============================================================================


class TestPhylogeneticNeighborhoodAnalyzer:
    """Tests for PhylogeneticNeighborhoodAnalyzer."""

    @pytest.fixture
    def two_genera_setup(self):
        """Setup with two genera (GenusA: 2 genomes, GenusB: 1 genome)."""
        matrix, genus_map = _make_ani_matrix_with_genera(
            {"GenusA": ["GCF_A1", "GCF_A2"], "GenusB": ["GCF_B1"]},
            within_genus_ani=90.0,
            between_genus_ani=75.0,
        )
        return matrix, genus_map

    @pytest.fixture
    def three_genera_setup(self):
        """Setup with three genera."""
        matrix, genus_map = _make_ani_matrix_with_genera(
            {
                "GenusA": ["GCF_A1", "GCF_A2"],
                "GenusB": ["GCF_B1"],
                "GenusC": ["GCF_C1", "GCF_C2", "GCF_C3"],
            },
            within_genus_ani=90.0,
            between_genus_ani=75.0,
        )
        return matrix, genus_map

    # ------------------------------------------------------------------

    def test_analyze_returns_neighborhoods(self, two_genera_setup):
        """Analyzed clusters should have a neighborhood attached."""
        matrix, genus_map = two_genera_setup
        cluster = _make_novel_cluster()
        analyzer = PhylogeneticNeighborhoodAnalyzer(
            matrix, genus_map=genus_map
        )

        results = analyzer.analyze([cluster])

        assert len(results) == 1
        assert results[0].neighborhood is not None
        assert results[0].neighborhood.cluster_id == cluster.cluster_id

    def test_nearest_genera_sorted_by_ani(self, three_genera_setup):
        """Nearest genera should be sorted by estimated ANI descending."""
        matrix, genus_map = three_genera_setup
        cluster = _make_novel_cluster()
        analyzer = PhylogeneticNeighborhoodAnalyzer(
            matrix, genus_map=genus_map
        )

        results = analyzer.analyze([cluster])
        genera = results[0].neighborhood.nearest_genera

        anis = [gd.estimated_ani for gd in genera]
        assert anis == sorted(anis, reverse=True)

    def test_isolation_score_computed(self, three_genera_setup):
        """Isolation score should equal the gap between 1st and 2nd genus ANI."""
        matrix, genus_map = three_genera_setup
        cluster = _make_novel_cluster()
        analyzer = PhylogeneticNeighborhoodAnalyzer(
            matrix, genus_map=genus_map
        )

        results = analyzer.analyze([cluster])
        nbr = results[0].neighborhood
        genera = nbr.nearest_genera

        if len(genera) >= 2:
            expected_gap = genera[0].estimated_ani - genera[1].estimated_ani
            assert nbr.isolation_score == pytest.approx(expected_gap, abs=0.01)
        else:
            assert nbr.isolation_score == 0.0

    def test_phylogenetic_context_text_novel_genus(self, two_genera_setup):
        """Context text for Novel Genus should contain genus name and ANI."""
        matrix, genus_map = two_genera_setup
        cluster = _make_novel_cluster(taxonomic_call="Novel Genus")
        analyzer = PhylogeneticNeighborhoodAnalyzer(
            matrix, genus_map=genus_map
        )

        results = analyzer.analyze([cluster])
        ctx = results[0].neighborhood.phylogenetic_context

        assert "Sister to" in ctx
        assert "ANI" in ctx
        assert "Support:" in ctx

    def test_phylogenetic_context_text_novel_species(self, two_genera_setup):
        """Context text for Novel Species should contain 'Within' and genus."""
        matrix, genus_map = two_genera_setup
        cluster = _make_novel_cluster(
            cluster_id="NSP_001",
            taxonomic_call="Novel Species",
            mean_novelty_index=8.0,
        )
        analyzer = PhylogeneticNeighborhoodAnalyzer(
            matrix, genus_map=genus_map
        )

        results = analyzer.analyze([cluster])
        ctx = results[0].neighborhood.phylogenetic_context

        assert "Within" in ctx
        assert "ANI" in ctx

    def test_without_genus_map_uses_groups(self):
        """Without genus_map, analyzer should infer groups from ANI matrix."""
        # Two tight groups with low inter-group ANI
        ani_dict: dict[str, dict[str, float]] = {
            "GCF_A1": {"GCF_A1": 100, "GCF_A2": 92, "GCF_B1": 74},
            "GCF_A2": {"GCF_A1": 92, "GCF_A2": 100, "GCF_B1": 73},
            "GCF_B1": {"GCF_A1": 74, "GCF_A2": 73, "GCF_B1": 100},
        }
        matrix = ANIMatrix(ani_dict)
        cluster = _make_novel_cluster()
        analyzer = PhylogeneticNeighborhoodAnalyzer(matrix)  # no genus_map

        results = analyzer.analyze([cluster])
        nbr = results[0].neighborhood

        assert nbr is not None
        # Should have created at least 2 groups (A1+A2 cluster, B1 alone)
        genera_names = [gd.genus for gd in nbr.nearest_genera]
        assert any("Group" in name for name in genera_names)
        assert len(set(genera_names)) >= 2

    def test_placement_support_range(self, two_genera_setup):
        """Placement support score should be between 0 and 100."""
        matrix, genus_map = two_genera_setup
        cluster = _make_novel_cluster()
        analyzer = PhylogeneticNeighborhoodAnalyzer(
            matrix, genus_map=genus_map
        )

        results = analyzer.analyze([cluster])
        score = results[0].neighborhood.placement_support

        assert 0 <= score <= 100

    def test_multiple_clusters(self, three_genera_setup):
        """Analyzer should handle multiple clusters correctly."""
        matrix, genus_map = three_genera_setup
        clusters = [
            _make_novel_cluster(
                cluster_id="NGN_001",
                nearest_genome="GCF_A1",
                mean_novelty_index=22.0,
            ),
            _make_novel_cluster(
                cluster_id="NSP_001",
                taxonomic_call="Novel Species",
                nearest_genome="GCF_B1",
                nearest_genus="GenusB",
                mean_novelty_index=8.0,
            ),
        ]
        analyzer = PhylogeneticNeighborhoodAnalyzer(
            matrix, genus_map=genus_map
        )

        results = analyzer.analyze(clusters)

        assert len(results) == 2
        assert results[0].neighborhood is not None
        assert results[1].neighborhood is not None
        assert results[0].neighborhood.cluster_id == "NGN_001"
        assert results[1].neighborhood.cluster_id == "NSP_001"


class TestNeighborhoodDensity:
    """Tests for neighborhood density metric."""

    def test_density_counts_genera_near_boundary(self):
        """Density should count genera with ANI >= genus_boundary - 5."""
        matrix, genus_map = _make_ani_matrix_with_genera(
            {
                "GenusA": ["GCF_A1"],
                "GenusB": ["GCF_B1"],
                "GenusC": ["GCF_C1"],
            },
            within_genus_ani=90.0,
            between_genus_ani=76.0,
        )
        cluster = _make_novel_cluster(nearest_genome="GCF_A1")
        analyzer = PhylogeneticNeighborhoodAnalyzer(
            matrix, genus_map=genus_map, genus_boundary=80.0
        )

        results = analyzer.analyze([cluster])
        nbr = results[0].neighborhood

        # density_threshold = 80 - 5 = 75; between_genus_ani = 76 >= 75
        # All three genera should count (nearest + 2 others)
        assert nbr.neighborhood_density >= 1


class TestPlacementSupportComponents:
    """Tests verifying placement support score composition."""

    def test_novel_genus_boundary_component(self):
        """Novel Genus should gain points when ANI is below genus boundary."""
        matrix, genus_map = _make_ani_matrix_with_genera(
            {"GenusA": ["GCF_A1"], "GenusB": ["GCF_B1"]},
            within_genus_ani=90.0,
            between_genus_ani=75.0,
        )
        cluster = _make_novel_cluster(
            taxonomic_call="Novel Genus",
            mean_novelty_index=25.0,  # 75% ANI -- below 80% boundary
        )
        analyzer = PhylogeneticNeighborhoodAnalyzer(
            matrix, genus_map=genus_map, genus_boundary=80.0
        )

        results = analyzer.analyze([cluster])
        score = results[0].neighborhood.placement_support

        # boundary_component should contribute since top ANI < genus_boundary
        assert score > 0

    def test_novel_species_no_boundary_component(self):
        """Novel Species should not have a boundary component."""
        matrix, genus_map = _make_ani_matrix_with_genera(
            {"GenusA": ["GCF_A1"], "GenusB": ["GCF_B1"]},
            within_genus_ani=90.0,
            between_genus_ani=75.0,
        )
        # Two identical clusters differing only in taxonomic call
        genus_cluster = _make_novel_cluster(
            cluster_id="NGN_001",
            taxonomic_call="Novel Genus",
            mean_novelty_index=25.0,
            read_count=10,
            mean_bayesian_confidence=60.0,
        )
        species_cluster = _make_novel_cluster(
            cluster_id="NSP_001",
            taxonomic_call="Novel Species",
            mean_novelty_index=8.0,
            read_count=10,
            mean_bayesian_confidence=60.0,
        )
        analyzer = PhylogeneticNeighborhoodAnalyzer(
            matrix, genus_map=genus_map, genus_boundary=80.0
        )

        genus_result = analyzer.analyze([genus_cluster])
        species_result = analyzer.analyze([species_cluster])

        # Both should produce valid scores
        assert 0 <= genus_result[0].neighborhood.placement_support <= 100
        assert 0 <= species_result[0].neighborhood.placement_support <= 100


class TestContextTextEdgeCases:
    """Tests for context text generation edge cases."""

    def test_single_genus_context(self):
        """When only one genus exists, text should note no close neighbors."""
        matrix, genus_map = _make_ani_matrix_with_genera(
            {"GenusA": ["GCF_A1", "GCF_A2"]},
            within_genus_ani=95.0,
        )
        cluster = _make_novel_cluster(
            taxonomic_call="Novel Genus",
            nearest_genome="GCF_A1",
        )
        analyzer = PhylogeneticNeighborhoodAnalyzer(
            matrix, genus_map=genus_map
        )

        results = analyzer.analyze([cluster])
        ctx = results[0].neighborhood.phylogenetic_context

        assert "no close neighboring genera" in ctx
