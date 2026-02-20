"""
Tests for the extended ANI matrix builder module.

Tests ANI estimation functions, matrix construction with novel clusters,
and reference selection algorithms.
"""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from metadarkmatter.core.novel_diversity import NovelCluster
from metadarkmatter.visualization.report.components.extended_matrix_builder import (
    build_extended_similarity_matrix,
    estimate_novel_to_novel_similarity,
    estimate_novel_to_reference_similarity,
    select_relevant_references,
)


@pytest.fixture
def sample_ani_matrix() -> pl.DataFrame:
    """Create a sample ANI matrix for testing."""
    genomes = ["GCF_001", "GCF_002", "GCF_003", "GCF_004"]
    # Create symmetric ANI matrix
    data = {
        "genome": genomes,
        "GCF_001": [100.0, 95.0, 85.0, 75.0],
        "GCF_002": [95.0, 100.0, 88.0, 76.0],
        "GCF_003": [85.0, 88.0, 100.0, 92.0],
        "GCF_004": [75.0, 76.0, 92.0, 100.0],
    }
    return pl.DataFrame(data)


@pytest.fixture
def sample_clusters() -> list[NovelCluster]:
    """Create sample novel clusters for testing."""
    cluster1 = NovelCluster(
        cluster_id="NSP_001",
        taxonomic_call="Novel Species",
        nearest_genome="GCF_001",
        nearest_species="Francisella tularensis",
        nearest_genus="Francisella",
        nearest_family="Francisellaceae",
        novelty_band=5,
        read_count=25,
        mean_novelty_index=7.5,
        novelty_min=5.0,
        novelty_max=10.0,
        mean_placement_uncertainty=1.5,
        mean_discovery_score=82.0,
        suggested_name="Francisella sp. MDM-001",
        confidence="High",
        phylogenetic_context="Within genus Francisella, novel species",
        read_ids=["read_1", "read_2"],
    )

    cluster2 = NovelCluster(
        cluster_id="NGN_001",
        taxonomic_call="Novel Genus",
        nearest_genome="GCF_003",
        nearest_species="Francisella philomiragia",
        nearest_genus="Francisella",
        nearest_family="Francisellaceae",
        novelty_band=20,
        read_count=12,
        mean_novelty_index=22.0,
        novelty_min=20.0,
        novelty_max=25.0,
        mean_placement_uncertainty=3.0,
        mean_discovery_score=65.0,
        suggested_name="Francisellaceae gen. nov. MDM-001",
        confidence="Medium",
        phylogenetic_context="Within family Francisellaceae, novel genus",
        read_ids=["read_3", "read_4"],
    )

    return [cluster1, cluster2]


@pytest.fixture
def genome_labels_map() -> dict[str, str]:
    """Create sample genome labels map."""
    return {
        "GCF_001": "GCF_001 (F. tularensis)",
        "GCF_002": "GCF_002 (F. novicida)",
        "GCF_003": "GCF_003 (F. philomiragia)",
        "GCF_004": "GCF_004 (F. orientalis)",
    }


class TestEstimateNovelToReferenceSimilarity:
    """Tests for estimate_novel_to_reference_similarity function."""

    def test_direct_reference_returns_estimated_ani(
        self, sample_clusters: list[NovelCluster]
    ):
        """When reference is the nearest genome, return cluster's estimated ANI."""
        cluster = sample_clusters[0]
        ani_dict: dict[tuple[str, str], float] = {}

        result = estimate_novel_to_reference_similarity(
            cluster, "GCF_001", ani_dict, default_value=70.0
        )

        # estimated_ani = 100 - mean_novelty_index = 100 - 7.5 = 92.5
        assert result == cluster.estimated_ani
        assert abs(result - 92.5) < 0.1

    def test_other_reference_uses_triangular_estimation(
        self, sample_clusters: list[NovelCluster]
    ):
        """For other references, use triangular estimation via nearest."""
        cluster = sample_clusters[0]  # nearest is GCF_001
        ani_dict = {
            ("GCF_001", "GCF_002"): 95.0,
            ("GCF_002", "GCF_001"): 95.0,
        }

        result = estimate_novel_to_reference_similarity(
            cluster, "GCF_002", ani_dict, default_value=70.0
        )

        # Should be less than the direct ANI due to additional divergence
        assert result < 95.0
        assert result > 70.0  # Should be above default

    def test_distant_reference_clamps_to_default(
        self, sample_clusters: list[NovelCluster]
    ):
        """Very distant references should return values clamped to default."""
        cluster = sample_clusters[1]  # Novel genus with high novelty
        ani_dict = {
            ("GCF_003", "GCF_004"): 72.0,  # Distant reference
            ("GCF_004", "GCF_003"): 72.0,
        }

        result = estimate_novel_to_reference_similarity(
            cluster, "GCF_004", ani_dict, default_value=70.0
        )

        # Should be at or near default for distant combinations
        assert result >= 70.0
        assert result <= 100.0


class TestEstimateNovelToNovelSimilarity:
    """Tests for estimate_novel_to_novel_similarity function."""

    def test_same_nearest_reference(self):
        """Clusters with same nearest reference should have higher similarity."""
        cluster1 = NovelCluster(
            cluster_id="NSP_001",
            taxonomic_call="Novel Species",
            nearest_genome="GCF_001",
            nearest_species="Species A",
            nearest_genus="Genus A",
            nearest_family="Family A",
            novelty_band=5,
            read_count=10,
            mean_novelty_index=7.0,
            novelty_min=5.0,
            novelty_max=10.0,
            mean_placement_uncertainty=1.0,
            suggested_name="Sp. nov. 1",
            confidence="High",
            phylogenetic_context="Context 1",
        )

        cluster2 = NovelCluster(
            cluster_id="NSP_002",
            taxonomic_call="Novel Species",
            nearest_genome="GCF_001",  # Same nearest
            nearest_species="Species A",
            nearest_genus="Genus A",
            nearest_family="Family A",
            novelty_band=10,
            read_count=8,
            mean_novelty_index=12.0,
            novelty_min=10.0,
            novelty_max=15.0,
            mean_placement_uncertainty=2.0,
            suggested_name="Sp. nov. 2",
            confidence="Medium",
            phylogenetic_context="Context 2",
        )

        ani_dict: dict[tuple[str, str], float] = {}

        result = estimate_novel_to_novel_similarity(cluster1, cluster2, ani_dict, 70.0)

        # Average of estimated ANIs minus penalty
        # (93 + 88) / 2 - 2 = 88.5
        assert result > 80.0
        assert result < 95.0

    def test_different_nearest_references(self, sample_clusters: list[NovelCluster]):
        """Clusters with different nearest references estimate through ref-to-ref ANI."""
        cluster1 = sample_clusters[0]  # nearest GCF_001
        cluster2 = sample_clusters[1]  # nearest GCF_003

        ani_dict = {
            ("GCF_001", "GCF_003"): 85.0,
            ("GCF_003", "GCF_001"): 85.0,
        }

        result = estimate_novel_to_novel_similarity(cluster1, cluster2, ani_dict, 70.0)

        # Should be lower than ref-to-ref due to cluster divergences
        assert result < 85.0
        assert result >= 70.0


class TestSelectRelevantReferences:
    """Tests for select_relevant_references function."""

    def test_includes_all_nearest_genomes(
        self,
        sample_clusters: list[NovelCluster],
    ):
        """All nearest genomes from clusters should be included."""
        all_genomes = ["GCF_001", "GCF_002", "GCF_003", "GCF_004", "GCF_005"]
        ani_dict: dict[tuple[str, str], float] = {}

        selected = select_relevant_references(
            ani_dict, sample_clusters, all_genomes, max_refs=10
        )

        # NSP_001 nearest is GCF_001, NGN_001 nearest is GCF_003
        assert "GCF_001" in selected
        assert "GCF_003" in selected

    def test_includes_similar_genomes(self):
        """Should include genomes with high ANI to nearest references."""
        cluster = NovelCluster(
            cluster_id="NSP_001",
            taxonomic_call="Novel Species",
            nearest_genome="GCF_001",
            nearest_species="Species A",
            nearest_genus="Genus A",
            nearest_family="Family A",
            novelty_band=5,
            read_count=10,
            mean_novelty_index=7.0,
            novelty_min=5.0,
            novelty_max=10.0,
            mean_placement_uncertainty=1.0,
            suggested_name="Sp. nov. 1",
            confidence="High",
            phylogenetic_context="Context 1",
        )

        all_genomes = ["GCF_001", "GCF_002", "GCF_003"]
        ani_dict = {
            ("GCF_001", "GCF_002"): 96.0,  # Very similar to nearest
            ("GCF_002", "GCF_001"): 96.0,
            ("GCF_001", "GCF_003"): 78.0,  # Less similar
            ("GCF_003", "GCF_001"): 78.0,
        }

        selected = select_relevant_references(ani_dict, [cluster], all_genomes, max_refs=5)

        assert "GCF_001" in selected
        assert "GCF_002" in selected  # Should be included due to high ANI

    def test_respects_max_refs_limit(self, sample_clusters: list[NovelCluster]):
        """Should not exceed max_refs limit."""
        all_genomes = [f"GCF_{i:03d}" for i in range(100)]
        ani_dict: dict[tuple[str, str], float] = {}

        # Set clusters' nearest genomes to be in the list
        sample_clusters[0] = sample_clusters[0].model_copy(
            update={"nearest_genome": "GCF_001"}
        )
        sample_clusters[1] = sample_clusters[1].model_copy(
            update={"nearest_genome": "GCF_002"}
        )

        selected = select_relevant_references(
            ani_dict, sample_clusters, all_genomes, max_refs=10
        )

        assert len(selected) <= 10

    def test_empty_clusters_returns_subset_of_all_genomes(self):
        """With no clusters, should return first max_refs genomes."""
        all_genomes = ["GCF_001", "GCF_002", "GCF_003", "GCF_004", "GCF_005"]
        ani_dict: dict[tuple[str, str], float] = {}

        selected = select_relevant_references(ani_dict, [], all_genomes, max_refs=3)

        assert len(selected) == 3


class TestBuildExtendedSimilarityMatrix:
    """Tests for build_extended_similarity_matrix function."""

    def test_basic_matrix_construction(
        self,
        sample_ani_matrix: pl.DataFrame,
        sample_clusters: list[NovelCluster],
        genome_labels_map: dict[str, str],
    ):
        """Test basic extended matrix construction."""
        matrix, labels, is_novel = build_extended_similarity_matrix(
            similarity_matrix=sample_ani_matrix,
            novel_clusters=sample_clusters,
            genome_labels_map=genome_labels_map,
            default_value=70.0,
            max_references=10,
            max_clusters=10,
        )

        # Should have references + clusters
        n_refs = 4
        n_clusters = 2
        expected_size = n_refs + n_clusters

        assert matrix.shape == (expected_size, expected_size)
        assert len(labels) == expected_size
        assert len(is_novel) == expected_size

        # First entries should be references (not novel)
        assert not is_novel[0]
        assert not is_novel[1]

        # Last entries should be novel clusters
        assert is_novel[-1]
        assert is_novel[-2]

    def test_diagonal_is_100(
        self,
        sample_ani_matrix: pl.DataFrame,
        sample_clusters: list[NovelCluster],
        genome_labels_map: dict[str, str],
    ):
        """Diagonal should be 100% (self-similarity)."""
        matrix, _, _ = build_extended_similarity_matrix(
            similarity_matrix=sample_ani_matrix,
            novel_clusters=sample_clusters,
            genome_labels_map=genome_labels_map,
        )

        for i in range(matrix.shape[0]):
            assert matrix[i, i] == 100.0

    def test_matrix_is_symmetric(
        self,
        sample_ani_matrix: pl.DataFrame,
        sample_clusters: list[NovelCluster],
        genome_labels_map: dict[str, str],
    ):
        """Matrix should be symmetric."""
        matrix, _, _ = build_extended_similarity_matrix(
            similarity_matrix=sample_ani_matrix,
            novel_clusters=sample_clusters,
            genome_labels_map=genome_labels_map,
        )

        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                assert abs(matrix[i, j] - matrix[j, i]) < 0.01

    def test_novel_cluster_labels_format(
        self,
        sample_ani_matrix: pl.DataFrame,
        sample_clusters: list[NovelCluster],
        genome_labels_map: dict[str, str],
    ):
        """Novel cluster labels should have [*] prefix."""
        _, labels, is_novel = build_extended_similarity_matrix(
            similarity_matrix=sample_ani_matrix,
            novel_clusters=sample_clusters,
            genome_labels_map=genome_labels_map,
        )

        for label, is_novel_flag in zip(labels, is_novel, strict=True):
            if is_novel_flag:
                assert label.startswith("[*]")
            else:
                assert not label.startswith("[*]")

    def test_respects_max_limits(
        self,
        sample_ani_matrix: pl.DataFrame,
        genome_labels_map: dict[str, str],
    ):
        """Should respect max_references and max_clusters limits."""
        # Create many clusters
        clusters = []
        for i in range(30):
            cluster = NovelCluster(
                cluster_id=f"NSP_{i:03d}",
                taxonomic_call="Novel Species",
                nearest_genome="GCF_001",
                nearest_species="Species A",
                nearest_genus="Genus A",
                nearest_family="Family A",
                novelty_band=5,
                read_count=10,
                mean_novelty_index=7.0,
                novelty_min=5.0,
                novelty_max=10.0,
                mean_placement_uncertainty=1.0,
                suggested_name=f"Sp. nov. {i}",
                confidence="High",
                phylogenetic_context="Context",
            )
            clusters.append(cluster)

        matrix, labels, is_novel = build_extended_similarity_matrix(
            similarity_matrix=sample_ani_matrix,
            novel_clusters=clusters,
            genome_labels_map=genome_labels_map,
            max_references=4,
            max_clusters=5,
        )

        n_clusters = sum(1 for flag in is_novel if flag)
        n_refs = sum(1 for flag in is_novel if not flag)

        assert n_clusters <= 5
        assert n_refs <= 4

    def test_empty_clusters_returns_references_only(
        self,
        sample_ani_matrix: pl.DataFrame,
        genome_labels_map: dict[str, str],
    ):
        """With no clusters, should return references only."""
        matrix, labels, is_novel = build_extended_similarity_matrix(
            similarity_matrix=sample_ani_matrix,
            novel_clusters=[],
            genome_labels_map=genome_labels_map,
        )

        assert all(not flag for flag in is_novel)
        assert len(labels) > 0
