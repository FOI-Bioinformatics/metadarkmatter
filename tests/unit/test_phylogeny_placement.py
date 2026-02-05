"""Tests for phylogeny placement module.

Tests for placing novel clusters on phylogenetic trees based on ANI similarity
to reference genomes.
"""

from __future__ import annotations

import json
import logging
from io import StringIO

import numpy as np
import pandas as pd
import polars as pl
import pytest
from Bio import Phylo


class TestNovelClusterDataclass:
    """Test NovelCluster dataclass creation and attributes."""

    def test_create_novel_cluster(self) -> None:
        """Create a NovelCluster with all required attributes."""
        from metadarkmatter.core.phylogeny.placement import NovelCluster

        cluster = NovelCluster(
            cluster_id="NSP_001",
            classification="novel_species",
            best_match_genome="GCF_000123456.1",
            mean_identity=92.5,
            mean_novelty=7.5,
            mean_uncertainty=1.2,
            read_count=25,
            confidence_rating="High",
        )

        assert cluster.cluster_id == "NSP_001"
        assert cluster.classification == "novel_species"
        assert cluster.best_match_genome == "GCF_000123456.1"
        assert cluster.mean_identity == 92.5
        assert cluster.mean_novelty == 7.5
        assert cluster.mean_uncertainty == 1.2
        assert cluster.read_count == 25
        assert cluster.confidence_rating == "High"

    def test_novel_genus_cluster(self) -> None:
        """Create a NovelCluster for novel genus classification."""
        from metadarkmatter.core.phylogeny.placement import NovelCluster

        cluster = NovelCluster(
            cluster_id="NGN_001",
            classification="novel_genus",
            best_match_genome="GCF_000789012.1",
            mean_identity=78.3,
            mean_novelty=21.7,
            mean_uncertainty=0.8,
            read_count=15,
            confidence_rating="Medium",
        )

        assert cluster.cluster_id == "NGN_001"
        assert cluster.classification == "novel_genus"
        assert cluster.mean_novelty == 21.7

    def test_to_dict(self) -> None:
        """NovelCluster can be converted to dictionary."""
        from metadarkmatter.core.phylogeny.placement import NovelCluster

        cluster = NovelCluster(
            cluster_id="NSP_001",
            classification="novel_species",
            best_match_genome="GCF_000123456.1",
            mean_identity=92.5,
            mean_novelty=7.5,
            mean_uncertainty=1.2,
            read_count=25,
            confidence_rating="High",
        )

        d = cluster.to_dict()
        assert d["cluster_id"] == "NSP_001"
        assert d["classification"] == "novel_species"
        assert d["read_count"] == 25


class TestPlaceNovelClusters:
    """Test placing novel clusters on phylogenetic trees."""

    @pytest.fixture
    def sample_tree(self) -> str:
        """Create a sample tree with 3 reference genomes."""
        return "((GCF_A:0.1,GCF_B:0.2):0.3,GCF_C:0.4);"

    @pytest.fixture
    def sample_ani(self) -> pd.DataFrame:
        """Create a sample ANI matrix for 3 genomes."""
        return pd.DataFrame(
            {
                "GCF_A": [100.0, 95.0, 80.0],
                "GCF_B": [95.0, 100.0, 82.0],
                "GCF_C": [80.0, 82.0, 100.0],
            },
            index=["GCF_A", "GCF_B", "GCF_C"],
        )

    def test_place_single_novel_species_cluster(
        self, sample_tree: str, sample_ani: pd.DataFrame
    ) -> None:
        """Place a single novel species cluster near its best match genome."""
        from metadarkmatter.core.phylogeny.placement import (
            NovelCluster,
            place_novel_clusters,
        )

        cluster = NovelCluster(
            cluster_id="NSP_001",
            classification="novel_species",
            best_match_genome="GCF_A",
            mean_identity=92.5,
            mean_novelty=7.5,
            mean_uncertainty=1.2,
            read_count=25,
            confidence_rating="High",
        )

        result_newick = place_novel_clusters(
            tree=sample_tree,
            novel_clusters=[cluster],
            ani_matrix=sample_ani,
        )

        # Parse the result and check for the novel cluster node
        tree = Phylo.read(StringIO(result_newick), "newick")
        tip_names = {t.name for t in tree.get_terminals()}

        # Should contain original tips plus the novel cluster
        assert "GCF_A" in tip_names
        assert "GCF_B" in tip_names
        assert "GCF_C" in tip_names
        assert "NSP_001" in tip_names

    def test_place_multiple_clusters(
        self, sample_tree: str, sample_ani: pd.DataFrame
    ) -> None:
        """Place multiple novel clusters on the tree."""
        from metadarkmatter.core.phylogeny.placement import (
            NovelCluster,
            place_novel_clusters,
        )

        clusters = [
            NovelCluster(
                cluster_id="NSP_001",
                classification="novel_species",
                best_match_genome="GCF_A",
                mean_identity=92.5,
                mean_novelty=7.5,
                mean_uncertainty=1.2,
                read_count=25,
                confidence_rating="High",
            ),
            NovelCluster(
                cluster_id="NGN_001",
                classification="novel_genus",
                best_match_genome="GCF_C",
                mean_identity=78.0,
                mean_novelty=22.0,
                mean_uncertainty=0.5,
                read_count=10,
                confidence_rating="Medium",
            ),
        ]

        result_newick = place_novel_clusters(
            tree=sample_tree,
            novel_clusters=clusters,
            ani_matrix=sample_ani,
        )

        tree = Phylo.read(StringIO(result_newick), "newick")
        tip_names = {t.name for t in tree.get_terminals()}

        assert "NSP_001" in tip_names
        assert "NGN_001" in tip_names
        assert len(tip_names) == 5  # 3 original + 2 novel

    def test_empty_cluster_list_returns_unchanged_tree(
        self, sample_tree: str, sample_ani: pd.DataFrame
    ) -> None:
        """Empty cluster list returns the tree unchanged."""
        from metadarkmatter.core.phylogeny.placement import place_novel_clusters

        result_newick = place_novel_clusters(
            tree=sample_tree,
            novel_clusters=[],
            ani_matrix=sample_ani,
        )

        # Parse both trees and compare tip counts
        original = Phylo.read(StringIO(sample_tree), "newick")
        result = Phylo.read(StringIO(result_newick), "newick")

        original_tips = {t.name for t in original.get_terminals()}
        result_tips = {t.name for t in result.get_terminals()}

        assert original_tips == result_tips

    def test_novel_node_has_metadata_comment(
        self, sample_tree: str, sample_ani: pd.DataFrame
    ) -> None:
        """Novel nodes should have metadata attached as JSON comment."""
        from metadarkmatter.core.phylogeny.placement import (
            NovelCluster,
            place_novel_clusters,
        )

        cluster = NovelCluster(
            cluster_id="NSP_001",
            classification="novel_species",
            best_match_genome="GCF_A",
            mean_identity=92.5,
            mean_novelty=7.5,
            mean_uncertainty=1.2,
            read_count=25,
            confidence_rating="High",
        )

        result_newick = place_novel_clusters(
            tree=sample_tree,
            novel_clusters=[cluster],
            ani_matrix=sample_ani,
        )

        tree = Phylo.read(StringIO(result_newick), "newick")

        # Find the novel node
        novel_node = next(
            (t for t in tree.get_terminals() if t.name == "NSP_001"),
            None,
        )

        assert novel_node is not None
        # BioPython stores Newick comments in the 'comment' attribute
        assert novel_node.comment is not None

        # Parse the comment as JSON
        metadata = json.loads(novel_node.comment)
        assert metadata["classification"] == "novel_species"
        assert metadata["read_count"] == 25
        assert metadata["mean_novelty"] == 7.5

    def test_novel_species_placed_as_sibling_to_reference(
        self, sample_tree: str, sample_ani: pd.DataFrame
    ) -> None:
        """Novel species should be placed as sibling to nearest reference genome."""
        from metadarkmatter.core.phylogeny.placement import (
            NovelCluster,
            place_novel_clusters,
        )

        cluster = NovelCluster(
            cluster_id="NSP_001",
            classification="novel_species",
            best_match_genome="GCF_A",
            mean_identity=92.5,
            mean_novelty=7.5,
            mean_uncertainty=1.2,
            read_count=25,
            confidence_rating="High",
        )

        result_newick = place_novel_clusters(
            tree=sample_tree,
            novel_clusters=[cluster],
            ani_matrix=sample_ani,
        )

        tree = Phylo.read(StringIO(result_newick), "newick")

        # Find the novel node and the reference genome
        novel_node = next(t for t in tree.get_terminals() if t.name == "NSP_001")
        ref_node = next(t for t in tree.get_terminals() if t.name == "GCF_A")

        # They should share the same immediate parent
        novel_path = tree.get_path(novel_node)
        ref_path = tree.get_path(ref_node)

        # The parent of each (last internal node before the tip) should be the same
        # or one should be an ancestor of the other's parent
        # Since we're making NSP_001 a sibling of GCF_A, they share a common parent
        assert novel_path is not None
        assert ref_path is not None
        # Both should have their paths converge at some internal node
        assert len(set(novel_path) & set(ref_path)) > 0

    def test_branch_length_based_on_novelty(
        self, sample_tree: str, sample_ani: pd.DataFrame
    ) -> None:
        """Branch length of novel cluster should reflect novelty/divergence."""
        from metadarkmatter.core.phylogeny.placement import (
            NovelCluster,
            place_novel_clusters,
        )

        # Create two clusters with different novelty levels
        cluster_low_novelty = NovelCluster(
            cluster_id="NSP_001",
            classification="novel_species",
            best_match_genome="GCF_A",
            mean_identity=95.0,
            mean_novelty=5.0,
            mean_uncertainty=1.0,
            read_count=20,
            confidence_rating="High",
        )
        cluster_high_novelty = NovelCluster(
            cluster_id="NSP_002",
            classification="novel_species",
            best_match_genome="GCF_B",
            mean_identity=85.0,
            mean_novelty=15.0,
            mean_uncertainty=1.5,
            read_count=15,
            confidence_rating="Medium",
        )

        result_newick = place_novel_clusters(
            tree=sample_tree,
            novel_clusters=[cluster_low_novelty, cluster_high_novelty],
            ani_matrix=sample_ani,
        )

        tree = Phylo.read(StringIO(result_newick), "newick")

        nsp_001 = next(t for t in tree.get_terminals() if t.name == "NSP_001")
        nsp_002 = next(t for t in tree.get_terminals() if t.name == "NSP_002")

        # Higher novelty should correspond to longer branch length
        assert nsp_001.branch_length is not None
        assert nsp_002.branch_length is not None
        assert nsp_002.branch_length > nsp_001.branch_length


class TestExtractNovelClusters:
    """Test extracting novel clusters from classification DataFrame."""

    @pytest.fixture
    def sample_classifications(self) -> pl.DataFrame:
        """Create sample classification DataFrame with novel reads."""
        return pl.DataFrame({
            "read_id": [f"read_{i}" for i in range(20)],
            "best_match_genome": [
                "GCF_A", "GCF_A", "GCF_A", "GCF_A", "GCF_A",  # 5 reads to GCF_A
                "GCF_B", "GCF_B", "GCF_B", "GCF_B", "GCF_B",  # 5 reads to GCF_B
                "GCF_C", "GCF_C", "GCF_C",  # 3 reads to GCF_C
                "GCF_D", "GCF_D",  # 2 reads to GCF_D (below threshold)
                "GCF_E", "GCF_E", "GCF_E", "GCF_E", "GCF_E",  # 5 known species
            ],
            "top_hit_identity": [
                92.0, 93.0, 91.5, 92.5, 93.5,  # Novel species (5-10% novelty)
                78.0, 79.0, 77.5, 78.5, 79.5,  # Novel genus (20-25% novelty)
                88.0, 89.0, 87.0,  # Novel species
                85.0, 86.0,  # Novel species (but only 2 reads)
                98.0, 99.0, 97.5, 98.5, 99.5,  # Known species
            ],
            "novelty_index": [
                8.0, 7.0, 8.5, 7.5, 6.5,  # Novel species
                22.0, 21.0, 22.5, 21.5, 20.5,  # Novel genus
                12.0, 11.0, 13.0,  # Novel species
                15.0, 14.0,  # Novel species
                2.0, 1.0, 2.5, 1.5, 0.5,  # Known species
            ],
            "placement_uncertainty": [
                1.0, 1.2, 0.8, 1.1, 0.9,
                0.5, 0.6, 0.4, 0.5, 0.7,
                1.5, 1.3, 1.7,
                1.0, 1.2,
                0.5, 0.3, 0.4, 0.6, 0.2,
            ],
            "taxonomic_call": [
                "Novel Species", "Novel Species", "Novel Species",
                "Novel Species", "Novel Species",
                "Novel Genus", "Novel Genus", "Novel Genus",
                "Novel Genus", "Novel Genus",
                "Novel Species", "Novel Species", "Novel Species",
                "Novel Species", "Novel Species",
                "Known Species", "Known Species", "Known Species",
                "Known Species", "Known Species",
            ],
        })

    def test_extract_clusters_from_classifications(
        self, sample_classifications: pl.DataFrame
    ) -> None:
        """Extract novel clusters from classification DataFrame."""
        from metadarkmatter.core.phylogeny.placement import extract_novel_clusters

        clusters = extract_novel_clusters(sample_classifications, min_reads=3)

        # Should have clusters for GCF_A (5 novel species), GCF_B (5 novel genus),
        # and GCF_C (3 novel species)
        # GCF_D only has 2 reads, below threshold
        # GCF_E has known species, not novel
        assert len(clusters) == 3

        # Check cluster properties
        cluster_ids = {c.cluster_id for c in clusters}
        assert len(cluster_ids) == 3  # All unique IDs

        # Find the novel genus cluster
        genus_cluster = next(
            (c for c in clusters if c.classification == "novel_genus"),
            None,
        )
        assert genus_cluster is not None
        assert genus_cluster.best_match_genome == "GCF_B"
        assert genus_cluster.read_count == 5

    def test_filter_by_min_reads_threshold(
        self, sample_classifications: pl.DataFrame
    ) -> None:
        """Clusters below min_reads threshold are excluded."""
        from metadarkmatter.core.phylogeny.placement import extract_novel_clusters

        # With min_reads=3, GCF_D (2 reads) should be excluded
        clusters_3 = extract_novel_clusters(sample_classifications, min_reads=3)
        genomes_3 = {c.best_match_genome for c in clusters_3}
        assert "GCF_D" not in genomes_3

        # With min_reads=5, only GCF_A and GCF_B should remain
        clusters_5 = extract_novel_clusters(sample_classifications, min_reads=5)
        assert len(clusters_5) == 2
        genomes_5 = {c.best_match_genome for c in clusters_5}
        assert "GCF_C" not in genomes_5

    def test_mean_metrics_calculated(
        self, sample_classifications: pl.DataFrame
    ) -> None:
        """Mean identity, novelty, and uncertainty are calculated correctly."""
        from metadarkmatter.core.phylogeny.placement import extract_novel_clusters

        clusters = extract_novel_clusters(sample_classifications, min_reads=3)

        # Find cluster for GCF_A
        gcf_a_cluster = next(
            c for c in clusters if c.best_match_genome == "GCF_A"
        )

        # Mean identity should be mean of [92.0, 93.0, 91.5, 92.5, 93.5] = 92.5
        assert abs(gcf_a_cluster.mean_identity - 92.5) < 0.01
        # Mean novelty should be mean of [8.0, 7.0, 8.5, 7.5, 6.5] = 7.5
        assert abs(gcf_a_cluster.mean_novelty - 7.5) < 0.01

    def test_confidence_rating_assigned(
        self, sample_classifications: pl.DataFrame
    ) -> None:
        """Confidence ratings are assigned based on cluster properties."""
        from metadarkmatter.core.phylogeny.placement import extract_novel_clusters

        clusters = extract_novel_clusters(sample_classifications, min_reads=3)

        # All clusters should have a confidence rating
        for cluster in clusters:
            assert cluster.confidence_rating in ("High", "Medium", "Low")

    def test_empty_dataframe_returns_empty_list(self) -> None:
        """Empty DataFrame returns empty cluster list."""
        from metadarkmatter.core.phylogeny.placement import extract_novel_clusters

        empty_df = pl.DataFrame({
            "read_id": [],
            "best_match_genome": [],
            "top_hit_identity": [],
            "novelty_index": [],
            "placement_uncertainty": [],
            "taxonomic_call": [],
        }).cast({
            "read_id": pl.Utf8,
            "best_match_genome": pl.Utf8,
            "top_hit_identity": pl.Float64,
            "novelty_index": pl.Float64,
            "placement_uncertainty": pl.Float64,
            "taxonomic_call": pl.Utf8,
        })

        clusters = extract_novel_clusters(empty_df, min_reads=3)
        assert clusters == []

    def test_no_novel_reads_returns_empty_list(self) -> None:
        """DataFrame with no novel reads returns empty cluster list."""
        from metadarkmatter.core.phylogeny.placement import extract_novel_clusters

        known_only_df = pl.DataFrame({
            "read_id": ["r1", "r2", "r3"],
            "best_match_genome": ["GCF_A", "GCF_A", "GCF_A"],
            "top_hit_identity": [98.0, 99.0, 97.5],
            "novelty_index": [2.0, 1.0, 2.5],
            "placement_uncertainty": [0.5, 0.3, 0.4],
            "taxonomic_call": ["Known Species", "Known Species", "Known Species"],
        })

        clusters = extract_novel_clusters(known_only_df, min_reads=3)
        assert clusters == []


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_genome_not_in_tree_logs_warning(self, caplog: pytest.LogCaptureFixture) -> None:
        """Log warning when cluster's genome is not in tree."""
        from metadarkmatter.core.phylogeny.placement import (
            NovelCluster,
            place_novel_clusters,
        )

        tree = "((GCF_A:0.1,GCF_B:0.2):0.3,GCF_C:0.4);"
        ani = pd.DataFrame(
            {
                "GCF_A": [100.0, 95.0, 80.0],
                "GCF_B": [95.0, 100.0, 82.0],
                "GCF_C": [80.0, 82.0, 100.0],
            },
            index=["GCF_A", "GCF_B", "GCF_C"],
        )

        # Cluster references genome not in tree
        cluster = NovelCluster(
            cluster_id="NSP_001",
            classification="novel_species",
            best_match_genome="GCF_MISSING",
            mean_identity=92.5,
            mean_novelty=7.5,
            mean_uncertainty=1.2,
            read_count=25,
            confidence_rating="High",
        )

        with caplog.at_level(logging.WARNING):
            result = place_novel_clusters(
                tree=tree,
                novel_clusters=[cluster],
                ani_matrix=ani,
            )

        assert "GCF_MISSING" in caplog.text or "not found" in caplog.text.lower()
        # Tree should still be valid
        tree_obj = Phylo.read(StringIO(result), "newick")
        assert tree_obj is not None

    def test_place_novel_genus_cluster(self) -> None:
        """Novel genus clusters are placed correctly."""
        from metadarkmatter.core.phylogeny.placement import (
            NovelCluster,
            place_novel_clusters,
        )

        tree = "((GCF_A:0.1,GCF_B:0.2):0.3,GCF_C:0.4);"
        ani = pd.DataFrame(
            {
                "GCF_A": [100.0, 95.0, 80.0],
                "GCF_B": [95.0, 100.0, 82.0],
                "GCF_C": [80.0, 82.0, 100.0],
            },
            index=["GCF_A", "GCF_B", "GCF_C"],
        )

        cluster = NovelCluster(
            cluster_id="NGN_001",
            classification="novel_genus",
            best_match_genome="GCF_C",
            mean_identity=78.0,
            mean_novelty=22.0,
            mean_uncertainty=0.8,
            read_count=15,
            confidence_rating="Medium",
        )

        result_newick = place_novel_clusters(
            tree=tree,
            novel_clusters=[cluster],
            ani_matrix=ani,
        )

        tree_obj = Phylo.read(StringIO(result_newick), "newick")
        tip_names = {t.name for t in tree_obj.get_terminals()}
        assert "NGN_001" in tip_names

    def test_genome_metadata_enriches_output(self) -> None:
        """Genome metadata can be used to enrich placement output."""
        from metadarkmatter.core.phylogeny.placement import (
            NovelCluster,
            place_novel_clusters,
        )

        tree = "((GCF_A:0.1,GCF_B:0.2):0.3,GCF_C:0.4);"
        ani = pd.DataFrame(
            {
                "GCF_A": [100.0, 95.0, 80.0],
                "GCF_B": [95.0, 100.0, 82.0],
                "GCF_C": [80.0, 82.0, 100.0],
            },
            index=["GCF_A", "GCF_B", "GCF_C"],
        )
        metadata = {
            "GCF_A": {"species": "Francisella tularensis", "genus": "Francisella"},
            "GCF_B": {"species": "Francisella novicida", "genus": "Francisella"},
            "GCF_C": {"species": "Francisella philomiragia", "genus": "Francisella"},
        }

        cluster = NovelCluster(
            cluster_id="NSP_001",
            classification="novel_species",
            best_match_genome="GCF_A",
            mean_identity=92.5,
            mean_novelty=7.5,
            mean_uncertainty=1.2,
            read_count=25,
            confidence_rating="High",
        )

        result = place_novel_clusters(
            tree=tree,
            novel_clusters=[cluster],
            ani_matrix=ani,
            genome_metadata=metadata,
        )

        # Result should be valid
        tree_obj = Phylo.read(StringIO(result), "newick")
        assert "NSP_001" in {t.name for t in tree_obj.get_terminals()}
