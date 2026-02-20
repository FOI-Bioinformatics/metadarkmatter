"""
Unit tests for ANIMatrix class in core.classification.ani_matrix.

Covers initialization edge cases, default ANI for missing values,
hierarchical clustering, genus cluster queries, and union-find
genus counting.
"""

from __future__ import annotations

import numpy as np
import pytest

from metadarkmatter.core.classification.ani_matrix import ANIMatrix


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def three_genome_dict() -> dict[str, dict[str, float]]:
    """Three genomes with known ANI relationships.

    GCF_A and GCF_B are same species (~96% ANI).
    GCF_C is same genus but different species (~82% ANI to both).
    """
    return {
        "GCF_A": {
            "GCF_A": 100.0,
            "GCF_B": 96.0,
            "GCF_C": 82.0,
        },
        "GCF_B": {
            "GCF_A": 96.0,
            "GCF_B": 100.0,
            "GCF_C": 83.0,
        },
        "GCF_C": {
            "GCF_A": 82.0,
            "GCF_B": 83.0,
            "GCF_C": 100.0,
        },
    }


@pytest.fixture
def two_genera_dict() -> dict[str, dict[str, float]]:
    """Five genomes spanning two genera.

    Genus 1: GCF_A, GCF_B, GCF_C (ANI >= 82% within)
    Genus 2: GCF_D, GCF_E (ANI ~90% between them, <75% to genus 1)
    """
    return {
        "GCF_A": {
            "GCF_A": 100.0, "GCF_B": 96.0, "GCF_C": 82.0,
            "GCF_D": 72.0, "GCF_E": 73.0,
        },
        "GCF_B": {
            "GCF_A": 96.0, "GCF_B": 100.0, "GCF_C": 84.0,
            "GCF_D": 71.0, "GCF_E": 74.0,
        },
        "GCF_C": {
            "GCF_A": 82.0, "GCF_B": 84.0, "GCF_C": 100.0,
            "GCF_D": 70.0, "GCF_E": 72.0,
        },
        "GCF_D": {
            "GCF_A": 72.0, "GCF_B": 71.0, "GCF_C": 70.0,
            "GCF_D": 100.0, "GCF_E": 90.0,
        },
        "GCF_E": {
            "GCF_A": 73.0, "GCF_B": 74.0, "GCF_C": 72.0,
            "GCF_D": 90.0, "GCF_E": 100.0,
        },
    }


@pytest.fixture
def sparse_dict_with_missing() -> dict[str, dict[str, float]]:
    """Matrix with some zero (missing) ANI values between existing genomes."""
    return {
        "GCF_X": {
            "GCF_X": 100.0,
            "GCF_Y": 0.0,  # Not computed by skani/fastANI
            "GCF_Z": 88.0,
        },
        "GCF_Y": {
            "GCF_X": 0.0,
            "GCF_Y": 100.0,
            "GCF_Z": 0.0,
        },
        "GCF_Z": {
            "GCF_X": 88.0,
            "GCF_Y": 0.0,
            "GCF_Z": 100.0,
        },
    }


@pytest.fixture
def dict_with_dangling_reference() -> dict[str, dict[str, float]]:
    """Inner dict references a genome not present as outer key.

    This exercises the `if genome2 in self._genome_to_idx` guard
    (line 80) where the referenced genome is absent from the matrix.
    """
    return {
        "GCF_A": {
            "GCF_A": 100.0,
            "GCF_B": 95.0,
            "GCF_MISSING": 78.0,  # Not an outer key
        },
        "GCF_B": {
            "GCF_A": 95.0,
            "GCF_B": 100.0,
        },
    }


# =============================================================================
# Initialization and basic operations
# =============================================================================


class TestANIMatrixInit:
    """Tests for ANIMatrix construction and basic queries."""

    def test_empty_matrix(self):
        """Empty dict should produce an empty matrix."""
        matrix = ANIMatrix({})
        assert len(matrix) == 0
        assert matrix.genomes == set()

    def test_single_genome(self):
        """Single genome matrix should have length 1 and self-ANI 100."""
        matrix = ANIMatrix({"GCF_001": {"GCF_001": 100.0}})
        assert len(matrix) == 1
        assert matrix.get_ani("GCF_001", "GCF_001") == 100.0

    def test_genomes_sorted(self, three_genome_dict):
        """Internal genome ordering should be sorted alphabetically."""
        matrix = ANIMatrix(three_genome_dict)
        assert matrix._genomes == tuple(sorted(three_genome_dict.keys()))

    def test_diagonal_is_100(self, three_genome_dict):
        """Diagonal entries should be 100.0 regardless of input."""
        matrix = ANIMatrix(three_genome_dict)
        for g in three_genome_dict:
            assert matrix.get_ani(g, g) == 100.0

    def test_dangling_reference_skipped(self, dict_with_dangling_reference):
        """Genomes referenced in inner dict but absent from keys should be ignored.

        This covers the `if genome2 in self._genome_to_idx` guard (line 80).
        """
        matrix = ANIMatrix(dict_with_dangling_reference)
        assert len(matrix) == 2
        assert not matrix.has_genome("GCF_MISSING")
        # The valid pair should still be populated
        assert matrix.get_ani("GCF_A", "GCF_B") == pytest.approx(95.0, abs=0.5)


# =============================================================================
# get_ani with default ANI for missing values (line 146)
# =============================================================================


class TestGetANIDefaultValue:
    """Tests that 0.0 entries return default_ani instead of zero."""

    def test_missing_pair_returns_default(self, sparse_dict_with_missing):
        """When stored ANI is 0.0 (uncomputed), default_ani should be returned."""
        matrix = ANIMatrix(sparse_dict_with_missing, default_ani=70.0)

        # GCF_X <-> GCF_Y stored as 0.0
        ani = matrix.get_ani("GCF_X", "GCF_Y")
        assert ani == 70.0

    def test_missing_pair_custom_default(self, sparse_dict_with_missing):
        """Custom default_ani should be used for uncomputed pairs."""
        matrix = ANIMatrix(sparse_dict_with_missing, default_ani=65.0)
        assert matrix.get_ani("GCF_Y", "GCF_Z") == 65.0

    def test_nonzero_pair_not_affected(self, sparse_dict_with_missing):
        """Computed (nonzero) ANI values should be returned directly."""
        matrix = ANIMatrix(sparse_dict_with_missing, default_ani=70.0)
        assert matrix.get_ani("GCF_X", "GCF_Z") == pytest.approx(88.0, abs=0.5)

    def test_unknown_genome_returns_zero_not_default(self, sparse_dict_with_missing):
        """If a genome is not in the matrix at all, return 0.0 (not default)."""
        matrix = ANIMatrix(sparse_dict_with_missing, default_ani=70.0)
        assert matrix.get_ani("GCF_X", "NONEXISTENT") == 0.0
        assert matrix.get_ani("NONEXISTENT", "GCF_X") == 0.0


# =============================================================================
# cluster_genomes_by_ani (lines 214-251)
# =============================================================================


class TestClusterGenomesByANI:
    """Tests for hierarchical clustering of genomes by ANI threshold."""

    def test_single_cluster_at_low_threshold(self, three_genome_dict):
        """All genomes should cluster together at a low genus threshold."""
        matrix = ANIMatrix(three_genome_dict)
        clusters = matrix.cluster_genomes_by_ani(genus_threshold=75.0)

        # All 3 genomes should be in one cluster
        assert len(clusters) == 3
        for genome in three_genome_dict:
            assert genome in clusters
            assert clusters[genome] == set(three_genome_dict.keys())

    def test_two_clusters_at_high_threshold(self, two_genera_dict):
        """At genus_threshold=80, two genera should form separate clusters."""
        matrix = ANIMatrix(two_genera_dict)
        clusters = matrix.cluster_genomes_by_ani(genus_threshold=80.0)

        genus1 = {"GCF_A", "GCF_B", "GCF_C"}
        genus2 = {"GCF_D", "GCF_E"}

        # Each genome maps to its genus
        assert clusters["GCF_A"] == genus1
        assert clusters["GCF_D"] == genus2
        # Cross-genus genomes should not share clusters
        assert clusters["GCF_A"] != clusters["GCF_D"]

    def test_species_level_clustering(self, two_genera_dict):
        """At species threshold (95%), only very close genomes cluster."""
        matrix = ANIMatrix(two_genera_dict)
        clusters = matrix.cluster_genomes_by_ani(
            species_threshold=95.0,
            genus_threshold=95.0,
        )

        # GCF_A and GCF_B are at 96% so they cluster together
        assert "GCF_B" in clusters["GCF_A"]
        assert "GCF_A" in clusters["GCF_B"]
        # GCF_C is 82-84% from A/B, so separate at 95%
        assert "GCF_C" not in clusters["GCF_A"]

    def test_empty_matrix_returns_empty(self):
        """Empty matrix should return empty dict."""
        matrix = ANIMatrix({})
        clusters = matrix.cluster_genomes_by_ani()
        assert clusters == {}

    def test_all_genomes_present_in_result(self, two_genera_dict):
        """Every genome in the matrix should have an entry in the result."""
        matrix = ANIMatrix(two_genera_dict)
        clusters = matrix.cluster_genomes_by_ani(genus_threshold=80.0)
        assert set(clusters.keys()) == set(two_genera_dict.keys())

    def test_cluster_values_are_sets(self, three_genome_dict):
        """Cluster values should be sets of genome identifiers."""
        matrix = ANIMatrix(three_genome_dict)
        clusters = matrix.cluster_genomes_by_ani()
        for genome, members in clusters.items():
            assert isinstance(members, set)
            assert genome in members  # Genome is always in its own cluster


# =============================================================================
# get_genus_cluster_for_genome (lines 270-279)
# =============================================================================


class TestGetGenusClusterForGenome:
    """Tests for single-genome genus cluster query."""

    def test_returns_related_genomes(self, two_genera_dict):
        """Should return all genomes above the threshold for a given genome."""
        matrix = ANIMatrix(two_genera_dict)
        cluster = matrix.get_genus_cluster_for_genome("GCF_A", genus_threshold=80.0)

        # GCF_A, GCF_B (96%), GCF_C (82%) are above 80%
        assert "GCF_A" in cluster
        assert "GCF_B" in cluster
        assert "GCF_C" in cluster
        # GCF_D (72%) and GCF_E (73%) are below 80%
        assert "GCF_D" not in cluster
        assert "GCF_E" not in cluster

    def test_unknown_genome_returns_singleton(self, three_genome_dict):
        """Unknown genome should return a singleton set containing itself."""
        matrix = ANIMatrix(three_genome_dict)
        cluster = matrix.get_genus_cluster_for_genome("UNKNOWN_GENOME")
        assert cluster == {"UNKNOWN_GENOME"}

    def test_high_threshold_isolates_genome(self, two_genera_dict):
        """Very high threshold should return only the genome itself (self=100)."""
        matrix = ANIMatrix(two_genera_dict)
        cluster = matrix.get_genus_cluster_for_genome("GCF_C", genus_threshold=99.0)

        # Only self-ANI (100%) exceeds 99%
        assert cluster == {"GCF_C"}

    def test_low_threshold_returns_all(self, two_genera_dict):
        """Very low threshold should return all genomes."""
        matrix = ANIMatrix(two_genera_dict)
        cluster = matrix.get_genus_cluster_for_genome("GCF_A", genus_threshold=60.0)

        assert cluster == set(two_genera_dict.keys())

    def test_includes_self(self, three_genome_dict):
        """Query genome should always be included in its own cluster."""
        matrix = ANIMatrix(three_genome_dict)
        for genome in three_genome_dict:
            cluster = matrix.get_genus_cluster_for_genome(genome)
            assert genome in cluster


# =============================================================================
# count_distinct_genera (lines 299-343)
# =============================================================================


class TestCountDistinctGenera:
    """Tests for genus counting via union-find."""

    def test_empty_list(self, three_genome_dict):
        """Empty genome list should return 0 genera."""
        matrix = ANIMatrix(three_genome_dict)
        count, groups = matrix.count_distinct_genera([])
        assert count == 0
        assert groups == []

    def test_single_genome(self, three_genome_dict):
        """Single genome should be one genus."""
        matrix = ANIMatrix(three_genome_dict)
        count, groups = matrix.count_distinct_genera(["GCF_A"])
        assert count == 1
        assert len(groups) == 1
        assert "GCF_A" in groups[0]

    def test_two_genera(self, two_genera_dict):
        """Should detect two distinct genera at threshold 80%."""
        matrix = ANIMatrix(two_genera_dict)
        genomes = list(two_genera_dict.keys())
        count, groups = matrix.count_distinct_genera(genomes, genus_threshold=80.0)

        assert count == 2
        genus1_members = {"GCF_A", "GCF_B", "GCF_C"}
        genus2_members = {"GCF_D", "GCF_E"}
        group_sets = [g for g in groups]
        assert genus1_members in group_sets
        assert genus2_members in group_sets

    def test_all_same_genus(self, three_genome_dict):
        """All genomes above threshold should form a single genus."""
        matrix = ANIMatrix(three_genome_dict)
        genomes = list(three_genome_dict.keys())
        count, groups = matrix.count_distinct_genera(genomes, genus_threshold=80.0)

        assert count == 1
        assert len(groups) == 1
        assert groups[0] == set(genomes)

    def test_high_threshold_each_separate(self, three_genome_dict):
        """Above species level, each genome becomes its own genus."""
        matrix = ANIMatrix(three_genome_dict)
        genomes = list(three_genome_dict.keys())
        # GCF_A<->GCF_B=96, GCF_A<->GCF_C=82, GCF_B<->GCF_C=83
        # At threshold 97%, no pairs exceed, so 3 genera
        count, groups = matrix.count_distinct_genera(genomes, genus_threshold=97.0)

        assert count == 3
        assert all(len(g) == 1 for g in groups)

    def test_missing_genomes_treated_as_separate(self, three_genome_dict):
        """Genomes not in the matrix should each count as a separate genus."""
        matrix = ANIMatrix(three_genome_dict)
        genomes = ["GCF_A", "GCF_B", "NOT_IN_MATRIX_1", "NOT_IN_MATRIX_2"]
        count, groups = matrix.count_distinct_genera(genomes, genus_threshold=80.0)

        # GCF_A + GCF_B cluster (96% > 80%), plus 2 singletons
        assert count == 3

        # Verify the two missing genomes are singletons
        singletons = [g for g in groups if len(g) == 1]
        singleton_names = {name for s in singletons for name in s}
        assert "NOT_IN_MATRIX_1" in singleton_names
        assert "NOT_IN_MATRIX_2" in singleton_names

    def test_all_genomes_missing(self, three_genome_dict):
        """When all genomes are absent from the matrix, each is separate."""
        matrix = ANIMatrix(three_genome_dict)
        genomes = ["MISSING_A", "MISSING_B"]
        count, groups = matrix.count_distinct_genera(genomes)

        assert count == 2
        assert all(len(g) == 1 for g in groups)

    def test_transitive_closure(self, two_genera_dict):
        """Union-find should handle transitive relationships.

        GCF_A<->GCF_C=82 and GCF_B<->GCF_C=84 at threshold 80 means
        GCF_A, GCF_B, GCF_C should all be in one genus even if
        GCF_A<->GCF_C is the weakest link.
        """
        matrix = ANIMatrix(two_genera_dict)
        genomes = ["GCF_A", "GCF_B", "GCF_C"]
        count, groups = matrix.count_distinct_genera(genomes, genus_threshold=82.0)

        # At 82%, GCF_A<->GCF_C=82 (passes), GCF_B<->GCF_C=84 (passes),
        # GCF_A<->GCF_B=96 (passes), so all three are one genus
        assert count == 1
        assert groups[0] == {"GCF_A", "GCF_B", "GCF_C"}

    def test_subset_of_matrix_genomes(self, two_genera_dict):
        """Counting genera for a subset should only include queried genomes."""
        matrix = ANIMatrix(two_genera_dict)
        # Query only genus 2 genomes
        genomes = ["GCF_D", "GCF_E"]
        count, groups = matrix.count_distinct_genera(genomes, genus_threshold=80.0)

        assert count == 1
        assert groups[0] == {"GCF_D", "GCF_E"}


# =============================================================================
# from_file class method
# =============================================================================


class TestFromFile:
    """Tests for loading ANIMatrix from file."""

    def test_from_csv_file(self, tmp_path, three_genome_dict):
        """Should load matrix correctly from a CSV file."""
        genomes = sorted(three_genome_dict.keys())
        lines = ["genome," + ",".join(genomes)]
        for g1 in genomes:
            row = [g1] + [str(three_genome_dict[g1].get(g2, 0.0)) for g2 in genomes]
            lines.append(",".join(row))

        csv_path = tmp_path / "ani.csv"
        csv_path.write_text("\n".join(lines))

        matrix = ANIMatrix.from_file(csv_path)
        assert len(matrix) == 3
        assert matrix.get_ani("GCF_A", "GCF_B") == pytest.approx(96.0, abs=0.5)

    def test_from_file_custom_default_ani(self, tmp_path):
        """Custom default_ani should propagate through from_file."""
        lines = [
            "genome,GCF_1,GCF_2",
            "GCF_1,100.0,0.0",
            "GCF_2,0.0,100.0",
        ]
        csv_path = tmp_path / "sparse.csv"
        csv_path.write_text("\n".join(lines))

        matrix = ANIMatrix.from_file(csv_path, default_ani=65.0)
        assert matrix.get_ani("GCF_1", "GCF_2") == 65.0


# =============================================================================
# memory_usage_bytes
# =============================================================================


class TestMemoryUsage:
    """Tests for memory usage estimation."""

    def test_memory_scales_with_size(self, three_genome_dict, two_genera_dict):
        """Larger matrix should report more memory usage."""
        small = ANIMatrix(three_genome_dict)
        large = ANIMatrix(two_genera_dict)
        assert large.memory_usage_bytes() > small.memory_usage_bytes()

    def test_memory_positive(self, three_genome_dict):
        """Memory usage should always be positive for nonempty matrix."""
        matrix = ANIMatrix(three_genome_dict)
        assert matrix.memory_usage_bytes() > 0
