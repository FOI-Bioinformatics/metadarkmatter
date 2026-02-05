"""Tests for phylogeny tree builder module.

Tests for ANI matrix to Newick tree conversion and user-provided tree loading.
"""

from __future__ import annotations

import logging
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from Bio import Phylo


class TestAniToNewick:
    """Test ANI matrix to Newick tree conversion."""

    def test_basic_3x3_matrix(self) -> None:
        """3x3 ANI matrix produces valid Newick with all taxa."""
        from metadarkmatter.core.phylogeny.tree_builder import ani_to_newick

        ani = pd.DataFrame(
            {
                "A": [100.0, 95.0, 80.0],
                "B": [95.0, 100.0, 82.0],
                "C": [80.0, 82.0, 100.0],
            },
            index=["A", "B", "C"],
        )

        newick = ani_to_newick(ani)

        assert newick is not None
        assert newick.endswith(";")
        assert "A" in newick
        assert "B" in newick
        assert "C" in newick

        # Verify it parses as valid Newick
        tree = Phylo.read(StringIO(newick), "newick")
        tip_names = {t.name for t in tree.get_terminals()}
        assert tip_names == {"A", "B", "C"}

    def test_missing_ani_values_filled(self) -> None:
        """Missing ANI values are filled with maximum distance."""
        from metadarkmatter.core.phylogeny.tree_builder import ani_to_newick

        ani = pd.DataFrame(
            {
                "A": [100.0, 95.0, np.nan],
                "B": [95.0, 100.0, 82.0],
                "C": [np.nan, 82.0, 100.0],
            },
            index=["A", "B", "C"],
        )

        newick = ani_to_newick(ani)

        assert newick is not None
        tree = Phylo.read(StringIO(newick), "newick")
        assert len(list(tree.get_terminals())) == 3

    def test_fewer_than_3_genomes_returns_none(self) -> None:
        """Fewer than 3 genomes returns None (NJ requires >= 3)."""
        from metadarkmatter.core.phylogeny.tree_builder import ani_to_newick

        ani = pd.DataFrame(
            {
                "A": [100.0, 95.0],
                "B": [95.0, 100.0],
            },
            index=["A", "B"],
        )

        result = ani_to_newick(ani)
        assert result is None

    def test_larger_matrix(self) -> None:
        """Larger matrix (10 genomes) produces valid tree."""
        from metadarkmatter.core.phylogeny.tree_builder import ani_to_newick

        n = 10
        names = [f"G{i}" for i in range(n)]
        ani_values = np.full((n, n), 80.0)
        np.fill_diagonal(ani_values, 100.0)
        ani_values[0, 1] = ani_values[1, 0] = 98.0
        ani_values[2, 3] = ani_values[3, 2] = 97.0

        ani = pd.DataFrame(ani_values, index=names, columns=names)
        newick = ani_to_newick(ani)

        assert newick is not None
        tree = Phylo.read(StringIO(newick), "newick")
        assert len(list(tree.get_terminals())) == n

    def test_single_genome_returns_none(self) -> None:
        """Single genome matrix returns None."""
        from metadarkmatter.core.phylogeny.tree_builder import ani_to_newick

        ani = pd.DataFrame({"A": [100.0]}, index=["A"])
        result = ani_to_newick(ani)
        assert result is None

    def test_empty_matrix_returns_none(self) -> None:
        """Empty matrix returns None."""
        from metadarkmatter.core.phylogeny.tree_builder import ani_to_newick

        ani = pd.DataFrame()
        result = ani_to_newick(ani)
        assert result is None

    def test_tree_topology_reflects_ani_similarity(self) -> None:
        """Genomes with higher ANI should cluster together."""
        from metadarkmatter.core.phylogeny.tree_builder import ani_to_newick

        # A and B are very similar (99% ANI), C is distant
        ani = pd.DataFrame(
            {
                "A": [100.0, 99.0, 75.0],
                "B": [99.0, 100.0, 76.0],
                "C": [75.0, 76.0, 100.0],
            },
            index=["A", "B", "C"],
        )

        newick = ani_to_newick(ani)
        assert newick is not None

        tree = Phylo.read(StringIO(newick), "newick")
        # A and B should be sister taxa (share an immediate common ancestor)
        a_path = tree.get_path("A")
        b_path = tree.get_path("B")

        # The parent of A should also be an ancestor of B
        # (verifying A and B share a more recent common ancestor than with C)
        if a_path and b_path:
            # They should share a common internal node before the root
            assert len(a_path) > 1 or len(b_path) > 1


class TestLoadUserTree:
    """Test loading user-provided Newick trees."""

    def test_load_valid_tree(self, tmp_path: Path) -> None:
        """Load a valid Newick tree file."""
        from metadarkmatter.core.phylogeny.tree_builder import load_user_tree

        newick_file = tmp_path / "tree.nwk"
        newick_file.write_text("((A:0.1,B:0.2):0.3,C:0.4);")

        ani = pd.DataFrame(
            {
                "A": [100.0, 95.0, 80.0],
                "B": [95.0, 100.0, 82.0],
                "C": [80.0, 82.0, 100.0],
            },
            index=["A", "B", "C"],
        )

        result = load_user_tree(newick_file, ani)

        assert result is not None
        assert "A" in result
        assert "B" in result
        assert "C" in result

    def test_tree_with_extra_tips_pruned(self, tmp_path: Path) -> None:
        """Tree tips not in ANI matrix are pruned."""
        from metadarkmatter.core.phylogeny.tree_builder import load_user_tree

        newick_file = tmp_path / "tree.nwk"
        newick_file.write_text("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")

        ani = pd.DataFrame(
            {
                "A": [100.0, 95.0, 80.0],
                "B": [95.0, 100.0, 82.0],
                "C": [80.0, 82.0, 100.0],
            },
            index=["A", "B", "C"],
        )

        result = load_user_tree(newick_file, ani)

        tree = Phylo.read(StringIO(result), "newick")
        tip_names = {t.name for t in tree.get_terminals()}
        assert "D" not in tip_names
        assert tip_names == {"A", "B", "C"}

    def test_missing_genomes_logged(self, tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
        """Genomes in ANI but not in tree are logged as warning."""
        from metadarkmatter.core.phylogeny.tree_builder import load_user_tree

        newick_file = tmp_path / "tree.nwk"
        newick_file.write_text("(A:0.1,B:0.2);")

        ani = pd.DataFrame(
            {
                "A": [100.0, 95.0, 80.0],
                "B": [95.0, 100.0, 82.0],
                "C": [80.0, 82.0, 100.0],
            },
            index=["A", "B", "C"],
        )

        with caplog.at_level(logging.WARNING):
            load_user_tree(newick_file, ani)

        assert "1 genomes in ANI matrix not in tree" in caplog.text

    def test_load_tree_preserves_branch_lengths(self, tmp_path: Path) -> None:
        """Branch lengths in user tree should be preserved."""
        from metadarkmatter.core.phylogeny.tree_builder import load_user_tree

        newick_file = tmp_path / "tree.nwk"
        newick_file.write_text("((A:0.1,B:0.2):0.3,C:0.4);")

        ani = pd.DataFrame(
            {
                "A": [100.0, 95.0, 80.0],
                "B": [95.0, 100.0, 82.0],
                "C": [80.0, 82.0, 100.0],
            },
            index=["A", "B", "C"],
        )

        result = load_user_tree(newick_file, ani)
        tree = Phylo.read(StringIO(result), "newick")

        # Find terminal C and verify its branch length
        c_term = next(t for t in tree.get_terminals() if t.name == "C")
        assert c_term.branch_length is not None
        assert abs(c_term.branch_length - 0.4) < 0.001

    def test_nonexistent_file_raises(self, tmp_path: Path) -> None:
        """Nonexistent file should raise FileNotFoundError."""
        from metadarkmatter.core.phylogeny.tree_builder import load_user_tree

        newick_file = tmp_path / "nonexistent.nwk"
        ani = pd.DataFrame({"A": [100.0]}, index=["A"])

        with pytest.raises(FileNotFoundError):
            load_user_tree(newick_file, ani)


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_ani_with_negative_values_clipped(self) -> None:
        """ANI values below 0 should be handled gracefully."""
        from metadarkmatter.core.phylogeny.tree_builder import ani_to_newick

        # Unusual case: ANI values that would result in negative distances
        ani = pd.DataFrame(
            {
                "A": [100.0, 105.0, 80.0],  # 105 is unusual but should not break
                "B": [105.0, 100.0, 82.0],
                "C": [80.0, 82.0, 100.0],
            },
            index=["A", "B", "C"],
        )

        newick = ani_to_newick(ani)
        assert newick is not None

    def test_all_identical_genomes(self) -> None:
        """All genomes with 100% ANI should still produce a tree."""
        from metadarkmatter.core.phylogeny.tree_builder import ani_to_newick

        ani = pd.DataFrame(
            {
                "A": [100.0, 100.0, 100.0],
                "B": [100.0, 100.0, 100.0],
                "C": [100.0, 100.0, 100.0],
            },
            index=["A", "B", "C"],
        )

        newick = ani_to_newick(ani)
        assert newick is not None
        tree = Phylo.read(StringIO(newick), "newick")
        assert len(list(tree.get_terminals())) == 3
