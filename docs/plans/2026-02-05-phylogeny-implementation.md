# Phylogenetic Tree Visualization Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add an interactive phylogenetic tree tab to the HTML report showing reference genomes and novel clusters in evolutionary context.

**Architecture:** Create a new `core/phylogeny/` module with tree building (ANI-to-Newick) and novel cluster placement. Integrate with the existing report generator by adding a new `_build_phylogeny_section()` method. Use Phylotree.js (embedded as inline JavaScript) for interactive visualization.

**Tech Stack:** BioPython (Bio.Phylo), Polars, Phylotree.js (D3-based), existing report HTML templates.

---

## Task 1: Create phylogeny module with tree builder

**Files:**
- Create: `src/metadarkmatter/core/phylogeny/__init__.py`
- Create: `src/metadarkmatter/core/phylogeny/tree_builder.py`
- Create: `tests/unit/test_phylogeny_tree_builder.py`

**Step 1: Create the test file with failing tests**

```python
# tests/unit/test_phylogeny_tree_builder.py
"""Tests for phylogeny tree builder module."""

from __future__ import annotations

from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from Bio import Phylo


class TestAniToNewick:
    """Test ANI matrix to Newick tree conversion."""

    def test_basic_3x3_matrix(self):
        """3x3 ANI matrix produces valid Newick with all taxa."""
        from metadarkmatter.core.phylogeny.tree_builder import ani_to_newick

        ani = pd.DataFrame({
            "A": [100.0, 95.0, 80.0],
            "B": [95.0, 100.0, 82.0],
            "C": [80.0, 82.0, 100.0],
        }, index=["A", "B", "C"])

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

    def test_missing_ani_values_filled(self):
        """Missing ANI values are filled with maximum distance."""
        from metadarkmatter.core.phylogeny.tree_builder import ani_to_newick

        ani = pd.DataFrame({
            "A": [100.0, 95.0, np.nan],
            "B": [95.0, 100.0, 82.0],
            "C": [np.nan, 82.0, 100.0],
        }, index=["A", "B", "C"])

        newick = ani_to_newick(ani)

        assert newick is not None
        tree = Phylo.read(StringIO(newick), "newick")
        assert len(list(tree.get_terminals())) == 3

    def test_fewer_than_3_genomes_returns_none(self):
        """Fewer than 3 genomes returns None (NJ requires >= 3)."""
        from metadarkmatter.core.phylogeny.tree_builder import ani_to_newick

        ani = pd.DataFrame({
            "A": [100.0, 95.0],
            "B": [95.0, 100.0],
        }, index=["A", "B"])

        result = ani_to_newick(ani)
        assert result is None

    def test_larger_matrix(self):
        """Larger matrix (10 genomes) produces valid tree."""
        from metadarkmatter.core.phylogeny.tree_builder import ani_to_newick

        n = 10
        names = [f"G{i}" for i in range(n)]
        # Create symmetric ANI matrix with realistic values
        ani_values = np.full((n, n), 80.0)
        np.fill_diagonal(ani_values, 100.0)
        # Make some pairs more similar
        ani_values[0, 1] = ani_values[1, 0] = 98.0
        ani_values[2, 3] = ani_values[3, 2] = 97.0

        ani = pd.DataFrame(ani_values, index=names, columns=names)
        newick = ani_to_newick(ani)

        assert newick is not None
        tree = Phylo.read(StringIO(newick), "newick")
        assert len(list(tree.get_terminals())) == n


class TestLoadUserTree:
    """Test loading user-provided Newick trees."""

    def test_load_valid_tree(self, tmp_path):
        """Load a valid Newick tree file."""
        from metadarkmatter.core.phylogeny.tree_builder import load_user_tree

        newick_file = tmp_path / "tree.nwk"
        newick_file.write_text("((A:0.1,B:0.2):0.3,C:0.4);")

        ani = pd.DataFrame({
            "A": [100.0, 95.0, 80.0],
            "B": [95.0, 100.0, 82.0],
            "C": [80.0, 82.0, 100.0],
        }, index=["A", "B", "C"])

        result = load_user_tree(newick_file, ani)

        assert result is not None
        assert "A" in result
        assert "B" in result
        assert "C" in result

    def test_tree_with_extra_tips_pruned(self, tmp_path):
        """Tree tips not in ANI matrix are pruned."""
        from metadarkmatter.core.phylogeny.tree_builder import load_user_tree

        newick_file = tmp_path / "tree.nwk"
        newick_file.write_text("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")

        # ANI matrix only has A, B, C (not D)
        ani = pd.DataFrame({
            "A": [100.0, 95.0, 80.0],
            "B": [95.0, 100.0, 82.0],
            "C": [80.0, 82.0, 100.0],
        }, index=["A", "B", "C"])

        result = load_user_tree(newick_file, ani)

        tree = Phylo.read(StringIO(result), "newick")
        tip_names = {t.name for t in tree.get_terminals()}
        assert "D" not in tip_names
        assert tip_names == {"A", "B", "C"}

    def test_missing_genomes_logged(self, tmp_path, caplog):
        """Genomes in ANI but not in tree are logged as warning."""
        from metadarkmatter.core.phylogeny.tree_builder import load_user_tree
        import logging

        newick_file = tmp_path / "tree.nwk"
        newick_file.write_text("(A:0.1,B:0.2);")

        # ANI matrix has A, B, C but tree only has A, B
        ani = pd.DataFrame({
            "A": [100.0, 95.0, 80.0],
            "B": [95.0, 100.0, 82.0],
            "C": [80.0, 82.0, 100.0],
        }, index=["A", "B", "C"])

        with caplog.at_level(logging.WARNING):
            load_user_tree(newick_file, ani)

        assert "1 genomes in ANI matrix not in tree" in caplog.text
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/unit/test_phylogeny_tree_builder.py -v`
Expected: FAIL with "ModuleNotFoundError: No module named 'metadarkmatter.core.phylogeny'"

**Step 3: Create the phylogeny module**

```python
# src/metadarkmatter/core/phylogeny/__init__.py
"""Phylogeny module for tree building and novel cluster placement."""

from metadarkmatter.core.phylogeny.tree_builder import (
    ani_to_newick,
    load_user_tree,
)

__all__ = [
    "ani_to_newick",
    "load_user_tree",
]
```

```python
# src/metadarkmatter/core/phylogeny/tree_builder.py
"""Build phylogenetic trees from ANI matrices."""

from __future__ import annotations

import logging
from io import StringIO
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    pass

logger = logging.getLogger(__name__)


def _to_lower_triangular(matrix: pd.DataFrame) -> list[list[float]]:
    """Convert symmetric matrix to lower triangular list format for BioPython."""
    n = len(matrix)
    result = []
    for i in range(n):
        row = []
        for j in range(i + 1):
            row.append(float(matrix.iloc[i, j]))
        result.append(row)
    return result


def ani_to_newick(ani_matrix: pd.DataFrame) -> str | None:
    """Convert ANI matrix to neighbor-joining tree in Newick format.

    Args:
        ani_matrix: Square DataFrame with ANI values (0-100 scale).
                   Index and columns should be genome identifiers.

    Returns:
        Newick-formatted tree string with branch lengths, or None if < 3 genomes.
    """
    from Bio import Phylo
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

    if len(ani_matrix) < 3:
        logger.warning(
            "Too few genomes for phylogenetic tree (need >= 3). "
            f"Got {len(ani_matrix)}."
        )
        return None

    # Convert ANI (similarity) to distance
    distance_matrix = 100 - ani_matrix
    distance_matrix = distance_matrix.clip(lower=0, upper=50)

    # Handle missing values
    missing_count = distance_matrix.isna().sum().sum() // 2
    if missing_count > 0:
        logger.warning(
            f"{missing_count} genome pairs lack ANI values; "
            "using maximum distance (50)"
        )
        distance_matrix = distance_matrix.fillna(50.0)

    # Build neighbor-joining tree
    names = list(ani_matrix.columns)
    matrix = _to_lower_triangular(distance_matrix)
    dm = DistanceMatrix(names, matrix)

    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    # Root at midpoint
    tree.root_at_midpoint()

    # Export as Newick
    output = StringIO()
    Phylo.write(tree, output, "newick")
    return output.getvalue().strip()


def _prune_tips(tree, keep: set[str]):
    """Prune tree to only keep specified tip names."""
    from Bio import Phylo
    from Bio.Phylo.BaseTree import Tree

    # Get tips to remove
    tips_to_remove = [
        tip for tip in tree.get_terminals()
        if tip.name not in keep
    ]

    for tip in tips_to_remove:
        tree.prune(tip)

    return tree


def load_user_tree(newick_path: Path, ani_matrix: pd.DataFrame) -> str:
    """Load and validate user-provided Newick tree.

    Args:
        newick_path: Path to Newick file.
        ani_matrix: ANI matrix for validation.

    Returns:
        Validated Newick string (possibly pruned).
    """
    from Bio import Phylo

    tree = Phylo.read(newick_path, "newick")

    # Validate tip names against ANI matrix
    tree_tips = {tip.name for tip in tree.get_terminals()}
    matrix_genomes = set(ani_matrix.columns)

    missing = matrix_genomes - tree_tips
    extra = tree_tips - matrix_genomes

    if missing:
        logger.warning(
            f"{len(missing)} genomes in ANI matrix not in tree - "
            "excluded from phylogeny"
        )

    if extra:
        tree = _prune_tips(tree, keep=matrix_genomes)
        logger.info(f"Pruned {len(extra)} unused tips from provided tree")

    # Export as Newick
    output = StringIO()
    Phylo.write(tree, output, "newick")
    return output.getvalue().strip()
```

**Step 4: Run tests to verify they pass**

Run: `pytest tests/unit/test_phylogeny_tree_builder.py -v`
Expected: All tests PASS

**Step 5: Commit**

```bash
git add src/metadarkmatter/core/phylogeny/ tests/unit/test_phylogeny_tree_builder.py
git commit -m "feat(phylogeny): add tree builder module for ANI-to-Newick conversion"
```

---

## Task 2: Add novel cluster placement

**Files:**
- Create: `src/metadarkmatter/core/phylogeny/placement.py`
- Modify: `src/metadarkmatter/core/phylogeny/__init__.py`
- Create: `tests/unit/test_phylogeny_placement.py`

**Step 1: Create the test file with failing tests**

```python
# tests/unit/test_phylogeny_placement.py
"""Tests for novel cluster placement on phylogenetic trees."""

from __future__ import annotations

from io import StringIO

import pandas as pd
import pytest
from Bio import Phylo


class TestNovelCluster:
    """Test NovelCluster dataclass."""

    def test_create_novel_cluster(self):
        """Create a NovelCluster instance."""
        from metadarkmatter.core.phylogeny.placement import NovelCluster

        cluster = NovelCluster(
            cluster_id="Cluster_1",
            classification="novel_species",
            best_match_genome="GCF_000001.1",
            mean_identity=88.5,
            mean_novelty=11.5,
            mean_uncertainty=1.2,
            read_count=150,
            confidence_rating="High",
        )

        assert cluster.cluster_id == "Cluster_1"
        assert cluster.classification == "novel_species"
        assert cluster.read_count == 150


class TestPlaceNovelClusters:
    """Test placing novel clusters on trees."""

    @pytest.fixture
    def simple_tree(self):
        """Create a simple 3-taxon tree."""
        newick = "((A:5,B:5):3,C:8);"
        return Phylo.read(StringIO(newick), "newick")

    @pytest.fixture
    def sample_ani_matrix(self):
        """Create sample ANI matrix."""
        return pd.DataFrame({
            "A": [100.0, 95.0, 80.0],
            "B": [95.0, 100.0, 82.0],
            "C": [80.0, 82.0, 100.0],
        }, index=["A", "B", "C"])

    def test_place_single_novel_species(self, simple_tree, sample_ani_matrix):
        """Novel species cluster is inserted near its best match."""
        from metadarkmatter.core.phylogeny.placement import (
            NovelCluster,
            place_novel_clusters,
        )

        cluster = NovelCluster(
            cluster_id="Novel_1",
            classification="novel_species",
            best_match_genome="A",
            mean_identity=88.0,
            mean_novelty=12.0,
            mean_uncertainty=1.5,
            read_count=100,
            confidence_rating="High",
        )

        result = place_novel_clusters(
            simple_tree, [cluster], sample_ani_matrix
        )

        tip_names = [t.name for t in result.get_terminals()]
        assert "[NOVEL] Novel_1" in tip_names
        assert len(tip_names) == 4  # Original 3 + 1 novel

    def test_place_multiple_clusters(self, simple_tree, sample_ani_matrix):
        """Multiple novel clusters are all inserted."""
        from metadarkmatter.core.phylogeny.placement import (
            NovelCluster,
            place_novel_clusters,
        )

        clusters = [
            NovelCluster(
                cluster_id="Novel_1",
                classification="novel_species",
                best_match_genome="A",
                mean_identity=88.0,
                mean_novelty=12.0,
                mean_uncertainty=1.5,
                read_count=100,
                confidence_rating="High",
            ),
            NovelCluster(
                cluster_id="Novel_2",
                classification="novel_genus",
                best_match_genome="C",
                mean_identity=78.0,
                mean_novelty=22.0,
                mean_uncertainty=1.8,
                read_count=50,
                confidence_rating="Medium",
            ),
        ]

        result = place_novel_clusters(
            simple_tree, clusters, sample_ani_matrix
        )

        tip_names = [t.name for t in result.get_terminals()]
        assert "[NOVEL] Novel_1" in tip_names
        assert "[NOVEL] Novel_2" in tip_names
        assert len(tip_names) == 5

    def test_novel_node_has_metadata(self, simple_tree, sample_ani_matrix):
        """Novel nodes have metadata attached as comment."""
        import json
        from metadarkmatter.core.phylogeny.placement import (
            NovelCluster,
            place_novel_clusters,
        )

        cluster = NovelCluster(
            cluster_id="Novel_1",
            classification="novel_species",
            best_match_genome="A",
            mean_identity=88.0,
            mean_novelty=12.0,
            mean_uncertainty=1.5,
            read_count=100,
            confidence_rating="High",
        )

        result = place_novel_clusters(
            simple_tree, [cluster], sample_ani_matrix
        )

        # Find the novel node
        novel_node = None
        for tip in result.get_terminals():
            if tip.name == "[NOVEL] Novel_1":
                novel_node = tip
                break

        assert novel_node is not None
        assert novel_node.comment is not None

        metadata = json.loads(novel_node.comment)
        assert metadata["type"] == "novel_species"
        assert metadata["read_count"] == 100
        assert metadata["is_novel"] is True

    def test_empty_cluster_list(self, simple_tree, sample_ani_matrix):
        """Empty cluster list returns tree unchanged."""
        from metadarkmatter.core.phylogeny.placement import place_novel_clusters

        result = place_novel_clusters(simple_tree, [], sample_ani_matrix)

        tip_names = [t.name for t in result.get_terminals()]
        assert len(tip_names) == 3  # Unchanged


class TestExtractNovelClustersFromClassifications:
    """Test extracting novel clusters from classification DataFrame."""

    def test_extract_from_classifications(self):
        """Extract novel clusters from classification results."""
        import polars as pl
        from metadarkmatter.core.phylogeny.placement import (
            extract_novel_clusters,
        )

        df = pl.DataFrame({
            "read_id": ["r1", "r2", "r3", "r4", "r5"],
            "best_match_genome": ["G1", "G1", "G1", "G2", "G2"],
            "top_hit_identity": [88.0, 87.5, 89.0, 78.0, 77.0],
            "novelty_index": [12.0, 12.5, 11.0, 22.0, 23.0],
            "placement_uncertainty": [1.5, 1.6, 1.4, 1.8, 1.9],
            "taxonomic_call": [
                "Novel Species", "Novel Species", "Novel Species",
                "Novel Genus", "Novel Genus"
            ],
        })

        clusters = extract_novel_clusters(df, min_reads=2)

        assert len(clusters) == 2

        # Find the Novel Species cluster
        species_cluster = next(c for c in clusters if c.classification == "novel_species")
        assert species_cluster.read_count == 3
        assert species_cluster.best_match_genome == "G1"

        # Find the Novel Genus cluster
        genus_cluster = next(c for c in clusters if c.classification == "novel_genus")
        assert genus_cluster.read_count == 2
        assert genus_cluster.best_match_genome == "G2"

    def test_filter_by_min_reads(self):
        """Clusters below min_reads threshold are excluded."""
        import polars as pl
        from metadarkmatter.core.phylogeny.placement import (
            extract_novel_clusters,
        )

        df = pl.DataFrame({
            "read_id": ["r1", "r2"],
            "best_match_genome": ["G1", "G2"],
            "top_hit_identity": [88.0, 78.0],
            "novelty_index": [12.0, 22.0],
            "placement_uncertainty": [1.5, 1.8],
            "taxonomic_call": ["Novel Species", "Novel Genus"],
        })

        clusters = extract_novel_clusters(df, min_reads=3)
        assert len(clusters) == 0
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/unit/test_phylogeny_placement.py -v`
Expected: FAIL with "ImportError"

**Step 3: Implement placement module**

```python
# src/metadarkmatter/core/phylogeny/placement.py
"""Place novel clusters on phylogenetic trees."""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from Bio.Phylo.BaseTree import Tree
    import polars as pl

logger = logging.getLogger(__name__)


@dataclass
class NovelCluster:
    """Represents a cluster of novel reads."""

    cluster_id: str
    classification: str  # "novel_species" or "novel_genus"
    best_match_genome: str
    mean_identity: float
    mean_novelty: float
    mean_uncertainty: float
    read_count: int
    confidence_rating: str  # "High", "Medium", "Low"


def _find_clade_by_name(tree, name: str):
    """Find a clade (tip) by its name."""
    for tip in tree.get_terminals():
        if tip.name == name:
            return tip
    return None


def _get_parent(tree, child):
    """Get the parent clade of a given clade."""
    path = tree.get_path(child)
    if len(path) >= 2:
        return path[-2]
    return tree.root


def _insert_as_sibling(tree, target, new_node, distance: float):
    """Insert new_node as a sibling to target.

    Creates a new internal node that becomes parent to both target and new_node.
    """
    from Bio.Phylo.BaseTree import Clade

    parent = _get_parent(tree, target)

    # Remove target from parent
    parent.clades.remove(target)

    # Create new internal node
    internal = Clade(branch_length=target.branch_length / 2)

    # Adjust target branch length
    target.branch_length = target.branch_length / 2

    # Set new node branch length based on distance
    new_node.branch_length = distance / 2

    # Add both as children of new internal node
    internal.clades = [target, new_node]

    # Add internal node to parent
    parent.clades.append(internal)


def place_novel_clusters(
    tree,
    novel_clusters: list[NovelCluster],
    ani_matrix: pd.DataFrame,
    genome_metadata: pd.DataFrame | None = None,
):
    """Insert novel clusters as leaf nodes at estimated positions.

    Args:
        tree: BioPython Phylo tree object.
        novel_clusters: List of NovelCluster objects to place.
        ani_matrix: ANI matrix (used for distance estimation).
        genome_metadata: Optional metadata for genus-level placement.

    Returns:
        Modified tree with novel clusters inserted.
    """
    from Bio.Phylo.BaseTree import Clade

    for cluster in novel_clusters:
        # Find the target clade (nearest reference)
        target = _find_clade_by_name(tree, cluster.best_match_genome)

        if target is None:
            logger.warning(
                f"Could not find genome {cluster.best_match_genome} in tree; "
                f"skipping cluster {cluster.cluster_id}"
            )
            continue

        # Estimate distance from identity
        estimated_distance = 100 - cluster.mean_identity

        # Create novel node
        novel_node = Clade(
            name=f"[NOVEL] {cluster.cluster_id}",
            branch_length=estimated_distance / 2,
        )

        # Attach metadata as JSON comment
        novel_node.comment = json.dumps({
            "type": cluster.classification,
            "read_count": cluster.read_count,
            "mean_novelty": cluster.mean_novelty,
            "mean_uncertainty": cluster.mean_uncertainty,
            "confidence": cluster.confidence_rating,
            "nearest_ref": cluster.best_match_genome,
            "est_ani": cluster.mean_identity,
            "is_novel": True,
        })

        # Insert as sibling to target
        _insert_as_sibling(tree, target, novel_node, estimated_distance)

    return tree


def _calculate_confidence_rating(
    read_count: int,
    mean_uncertainty: float,
    mean_novelty: float,
) -> str:
    """Calculate confidence rating for a cluster."""
    if read_count >= 10 and mean_uncertainty < 5.0:
        return "High"
    elif read_count >= 5 and mean_uncertainty < 10.0:
        return "Medium"
    return "Low"


def extract_novel_clusters(
    classifications: "pl.DataFrame",
    min_reads: int = 3,
) -> list[NovelCluster]:
    """Extract novel clusters from classification results.

    Groups reads by best_match_genome and taxonomic_call to form clusters.

    Args:
        classifications: Polars DataFrame with classification results.
        min_reads: Minimum reads to form a cluster.

    Returns:
        List of NovelCluster objects.
    """
    import polars as pl

    # Filter to novel reads only
    novel_df = classifications.filter(
        pl.col("taxonomic_call").is_in(["Novel Species", "Novel Genus"])
    )

    if len(novel_df) == 0:
        return []

    # Group by genome and classification type
    clusters = (
        novel_df.group_by(["best_match_genome", "taxonomic_call"])
        .agg([
            pl.len().alias("read_count"),
            pl.col("top_hit_identity").mean().alias("mean_identity"),
            pl.col("novelty_index").mean().alias("mean_novelty"),
            pl.col("placement_uncertainty").mean().alias("mean_uncertainty"),
        ])
        .filter(pl.col("read_count") >= min_reads)
        .sort("read_count", descending=True)
    )

    result = []
    cluster_num = 1

    for row in clusters.iter_rows(named=True):
        classification = (
            "novel_species" if row["taxonomic_call"] == "Novel Species"
            else "novel_genus"
        )

        confidence = _calculate_confidence_rating(
            row["read_count"],
            row["mean_uncertainty"],
            row["mean_novelty"],
        )

        result.append(NovelCluster(
            cluster_id=f"Cluster_{cluster_num}",
            classification=classification,
            best_match_genome=row["best_match_genome"],
            mean_identity=row["mean_identity"],
            mean_novelty=row["mean_novelty"],
            mean_uncertainty=row["mean_uncertainty"],
            read_count=row["read_count"],
            confidence_rating=confidence,
        ))
        cluster_num += 1

    return result
```

**Step 4: Update __init__.py**

```python
# src/metadarkmatter/core/phylogeny/__init__.py
"""Phylogeny module for tree building and novel cluster placement."""

from metadarkmatter.core.phylogeny.tree_builder import (
    ani_to_newick,
    load_user_tree,
)
from metadarkmatter.core.phylogeny.placement import (
    NovelCluster,
    extract_novel_clusters,
    place_novel_clusters,
)

__all__ = [
    "ani_to_newick",
    "load_user_tree",
    "NovelCluster",
    "extract_novel_clusters",
    "place_novel_clusters",
]
```

**Step 5: Run tests to verify they pass**

Run: `pytest tests/unit/test_phylogeny_placement.py tests/unit/test_phylogeny_tree_builder.py -v`
Expected: All tests PASS

**Step 6: Commit**

```bash
git add src/metadarkmatter/core/phylogeny/ tests/unit/test_phylogeny_placement.py
git commit -m "feat(phylogeny): add novel cluster placement on trees"
```

---

## Task 3: Add Phylotree.js template and styles

**Files:**
- Modify: `src/metadarkmatter/visualization/report/templates.py`
- Modify: `src/metadarkmatter/visualization/report/styles.py`

**Step 1: Create test for phylogeny templates**

```python
# Add to tests/unit/test_report_generator.py

class TestPhylogenyTemplates:
    """Test phylogeny-related templates."""

    def test_phylogeny_section_template_exists(self):
        """Phylogeny section template is defined."""
        from metadarkmatter.visualization.report.templates import (
            PHYLOGENY_SECTION_TEMPLATE,
        )
        assert "phylogeny-container" in PHYLOGENY_SECTION_TEMPLATE
        assert "TREE_DATA" in PHYLOGENY_SECTION_TEMPLATE

    def test_phylotree_js_template_exists(self):
        """Phylotree.js code template is defined."""
        from metadarkmatter.visualization.report.templates import (
            PHYLOTREE_JS_TEMPLATE,
        )
        assert "Phylotree" in PHYLOTREE_JS_TEMPLATE or "d3" in PHYLOTREE_JS_TEMPLATE
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_report_generator.py::TestPhylogenyTemplates -v`
Expected: FAIL with "cannot import name 'PHYLOGENY_SECTION_TEMPLATE'"

**Step 3: Add templates to templates.py**

Add to the end of `src/metadarkmatter/visualization/report/templates.py`:

```python
# Phylogeny section template
PHYLOGENY_SECTION_TEMPLATE = """
<div class="phylogeny-section">
    <div class="phylogeny-header">
        <h3>Phylogenetic Tree</h3>
        <p class="section-description">
            Interactive phylogenetic tree showing reference genomes and novel clusters.
            {tree_source_note}
        </p>
    </div>

    <div class="phylogeny-controls">
        <button id="layout-toggle" class="btn btn-secondary">
            Toggle Radial/Rectangular
        </button>
        <button id="expand-all" class="btn btn-secondary">
            Expand All
        </button>
        <button id="collapse-genera" class="btn btn-secondary">
            Collapse to Genera
        </button>
    </div>

    <div class="phylogeny-legend">
        <div class="legend-item">
            <span class="legend-circle reference"></span>
            Reference genome
        </div>
        <div class="legend-item">
            <span class="legend-circle novel-species"></span>
            Novel species cluster
        </div>
        <div class="legend-item">
            <span class="legend-circle novel-genus"></span>
            Novel genus cluster
        </div>
        <div class="legend-item">
            <span class="confidence-badge">H</span> High /
            <span class="confidence-badge medium">M</span> Medium /
            <span class="confidence-badge low">L</span> Low confidence
        </div>
    </div>

    <div id="phylogeny-container" class="phylogeny-tree-container">
        <!-- Tree will be rendered here by JavaScript -->
    </div>

    <div id="phylogeny-tooltip" class="phylogeny-tooltip" style="display: none;"></div>
</div>

<script>
const TREE_DATA = {tree_data_json};
</script>
"""

# D3-based tree visualization (simplified Phylotree-like implementation)
PHYLOTREE_JS_TEMPLATE = """
<script src="https://d3js.org/d3.v7.min.js"></script>
<script>
(function() {
    const treeData = TREE_DATA;
    if (!treeData || !treeData.newick) {
        document.getElementById('phylogeny-container').innerHTML =
            '<p class="warning">No phylogenetic tree data available.</p>';
        return;
    }

    // Parse Newick string into hierarchical structure
    function parseNewick(newick) {
        const ancestors = [];
        let tree = {};
        const tokens = newick.split(/\\s*(;|\\(|\\)|,|:)\\s*/);
        let subtree;

        for (let i = 0; i < tokens.length; i++) {
            const token = tokens[i];
            switch (token) {
                case '(':
                    subtree = {};
                    tree.children = [subtree];
                    ancestors.push(tree);
                    tree = subtree;
                    break;
                case ',':
                    subtree = {};
                    ancestors[ancestors.length - 1].children.push(subtree);
                    tree = subtree;
                    break;
                case ')':
                    tree = ancestors.pop();
                    break;
                case ':':
                    break;
                default:
                    const x = tokens[i - 1];
                    if (x === ')' || x === '(' || x === ',') {
                        tree.name = token;
                    } else if (x === ':') {
                        tree.length = parseFloat(token);
                    }
            }
        }
        return tree;
    }

    // Build tree visualization
    const container = document.getElementById('phylogeny-container');
    const width = container.clientWidth || 900;
    const tipCount = treeData.tip_count || 20;
    const height = Math.max(400, tipCount * 20);

    const svg = d3.select('#phylogeny-container')
        .append('svg')
        .attr('width', width)
        .attr('height', height)
        .append('g')
        .attr('transform', 'translate(40, 20)');

    const root = d3.hierarchy(parseNewick(treeData.newick));
    const treeLayout = d3.cluster().size([height - 40, width - 200]);
    treeLayout(root);

    // Draw links
    svg.selectAll('.link')
        .data(root.links())
        .enter()
        .append('path')
        .attr('class', 'link')
        .attr('d', d => {
            return `M${d.source.y},${d.source.x}
                    L${d.target.y},${d.source.x}
                    L${d.target.y},${d.target.x}`;
        })
        .style('fill', 'none')
        .style('stroke', d => {
            const name = d.target.data.name || '';
            return name.startsWith('[NOVEL]') ? '#f59e0b' : '#999';
        })
        .style('stroke-width', d => {
            const name = d.target.data.name || '';
            return name.startsWith('[NOVEL]') ? 2 : 1;
        })
        .style('stroke-dasharray', d => {
            const name = d.target.data.name || '';
            return name.startsWith('[NOVEL]') ? '4,2' : 'none';
        });

    // Draw nodes
    const nodes = svg.selectAll('.node')
        .data(root.descendants().filter(d => d.data.name))
        .enter()
        .append('g')
        .attr('class', 'node')
        .attr('transform', d => `translate(${d.y},${d.x})`);

    nodes.append('circle')
        .attr('r', d => {
            const name = d.data.name || '';
            const meta = treeData.annotations[name];
            if (meta && meta.is_novel) {
                return Math.sqrt(meta.read_count / 10) + 4;
            }
            return 3;
        })
        .style('fill', d => {
            const name = d.data.name || '';
            const meta = treeData.annotations[name];
            if (meta && meta.is_novel) {
                return meta.type === 'novel_species' ? '#f59e0b' : '#ef4444';
            }
            return '#6b7280';
        })
        .style('stroke', '#fff')
        .style('stroke-width', d => {
            const name = d.data.name || '';
            const meta = treeData.annotations[name];
            return (meta && meta.is_novel) ? 2 : 1;
        });

    // Add labels for tips only
    nodes.filter(d => !d.children)
        .append('text')
        .attr('dx', 8)
        .attr('dy', 3)
        .style('font-size', '11px')
        .style('fill', d => {
            const name = d.data.name || '';
            return name.startsWith('[NOVEL]') ? '#f59e0b' : '#333';
        })
        .text(d => {
            const name = d.data.name || '';
            // Truncate long names
            return name.length > 30 ? name.substring(0, 27) + '...' : name;
        });

    // Add confidence badges for novel clusters
    nodes.filter(d => {
        const meta = treeData.annotations[d.data.name];
        return meta && meta.is_novel;
    })
    .append('text')
    .attr('class', 'confidence-badge-svg')
    .attr('dx', d => {
        const meta = treeData.annotations[d.data.name];
        return Math.sqrt(meta.read_count / 10) + 6;
    })
    .attr('dy', -8)
    .style('font-size', '10px')
    .style('font-weight', 'bold')
    .style('fill', d => {
        const meta = treeData.annotations[d.data.name];
        if (meta.confidence === 'High') return '#16a34a';
        if (meta.confidence === 'Medium') return '#f59e0b';
        return '#6b7280';
    })
    .text(d => {
        const meta = treeData.annotations[d.data.name];
        return meta.confidence[0];
    });

    // Tooltip
    const tooltip = d3.select('#phylogeny-tooltip');

    nodes.on('mouseover', function(event, d) {
        const meta = treeData.annotations[d.data.name];
        if (meta && meta.is_novel) {
            tooltip.style('display', 'block')
                .style('left', (event.pageX + 10) + 'px')
                .style('top', (event.pageY - 10) + 'px')
                .html(`
                    <strong>${d.data.name}</strong><br>
                    Classification: ${meta.type.replace('_', ' ')}<br>
                    Reads: ${meta.read_count.toLocaleString()}<br>
                    Mean Novelty: ${meta.mean_novelty.toFixed(1)}%<br>
                    Mean Uncertainty: ${meta.mean_uncertainty.toFixed(1)}%<br>
                    Confidence: ${meta.confidence}<br>
                    Nearest Reference: ${meta.nearest_ref}<br>
                    Est. ANI: ${meta.est_ani.toFixed(1)}%
                `);
        }
    })
    .on('mouseout', function() {
        tooltip.style('display', 'none');
    });

    // Layout toggle (simplified - just logs for now)
    document.getElementById('layout-toggle').addEventListener('click', function() {
        console.log('Layout toggle clicked - radial layout not yet implemented');
    });

    document.getElementById('expand-all').addEventListener('click', function() {
        console.log('Expand all clicked');
    });

    document.getElementById('collapse-genera').addEventListener('click', function() {
        console.log('Collapse to genera clicked');
    });
})();
</script>
"""

PHYLOGENY_NOT_PROVIDED_MESSAGE = """
<div class="info-message">
    <p><strong>Phylogenetic tree not available.</strong></p>
    <p>To enable the phylogeny tab, provide an ANI matrix with at least 3 genomes.</p>
</div>
"""
```

**Step 4: Add CSS styles to styles.py**

Add to `src/metadarkmatter/visualization/report/styles.py` in the `get_css_styles()` function:

```python
# Add these styles to the CSS string (inside get_css_styles function)

PHYLOGENY_STYLES = """
/* Phylogeny Tab Styles */
.phylogeny-section {
    padding: 20px;
}

.phylogeny-header {
    margin-bottom: 20px;
}

.phylogeny-header h3 {
    margin: 0 0 10px 0;
    color: var(--text-primary);
}

.phylogeny-controls {
    display: flex;
    gap: 10px;
    margin-bottom: 15px;
}

.phylogeny-controls .btn {
    padding: 8px 16px;
    border: 1px solid var(--border-color);
    border-radius: 4px;
    background: var(--bg-secondary);
    color: var(--text-primary);
    cursor: pointer;
    font-size: 13px;
}

.phylogeny-controls .btn:hover {
    background: var(--bg-tertiary);
}

.phylogeny-legend {
    display: flex;
    flex-wrap: wrap;
    gap: 20px;
    margin-bottom: 20px;
    padding: 10px 15px;
    background: var(--bg-secondary);
    border-radius: 6px;
    font-size: 13px;
}

.legend-item {
    display: flex;
    align-items: center;
    gap: 6px;
}

.legend-circle {
    width: 12px;
    height: 12px;
    border-radius: 50%;
    border: 2px solid #fff;
    box-shadow: 0 0 0 1px var(--border-color);
}

.legend-circle.reference {
    background: #6b7280;
}

.legend-circle.novel-species {
    background: #f59e0b;
}

.legend-circle.novel-genus {
    background: #ef4444;
}

.confidence-badge {
    display: inline-block;
    padding: 2px 6px;
    border-radius: 3px;
    font-size: 11px;
    font-weight: bold;
    background: #16a34a;
    color: white;
}

.confidence-badge.medium {
    background: #f59e0b;
}

.confidence-badge.low {
    background: #6b7280;
}

.phylogeny-tree-container {
    border: 1px solid var(--border-color);
    border-radius: 6px;
    background: var(--bg-primary);
    overflow: auto;
    min-height: 400px;
}

.phylogeny-tree-container svg {
    display: block;
}

.phylogeny-tooltip {
    position: absolute;
    padding: 10px 14px;
    background: var(--bg-primary);
    border: 1px solid var(--border-color);
    border-radius: 6px;
    box-shadow: 0 4px 12px rgba(0,0,0,0.15);
    font-size: 12px;
    line-height: 1.5;
    z-index: 1000;
    max-width: 300px;
}

.phylogeny-tooltip strong {
    color: var(--text-primary);
}
"""
```

**Step 5: Run tests**

Run: `pytest tests/unit/test_report_generator.py::TestPhylogenyTemplates -v`
Expected: PASS

**Step 6: Commit**

```bash
git add src/metadarkmatter/visualization/report/templates.py src/metadarkmatter/visualization/report/styles.py tests/unit/test_report_generator.py
git commit -m "feat(report): add phylogeny tab templates and styles"
```

---

## Task 4: Integrate phylogeny section into report generator

**Files:**
- Modify: `src/metadarkmatter/visualization/report/generator.py`
- Add tests to: `tests/unit/test_report_generator.py`

**Step 1: Add integration test**

Add to `tests/unit/test_report_generator.py`:

```python
class TestPhylogenySection:
    """Test phylogeny section generation."""

    def test_build_phylogeny_section_with_ani(self, sample_classifications, sample_ani_matrix):
        """Phylogeny section is built when ANI matrix provided."""
        from metadarkmatter.visualization.report.generator import (
            ReportConfig,
            ReportGenerator,
        )
        from metadarkmatter.visualization.plots.base import PlotConfig, ThresholdConfig

        config = ReportConfig(
            sample_name="Test",
            title="Test Report",
            plot_config=PlotConfig(),
            thresholds=ThresholdConfig(),
        )

        generator = ReportGenerator(
            classifications=sample_classifications,
            config=config,
            ani_matrix=sample_ani_matrix,
        )

        # Convert Polars to Pandas for the phylogeny module
        import pandas as pd
        ani_pd = sample_ani_matrix.to_pandas()
        ani_pd = ani_pd.set_index("genome")

        section = generator._build_phylogeny_section(ani_pd)

        assert section is not None
        assert "phylogeny-container" in section
        assert "TREE_DATA" in section

    def test_phylogeny_section_skipped_without_ani(self, sample_classifications):
        """Phylogeny section returns None without ANI matrix."""
        from metadarkmatter.visualization.report.generator import (
            ReportConfig,
            ReportGenerator,
        )
        from metadarkmatter.visualization.plots.base import PlotConfig, ThresholdConfig

        config = ReportConfig(
            sample_name="Test",
            title="Test Report",
            plot_config=PlotConfig(),
            thresholds=ThresholdConfig(),
        )

        generator = ReportGenerator(
            classifications=sample_classifications,
            config=config,
            ani_matrix=None,
        )

        section = generator._build_phylogeny_section(None)
        assert section is None

    def test_phylogeny_section_with_novel_clusters(self, sample_ani_matrix):
        """Phylogeny section includes novel clusters when present."""
        import polars as pl
        from metadarkmatter.visualization.report.generator import (
            ReportConfig,
            ReportGenerator,
        )
        from metadarkmatter.visualization.plots.base import PlotConfig, ThresholdConfig

        # Create classifications with novel reads
        classifications = pl.DataFrame({
            "read_id": [f"read_{i}" for i in range(100)],
            "best_match_genome": ["GCF_000001.1"] * 50 + ["GCF_000002.1"] * 50,
            "top_hit_identity": [99.5] * 50 + [88.0] * 50,
            "novelty_index": [0.5] * 50 + [12.0] * 50,
            "placement_uncertainty": [0.2] * 50 + [1.5] * 50,
            "num_ambiguous_hits": [1] * 100,
            "taxonomic_call": ["Known Species"] * 50 + ["Novel Species"] * 50,
            "is_novel": [False] * 50 + [True] * 50,
        })

        config = ReportConfig(
            sample_name="Test",
            title="Test Report",
            plot_config=PlotConfig(),
            thresholds=ThresholdConfig(),
        )

        generator = ReportGenerator(
            classifications=classifications,
            config=config,
            ani_matrix=sample_ani_matrix,
        )

        import pandas as pd
        ani_pd = sample_ani_matrix.to_pandas()
        ani_pd = ani_pd.set_index("genome")

        section = generator._build_phylogeny_section(ani_pd)

        assert "[NOVEL]" in section
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_report_generator.py::TestPhylogenySection -v`
Expected: FAIL

**Step 3: Add _build_phylogeny_section to generator.py**

Add this method to the `ReportGenerator` class in `src/metadarkmatter/visualization/report/generator.py`:

```python
def _build_phylogeny_section(
    self,
    ani_matrix_pd: "pd.DataFrame | None",
    user_tree_path: Path | None = None,
) -> str | None:
    """Build interactive phylogenetic tree section.

    Args:
        ani_matrix_pd: Pandas DataFrame with ANI values (genome x genome).
        user_tree_path: Optional path to user-provided Newick tree.

    Returns:
        HTML string for Phylogeny tab content, or None if not available.
    """
    import json

    from metadarkmatter.visualization.report.templates import (
        PHYLOGENY_NOT_PROVIDED_MESSAGE,
        PHYLOGENY_SECTION_TEMPLATE,
        PHYLOTREE_JS_TEMPLATE,
    )

    if ani_matrix_pd is None or len(ani_matrix_pd) < 3:
        return None

    try:
        from metadarkmatter.core.phylogeny import (
            ani_to_newick,
            extract_novel_clusters,
            load_user_tree,
            place_novel_clusters,
        )
        from Bio import Phylo
        from io import StringIO

        # Build or load tree
        if user_tree_path is not None:
            newick = load_user_tree(user_tree_path, ani_matrix_pd)
            tree_source_note = "Tree source: user-provided Newick file."
        else:
            newick = ani_to_newick(ani_matrix_pd)
            tree_source_note = "Tree source: neighbor-joining from ANI matrix."

        if newick is None:
            return None

        # Extract novel clusters from classifications
        novel_clusters = extract_novel_clusters(
            self.classifications,
            min_reads=3,
        )

        # Place novel clusters on tree
        tree = Phylo.read(StringIO(newick), "newick")
        if novel_clusters:
            tree = place_novel_clusters(
                tree, novel_clusters, ani_matrix_pd
            )

        # Export modified tree as Newick
        output = StringIO()
        Phylo.write(tree, output, "newick")
        final_newick = output.getvalue().strip()

        # Build annotations dict for JavaScript
        annotations = {}

        # Add annotations for novel clusters
        for cluster in novel_clusters:
            node_name = f"[NOVEL] {cluster.cluster_id}"
            annotations[node_name] = {
                "type": cluster.classification,
                "read_count": cluster.read_count,
                "mean_novelty": cluster.mean_novelty,
                "mean_uncertainty": cluster.mean_uncertainty,
                "confidence": cluster.confidence_rating,
                "nearest_ref": cluster.best_match_genome,
                "est_ani": cluster.mean_identity,
                "is_novel": True,
            }

        # Count tips
        tip_count = len(list(tree.get_terminals()))

        # Build tree data JSON
        tree_data = {
            "newick": final_newick,
            "annotations": annotations,
            "tip_count": tip_count,
            "novel_count": len(novel_clusters),
        }

        # Render HTML
        section_html = PHYLOGENY_SECTION_TEMPLATE.format(
            tree_source_note=tree_source_note,
            tree_data_json=json.dumps(tree_data),
        )

        return section_html + PHYLOTREE_JS_TEMPLATE

    except ImportError as e:
        logger.warning(f"Could not build phylogeny section: {e}")
        return None
    except Exception as e:
        logger.warning(f"Error building phylogeny section: {e}")
        return None
```

**Step 4: Integrate into generate() method**

Find the `generate()` method in `ReportGenerator` and add the phylogeny tab. Look for where tabs are being built and add:

```python
# In the generate() method, find where tabs are assembled and add:

# Build phylogeny section if ANI matrix available
phylogeny_html = None
if self.ani_matrix is not None:
    # Convert Polars to Pandas for phylogeny module
    ani_pd = self.ani_matrix.to_pandas()
    if "genome" in ani_pd.columns:
        ani_pd = ani_pd.set_index("genome")
    phylogeny_html = self._build_phylogeny_section(ani_pd)

# Add to tabs list (find existing tab assembly code and add):
if phylogeny_html:
    tabs.append(("Phylogeny", "phylogeny", phylogeny_html))
```

**Step 5: Run tests**

Run: `pytest tests/unit/test_report_generator.py::TestPhylogenySection -v`
Expected: PASS

**Step 6: Run full test suite**

Run: `pytest tests/unit/test_report_generator.py -v`
Expected: All PASS

**Step 7: Commit**

```bash
git add src/metadarkmatter/visualization/report/generator.py tests/unit/test_report_generator.py
git commit -m "feat(report): integrate phylogeny section into report generator"
```

---

## Task 5: Add CLI options for phylogeny

**Files:**
- Modify: `src/metadarkmatter/cli/report.py`
- Add tests to: `tests/unit/test_cli_report.py`

**Step 1: Add CLI test**

Add to `tests/unit/test_cli_report.py`:

```python
class TestPhylogenyCLI:
    """Test phylogeny-related CLI options."""

    def test_tree_option_exists(self):
        """--tree option is available in generate command."""
        from typer.testing import CliRunner
        from metadarkmatter.cli.report import app

        runner = CliRunner()
        result = runner.invoke(app, ["generate", "--help"])

        assert result.exit_code == 0
        assert "--tree" in result.output

    def test_no_phylogeny_option_exists(self):
        """--no-phylogeny option is available in generate command."""
        from typer.testing import CliRunner
        from metadarkmatter.cli.report import app

        runner = CliRunner()
        result = runner.invoke(app, ["generate", "--help"])

        assert result.exit_code == 0
        assert "--no-phylogeny" in result.output
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_cli_report.py::TestPhylogenyCLI -v`
Expected: FAIL (options don't exist yet)

**Step 3: Add CLI options to report.py**

In `src/metadarkmatter/cli/report.py`, add these parameters to `generate_report()`:

```python
@app.command(name="generate")
def generate_report(
    # ... existing parameters ...

    tree: Path | None = typer.Option(
        None,
        "--tree",
        help="Newick tree file for phylogenetic context. "
             "If not provided, neighbor-joining tree is built from ANI matrix.",
        exists=True,
        dir_okay=False,
    ),
    no_phylogeny: bool = typer.Option(
        False,
        "--no-phylogeny",
        help="Skip phylogeny tab generation (faster for large datasets).",
    ),

    # ... rest of parameters ...
) -> None:
```

Then pass these to the generator in the function body:

```python
# In the ReportConfig creation, add:
config = ReportConfig(
    # ... existing params ...
    skip_phylogeny=no_phylogeny,
    user_tree_path=tree,
)
```

**Step 4: Update ReportConfig dataclass**

In `src/metadarkmatter/visualization/report/generator.py`, add to `ReportConfig`:

```python
@dataclass
class ReportConfig:
    # ... existing fields ...
    skip_phylogeny: bool = False
    user_tree_path: Path | None = None
```

**Step 5: Use config in generate()**

Update the phylogeny section building to respect the config:

```python
# In generate() method:
phylogeny_html = None
if not self.config.skip_phylogeny and self.ani_matrix is not None:
    ani_pd = self.ani_matrix.to_pandas()
    if "genome" in ani_pd.columns:
        ani_pd = ani_pd.set_index("genome")
    phylogeny_html = self._build_phylogeny_section(
        ani_pd,
        user_tree_path=self.config.user_tree_path,
    )
```

**Step 6: Run tests**

Run: `pytest tests/unit/test_cli_report.py::TestPhylogenyCLI -v`
Expected: PASS

**Step 7: Commit**

```bash
git add src/metadarkmatter/cli/report.py src/metadarkmatter/visualization/report/generator.py tests/unit/test_cli_report.py
git commit -m "feat(cli): add --tree and --no-phylogeny options to report generate"
```

---

## Task 6: Integration test with full report generation

**Files:**
- Create: `tests/integration/test_report_phylogeny.py`

**Step 1: Create integration test file**

```python
# tests/integration/test_report_phylogeny.py
"""Integration tests for phylogeny in report generation."""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.report import app

runner = CliRunner()


@pytest.fixture
def sample_data(tmp_path):
    """Create sample data files for testing."""
    # Classifications with novel reads
    classifications = pl.DataFrame({
        "read_id": [f"read_{i}" for i in range(100)],
        "best_match_genome": ["GCF_000001.1"] * 40 + ["GCF_000002.1"] * 30 + ["GCF_000003.1"] * 30,
        "top_hit_identity": [99.5] * 40 + [88.0] * 30 + [78.0] * 30,
        "novelty_index": [0.5] * 40 + [12.0] * 30 + [22.0] * 30,
        "placement_uncertainty": [0.2] * 40 + [1.5] * 30 + [1.8] * 30,
        "num_ambiguous_hits": [1] * 40 + [2] * 30 + [3] * 30,
        "taxonomic_call": ["Known Species"] * 40 + ["Novel Species"] * 30 + ["Novel Genus"] * 30,
        "is_novel": [False] * 40 + [True] * 60,
    })

    # ANI matrix
    ani_matrix = pl.DataFrame({
        "genome": ["GCF_000001.1", "GCF_000002.1", "GCF_000003.1"],
        "GCF_000001.1": [100.0, 85.5, 78.2],
        "GCF_000002.1": [85.5, 100.0, 82.1],
        "GCF_000003.1": [78.2, 82.1, 100.0],
    })

    # Save files
    class Data:
        pass

    data = Data()
    data.classifications = tmp_path / "classifications.csv"
    data.ani = tmp_path / "ani_matrix.csv"
    data.output = tmp_path / "report.html"

    classifications.write_csv(data.classifications)
    ani_matrix.write_csv(data.ani)

    return data


class TestReportPhylogenyIntegration:
    """Integration tests for phylogeny tab in reports."""

    def test_report_with_phylogeny_tab(self, sample_data):
        """Report includes phylogeny tab when ANI matrix provided."""
        result = runner.invoke(app, [
            "generate",
            "--classifications", str(sample_data.classifications),
            "--ani", str(sample_data.ani),
            "--output", str(sample_data.output),
        ])

        assert result.exit_code == 0
        html = sample_data.output.read_text()

        # Check phylogeny tab exists
        assert "Phylogeny" in html
        assert "phylogeny-container" in html
        assert "TREE_DATA" in html

    def test_report_phylogeny_includes_novel_clusters(self, sample_data):
        """Phylogeny tab shows novel clusters."""
        result = runner.invoke(app, [
            "generate",
            "--classifications", str(sample_data.classifications),
            "--ani", str(sample_data.ani),
            "--output", str(sample_data.output),
        ])

        assert result.exit_code == 0
        html = sample_data.output.read_text()

        # Check novel clusters are included
        assert "[NOVEL]" in html

    def test_report_no_phylogeny_flag(self, sample_data):
        """--no-phylogeny skips phylogeny tab."""
        result = runner.invoke(app, [
            "generate",
            "--classifications", str(sample_data.classifications),
            "--ani", str(sample_data.ani),
            "--no-phylogeny",
            "--output", str(sample_data.output),
        ])

        assert result.exit_code == 0
        html = sample_data.output.read_text()

        # Check phylogeny tab is NOT present
        assert "phylogeny-container" not in html

    def test_report_without_ani(self, sample_data, tmp_path):
        """Report without ANI matrix has no phylogeny tab."""
        output = tmp_path / "report_no_ani.html"

        result = runner.invoke(app, [
            "generate",
            "--classifications", str(sample_data.classifications),
            "--output", str(output),
        ])

        assert result.exit_code == 0
        html = output.read_text()

        # Check phylogeny tab is NOT present (no ANI = no tree)
        assert "phylogeny-container" not in html

    def test_report_with_user_tree(self, sample_data, tmp_path):
        """Report uses user-provided tree file."""
        tree_file = tmp_path / "tree.nwk"
        tree_file.write_text("((GCF_000001.1:0.1,GCF_000002.1:0.2):0.3,GCF_000003.1:0.4);")

        result = runner.invoke(app, [
            "generate",
            "--classifications", str(sample_data.classifications),
            "--ani", str(sample_data.ani),
            "--tree", str(tree_file),
            "--output", str(sample_data.output),
        ])

        assert result.exit_code == 0
        html = sample_data.output.read_text()

        # Check tree was used
        assert "phylogeny-container" in html
        assert "user-provided" in html
```

**Step 2: Run integration tests**

Run: `pytest tests/integration/test_report_phylogeny.py -v`
Expected: All PASS

**Step 3: Commit**

```bash
git add tests/integration/test_report_phylogeny.py
git commit -m "test: add integration tests for phylogeny report generation"
```

---

## Task 7: Final verification and documentation

**Files:**
- Verify all tests pass
- Update CLAUDE.md if needed

**Step 1: Run full test suite**

Run: `pytest tests/ -v --ignore=tests/integration/test_mmseqs2_integration.py --ignore=tests/e2e/`
Expected: All tests PASS

**Step 2: Manual verification**

Generate a test report with phylogeny:

```bash
# Create test data if needed, then:
metadarkmatter report generate \
    --classifications test_classifications.csv \
    --ani test_ani_matrix.csv \
    --output test_report.html
```

Verify in browser:
- [ ] Phylogeny tab appears
- [ ] Tree renders correctly
- [ ] Novel clusters are colored and labeled
- [ ] Hover tooltips work
- [ ] Legend is clear

**Step 3: Final commit**

```bash
git add -A
git commit -m "feat(phylogeny): complete phylogenetic tree visualization feature

Adds interactive phylogenetic tree tab to HTML reports:
- Auto-generates neighbor-joining tree from ANI matrix
- Accepts user-provided Newick trees via --tree option
- Places novel clusters as leaf nodes at estimated positions
- Color-coded by classification type (species/genus)
- Shows confidence badges and hover tooltips
- Supports --no-phylogeny flag to skip for performance

Closes design: docs/plans/2026-02-05-phylogeny-visualization-design.md"
```

---

## Summary

| Task | Description | Files Changed |
|------|-------------|---------------|
| 1 | Tree builder module | `core/phylogeny/tree_builder.py`, tests |
| 2 | Novel cluster placement | `core/phylogeny/placement.py`, tests |
| 3 | Templates and styles | `templates.py`, `styles.py` |
| 4 | Report generator integration | `generator.py`, tests |
| 5 | CLI options | `cli/report.py`, tests |
| 6 | Integration tests | `test_report_phylogeny.py` |
| 7 | Final verification | Documentation, manual testing |

**Total estimated tasks:** 7 main tasks with ~25 individual steps

**Dependencies:** BioPython (already installed), D3.js (loaded from CDN)
