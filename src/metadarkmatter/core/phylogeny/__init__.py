"""Phylogeny module for tree building and novel cluster placement.

Provides functions to convert ANI matrices to phylogenetic trees using
neighbor-joining, load and validate user-provided Newick trees, and place
novel taxon clusters on trees at estimated phylogenetic positions.
"""

from metadarkmatter.core.phylogeny.placement import (
    NovelCluster,
    extract_novel_clusters,
    place_novel_clusters,
)
from metadarkmatter.core.phylogeny.tree_builder import (
    ani_to_newick,
    load_user_tree,
)

__all__ = [
    "ani_to_newick",
    "load_user_tree",
    "NovelCluster",
    "extract_novel_clusters",
    "place_novel_clusters",
]
