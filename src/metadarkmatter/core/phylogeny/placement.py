"""Place novel clusters on phylogenetic trees.

This module provides functionality to insert novel taxon clusters as leaf nodes
on phylogenetic trees at estimated positions based on their divergence from
reference genomes. Novel species are placed as siblings to their nearest
reference, while novel genera are placed at greater distances.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from io import StringIO
from typing import TYPE_CHECKING, Any, Literal

import polars as pl

if TYPE_CHECKING:
    import pandas as pd
    from Bio.Phylo.BaseTree import Tree

logger = logging.getLogger(__name__)


@dataclass
class NovelCluster:
    """A cluster of novel reads for phylogenetic tree placement.

    This is a lightweight dataclass specifically for tree placement operations.
    For the full-featured cluster model with GTDB-style naming and phylogenetic
    context, see metadarkmatter.core.novel_diversity.models.NovelCluster.

    Attributes:
        cluster_id: Unique identifier (e.g., NSP_001 for Novel Species).
        classification: Either "novel_species" or "novel_genus".
        best_match_genome: Accession of the nearest reference genome.
        mean_identity: Mean percent identity of reads in the cluster.
        mean_novelty: Mean novelty index (100 - identity).
        mean_uncertainty: Mean placement uncertainty of reads.
        read_count: Number of reads in this cluster.
        confidence_rating: Quality rating ("High", "Medium", or "Low").
    """

    cluster_id: str
    classification: Literal["novel_species", "novel_genus"]
    best_match_genome: str
    mean_identity: float
    mean_novelty: float
    mean_uncertainty: float
    read_count: int
    confidence_rating: Literal["High", "Medium", "Low"]

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for serialization.

        Returns:
            Dictionary with all cluster attributes.
        """
        return {
            "cluster_id": self.cluster_id,
            "classification": self.classification,
            "best_match_genome": self.best_match_genome,
            "mean_identity": self.mean_identity,
            "mean_novelty": self.mean_novelty,
            "mean_uncertainty": self.mean_uncertainty,
            "read_count": self.read_count,
            "confidence_rating": self.confidence_rating,
        }


def _estimate_branch_length(novelty: float, classification: str) -> float:
    """Estimate branch length based on novelty level.

    Uses a simple linear mapping from novelty index to branch length.
    Novel genus clusters get longer branches than novel species.

    Args:
        novelty: Mean novelty index (0-100).
        classification: "novel_species" or "novel_genus".

    Returns:
        Estimated branch length for tree insertion.
    """
    # Base branch length scales with novelty
    # Novelty of 5% -> ~0.05, Novelty of 25% -> ~0.25
    base_length = novelty / 100.0

    # Novel genus gets additional length to reflect greater divergence
    if classification == "novel_genus":
        base_length *= 1.5

    return round(base_length, 4)


def _find_node_by_name(tree: Tree, name: str) -> Any:
    """Find a node in the tree by its name.

    Args:
        tree: BioPython Tree object.
        name: Name of the node to find.

    Returns:
        The node with the matching name, or None if not found.
    """
    for node in tree.find_clades():
        if node.name == name:
            return node
    return None


def _create_sibling_to_node(
    tree: Tree,
    target_name: str,
    new_name: str,
    new_branch_length: float,
    comment: str | None = None,
) -> bool:
    """Create a new node as a sibling to an existing node.

    Inserts a new leaf node as a sibling to the target node by:
    1. Finding the parent of the target node
    2. Creating a new internal node between parent and target
    3. Adding the new leaf as a sibling to target under the new internal node

    Args:
        tree: BioPython Tree object to modify.
        target_name: Name of the node to become sibling to.
        new_name: Name for the new node.
        new_branch_length: Branch length for the new node.
        comment: Optional JSON comment to attach to the new node.

    Returns:
        True if insertion succeeded, False otherwise.
    """
    from Bio.Phylo.BaseTree import Clade

    target = _find_node_by_name(tree, target_name)
    if target is None:
        return False

    # Find parent of target
    parent = None
    for clade in tree.find_clades():
        if target in clade.clades:
            parent = clade
            break

    if parent is None:
        # Target might be at root level
        logger.warning(f"Could not find parent for {target_name}")
        return False

    # Store target's current branch length
    target_original_length = target.branch_length or 0.1

    # Remove target from parent temporarily
    parent.clades.remove(target)

    # Create new internal node
    # The internal node gets half the original branch length
    internal_length = target_original_length / 2
    internal_node = Clade(branch_length=internal_length)

    # Target becomes child of internal node with remaining length
    target.branch_length = target_original_length / 2

    # Create new leaf node
    new_leaf = Clade(
        name=new_name,
        branch_length=new_branch_length,
    )
    if comment:
        new_leaf.comment = comment

    # Add both as children of the new internal node
    internal_node.clades = [target, new_leaf]

    # Add internal node to parent
    parent.clades.append(internal_node)

    return True


def place_novel_clusters(
    tree: str,
    novel_clusters: list[NovelCluster],
    ani_matrix: pd.DataFrame,
    genome_metadata: dict[str, dict[str, str]] | None = None,
) -> str:
    """Place novel clusters on a phylogenetic tree.

    Inserts novel cluster nodes at estimated positions based on their
    relationship to reference genomes. Novel species are placed as siblings
    to their nearest reference, while novel genera are placed at greater
    distances reflecting their divergence.

    Args:
        tree: Newick format tree string.
        novel_clusters: List of NovelCluster objects to place.
        ani_matrix: ANI matrix DataFrame for reference (not currently used
            for placement but available for future distance-based placement).
        genome_metadata: Optional dictionary mapping genome accessions to
            metadata dicts with species/genus information.

    Returns:
        Modified Newick tree string with novel clusters inserted.

    Note:
        Clusters whose best_match_genome is not found in the tree are skipped
        with a warning logged.
    """
    from Bio import Phylo

    if not novel_clusters:
        # Return tree unchanged
        return tree

    # Parse the tree
    tree_obj = Phylo.read(StringIO(tree), "newick")

    # Get existing tip names
    existing_tips = {t.name for t in tree_obj.get_terminals()}

    # Place each cluster
    for cluster in novel_clusters:
        if cluster.best_match_genome not in existing_tips:
            logger.warning(
                f"Genome {cluster.best_match_genome} not found in tree; "
                f"skipping placement of cluster {cluster.cluster_id}"
            )
            continue

        # Estimate branch length based on novelty
        branch_length = _estimate_branch_length(
            cluster.mean_novelty,
            cluster.classification,
        )

        # Create metadata comment as JSON
        comment_data: dict[str, str | int | float] = {
            "classification": cluster.classification,
            "read_count": cluster.read_count,
            "mean_novelty": cluster.mean_novelty,
            "mean_uncertainty": cluster.mean_uncertainty,
            "confidence": cluster.confidence_rating,
        }

        # Add genome metadata if available
        if genome_metadata and cluster.best_match_genome in genome_metadata:
            meta = genome_metadata[cluster.best_match_genome]
            comment_data["nearest_species"] = meta.get("species", "Unknown")
            comment_data["nearest_genus"] = meta.get("genus", "Unknown")

        comment = json.dumps(comment_data)

        # Insert as sibling to best match genome
        success = _create_sibling_to_node(
            tree_obj,
            cluster.best_match_genome,
            cluster.cluster_id,
            branch_length,
            comment,
        )

        if success:
            logger.debug(f"Placed {cluster.cluster_id} as sibling to {cluster.best_match_genome}")
        else:
            logger.warning(f"Failed to place cluster {cluster.cluster_id}")

    # Convert back to Newick
    output = StringIO()
    Phylo.write(tree_obj, output, "newick")
    return output.getvalue().strip()


def _assign_confidence_rating(
    read_count: int,
    mean_uncertainty: float,
) -> Literal["High", "Medium", "Low"]:
    """Assign a confidence rating to a cluster based on its properties.

    Criteria:
        High: read_count >= 10 and mean_uncertainty < 2%
        Medium: read_count >= 5 and mean_uncertainty < 5%
        Low: All other clusters

    Args:
        read_count: Number of reads in the cluster.
        mean_uncertainty: Mean placement uncertainty.

    Returns:
        Confidence rating string.
    """
    if read_count >= 10 and mean_uncertainty < 2.0:
        return "High"
    elif read_count >= 5 and mean_uncertainty < 5.0:
        return "Medium"
    return "Low"


def extract_novel_clusters(
    classifications: pl.DataFrame,
    min_reads: int = 3,
) -> list[NovelCluster]:
    """Extract novel clusters from a classification DataFrame.

    Groups novel reads (Novel Species and Novel Genus) by their best match
    genome and taxonomic call, then creates NovelCluster objects for groups
    meeting the minimum read threshold.

    Args:
        classifications: Polars DataFrame with classification results.
            Required columns: read_id, best_match_genome, top_hit_identity,
            novelty_index, placement_uncertainty, taxonomic_call.
        min_reads: Minimum reads required to form a cluster (default 3).

    Returns:
        List of NovelCluster objects sorted by read count (descending).

    Example:
        >>> clusters = extract_novel_clusters(classifications_df, min_reads=5)
        >>> for c in clusters:
        ...     print(f"{c.cluster_id}: {c.read_count} reads")
    """
    # Filter to novel reads only
    novel_df = classifications.filter(
        pl.col("taxonomic_call").is_in(["Novel Species", "Novel Genus"])
    )

    if len(novel_df) == 0:
        return []

    # Group by best_match_genome and taxonomic_call
    grouped = (
        novel_df
        .group_by(["best_match_genome", "taxonomic_call"])
        .agg([
            pl.len().alias("read_count"),
            pl.col("top_hit_identity").mean().alias("mean_identity"),
            pl.col("novelty_index").mean().alias("mean_novelty"),
            pl.col("placement_uncertainty").mean().alias("mean_uncertainty"),
        ])
    )

    # Filter by minimum cluster size
    grouped = grouped.filter(pl.col("read_count") >= min_reads)

    if len(grouped) == 0:
        return []

    # Sort by read count descending
    grouped = grouped.sort("read_count", descending=True)

    # Build NovelCluster objects
    clusters: list[NovelCluster] = []
    nsp_counter = 0
    ngn_counter = 0

    for row in grouped.iter_rows(named=True):
        tax_call = row["taxonomic_call"]

        # Generate cluster ID
        if tax_call == "Novel Species":
            nsp_counter += 1
            cluster_id = f"NSP_{nsp_counter:03d}"
            classification: Literal["novel_species", "novel_genus"] = "novel_species"
        else:
            ngn_counter += 1
            cluster_id = f"NGN_{ngn_counter:03d}"
            classification = "novel_genus"

        # Assign confidence rating
        confidence = _assign_confidence_rating(
            row["read_count"],
            row["mean_uncertainty"],
        )

        cluster = NovelCluster(
            cluster_id=cluster_id,
            classification=classification,
            best_match_genome=row["best_match_genome"],
            mean_identity=round(row["mean_identity"], 2),
            mean_novelty=round(row["mean_novelty"], 2),
            mean_uncertainty=round(row["mean_uncertainty"], 2),
            read_count=row["read_count"],
            confidence_rating=confidence,
        )
        clusters.append(cluster)

    logger.info(f"Extracted {len(clusters)} novel clusters from classifications")
    return clusters
