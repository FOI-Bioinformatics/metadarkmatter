"""
Extended similarity matrix builder for novel diversity visualization.

This module provides functions for building ANI/AAI matrices that include
both reference genomes and novel clusters, enabling phylogenetic context
visualization where novel diversity can be seen alongside known species.

Supports both:
- ANI (Average Nucleotide Identity) for species-level comparisons
- AAI (Average Amino Acid Identity) for genus-level comparisons
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import polars as pl
    from metadarkmatter.core.novel_diversity import NovelCluster

logger = logging.getLogger(__name__)


def build_extended_similarity_matrix(
    similarity_matrix: pl.DataFrame,
    novel_clusters: list[NovelCluster],
    genome_labels_map: dict[str, str],
    similarity_type: str = "ANI",
    default_value: float | None = None,
    max_references: int = 50,
    max_clusters: int = 20,
) -> tuple[np.ndarray, list[str], list[bool]]:
    """
    Build extended similarity matrix (ANI or AAI) including novel clusters.

    Creates a combined matrix with reference genomes and novel cluster
    representatives, estimating similarity values between novel clusters
    and references based on classification data.

    Args:
        similarity_matrix: Similarity matrix DataFrame (wide format with 'genome' column)
        novel_clusters: List of NovelCluster objects to include
        genome_labels_map: Mapping from accession to display label
        similarity_type: "ANI" for nucleotide or "AAI" for protein comparisons
        default_value: Default value for distant/unknown pairs
            (defaults to 70 for ANI, 40 for AAI)
        max_references: Maximum number of reference genomes to include
        max_clusters: Maximum number of novel clusters to include

    Returns:
        Tuple of (extended_matrix, labels, is_novel_mask)
        - extended_matrix: NxN numpy array with similarity values
        - labels: Ordered list of labels for rows/columns
        - is_novel_mask: Boolean list indicating novel cluster entries
    """
    import polars as pl

    # Set default based on similarity type
    is_aai = similarity_type.upper() == "AAI"
    if default_value is None:
        default_value = 40.0 if is_aai else 70.0

    # Extract genome accessions from similarity matrix
    matrix_cols = similarity_matrix.columns
    if "genome" in matrix_cols:
        all_genomes = similarity_matrix["genome"].to_list()
        z_data = similarity_matrix.drop("genome").to_numpy()
    else:
        all_genomes = list(matrix_cols)
        z_data = similarity_matrix.to_numpy()

    # Build similarity lookup dictionary for fast access
    sim_dict = _build_similarity_dict(all_genomes, z_data, default_value)

    # Select relevant references and clusters
    selected_refs = select_relevant_references(
        sim_dict, novel_clusters, all_genomes, max_references
    )
    selected_clusters = novel_clusters[:max_clusters]

    # Build combined entity list
    n_refs = len(selected_refs)
    n_clusters = len(selected_clusters)
    n_total = n_refs + n_clusters

    if n_total == 0:
        return np.array([[]]), [], []

    # Initialize extended matrix
    extended = np.full((n_total, n_total), default_value, dtype=np.float64)

    # Fill ref-to-ref values from similarity matrix
    for i, ref_i in enumerate(selected_refs):
        for j, ref_j in enumerate(selected_refs):
            if i == j:
                extended[i, j] = 100.0
            else:
                extended[i, j] = sim_dict.get((ref_i, ref_j), default_value)

    # Fill ref-to-novel and novel-to-ref
    for c_idx, cluster in enumerate(selected_clusters):
        matrix_idx = n_refs + c_idx

        # Self-similarity for novel cluster (identity)
        extended[matrix_idx, matrix_idx] = 100.0

        for r_idx, ref_genome in enumerate(selected_refs):
            sim_val = estimate_novel_to_reference_similarity(
                cluster, ref_genome, sim_dict, default_value, is_aai
            )
            # Symmetric matrix
            extended[r_idx, matrix_idx] = sim_val
            extended[matrix_idx, r_idx] = sim_val

    # Fill novel-to-novel
    for i, c1 in enumerate(selected_clusters):
        for j, c2 in enumerate(selected_clusters):
            if i == j:
                continue
            matrix_i = n_refs + i
            matrix_j = n_refs + j
            sim_val = estimate_novel_to_novel_similarity(
                c1, c2, sim_dict, default_value, is_aai
            )
            extended[matrix_i, matrix_j] = sim_val
            extended[matrix_j, matrix_i] = sim_val

    # Build labels list
    labels = []
    is_novel_mask = []

    for ref in selected_refs:
        label = genome_labels_map.get(ref, ref)
        labels.append(label)
        is_novel_mask.append(False)

    for cluster in selected_clusters:
        # Novel cluster label format: [*] NSP_001: Name (N reads)
        label = f"[*] {cluster.cluster_id}: {cluster.suggested_name} ({cluster.read_count} reads)"
        labels.append(label)
        is_novel_mask.append(True)

    logger.debug(
        f"Built extended {similarity_type} matrix: {n_refs} references + {n_clusters} novel clusters"
    )

    return extended, labels, is_novel_mask


def _build_similarity_dict(
    genomes: list[str],
    matrix: np.ndarray,
    default_value: float,
) -> dict[tuple[str, str], float]:
    """
    Build a dictionary for fast similarity lookups.

    Args:
        genomes: List of genome accessions
        matrix: Similarity matrix as numpy array
        default_value: Default value for missing entries

    Returns:
        Dictionary mapping (genome1, genome2) -> similarity value
    """
    sim_dict: dict[tuple[str, str], float] = {}

    for i, g1 in enumerate(genomes):
        for j, g2 in enumerate(genomes):
            val = matrix[i, j]
            # Skip zeros/default values, store actual similarity values
            if val > default_value:
                sim_dict[(g1, g2)] = float(val)

    return sim_dict


def estimate_novel_to_reference_similarity(
    cluster: NovelCluster,
    ref_genome: str,
    sim_dict: dict[tuple[str, str], float],
    default_value: float = 70.0,
    is_aai: bool = False,
) -> float:
    """
    Estimate similarity (ANI or AAI) between a novel cluster and a reference genome.

    Uses the cluster's estimated similarity for the nearest reference, and
    triangular estimation for other references based on their similarity
    to the nearest genome.

    The cluster's estimated_ani (100 - novelty_index) represents the similarity
    to the nearest reference based on whatever alignment mode was used:
    - Nucleotide mode (BLASTN): estimated_ani reflects nucleotide identity
    - Protein mode (BLASTX): estimated_ani reflects protein identity (AAI-like)

    Args:
        cluster: NovelCluster object
        ref_genome: Reference genome accession
        sim_dict: Pre-computed similarity lookup dictionary
        default_value: Default value for unknown pairs
        is_aai: True if using AAI mode (protein), False for ANI (nucleotide)

    Returns:
        Estimated similarity value (0-100)
    """
    nearest = cluster.nearest_genome

    # Use estimated_ani directly - it represents 100 - novelty_index
    # which is the similarity to the nearest reference from the BLAST alignment.
    # This value already reflects the alignment mode used (nucleotide or protein).
    estimated_to_nearest = cluster.estimated_ani

    if ref_genome == nearest:
        # Direct relationship: use cluster's estimated similarity
        return estimated_to_nearest

    # Triangular estimation via nearest reference
    # Novel -> Nearest: estimated_to_nearest
    # Nearest -> Other Ref: lookup in similarity matrix
    nearest_to_ref = sim_dict.get((nearest, ref_genome), default_value)

    # Estimate: combine the divergences with appropriate scaling
    divergence_from_nearest = 100.0 - estimated_to_nearest
    divergence_nearest_to_ref = 100.0 - nearest_to_ref

    # Combined divergence with diminishing returns
    combined_divergence = divergence_from_nearest + (divergence_nearest_to_ref * 0.7)
    estimated = 100.0 - combined_divergence

    # Clamp to valid range
    return max(default_value, min(100.0, estimated))


def estimate_novel_to_novel_similarity(
    cluster1: NovelCluster,
    cluster2: NovelCluster,
    sim_dict: dict[tuple[str, str], float],
    default_value: float = 70.0,
    is_aai: bool = False,
) -> float:
    """
    Estimate similarity (ANI or AAI) between two novel clusters.

    Uses the relationship through their nearest reference genomes
    to estimate similarity.

    Args:
        cluster1: First NovelCluster
        cluster2: Second NovelCluster
        sim_dict: Pre-computed similarity lookup dictionary
        default_value: Default value for unknown pairs
        is_aai: True if using AAI mode (protein), False for ANI (nucleotide)

    Returns:
        Estimated similarity value (0-100)
    """
    nearest1 = cluster1.nearest_genome
    nearest2 = cluster2.nearest_genome

    # Use estimated_ani directly - it reflects the alignment mode used
    est1 = cluster1.estimated_ani
    est2 = cluster2.estimated_ani

    if nearest1 == nearest2:
        # Same nearest reference: clusters are in similar phylogenetic space
        avg_sim = (est1 + est2) / 2
        # Subtract small penalty for being distinct clusters
        return max(default_value, avg_sim - 2.0)

    # Different nearest references: go through similarity matrix
    ref_to_ref_sim = sim_dict.get((nearest1, nearest2), default_value)

    # Calculate total divergence through the path
    div1 = 100.0 - est1
    div2 = 100.0 - est2
    div_refs = 100.0 - ref_to_ref_sim

    # Combined divergence (not simply additive)
    combined = div_refs + (div1 + div2) * 0.5
    estimated = 100.0 - combined

    return max(default_value, min(100.0, estimated))


def select_relevant_references(
    ani_dict: dict[tuple[str, str], float],
    clusters: list[NovelCluster],
    all_genomes: list[str],
    max_refs: int = 50,
) -> list[str]:
    """
    Select the most relevant reference genomes for visualization.

    Prioritizes:
    1. All nearest genomes from novel clusters
    2. Genomes with high ANI to nearest genomes (phylogenetic context)
    3. Diverse representation across the ANI matrix

    Args:
        ani_dict: Pre-computed ANI lookup dictionary
        clusters: List of NovelCluster objects
        all_genomes: All available genome accessions
        max_refs: Maximum number of references to select

    Returns:
        List of selected genome accessions
    """
    if not clusters or not all_genomes:
        return all_genomes[:max_refs]

    selected: set[str] = set()

    # Priority 1: Include all nearest genomes from clusters
    for cluster in clusters:
        if cluster.nearest_genome in all_genomes:
            selected.add(cluster.nearest_genome)

    # Priority 2: Add genomes with high ANI to the nearest genomes
    nearest_genomes = list(selected)
    for nearest in nearest_genomes:
        # Find genomes similar to this nearest reference
        similar = []
        for genome in all_genomes:
            if genome in selected:
                continue
            ani = ani_dict.get((nearest, genome), 0.0)
            if ani > 80.0:  # Include if reasonably similar
                similar.append((genome, ani))

        # Sort by ANI and add top candidates
        similar.sort(key=lambda x: x[1], reverse=True)
        for genome, _ in similar[:5]:  # Add up to 5 per nearest
            selected.add(genome)
            if len(selected) >= max_refs:
                break

        if len(selected) >= max_refs:
            break

    # Priority 3: Fill remaining slots with diverse genomes
    # Select genomes that maximize coverage of the ANI space
    if len(selected) < max_refs:
        remaining = [g for g in all_genomes if g not in selected]

        # Add remaining genomes to reach max_refs
        for genome in remaining:
            selected.add(genome)
            if len(selected) >= max_refs:
                break

    # Return in original order for consistency
    return [g for g in all_genomes if g in selected]
