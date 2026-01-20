"""
Hierarchical clustering utilities for ANI/AAI heatmap visualization.

This module provides hierarchical clustering functionality for reordering
genome similarity matrices to create visually grouped heatmaps.
"""

from __future__ import annotations

import logging

import numpy as np

logger = logging.getLogger(__name__)


def perform_hierarchical_clustering(
    matrix: np.ndarray,
    labels: list[str],
    default_value: float,
) -> tuple[np.ndarray, list[str], bool]:
    """
    Perform hierarchical clustering on a similarity matrix (ANI or AAI).

    Uses scipy hierarchical clustering with average linkage (UPGMA) to reorder
    genomes by similarity, creating visual blocks of related taxa in heatmaps.

    Args:
        matrix: Square similarity matrix (ANI or AAI values)
        labels: Genome labels corresponding to matrix rows/columns
        default_value: Fill value for missing data (e.g., 70 for ANI, 40 for AAI)

    Returns:
        Tuple of (clustered_matrix, ordered_labels, clustering_succeeded)
        - clustered_matrix: Matrix reordered by clustering
        - ordered_labels: Labels reordered to match clustered matrix
        - clustering_succeeded: True if scipy available and clustering worked

    Note:
        If scipy is unavailable, returns original matrix/labels with success=False.
        Logs warning when scipy import fails.
    """
    try:
        from scipy.cluster.hierarchy import leaves_list, linkage
        from scipy.spatial.distance import squareform

        # Fill diagonal (self-comparison should be 100%)
        np.fill_diagonal(matrix, 100.0)

        # Convert similarity to distance (100 - similarity)
        dist_matrix = 100.0 - matrix

        # Make symmetric (average upper and lower triangles)
        dist_matrix = (dist_matrix + dist_matrix.T) / 2
        np.fill_diagonal(dist_matrix, 0.0)

        # Perform hierarchical clustering with average linkage (UPGMA)
        condensed_dist = squareform(dist_matrix)
        linkage_matrix = linkage(condensed_dist, method="average")

        # Get optimal leaf ordering
        order = leaves_list(linkage_matrix)

        # Reorder matrix and labels
        clustered_matrix = matrix[order][:, order]
        ordered_labels = [labels[i] for i in order]

        logger.debug(
            f"Hierarchical clustering successful ({len(labels)} genomes)"
        )
        return (clustered_matrix, ordered_labels, True)

    except ImportError:
        logger.warning(
            "scipy not available - hierarchical clustering skipped. "
            "Install scipy for clustered heatmaps: pip install scipy>=1.11.0"
        )
        return (matrix, labels, False)
