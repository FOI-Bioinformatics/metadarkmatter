"""
Heatmap generation utilities for ANI/AAI similarity matrices.

This module provides functions for building interactive heatmaps with
statistics cards for genome similarity matrices (ANI and AAI).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import plotly.graph_objects as go

from metadarkmatter.visualization.report.components.clustering import (
    perform_hierarchical_clustering,
)

if TYPE_CHECKING:
    import polars as pl


@dataclass
class SimilarityStats:
    """Statistics computed from a similarity matrix."""

    num_genomes: int
    min_value: float
    max_value: float
    mean_value: float
    same_category_count: int  # e.g., same species for ANI
    boundary_count: int
    different_category_count: int


def compute_similarity_stats(
    matrix: np.ndarray,
    default_value: float,
    same_threshold: float,
    boundary_lower: float,
) -> SimilarityStats:
    """
    Compute statistics from a similarity matrix.

    Args:
        matrix: Square similarity matrix (clustered)
        default_value: Default value for missing data
        same_threshold: Threshold for "same category" (e.g., 95 for species, 65 for genus)
        boundary_lower: Lower boundary for boundary zone

    Returns:
        SimilarityStats with computed values
    """
    # Get non-diagonal values
    non_diag = matrix[~np.eye(matrix.shape[0], dtype=bool)]
    non_zero = non_diag[non_diag > default_value]

    if len(non_zero) > 0:
        min_val = float(np.min(non_zero))
        max_val = float(np.max(non_zero))
        mean_val = float(np.mean(non_zero))
    else:
        min_val = max_val = mean_val = 0.0

    # Category counts
    same_count = int(np.sum(non_zero >= same_threshold))
    diff_count = int(np.sum(non_zero < boundary_lower))
    boundary_count = int(np.sum((non_zero >= boundary_lower) & (non_zero < same_threshold)))

    return SimilarityStats(
        num_genomes=matrix.shape[0],
        min_value=min_val,
        max_value=max_val,
        mean_value=mean_val,
        same_category_count=same_count,
        boundary_count=boundary_count,
        different_category_count=diff_count,
    )


def build_ani_heatmap(
    ani_matrix: pl.DataFrame,
    genome_labels_map: dict[str, str],
    default_ani: float = 70.0,
) -> tuple[go.Figure, SimilarityStats, bool]:
    """
    Build ANI heatmap with hierarchical clustering.

    Args:
        ani_matrix: ANI matrix DataFrame (wide format)
        genome_labels_map: Mapping from accession to display label
        default_ani: Default ANI value for missing data

    Returns:
        Tuple of (plotly_figure, statistics, clustering_succeeded)
    """
    # Extract genome names and matrix data
    matrix_cols = ani_matrix.columns
    if "genome" in matrix_cols:
        genome_accessions = ani_matrix["genome"].to_list()
        z_data = ani_matrix.drop("genome").to_numpy()
    else:
        genome_accessions = list(matrix_cols)
        z_data = ani_matrix.to_numpy()

    # Fill missing values
    z_filled = np.where(z_data == 0.0, default_ani, z_data)

    # Hierarchical clustering
    z_clustered, genome_accessions_ordered, clustering_succeeded = (
        perform_hierarchical_clustering(z_filled, genome_accessions, default_ani)
    )

    # Create labels
    genome_labels = [
        genome_labels_map.get(acc, acc) for acc in genome_accessions_ordered
    ]

    # Compute statistics
    stats = compute_similarity_stats(
        z_clustered,
        default_value=default_ani,
        same_threshold=95.0,  # Species boundary
        boundary_lower=93.0,
    )

    # Convert to Python list for JSON serialization
    z_list = [[float(v) for v in row] for row in z_clustered]

    # Inverted colorscale: red = high ANI (similar), blue = low ANI (distant)
    ani_colorscale = [
        [0.0, "#4575b4"],    # 70% - Blue (distant/different genera)
        [0.167, "#74add1"],  # 75%
        [0.333, "#abd9e9"],  # 80%
        [0.5, "#ffffbf"],    # 85% - Yellow (genus boundary)
        [0.667, "#fee090"],  # 90%
        [0.833, "#fdae61"],  # 95% - Orange (species boundary)
        [1.0, "#d73027"],    # 100% - Red (identical/same species)
    ]

    # Create heatmap
    fig = go.Figure(
        data=go.Heatmap(
            z=z_list,
            x=genome_labels,
            y=genome_labels,
            colorscale=ani_colorscale,
            zmin=70,
            zmax=100,
            colorbar=dict(
                title="ANI (%)",
                tickvals=[70, 75, 80, 85, 90, 95, 100],
                ticktext=["70", "75", "80", "85 (genus)", "90", "95 (species)", "100"],
            ),
            hovertemplate=(
                "Genome 1: %{y}<br>"
                "Genome 2: %{x}<br>"
                "ANI: %{z:.1f}%<extra></extra>"
            ),
        )
    )

    # Conditional title
    title = "ANI Matrix Heatmap"
    if clustering_succeeded:
        title += " (Hierarchically Clustered)"

    fig.update_layout(
        title=title,
        width=1000,
        height=900,
        template="plotly_white",
        xaxis=dict(tickangle=45, tickfont=dict(size=8)),
        yaxis=dict(tickfont=dict(size=8)),
    )

    return fig, stats, clustering_succeeded


def build_aai_heatmap(
    aai_matrix: pl.DataFrame,
    genome_labels_map: dict[str, str],
    default_aai: float = 40.0,
) -> tuple[go.Figure, SimilarityStats, bool]:
    """
    Build AAI heatmap with hierarchical clustering.

    Args:
        aai_matrix: AAI matrix DataFrame (wide format)
        genome_labels_map: Mapping from accession to display label
        default_aai: Default AAI value for missing data

    Returns:
        Tuple of (plotly_figure, statistics, clustering_succeeded)
    """
    # Extract genome names and matrix data
    matrix_cols = aai_matrix.columns
    if "genome" in matrix_cols:
        genome_accessions = aai_matrix["genome"].to_list()
        z_data = aai_matrix.drop("genome").to_numpy()
    else:
        genome_accessions = list(matrix_cols)
        z_data = aai_matrix.to_numpy()

    # Fill missing values
    z_filled = np.where(z_data == 0.0, default_aai, z_data)

    # Hierarchical clustering
    z_clustered, genome_accessions_ordered, clustering_succeeded = (
        perform_hierarchical_clustering(z_filled, genome_accessions, default_aai)
    )

    # Create labels
    genome_labels = [
        genome_labels_map.get(acc, acc) for acc in genome_accessions_ordered
    ]

    # Compute statistics
    stats = compute_similarity_stats(
        z_clustered,
        default_value=default_aai,
        same_threshold=65.0,  # Genus boundary
        boundary_lower=58.0,
    )

    # Convert to Python list for JSON serialization
    z_list = [[float(v) for v in row] for row in z_clustered]

    # AAI colorscale: 40-100% range for protein-level comparisons
    # Matches the phylogenetic context heatmap scale
    aai_colorscale = [
        [0.0, "#313695"],     # 40% - Dark blue (very distant)
        [0.167, "#4575b4"],   # 50%
        [0.333, "#abd9e9"],   # 60%
        [0.417, "#ffffbf"],   # 65% - Yellow (genus boundary)
        [0.5, "#fee090"],     # 70%
        [0.583, "#fdae61"],   # 75%
        [0.75, "#f46d43"],    # 85%
        [1.0, "#a50026"],     # 100% - Dark red (identical)
    ]

    # Create heatmap
    fig = go.Figure(
        data=go.Heatmap(
            z=z_list,
            x=genome_labels,
            y=genome_labels,
            colorscale=aai_colorscale,
            zmin=40,
            zmax=100,
            colorbar=dict(
                title="AAI (%)",
                tickvals=[40, 50, 60, 65, 70, 80, 90, 100],
                ticktext=["40", "50", "60", "65 (genus)", "70", "80", "90", "100"],
            ),
            hovertemplate=(
                "Genome 1: %{y}<br>"
                "Genome 2: %{x}<br>"
                "AAI: %{z:.1f}%<extra></extra>"
            ),
        )
    )

    # Conditional title
    title = "AAI Matrix Heatmap"
    if clustering_succeeded:
        title += " (Hierarchically Clustered)"

    fig.update_layout(
        title=title,
        width=1000,
        height=900,
        template="plotly_white",
        xaxis=dict(tickangle=45, tickfont=dict(size=8)),
        yaxis=dict(tickfont=dict(size=8)),
    )

    return fig, stats, clustering_succeeded


def build_ani_stats_cards(stats: SimilarityStats) -> str:
    """
    Generate HTML for ANI statistics cards.

    Args:
        stats: Computed similarity statistics

    Returns:
        HTML string with metric cards
    """
    return f"""
    <div class="metric-cards">
        <div class="metric-card">
            <div class="metric-value">{stats.num_genomes}</div>
            <div class="metric-label">Genomes</div>
            <div class="metric-subtext">In ANI matrix</div>
        </div>
        <div class="metric-card">
            <div class="metric-value">{stats.mean_value:.1f}%</div>
            <div class="metric-label">Mean ANI</div>
            <div class="metric-subtext">Range: {stats.min_value:.1f}% - {stats.max_value:.1f}%</div>
        </div>
        <div class="metric-card success">
            <div class="metric-value">{stats.same_category_count}</div>
            <div class="metric-label">Genome Pairs</div>
            <div class="metric-subtext">Same species (ANI >= 95%)</div>
        </div>
        <div class="metric-card warning">
            <div class="metric-value">{stats.boundary_count}</div>
            <div class="metric-label">Genome Pairs</div>
            <div class="metric-subtext">Boundary zone (ANI 93-95%)</div>
        </div>
        <div class="metric-card danger">
            <div class="metric-value">{stats.different_category_count}</div>
            <div class="metric-label">Genome Pairs</div>
            <div class="metric-subtext">Different species (ANI < 93%)</div>
        </div>
    </div>
    """


def build_aai_stats_cards(stats: SimilarityStats) -> str:
    """
    Generate HTML for AAI statistics cards.

    Args:
        stats: Computed similarity statistics

    Returns:
        HTML string with metric cards
    """
    return f"""
    <div class="metric-cards">
        <div class="metric-card">
            <div class="metric-value">{stats.num_genomes}</div>
            <div class="metric-label">Genomes</div>
            <div class="metric-subtext">In AAI matrix</div>
        </div>
        <div class="metric-card">
            <div class="metric-value">{stats.mean_value:.1f}%</div>
            <div class="metric-label">Mean AAI</div>
            <div class="metric-subtext">Range: {stats.min_value:.1f}% - {stats.max_value:.1f}%</div>
        </div>
        <div class="metric-card success">
            <div class="metric-value">{stats.same_category_count}</div>
            <div class="metric-label">Genome Pairs</div>
            <div class="metric-subtext">Same genus (AAI >= 65%)</div>
        </div>
        <div class="metric-card warning">
            <div class="metric-value">{stats.boundary_count}</div>
            <div class="metric-label">Genome Pairs</div>
            <div class="metric-subtext">Boundary zone (AAI 58-65%)</div>
        </div>
        <div class="metric-card danger">
            <div class="metric-value">{stats.different_category_count}</div>
            <div class="metric-label">Genome Pairs</div>
            <div class="metric-subtext">Different genera (AAI < 58%)</div>
        </div>
    </div>
    """


def build_phylogenetic_context_heatmap(
    similarity_matrix: pl.DataFrame,
    novel_clusters: list,
    genome_labels_map: dict[str, str],
    similarity_type: str = "ANI",
    default_value: float | None = None,
    max_references: int = 50,
    max_clusters: int = 20,
) -> tuple[go.Figure, dict, bool]:
    """
    Build heatmap showing novel clusters in phylogenetic context.

    Creates an extended similarity heatmap (ANI or AAI) that includes both
    reference genomes and novel clusters, positioning novel taxa near their
    closest relatives.

    Args:
        similarity_matrix: Similarity matrix DataFrame (wide format) - ANI or AAI
        novel_clusters: List of NovelCluster objects
        genome_labels_map: Mapping from accession to display label
        similarity_type: "ANI" for nucleotide or "AAI" for protein comparisons
        default_value: Default value for missing/distant pairs
            (defaults to 70 for ANI, 40 for AAI)
        max_references: Maximum reference genomes to include
        max_clusters: Maximum novel clusters to include

    Returns:
        Tuple of (plotly_figure, metadata_dict, clustering_succeeded)
        - plotly_figure: Plotly heatmap figure
        - metadata_dict: Dictionary with counts and cluster info
        - clustering_succeeded: Whether hierarchical clustering worked
    """
    from metadarkmatter.visualization.report.components.extended_matrix_builder import (
        build_extended_similarity_matrix,
    )

    # Set defaults based on similarity type
    is_aai = similarity_type.upper() == "AAI"
    if default_value is None:
        default_value = 40.0 if is_aai else 70.0

    # Build extended matrix
    extended_matrix, labels, is_novel_mask = build_extended_similarity_matrix(
        similarity_matrix=similarity_matrix,
        novel_clusters=novel_clusters,
        genome_labels_map=genome_labels_map,
        similarity_type=similarity_type,
        default_value=default_value,
        max_references=max_references,
        max_clusters=max_clusters,
    )

    if len(labels) == 0:
        # Empty result - return empty figure
        fig = go.Figure()
        fig.update_layout(
            title=f"No data available for {similarity_type} phylogenetic context heatmap",
            template="plotly_white",
        )
        return fig, {"n_references": 0, "n_clusters": 0, "similarity_type": similarity_type}, False

    # Count entities
    n_references = sum(1 for is_novel in is_novel_mask if not is_novel)
    n_clusters = sum(1 for is_novel in is_novel_mask if is_novel)

    # Attempt hierarchical clustering
    clustering_succeeded = False
    try:
        z_clustered, labels_ordered, is_novel_ordered, clustering_succeeded = (
            _cluster_extended_matrix(extended_matrix, labels, is_novel_mask, default_value)
        )
    except Exception:
        # Fallback to unclustered
        z_clustered = extended_matrix
        labels_ordered = labels
        is_novel_ordered = is_novel_mask

    # Convert to Python list for JSON serialization
    z_list = [[float(v) for v in row] for row in z_clustered]

    # Colorscale and range based on similarity type
    if is_aai:
        # AAI colorscale: 40-100% range (full range for proper visualization)
        # Key threshold: ~65% AAI = genus boundary
        colorscale = [
            [0.0, "#313695"],     # 40% - Dark blue (very distant)
            [0.167, "#4575b4"],   # 50%
            [0.333, "#abd9e9"],   # 60%
            [0.417, "#ffffbf"],   # 65% - Yellow (genus boundary)
            [0.5, "#fee090"],     # 70%
            [0.583, "#fdae61"],   # 75%
            [0.75, "#f46d43"],    # 85%
            [1.0, "#a50026"],     # 100% - Dark red (identical)
        ]
        zmin, zmax = 40, 100
        colorbar_tickvals = [40, 50, 60, 65, 70, 80, 90, 100]
        colorbar_ticktext = ["40", "50", "60", "65 (genus)", "70", "80", "90", "100"]
        boundary_label = "genus"
    else:
        # ANI colorscale: 70-100% range (species-level comparisons)
        colorscale = [
            [0.0, "#4575b4"],    # 70% - Blue (distant/different genera)
            [0.167, "#74add1"],  # 75%
            [0.333, "#abd9e9"],  # 80%
            [0.5, "#ffffbf"],    # 85% - Yellow (genus boundary)
            [0.667, "#fee090"],  # 90%
            [0.833, "#fdae61"],  # 95% - Orange (species boundary)
            [1.0, "#d73027"],    # 100% - Red (identical/same species)
        ]
        zmin, zmax = 70, 100
        colorbar_tickvals = [70, 75, 80, 85, 90, 95, 100]
        colorbar_ticktext = ["70", "75", "80", "85 (genus)", "90", "95 (species)", "100"]
        boundary_label = "species"

    # Create custom hover text to indicate novel vs reference
    hover_text = []
    for i, label_i in enumerate(labels_ordered):
        row_hover = []
        is_novel_i = is_novel_ordered[i]
        type_i = "Novel cluster" if is_novel_i else "Reference genome"
        for j, label_j in enumerate(labels_ordered):
            is_novel_j = is_novel_ordered[j]
            type_j = "Novel cluster" if is_novel_j else "Reference genome"
            sim_val = z_clustered[i, j]

            # Determine if this is an estimated value
            is_estimated = is_novel_i or is_novel_j
            estimate_note = " (estimated)" if is_estimated and i != j else ""

            hover_text_val = (
                f"{type_i}: {label_i}<br>"
                f"{type_j}: {label_j}<br>"
                f"{similarity_type}: {sim_val:.1f}%{estimate_note}"
            )
            row_hover.append(hover_text_val)
        hover_text.append(row_hover)

    # Create heatmap
    fig = go.Figure(
        data=go.Heatmap(
            z=z_list,
            x=labels_ordered,
            y=labels_ordered,
            colorscale=colorscale,
            zmin=zmin,
            zmax=zmax,
            colorbar=dict(
                title=f"{similarity_type} (%)",
                tickvals=colorbar_tickvals,
                ticktext=colorbar_ticktext,
            ),
            hovertext=hover_text,
            hovertemplate="%{hovertext}<extra></extra>",
        )
    )

    # Build title
    title = f"Phylogenetic Context Heatmap ({similarity_type})"
    if clustering_succeeded:
        title += " - Hierarchically Clustered"

    # Calculate appropriate dimensions based on number of entities
    n_total = len(labels_ordered)
    height = max(600, min(1200, 50 + n_total * 15))

    # Dynamic width based on matrix size for better label readability
    if n_total > 50:
        width = 1400
    elif n_total > 30:
        width = 1200
    else:
        width = 1000

    fig.update_layout(
        title=title,
        width=width,
        height=height,
        template="plotly_white",
        xaxis=dict(tickangle=45, tickfont=dict(size=8)),
        yaxis=dict(tickfont=dict(size=8)),
    )

    # Metadata for reporting
    metadata = {
        "n_references": n_references,
        "n_clusters": n_clusters,
        "n_total": n_total,
        "clustering_succeeded": clustering_succeeded,
        "similarity_type": similarity_type,
    }

    return fig, metadata, clustering_succeeded


def _cluster_extended_matrix(
    matrix: np.ndarray,
    labels: list[str],
    is_novel_mask: list[bool],
    default_value: float,
) -> tuple[np.ndarray, list[str], list[bool], bool]:
    """
    Perform hierarchical clustering on extended matrix.

    Args:
        matrix: Extended ANI matrix
        labels: Labels for matrix entries
        is_novel_mask: Boolean mask indicating novel clusters
        default_value: Default ANI value

    Returns:
        Tuple of (clustered_matrix, ordered_labels, ordered_mask, success)
    """
    try:
        from scipy.cluster.hierarchy import leaves_list, linkage
        from scipy.spatial.distance import squareform

        # Fill diagonal
        np.fill_diagonal(matrix, 100.0)

        # Convert similarity to distance
        dist_matrix = 100.0 - matrix
        dist_matrix = (dist_matrix + dist_matrix.T) / 2
        np.fill_diagonal(dist_matrix, 0.0)

        # Hierarchical clustering
        condensed_dist = squareform(dist_matrix)
        linkage_matrix = linkage(condensed_dist, method="average")
        order = leaves_list(linkage_matrix)

        # Reorder everything
        clustered = matrix[order][:, order]
        ordered_labels = [labels[i] for i in order]
        ordered_mask = [is_novel_mask[i] for i in order]

        return clustered, ordered_labels, ordered_mask, True

    except ImportError:
        return matrix, labels, is_novel_mask, False
