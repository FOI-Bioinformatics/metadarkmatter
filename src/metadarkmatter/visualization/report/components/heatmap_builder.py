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

    # AAI colorscale: adjusted range for protein-level comparisons
    aai_colorscale = [
        [0.0, "#313695"],    # 40% - Dark blue (very distant)
        [0.15, "#4575b4"],   # 45%
        [0.3, "#74add1"],    # 50%
        [0.45, "#abd9e9"],   # 55%
        [0.6, "#ffffbf"],    # 60% - Yellow (approaching genus)
        [0.75, "#fdae61"],   # 65% - Orange (genus boundary)
        [0.9, "#f46d43"],    # 70%
        [1.0, "#a50026"],    # 75%+ - Red (same genus)
    ]

    # Create heatmap
    fig = go.Figure(
        data=go.Heatmap(
            z=z_list,
            x=genome_labels,
            y=genome_labels,
            colorscale=aai_colorscale,
            zmin=40,
            zmax=75,
            colorbar=dict(
                title="AAI (%)",
                tickvals=[40, 45, 50, 55, 60, 65, 70, 75],
                ticktext=["40", "45", "50", "55", "60", "65 (genus)", "70", "75+"],
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
