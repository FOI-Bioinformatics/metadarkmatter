"""
Report generation components for modular heatmap and statistics building.

This package provides reusable components for building interactive
visualizations and statistics cards in HTML reports.
"""

from metadarkmatter.visualization.report.components.clustering import (
    perform_hierarchical_clustering,
)
from metadarkmatter.visualization.report.components.heatmap_builder import (
    SimilarityStats,
    build_aai_heatmap,
    build_aai_stats_cards,
    build_ani_heatmap,
    build_ani_stats_cards,
    compute_similarity_stats,
)

__all__ = [
    "SimilarityStats",
    "build_aai_heatmap",
    "build_aai_stats_cards",
    "build_ani_heatmap",
    "build_ani_stats_cards",
    "compute_similarity_stats",
    "perform_hierarchical_clustering",
]
