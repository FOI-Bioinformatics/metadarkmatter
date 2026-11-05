"""
Report generation components for modular heatmap and statistics building.

This package provides reusable components for building interactive
visualizations and statistics cards in HTML reports.
"""

from metadarkmatter.visualization.report.components.clustering import (
    perform_hierarchical_clustering,
)
from metadarkmatter.visualization.report.components.extended_matrix_builder import (
    build_extended_similarity_matrix,
    estimate_novel_to_novel_similarity,
    estimate_novel_to_reference_similarity,
    select_relevant_references,
)
from metadarkmatter.visualization.report.components.heatmap_builder import (
    SimilarityStats,
    build_aai_heatmap,
    build_aai_stats_cards,
    build_ani_heatmap,
    build_ani_stats_cards,
    build_missing_pairs_note,
    build_phylogenetic_context_heatmap,
    compute_similarity_stats,
    count_missing_pairs,
)

__all__ = [
    "SimilarityStats",
    "build_aai_heatmap",
    "build_aai_stats_cards",
    "build_ani_heatmap",
    "build_ani_stats_cards",
    "build_extended_similarity_matrix",
    "build_missing_pairs_note",
    "build_phylogenetic_context_heatmap",
    "compute_similarity_stats",
    "count_missing_pairs",
    "estimate_novel_to_novel_similarity",
    "estimate_novel_to_reference_similarity",
    "perform_hierarchical_clustering",
    "select_relevant_references",
]
