"""
Visualization module for metadarkmatter.

Provides Plotly-based visualization tools for recruitment plots
and taxonomic classification results.
"""

from metadarkmatter.visualization.recruitment_plots import (
    IdentityBand,
    RecruitmentPlotGenerator,
    create_recruitment_plot,
)

__all__ = [
    "IdentityBand",
    "RecruitmentPlotGenerator",
    "create_recruitment_plot",
]
