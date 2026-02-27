"""
Novel diversity analysis module.

Provides tools for clustering novel reads into putative taxa
and analyzing the phylogenetic context of novel diversity.
"""

from metadarkmatter.core.novel_diversity.clustering import NovelDiversityAnalyzer
from metadarkmatter.core.novel_diversity.models import (
    GenusDistance,
    NovelCluster,
    NovelDiversitySummary,
    PhylogeneticNeighborhood,
)
from metadarkmatter.core.novel_diversity.neighborhood import (
    PhylogeneticNeighborhoodAnalyzer,
)

__all__ = [
    "GenusDistance",
    "NovelCluster",
    "NovelDiversityAnalyzer",
    "NovelDiversitySummary",
    "PhylogeneticNeighborhood",
    "PhylogeneticNeighborhoodAnalyzer",
]
