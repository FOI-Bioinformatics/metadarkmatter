"""
Novel diversity analysis module.

Provides tools for clustering novel reads into putative taxa
and analyzing the phylogenetic context of novel diversity.
"""

from metadarkmatter.core.novel_diversity.clustering import NovelDiversityAnalyzer
from metadarkmatter.core.novel_diversity.models import (
    NovelCluster,
    NovelDiversitySummary,
)

__all__ = [
    "NovelCluster",
    "NovelDiversityAnalyzer",
    "NovelDiversitySummary",
]
