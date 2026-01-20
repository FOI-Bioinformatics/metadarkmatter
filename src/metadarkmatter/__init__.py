"""
Metadarkmatter: ANI-weighted placement for detecting novel microbial diversity.

A tool for analyzing environmental DNA (eDNA) metagenomic data from air filters,
water samples, and other environmental sources to characterize microbial diversity
and detect novel bacterial taxa using whole-genome competitive read recruitment.
"""

__version__ = "0.1.0"
__author__ = "Metadarkmatter Team"

from metadarkmatter.core.ani_placement import ANIWeightedClassifier
from metadarkmatter.models.classification import (
    ReadClassification,
    TaxonomicCall,
    TaxonomicSummary,
)

__all__ = [
    "ANIWeightedClassifier",
    "ReadClassification",
    "TaxonomicCall",
    "TaxonomicSummary",
    "__version__",
]
