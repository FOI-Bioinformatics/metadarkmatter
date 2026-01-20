"""
ANI-weighted placement uncertainty algorithm for detecting novel microbial diversity.

This module re-exports classification classes from the refactored classification package
for backward compatibility. New code should import directly from:
- metadarkmatter.core.classification.ani_matrix.ANIMatrix
- metadarkmatter.core.classification.classifiers.ANIWeightedClassifier
- metadarkmatter.core.classification.classifiers.ParallelClassifier
- metadarkmatter.core.classification.classifiers.VectorizedClassifier
- metadarkmatter.core.classification.sparse_ani_matrix.SparseANIMatrix
"""

from __future__ import annotations

# Re-export all classes for backward compatibility
from metadarkmatter.core.classification import (
    ANIMatrix,
    ANIWeightedClassifier,
    ParallelClassifier,
    SparseANIMatrix,
    VectorizedClassifier,
)

__all__ = [
    "ANIMatrix",
    "ANIWeightedClassifier",
    "ParallelClassifier",
    "SparseANIMatrix",
    "VectorizedClassifier",
]
