"""
Sparse ANI matrix implementation for very large genome datasets.

This module provides SparseANIMatrix, an alternative to ANIMatrix optimized
for datasets with 10K+ genomes where most genome pairs have no precomputed
ANI values. Stores only non-zero values to minimize memory usage.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from metadarkmatter.core.parsers import ANIMatrixParser


class SparseANIMatrix:
    """
    Sparse ANI matrix for very large genome sets (10K+ genomes).

    For genome databases with 10K+ genomes, a dense matrix becomes
    impractical (10K x 10K = 100M cells = 400MB). This sparse
    representation stores only non-zero ANI values.

    Memory comparison (10K genomes):
    - Dense: 400 MB
    - Sparse (5% density): ~20 MB

    Falls back to default ANI value (70.0) for missing pairs,
    which is typical for distantly related genomes.
    """

    __slots__ = ("_ani_sparse", "_default_ani", "_genome_to_idx", "_genomes")

    def __init__(
        self,
        ani_dict: dict[str, dict[str, float]],
        default_ani: float = 70.0,
        min_ani: float = 75.0,
    ) -> None:
        """
        Initialize sparse ANI matrix.

        Args:
            ani_dict: Nested dict {genome1: {genome2: ani_value}}
            default_ani: Default ANI for missing pairs (default: 70.0)
            min_ani: Minimum ANI to store; below this, use default (default: 75.0)
        """
        self._genomes: tuple[str, ...] = tuple(sorted(ani_dict.keys()))
        self._genome_to_idx: dict[str, int] = {
            g: i for i, g in enumerate(self._genomes)
        }
        self._default_ani = default_ani

        # Store only significant ANI values as (i, j) -> ani
        self._ani_sparse: dict[tuple[int, int], float] = {}

        for genome1, inner_dict in ani_dict.items():
            i = self._genome_to_idx[genome1]
            for genome2, ani_value in inner_dict.items():
                if genome2 in self._genome_to_idx and ani_value >= min_ani:
                    j = self._genome_to_idx[genome2]
                    # Store only upper triangle (symmetric)
                    key = (min(i, j), max(i, j))
                    self._ani_sparse[key] = ani_value

    @classmethod
    def from_file(cls, path: Path, **kwargs: Any) -> SparseANIMatrix:
        """Load sparse ANI matrix from file."""
        parser = ANIMatrixParser(path)
        ani_dict = parser.to_dict()
        return cls(ani_dict, **kwargs)

    def __len__(self) -> int:
        """Number of genomes."""
        return len(self._genomes)

    @property
    def genomes(self) -> set[str]:
        """Set of genome names."""
        return set(self._genomes)

    def get_ani(self, genome1: str, genome2: str) -> float:
        """Get ANI between two genomes."""
        if genome1 == genome2:
            return 100.0

        i = self._genome_to_idx.get(genome1)
        j = self._genome_to_idx.get(genome2)

        if i is None or j is None:
            return self._default_ani

        key = (min(i, j), max(i, j))
        return self._ani_sparse.get(key, self._default_ani)

    def has_genome(self, genome: str) -> bool:
        """Check if genome exists."""
        return genome in self._genome_to_idx

    def memory_usage_bytes(self) -> int:
        """Estimate memory usage."""
        # Each entry: 2 ints (16 bytes) + 1 float (8 bytes) + dict overhead (~50 bytes)
        sparse_bytes = len(self._ani_sparse) * 74
        # Genome index dict
        index_bytes = len(self._genomes) * 100
        return sparse_bytes + index_bytes

    def density(self) -> float:
        """Calculate matrix density (fraction of non-default values)."""
        n = len(self._genomes)
        total_possible = n * (n - 1) // 2  # Upper triangle
        if total_possible == 0:
            return 0.0
        return len(self._ani_sparse) / total_possible
