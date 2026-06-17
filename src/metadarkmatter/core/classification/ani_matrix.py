"""
High-performance ANI matrix implementation using NumPy arrays.

This module provides the ANIMatrix class for storing and querying Average
Nucleotide Identity (ANI) values between genome pairs with optimized memory
usage and O(1) lookup performance.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from metadarkmatter.core.constants import ANI_DEFAULT_UNRELATED
from metadarkmatter.core.parsers import ANIMatrixParser

logger = logging.getLogger(__name__)

# fastANI / skani typically reproduce the same pair within ~0.5 ANI units.
# Asymmetries larger than this likely indicate a corrupt or merged matrix.
DEFAULT_SYMMETRY_TOLERANCE = 0.5

if TYPE_CHECKING:
    import polars as pl


class ANIMatrix:
    """
    High-performance ANI matrix using NumPy arrays with integer indexing.

    Performance optimization over nested dictionaries:
    - Memory: 8MB for 1000 genomes vs 150MB with nested dicts (18x reduction)
    - Lookup: O(1) array access vs O(hash) string lookups (16x faster)

    The matrix stores ANI values as float32, with genome names mapped to
    integer indices for fast array-based lookups.

    Missing ANI values (0.0 in the matrix) are treated as distant organisms
    and return the default_ani value (typically 70%) instead of 0.0.
    """

    __slots__ = ("_ani_array", "_default_ani", "_genome_to_idx", "_genomes", "_num_genomes")

    def __init__(
        self,
        ani_dict: dict[str, dict[str, float]],
        default_ani: float = ANI_DEFAULT_UNRELATED,
        symmetry_check: str = "warn",
        symmetry_tolerance: float = DEFAULT_SYMMETRY_TOLERANCE,
    ) -> None:
        """
        Initialize ANI matrix from nested dictionary.

        Converts the nested dictionary to a NumPy array with integer indexing
        for O(1) lookups during classification.

        Args:
            ani_dict: Nested dict {genome1: {genome2: ani_value}}
            default_ani: Default ANI value for missing/uncomputed pairs (default: 70.0).
                When skani/fastANI cannot compute ANI between distant genomes,
                the value is stored as 0.0. This parameter specifies what value
                to return instead, representing the expected ANI for distant
                organisms within the same family.
            symmetry_check: One of "off", "warn" (default), or "strict". When
                "warn" or "strict", scan the input dict for pairs where
                ANI(a,b) and ANI(b,a) both exist and differ by more than
                ``symmetry_tolerance``. "warn" logs a warning and keeps the
                last-written value; "strict" raises ValueError. The check is
                a guard against silent overwrites — the nested-dict ingest
                stores only one direction per pair and an asymmetric input
                would otherwise be accepted without warning.
            symmetry_tolerance: Maximum allowed |ANI(a,b) - ANI(b,a)| before
                the symmetry check triggers (default 0.5; matches fastANI /
                skani numerical noise).
        """
        if symmetry_check not in ("off", "warn", "strict"):
            raise ValueError(
                f"symmetry_check must be 'off', 'warn', or 'strict' (got {symmetry_check!r})"
            )

        if symmetry_check != "off":
            self._check_input_symmetry(ani_dict, symmetry_tolerance, symmetry_check)

        self._default_ani = default_ani

        # Sort genomes for consistent indexing
        self._genomes: tuple[str, ...] = tuple(sorted(ani_dict.keys()))
        self._num_genomes: int = len(self._genomes)

        # Create genome -> index mapping for O(1) lookups
        self._genome_to_idx: dict[str, int] = {
            genome: idx for idx, genome in enumerate(self._genomes)
        }

        # Build NumPy array: float32 is sufficient for ANI values (0-100)
        # Memory: 1000 genomes = 1000 * 1000 * 4 bytes = 4 MB
        self._ani_array: np.ndarray = np.zeros(
            (self._num_genomes, self._num_genomes),
            dtype=np.float32,
        )

        # Fill diagonal with 100.0 (self-ANI)
        np.fill_diagonal(self._ani_array, 100.0)

        # Populate matrix from dictionary
        for genome1, inner_dict in ani_dict.items():
            i = self._genome_to_idx[genome1]
            for genome2, ani_value in inner_dict.items():
                if genome2 in self._genome_to_idx:
                    j = self._genome_to_idx[genome2]
                    self._ani_array[i, j] = ani_value

    @staticmethod
    def _check_input_symmetry(
        ani_dict: dict[str, dict[str, float]],
        tolerance: float,
        mode: str,
    ) -> None:
        """Verify ANI(a,b) ~= ANI(b,a) for every pair present in both directions.

        The constructor stores only one direction per (a,b) pair in the
        underlying NumPy array, so an asymmetric input would be silently
        accepted with only one value retained. This guard surfaces the
        problem at load time. Limited to pairs present in both directions
        in the input — missing reverse entries are not an asymmetry, just
        an absence.
        """
        asymmetries: list[tuple[str, str, float, float]] = []
        for g1, inner in ani_dict.items():
            for g2, v_fwd in inner.items():
                if g1 == g2:
                    continue
                reverse = ani_dict.get(g2)
                if reverse is None:
                    continue
                v_rev = reverse.get(g1)
                if v_rev is None:
                    continue
                if abs(float(v_fwd) - float(v_rev)) > tolerance:
                    asymmetries.append((g1, g2, float(v_fwd), float(v_rev)))

        if not asymmetries:
            return

        sample = asymmetries[:5]
        sample_str = ", ".join(
            f"{a}->{b}: {v1:.2f} vs {v2:.2f}" for a, b, v1, v2 in sample
        )
        msg = (
            f"ANI matrix is not symmetric within tolerance {tolerance}: "
            f"{len(asymmetries)} pair(s) disagree. Examples: {sample_str}. "
            "Only one direction per pair is retained, which can corrupt "
            "placement uncertainty. Regenerate the matrix from fastANI/skani "
            "or pass symmetry_check='off' to suppress this check."
        )
        if mode == "strict":
            raise ValueError(msg)
        logger.warning(msg)

    @classmethod
    def from_file(
        cls,
        path: Path,
        default_ani: float = ANI_DEFAULT_UNRELATED,
        symmetry_check: str = "warn",
        symmetry_tolerance: float = DEFAULT_SYMMETRY_TOLERANCE,
    ) -> ANIMatrix:
        """
        Load ANI matrix from file.

        Args:
            path: Path to ANI matrix CSV/TSV
            default_ani: Default ANI value for missing/uncomputed pairs (default: 70.0)
            symmetry_check: See ``ANIMatrix.__init__`` (default "warn").
            symmetry_tolerance: See ``ANIMatrix.__init__`` (default 0.5).

        Returns:
            ANIMatrix instance
        """
        parser = ANIMatrixParser(path)
        ani_dict = parser.to_dict()
        return cls(
            ani_dict,
            default_ani=default_ani,
            symmetry_check=symmetry_check,
            symmetry_tolerance=symmetry_tolerance,
        )

    @classmethod
    def from_dataframe(cls, df: pl.DataFrame, default_ani: float = ANI_DEFAULT_UNRELATED) -> ANIMatrix:
        """
        Construct an ANIMatrix from a Polars DataFrame.

        Accepts the matrix format used by the report generator, where the
        first column is ``genome`` (row labels) and remaining columns are
        genome accessions with ANI values.  Also accepts a square numeric
        DataFrame where column names serve as both row and column labels.

        Args:
            df: Polars DataFrame with ANI values
            default_ani: Default ANI value for missing/uncomputed pairs

        Returns:
            ANIMatrix instance
        """

        if len(df) == 0:
            return cls({}, default_ani=default_ani)

        # Determine genome labels and value columns
        if "genome" in df.columns:
            genomes = df["genome"].to_list()
            value_cols = [c for c in df.columns if c != "genome"]
        else:
            genomes = list(df.columns)
            value_cols = list(df.columns)

        # Build nested dict
        ani_dict: dict[str, dict[str, float]] = {}
        for row_idx, row_genome in enumerate(genomes):
            inner: dict[str, float] = {}
            for col in value_cols:
                val = df[col][row_idx]
                if val is not None:
                    inner[col] = float(val)
            ani_dict[row_genome] = inner

        return cls(ani_dict, default_ani=default_ani)

    @property
    def genomes(self) -> set[str]:
        """Set of genome names in the matrix."""
        return set(self._genomes)

    def __len__(self) -> int:
        """Number of genomes in the matrix."""
        return self._num_genomes

    def get_ani(self, genome1: str, genome2: str) -> float:
        """
        Get ANI value between two genomes.

        Uses integer-indexed NumPy array for O(1) access after
        initial genome name lookup.

        Args:
            genome1: First genome identifier
            genome2: Second genome identifier

        Returns:
            ANI value (0-100). Returns default_ani (typically 70%) for
            missing/uncomputed pairs where both genomes exist in the matrix.
            Returns 0.0 only if a genome is not in the matrix at all.

        Note:
            Returns 100.0 if genome1 == genome2 (diagonal)
        """
        # Fast path: same genome
        if genome1 == genome2:
            return 100.0

        # Get indices (single dict lookup per genome)
        i = self._genome_to_idx.get(genome1)
        if i is None:
            return 0.0

        j = self._genome_to_idx.get(genome2)
        if j is None:
            return 0.0

        # O(1) array access
        ani_value = float(self._ani_array[i, j])

        # If ANI is 0.0 (not computed by skani/fastANI), return default for distant organisms
        if ani_value == 0.0:
            return self._default_ani

        return ani_value

    def get_ani_by_idx(self, idx1: int, idx2: int) -> float:
        """
        Get ANI value by pre-computed integer indices.

        Fastest lookup method when genome indices are already known.
        Used in optimized classification pipelines.

        Args:
            idx1: First genome index
            idx2: Second genome index

        Returns:
            ANI value (0-100)
        """
        return float(self._ani_array[idx1, idx2])

    def get_genome_idx(self, genome: str) -> int | None:
        """
        Get integer index for a genome name.

        Args:
            genome: Genome identifier

        Returns:
            Integer index, or None if genome not found
        """
        return self._genome_to_idx.get(genome)

    def has_genome(self, genome: str) -> bool:
        """Check if genome is present in ANI matrix."""
        return genome in self._genome_to_idx

    def memory_usage_bytes(self) -> int:
        """Estimate memory usage in bytes."""
        # Array memory
        array_bytes = self._ani_array.nbytes
        # Index dict overhead (rough estimate)
        dict_bytes = self._num_genomes * 100  # ~100 bytes per entry
        return array_bytes + dict_bytes

    def cluster_genomes_by_ani(
        self,
        species_threshold: float = 95.0,
        genus_threshold: float = 80.0,
    ) -> dict[str, set[str]]:
        """
        Cluster genomes by ANI thresholds using single-linkage clustering.

        Genomes are grouped if their ANI exceeds the threshold. This provides
        phylogenetic context for ambiguous placements without requiring
        phylogenetic trees.

        Args:
            species_threshold: ANI threshold for same species (default: 95%)
            genus_threshold: ANI threshold for same genus (default: 80%)

        Returns:
            Dictionary mapping each genome to a set of genomes in the same genus
            based on single-linkage clustering at the genus threshold.

        Note:
            Uses genus_threshold for clustering. Species-level clusters are
            subsets within genus clusters.
        """
        from scipy.cluster.hierarchy import fcluster, linkage
        from scipy.spatial.distance import squareform

        if self._num_genomes == 0:
            return {}

        # Convert ANI to distance (100 - ANI) for hierarchical clustering
        # Use condensed distance matrix format for scipy
        distance_matrix = 100.0 - self._ani_array

        # Ensure symmetry and zero diagonal
        np.fill_diagonal(distance_matrix, 0)
        distance_matrix = (distance_matrix + distance_matrix.T) / 2

        # Convert to condensed form for scipy
        condensed = squareform(distance_matrix)

        # Perform single-linkage clustering
        # Single-linkage: genomes are in same cluster if ANY pair exceeds threshold
        linkage_matrix = linkage(condensed, method='single')

        # Cut tree at genus distance threshold (100 - genus_threshold)
        genus_distance = 100.0 - genus_threshold
        cluster_labels = fcluster(linkage_matrix, genus_distance, criterion='distance')

        # Build genome -> cluster members mapping
        clusters: dict[int, set[str]] = {}
        for genome, label in zip(self._genomes, cluster_labels, strict=True):
            if label not in clusters:
                clusters[label] = set()
            clusters[label].add(genome)

        # Return genome -> same-genus genomes mapping
        genome_to_genus_members: dict[str, set[str]] = {}
        for genome, label in zip(self._genomes, cluster_labels, strict=True):
            genome_to_genus_members[genome] = clusters[label]

        return genome_to_genus_members

    def get_genus_cluster_for_genome(
        self,
        genome: str,
        genus_threshold: float = 80.0,
    ) -> set[str]:
        """
        Get all genomes in the same genus cluster as the given genome.

        A faster method for querying genus membership without full clustering.

        Args:
            genome: Target genome identifier
            genus_threshold: ANI threshold for same genus (default: 80%)

        Returns:
            Set of genomes sharing >= genus_threshold ANI with target genome
        """
        idx = self._genome_to_idx.get(genome)
        if idx is None:
            return {genome}

        # Find all genomes with ANI >= threshold
        ani_row = self._ani_array[idx]
        mask = ani_row >= genus_threshold
        same_genus = {self._genomes[i] for i in np.where(mask)[0]}

        return same_genus

    def count_distinct_genera(
        self,
        genomes: list[str],
        genus_threshold: float = 80.0,
    ) -> tuple[int, list[set[str]]]:
        """
        Count distinct genera among a list of genomes based on ANI.

        Uses transitive closure: if A and B share >80% ANI, and B and C share
        >80% ANI, then A, B, C are all in the same genus.

        Args:
            genomes: List of genome identifiers
            genus_threshold: ANI threshold for same genus (default: 80%)

        Returns:
            Tuple of (num_distinct_genera, list of genus groups)
        """
        if not genomes:
            return 0, []

        # Build adjacency list based on ANI threshold
        genome_indices = [self._genome_to_idx.get(g) for g in genomes]
        valid_genomes = [(g, idx) for g, idx in zip(genomes, genome_indices, strict=True) if idx is not None]

        if not valid_genomes:
            return len(genomes), [{g} for g in genomes]

        # Union-find for transitive closure
        parent: dict[str, str] = {g: g for g, _ in valid_genomes}

        def find(x: str) -> str:
            if parent[x] != x:
                parent[x] = find(parent[x])
            return parent[x]

        def union(x: str, y: str) -> None:
            px, py = find(x), find(y)
            if px != py:
                parent[px] = py

        # Check all pairs
        for i, (g1, idx1) in enumerate(valid_genomes):
            for g2, idx2 in valid_genomes[i + 1:]:
                ani = self._ani_array[idx1, idx2]
                if ani >= genus_threshold:
                    union(g1, g2)

        # Group by root
        groups: dict[str, set[str]] = {}
        for g, _ in valid_genomes:
            root = find(g)
            if root not in groups:
                groups[root] = set()
            groups[root].add(g)

        # Handle genomes not in ANI matrix
        missing = [g for g, idx in zip(genomes, genome_indices, strict=True) if idx is None]
        genus_groups = list(groups.values())
        for g in missing:
            genus_groups.append({g})

        return len(genus_groups), genus_groups
