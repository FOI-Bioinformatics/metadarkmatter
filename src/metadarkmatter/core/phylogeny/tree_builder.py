"""Build phylogenetic trees from ANI matrices or genome assemblies.

This module provides functions to convert ANI (Average Nucleotide Identity)
matrices into neighbor-joining or UPGMA phylogenetic trees in Newick format,
to build trees from genome FASTA files using Mashtree, and to load and
validate user-provided phylogenetic trees.
"""

from __future__ import annotations

import logging
from enum import StrEnum
from io import StringIO
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from Bio.Phylo.BaseTree import Tree

logger = logging.getLogger(__name__)


class TreeMethod(StrEnum):
    """Phylogenetic tree building method."""

    NJ = "nj"
    UPGMA = "upgma"
    MASHTREE = "mashtree"


MISSING_ANI_FALLBACK = 70.0
"""ANI (percent) assumed for pairs that were not computed (0 or NaN).

Matches the classifier's default for missing/uncomputed pairs
(``ANIMatrix.default_ani``), so the tree and the classification treat an
undetected pair the same way (distance 30).
"""


def count_missing_pairs(ani_matrix: pd.DataFrame) -> int:
    """Count off-diagonal genome pairs whose ANI is missing (NaN or 0)."""
    ani_matrix = _align_rows_to_columns(ani_matrix)
    values = ani_matrix.to_numpy(dtype=float)
    n = min(values.shape)
    sub = values[:n, :n]
    sub = np.where(sub > 0, sub, sub.T)
    missing = ~(sub > 0)
    for i in range(n):
        missing[i, i] = False
    return int(missing.sum() // 2)


def _ani_to_distance(ani_matrix: pd.DataFrame) -> pd.DataFrame:
    """Convert ANI to distance, treating 0/NaN as 'not computed' (70% ANI)."""
    ani = ani_matrix.astype(float).copy()
    if ani.shape[0] == ani.shape[1]:
        # Mirror one-directional pairs (a 0/NaN cell with a value in the
        # opposite direction is not missing).
        mirrored = ani.T.set_axis(ani.index, axis=0).set_axis(ani.columns, axis=1)
        ani = ani.where(ani > 0, mirrored)
    missing = ~(ani > 0)
    n_missing = count_missing_pairs(ani)
    if n_missing > 0:
        logger.warning(
            f"{n_missing} genome pairs lack ANI values (0 or NaN); "
            f"assuming {MISSING_ANI_FALLBACK:.0f}% ANI"
        )
    ani = ani.mask(missing, MISSING_ANI_FALLBACK)
    # Diagonal is always identical to itself.
    for i in range(min(ani.shape)):
        ani.iloc[i, i] = 100.0
    return (100 - ani).clip(lower=0, upper=50)


def _align_rows_to_columns(ani_matrix: pd.DataFrame) -> pd.DataFrame:
    """Reorder rows to match column order (label-based, not positional).

    When the row index carries the same genome labels as the columns, the
    matrix is re-indexed as ``ani.loc[cols, cols]`` so a different row order
    cannot silently pair the wrong distances with the wrong names. Matrices
    without a matching row index are returned unchanged.
    """
    cols = list(ani_matrix.columns)
    if set(ani_matrix.index) == set(cols) and len(ani_matrix.index) == len(cols):
        return ani_matrix.loc[cols, cols]
    return ani_matrix


def _to_lower_triangular(matrix: pd.DataFrame) -> list[list[float]]:
    """Convert symmetric matrix to lower triangular list format for BioPython.

    BioPython's DistanceMatrix requires a lower triangular matrix format
    where row i contains values for columns 0 to i (inclusive).

    Args:
        matrix: Symmetric distance matrix as DataFrame.

    Returns:
        Lower triangular matrix as nested lists.
    """
    n = len(matrix)
    result: list[list[float]] = []
    for i in range(n):
        row: list[float] = []
        for j in range(i + 1):
            row.append(float(matrix.iloc[i, j]))
        result.append(row)
    return result


def ani_to_newick(ani_matrix: pd.DataFrame) -> str | None:
    """Convert ANI matrix to neighbor-joining tree in Newick format.

    Uses BioPython's neighbor-joining implementation to construct a
    phylogenetic tree from ANI values. The tree is rooted at its midpoint.

    Args:
        ani_matrix: Square DataFrame with ANI values (0-100 scale) where
            rows and columns are genome identifiers. Diagonal should be 100.

    Returns:
        Newick format string if tree construction succeeds, None if the
        matrix has fewer than 3 genomes (minimum required for NJ algorithm).

    Note:
        Missing ANI values (NaN or 0, i.e. not computed) are treated as 70%
        ANI (distance 30), the same fallback the classifier uses.
    """
    from Bio import Phylo
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

    if len(ani_matrix) < 3:
        if len(ani_matrix) > 0:
            logger.warning(
                f"Too few genomes for phylogenetic tree (need >= 3). Got {len(ani_matrix)}."
            )
        return None

    ani_matrix = _align_rows_to_columns(ani_matrix)

    # Convert ANI (similarity) to distance (100 - ANI); 0/NaN = not computed.
    distance_matrix = _ani_to_distance(ani_matrix)

    # Build BioPython DistanceMatrix
    names = list(ani_matrix.columns)
    matrix = _to_lower_triangular(distance_matrix)
    dm = DistanceMatrix(names, matrix)

    # Construct neighbor-joining tree
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    # Root at midpoint if there is variation in distances
    # BioPython's root_at_midpoint fails when all distances are zero
    try:
        tree.root_at_midpoint()
    except (UnboundLocalError, ValueError):
        # All distances are identical (e.g., all 100% ANI)
        # Tree remains unrooted or arbitrarily rooted
        logger.debug("Could not root at midpoint (all distances may be identical)")

    # Convert to Newick string
    output = StringIO()
    Phylo.write(tree, output, "newick")
    return output.getvalue().strip()


def ani_to_upgma(ani_matrix: pd.DataFrame) -> str | None:
    """Convert ANI matrix to UPGMA tree in Newick format.

    Uses BioPython's UPGMA implementation. Same distance conversion
    as ani_to_newick() but produces an ultrametric tree.

    Args:
        ani_matrix: Square DataFrame with ANI values (0-100 scale) where
            rows and columns are genome identifiers. Diagonal should be 100.

    Returns:
        Newick format string if tree construction succeeds, None if the
        matrix has fewer than 3 genomes.

    Note:
        Missing ANI values (NaN or 0) are treated as 70% ANI, as in the classifier.
    """
    from Bio import Phylo
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

    if len(ani_matrix) < 3:
        if len(ani_matrix) > 0:
            logger.warning(
                f"Too few genomes for phylogenetic tree (need >= 3). Got {len(ani_matrix)}."
            )
        return None

    ani_matrix = _align_rows_to_columns(ani_matrix)
    distance_matrix = _ani_to_distance(ani_matrix)

    names = list(ani_matrix.columns)
    matrix = _to_lower_triangular(distance_matrix)
    dm = DistanceMatrix(names, matrix)

    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)

    output = StringIO()
    Phylo.write(tree, output, "newick")
    return output.getvalue().strip()


def mashtree_to_newick(
    genome_dir: Path,
    genome_pattern: str = "*.fna",
    threads: int = 4,
) -> str:
    """Build tree from genome FASTA files using Mashtree.

    Mashtree computes Mash (MinHash) distances between genome assemblies
    and constructs a neighbor-joining tree.

    Args:
        genome_dir: Directory containing genome FASTA files.
        genome_pattern: Glob pattern for genome files (default: "*.fna").
        threads: Number of threads for mashtree.

    Returns:
        Newick format string.

    Raises:
        FileNotFoundError: If no genome files found matching the pattern.
        ValueError: If fewer than 3 genome files found.
    """
    from metadarkmatter.external.mashtree import Mashtree

    # Find genome files, trying multiple patterns if primary fails
    genome_files = sorted(genome_dir.glob(genome_pattern))
    if not genome_files:
        for alt in ("*.fna", "*.fa", "*.fasta", "*.fna.gz", "*.fa.gz"):
            if alt != genome_pattern:
                genome_files = sorted(genome_dir.glob(alt))
                if genome_files:
                    break

    if not genome_files:
        raise FileNotFoundError(
            f"No genome files found in {genome_dir} matching '{genome_pattern}'"
        )

    if len(genome_files) < 3:
        raise ValueError(
            f"Too few genome files for tree building (need >= 3, got {len(genome_files)})"
        )

    mt = Mashtree()
    result = mt.run_or_raise(genomes=genome_files, threads=threads)
    newick = result.stdout.strip()

    if not newick:
        raise RuntimeError("Mashtree produced no output")

    logger.info(f"Mashtree tree built from {len(genome_files)} genomes")
    return newick


def build_tree(
    method: TreeMethod,
    *,
    ani_matrix: pd.DataFrame | None = None,
    genome_dir: Path | None = None,
    genome_pattern: str = "*.fna",
    threads: int = 4,
) -> str | None:
    """Build a phylogenetic tree using the specified method.

    Dispatches to the appropriate tree builder based on the method.

    Args:
        method: Tree building method (NJ, UPGMA, or MASHTREE).
        ani_matrix: ANI matrix DataFrame. Required for NJ and UPGMA.
        genome_dir: Directory of genome FASTAs. Required for MASHTREE.
        genome_pattern: Glob pattern for genome files (MASHTREE only).
        threads: Number of threads (MASHTREE only).

    Returns:
        Newick format string, or None if NJ/UPGMA with < 3 genomes.

    Raises:
        ValueError: If required inputs are missing for the chosen method.
    """
    if method in (TreeMethod.NJ, TreeMethod.UPGMA):
        if ani_matrix is None:
            raise ValueError(f"ANI matrix required for {method.value} method")
        if method == TreeMethod.NJ:
            return ani_to_newick(ani_matrix)
        return ani_to_upgma(ani_matrix)

    if method == TreeMethod.MASHTREE:
        if genome_dir is None:
            raise ValueError("genome_dir required for mashtree method")
        return mashtree_to_newick(genome_dir, genome_pattern, threads)

    raise ValueError(f"Unknown method: {method}")


def _prune_tips(tree: Tree, keep: set[str]) -> Tree:
    """Prune tree to only keep specified tip names.

    Removes terminal nodes whose names are not in the keep set.

    Args:
        tree: BioPython Tree object to prune.
        keep: Set of tip names to retain.

    Returns:
        The pruned tree (modified in place).
    """
    tips_to_remove = [tip for tip in tree.get_terminals() if tip.name not in keep]
    for tip in tips_to_remove:
        tree.prune(tip)
    return tree


def load_user_tree(newick_path: Path, ani_matrix: pd.DataFrame) -> str:
    """Load and validate user-provided Newick tree.

    Reads a Newick tree file and validates it against the ANI matrix genomes.
    Extra tips in the tree (not in ANI matrix) are pruned. Genomes in the
    ANI matrix but missing from the tree are logged as warnings.

    Args:
        newick_path: Path to Newick format tree file.
        ani_matrix: ANI matrix DataFrame for validation. Column names are
            used to determine which genomes should be in the tree.

    Returns:
        Newick format string of the validated (and possibly pruned) tree.

    Raises:
        FileNotFoundError: If the newick_path does not exist.
    """
    from Bio import Phylo

    if not newick_path.exists():
        raise FileNotFoundError(f"Tree file not found: {newick_path}")

    tree = Phylo.read(newick_path, "newick")
    tree_tips = {tip.name for tip in tree.get_terminals()}
    matrix_genomes = set(ani_matrix.columns)

    missing = matrix_genomes - tree_tips
    extra = tree_tips - matrix_genomes

    if missing:
        logger.warning(
            f"{len(missing)} genomes in ANI matrix not in tree - excluded from phylogeny"
        )

    if extra:
        tree = _prune_tips(tree, keep=matrix_genomes)
        logger.info(f"Pruned {len(extra)} unused tips from provided tree")

    output = StringIO()
    Phylo.write(tree, output, "newick")
    return output.getvalue().strip()
