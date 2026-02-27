# Tree Build Command Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a `metadarkmatter tree build` CLI command that builds phylogenetic trees using NJ, UPGMA, or Mashtree.

**Architecture:** Strategy-based dispatch in `core/phylogeny/tree_builder.py` with a new Mashtree external tool wrapper. A new `cli/tree.py` module provides the Typer CLI. Existing `ani_to_newick()` is preserved unchanged for backward compatibility.

**Tech Stack:** BioPython (DistanceTreeConstructor), Typer (CLI), Rich (progress), Mashtree (external binary, optional)

---

### Task 1: Add TreeMethod Enum and UPGMA Builder

**Files:**
- Modify: `src/metadarkmatter/core/phylogeny/tree_builder.py`
- Test: `tests/unit/test_phylogeny_tree_builder.py`

**Step 1: Write the failing tests**

Add to `tests/unit/test_phylogeny_tree_builder.py`:

```python
class TestTreeMethod:
    """Test TreeMethod enum."""

    def test_enum_values(self) -> None:
        """TreeMethod has nj, upgma, mashtree."""
        from metadarkmatter.core.phylogeny.tree_builder import TreeMethod

        assert TreeMethod.NJ == "nj"
        assert TreeMethod.UPGMA == "upgma"
        assert TreeMethod.MASHTREE == "mashtree"


class TestAniToUpgma:
    """Test ANI matrix to UPGMA tree conversion."""

    def test_basic_3x3_matrix(self) -> None:
        """3x3 ANI matrix produces valid UPGMA Newick with all taxa."""
        from metadarkmatter.core.phylogeny.tree_builder import ani_to_upgma

        ani = pd.DataFrame(
            {
                "A": [100.0, 95.0, 80.0],
                "B": [95.0, 100.0, 82.0],
                "C": [80.0, 82.0, 100.0],
            },
            index=["A", "B", "C"],
        )

        newick = ani_to_upgma(ani)

        assert newick is not None
        assert newick.endswith(";")
        tree = Phylo.read(StringIO(newick), "newick")
        tip_names = {t.name for t in tree.get_terminals()}
        assert tip_names == {"A", "B", "C"}

    def test_fewer_than_3_genomes_returns_none(self) -> None:
        """Fewer than 3 genomes returns None."""
        from metadarkmatter.core.phylogeny.tree_builder import ani_to_upgma

        ani = pd.DataFrame(
            {"A": [100.0, 95.0], "B": [95.0, 100.0]},
            index=["A", "B"],
        )
        assert ani_to_upgma(ani) is None

    def test_empty_matrix_returns_none(self) -> None:
        """Empty matrix returns None."""
        from metadarkmatter.core.phylogeny.tree_builder import ani_to_upgma

        assert ani_to_upgma(pd.DataFrame()) is None

    def test_upgma_produces_ultrametric_tree(self) -> None:
        """UPGMA should produce a tree where root-to-tip distances are equal."""
        from metadarkmatter.core.phylogeny.tree_builder import ani_to_upgma

        ani = pd.DataFrame(
            {
                "A": [100.0, 99.0, 80.0],
                "B": [99.0, 100.0, 81.0],
                "C": [80.0, 81.0, 100.0],
            },
            index=["A", "B", "C"],
        )

        newick = ani_to_upgma(ani)
        assert newick is not None

        tree = Phylo.read(StringIO(newick), "newick")
        # UPGMA trees are ultrametric: all root-to-tip distances should be approx equal
        distances = [tree.distance(t) for t in tree.get_terminals()]
        assert max(distances) - min(distances) < 1.0  # Within 1 unit tolerance

    def test_missing_values_handled(self) -> None:
        """Missing ANI values should be filled and produce valid tree."""
        from metadarkmatter.core.phylogeny.tree_builder import ani_to_upgma

        ani = pd.DataFrame(
            {
                "A": [100.0, 95.0, np.nan],
                "B": [95.0, 100.0, 82.0],
                "C": [np.nan, 82.0, 100.0],
            },
            index=["A", "B", "C"],
        )

        newick = ani_to_upgma(ani)
        assert newick is not None
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/unit/test_phylogeny_tree_builder.py::TestTreeMethod -v`
Expected: FAIL with `ImportError` (TreeMethod not defined)

Run: `python -m pytest tests/unit/test_phylogeny_tree_builder.py::TestAniToUpgma -v`
Expected: FAIL with `ImportError` (ani_to_upgma not defined)

**Step 3: Implement TreeMethod enum and ani_to_upgma**

In `src/metadarkmatter/core/phylogeny/tree_builder.py`, add after the imports:

```python
from enum import Enum

class TreeMethod(str, Enum):
    """Phylogenetic tree building method."""
    NJ = "nj"
    UPGMA = "upgma"
    MASHTREE = "mashtree"
```

Add the `ani_to_upgma` function after `ani_to_newick`:

```python
def ani_to_upgma(ani_matrix: pd.DataFrame) -> str | None:
    """Convert ANI matrix to UPGMA tree in Newick format.

    Uses BioPython's UPGMA implementation to construct an ultrametric
    phylogenetic tree from ANI values. UPGMA assumes a constant rate
    of evolution (molecular clock hypothesis).

    Args:
        ani_matrix: Square DataFrame with ANI values (0-100 scale) where
            rows and columns are genome identifiers. Diagonal should be 100.

    Returns:
        Newick format string if tree construction succeeds, None if the
        matrix has fewer than 3 genomes.

    Note:
        Missing ANI values (NaN) are filled with maximum distance (50).
    """
    from Bio import Phylo
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

    if len(ani_matrix) < 3:
        if len(ani_matrix) > 0:
            logger.warning(
                f"Too few genomes for phylogenetic tree (need >= 3). Got {len(ani_matrix)}."
            )
        return None

    # Convert ANI (similarity) to distance: distance = 100 - ANI
    distance_matrix = 100 - ani_matrix
    distance_matrix = distance_matrix.clip(lower=0, upper=50)

    # Handle missing values
    missing_count = distance_matrix.isna().sum().sum() // 2
    if missing_count > 0:
        logger.warning(f"{missing_count} genome pairs lack ANI values; using maximum distance (50)")
        distance_matrix = distance_matrix.fillna(50.0)

    # Build BioPython DistanceMatrix
    names = list(ani_matrix.columns)
    matrix = _to_lower_triangular(distance_matrix)
    dm = DistanceMatrix(names, matrix)

    # Construct UPGMA tree
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)

    # Convert to Newick string
    output = StringIO()
    Phylo.write(tree, output, "newick")
    return output.getvalue().strip()
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/unit/test_phylogeny_tree_builder.py::TestTreeMethod tests/unit/test_phylogeny_tree_builder.py::TestAniToUpgma -v`
Expected: All PASS

**Step 5: Run existing tree builder tests to verify no regressions**

Run: `python -m pytest tests/unit/test_phylogeny_tree_builder.py -v`
Expected: All PASS (existing tests + new tests)

**Step 6: Commit**

```bash
git add src/metadarkmatter/core/phylogeny/tree_builder.py tests/unit/test_phylogeny_tree_builder.py
git commit -m "feat: add TreeMethod enum and UPGMA tree builder"
```

---

### Task 2: Add Mashtree External Tool Wrapper

**Files:**
- Create: `src/metadarkmatter/external/mashtree.py`
- Modify: `src/metadarkmatter/external/__init__.py`
- Test: `tests/unit/test_mashtree_wrapper.py`

**Step 1: Write the failing tests**

Create `tests/unit/test_mashtree_wrapper.py`:

```python
"""Tests for Mashtree external tool wrapper."""
from __future__ import annotations

from pathlib import Path

import pytest

from metadarkmatter.external.mashtree import Mashtree


class TestMashtreeCommand:
    """Test Mashtree command building."""

    @pytest.fixture(autouse=True)
    def mock_executable(self) -> None:
        """Mock the executable resolver so tests don't need mashtree installed."""
        Mashtree.set_executable_resolver(lambda name: f"/usr/bin/{name}")
        Mashtree.clear_cache()
        yield
        Mashtree.reset_executable_resolver()

    def test_build_command_basic(self, tmp_path: Path) -> None:
        """Basic command builds correctly."""
        g1 = tmp_path / "genome1.fna"
        g2 = tmp_path / "genome2.fna"
        g1.touch()
        g2.touch()

        mt = Mashtree()
        cmd = mt.build_command(genomes=[g1, g2], threads=4)

        assert cmd[0] == "/usr/bin/mashtree"
        assert "--numcpus" in cmd
        assert "4" in cmd
        assert str(g1) in cmd
        assert str(g2) in cmd

    def test_build_command_custom_threads(self, tmp_path: Path) -> None:
        """Threads parameter is passed correctly."""
        g1 = tmp_path / "genome1.fna"
        g1.touch()

        mt = Mashtree()
        cmd = mt.build_command(genomes=[g1], threads=16)

        idx = cmd.index("--numcpus")
        assert cmd[idx + 1] == "16"

    def test_tool_name(self) -> None:
        """Tool name is 'mashtree'."""
        assert Mashtree.TOOL_NAME == "mashtree"

    def test_install_hint(self) -> None:
        """Install hint includes conda."""
        assert "conda" in Mashtree.INSTALL_HINT

    def test_check_available_with_mock(self) -> None:
        """check_available returns True when mock resolver finds it."""
        assert Mashtree.check_available() is True

    def test_check_available_not_installed(self) -> None:
        """check_available returns False when not installed."""
        Mashtree.set_executable_resolver(lambda name: None)
        Mashtree.clear_cache()
        assert Mashtree.check_available() is False
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/unit/test_mashtree_wrapper.py -v`
Expected: FAIL with `ModuleNotFoundError` (mashtree module doesn't exist)

**Step 3: Implement the Mashtree wrapper**

Create `src/metadarkmatter/external/mashtree.py`:

```python
"""
Mashtree wrapper for building phylogenetic trees from genome assemblies.

Mashtree uses Mash (MinHash) distances between genomes to construct a
neighbor-joining tree. It is fast and suitable for quick phylogenetic
overviews of genome collections.
"""
from __future__ import annotations

from pathlib import Path
from typing import ClassVar

from metadarkmatter.external.base import ExternalTool


class Mashtree(ExternalTool):
    """Wrapper for mashtree phylogenetic tree builder.

    Mashtree builds a tree from genome assemblies using Mash distances.
    It takes FASTA files as input and writes a Newick tree to stdout.

    Example:
        >>> mt = Mashtree()
        >>> result = mt.run(genomes=[Path("g1.fna"), Path("g2.fna")], threads=4)
        >>> newick = result.stdout.strip()
    """

    TOOL_NAME: ClassVar[str] = "mashtree"
    INSTALL_HINT: ClassVar[str] = "conda install -c bioconda mashtree"

    def build_command(
        self,
        *,
        genomes: list[Path],
        threads: int = 4,
    ) -> list[str]:
        """Build mashtree command.

        Mashtree reads genome FASTA files and writes a Newick tree to stdout.

        Args:
            genomes: List of paths to genome FASTA files.
            threads: Number of CPU threads for parallel execution.

        Returns:
            Command as list of strings.
        """
        exe = str(self.get_executable())
        cmd = [exe, "--numcpus", str(threads)]
        cmd.extend(str(g) for g in genomes)
        return cmd
```

**Step 4: Update `__init__.py` exports**

In `src/metadarkmatter/external/__init__.py`, add the import and export:

```python
from metadarkmatter.external.mashtree import Mashtree
```

And add `"Mashtree"` to the `__all__` list.

**Step 5: Run tests to verify they pass**

Run: `python -m pytest tests/unit/test_mashtree_wrapper.py -v`
Expected: All PASS

**Step 6: Commit**

```bash
git add src/metadarkmatter/external/mashtree.py src/metadarkmatter/external/__init__.py tests/unit/test_mashtree_wrapper.py
git commit -m "feat: add Mashtree external tool wrapper"
```

---

### Task 3: Add mashtree_to_newick and build_tree Dispatcher

**Files:**
- Modify: `src/metadarkmatter/core/phylogeny/tree_builder.py`
- Test: `tests/unit/test_phylogeny_tree_builder.py`

**Step 1: Write the failing tests**

Add to `tests/unit/test_phylogeny_tree_builder.py`:

```python
class TestMashtreeToNewick:
    """Test mashtree_to_newick function."""

    def test_no_genome_files_raises(self, tmp_path: Path) -> None:
        """Empty directory should raise FileNotFoundError."""
        from metadarkmatter.core.phylogeny.tree_builder import mashtree_to_newick

        with pytest.raises(FileNotFoundError):
            mashtree_to_newick(tmp_path)

    def test_fewer_than_3_genomes_raises(self, tmp_path: Path) -> None:
        """Fewer than 3 genome files should raise ValueError."""
        from metadarkmatter.core.phylogeny.tree_builder import mashtree_to_newick

        (tmp_path / "g1.fna").write_text(">seq\nACGT")
        (tmp_path / "g2.fna").write_text(">seq\nACGT")

        with pytest.raises(ValueError, match="need >= 3"):
            mashtree_to_newick(tmp_path)


class TestBuildTree:
    """Test build_tree dispatcher."""

    def test_nj_method(self) -> None:
        """NJ method produces valid tree from ANI matrix."""
        from metadarkmatter.core.phylogeny.tree_builder import TreeMethod, build_tree

        ani = pd.DataFrame(
            {
                "A": [100.0, 95.0, 80.0],
                "B": [95.0, 100.0, 82.0],
                "C": [80.0, 82.0, 100.0],
            },
            index=["A", "B", "C"],
        )

        newick = build_tree(TreeMethod.NJ, ani_matrix=ani)
        assert newick is not None
        assert "A" in newick
        tree = Phylo.read(StringIO(newick), "newick")
        assert len(list(tree.get_terminals())) == 3

    def test_upgma_method(self) -> None:
        """UPGMA method produces valid tree from ANI matrix."""
        from metadarkmatter.core.phylogeny.tree_builder import TreeMethod, build_tree

        ani = pd.DataFrame(
            {
                "A": [100.0, 95.0, 80.0],
                "B": [95.0, 100.0, 82.0],
                "C": [80.0, 82.0, 100.0],
            },
            index=["A", "B", "C"],
        )

        newick = build_tree(TreeMethod.UPGMA, ani_matrix=ani)
        assert newick is not None
        tree = Phylo.read(StringIO(newick), "newick")
        assert len(list(tree.get_terminals())) == 3

    def test_nj_without_ani_raises(self) -> None:
        """NJ method without ANI matrix should raise ValueError."""
        from metadarkmatter.core.phylogeny.tree_builder import TreeMethod, build_tree

        with pytest.raises(ValueError, match="ANI matrix required"):
            build_tree(TreeMethod.NJ)

    def test_upgma_without_ani_raises(self) -> None:
        """UPGMA method without ANI matrix should raise ValueError."""
        from metadarkmatter.core.phylogeny.tree_builder import TreeMethod, build_tree

        with pytest.raises(ValueError, match="ANI matrix required"):
            build_tree(TreeMethod.UPGMA)

    def test_mashtree_without_genomes_raises(self) -> None:
        """Mashtree method without genome_dir should raise ValueError."""
        from metadarkmatter.core.phylogeny.tree_builder import TreeMethod, build_tree

        with pytest.raises(ValueError, match="genome_dir required"):
            build_tree(TreeMethod.MASHTREE)
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/unit/test_phylogeny_tree_builder.py::TestMashtreeToNewick tests/unit/test_phylogeny_tree_builder.py::TestBuildTree -v`
Expected: FAIL with `ImportError`

**Step 3: Implement mashtree_to_newick and build_tree**

Add to `src/metadarkmatter/core/phylogeny/tree_builder.py`:

```python
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
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/unit/test_phylogeny_tree_builder.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
git add src/metadarkmatter/core/phylogeny/tree_builder.py tests/unit/test_phylogeny_tree_builder.py
git commit -m "feat: add mashtree_to_newick and build_tree dispatcher"
```

---

### Task 4: Create CLI Command

**Files:**
- Create: `src/metadarkmatter/cli/tree.py`
- Modify: `src/metadarkmatter/cli/main.py`
- Test: `tests/unit/test_cli_tree.py`

**Step 1: Write the failing tests**

Create `tests/unit/test_cli_tree.py`:

```python
"""Tests for the tree CLI command."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.tree import app

runner = CliRunner()


@pytest.fixture()
def ani_csv(tmp_path: Path) -> Path:
    """Create a minimal ANI matrix CSV."""
    ani = pd.DataFrame(
        {
            "A": [100.0, 95.0, 80.0],
            "B": [95.0, 100.0, 82.0],
            "C": [80.0, 82.0, 100.0],
        },
        index=["A", "B", "C"],
    )
    path = tmp_path / "ani_matrix.csv"
    # Write with index as first column (metadarkmatter format)
    ani.to_csv(path)
    return path


class TestTreeBuildNJ:
    """Test tree build with NJ method."""

    def test_nj_produces_newick(self, ani_csv: Path, tmp_path: Path) -> None:
        """NJ method produces valid Newick output file."""
        output = tmp_path / "tree.nwk"
        result = runner.invoke(app, [
            "build",
            "--method", "nj",
            "--ani", str(ani_csv),
            "--output", str(output),
        ])
        assert result.exit_code == 0, result.output
        assert output.exists()
        content = output.read_text().strip()
        assert content.endswith(";")
        assert "A" in content

    def test_nj_without_ani_fails(self, tmp_path: Path) -> None:
        """NJ without --ani should fail."""
        output = tmp_path / "tree.nwk"
        result = runner.invoke(app, [
            "build",
            "--method", "nj",
            "--output", str(output),
        ])
        assert result.exit_code != 0


class TestTreeBuildUPGMA:
    """Test tree build with UPGMA method."""

    def test_upgma_produces_newick(self, ani_csv: Path, tmp_path: Path) -> None:
        """UPGMA method produces valid Newick output file."""
        output = tmp_path / "tree.nwk"
        result = runner.invoke(app, [
            "build",
            "--method", "upgma",
            "--ani", str(ani_csv),
            "--output", str(output),
        ])
        assert result.exit_code == 0, result.output
        assert output.exists()
        content = output.read_text().strip()
        assert content.endswith(";")


class TestTreeBuildMashtree:
    """Test tree build with Mashtree method."""

    def test_mashtree_without_genomes_fails(self, tmp_path: Path) -> None:
        """Mashtree without --genomes should fail."""
        output = tmp_path / "tree.nwk"
        result = runner.invoke(app, [
            "build",
            "--method", "mashtree",
            "--output", str(output),
        ])
        assert result.exit_code != 0


class TestInputValidation:
    """Test CLI input validation."""

    def test_nonexistent_ani_fails(self, tmp_path: Path) -> None:
        """Non-existent ANI file should fail."""
        output = tmp_path / "tree.nwk"
        result = runner.invoke(app, [
            "build",
            "--method", "nj",
            "--ani", str(tmp_path / "missing.csv"),
            "--output", str(output),
        ])
        assert result.exit_code != 0

    def test_too_few_genomes_in_ani(self, tmp_path: Path) -> None:
        """ANI matrix with < 3 genomes should fail gracefully."""
        ani = pd.DataFrame(
            {"A": [100.0, 95.0], "B": [95.0, 100.0]},
            index=["A", "B"],
        )
        ani_path = tmp_path / "small_ani.csv"
        ani.to_csv(ani_path)

        output = tmp_path / "tree.nwk"
        result = runner.invoke(app, [
            "build",
            "--method", "nj",
            "--ani", str(ani_path),
            "--output", str(output),
        ])
        assert result.exit_code != 0
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/unit/test_cli_tree.py -v`
Expected: FAIL with `ModuleNotFoundError` (cli/tree.py doesn't exist)

**Step 3: Implement the CLI command**

Create `src/metadarkmatter/cli/tree.py`:

```python
"""
Tree command for building phylogenetic trees.

Provides subcommands:
- build: Build a phylogenetic tree from ANI matrix or genome assemblies
"""
from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console

from metadarkmatter.cli.utils import QuietConsole, spinner_progress
from metadarkmatter.core.phylogeny.tree_builder import TreeMethod

app = typer.Typer(
    name="tree",
    help="Build phylogenetic trees from genomes or ANI matrices",
    no_args_is_help=True,
)

console = Console()


@app.command(name="build")
def build(
    method: TreeMethod = typer.Option(
        ...,
        "--method",
        "-m",
        help="Tree building method: nj, upgma, or mashtree",
    ),
    ani: Path | None = typer.Option(
        None,
        "--ani",
        "-a",
        help="ANI matrix CSV file (required for nj/upgma methods)",
        exists=True,
        dir_okay=False,
    ),
    genomes: Path | None = typer.Option(
        None,
        "--genomes",
        "-g",
        help="Directory containing genome FASTA files (required for mashtree method)",
        exists=True,
        file_okay=False,
    ),
    output: Path = typer.Option(
        ...,
        "--output",
        "-o",
        help="Output Newick tree file",
    ),
    threads: int = typer.Option(
        4,
        "--threads",
        "-t",
        help="Number of threads (mashtree only)",
        min=1,
    ),
    genome_pattern: str = typer.Option(
        "*.fna",
        "--genome-pattern",
        "-p",
        help="Glob pattern for genome files (mashtree only)",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Enable verbose output",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet",
        "-q",
        help="Suppress progress output",
    ),
) -> None:
    """
    Build a phylogenetic tree.

    Supports three methods:

    - nj: Neighbor-joining from ANI distance matrix (BioPython)

    - upgma: UPGMA from ANI distance matrix (BioPython, ultrametric)

    - mashtree: Mash distance-based NJ from genome assemblies (requires mashtree)

    Examples:

        # NJ tree from ANI matrix
        metadarkmatter tree build --method nj --ani ani_matrix.csv --output tree.nwk

        # UPGMA tree
        metadarkmatter tree build --method upgma --ani ani_matrix.csv --output tree.nwk

        # Mashtree from genome files
        metadarkmatter tree build --method mashtree --genomes genomes/ --output tree.nwk -t 16
    """
    import polars as pl

    from metadarkmatter.core.phylogeny.tree_builder import build_tree

    out = QuietConsole(console, quiet=quiet)

    out.print("\n[bold blue]Metadarkmatter Tree Builder[/bold blue]\n")
    out.print(f"[bold]Method:[/bold] {method.value}")

    # Validate inputs
    if method in (TreeMethod.NJ, TreeMethod.UPGMA):
        if ani is None:
            console.print("[red]Error: --ani is required for nj/upgma methods[/red]")
            raise typer.Exit(code=1) from None

        # Load ANI matrix
        out.print(f"[bold]ANI matrix:[/bold] {ani}")
        try:
            ani_df = pl.read_csv(ani)
            # First column is the genome names (index)
            genome_names = ani_df.get_column(ani_df.columns[0]).to_list()
            ani_values = ani_df.select(ani_df.columns[1:])
            import pandas as pd
            ani_matrix = pd.DataFrame(
                ani_values.to_numpy(),
                index=genome_names,
                columns=genome_names,
            )
        except Exception as e:
            console.print(f"[red]Error loading ANI matrix: {e}[/red]")
            raise typer.Exit(code=1) from None

        n_genomes = len(ani_matrix)
        out.print(f"[bold]Genomes:[/bold] {n_genomes}")

        if n_genomes < 3:
            console.print(
                f"[red]Error: Need >= 3 genomes for tree building (got {n_genomes})[/red]"
            )
            raise typer.Exit(code=1) from None

        with spinner_progress(
            f"Building {method.value.upper()} tree from {n_genomes} genomes...",
            console,
            quiet,
        ):
            newick = build_tree(method, ani_matrix=ani_matrix)

    elif method == TreeMethod.MASHTREE:
        if genomes is None:
            console.print("[red]Error: --genomes is required for mashtree method[/red]")
            raise typer.Exit(code=1) from None

        out.print(f"[bold]Genomes:[/bold] {genomes}")
        out.print(f"[bold]Threads:[/bold] {threads}")

        with spinner_progress(
            f"Building Mashtree tree from {genomes}...",
            console,
            quiet,
        ):
            newick = build_tree(
                method,
                genome_dir=genomes,
                genome_pattern=genome_pattern,
                threads=threads,
            )
    else:
        console.print(f"[red]Error: Unknown method '{method}'[/red]")
        raise typer.Exit(code=1) from None

    if newick is None:
        console.print("[red]Error: Tree construction failed[/red]")
        raise typer.Exit(code=1) from None

    # Write output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(newick + "\n")

    out.print(f"\n[bold green]Tree built successfully![/bold green]")
    out.print(f"[bold]Output:[/bold] {output}")
    out.print(f"\n[dim]Use with: metadarkmatter report generate --tree {output}[/dim]")
    out.print()
```

**Step 4: Register in main.py**

In `src/metadarkmatter/cli/main.py`, add the import and registration:

After line 58 (the existing imports line), add `tree` to the import:
```python
from metadarkmatter.cli import aai, ani, blast, blastx, download, kraken2, mapping, mmseqs2, proteins, report, score, tree, visualize
```

After line 74 (the last `app.add_typer` line), add:
```python
app.add_typer(tree.app, name="tree")
```

**Step 5: Run tests to verify they pass**

Run: `python -m pytest tests/unit/test_cli_tree.py -v`
Expected: All PASS

**Step 6: Run full test suite for regressions**

Run: `python -m pytest tests/ -x -q`
Expected: All pass

**Step 7: Commit**

```bash
git add src/metadarkmatter/cli/tree.py src/metadarkmatter/cli/main.py tests/unit/test_cli_tree.py
git commit -m "feat: add metadarkmatter tree build CLI command"
```

---

### Task 5: End-to-End Verification

**Step 1: Run the full test suite**

Run: `python -m pytest tests/ -x -q`
Expected: All pass with 0 failures

**Step 2: Verify CLI help output**

Run: `python -m metadarkmatter tree --help`
Expected: Shows "Build phylogenetic trees from genomes or ANI matrices" and the `build` subcommand

Run: `python -m metadarkmatter tree build --help`
Expected: Shows all options (--method, --ani, --genomes, --output, --threads, etc.)

**Step 3: Test NJ tree build from CLI (functional)**

If an ANI matrix CSV exists in the project's test fixtures or can be created quickly:

```bash
# Create a small test matrix
python -c "
import pandas as pd
ani = pd.DataFrame({
    'A': [100.0, 95.0, 80.0, 75.0],
    'B': [95.0, 100.0, 82.0, 77.0],
    'C': [80.0, 82.0, 100.0, 90.0],
    'D': [75.0, 77.0, 90.0, 100.0],
}, index=['A', 'B', 'C', 'D'])
ani.to_csv('/tmp/test_ani.csv')
"

python -m metadarkmatter tree build --method nj --ani /tmp/test_ani.csv --output /tmp/test_tree.nwk
cat /tmp/test_tree.nwk
```

Expected: Valid Newick tree printed, file created

**Step 4: Test UPGMA tree build**

```bash
python -m metadarkmatter tree build --method upgma --ani /tmp/test_ani.csv --output /tmp/test_upgma.nwk
cat /tmp/test_upgma.nwk
```

Expected: Valid Newick tree (may differ from NJ)

**Step 5: Commit (if any fixes needed)**

Only if fixes were made during verification.
