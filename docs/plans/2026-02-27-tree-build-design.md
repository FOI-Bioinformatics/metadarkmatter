# Phylogenetic Tree Build Command

**Date:** 2026-02-27
**Status:** Approved

## Problem

Phylogenetic tree construction is currently embedded in the report pipeline (`ani_to_newick()` in `tree_builder.py`), supporting only Neighbor-Joining from ANI matrices. Users cannot build standalone trees, choose alternative methods, or use Mash-based distance estimation for rapid tree construction from genome assemblies.

## Design Goals

- Expose tree building as a standalone CLI command (`metadarkmatter tree build`)
- Support three methods: NJ (ANI), UPGMA (ANI), and Mashtree (genome FASTAs)
- Follow existing external tool wrapper patterns for Mashtree
- Maintain backward compatibility with the report pipeline
- Clean strategy-based architecture for extensibility

## CLI Interface

```bash
# Neighbor-joining from ANI matrix (default)
metadarkmatter tree build --method nj --ani ani_matrix.csv --output tree.nwk

# UPGMA from ANI matrix
metadarkmatter tree build --method upgma --ani ani_matrix.csv --output tree.nwk

# Mashtree from genome FASTA files
metadarkmatter tree build --method mashtree --genomes genomes/ --output tree.nwk --threads 16

# Use with report
metadarkmatter report generate --classifications results.csv --ani ani.csv --tree tree.nwk --output report.html
```

### Parameters

| Parameter | Required | Methods | Description |
|-----------|----------|---------|-------------|
| `--method` | Yes | all | Tree-building method: `nj`, `upgma`, `mashtree` |
| `--ani` | NJ/UPGMA | nj, upgma | Path to ANI matrix CSV |
| `--genomes` | Mashtree | mashtree | Directory containing genome FASTA files |
| `--output` | Yes | all | Output Newick file path |
| `--threads` | No | mashtree | Number of threads (default: 4) |
| `--genome-pattern` | No | mashtree | Glob pattern for genome files (default: `*.fna`) |
| `--verbose` | No | all | Enable verbose output |
| `--quiet` | No | all | Suppress progress output |

### Input Validation

- NJ/UPGMA: `--ani` required, `--genomes` ignored with warning
- Mashtree: `--genomes` required, `--ani` ignored with warning
- NJ/UPGMA: ANI matrix must have >= 3 genomes
- Mashtree: genome directory must contain >= 3 genome files

## Architecture

### Component 1: TreeMethod Enum and Dispatcher

**File:** `src/metadarkmatter/core/phylogeny/tree_builder.py`

Extends the existing module with:

```python
class TreeMethod(str, Enum):
    NJ = "nj"
    UPGMA = "upgma"
    MASHTREE = "mashtree"

def build_tree(
    method: TreeMethod,
    *,
    ani_matrix: pd.DataFrame | None = None,
    genome_dir: Path | None = None,
    genome_pattern: str = "*.fna",
    threads: int = 4,
) -> str:
    """Build a phylogenetic tree using the specified method.

    Returns Newick format string.
    """
```

The dispatcher validates inputs per method and calls the appropriate builder function.

### Component 2: UPGMA Builder

**File:** `src/metadarkmatter/core/phylogeny/tree_builder.py`

```python
def ani_to_upgma(ani_matrix: pd.DataFrame) -> str | None:
    """Convert ANI matrix to UPGMA tree in Newick format.

    Uses BioPython's UPGMA implementation. Same distance conversion
    as ani_to_newick() but produces an ultrametric tree.
    """
```

Shares the existing `_to_lower_triangular()` helper and ANI-to-distance conversion logic with `ani_to_newick()`.

### Component 3: Mashtree External Tool Wrapper

**File:** `src/metadarkmatter/external/mashtree.py` (new)

```python
class Mashtree(ExternalTool):
    TOOL_NAME = "mashtree"
    INSTALL_HINT = "conda install -c bioconda mashtree"

    def build_command(
        self,
        *,
        genomes: list[Path],
        threads: int = 4,
    ) -> list[str]:
        """Build mashtree command.

        mashtree --numcpus {threads} {genome_files...}

        Mashtree writes Newick to stdout.
        """
```

Follows the exact same pattern as `FastANI` and `Skani` wrappers.

### Component 4: Mashtree Integration in tree_builder

**File:** `src/metadarkmatter/core/phylogeny/tree_builder.py`

```python
def mashtree_to_newick(
    genome_dir: Path,
    genome_pattern: str = "*.fna",
    threads: int = 4,
) -> str:
    """Build tree from genome FASTAs using Mashtree.

    Returns Newick format string from mashtree stdout.
    """
```

### Component 5: CLI Command

**File:** `src/metadarkmatter/cli/tree.py` (new)

```python
app = typer.Typer(
    name="tree",
    help="Build phylogenetic trees from genomes or ANI matrices",
    no_args_is_help=True,
)

@app.command(name="build")
def build(
    method: TreeMethod,
    ani: Path | None,
    genomes: Path | None,
    output: Path,
    threads: int = 4,
    genome_pattern: str = "*.fna",
    verbose: bool = False,
    quiet: bool = False,
) -> None:
    """Build a phylogenetic tree."""
```

Rich progress spinner during computation. Summary output with genome count, method, output path.

### Component 6: Registration

**File:** `src/metadarkmatter/cli/main.py`

```python
from metadarkmatter.cli import tree
app.add_typer(tree.app, name="tree")
```

## What Stays Unchanged

- `ani_to_newick()` function signature and behavior (report pipeline calls it)
- `load_user_tree()` and `_prune_tips()`
- Novel cluster placement logic (`placement.py`)
- Report `--tree` flag behavior (accepts any Newick file)
- All existing CLI commands

## Files to Create/Modify

| File | Action |
|------|--------|
| `external/mashtree.py` | **New** - Mashtree ExternalTool wrapper |
| `external/__init__.py` | Add Mashtree export |
| `core/phylogeny/tree_builder.py` | Add `TreeMethod`, `build_tree()`, `ani_to_upgma()`, `mashtree_to_newick()` |
| `cli/tree.py` | **New** - `tree build` command |
| `cli/main.py` | Register tree app |
| Tests | Unit tests for UPGMA, Mashtree wrapper, build_tree dispatcher, CLI |
