# Codebase Simplification Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Remove ~3,000 lines of dead code and duplication, consolidate to a single classifier implementation.

**Architecture:** Four phases: (1) safe dead code removal, (2) classifier consolidation routing all CLI through VectorizedClassifier, (3) duplication reduction across heatmaps/thresholds/streaming, (4) structural cleanup. Each phase ends with full test suite verification.

**Tech Stack:** Python 3.11+, Polars, Pydantic v2, Typer, Plotly, pytest

---

### Task 1: Delete backup file and unused constants

**Files:**
- Delete: `src/metadarkmatter/core/ani_placement.py.backup`
- Modify: `src/metadarkmatter/core/constants.py`

**Step 1: Delete the backup file**

```bash
rm src/metadarkmatter/core/ani_placement.py.backup
```

**Step 2: Remove unused constants from constants.py**

Remove these constants that are never imported outside the file (lines 18-21, 32, 35-36, 42, 57, 60, 66, 117-118):

```python
# DELETE these lines:
BLAST_OUTFMT_12COL = ...        # line ~18
ANI_GENUS_BOUNDARY = 75.0       # line 32
ANI_SPECIES_BOUNDARY_LOW = 95.0  # line 35
ANI_SPECIES_BOUNDARY_HIGH = 96.0 # line 36
ANI_MIN_RELATED = ...           # line 42
AAI_GENUS_BOUNDARY_HIGH = 65.0  # line 57
AAI_GENUS_BOUNDARY_LOW = 58.0   # line 60
AAI_MIN_RELATED = ...           # line 66
UNCERTAINTY_AMBIGUOUS_MIN = ... # line 117
UNCERTAINTY_AMBIGUOUS_MAX = ... # line 118
```

Keep `ANI_DEFAULT_UNRELATED` (line 39) and `AAI_DEFAULT_UNRELATED` (line 63) - they ARE used.
Keep `CATEGORY_OFF_TARGET` (line 81) - used by family validation.

**Step 3: Run tests**

Run: `pytest tests/ --ignore=tests/unit/test_phylogeny_placement.py --ignore=tests/unit/test_phylogeny_tree_builder.py -x -q`
Expected: All pass (these constants are never imported)

**Step 4: Commit**

```bash
git add -A
git commit -m "chore: remove backup file and unused constants"
```

---

### Task 2: Delete ParallelClassifier and SparseANIMatrix

**Files:**
- Delete: `src/metadarkmatter/core/classification/classifiers/parallel.py`
- Delete: `src/metadarkmatter/core/classification/sparse_ani_matrix.py`
- Modify: `src/metadarkmatter/core/classification/classifiers/__init__.py`
- Modify: `src/metadarkmatter/core/classification/__init__.py`
- Modify: `src/metadarkmatter/core/ani_placement.py`
- Modify: `src/metadarkmatter/core/__init__.py`
- Modify: `src/metadarkmatter/__init__.py`
- Modify: `tests/unit/test_classifier_equivalence.py`
- Modify: `tests/unit/test_ani_placement.py`

**Step 1: Delete the files**

```bash
rm src/metadarkmatter/core/classification/classifiers/parallel.py
rm src/metadarkmatter/core/classification/sparse_ani_matrix.py
```

**Step 2: Update classifiers/__init__.py**

Remove `ParallelClassifier` from imports and `__all__`:

```python
# Keep only:
from metadarkmatter.core.classification.classifiers.base import ANIWeightedClassifier
from metadarkmatter.core.classification.classifiers.vectorized import VectorizedClassifier

__all__ = [
    "ANIWeightedClassifier",
    "VectorizedClassifier",
]
```

**Step 3: Update classification/__init__.py**

Remove `SparseANIMatrix` and `ParallelClassifier` from imports and `__all__`.

**Step 4: Update core/ani_placement.py**

Remove `ParallelClassifier` and `SparseANIMatrix` from imports and `__all__`.

**Step 5: Update core/__init__.py**

Remove `ParallelClassifier` and `SparseANIMatrix` from imports and `__all__`.

**Step 6: Update root __init__.py**

Remove `ParallelClassifier` and `SparseANIMatrix` if present in imports and `__all__`.

**Step 7: Update test_classifier_equivalence.py**

Remove or skip the `ParallelClassifier` tests. The test file tests equivalence between classifiers - remove the parallel classifier portion. If the entire file only tests ParallelClassifier equivalence, delete it.

**Step 8: Update test_ani_placement.py**

Remove `TestSparseANIMatrix` class (lines 147-201) and the `SparseANIMatrix` import (line 22).

**Step 9: Run tests**

Run: `pytest tests/ --ignore=tests/unit/test_phylogeny_placement.py --ignore=tests/unit/test_phylogeny_tree_builder.py -x -q`
Expected: All pass

**Step 10: Commit**

```bash
git add -A
git commit -m "chore: remove ParallelClassifier and SparseANIMatrix (unused from CLI)"
```

---

### Task 3: Remove unused imports and add missing imports

**Files:**
- Modify: `src/metadarkmatter/core/classification/classifiers/vectorized.py`
- Modify: `src/metadarkmatter/core/classification/classifiers/base.py`

**Step 1: Remove unused imports from vectorized.py**

Remove these unused constant imports (lines ~18-22):
- `NOVELTY_KNOWN_MAX`
- `NOVELTY_NOVEL_GENUS_MIN`
- `NOVELTY_NOVEL_GENUS_MAX`
- `NOVELTY_NOVEL_SPECIES_MAX`
- `NOVELTY_NOVEL_SPECIES_MIN`

Keep `calculate_confidence_score` and any other actually-used imports.

**Step 2: Add missing Callable import to base.py**

Add `from collections.abc import Callable` at the top of `base.py` (near other imports). The file uses `Callable` in type hints but never imports it.

**Step 3: Run tests**

Run: `pytest tests/ --ignore=tests/unit/test_phylogeny_placement.py --ignore=tests/unit/test_phylogeny_tree_builder.py -x -q`
Expected: All pass

**Step 4: Commit**

```bash
git add -A
git commit -m "chore: clean up imports in classifier files"
```

---

### Task 4: Remove backward-compat aliases from extended_matrix_builder

**Files:**
- Modify: `src/metadarkmatter/visualization/report/components/extended_matrix_builder.py`
- Modify: `src/metadarkmatter/visualization/report/components/__init__.py`
- Modify: `tests/unit/test_extended_matrix_builder.py`

**Step 1: Remove alias functions from extended_matrix_builder.py**

Remove these functions that just wrap the generic `*_similarity_*` counterparts:
- `build_extended_ani_matrix()` (wraps `build_extended_similarity_matrix`)
- `_build_ani_dict()` (wraps `_build_similarity_dict`)
- `estimate_novel_to_reference_ani()` (wraps `estimate_novel_to_reference_similarity`)
- `estimate_novel_to_novel_ani()` (wraps `estimate_novel_to_novel_similarity`)
- `get_cluster_label()` (never called from production)

**Step 2: Update components/__init__.py**

Remove the alias imports:
```python
# Remove these from imports:
build_extended_ani_matrix,
estimate_novel_to_novel_ani,
estimate_novel_to_reference_ani,
get_cluster_label,

# Remove from __all__:
"build_extended_ani_matrix",
"estimate_novel_to_novel_ani",
"estimate_novel_to_reference_ani",
"get_cluster_label",
```

**Step 3: Update test_extended_matrix_builder.py**

Update tests to call the generic `build_extended_similarity_matrix()` / `estimate_novel_to_reference_similarity()` / `estimate_novel_to_novel_similarity()` instead of the removed ANI aliases.

**Step 4: Run tests**

Run: `pytest tests/unit/test_extended_matrix_builder.py -v`
Expected: All pass

**Step 5: Commit**

```bash
git add -A
git commit -m "chore: remove backward-compat alias functions in extended_matrix_builder"
```

---

### Task 5: Remove GlobalConfig and other dead code

**Files:**
- Modify: `src/metadarkmatter/models/config.py`
- Modify: `src/metadarkmatter/models/__init__.py`
- Modify: `src/metadarkmatter/core/protein_constants.py`
- Modify: `src/metadarkmatter/models/blast.py`

**Step 1: Remove GlobalConfig from config.py**

Remove the `GlobalConfig` class (line 686 to end of file). It is never used by any CLI command.

Also remove any sub-config classes it references that are not used independently (check if `Bowtie2Config`, `MapConfig`, `KrakenConfig`, `ExtractConfig`, `VisualizeConfig` are used outside GlobalConfig - if they are used by CLI commands independently, keep them).

**Step 2: Update models/__init__.py**

Remove `GlobalConfig` from imports and `__all__`.

**Step 3: Remove calculate_protein_confidence_score from protein_constants.py**

Remove the thin wrapper function `calculate_protein_confidence_score()` (line ~97). It only delegates to `calculate_confidence_score()` with different defaults and is never called from production code.

**Step 4: Remove sorted_by_bitscore from blast.py**

Remove `BlastResult.sorted_by_bitscore()` method - it just returns `self` and is never meaningfully called.

**Step 5: Run tests**

Run: `pytest tests/ --ignore=tests/unit/test_phylogeny_placement.py --ignore=tests/unit/test_phylogeny_tree_builder.py -x -q`

If `test_protein_classification.py` fails due to removed function, update those tests to call `calculate_confidence_score()` directly with protein thresholds.

**Step 6: Commit**

```bash
git add -A
git commit -m "chore: remove GlobalConfig, protein confidence wrapper, and dead methods"
```

---

### Task 6: Phase 1 verification

**Step 1: Run full test suite**

Run: `pytest tests/ --ignore=tests/unit/test_phylogeny_placement.py --ignore=tests/unit/test_phylogeny_tree_builder.py -v --tb=short 2>&1 | tail -30`
Expected: Same pass count as before Phase 1 (minus removed tests for deleted code)

**Step 2: Count lines removed**

```bash
git diff --stat main
```

Expected: Net deletion of ~3,000+ lines

---

### Task 7: Route all CLI paths through VectorizedClassifier

This is the largest and most important task. The `--fast` and default (standard) modes in `cli/score.py` currently use `ANIWeightedClassifier`. We route them through `VectorizedClassifier` instead.

**Files:**
- Modify: `src/metadarkmatter/cli/score.py`

**Step 1: Simplify ProcessingMode enum**

Replace the 4-mode enum with a simpler flag:

```python
# Remove ProcessingMode enum entirely (lines 106-112)
# Remove validate_processing_modes() function (lines 115-154)
```

**Step 2: Simplify CLI classify function parameters**

In the `classify()` function:
- Remove `--fast` parameter (line ~446-449)
- Remove `--parallel` parameter (line ~451-455)
- Keep `--streaming` parameter
- Remove `validate_processing_modes()` call

**Step 3: Consolidate classification logic**

Replace the if/elif chain (lines ~945-1060) with unified VectorizedClassifier usage:

```python
# Old code had:
#   if streaming: VectorizedClassifier.stream_to_file()
#   elif parallel: VectorizedClassifier.classify_file()
#   elif fast: ANIWeightedClassifier.classify_to_dataframe_fast()
#   else: ANIWeightedClassifier.classify_to_dataframe()

# New code:
vectorized = VectorizedClassifier(ani_matrix=ani_matrix, aai_matrix=aai_matrix, config=config)

if streaming:
    # streaming path (unchanged)
    num_classified = vectorized.stream_to_file(...)
else:
    # unified path - VectorizedClassifier handles everything
    result = vectorized.classify_file(alignment)
    classification_df = result
    num_classified = _finalize_classification(...)
```

**Step 4: Update batch classify function**

Similarly update the batch `classify_batch()` function (lines ~1467-1495) to always use `VectorizedClassifier`.

**Step 5: Remove ANIWeightedClassifier import**

Remove the import of `ANIWeightedClassifier` from `cli/score.py` (line 35) if no longer needed.

**Step 6: Run tests**

Run: `pytest tests/ --ignore=tests/unit/test_phylogeny_placement.py --ignore=tests/unit/test_phylogeny_tree_builder.py -x -q`
Expected: All pass. CLI tests may need flag updates (remove --fast/--parallel test cases).

Check CLI help:
```bash
metadarkmatter score classify --help
```
Expected: No `--fast` or `--parallel` flags shown

**Step 7: Commit**

```bash
git add -A
git commit -m "refactor: route all CLI classification through VectorizedClassifier"
```

---

### Task 8: Remove dead code from base.py

Now that CLI no longer calls the fast/parallel paths, remove the dead code.

**Files:**
- Modify: `src/metadarkmatter/core/classification/classifiers/base.py`
- Modify: `tests/unit/test_ani_placement.py`

**Step 1: Remove module-level worker functions from base.py**

Delete these standalone functions (lines 1080-1636):
- `_classify_chunk_worker()` (~320 lines)
- `_calculate_inferred_uncertainty_inline()` (~45 lines)
- `_calculate_confidence_score_inline()` (~95 lines)
- `_calculate_identity_confidence_inline()` (~25 lines)
- `_calculate_placement_confidence_inline()` (~40 lines)
- `_calculate_discovery_score_inline()` (~30 lines)

**Step 2: Remove fast/batch methods from ANIWeightedClassifier**

Delete these methods from the class:
- `classify_read_fast()` (line 705)
- `classify_blast_file_fast()` (line 899)
- `classify_to_dataframe_fast()` (line 924)
- `stream_to_file_fast()` (line 978)
- `_write_chunk()` (line 1063)
- `classify_blast_file()` (line 617)
- `classify_to_dataframe()` (line 639)
- `write_classifications()` (line 677)

Keep:
- `__init__()` (line 55)
- `classify_read()` (line 81) - public programmatic API
- `_calculate_placement_uncertainty()` (line 298)
- `_calculate_genus_uncertainty()` (line 359)
- `_calculate_aai_uncertainty()` (line 410)
- `_classify_by_thresholds()` (line 459)

**Step 3: Remove unused imports from base.py**

Remove imports that are only needed by deleted code: `BlastResultFast`, `BlastHitFast`, multiprocessing imports, etc.

**Step 4: Update test_ani_placement.py**

Remove tests for deleted methods:
- `test_classify_blast_file` (line 408)
- `test_classify_to_dataframe_fast` (line 426)
- `test_classify_blast_file_fast` (line 740)

Keep tests for `classify_read()` as it remains the programmatic API.

**Step 5: Run tests**

Run: `pytest tests/ --ignore=tests/unit/test_phylogeny_placement.py --ignore=tests/unit/test_phylogeny_tree_builder.py -x -q`
Expected: All pass

**Step 6: Commit**

```bash
git add -A
git commit -m "refactor: remove dead classifier code from base.py (~1000 lines)"
```

---

### Task 9: Update re-export shim and Phase 2 verification

**Files:**
- Modify: `src/metadarkmatter/core/ani_placement.py`
- Modify: `src/metadarkmatter/cli/score.py` (if still importing from ani_placement)

**Step 1: Simplify ani_placement.py**

Update to only re-export what still exists:

```python
"""Backward compatibility re-exports."""
from __future__ import annotations

from metadarkmatter.core.classification import (
    ANIMatrix,
    ANIWeightedClassifier,
    VectorizedClassifier,
)

__all__ = [
    "ANIMatrix",
    "ANIWeightedClassifier",
    "VectorizedClassifier",
]
```

**Step 2: Run full test suite**

Run: `pytest tests/ --ignore=tests/unit/test_phylogeny_placement.py --ignore=tests/unit/test_phylogeny_tree_builder.py -v --tb=short 2>&1 | tail -30`
Expected: All pass

**Step 3: Commit**

```bash
git add -A
git commit -m "refactor: simplify re-export shim after classifier consolidation"
```

---

### Task 10: Unify ANI/AAI heatmap builders

**Files:**
- Modify: `src/metadarkmatter/visualization/report/components/heatmap_builder.py`

**Step 1: Create a shared _build_similarity_heatmap function**

The two functions `build_ani_heatmap()` (line 82) and `build_aai_heatmap()` (line 181) are ~100 lines each and differ only in:
- `default_value`: 70.0 vs 40.0
- `same_threshold`: 95.0 vs 65.0
- `boundary_lower`: 93.0 vs 58.0
- `colorscale`: different color arrays
- `zmin`: 70 vs 40
- `tick_vals/tick_text`: different labels
- `title`: "ANI" vs "AAI"

Create a dataclass and shared function:

```python
@dataclass
class HeatmapConfig:
    metric_name: str
    default_value: float
    same_threshold: float
    boundary_lower: float
    colorscale: list
    zmin: float
    zmax: float = 100.0
    tick_vals: list[float] = field(default_factory=list)
    tick_text: list[str] = field(default_factory=list)

ANI_HEATMAP_CONFIG = HeatmapConfig(
    metric_name="ANI",
    default_value=70.0,
    same_threshold=95.0,
    boundary_lower=93.0,
    colorscale=[...],  # existing ANI colorscale
    zmin=70,
    tick_vals=[70, 75, 80, 85, 90, 95, 100],
    tick_text=["70", "75", "80", "85 (genus)", "90", "95 (species)", "100"],
)

AAI_HEATMAP_CONFIG = HeatmapConfig(
    metric_name="AAI",
    default_value=40.0,
    same_threshold=65.0,
    boundary_lower=58.0,
    colorscale=[...],  # existing AAI colorscale
    zmin=40,
    tick_vals=[40, 50, 60, 65, 70, 80, 90, 100],
    tick_text=["40", "50", "60", "65 (genus)", "70", "80", "90", "100"],
)

def _build_similarity_heatmap(
    matrix: pl.DataFrame,
    genome_labels_map: dict[str, str],
    config: HeatmapConfig,
) -> tuple[go.Figure, SimilarityStats, bool]:
    # Single implementation using config values
    ...
```

Make `build_ani_heatmap` and `build_aai_heatmap` thin wrappers:
```python
def build_ani_heatmap(ani_matrix, genome_labels_map, default_ani=70.0):
    return _build_similarity_heatmap(ani_matrix, genome_labels_map, ANI_HEATMAP_CONFIG)

def build_aai_heatmap(aai_matrix, genome_labels_map, default_aai=40.0):
    return _build_similarity_heatmap(aai_matrix, genome_labels_map, AAI_HEATMAP_CONFIG)
```

**Step 2: Do the same for stats cards**

`build_ani_stats_cards()` (line 282) and `build_aai_stats_cards()` (line 323) generate nearly identical HTML. Create a shared `_build_stats_cards(stats, metric_name, boundary_desc)` function.

**Step 3: Run tests**

Run: `pytest tests/unit/test_report_generator.py -v --tb=short`
Expected: All pass

**Step 4: Commit**

```bash
git add -A
git commit -m "refactor: unify ANI/AAI heatmap builders into parameterized function"
```

---

### Task 11: Consolidate _classify_partition into classify_file

**Files:**
- Modify: `src/metadarkmatter/core/classification/classifiers/vectorized.py`

**Step 1: Refactor stream_to_file to use classify_file**

`_classify_partition()` (line 1010, ~260 lines) duplicates `classify_file()`. Refactor `stream_to_file()` to call `classify_file()` per partition instead.

In `stream_to_file()`, replace:
```python
result = self._classify_partition(partition_df)
```
with:
```python
result = self.classify_file(partition_df)
```

Note: `classify_file()` currently accepts a `Path` argument. If `_classify_partition()` accepts a DataFrame, you may need to make `classify_file()` accept either a Path or DataFrame. Check the signature and adjust accordingly.

**Step 2: Delete _classify_partition method**

Remove the entire `_classify_partition()` method (~260 lines starting at line 1010).

**Step 3: Run tests**

Run: `pytest tests/unit/test_family_validation.py tests/unit/test_classification.py -v --tb=short`
Expected: All pass

**Step 4: Commit**

```bash
git add -A
git commit -m "refactor: remove _classify_partition duplication, reuse classify_file"
```

---

### Task 12: Unify _effective_thresholds

**Files:**
- Modify: `src/metadarkmatter/core/classification/thresholds.py`

**Step 1: Replace _effective_thresholds with config.get_effective_thresholds**

In `thresholds.py`, the `_effective_thresholds()` function (line 163-203) duplicates `ScoringConfig.get_effective_thresholds()`. Replace:

```python
# Old:
eff = _effective_thresholds(config)

# New:
eff = config.get_effective_thresholds()
```

Then delete the `_effective_thresholds()` function entirely.

Note: Check that the keys used by `apply_classification_thresholds()` match the dict keys returned by `config.get_effective_thresholds()`. Both should have `novelty_known_max`, `novelty_novel_species_min`, etc.

**Step 2: Run tests**

Run: `pytest tests/unit/test_sensitivity.py tests/unit/test_classification.py -v --tb=short`
Expected: All pass

**Step 3: Commit**

```bash
git add -A
git commit -m "refactor: unify threshold resolution to single implementation"
```

---

### Task 13: Consolidate clustering and fix ThresholdConfig

**Files:**
- Modify: `src/metadarkmatter/visualization/report/components/heatmap_builder.py`
- Modify: `src/metadarkmatter/visualization/plots/base.py`

**Step 1: Consolidate clustering in heatmap_builder.py**

Find `_cluster_extended_matrix()` (line ~558 in heatmap_builder.py). It duplicates `perform_hierarchical_clustering()` from `clustering.py` but also reorders an `is_novel_mask`. Refactor to call the canonical function:

```python
def _cluster_extended_matrix(z_filled, genome_accessions, default_value, is_novel_mask):
    z_clustered, ordered_accessions, success = perform_hierarchical_clustering(
        z_filled, genome_accessions, default_value
    )
    if success:
        # Reorder mask to match clustered order
        acc_to_idx = {acc: i for i, acc in enumerate(genome_accessions)}
        ordered_mask = [is_novel_mask[acc_to_idx[acc]] for acc in ordered_accessions]
        return z_clustered, ordered_accessions, ordered_mask, success
    return z_filled, genome_accessions, is_novel_mask, success
```

**Step 2: Fix ThresholdConfig defaults in plots/base.py**

`ThresholdConfig` (line 123) has `novelty_known_max=2.0` but `ScoringConfig` uses `4.0`. Fix:

```python
@dataclass
class ThresholdConfig:
    """Classification threshold configuration."""
    novelty_known_max: float = 4.0          # was 2.0 - bug fix
    novelty_novel_species_min: float = 4.0  # was 5.0 - align with ScoringConfig
    novelty_novel_species_max: float = 20.0
    novelty_novel_genus_min: float = 20.0
    novelty_novel_genus_max: float = 25.0
    uncertainty_known_max: float = 1.5      # was 0.5 - align with ScoringConfig
    uncertainty_novel_species_max: float = 1.5  # was 0.5
    uncertainty_novel_genus_max: float = 1.5    # was 2.0
    uncertainty_conserved_min: float = 5.0
```

Check `ScoringConfig` defaults (config.py lines ~94-144) for the correct values and match them exactly.

**Step 3: Run tests**

Run: `pytest tests/unit/test_visualization_plots.py tests/unit/test_report_generator.py -v --tb=short`
Expected: All pass (visualization tests may need threshold value updates)

**Step 4: Commit**

```bash
git add -A
git commit -m "refactor: consolidate clustering logic and fix ThresholdConfig defaults"
```

---

### Task 14: Phase 3 verification

**Step 1: Run full test suite**

Run: `pytest tests/ --ignore=tests/unit/test_phylogeny_placement.py --ignore=tests/unit/test_phylogeny_tree_builder.py -v --tb=short 2>&1 | tail -30`
Expected: All pass

**Step 2: Commit if any fixes needed**

---

### Task 15: Rename cli/util.py to cli/mapping.py

**Files:**
- Rename: `src/metadarkmatter/cli/util.py` -> `src/metadarkmatter/cli/mapping.py`
- Modify: `src/metadarkmatter/cli/main.py` (update import of util subcommand)

**Step 1: Rename the file**

```bash
git mv src/metadarkmatter/cli/util.py src/metadarkmatter/cli/mapping.py
```

**Step 2: Update main.py import**

Find where `util` is imported/registered as a Typer subcommand and update to `mapping`:

```python
# Old:
from metadarkmatter.cli.util import app as util_app
# New:
from metadarkmatter.cli.mapping import app as util_app
```

Note: Keep the CLI command name as `util` for backward compatibility - only the file name changes.

**Step 3: Run tests**

Run: `pytest tests/ --ignore=tests/unit/test_phylogeny_placement.py --ignore=tests/unit/test_phylogeny_tree_builder.py -x -q`
Expected: All pass

**Step 4: Commit**

```bash
git add -A
git commit -m "refactor: rename cli/util.py to cli/mapping.py for clarity"
```

---

### Task 16: Final cleanup and verification

**Files:**
- Modify: `src/metadarkmatter/__init__.py`
- Modify: `src/metadarkmatter/core/__init__.py`
- Modify: `CLAUDE.md`

**Step 1: Update root __init__.py re-exports**

Verify all re-exports reference existing classes only. Remove any references to `ParallelClassifier`, `SparseANIMatrix`, or `GlobalConfig`.

**Step 2: Update core/__init__.py**

Same cleanup.

**Step 3: Update CLAUDE.md**

Remove references to `ParallelClassifier`, `--fast` flag, `--parallel` flag. Update the Key Architecture section to reflect simplified classifier structure. Note that all classification now uses VectorizedClassifier.

**Step 4: Run full test suite**

Run: `pytest tests/ --ignore=tests/unit/test_phylogeny_placement.py --ignore=tests/unit/test_phylogeny_tree_builder.py -v --tb=short`
Expected: All pass

**Step 5: Check CLI help**

```bash
metadarkmatter score classify --help
```
Expected: No `--fast` or `--parallel` flags. `--streaming` still present.

**Step 6: Count total impact**

```bash
git diff --stat main
```

**Step 7: Commit**

```bash
git add -A
git commit -m "chore: final cleanup and documentation updates for codebase simplification"
```
