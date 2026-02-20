# Codebase Simplification Design

**Date:** 2026-02-20
**Status:** Design

## Problem

The metadarkmatter codebase has accumulated substantial duplication and dead code through iterative development. The classification algorithm is implemented independently in 4-5 locations, scoring functions exist in 3 copies each, and approximately 3,000 lines of dead code remain in the repository. This increases maintenance burden and risk of logic drift between implementations.

## Design

### Overview

Simplify the codebase in four phases: dead code removal, classifier consolidation, duplication reduction, and structural cleanup. The primary architectural decision is to make `VectorizedClassifier` the single classifier implementation, removing `ParallelClassifier` and the iterator-based paths in `ANIWeightedClassifier`.

### Phase 1: Dead Code Removal

Remove code that is defined but never used from any production path.

| Item | File | Lines Removed |
|------|------|---------------|
| `ani_placement.py.backup` | `core/ani_placement.py.backup` | ~2,762 |
| `ParallelClassifier` | `classifiers/parallel.py` | ~255 |
| `SparseANIMatrix` | `classification/sparse_ani_matrix.py` | ~116 |
| `GlobalConfig` + unused sub-configs | `models/config.py` | ~100 |
| 11 unused constants | `core/constants.py` | ~15 |
| 5 unused imports | `classifiers/vectorized.py` | 5 |
| Backward-compat aliases | `components/extended_matrix_builder.py` | ~50 |
| `calculate_protein_confidence_score` | `core/protein_constants.py` | ~20 |
| `BlastResult.sorted_by_bitscore` | `models/blast.py` | ~10 |

Also: add missing `Callable` import to classifier files.

### Phase 2: Classifier Consolidation

Route all CLI paths through `VectorizedClassifier` and remove redundant code from `base.py`.

**Current CLI flow:**
- `--parallel` / `--streaming` -> `VectorizedClassifier.classify_file()` / `stream_to_file()`
- `--fast` -> `ANIWeightedClassifier.classify_to_dataframe_fast()`
- default -> `ANIWeightedClassifier.classify_to_dataframe()`

**New CLI flow:**
- All modes -> `VectorizedClassifier.classify_file()` or `stream_to_file()`

**Code to remove from `base.py`:**
- `_classify_chunk_worker()` (~320 lines)
- Five `_*_inline()` scoring functions (~160 lines)
- `classify_read_fast()` and `classify_to_dataframe_fast()` (~200 lines)
- `classify_blast_file()` and `classify_to_dataframe()` (~80 lines)

**Code to keep in `base.py`:**
- `ANIWeightedClassifier` class with `classify_read()` as a public programmatic API
- Core methods like `_classify_by_thresholds()` for the API path
- Constructor, config handling

**CLI changes (`cli/score.py`):**
- Remove `--fast` flag (VectorizedClassifier handles all cases)
- Remove `--parallel` flag (Polars is inherently parallel)
- Route all classification through VectorizedClassifier
- Keep `--streaming` flag (controls `stream_to_file` vs `classify_file`)

**Update `ani_placement.py` re-export shim** to import directly from classification package. Update CLI imports accordingly.

### Phase 3: Duplication Reduction

**3a. Unify ANI/AAI heatmap builders** (`heatmap_builder.py`)

Replace `build_ani_heatmap()` and `build_aai_heatmap()` with a single parameterized `_build_similarity_heatmap()`. The ANI/AAI functions become thin wrappers passing config. Same for `build_ani_stats_cards()` / `build_aai_stats_cards()`.

**3b. Consolidate `_classify_partition`** (`vectorized.py`)

Refactor `stream_to_file()` to call `classify_file()` per partition instead of using the duplicated `_classify_partition()` method. Remove `_classify_partition()`.

**3c. Unify `_effective_thresholds`** (`thresholds.py`)

Remove `_effective_thresholds()` from `thresholds.py`. Have `apply_classification_thresholds()` call `config.get_effective_thresholds()` instead.

**3d. Consolidate clustering** (`heatmap_builder.py`)

Have `_cluster_extended_matrix()` call `perform_hierarchical_clustering()` from `clustering.py` and handle mask reordering separately.

**3e. Fix ThresholdConfig defaults** (`plots/base.py`)

Align `ThresholdConfig` default values with `ScoringConfig` defaults, or derive from `ScoringConfig`.

### Phase 4: Structural Cleanup

**4a. Rename `cli/util.py`** to `cli/mapping.py` to disambiguate from `cli/utils.py`.

**4b. Update `__init__.py` re-exports** to remove references to deleted classes.

**4c. Update tests** to use VectorizedClassifier directly where they previously used ParallelClassifier or the slow base classifier paths.

### Backward Compatibility

- The `ANIWeightedClassifier.classify_read()` method remains available as a programmatic API
- CLI output format is unchanged
- `--parallel` and `--fast` flags are removed (VectorizedClassifier handles everything)
- All external tool wrappers remain unchanged

### Testing Strategy

- Run full test suite after each phase
- Update tests that reference removed classes/functions
- Verify classification output is identical before/after for a reference dataset
- Remove tests for deleted code, update tests for modified code

### Estimated Impact

| Metric | Before | After |
|--------|--------|-------|
| Lines in classifiers/ | ~3,300 | ~1,600 |
| Classification implementations | 4-5 | 1 (Polars) + 1 (Python API) |
| Scoring function copies | 3 each | 1 each |
| Dead code files | 1 backup | 0 |
