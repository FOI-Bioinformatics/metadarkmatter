# Alignment Length/Coverage Integration Plan

## Executive Summary

This plan addresses the incorporation of alignment length and coverage into the metadarkmatter classification system to prioritize hits that span larger portions of reads, resulting in cleaner classification results.

**Current Problem**: Classification uses pure bitscore ranking, which can select short high-identity hits (e.g., conserved domains) over longer moderate-identity hits that better represent the full read sequence.

**Recommended Solution**: Implement coverage-weighted hit selection (Approach 2) that applies a composite scoring function combining bitscore and alignment coverage, with configurable weighting modes and backward compatibility.

**Implementation Timeline**: 5-7 days across 4 phases

---

## Problem Statement

### Current Behavior

The classification algorithm selects the "best hit" purely by bitscore:

```python
# From models/blast.py:200-202
@property
def best_hit(self) -> BlastHit | None:
    """Get the highest-scoring BLAST hit."""
    return self.hits[0] if self.hits else None  # Sorted by bitscore descending
```

**Issue**: A read with two hits:
- Hit A: 98% identity, 50bp alignment (conserved domain), bitscore 95
- Hit B: 92% identity, 280bp alignment (full read), bitscore 90

Current system selects Hit A (higher bitscore), classifying as "Known Species" based on 98% identity. However, Hit B provides better evidence for the read's origin because it covers 93% of the read length versus 17% for Hit A.

### Biological Rationale

Short high-identity alignments often represent:
- Conserved domains (16S rRNA, housekeeping genes)
- Repetitive sequences
- Horizontal gene transfer elements

Longer alignments spanning most of the read provide:
- More phylogenetic information
- Better species-level placement
- Reduced false positive classifications

This aligns with GTDB's requirement for ≥50% alignment fraction (AF) for genome-based species delineation.

---

## Recommended Approach: Coverage-Weighted Hit Selection

### Overview

Implement a composite scoring system that weights bitscore by alignment coverage, with three configurable weighting modes:

**Composite Score Formula**:
```
weighted_score = bitscore × coverage_weight_function(coverage)
```

Where coverage = `(qend - qstart + 1) / read_length_proxy`

### Three Weighting Modes

#### 1. Linear Mode (Default)
```python
weight = min_weight + (max_weight - min_weight) × coverage
# Example at strength=0.5: weight = 0.5 + 0.5 × coverage
```

- 50% coverage → weight = 0.75
- 90% coverage → weight = 0.95
- Gradual penalty for low coverage

#### 2. Logarithmic Mode (Gentle)
```python
weight = min_weight + (max_weight - min_weight) × log(1 + coverage × 9) / log(10)
```

- Reduces penalty for moderate coverage (50-70%)
- Still penalizes very short hits (<30%)

#### 3. Sigmoid Mode (Threshold-like)
```python
weight = min_weight + (max_weight - min_weight) / (1 + exp(-10 × (coverage - 0.6)))
```

- Sharp transition around 60% coverage
- Mimics binary filter but smoother
- Good for enforcing stricter standards

### Configuration Parameters

Add to `ScoringConfig` (models/config.py):

```python
coverage_weight_mode: Literal["none", "linear", "log", "sigmoid"] = Field(
    default="none",
    description="Mode for weighting hits by coverage (none = disabled, backward compatible)",
)

coverage_weight_strength: float = Field(
    default=0.5,
    ge=0.0,
    le=1.0,
    description=(
        "Strength of coverage weighting (0.0 = no effect, 1.0 = full effect). "
        "At 0.5, a 50% coverage hit gets 0.75× weight in linear mode."
    ),
)
```

### Why This Approach?

**Pros**:
- Graduated weighting (not binary cutoff)
- Tunable via presets
- Biological intuition: longer alignments = more confidence
- Can be tested incrementally
- Backward compatible (default mode="none")

**Cons**:
- Adds complexity to hit selection
- Requires validation across three classifier pathways
- Parameter tuning needed for optimal defaults

---

## Alternative Approaches (For Reference)

### Approach 1: Enhanced Filtering (Conservative)

Extend existing `min_alignment_fraction` filter with stricter defaults in new presets. Trivial to implement but binary cutoff loses borderline reads.

### Approach 3: Adaptive Coverage-Quality Score (Advanced)

Unified metric `AQS = bitscore^0.7 × pident^0.2 × coverage^0.5` for both hit selection AND confidence scoring. More comprehensive but higher complexity; consider for future version.

---

## Critical Files to Modify

### 1. `/src/metadarkmatter/models/blast.py`
**Purpose**: Add coverage calculation and weighted scoring methods to BlastHit model

**Changes**:
- Add `calculate_coverage(read_length)` method (after line 95)
- Add `calculate_weighted_score(mode, strength, read_length)` method
- Update `BlastResult.get_best_hit()` to support coverage weighting
- Update `BlastResult.iter_ambiguous_hits()` for weighted threshold

### 2. `/src/metadarkmatter/models/config.py`
**Purpose**: Add configuration parameters

**Changes** (after line 209):
- Add `coverage_weight_mode: Literal["none", "linear", "log", "sigmoid"]`
- Add `coverage_weight_strength: float` (0.0 to 1.0 range)

### 3. `/src/metadarkmatter/core/classification/classifiers/base.py`
**Purpose**: Update canonical classifier to use coverage weighting

**Changes** (lines 74-89):
- Call `blast_result.get_best_hit(mode, strength)` instead of `.best_hit`
- Pass coverage parameters to `iter_ambiguous_hits()`

### 4. `/src/metadarkmatter/core/classification/classifiers/vectorized.py`
**Purpose**: Add Polars expressions for coverage weighting

**Changes** (lines 190-235):
- Add coverage calculation column
- Add weight calculation column based on mode
- Add `weighted_score` column
- Sort by `weighted_score` instead of `bitscore` when enabled

### 5. `/src/metadarkmatter/cli/score.py`
**Purpose**: Add CLI parameters and presets

**Changes**:
- Add presets: `coverage-linear`, `coverage-strict`, `coverage-gentle`, `gtdb-coverage`
- Add CLI parameters: `--coverage-weight-mode`, `--coverage-weight-strength`

---

## Implementation Phases

### Phase 1: Core Models (Days 1-2)
- Add coverage calculation to `BlastHit`
- Add weighted score calculation with three modes (linear, log, sigmoid)
- Update `BlastResult` hit selection methods
- Write comprehensive unit tests

**Testing**: Verify that longer hits are preferred when weighting enabled, mode="none" preserves original behavior

### Phase 2: Configuration (Day 2)
- Add config parameters to `ScoringConfig`
- Validate parameter ranges with Pydantic
- Write config validation tests

### Phase 3: Classifiers (Days 3-5)
- Update base classifier to use weighted hit selection
- Update vectorized classifier with Polars expressions for batch weighting
- Verify parallel classifier compatibility
- Write integration tests comparing all three classifiers

**Critical Test**: Ensure vectorized and parallel classifiers produce identical results to base classifier with coverage weighting enabled

### Phase 4: CLI and Presets (Days 6-7)
- Add CLI parameters with typer.Option
- Create four presets with different coverage strategies
- Update documentation (REFERENCE.md, ALGORITHM_DETAILED.md, USER_GUIDE.md)
- Write end-to-end CLI tests

---

## Testing Strategy

### Unit Tests (`tests/models/test_blast_coverage.py`)
- Coverage calculation with various read lengths
- Weighted score calculation for all three modes
- Strength parameter validation (0.0 to 1.0)
- Best hit selection changes with weighting
- Ambiguous hit detection with weighting

### Integration Tests (`tests/integration/test_coverage_pipeline.py`)
- End-to-end pipeline with coverage weighting
- Vectorized vs parallel classifier consistency
- Backward compatibility (mode="none" produces identical results)

### Validation Scenarios

**Scenario 1**: Short conserved domain vs long divergent hit
- Read: 300bp
- Hit A: 98% identity, 50bp, bitscore 95 (conserved domain)
- Hit B: 88% identity, 280bp, bitscore 88 (nearly full read)
- **Expected**: Mode "none" selects A; mode "linear" selects B

**Scenario 2**: Multiple hits at species boundary
- Verify that coverage weighting affects placement uncertainty calculation correctly

---

## Backward Compatibility

All new parameters default to disabled:
```python
coverage_weight_mode: str = "none"  # No weighting by default
coverage_weight_strength: float = 0.5  # Irrelevant when mode="none"
```

**Guarantees**:
1. Default configuration unchanged (existing pipelines produce identical results)
2. Existing presets unchanged (`default`, `gtdb-strict` maintain original behavior)
3. Python API compatible (existing calls work without modification)
4. CLI compatible (existing commands produce same output)

---

## Verification Steps

### After Implementation

1. **Run unit tests**: All new tests pass with >95% coverage
   ```bash
   pytest tests/models/test_blast_coverage.py -v
   ```

2. **Run integration tests**: Verify classifier consistency
   ```bash
   pytest tests/integration/test_coverage_pipeline.py -v
   ```

3. **Backward compatibility test**: Confirm mode="none" produces identical results
   ```bash
   # Run same classification with and without explicit mode="none"
   diff <(metadarkmatter score classify --alignment test.tsv --ani ani.csv) \
        <(metadarkmatter score classify --alignment test.tsv --ani ani.csv --coverage-weight-mode none)
   ```

4. **Benchmark performance**: Measure overhead (<5% expected)
   ```bash
   time metadarkmatter score classify --alignment large.tsv --ani ani.csv  # Baseline
   time metadarkmatter score classify --alignment large.tsv --ani ani.csv --preset coverage-linear  # Weighted
   ```

5. **Visual validation**: Generate reports with both modes and compare classification distributions
   ```bash
   metadarkmatter report generate --classifications default.csv --output report_default.html
   metadarkmatter report generate --classifications weighted.csv --output report_weighted.html
   # Expect: 5-15% reduction in "Known Species", increase in "Conserved Region"
   ```

---

## Success Metrics

### Quantitative
- Classification rate change: 5-15% reduction in "Known Species" as conserved domains deprioritized
- Confidence score improvement: +5-10 points average for novel taxa
- Performance overhead: <5% runtime increase
- Test coverage: >95% for new code

### Qualitative
- Longer alignments provide better phylogenetic signal (manual review)
- Users can choose appropriate preset without consulting developers
- Documentation includes clear worked examples

---

## Documentation Updates

### Files to Update
1. **docs/REFERENCE.md**: Add "Alignment Coverage Weighting" section with formulas and examples
2. **docs/ALGORITHM_DETAILED.md**: Update "Step 2: Identify Top Hit" with coverage weighting logic
3. **docs/CLASSIFICATION_STATISTICS.md**: Add coverage weighting to "Additional Quality Metrics"
4. **docs/USER_GUIDE.md**: Add "Coverage Weighting" workflow guide with preset recommendations
5. **README.md**: Add coverage weighting to feature list

### Example Documentation Snippet

```markdown
## Coverage Weighting

By default, metadarkmatter ranks hits by bitscore. Enable coverage weighting to prioritize hits spanning larger portions of reads:

# Linear coverage weighting (balanced)
metadarkmatter score classify \
  --alignment results.tsv \
  --ani ani.csv \
  --preset coverage-linear \
  --output classifications.csv

# Strict coverage requirement (sigmoid mode, ~60% threshold)
metadarkmatter score classify --preset coverage-strict ...

# Combine with GTDB alignment fraction filter
metadarkmatter score classify --preset gtdb-coverage ...
```

---

## Risk Assessment

### High Risk
**Bitscore interpretation change**: Coverage weighting changes which hit is "best", affecting placement uncertainty calculations.

**Mitigation**: Extensive validation with real datasets, comparison to original results, user testing before default change.

### Medium Risk
**Performance overhead**: Calculating weighted scores adds computation.

**Mitigation**: Benchmark with large datasets, optimize Polars expressions, profile hot paths.

### Low Risk
**Backward compatibility**: New parameters might conflict with existing workflows.

**Mitigation**: Default to disabled, maintain existing presets unchanged, comprehensive documentation.

---

## Future Enhancements

1. **Read length tracking**: Track actual read length from FASTQ instead of using `qend` proxy
2. **Coverage-dependent thresholds**: Adjust novelty/uncertainty thresholds based on coverage
3. **Multi-hit aggregation**: Combine evidence from multiple alignments to same genome
4. **Approach 3 implementation**: If Approach 2 proves successful, implement unified AQS metric

---

## Summary

This plan implements coverage-weighted hit selection to address short high-identity hits dominating classification. The solution:

1. Adds coverage calculation to BlastHit model
2. Implements three weighting modes (linear, log, sigmoid) with tunable strength
3. Updates all three classifier pathways (base, vectorized, parallel)
4. Maintains backward compatibility via default mode="none"
5. Provides user-friendly presets

**Estimated timeline**: 5-7 days
**Success criteria**: Improved classification accuracy, <5% performance overhead, positive user feedback
**Next steps**: Await approval, then begin Phase 1 (Core Models)
