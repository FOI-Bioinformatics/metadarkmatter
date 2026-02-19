# Family Validation for Broad-Database Classification

**Date:** 2026-02-19
**Status:** Design

## Problem

When users run BLAST or MMseqs2 against a broad database (e.g., all bacteria rather than a single family), reads produce hits spanning multiple taxonomic levels. The current classifier treats all hits as peers within a single family, which leads to:

1. **Kraken2 misclassifications appearing as novel diversity**: Reads that don't belong to the target family are forced to match family genomes and appear as "highly novel."
2. **No way to validate family assignment**: Users cannot distinguish genuine novel diversity from taxonomic misassignment without a separate validation step.
3. **Inflated ambiguity counts**: Out-of-family hits compete with in-family hits, increasing placement uncertainty for reads that are actually well-placed within the family.

## Design

### Overview

Integrate taxonomic context from broad-database alignments into the classification pipeline. The core mechanism: partition BLAST hits by ANI matrix membership (in-family vs external), compute family affinity metrics, and add an "Off-target" classification category for reads better explained by an external genome.

### Hit Partitioning

Every hit is classified as in-family or external based on whether the subject genome appears in the ANI matrix. The ANI matrix defines the family scope.

```
Read -> All hits
          |-- In-family hits (genome in ANI matrix) -> existing classifier
          |-- External hits (genome not in ANI matrix) -> family context metrics
```

This partition happens after alignment quality filtering and ID mapping, but before bitscore threshold and ambiguous hit identification.

### Family Determination

The target family is determined by (in priority order):

1. `--target-family "f__Francisellaceae"` explicit CLI flag
2. Inferred from `--metadata` genome_metadata.tsv (most common family)
3. If neither available, family validation is disabled; classification runs as today

### Family Context Metrics

For each read, we compute:

| Metric | Formula | Interpretation |
|--------|---------|----------------|
| `family_bitscore_ratio` | `best_in_family_bitscore / best_overall_bitscore` | 1.0 = best hit is in-family. <0.8 = likely off-target. |
| `external_best_identity` | Top pident from external hits | Shows how well the read matches outside the family |
| `family_identity_gap` | `best_in_family_pident - best_external_pident` | Positive = family hit wins. Negative = external hit is better match. |
| `in_family_hit_fraction` | `count(in_family_hits) / count(all_hits)` | Fraction of competitive hits from within-family genomes |
| `external_best_genome` | string | Best-matching genome outside the family |
| `external_best_family` | string | Family of best external hit (from metadata if available) |

### Off-target Decision Rule

A read is classified as "Off-target" when:

1. `family_bitscore_ratio < threshold` (default 0.8, configurable via `--family-ratio-threshold`), OR
2. No in-family hits exist at all

The default threshold of 0.8 means: if the best external hit has a bitscore more than 25% better than the best in-family hit, the read probably does not belong to this family.

### Classification Integration

Modified flow in `VectorizedClassifier.classify_file()`:

```
1. Load BLAST data (existing)
2. Apply ID mapping if provided (existing)
3. Apply alignment quality filters (existing)
4. Extract genome names (existing)
5. NEW: Partition hits by ANI matrix membership
6. NEW: Compute family context metrics per read
7. NEW: Apply off-target rule - reads failing get "Off-target" category
8. For remaining reads: run existing classification on IN-FAMILY HITS ONLY
9. NEW: Adjust confidence score using family context
10. Combine off-target + classified results
11. Output with new columns
```

Key change: when a read has both in-family and external hits, only the in-family hits participate in best hit selection, ambiguous hit identification, and ANI-based uncertainty calculation. External hits are excluded from ANI lookup since they are not in the matrix.

### Confidence Score Adjustment

For in-family reads, a family context component (0-10 points) replaces half of the current phylogenetic context allocation:

- ratio >= 0.95 and positive identity gap: +10
- ratio >= 0.90: +7
- ratio >= 0.85: +3
- ratio < 0.85 (borderline): +0

The remaining 10 points come from within-family ambiguity scope (as today).

### New Output Columns

When family validation is active:

| Column | Type | Description |
|--------|------|-------------|
| `family_bitscore_ratio` | float | Best in-family / best overall bitscore |
| `family_identity_gap` | float | Best in-family pident - best external pident |
| `in_family_hit_fraction` | float | Fraction of competitive hits from in-family genomes |
| `external_best_genome` | string | Best-matching genome outside the family |
| `external_best_identity` | float | Percent identity of best external hit |
| `external_best_family` | string | Family of best external hit (from metadata) |

New taxonomic_call value: `"Off-target"` joins the existing set.

### CLI Interface

New flags on `score classify`:

```bash
metadarkmatter score classify \
    --alignment results.tsv.gz \
    --ani ani_matrix.csv \
    --target-family "f__Francisellaceae" \
    --family-ratio-threshold 0.8 \
    --output classifications.csv \
    --parallel
```

- `--target-family`: Optional. Enables family validation. Without this flag (and without metadata), classification runs identically to current behavior.
- `--family-ratio-threshold`: Default 0.8. Off-target threshold for family bitscore ratio.

### Report Changes

When family validation is active, add a "Family Validation" tab to the HTML report:

- Pie chart: in-family vs off-target read proportions
- Histogram: `family_bitscore_ratio` distribution across all reads
- Table: top external families receiving off-target reads
- Summary statistic: percentage of reads validated as belonging to target family

### Summary JSON Changes

Add to existing summary output:

```json
{
  "off_target_count": 1234,
  "off_target_pct": 12.3,
  "family_validated_pct": 87.7,
  "target_family": "f__Francisellaceae",
  "top_external_families": {
    "f__Piscirickettsiaceae": 456,
    "f__Enterobacteriaceae": 234
  }
}
```

### Backward Compatibility

- All new columns only appear when family validation is active
- Without `--target-family` or `--metadata`, output is identical to current behavior
- No breaking changes to existing CLI flags, output formats, or report sections
- The Off-target category is added to the TAXONOMIC_TO_DIVERSITY mapping

### Key Files to Modify

| File | Changes |
|------|---------|
| `core/classification/classifiers/vectorized.py` | Hit partitioning, family metrics, off-target rule, in-family classification |
| `core/classification/thresholds.py` | Add Off-target to threshold application logic |
| `models/config.py` | Add `target_family`, `family_ratio_threshold` to ScoringConfig |
| `models/classification.py` | Add "Off-target" to TaxonomicCall and TAXONOMIC_TO_DIVERSITY |
| `cli/score.py` | Add `--target-family` and `--family-ratio-threshold` flags |
| `core/metadata.py` | Add `infer_target_family()` method |
| `visualization/report/generator.py` | Add Family Validation report tab |
| `visualization/report/templates.py` | Add Family Validation HTML/JS templates |
| `core/constants.py` | Add CATEGORY_OFF_TARGET constant |

### Testing Strategy

1. **Unit tests**: Family context metric calculations with known input
2. **Unit tests**: Off-target decision rule with various ratio thresholds
3. **Unit tests**: In-family-only classification produces same results as current classifier when no external hits exist
4. **Integration test**: Full pipeline with broad-database alignment file
5. **Backward compatibility test**: Verify identical output when family validation is not enabled
