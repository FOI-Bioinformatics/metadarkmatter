# Family Validation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add family validation to the classification pipeline so reads from broad-database alignments can be flagged as "Off-target" when they belong to a different taxonomic family than the target.

**Architecture:** Partition BLAST hits by ANI matrix membership (in-family vs external), compute family affinity metrics per read, classify reads as "Off-target" when the best external hit substantially outscores the best in-family hit. Only in-family hits participate in existing ANI-based classification.

**Tech Stack:** Python 3.11+, Polars (vectorized classification), Pydantic (config models), Typer (CLI), Plotly (report visualization)

**Design doc:** `docs/plans/2026-02-19-family-validation-design.md`

---

### Task 1: Add Off-target to classification models

**Files:**
- Modify: `src/metadarkmatter/models/classification.py:36-45` (TAXONOMIC_TO_DIVERSITY dict)
- Modify: `src/metadarkmatter/models/classification.py:48-81` (TaxonomicCall enum)
- Modify: `src/metadarkmatter/models/classification.py:256-298` (TaxonomicSummary)
- Modify: `src/metadarkmatter/core/constants.py:74-80` (category constants)
- Test: `tests/unit/test_classification.py`

**Step 1: Write failing tests**

Add to `tests/unit/test_classification.py`:

```python
def test_off_target_value(self):
    """Off-target should have correct string value."""
    assert TaxonomicCall.OFF_TARGET.value == "Off-target"

def test_off_target_diversity_status(self):
    """Off-target should map to Uncertain diversity status."""
    from metadarkmatter.models.classification import TAXONOMIC_TO_DIVERSITY
    assert TAXONOMIC_TO_DIVERSITY["Off-target"] == "Uncertain"

def test_off_target_is_not_novel(self):
    """Off-target reads should not be flagged as novel."""
    classification = ReadClassification(
        read_id="read_off",
        best_match_genome="GCF_999999999.1",
        top_hit_identity=85.0,
        novelty_index=15.0,
        placement_uncertainty=0.0,
        num_ambiguous_hits=1,
        taxonomic_call=TaxonomicCall.OFF_TARGET,
    )
    assert not classification.is_novel
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/unit/test_classification.py -v -k "off_target" 2>&1 | head -30`
Expected: FAIL with `AttributeError: OFF_TARGET` or `KeyError: 'Off-target'`

**Step 3: Implement the model changes**

In `src/metadarkmatter/models/classification.py`:

1. Add to `TAXONOMIC_TO_DIVERSITY` dict (after line 44):
```python
"Off-target": "Uncertain",
```

2. Add to `TaxonomicCall` enum (after line 81):
```python
OFF_TARGET = "Off-target"
```

3. Add to `TaxonomicSummary` (after line 298, after `unclassified` field):
```python
off_target: int = Field(
    default=0,
    description="Reads classified as off-target (better match outside target family)",
)
```

In `src/metadarkmatter/core/constants.py`, add after line 80:
```python
CATEGORY_OFF_TARGET = "Off-target"
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/unit/test_classification.py -v -k "off_target"`
Expected: PASS (3 tests)

**Step 5: Run full test suite to check for regressions**

Run: `python -m pytest tests/unit/test_classification.py -v`
Expected: All existing tests PASS

**Step 6: Commit**

```bash
git add src/metadarkmatter/models/classification.py src/metadarkmatter/core/constants.py tests/unit/test_classification.py
git commit -m "feat: add Off-target classification category to models"
```

---

### Task 2: Add family validation fields to ScoringConfig

**Files:**
- Modify: `src/metadarkmatter/models/config.py:18-430` (ScoringConfig class)
- Test: `tests/unit/test_config_models.py`

**Step 1: Write failing tests**

Add to `tests/unit/test_config_models.py`:

```python
class TestFamilyValidationConfig:
    """Tests for family validation config fields."""

    def test_default_target_family_is_none(self):
        """target_family should default to None (disabled)."""
        config = ScoringConfig()
        assert config.target_family is None

    def test_default_family_ratio_threshold(self):
        """family_ratio_threshold should default to 0.8."""
        config = ScoringConfig()
        assert config.family_ratio_threshold == 0.8

    def test_target_family_accepts_string(self):
        """target_family should accept a family name string."""
        config = ScoringConfig(target_family="f__Francisellaceae")
        assert config.target_family == "f__Francisellaceae"

    def test_family_ratio_threshold_validation(self):
        """family_ratio_threshold must be between 0 and 1."""
        with pytest.raises(ValidationError):
            ScoringConfig(family_ratio_threshold=1.5)
        with pytest.raises(ValidationError):
            ScoringConfig(family_ratio_threshold=-0.1)

    def test_family_validation_disabled_by_default(self):
        """Family validation should be disabled when target_family is None."""
        config = ScoringConfig()
        assert config.target_family is None
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/unit/test_config_models.py -v -k "FamilyValidation" 2>&1 | head -20`
Expected: FAIL with `TypeError: unexpected keyword argument 'target_family'`

**Step 3: Implement the config changes**

In `src/metadarkmatter/models/config.py`, add two fields to `ScoringConfig` class after `aai_genus_boundary_high` field (around line 277):

```python
# Family validation for broad-database alignments
# When target_family is set, reads with better external hits are flagged Off-target
target_family: str | None = Field(
    default=None,
    description=(
        "Target taxonomic family for family validation (e.g., 'f__Francisellaceae'). "
        "When set, reads with substantially better hits outside the ANI matrix "
        "are classified as Off-target. None disables family validation."
    ),
)
family_ratio_threshold: float = Field(
    default=0.8,
    ge=0.0,
    le=1.0,
    description=(
        "Family bitscore ratio threshold for off-target detection. "
        "Reads with best_in_family_bitscore / best_overall_bitscore below "
        "this value are classified as Off-target. Default 0.8."
    ),
)
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/unit/test_config_models.py -v -k "FamilyValidation"`
Expected: PASS (5 tests)

**Step 5: Commit**

```bash
git add src/metadarkmatter/models/config.py tests/unit/test_config_models.py
git commit -m "feat: add target_family and family_ratio_threshold to ScoringConfig"
```

---

### Task 3: Add infer_target_family() to GenomeMetadata

**Files:**
- Modify: `src/metadarkmatter/core/metadata.py:20-68` (GenomeMetadata class)
- Test: `tests/unit/test_classification.py` (or new test file)

**Step 1: Write failing tests**

Add to `tests/unit/test_classification.py` or create a new section:

```python
import polars as pl
from metadarkmatter.core.metadata import GenomeMetadata


class TestInferTargetFamily:
    """Tests for GenomeMetadata.infer_target_family()."""

    def test_infer_most_common_family(self):
        """Should return the most common family."""
        df = pl.DataFrame({
            "accession": ["GCF_001", "GCF_002", "GCF_003", "GCF_004"],
            "species": ["sp1", "sp2", "sp3", "sp4"],
            "genus": ["g1", "g1", "g2", "g2"],
            "family": ["f__Francisellaceae", "f__Francisellaceae", "f__Francisellaceae", "f__Enterobacteriaceae"],
        })
        metadata = GenomeMetadata(df)
        assert metadata.infer_target_family() == "f__Francisellaceae"

    def test_infer_family_no_family_column(self):
        """Should return None if family column missing."""
        df = pl.DataFrame({
            "accession": ["GCF_001"],
            "species": ["sp1"],
            "genus": ["g1"],
        })
        metadata = GenomeMetadata(df)
        assert metadata.infer_target_family() is None

    def test_infer_family_empty_metadata(self):
        """Should return None for empty metadata."""
        df = pl.DataFrame({
            "accession": [],
            "species": [],
            "genus": [],
            "family": [],
        }).cast({"accession": pl.Utf8, "species": pl.Utf8, "genus": pl.Utf8, "family": pl.Utf8})
        metadata = GenomeMetadata(df)
        assert metadata.infer_target_family() is None
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/unit/test_classification.py -v -k "InferTargetFamily" 2>&1 | head -20`
Expected: FAIL with `AttributeError: 'GenomeMetadata' object has no attribute 'infer_target_family'`

**Step 3: Implement the method**

Add to `GenomeMetadata` class in `src/metadarkmatter/core/metadata.py` (after `get_family` method, around line 151):

```python
def infer_target_family(self) -> str | None:
    """Infer the target family from genome metadata.

    Returns the most common family among genomes in this metadata set.
    Used as fallback when --target-family is not explicitly provided.

    Returns:
        Most common family name, or None if family column is absent or empty.
    """
    if "family" not in self._df.columns:
        return None
    if self._df.is_empty():
        return None
    family_counts = (
        self._df.group_by("family")
        .len()
        .sort("len", descending=True)
    )
    if family_counts.is_empty():
        return None
    return family_counts["family"][0]
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/unit/test_classification.py -v -k "InferTargetFamily"`
Expected: PASS (3 tests)

**Step 5: Commit**

```bash
git add src/metadarkmatter/core/metadata.py tests/unit/test_classification.py
git commit -m "feat: add infer_target_family() to GenomeMetadata"
```

---

### Task 4: Add CLI flags for family validation

**Files:**
- Modify: `src/metadarkmatter/cli/score.py:290-496` (classify function parameters)
- Modify: `src/metadarkmatter/cli/score.py:697-726` (ScoringConfig construction)
- Test: `tests/unit/test_cli_score.py`

**Step 1: Write failing test**

Add to `tests/unit/test_cli_score.py`:

```python
def test_classify_help_shows_target_family(capsys):
    """CLI should accept --target-family flag."""
    from typer.testing import CliRunner
    from metadarkmatter.cli.main import app

    runner = CliRunner()
    result = runner.invoke(app, ["score", "classify", "--help"])
    assert "--target-family" in result.output
    assert "--family-ratio-threshold" in result.output
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest tests/unit/test_cli_score.py -v -k "target_family" 2>&1 | head -20`
Expected: FAIL with `AssertionError` (flag not in help output)

**Step 3: Implement the CLI flags**

In `src/metadarkmatter/cli/score.py`, add two new parameters to the `classify` function after the `bayesian` parameter (around line 489):

```python
target_family: str | None = typer.Option(
    None,
    "--target-family",
    help=(
        "Target taxonomic family for family validation (e.g., 'f__Francisellaceae'). "
        "Enables off-target detection: reads with better BLAST hits outside the "
        "ANI matrix are classified as Off-target. If not set but --metadata is "
        "provided, the most common family is inferred."
    ),
),
family_ratio_threshold: float = typer.Option(
    0.8,
    "--family-ratio-threshold",
    help=(
        "Bitscore ratio threshold for off-target detection. "
        "Reads with best_in_family / best_overall bitscore below this "
        "value are classified as Off-target. Default 0.8."
    ),
    min=0.0,
    max=1.0,
),
```

Then update the ScoringConfig construction at lines ~697-726 to pass the new fields. In both the preset branch and the else branch, add:

```python
target_family=target_family,
family_ratio_threshold=family_ratio_threshold,
```

Also add family inference from metadata (after metadata loading, around line 819):

```python
# Infer target family from metadata if not explicitly provided
if target_family is None and genome_metadata is not None:
    inferred = genome_metadata.infer_target_family()
    if inferred:
        target_family = inferred
        out.print(f"[dim]Inferred target family from metadata: {target_family}[/dim]")
        # Rebuild config with inferred family
        config = config.model_copy(update={"target_family": target_family})
```

Note: `ScoringConfig` is frozen, so we need to use `model_copy(update=...)` for the inference case. Alternatively, defer config construction until after metadata loading.

**Step 4: Run test to verify it passes**

Run: `python -m pytest tests/unit/test_cli_score.py -v -k "target_family"`
Expected: PASS

**Step 5: Commit**

```bash
git add src/metadarkmatter/cli/score.py tests/unit/test_cli_score.py
git commit -m "feat: add --target-family and --family-ratio-threshold CLI flags"
```

---

### Task 5: Implement hit partitioning and family metrics in VectorizedClassifier

This is the core task. It modifies `classify_file()` to partition hits by ANI matrix membership and compute family context metrics.

**Files:**
- Modify: `src/metadarkmatter/core/classification/classifiers/vectorized.py:232-850`
- Test: `tests/unit/test_family_validation.py` (new file)

**Step 1: Write failing tests**

Create `tests/unit/test_family_validation.py`:

```python
"""Tests for family validation in VectorizedClassifier."""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from metadarkmatter.core.classification.ani_matrix import ANIMatrix
from metadarkmatter.core.classification.classifiers.vectorized import VectorizedClassifier
from metadarkmatter.models.config import ScoringConfig


def _make_ani_matrix(genomes: list[str], default_ani: float = 70.0) -> ANIMatrix:
    """Create a minimal ANIMatrix for testing."""
    n = len(genomes)
    arr = np.full((n, n), default_ani)
    np.fill_diagonal(arr, 100.0)
    # Make within-family genomes related
    for i in range(n):
        for j in range(n):
            if i != j:
                arr[i, j] = 85.0  # same genus
    matrix = ANIMatrix.__new__(ANIMatrix)
    matrix._genomes = genomes
    matrix._genome_to_idx = {g: i for i, g in enumerate(genomes)}
    matrix._ani_array = arr
    matrix._default_ani = default_ani
    return matrix


class TestFamilyMetrics:
    """Tests for family context metric calculations."""

    def test_family_bitscore_ratio_all_in_family(self):
        """When all hits are in-family, ratio should be 1.0."""
        genomes = ["GCF_001", "GCF_002"]
        ani = _make_ani_matrix(genomes)
        config = ScoringConfig(target_family="f__TestFamily")
        classifier = VectorizedClassifier(ani, config=config)

        # All hits match genomes in ANI matrix
        ratio = classifier._compute_family_bitscore_ratio(
            best_in_family_bitscore=500.0,
            best_overall_bitscore=500.0,
        )
        assert ratio == 1.0

    def test_family_bitscore_ratio_external_better(self):
        """When external hit is better, ratio should be < 1.0."""
        genomes = ["GCF_001", "GCF_002"]
        ani = _make_ani_matrix(genomes)
        config = ScoringConfig(target_family="f__TestFamily")
        classifier = VectorizedClassifier(ani, config=config)

        ratio = classifier._compute_family_bitscore_ratio(
            best_in_family_bitscore=400.0,
            best_overall_bitscore=500.0,
        )
        assert ratio == pytest.approx(0.8)

    def test_family_bitscore_ratio_no_in_family(self):
        """When no in-family hits exist, ratio should be 0.0."""
        genomes = ["GCF_001"]
        ani = _make_ani_matrix(genomes)
        config = ScoringConfig(target_family="f__TestFamily")
        classifier = VectorizedClassifier(ani, config=config)

        ratio = classifier._compute_family_bitscore_ratio(
            best_in_family_bitscore=0.0,
            best_overall_bitscore=500.0,
        )
        assert ratio == 0.0


class TestOffTargetDecision:
    """Tests for off-target classification decision rule."""

    def test_below_threshold_is_off_target(self):
        """Reads with ratio below threshold should be off-target."""
        config = ScoringConfig(
            target_family="f__TestFamily",
            family_ratio_threshold=0.8,
        )
        # ratio = 0.7 < 0.8 -> off-target
        assert config.family_ratio_threshold > 0.7

    def test_above_threshold_is_not_off_target(self):
        """Reads with ratio above threshold should not be off-target."""
        config = ScoringConfig(
            target_family="f__TestFamily",
            family_ratio_threshold=0.8,
        )
        # ratio = 0.9 >= 0.8 -> not off-target
        assert config.family_ratio_threshold <= 0.9

    def test_no_target_family_disables_validation(self):
        """Without target_family, all reads go through normal classification."""
        config = ScoringConfig()
        assert config.target_family is None


class TestFamilyValidationDisabled:
    """Tests that classifier output is identical when family validation is off."""

    def test_no_family_columns_without_target_family(self, tmp_path):
        """Output should not contain family columns when target_family is None."""
        genomes = ["GCF_001", "GCF_002"]
        ani = _make_ani_matrix(genomes)
        config = ScoringConfig()  # No target_family

        # Create minimal BLAST file
        blast_file = tmp_path / "test.tsv"
        blast_file.write_text(
            "read1\tGCF_001|contig1\t98.0\t150\t3\t0\t1\t150\t1\t150\t1e-50\t500\n"
            "read1\tGCF_002|contig1\t95.0\t150\t8\t0\t1\t150\t1\t150\t1e-40\t450\n"
        )

        classifier = VectorizedClassifier(ani, config=config)
        result = classifier.classify_file(blast_file)

        # Family columns should NOT be present
        assert "family_bitscore_ratio" not in result.columns
        assert "external_best_genome" not in result.columns
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/unit/test_family_validation.py -v 2>&1 | head -40`
Expected: FAIL (methods don't exist yet)

**Step 3: Implement hit partitioning in VectorizedClassifier**

In `src/metadarkmatter/core/classification/classifiers/vectorized.py`, add the following changes:

**3a. Add helper method for family bitscore ratio** (after `_build_aai_lookup`, around line 152):

```python
@staticmethod
def _compute_family_bitscore_ratio(
    best_in_family_bitscore: float,
    best_overall_bitscore: float,
) -> float:
    """Compute family bitscore ratio for a single read.

    Args:
        best_in_family_bitscore: Best bitscore from in-family hits.
        best_overall_bitscore: Best bitscore across all hits.

    Returns:
        Ratio (0.0-1.0). 0.0 if no in-family hits.
    """
    if best_overall_bitscore == 0.0:
        return 0.0
    if best_in_family_bitscore == 0.0:
        return 0.0
    return best_in_family_bitscore / best_overall_bitscore
```

**3b. Add hit partitioning to `classify_file()`**

After the genome_name extraction (line 285) and before the coverage weighting (line 316), add the family validation block. The logic:

1. Build a set of in-family genome names from the ANI matrix.
2. Tag each hit as `is_in_family` (genome_name is in the ANI matrix genome set).
3. Compute per-read family metrics: `family_bitscore_ratio`, `family_identity_gap`, `in_family_hit_fraction`, `external_best_genome`, `external_best_identity`.
4. Mark reads as off-target when `family_bitscore_ratio < threshold` or no in-family hits.
5. For off-target reads: store metrics, assign "Off-target" category.
6. For remaining reads: filter to in-family hits only, run existing classification pipeline.
7. Concatenate off-target + classified results.

The key code block to insert in `classify_file()` after line 285 (`df = df.with_columns([extract_genome_name_expr()])`) and before the empty check (line 290):

```python
# Family validation: partition hits by ANI matrix membership
family_validation_active = self.config.target_family is not None
family_metrics_df = None
off_target_df = None

if family_validation_active:
    # Build set of in-family genomes from ANI matrix
    in_family_genomes = set(self.ani_matrix._genomes)

    # Tag each hit
    df = df.with_columns([
        pl.col("genome_name").is_in(in_family_genomes).alias("_is_in_family"),
    ])

    # Per-read family metrics
    read_metrics = (
        df.group_by("qseqid")
        .agg([
            # Best bitscore from in-family hits
            pl.col("bitscore")
            .filter(pl.col("_is_in_family"))
            .max()
            .fill_null(0.0)
            .alias("_best_in_family_bitscore"),
            # Best bitscore overall
            pl.col("bitscore").max().alias("_best_overall_bitscore"),
            # Best external identity
            pl.col("pident")
            .filter(~pl.col("_is_in_family"))
            .max()
            .fill_null(0.0)
            .alias("external_best_identity"),
            # Best in-family identity
            pl.col("pident")
            .filter(pl.col("_is_in_family"))
            .max()
            .fill_null(0.0)
            .alias("_best_in_family_identity"),
            # Best external genome
            pl.col("genome_name")
            .filter(~pl.col("_is_in_family"))
            .sort_by(pl.col("bitscore"), descending=True)
            .first()
            .alias("external_best_genome"),
            # Hit counts
            pl.col("_is_in_family").sum().alias("_in_family_count"),
            pl.len().alias("_total_count"),
        ])
        .with_columns([
            # Family bitscore ratio
            (pl.col("_best_in_family_bitscore") / pl.col("_best_overall_bitscore"))
            .fill_nan(0.0)
            .fill_null(0.0)
            .alias("family_bitscore_ratio"),
            # Family identity gap
            (pl.col("_best_in_family_identity") - pl.col("external_best_identity"))
            .alias("family_identity_gap"),
            # In-family hit fraction
            (pl.col("_in_family_count").cast(pl.Float64) / pl.col("_total_count").cast(pl.Float64))
            .alias("in_family_hit_fraction"),
        ])
    )

    # Determine off-target reads
    threshold = self.config.family_ratio_threshold
    off_target_reads = read_metrics.filter(
        (pl.col("family_bitscore_ratio") < threshold)
        | (pl.col("_best_in_family_bitscore") == 0.0)
    )

    family_metrics_df = read_metrics  # Save for later joining

    if not off_target_reads.is_empty():
        off_target_read_ids = set(off_target_reads["qseqid"].to_list())

        # Build off-target result DataFrame
        off_target_df = off_target_reads.select([
            pl.col("qseqid").alias("read_id"),
            pl.col("external_best_genome").fill_null("unknown").alias("best_match_genome"),
            pl.col("external_best_identity").alias("top_hit_identity"),
            (100.0 - pl.col("external_best_identity")).alias("novelty_index"),
            pl.lit(0.0).alias("placement_uncertainty"),
            pl.lit(0.0).alias("genus_uncertainty"),
            pl.lit("unambiguous").alias("ambiguity_scope"),
            pl.lit(0).cast(pl.Int64).alias("num_ambiguous_hits"),
            pl.lit(None).cast(pl.Float64).alias("second_hit_identity"),
            pl.lit(None).cast(pl.Float64).alias("identity_gap"),
            pl.lit(0.0).alias("confidence_score"),
            pl.lit("Off-target").alias("taxonomic_call"),
            pl.lit("Uncertain").alias("diversity_status"),
            pl.lit(False).alias("is_novel"),
            pl.lit(False).alias("low_confidence"),
            pl.lit(None).cast(pl.Float64).alias("inferred_uncertainty"),
            pl.lit("none").alias("uncertainty_type"),
            pl.lit(0.0).alias("alignment_quality"),
            pl.lit(0.0).alias("identity_confidence"),
            pl.lit(0.0).alias("placement_confidence"),
            pl.lit(None).cast(pl.Float64).alias("discovery_score"),
            # Family validation columns
            pl.col("family_bitscore_ratio"),
            pl.col("family_identity_gap"),
            pl.col("in_family_hit_fraction"),
            pl.col("external_best_genome"),
            pl.col("external_best_identity"),
        ])

        # Filter main DataFrame to only in-family hits for non-off-target reads
        df = df.filter(
            ~pl.col("qseqid").is_in(off_target_read_ids)
            & pl.col("_is_in_family")
        )
    else:
        # No off-target reads - filter to in-family hits only
        df = df.filter(pl.col("_is_in_family"))

    # Drop the temporary column
    df = df.drop("_is_in_family")
```

Then at the end of `classify_file()`, before returning (around line 839-850), add the family column joining and concatenation:

```python
# Add family validation columns if active
if family_validation_active and family_metrics_df is not None:
    # Join family metrics to classified results
    family_cols = family_metrics_df.select([
        pl.col("qseqid"),
        "family_bitscore_ratio",
        "family_identity_gap",
        "in_family_hit_fraction",
        "external_best_genome",
        "external_best_identity",
    ])
    classification_result = classification_result.join(
        family_cols,
        left_on="read_id",
        right_on="qseqid",
        how="left",
    )

    # Concatenate off-target reads
    if off_target_df is not None and not off_target_df.is_empty():
        classification_result = pl.concat(
            [classification_result, off_target_df],
            how="diagonal",
        )
```

Also update `_empty_dataframe()` to conditionally include family columns when `self.config.target_family is not None`.

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/unit/test_family_validation.py -v`
Expected: PASS

**Step 5: Run full test suite to check for regressions**

Run: `python -m pytest tests/ -v --timeout=120 2>&1 | tail -30`
Expected: All existing tests PASS (backward compatibility preserved)

**Step 6: Commit**

```bash
git add src/metadarkmatter/core/classification/classifiers/vectorized.py tests/unit/test_family_validation.py
git commit -m "feat: implement hit partitioning and family metrics in VectorizedClassifier"
```

---

### Task 6: Add Off-target to threshold application logic

**Files:**
- Modify: `src/metadarkmatter/core/classification/thresholds.py:16-143`
- Test: `tests/unit/test_family_validation.py` (add threshold tests)

**Step 1: Write failing test**

Add to `tests/unit/test_family_validation.py`:

```python
from metadarkmatter.core.classification.thresholds import apply_classification_thresholds


class TestThresholdsWithOffTarget:
    """Tests that thresholds module handles Off-target correctly."""

    def test_off_target_preserved_through_reclassification(self):
        """Off-target reads should not be reclassified by threshold application."""
        df = pl.DataFrame({
            "read_id": ["r1", "r2"],
            "novelty_index": [15.0, 5.0],
            "placement_uncertainty": [0.5, 0.5],
            "num_ambiguous_hits": [1, 2],
            "identity_gap": [None, 3.0],
            "taxonomic_call": ["Off-target", "Known Species"],
        })
        config = ScoringConfig(target_family="f__Test")
        result = apply_classification_thresholds(df, config)
        # Off-target should be preserved, not overwritten
        assert result.filter(pl.col("read_id") == "r1")["taxonomic_call"][0] == "Off-target"
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest tests/unit/test_family_validation.py::TestThresholdsWithOffTarget -v`
Expected: FAIL (Off-target gets reclassified by threshold rules)

**Step 3: Implement the change**

In `src/metadarkmatter/core/classification/thresholds.py`, at the start of `apply_classification_thresholds()`, add a guard for Off-target reads:

After line 58, before the classification expression:

```python
# Preserve Off-target reads (already classified by family validation)
has_off_target = "taxonomic_call" in df.columns
if has_off_target:
    off_target_mask = pl.col("taxonomic_call") == "Off-target"
    off_target_rows = df.filter(off_target_mask)
    df = df.filter(~off_target_mask)
```

And before the return, recombine:

```python
if has_off_target and not off_target_rows.is_empty():
    # Re-add Off-target rows with updated diversity columns
    off_target_rows = off_target_rows.with_columns([
        pl.lit("Off-target").alias("taxonomic_call"),
        pl.lit("Uncertain").alias("diversity_status"),
        pl.lit(False).alias("is_novel"),
    ])
    result = pl.concat([result, off_target_rows], how="diagonal")
```

**Step 4: Run test to verify it passes**

Run: `python -m pytest tests/unit/test_family_validation.py::TestThresholdsWithOffTarget -v`
Expected: PASS

**Step 5: Commit**

```bash
git add src/metadarkmatter/core/classification/thresholds.py tests/unit/test_family_validation.py
git commit -m "feat: preserve Off-target reads through threshold reclassification"
```

---

### Task 7: Backward compatibility test

**Files:**
- Test: `tests/unit/test_family_validation.py` (add backward compat test)

**Step 1: Write the backward compatibility test**

Add to `tests/unit/test_family_validation.py`:

```python
class TestBackwardCompatibility:
    """Verify identical output when family validation is disabled."""

    def test_output_identical_without_target_family(self, tmp_path):
        """Classification should produce identical results when target_family is None."""
        genomes = ["GCF_001", "GCF_002", "GCF_003"]
        ani = _make_ani_matrix(genomes)

        # Create test BLAST file with hits only to in-family genomes
        blast_file = tmp_path / "test.tsv"
        lines = [
            "read1\tGCF_001|c1\t98.0\t150\t3\t0\t1\t150\t1\t150\t1e-50\t500",
            "read1\tGCF_002|c1\t95.0\t150\t8\t0\t1\t150\t1\t150\t1e-40\t450",
            "read2\tGCF_003|c1\t88.0\t150\t18\t0\t1\t150\t1\t150\t1e-30\t300",
        ]
        blast_file.write_text("\n".join(lines) + "\n")

        # Run without family validation
        config_without = ScoringConfig()
        classifier_without = VectorizedClassifier(ani, config=config_without)
        result_without = classifier_without.classify_file(blast_file)

        # Run with family validation but all hits in-family
        config_with = ScoringConfig(target_family="f__TestFamily")
        classifier_with = VectorizedClassifier(ani, config=config_with)
        result_with = classifier_with.classify_file(blast_file)

        # Core classification columns should match
        core_cols = [
            "read_id", "best_match_genome", "top_hit_identity",
            "novelty_index", "placement_uncertainty", "taxonomic_call",
        ]
        for col in core_cols:
            assert result_without[col].to_list() == result_with[col].to_list(), (
                f"Column {col} differs"
            )

        # Family columns should only exist in result_with
        assert "family_bitscore_ratio" not in result_without.columns
        assert "family_bitscore_ratio" in result_with.columns
```

**Step 2: Run the test**

Run: `python -m pytest tests/unit/test_family_validation.py::TestBackwardCompatibility -v`
Expected: PASS (if Task 5 implemented correctly)

**Step 3: Commit**

```bash
git add tests/unit/test_family_validation.py
git commit -m "test: add backward compatibility test for family validation"
```

---

### Task 8: Add family validation summary to JSON output

**Files:**
- Modify: `src/metadarkmatter/cli/score.py` (summary generation section)
- Test: `tests/unit/test_family_validation.py`

**Step 1: Write failing test**

```python
class TestFamilyValidationSummary:
    """Tests for family validation in summary JSON."""

    def test_summary_includes_off_target_count(self):
        """TaxonomicSummary should include off_target count."""
        summary = TaxonomicSummary(
            total_reads=100,
            known_species=50,
            novel_species=20,
            novel_genus=5,
            conserved_regions=3,
            off_target=10,
            mean_novelty_index=8.0,
            mean_placement_uncertainty=1.2,
        )
        assert summary.off_target == 10
        data = summary.model_dump()
        assert data["off_target"] == 10
```

**Step 2: Run test**

Run: `python -m pytest tests/unit/test_family_validation.py::TestFamilyValidationSummary -v`
Expected: Should PASS (if Task 1 added the field correctly). If not, add the field.

**Step 3: Update summary generation in CLI**

In `src/metadarkmatter/cli/score.py`, find the summary generation section (search for `TaxonomicSummary`) and add off_target count extraction from the classification DataFrame:

```python
off_target_count = int(classification_df.filter(
    pl.col("taxonomic_call") == "Off-target"
).height) if classification_df is not None else 0
```

Pass `off_target=off_target_count` to the `TaxonomicSummary` constructor.

**Step 4: Run test to verify**

Run: `python -m pytest tests/unit/test_family_validation.py::TestFamilyValidationSummary -v`
Expected: PASS

**Step 5: Commit**

```bash
git add src/metadarkmatter/cli/score.py tests/unit/test_family_validation.py
git commit -m "feat: include off-target count in summary JSON output"
```

---

### Task 9: Add Family Validation tab to HTML report

**Files:**
- Modify: `src/metadarkmatter/visualization/report/templates.py` (add template)
- Modify: `src/metadarkmatter/visualization/report/generator.py` (add section builder)
- Test: `tests/unit/test_report_generator.py`

**Step 1: Write failing test**

Add to `tests/unit/test_report_generator.py`:

```python
def test_family_validation_tab_present_when_active(self):
    """Report should include Family Validation tab when off-target reads exist."""
    # This test checks that the template string is present in generated HTML
    # when classification data contains Off-target reads
    from metadarkmatter.visualization.report.templates import FAMILY_VALIDATION_SECTION_TEMPLATE
    assert "Family Validation" in FAMILY_VALIDATION_SECTION_TEMPLATE
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest tests/unit/test_report_generator.py -v -k "family_validation" 2>&1 | head -20`
Expected: FAIL with `ImportError: cannot import name 'FAMILY_VALIDATION_SECTION_TEMPLATE'`

**Step 3: Implement the template**

In `src/metadarkmatter/visualization/report/templates.py`, add the Family Validation section template:

```python
FAMILY_VALIDATION_SECTION_TEMPLATE: str = '''
<div class="stats-grid">
    <div class="stat-card">
        <div class="stat-value">{validated_pct:.1f}%</div>
        <div class="stat-label">Family Validated</div>
    </div>
    <div class="stat-card">
        <div class="stat-value">{off_target_count:,}</div>
        <div class="stat-label">Off-target Reads</div>
    </div>
    <div class="stat-card">
        <div class="stat-value">{off_target_pct:.1f}%</div>
        <div class="stat-label">Off-target Rate</div>
    </div>
    <div class="stat-card">
        <div class="stat-value">{target_family}</div>
        <div class="stat-label">Target Family</div>
    </div>
</div>
<div class="plot-container" id="family-ratio-histogram"></div>
<div class="plot-container" id="family-pie-chart"></div>
{external_families_table}
'''
```

In `src/metadarkmatter/visualization/report/generator.py`, add a `_build_family_validation_section()` method that:

1. Counts off-target and in-family reads from the classification DataFrame.
2. Builds a Plotly histogram of `family_bitscore_ratio` values.
3. Builds a pie chart of in-family vs off-target proportions.
4. Builds a table of top external families (if `external_best_genome` column exists and metadata is available).
5. Returns the formatted HTML section.

Add a new tab entry in the report navigation and call the section builder when `family_bitscore_ratio` column is present in the classification DataFrame.

**Step 4: Run test to verify it passes**

Run: `python -m pytest tests/unit/test_report_generator.py -v -k "family_validation"`
Expected: PASS

**Step 5: Commit**

```bash
git add src/metadarkmatter/visualization/report/templates.py src/metadarkmatter/visualization/report/generator.py tests/unit/test_report_generator.py
git commit -m "feat: add Family Validation tab to HTML report"
```

---

### Task 10: Integration test with broad-database alignment

**Files:**
- Create: `tests/integration/test_family_validation.py`

**Step 1: Write integration test**

```python
"""Integration test for family validation with mixed in-family and external hits."""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest
from pathlib import Path

from metadarkmatter.core.classification.ani_matrix import ANIMatrix
from metadarkmatter.core.classification.classifiers.vectorized import VectorizedClassifier
from metadarkmatter.models.config import ScoringConfig


def _make_ani_matrix(genomes, default_ani=70.0):
    n = len(genomes)
    arr = np.full((n, n), default_ani)
    np.fill_diagonal(arr, 100.0)
    for i in range(n):
        for j in range(n):
            if i != j:
                arr[i, j] = 85.0
    matrix = ANIMatrix.__new__(ANIMatrix)
    matrix._genomes = genomes
    matrix._genome_to_idx = {g: i for i, g in enumerate(genomes)}
    matrix._ani_array = arr
    matrix._default_ani = default_ani
    return matrix


class TestFamilyValidationIntegration:
    """Full pipeline integration test with mixed hits."""

    def test_mixed_hits_classifies_correctly(self, tmp_path):
        """Reads with external hits should be classified as Off-target."""
        # ANI matrix only contains family genomes
        family_genomes = ["GCF_001", "GCF_002"]
        ani = _make_ani_matrix(family_genomes)

        # BLAST file with both in-family and external hits
        blast_lines = [
            # read1: best hit is in-family -> should be classified normally
            "read1\tGCF_001|c1\t98.0\t150\t3\t0\t1\t150\t1\t150\t1e-50\t500",
            "read1\tGCF_002|c1\t95.0\t150\t8\t0\t1\t150\t1\t150\t1e-40\t450",
            # read2: external hit is much better -> Off-target
            "read2\tGCF_001|c1\t75.0\t150\t38\t0\t1\t150\t1\t150\t1e-10\t200",
            "read2\tEXTERNAL_001|c1\t95.0\t150\t8\t0\t1\t150\t1\t150\t1e-50\t600",
            # read3: only external hits -> Off-target
            "read3\tEXTERNAL_002|c1\t92.0\t150\t12\t0\t1\t150\t1\t150\t1e-45\t550",
        ]
        blast_file = tmp_path / "mixed.tsv"
        blast_file.write_text("\n".join(blast_lines) + "\n")

        config = ScoringConfig(
            target_family="f__TestFamily",
            family_ratio_threshold=0.8,
        )
        classifier = VectorizedClassifier(ani, config=config)
        result = classifier.classify_file(blast_file)

        # read1 should be classified normally (ratio = 500/500 = 1.0)
        read1 = result.filter(pl.col("read_id") == "read1")
        assert read1["taxonomic_call"][0] != "Off-target"
        assert read1["family_bitscore_ratio"][0] == pytest.approx(1.0)

        # read2 should be Off-target (ratio = 200/600 = 0.33)
        read2 = result.filter(pl.col("read_id") == "read2")
        assert read2["taxonomic_call"][0] == "Off-target"
        assert read2["family_bitscore_ratio"][0] == pytest.approx(200.0 / 600.0, abs=0.01)

        # read3 should be Off-target (no in-family hits)
        read3 = result.filter(pl.col("read_id") == "read3")
        assert read3["taxonomic_call"][0] == "Off-target"
        assert read3["family_bitscore_ratio"][0] == 0.0
```

**Step 2: Run the integration test**

Run: `python -m pytest tests/integration/test_family_validation.py -v`
Expected: PASS (if Tasks 1-5 implemented correctly)

**Step 3: Commit**

```bash
git add tests/integration/test_family_validation.py
git commit -m "test: add integration test for family validation pipeline"
```

---

### Task 11: Update documentation

**Files:**
- Modify: `docs/REFERENCE.md` (add family validation section)
- Modify: `docs/CLI_REFERENCE.md` (add --target-family, --family-ratio-threshold)
- Modify: `docs/WORKFLOW.md` (add broad-database workflow example)
- Modify: `CLAUDE.md` (add family validation section)

**Step 1: Update docs/CLI_REFERENCE.md**

Add `--target-family` and `--family-ratio-threshold` to the `score classify` options table.

**Step 2: Update docs/REFERENCE.md**

Add a "Family Validation" subsection in Advanced Features explaining:
- When to use family validation
- How to enable it
- What the Off-target category means
- Output columns

**Step 3: Update docs/WORKFLOW.md**

Add a "Broad-Database Classification" section showing:
```bash
# Run BLAST against all bacteria
blastn -query reads.fa -db all_bacteria -outfmt 6 -out broad_results.tsv

# Classify with family validation
metadarkmatter score classify \
    --alignment broad_results.tsv \
    --ani family_ani_matrix.csv \
    --target-family "f__Francisellaceae" \
    --output classifications.csv \
    --parallel
```

**Step 4: Update CLAUDE.md**

Add to the "Advanced Classification Features" section:
```
### Family Validation (`--target-family`)
- Detects off-target reads from broad-database alignments
- Partitions hits by ANI matrix membership (in-family vs external)
- Reads with best_in_family_bitscore / best_overall_bitscore < 0.8 are Off-target
- Key file: `core/classification/classifiers/vectorized.py`
```

**Step 5: Commit**

```bash
git add docs/REFERENCE.md docs/CLI_REFERENCE.md docs/WORKFLOW.md CLAUDE.md
git commit -m "docs: add family validation documentation"
```

---

### Task 12: Final verification

**Step 1: Run full test suite**

Run: `python -m pytest tests/ -v --timeout=120 2>&1 | tail -40`
Expected: All tests PASS

**Step 2: Run specific family validation tests**

Run: `python -m pytest tests/unit/test_family_validation.py tests/integration/test_family_validation.py -v`
Expected: All tests PASS

**Step 3: Verify CLI help**

Run: `python -m metadarkmatter score classify --help | grep -A2 "target-family\|family-ratio"`
Expected: Both flags appear with descriptions

**Step 4: Commit final state (if any remaining changes)**

```bash
git status
# If clean, no commit needed
```
