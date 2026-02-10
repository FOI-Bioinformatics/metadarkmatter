# Classification Pipeline Improvements: Design Document

> **Date:** 2026-02-10
> **Status:** Draft - Pending Review
> **Scope:** Five interconnected improvements to the classification stage, plus pipeline-wide QC

---

## Executive Summary

The metadarkmatter classification algorithm is well-implemented and performs reliably for typical use cases. However, three interconnected limitations reduce accuracy and interpretability for challenging environmental samples:

1. **Fixed species boundaries** (95% ANI) ignore taxon-specific variation
2. **Single-hit reads** (~70% of environmental data) produce overconfident classifications
3. **Conserved gene fragments** inflate identity scores, biasing toward "Known Species"

This document proposes five improvements that address these issues while preserving backward compatibility. They are ordered from quick wins to larger architectural changes:

| # | Improvement | Effort | Impact |
|---|------------|--------|--------|
| 1 | Automated QC checks | Days | High |
| 2 | Coverage weighting discoverability | Hours | Medium |
| 3 | Threshold sensitivity tab in report | ~1 week | Medium |
| 4 | Adaptive thresholds from ANI distribution | 1-2 weeks | High |
| 5 | Bayesian confidence framework | 2-4 weeks | Very high |

---

## Current State Analysis

### Classification Decision Flow

The current algorithm evaluates reads through a deterministic decision tree
(`vectorized.py` lines 611-692):

```
Read → Best hit selection → Metric computation → Hard threshold → Category
```

Key metrics computed per read:
- `novelty_index = 100 - top_hit_identity`
- `placement_uncertainty = 100 - max(ANI between competing genomes)`
- `identity_gap = top_identity - second_identity`
- `num_ambiguous_hits` (within 95% of top bitscore)
- `alignment_quality` (length, coverage, e-value)

### Three Core Problems

**Problem 1: Fixed Thresholds**

The species boundary sits at 95% ANI (novelty_index = 5%). This value comes
from Jain et al. (2018) and works for most prokaryotes, but:

- Some genera cluster tightly (e.g., Bacillus cereus group: 97-99.9% ANI
  within species, gaps at 94-95%)
- Others are more diffuse (e.g., Prochlorococcus: species separated by
  only 85-90% ANI)
- The ANI matrix already contains the information needed to detect these
  taxon-specific boundaries

**Problem 2: Single-Hit Reads**

For reads with exactly one BLAST hit:
- `placement_uncertainty = 0%` (no competing genomes to compare)
- The read is classified solely by novelty_index
- A read at 92% identity is called "Novel Species" with maximum confidence
- But we have *no evidence* that the placement is unambiguous - there
  simply were no other genomes close enough to score

Current workaround: `--infer-single-hit-uncertainty` estimates uncertainty
from the novelty level, but this is heuristic and disabled by default.

**Problem 3: Conserved Gene Bias**

Highly conserved genes (16S rRNA, 23S rRNA, recA, rpoB) are present in all
bacteria and share >95% identity even across phyla. When a read from a novel
organism aligns to a conserved gene in the reference database:
- Alignment length: ~200-500bp (short fragment of conserved gene)
- Percent identity: 95-99% (gene conservation, not genomic similarity)
- Classification: "Known Species" (false positive)

Coverage weighting (`--coverage-weight-mode`) addresses this by penalizing
short alignments, but it defaults to `none` for backward compatibility
and most users are not aware it exists.

---

## Improvement 1: Automated QC Checks

### Motivation

The pipeline currently provides no warning when input data quality is poor
or results are suspicious. Users discover problems only when examining
the final report.

### Design

Add validation gates at three points in the pipeline:

#### Gate 1: Pre-Classification Input QC

**Location:** `VectorizedClassifier.classify_file()` after loading data,
before classification (~line 256)

**Checks:**

| Check | Condition | Severity | Message |
|-------|-----------|----------|---------|
| Low hit rate | <20% of reads have BLAST hits | Warning | "Only {pct}% of reads aligned. Consider a broader reference database." |
| Short alignments | Median alignment length <200bp | Warning | "Median alignment is {len}bp. Short conserved gene fragments may dominate. Consider --coverage-weight-mode linear." |
| Extreme e-values | >10% of hits with e-value >1e-5 | Warning | "Many weak hits detected. Consider stricter BLAST parameters." |
| Empty input | 0 alignments parsed | Error | "No valid alignments found in {file}." |
| Missing ANI pairs | >10% of genome pairs missing | Warning | "ANI matrix is {pct}% incomplete. Uncertainty estimates may be unreliable." |

#### Gate 2: Post-Classification Consistency QC

**Location:** After classification assignment (~line 692)

**Checks added as columns:**

```python
result = result.with_columns([
    # Flag: very high novelty classified as Known Species
    ((pl.col("novelty_index") > 10)
     & (pl.col("taxonomic_call") == "Known Species"))
    .alias("qc_inconsistent_novelty"),

    # Flag: single-hit novel call (no placement evidence)
    ((pl.col("num_ambiguous_hits") <= 1)
     & (pl.col("taxonomic_call").is_in(["Novel Species", "Novel Genus"])))
    .alias("qc_single_hit_novel"),

    # Flag: very short alignment for confident classification
    ((pl.col("alignment_quality") < 0.3)
     & (pl.col("taxonomic_call") != "Ambiguous"))
    .alias("qc_short_alignment"),
])
```

#### Gate 3: Summary-Level QC

**Location:** After all reads classified, in CLI output and report.

| Check | Condition | Action |
|-------|-----------|--------|
| Single-hit dominance | >70% reads with 1 hit | Recommend `--infer-single-hit-uncertainty` |
| High ambiguous rate | >40% reads "Ambiguous" | Suggest reviewing reference set completeness |
| Conserved gene signal | >30% reads with alignment <300bp and identity >95% | Recommend `--coverage-weight-mode linear` |
| Classification imbalance | >90% one category | Note potential database bias |

### Implementation

**New file:** `src/metadarkmatter/core/classification/qc.py`

```python
@dataclass
class ClassificationQC:
    """Quality control metrics for a classification run."""

    total_reads: int
    alignment_rate: float           # Fraction of input reads with hits
    median_alignment_length: float
    single_hit_fraction: float
    short_alignment_fraction: float # <300bp
    ambiguous_fraction: float
    qc_flags: list[QCFlag]          # List of triggered warnings

    def summary(self) -> str:
        """Human-readable QC summary."""

    def to_dict(self) -> dict[str, Any]:
        """Serializable dict for report integration."""
```

**Integration in CLI** (`cli/score.py`): Display QC summary after classification
table. Include in report if `--output-qc` flag provided.

**Integration in report** (`visualization/report/generator.py`): Add QC
metrics to the Overview tab as a collapsible panel.

### Backward Compatibility

QC checks are purely additive. No existing behavior changes. QC columns
are added to the output DataFrame but do not affect classification logic.

---

## Improvement 2: Coverage Weighting Discoverability

### Motivation

Coverage weighting already solves the conserved gene bias problem but most
users do not know it exists. The feature is buried in CLI options and
defaults to `none`.

### Design

Three changes to surface the feature:

#### 2a: Proactive Recommendation

When the post-classification QC detects a conserved gene signal (>30% of
top hits are short alignments with high identity), print:

```
Note: 42% of top hits are short conserved gene fragments (<300bp, >95% identity).
This may inflate "Known Species" counts. To prioritize full-length alignments:

  metadarkmatter score classify ... --coverage-weight-mode linear

See docs/REFERENCE.md#coverage-weighting for details.
```

#### 2b: Add to Report

Include alignment length distribution histogram in the Distributions tab.
Annotate the 300bp mark with a vertical line and note.

#### 2c: Consider Changing Default

**Longer term:** After gathering user feedback, consider changing the default
from `none` to `linear` with `strength=0.3` (mild preference for longer
alignments). This would require a deprecation cycle:

1. Release N: Add deprecation warning when `none` is used implicitly
2. Release N+1: Change default to `linear`

This is a discussion point, not a firm recommendation.

---

## Improvement 3: Threshold Sensitivity Tab in Report

### Motivation

When a user gets "35% Novel Species," they need to know: is this robust or
does it swing to 10% if the threshold moves by 1%? Currently there is
no way to assess this without re-running the classification manually.

### Design

Add a "Sensitivity" tab to the HTML report showing how classification
counts change across a range of threshold values.

#### Visualization

**Primary plot:** Stacked area chart with novelty threshold on x-axis
(90-100% identity in 0.5% steps) and read counts on y-axis, colored by
classification category.

**Secondary plot:** Line chart showing the number of novel species calls
as a function of the species boundary threshold.

#### Computation

Re-classify all reads at each threshold without re-running the full
pipeline. Since novelty_index and placement_uncertainty are already
computed, only the threshold comparison changes:

```python
def compute_sensitivity(
    df: pl.DataFrame,
    thresholds: list[float],
) -> pl.DataFrame:
    """Re-classify reads at multiple threshold values.

    Only varies novelty_known_max (species boundary).
    Uncertainty thresholds remain fixed.

    Returns DataFrame with columns: threshold, known, novel_species,
    novel_genus, ambiguous, conserved.
    """
    results = []
    for threshold in thresholds:
        counts = _classify_at_threshold(df, threshold)
        results.append({"threshold": threshold, **counts})
    return pl.DataFrame(results)
```

#### Report Integration

Add as a sub-panel in the Overview or Distributions tab. Include brief
interpretation text:

> "The threshold sensitivity plot shows how classification proportions
> change as the species boundary moves. A steep slope near the current
> threshold (dashed line) indicates many reads near the boundary."

### Implementation

**New file:** `src/metadarkmatter/core/classification/sensitivity.py`
**Modified:** `visualization/report/generator.py` (new section builder)
**Modified:** `visualization/report/templates.py` (new template)

---

## Improvement 4: Adaptive Thresholds from ANI Distribution

### Motivation

The 95% ANI species boundary is a useful default but does not reflect the
actual taxonomic structure in every reference set. The ANI matrix already
encodes this structure: within-species pairs cluster at high ANI values,
between-species pairs cluster lower, and the species boundary sits in
the gap between them.

### Approach: ANI Gap Detection

#### Step 1: Extract Off-Diagonal ANI Values

From the symmetric ANI matrix, collect all unique pairwise values
(excluding the diagonal):

```python
import numpy as np
from scipy import stats

values = ani_matrix.values[np.triu_indices(n, k=1)]
values = values[~np.isnan(values)]
```

#### Step 2: Fit Gaussian Mixture Model

Fit a 2-component GMM to the ANI distribution:

```python
from sklearn.mixture import GaussianMixture

gmm = GaussianMixture(n_components=2, random_state=42)
gmm.fit(values.reshape(-1, 1))

# Component with higher mean = within-species
# Component with lower mean = between-species
means = sorted(gmm.means_.flatten())
between_species_mean = means[0]
within_species_mean = means[1]
```

#### Step 3: Find Crossing Point

The adaptive species boundary is where the two Gaussian densities cross:

```python
def find_crossing_point(gmm, low=70, high=100, n_points=1000):
    """Find ANI value where within-species and between-species densities cross."""
    x = np.linspace(low, high, n_points).reshape(-1, 1)
    log_probs = gmm.predict_proba(x)

    # Find where component probabilities are closest to equal
    diff = np.abs(log_probs[:, 0] - log_probs[:, 1])
    crossing_idx = np.argmin(diff)
    return x[crossing_idx, 0]
```

#### Step 4: Validate and Apply

```python
def detect_species_boundary(ani_matrix: pd.DataFrame) -> AdaptiveBoundary:
    """Detect the species boundary from ANI distribution.

    Returns:
        AdaptiveBoundary with:
        - boundary: ANI value at species boundary (e.g., 95.5)
        - confidence: How clear the bimodal separation is (0-1)
        - within_species_mean: Mean ANI for within-species pairs
        - between_species_mean: Mean ANI for between-species pairs
        - method: "gmm" or "fallback"
    """
    # ... GMM fitting ...

    # Validate: is there a clear bimodal gap?
    separation = within_species_mean - between_species_mean
    if separation < 5.0:
        # No clear bimodal structure - fall back to 95%
        return AdaptiveBoundary(
            boundary=95.0,
            confidence=0.0,
            method="fallback",
            note="No clear species gap detected in ANI distribution.",
        )

    # Compute boundary confidence based on separation and sample size
    confidence = min(1.0, separation / 15.0)  # Normalized by expected separation

    return AdaptiveBoundary(
        boundary=crossing_point,
        confidence=confidence,
        within_species_mean=within_species_mean,
        between_species_mean=between_species_mean,
        method="gmm",
    )
```

#### Integration

**CLI flag:** `--adaptive-thresholds` (opt-in initially)

When enabled:
1. Analyze ANI matrix before classification
2. Detect species boundary
3. Adjust `novelty_known_max` and related thresholds accordingly
4. Report detected boundary in output and report

**Report integration:** Show ANI distribution histogram with detected
boundary as a vertical line, GMM component curves overlaid.

### Edge Cases

| Scenario | Behavior |
|----------|----------|
| <10 genome pairs | Fall back to 95%, warn "too few genomes for adaptive thresholds" |
| Unimodal distribution | Fall back to 95%, note "no clear species gap" |
| >2 modes | Use 2-component GMM, report uncertainty |
| Very tight gap (e.g., 96.5%) | Use detected value, confidence reflects narrowness |
| Very wide gap (e.g., 88%) | Use detected value - this is a genuine loose-species taxon |

### New Dependencies

- `scikit-learn>=1.3.0` (GaussianMixture) - or implement simple GMM from
  scratch to avoid the dependency

### Backward Compatibility

Disabled by default. When `--adaptive-thresholds` is not specified,
behavior is identical to current version.

---

## Improvement 5: Bayesian Confidence Framework

### Motivation

The current classification uses hard thresholds that produce binary
decisions. A read at 94.9% identity is "Novel Species" while one at 95.1%
is "Known Species" - despite the measurements being essentially identical.
Bayesian inference replaces this knife-edge with a continuous probability
distribution that naturally accounts for uncertainty.

### Statistical Framework

#### The Question

For each read, we want to compute:

```
P(category | identity, uncertainty, alignment_quality, n_hits, ANI_landscape)
```

Where category is one of: Known Species, Novel Species, Novel Genus,
Ambiguous, Conserved Region.

#### Likelihood Model

The identity of a read depends on whether it truly belongs to a known
species or represents something novel. We model the identity distribution
for each category as a truncated normal:

```
P(identity | known_species) ~ TruncNorm(mu=97, sigma=2, low=95, high=100)
P(identity | novel_species) ~ TruncNorm(mu=90, sigma=4, low=80, high=95)
P(identity | novel_genus)   ~ TruncNorm(mu=82, sigma=3, low=75, high=85)
```

The parameters can be estimated from the data itself: if we have a
reference set with known within-species ANI values, the empirical
distribution of those values gives us the likelihood directly.

#### Prior Model

The prior probability of novelty depends on:

1. **Database completeness**: If the reference set contains all known
   species in the family, the prior for "Novel Species" is lower than
   if the database is sparse.

2. **ANI landscape**: The density of reference genomes near the read's
   identity range. A read at 92% in a region with many references at
   91-93% is more likely to be a known species variant than the same
   read in a region with no references between 85-95%.

```python
def compute_novelty_prior(
    novelty_index: float,
    ani_matrix: np.ndarray,
    genome_count: int,
) -> float:
    """Estimate prior probability of novelty.

    Based on database density in the novelty range. More reference
    genomes near the identity level → lower novelty prior.
    """
    # Count reference pairs near this ANI level
    target_ani = 100 - novelty_index
    nearby = np.sum(np.abs(ani_matrix - target_ani) < 2.0)
    density = nearby / (genome_count * (genome_count - 1) / 2)

    # Higher density → lower novelty prior
    # Base rate ~20% novel in environmental samples (adjustable)
    base_rate = 0.20
    adjusted = base_rate * (1 - density * 5)  # Density penalizes novelty
    return np.clip(adjusted, 0.01, 0.95)
```

#### Posterior Computation

```python
def compute_posterior(
    identity: float,
    uncertainty: float,
    alignment_quality: float,
    n_hits: int,
    prior_novel: float,
) -> dict[str, float]:
    """Compute posterior probability for each classification category.

    Returns dict with keys: known_species, novel_species, novel_genus,
    ambiguous, conserved_region. Values sum to 1.0.
    """
    # Compute likelihoods
    L_known = likelihood_known(identity, alignment_quality)
    L_novel_sp = likelihood_novel_species(identity, alignment_quality)
    L_novel_gn = likelihood_novel_genus(identity, alignment_quality)

    # Adjust for evidence quality
    evidence_weight = evidence_quality(n_hits, uncertainty, alignment_quality)

    # Single-hit penalty: wider posterior (less confident)
    if n_hits <= 1:
        # Spread probability mass toward ambiguous
        spread = 0.3 * (1 - evidence_weight)
    else:
        spread = 0.0

    # Compute unnormalized posteriors
    posteriors = {
        "known_species": L_known * (1 - prior_novel),
        "novel_species": L_novel_sp * prior_novel * 0.7,
        "novel_genus": L_novel_gn * prior_novel * 0.3,
    }

    # Add ambiguous mass
    posteriors["ambiguous"] = spread * sum(posteriors.values())

    # Conserved region: high uncertainty signal
    if uncertainty > 3.0:
        conserved_mass = uncertainty / 20.0
        posteriors["conserved_region"] = conserved_mass * sum(posteriors.values())

    # Normalize
    total = sum(posteriors.values())
    return {k: v / total for k, v in posteriors.items()}
```

#### Classification from Posterior

```python
def classify_from_posterior(
    posteriors: dict[str, float],
    confidence_threshold: float = 0.6,
) -> tuple[str, float]:
    """Assign classification from posterior probabilities.

    Args:
        posteriors: Category probabilities (sum to 1.0).
        confidence_threshold: Minimum posterior for confident call.

    Returns:
        (category, confidence) tuple.
    """
    best_category = max(posteriors, key=posteriors.get)
    best_prob = posteriors[best_category]

    if best_prob < confidence_threshold:
        return ("Ambiguous", best_prob)

    return (best_category, best_prob)
```

### Integration with Existing Architecture

The Bayesian framework would be an alternative scoring path, not a
replacement. Users select it with `--scoring-mode bayesian`:

```
--scoring-mode threshold    (default, current behavior)
--scoring-mode bayesian     (posterior probability classification)
```

#### Implementation Location

**New file:** `src/metadarkmatter/core/classification/bayesian.py`

Contains:
- `BayesianScorer` class
- Likelihood functions
- Prior computation
- Posterior calculation
- Parameter estimation from ANI matrix

**Modified:** `core/classification/classifiers/vectorized.py`
- After metric computation (line 550), branch based on scoring mode
- Threshold mode: existing logic (lines 611-692)
- Bayesian mode: call `BayesianScorer.score(metrics_df)` then
  `classify_from_posterior()`

**Output additions:**
- `posterior_known`: P(Known Species | data)
- `posterior_novel_species`: P(Novel Species | data)
- `posterior_novel_genus`: P(Novel Genus | data)
- `posterior_max`: max(posteriors) = classification confidence
- `posterior_category`: argmax(posteriors) = final call

### Report Integration

The Bayesian output enables richer visualization:
- Posterior probability distributions per read (violin plots by category)
- Confidence calibration plot (predicted vs observed novelty rate)
- Uncertainty decomposition (how much uncertainty from single-hit vs.
  low alignment quality vs. genuinely ambiguous placement)

### Calibration

The Bayesian model has parameters (likelihood distributions, prior rates)
that should be calibrated. Two approaches:

1. **Empirical Bayes:** Estimate parameters from the current dataset's ANI
   matrix and alignment distribution. Self-contained, no external data.

2. **Reference calibration:** Use well-studied families with known ground
   truth (e.g., Enterobacteriaceae) to set baseline parameters.

The empirical Bayes approach is recommended for the initial implementation
as it requires no external data.

### Backward Compatibility

Entirely opt-in. Default `--scoring-mode threshold` produces identical
results to current version.

---

## Implementation Order and Dependencies

```
Improvement 1: QC Checks
    └── No dependencies, can start immediately
    └── Provides data for Improvement 2 (proactive recommendations)

Improvement 2: Coverage Weighting Discoverability
    └── Depends on: QC Checks (for proactive recommendations)
    └── Simple changes to CLI output and report

Improvement 3: Threshold Sensitivity Tab
    └── No dependencies on other improvements
    └── Provides visualization foundation for Improvement 4

Improvement 4: Adaptive Thresholds
    └── Depends on: scipy (already installed), optionally scikit-learn
    └── Benefits from: Sensitivity tab (shows impact of adapted thresholds)

Improvement 5: Bayesian Confidence
    └── Can use adaptive thresholds as informative priors
    └── Benefits from: QC data as input features
    └── Most impactful but largest effort
```

Recommended implementation order: 1 → 2 → 3 → 4 → 5

---

## Testing Strategy

### Improvement 1 (QC Checks)
- Unit tests: Each QC check with synthetic data triggering/not triggering
- Integration: Full classification with poor-quality inputs, verify warnings
- Regression: Existing tests must pass unchanged (QC is additive)

### Improvement 2 (Coverage Discoverability)
- Unit tests: Recommendation logic triggers at correct thresholds
- Integration: CLI output contains recommendations when appropriate

### Improvement 3 (Sensitivity Tab)
- Unit tests: Re-classification at different thresholds produces correct counts
- Visual: Manual review of generated plots

### Improvement 4 (Adaptive Thresholds)
- Unit tests: GMM fitting on synthetic bimodal distributions
- Edge cases: Unimodal, trimodal, small sample, perfect separation
- Integration: Full classification with adaptive vs fixed, compare results
- Validation: Known taxa (e.g., E. coli vs Shigella at ~96% ANI)

### Improvement 5 (Bayesian Confidence)
- Unit tests: Posterior computation, likelihood functions, prior model
- Calibration tests: Posterior probabilities sum to 1.0 for all inputs
- Comparison tests: Bayesian and threshold modes agree on clear-cut reads
- Sensitivity: Bayesian mode correctly produces lower confidence for
  single-hit reads vs. multi-hit reads

### Benchmark Dataset

A synthetic benchmark should be created with:
- Known ground truth (simulated reads from known genomes)
- Controlled novel organisms (known divergence from reference)
- Single-hit reads (remove some reference genomes to create gaps)
- Conserved gene fragments (extract 16S sequences)

This benchmark enables validation of all five improvements.

---

## Open Questions

1. **scikit-learn dependency:** The adaptive threshold GMM requires
   scikit-learn. Should we implement a minimal GMM from scratch to
   avoid the dependency, or accept it as an optional dependency?

2. **Default coverage weighting:** Should the default change from `none`
   to `linear` in a future release? This affects backward compatibility.

3. **Bayesian prior calibration:** Empirical Bayes (self-calibrating)
   vs. pre-trained parameters from reference datasets?

4. **QC column naming:** Should QC flags be prefixed (`qc_`) or kept
   in a separate output file to avoid cluttering the main output?

5. **Report size:** The sensitivity analysis and Bayesian posteriors add
   data to the report. Set a size budget or let it grow?
