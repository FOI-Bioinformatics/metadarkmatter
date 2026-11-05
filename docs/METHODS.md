# Computational Methods

This document describes the computational methods used by metadarkmatter for detecting novel bacterial diversity from metagenomic sequencing data. All calculations are presented with their mathematical formulations and biological interpretations.

**Version:** Corresponds to metadarkmatter v0.2
**Last updated:** 2026-11-02 (re-checked against the classifier code after PR #13)

---

## Table of Contents

1. [Overview](#1-overview)
2. [Read Classification Algorithm](#2-read-classification-algorithm)
   - 2.1 [Novelty Index](#21-novelty-index-n)
   - 2.2 [Placement Uncertainty](#22-placement-uncertainty-u)
   - 2.3 [Classification Decision Tree](#23-classification-decision-tree)
   - 2.4 [Threshold Justification](#24-threshold-justification)
   - 2.5 [Hard Thresholds versus Bayesian Decision Boundaries](#25-hard-thresholds-versus-bayesian-decision-boundaries)
3. [Enhanced Scoring Metrics](#3-enhanced-scoring-metrics)
   - 3.1 [Confidence Score](#31-confidence-score)
   - 3.2 [Inferred Uncertainty](#32-inferred-uncertainty-for-single-hit-reads)
   - 3.3 [Alignment Quality](#33-alignment-quality-score)
   - 3.4 [Identity and Placement Confidence](#34-identity-and-placement-confidence)
   - 3.5 [Discovery Score](#35-discovery-score)
4. [Novel Diversity Clustering](#4-novel-diversity-clustering)
   - 4.1 [Clustering Algorithm](#41-clustering-algorithm)
   - 4.2 [Cluster Confidence Ratings](#42-cluster-confidence-ratings)
   - 4.3 [Estimated ANI Calculation](#43-estimated-ani-calculation)
5. [Coverage Weighting](#5-coverage-weighting)
6. [Protein Mode Thresholds](#6-protein-mode-thresholds)
7. [Quality Control Metrics](#7-quality-control-metrics)
   - 7.1 [Pre-Classification QC](#71-pre-classification-qc)
   - 7.2 [Post-Classification QC](#72-post-classification-qc)
8. [Threshold Sensitivity Analysis](#8-threshold-sensitivity-analysis)
9. [Adaptive Thresholds](#9-adaptive-thresholds-from-ani-distribution)
10. [Bayesian Classification Framework](#10-bayesian-classification-framework)
    - 10.1 [Likelihood Model](#101-likelihood-model)
    - 10.2 [Posterior Computation](#102-posterior-computation)
    - 10.3 [Stage 2 Refinement](#103-stage-2-refinement)
    - 10.4 [Shannon Entropy and Confidence](#104-shannon-entropy-and-confidence)
    - 10.5 [Clade-Novelty Rule](#105-clade-novelty-rule)
    - 10.6 [Interpretation](#106-interpretation)
11. [References](#11-references)

---

## 1. Overview

metadarkmatter classifies metagenomic reads by combining sequence alignment scores with pre-computed genome similarity matrices. The core innovation is **ANI-weighted placement uncertainty**, which distinguishes confident identification of novel taxa from ambiguous placement within known taxa.

**Input data:**
- BLAST/MMseqs2 alignment results (tabular format)
- ANI matrix (Average Nucleotide Identity between reference genomes)
- Optional: AAI matrix (Average Amino Acid Identity) for protein-level analysis

**Output:**
- Per-read classification with novelty and uncertainty metrics
- Confidence scores and discovery prioritization
- Clustered novel diversity analysis

---

## 2. Read Classification Algorithm

### 2.1 Novelty Index (N)

The Novelty Index quantifies sequence divergence from the closest reference genome.

**Definition:**

```
N = 100 - pident_top
```

Where:
- `pident_top` = percent identity of the best hit, i.e. the hit with the highest coverage-weighted bitscore (Section 5; raw bitscore when `coverage_weight_mode="none"`). Ties are broken by higher `pident`, then by genome accession (ascending), so the choice is deterministic.
- `pident_top` is clipped to [0, 100]; N is expressed as a percentage (0-100%)

Before the best hit is chosen, hits are filtered (`_apply_alignment_filters` in `core/classification/classifiers/vectorized.py`):
- Reverse-strand coordinates (`qstart > qend`, as reported by MMseqs2/DIAMOND) are normalised to `qstart = min(qstart, qend)`, `qend = max(qstart, qend)`.
- `length >= min_alignment_length` (default 100 bp). In protein mode the query span in nucleotides (`qend - qstart + 1`) is compared instead, because BLASTX/DIAMOND report `length` in amino acids.
- Alignment fraction `(qend - qstart + 1) / read_length >= min_alignment_fraction` (default 0.3), where `read_length = qlen` when the alignment file provides it, otherwise `qend` is used as a proxy (12-column input).
- Optional E-value, identity, bitscore, read-length and query-coverage filters (all disabled by default).

**Interpretation:**

| Novelty Index | BLAST Identity | Biological Interpretation |
|---------------|----------------|---------------------------|
| N < 4% | pident > 96% | At or above species boundary |
| 4% ≤ N < 20% | 80% < pident ≤ 96% | Below species boundary (novel species candidate) |
| 20% ≤ N ≤ 25% | 75% ≤ pident ≤ 80% | Genus-level divergence (novel genus candidate) |
| N > 25% | pident < 75% | Very high divergence |

**Note:** Read-level BLAST identity can be 10-20% lower than genome-level ANI due to partial alignments. The 20% threshold for novel genus accounts for this systematic offset.

### 2.2 Placement Uncertainty (U)

Placement Uncertainty measures classification ambiguity when a read aligns to multiple reference genomes with similar scores.

**Step 1: Identify competing hits**

```
threshold = (bitscore_threshold_pct / 100) × weighted_bitscore_top      # default 0.95
ambiguous_hits = {hit : hit.weighted_bitscore ≥ threshold}
num_ambiguous_hits = |ambiguous_hits|                                   # counts hits, not genomes
```

Each hit genome is mapped to its species representative (from `genome_metadata.tsv`; every genome maps to itself when no `representative` column is present). ANI lookups use representatives.

**Step 2: Calculate uncertainty from ANI**

```
secondary = representatives of ambiguous_hits, excluding rep(G_top)

if num_ambiguous_hits ≤ 1 or |secondary| = 0:
    U = 0                                    (no measurable competition)
elif uncertainty_mode == "second":           (default)
    U = 100 - ANI(rep(G_top), G_2)
elif uncertainty_mode == "max":
    U = 100 - max(ANI(rep(G_top), G_i)) for G_i in secondary
```

Where:
- `G_top` = genome of the best hit (Section 2.1)
- `G_2` = the first secondary representative under the same total order used for best-hit selection (weighted bitscore descending, then `pident` descending, then accession ascending). This tie-breaking was made explicit in PR #13; previously ties were left to an unstable sort, so `U` could differ between streaming and in-memory runs.
- `ANI(·,·)` = Average Nucleotide Identity from the pre-computed matrix. Pairs present in the matrix but not computed by the ANI tool (stored as 0) are looked up as `ANI_DEFAULT_UNRELATED` = 70% (so U = 30%). A secondary genome absent from the matrix yields no ANI value and U = 100%.
- The ANI matrix must be on a 0-100 scale (a matrix whose largest off-diagonal value is ≤ 1.0 is rejected), is stored as float64, and pairs supplied in only one direction are mirrored.

**Related per-read quantities** (used by the Bayesian Stage 2 refinement and the legacy cascade):

```
identity_gap      = pident_top - max(pident of any filtered hit to a different representative)
                    (null when no such hit exists; not restricted to the 95% window)
genus_uncertainty = 0                                  if |secondary| = 0
                    100 - min(ANI(rep(G_top), G_i))    otherwise (100 if no ANI value)
ambiguity_scope   = "unambiguous"     if |secondary| = 0
                    "within_species"  if all secondary ANI ≥ same_species_ani_threshold (95)
                    "within_genus"    if all secondary ANI ≥ same_genus_ani_threshold (80)
                    "across_genera"   otherwise
```

In the CLI classifier, `genus_uncertainty` is computed from the same 95% window as `U`; the `genus_bitscore_threshold_pct` (90%) window is used only by the programmatic `ANIWeightedClassifier`.

**Known issue (ambiguity scope):** the same-species and same-genus counts are taken over secondary *hits*, while the comparison is against the number of distinct secondary *representatives*. When one secondary representative contributes more than one hit inside the window (e.g. two contigs or two HSPs), the equality test fails and the read is labelled `across_genera` even if the competing genome shares > 95% ANI with the best hit. This can turn an `Ambiguous` MAP call into `Conserved Region` (Section 10.3). The behaviour is described here as it currently stands; it has not been changed.

**Interpretation:**

| Uncertainty | ANI Between Competitors | Biological Interpretation |
|-------------|-------------------------|---------------------------|
| U < 1.5% | ANI > 98.5% | Same species (confident placement) |
| 1.5% ≤ U < 5% | 95% < ANI ≤ 98.5% | Species boundary zone (ambiguous) |
| U ≥ 5% | ANI ≤ 95% | Different species (conserved gene region) |

**Biological rationale:** When competing genomes share high ANI (>98%), they represent the same species and placement ambiguity reflects strain-level variation. When competing genomes share low ANI (<95%), the read likely originates from a conserved gene present across multiple species.

### 2.3 Classification Decision Tree

`score classify` reports two calls per read:

- **`taxonomic_call`** (primary): the maximum a posteriori (MAP) category of the six-category Bayesian model, followed by the Stage 2 refinement (Section 10). This is the column used by the report, the novel-diversity clustering and `is_novel`.
- **`legacy_call`**: the hard-threshold cascade below, kept for side-by-side comparison (`core/classification/thresholds.py`). `score sensitivity` also re-classifies with this cascade.

The thresholds below are applied exactly only in `legacy_call`. In `taxonomic_call` they set the Gaussian means and widths, and the resulting decision boundaries are not the same (Section 2.5).

**Legacy cascade** (rules evaluated in order; first match wins; nucleotide defaults):

```
Rule 0   identity_gap < 2% AND num_ambiguous_hits > 1          -> Ambiguous
Rule 1   U ≥ 5% AND ambiguity_scope == "across_genera"         -> Conserved Region
Rule 2   U ≥ 5%                                                -> Ambiguous
Rule 3   U ≥ 1.5%                                              -> Species Boundary
Rule 4   N < 4% AND U < 1.5%                                   -> Known Species
Rule 5   num_ambiguous_hits ≤ 1 AND N ≥ 4% AND
         inferred_U(N) ≥ single_hit_uncertainty_threshold (10%) -> Ambiguous
Rule 6   4% ≤ N < 20% AND U < 1.5%                             -> Novel Species
Rule 7   20% ≤ N ≤ 25% AND U < 1.5% AND
         (genus_uncertainty < 10% OR no secondary genomes)      -> Novel Genus
else                                                           -> Unclassified
```

With the default inferred-uncertainty function (Section 3.2), Rule 5 makes every single-hit read with N ≥ 7% `Ambiguous` in `legacy_call`. A Novel Genus candidate that fails the genus-uncertainty condition of Rule 7 falls through to `Unclassified` in the legacy cascade (the CLI cascade has no `Ambiguous Within Genus` rule; that category comes only from the Bayesian Stage 2).

**Classification categories:**

| Category | Legacy criteria (`legacy_call`) | Biological Interpretation |
|----------|----------|---------------------------|
| Known Species | N < 4% AND U < 1.5% | Matches characterized species |
| Novel Species | 4% ≤ N < 20% AND U < 1.5% | Significant divergence, confident placement |
| Novel Genus | 20% ≤ N ≤ 25% AND U < 1.5% | Genus-level divergence, confident placement |
| Species Boundary | 1.5% ≤ U < 5% | Matches multiple closely related species |
| Ambiguous | U ≥ 5% (not across genera), identity gap < 2% with multiple hits, or single hit with N ≥ 7% | Placement cannot be resolved |
| Conserved Region | U ≥ 5% AND hits span genera | Conserved gene shared across genera |
| Ambiguous Within Genus | (Bayesian Stage 2 only) | Novel Genus candidate with high genus-level uncertainty |
| Unclassified | N > 25% or no rule matched | Very high divergence, manual review needed |
| Off-target | Family validation (`--target-family`) | Better hit outside the target family |

### 2.4 Threshold Justification

**Species boundary (96% ANI default):**
- Jain et al. (2018) analyzed 90,000+ prokaryotic genomes and established 95-96% ANI as the species boundary, corresponding to 70% DNA-DNA hybridization
- Goris et al. (2007) demonstrated ANI >95% correlates with >70% DDH
- Konstantinidis & Tiedje (2005) proposed ANI as the gold standard for species delineation
- metadarkmatter defaults to 96% ANI (N < 4%) for a slightly more conservative boundary

**Confident placement (>98.5% ANI):**
- When competing genomes share >98.5% ANI, they represent the same species
- Placement ambiguity at this level reflects strain variation, not taxonomic uncertainty
- Tightened from the previous 98% threshold to reduce false-confident assignments

**Adaptive thresholds:** When the `--adaptive-thresholds` option is enabled, the species boundary is detected from the ANI matrix distribution using a Gaussian Mixture Model rather than assuming the default 96% for all taxa. See [Section 9](#9-adaptive-thresholds-from-ani-distribution).

**Genus boundary (~75-80% ANI):**
- Qin et al. (2014) found genus boundaries vary from 65-83% ANI
- We use 75% (N=25%) as a conservative genus boundary
- For nucleotide mode, 80% identity (N=20%) marks the minimum for novel species

### 2.5 Hard Thresholds versus Bayesian Decision Boundaries

The Gaussian parameters of the Bayesian model are derived from the thresholds in Section 2.3 (Section 10.1), but the MAP decision boundaries that result from overlapping Gaussians do not coincide with those thresholds. The table lists the boundaries of the default model, obtained by scanning N (or U) in steps of 0.1 with a multi-hit read without prior modulation (`num_ambiguous_hits = 2`, `identity_gap = 10`):

| Transition (nucleotide mode) | `legacy_call` threshold | `taxonomic_call` boundary (approx.) |
|---|---|---|
| Known Species -> Novel Species (U = 0) | N = 4.0% | N ≈ 6.2% |
| Novel Species -> Novel Genus (U = 0) | N = 20.0% | N ≈ 18.6% |
| Novel Genus -> Unclassified (U = 0) | N > 25.0% | N ≈ 27.6% |
| Known Species -> Species Boundary (N = 1%) | U = 1.5% | U ≈ 3.8% |
| Species Boundary -> Ambiguous (N = 1%) | U = 5.0% | U ≈ 5.5% |
| Novel Species -> Species Boundary (N = 10%) | U = 1.5% | U ≈ 2.1% |
| Species Boundary -> Ambiguous (N = 10%) | U = 5.0% | U ≈ 6.4% |

| Transition (protein mode) | `legacy_call` threshold | `taxonomic_call` boundary (approx.) |
|---|---|---|
| Known Species -> Novel Species (U = 0) | N = 10.0% | N ≈ 11.3% |
| Novel Species -> Novel Genus (U = 0) | N = 25.0% | N ≈ 25.1% |
| Novel Genus -> Unclassified (U = 0) | N > 40.0% | N ≈ 40.4% |
| Along U at N = 1% | Species Boundary at U = 5%, Ambiguous at U = 10% | Ambiguous from U ≈ 8.6%; Species Boundary is not the MAP at this N |

At U = 0 in nucleotide mode the MAP becomes `Ambiguous` again only at N ≈ 48%, far beyond the Unclassified Gaussian.

**Single-hit reads.** For reads with `num_ambiguous_hits ≤ 1` (or whose hits all collapse to one representative), `placement_uncertainty` is 0 and this value is what the Bayesian model sees. The inferred uncertainty (Section 3.2) is reported in the `inferred_uncertainty` column and is used by legacy Rule 5, but it is not an input to `taxonomic_call`; changing `single_hit_uncertainty_threshold` therefore has no effect on `taxonomic_call`. The single-hit prior boost (1.5x on Ambiguous) does not move the boundaries above: with the default parameters a single-hit read is called Known Species for N < 6.2%, Novel Species for 6.2% ≤ N < 18.6% and Novel Genus for 18.6% ≤ N < 27.6%, whereas `legacy_call` calls every single-hit read with N ≥ 7% `Ambiguous`. Similarly, the 2.0x identity-gap boost does not change the MAP call at low U, while legacy Rule 0 always returns `Ambiguous`.

**Consequences.** Category counts from `taxonomic_call` and `legacy_call` can differ substantially, in particular for single-hit reads and for reads with small identity gaps. The threshold values quoted in this document and in `CLAUDE.md` describe `legacy_call` exactly and `taxonomic_call` only approximately. This is a description of current behaviour and a known limitation of the hand-tuned parameterisation; the defaults have not been changed (see `docs/CALIBRATION_RESULTS.md` for why fitted parameters are not shipped).

The boundaries can be reproduced by evaluating `BayesianClassifier(ScoringConfig()).classify_primary()` and `apply_legacy_thresholds()` on a grid of (N, U) values.

---

## 3. Enhanced Scoring Metrics

### 3.1 Confidence Score

The `confidence_score` column (0-100) is derived from the posterior entropy of the Bayesian model (Section 10.4):

```
confidence_score = clip((1 - H / log2(6)) × 100, 0, 100)      # rounded to 1 decimal
low_confidence   = H > bayesian.entropy_threshold              # default 1.5, i.e. confidence_score < ~42
```

When `bayesian.entropy_calibration` is set in a YAML config, the linear map is replaced by piecewise-linear interpolation between the supplied `[entropy, confidence]` knots. Off-target reads receive a fixed `confidence_score = 10.0` and `low_confidence = True`. `ScoringConfig.confidence_threshold` (default 50; 60 in the `literature-strict` preset) is not used by the CLI classifier.

**Multi-factor score (not used for output).** `core/constants.py` still contains `calculate_confidence_score()`, the earlier multi-factor score described below. No classifier calls it; it is documented here because the protein-mode scaling constants in Section 6 refer to it.

```
confidence_score = margin_score + placement_score + alignment_score
```

**Component 1: Margin from threshold boundaries (0-40 points)**

```python
if taxonomic_call == "Known Species":
    margin = novelty_known_max - novelty_index  # Distance from 4% boundary
    margin_score = min(40, (margin / 5.0) × 40)

elif taxonomic_call == "Novel Species":
    margin_lower = novelty_index - 4.0   # Distance from lower boundary (novelty_known_max)
    margin_upper = 20.0 - novelty_index  # Distance from upper boundary
    margin_score = min(40, (min(margin_lower, margin_upper) / 7.5) × 40)

elif taxonomic_call == "Novel Genus":
    margin_lower = novelty_index - 20.0
    margin_upper = 25.0 - novelty_index
    margin_score = min(40, (min(margin_lower, margin_upper) / 5.0) × 40)

else:
    margin_score = 10  # Base score for ambiguous/unclassified
```

**Uncertainty margin bonus (up to +10 points):**

```python
if placement_uncertainty < 1.5:              # uncertainty_confident_max
    uncertainty_margin = 1.5 - placement_uncertainty
    margin_score = min(40, margin_score + (uncertainty_margin / 1.5) × 10)
```

**Component 2: Placement certainty (0-40 points)**

Hit count factor (0-20 points):
| Ambiguous Hits | Points |
|----------------|--------|
| 1 | 20 |
| 2-3 | 15 |
| 4-5 | 10 |
| 6-10 | 5 |
| >10 | 0 |

Identity gap factor (0-20 points):
| Identity Gap | Points |
|--------------|--------|
| ≥ 5% | 20 |
| ≥ 2% | 15 |
| ≥ 1% | 10 |
| ≥ 0.5% | 5 |
| < 0.5% | 0 |
| no secondary hit, ≤ 1 ambiguous hit | 20 |

**Component 3: Alignment quality proxy (0-20 points)**

```python
alignment_score = min(20, max(0, (top_hit_identity - 70) / 25 × 20))
```

The score-range interpretation previously given for the multi-factor score does not apply to the entropy-based `confidence_score`; for the latter, the `low_confidence` flag (entropy > 1.5) is the only built-in cut-off.

### 3.2 Inferred Uncertainty for Single-Hit Reads

For reads with only one hit in the bitscore window (`num_ambiguous_hits ≤ 1`), or whose hits all map to the same species representative, placement uncertainty cannot be measured from ANI between competing genomes (`placement_uncertainty` is 0 for these reads). Instead, an uncertainty is inferred from the novelty level and reported in the `inferred_uncertainty` column.

**Rationale:** Higher novelty suggests either novel diversity OR database incompleteness, both of which increase placement uncertainty.

**Formula:**

```python
if novelty_index < 4%:
    # High identity: database likely complete for this species
    inferred_U = 5.0 + novelty_index × 0.5  # Range: 5.0-7.0%

elif novelty_index < 20%:
    # Novel species range: uncertain if truly novel or database gap
    inferred_U = 7.0 + (novelty_index - 4.0) × 1.0  # Range: 7.0-23.0%

elif novelty_index < 25%:
    # Novel genus range: high uncertainty
    inferred_U = 23.0 + (novelty_index - 20.0) × 1.5  # Range: 23.0-30.5%

else:
    # Very high divergence: maximum uncertainty
    inferred_U = 35.0
```

The function is continuous at N = 4% and N = 20% but steps from 30.5% to 35.0% at N = 25%. The breakpoints follow the effective novelty thresholds; in protein mode (10/25/40%) the same slopes give 10.0% at N = 10%, 25.0% at N = 25%, up to 47.5% just below N = 40%, and then 35.0%, so the protein-mode function is not monotone. (The programmatic `ANIWeightedClassifier` uses `calculate_inferred_uncertainty()` with the nucleotide breakpoints in both modes.)

**Use:** `inferred_uncertainty` feeds legacy Rule 5 (Section 2.3) and the legacy `placement_confidence` sub-score. It is not an input to the Bayesian `taxonomic_call` or to `confidence_score` (Section 2.5).

**Output field:** `uncertainty_type` = "inferred" or "measured"

### 3.3 Alignment Quality Score

Sections 3.3-3.5 describe the legacy sub-scores. In `score classify` they are written only with `--include-legacy-scores`; the programmatic `ANIWeightedClassifier` always computes them. The two implementations differ where noted.

**CLI (`VectorizedClassifier._compute_legacy_scores`):**

```
alignment_quality = clip(alignment_fraction × 50 + 25, 0, 100)
alignment_fraction = (qend - qstart + 1) / read_length   of the best hit
```

**Programmatic API (`calculate_alignment_quality` in `core/constants.py`):** incorporates BLAST statistics beyond percent identity.

**Formula:**

```
alignment_quality = mismatch_score + gap_score + coverage_score + evalue_score
```

**Components (each 0-25 points):**

```python
# Mismatch density penalty
mismatch_rate = mismatch / alignment_length
mismatch_score = max(0, 25 - mismatch_rate × 50)

# Gap complexity penalty
gap_rate = gapopen / alignment_length × 100
gap_score = max(0, 25 - gap_rate × 25)

# Coverage bonus
coverage = (qend - qstart + 1) / read_length
coverage_score = min(25, coverage × 25)

# E-value significance
if evalue ≤ 1e-50: evalue_score = 25
elif evalue ≤ 1e-20: evalue_score = 20
elif evalue ≤ 1e-10: evalue_score = 15
elif evalue ≤ 1e-5: evalue_score = 10
else: evalue_score = 5
```

### 3.4 Identity and Placement Confidence

These orthogonal metrics separate confidence in what a sequence is from confidence in where it belongs.

**Identity Confidence** - How reliable is the identity measurement?

```python
# Base score from novelty (inverted)
if novelty_index < 5%: base = 80
elif novelty_index < 15%: base = 60 - (novelty_index - 5)
elif novelty_index < 25%: base = 50 - (novelty_index - 15) × 1.5
else: base = 30

# Alignment quality contribution
identity_confidence = min(100, base + alignment_quality × 0.2)
```

**Placement Confidence** - How confident is the genome assignment?

```python
# Base score from uncertainty level
if uncertainty < 2%: base = 80
elif uncertainty < 5%: base = 60
elif uncertainty < 10%: base = 40
elif uncertainty < 20%: base = 25
else: base = 10

# Penalty for inferred (unmeasured) uncertainty
if uncertainty_type == "inferred":
    base -= 15

# Identity gap bonus (for multi-hit reads)
if identity_gap is not None and num_ambiguous_hits > 1:
    if identity_gap ≥ 5%: gap_bonus = 20
    elif identity_gap ≥ 2%: gap_bonus = 10
    elif identity_gap ≥ 1%: gap_bonus = 5
    else: gap_bonus = 0

placement_confidence = max(0, min(100, base + gap_bonus))
```

In the CLI, reads with `num_ambiguous_hits ≤ 1` use `inferred_uncertainty` with the 15-point penalty already applied (65/45/25/10/0) and no gap bonus; all other reads use the measured `placement_uncertainty` and the gap bonus.

### 3.5 Discovery Score

The discovery score prioritizes novel findings for experimental validation. Only calculated for reads whose `taxonomic_call` is Novel Species or Novel Genus.

**Formula:**

```
discovery_score = novelty_pts + quality_pts + confidence_pts
```

**Components:**

```python
# Novelty component (0-40 points)
if taxonomic_call == "Novel Species":
    # novelty_novel_species_min = 4 (nucleotide)
    novelty_pts = 15 + (novelty_index - 4) × 1.5

elif taxonomic_call == "Novel Genus":
    # novelty_novel_genus_min = 20 (nucleotide)
    novelty_pts = 30 + (novelty_index - 20) × 2.0

# CLI: Novel Species points clipped to [15, 30], Novel Genus to [30, 40]
# Programmatic API: novelty_pts = min(40, novelty_pts), no lower bound

# Quality component (0-30 points)
quality_pts = alignment_quality × 0.3

# Confidence component (0-30 points)
confidence_pts = (identity_confidence + placement_confidence) / 2 × 0.3
```

**Score interpretation:**

| Discovery Score | Priority | Recommended Action |
|-----------------|----------|-------------------|
| 75-100 | High | Prioritize for experimental validation |
| 50-74 | Moderate | Include in candidate list |
| 25-49 | Low | Needs additional evidence |
| < 25 | Unreliable | Likely artifact |

---

## 4. Novel Diversity Clustering

### 4.1 Clustering Algorithm

Reads classified as Novel Species or Novel Genus are grouped into clusters representing putative novel organisms.

**Clustering key:**

```
cluster_key = (representative_genome, novelty_band, taxonomic_call)
```

Where:
- `representative_genome` = neighborhood representative of `best_match_genome`. When an ANI matrix is available, genomes whose pairwise ANI is ≥ 80% (`genome_neighborhood_threshold`) are joined by union-find; the representative of each connected component is the genome with the most novel reads (ties broken by accession). Without an ANI matrix each genome is its own representative.
- `novelty_band` = floor(novelty_index / 5) × 5 (5% intervals: 5-10%, 10-15%, etc.)
- `taxonomic_call` = "Novel Species" or "Novel Genus" (the Bayesian `taxonomic_call`)

**Algorithm steps** (`core/novel_diversity/clustering.py`):

1. Filter to novel reads only (`taxonomic_call ∈ {"Novel Species", "Novel Genus"}`)
2. Compute novelty band and effective uncertainty (Section 4.3) for each read
3. Map genomes to neighborhood representatives
4. Group by clustering key
5. Merge adjacent bands with the same (representative, call) when `novelty_max` of the lower group ≥ `novelty_min` of the next group minus `merge_band_tolerance` (default 2.0); chained merges are allowed
6. Filter clusters with `read_count ≥ min_cluster_size` (default: 3); merging precedes filtering so that sub-threshold groups in neighbouring bands are not discarded
7. Sort by read count (descending; ties by representative, then call) and generate cluster IDs: `NSP_001`, `NSP_002`... for Novel Species; `NGN_001`... for Novel Genus
8. Assign confidence ratings and compute cluster statistics

**Cluster statistics:**

| Metric | Calculation |
|--------|-------------|
| read_count | Number of reads in cluster |
| mean_novelty_index | Mean of read novelty values |
| novelty_min, novelty_max | Range of novelty in cluster |
| mean_placement_uncertainty | Mean uncertainty (using effective uncertainty, Section 4.3) |
| mean_discovery_score | Mean discovery score (only if the input has a `discovery_score` column, i.e. `--include-legacy-scores`) |

### 4.2 Cluster Confidence Ratings

Clusters are assigned confidence ratings based on multiple criteria.

**Rating criteria when `confidence_score` is present** (always the case for `score classify` output):

| Rating | Criteria |
|--------|----------|
| High | mean confidence_score ≥ 70 AND read_count ≥ 10 |
| Medium | mean confidence_score ≥ 50 AND read_count ≥ 5 |
| Low | All other clusters meeting minimum size |

**Legacy rating criteria** (used only when the input has no `confidence_score` column):

| Rating | Criteria | Interpretation |
|--------|----------|----------------|
| High | read_count ≥ 10 AND mean_uncertainty < 5% AND mean_discovery ≥ 75 | Prioritize for experimental validation |
| Medium | read_count ≥ 5 AND mean_uncertainty < 10% AND mean_discovery ≥ 50 | Include in candidate list |
| Low | All other clusters meeting minimum size | May need additional evidence |

**Simplified criteria (without discovery score):**

| Rating | Criteria |
|--------|----------|
| High | read_count ≥ 10 AND mean_uncertainty < 5% |
| Medium | read_count ≥ 5 AND mean_uncertainty < 10% |
| Low | All other clusters |

### 4.3 Estimated ANI Calculation

For novel clusters, the estimated ANI to the nearest reference is calculated as:

```
estimated_ANI = 100 - mean_novelty_index
```

This provides an approximate measure of genomic similarity, though read-level identity may differ from genome-level ANI by 10-20%.

**Effective uncertainty for single-hit reads:**

The effective uncertainty used for cluster statistics (`mean_placement_uncertainty`) is:

```python
if num_ambiguous_hits <= 1:      # when the column is absent: placement_uncertainty < 0.1
    effective_uncertainty = min(novelty_index × 0.4, 15.0)
else:
    effective_uncertainty = placement_uncertainty
```

This is a separate approximation from the `inferred_uncertainty` column of Section 3.2.

---

## 5. Coverage Weighting

Coverage weighting adjusts hit selection to prioritize longer alignments over short conserved domains. This is enabled by default (`coverage_weight_mode="linear"`) to reduce conserved-gene bias in classification.

**Weighted score calculation:**

```python
coverage = clip((qend - qstart + 1) / read_length, 0, 1)   # after reverse-strand normalisation
weight = calculate_weight(coverage, mode, strength)
weighted_score = bitscore × weight
```

`read_length` is `qlen` when the alignment file has a 13th `qlen` column, otherwise `qend` (an approximation that depends on where the alignment lies on the read). The weighted score is used both to select the best hit (Section 2.1) and to define the 95% window of competing hits (Section 2.2). Before PR #13, minus-strand hits with `qstart > qend` received a negative coverage and were down-weighted or filtered out; coordinates are now normalised first.

**Weight calculation by mode:**

```python
min_weight = 1.0 - strength  # Default strength = 0.5 → min = 0.5
max_weight = 1.0 + strength  # Default strength = 0.5 → max = 1.5
weight_range = max_weight - min_weight

if mode == "linear":
    normalized = coverage

elif mode == "log":
    normalized = log(1 + 9 × coverage) / log(10)

elif mode == "sigmoid":
    normalized = 1 / (1 + exp(-10 × (coverage - 0.6)))

weight = min_weight + weight_range × normalized
```

**Mode characteristics:**

| Mode | Behavior | Use Case |
|------|----------|----------|
| none | Raw bitscore | Disables coverage weighting |
| linear | Proportional to coverage | Default; balanced weighting |
| log | Diminishing returns | Gentle preference for coverage |
| sigmoid | Step function at 60% | Strict coverage requirement |

---

## 6. Protein Mode Thresholds

When using protein-level alignment (BLASTX + AAI), wider thresholds are applied to account for slower protein evolution.

**Protein mode thresholds:**

| Parameter | Nucleotide Mode | Protein Mode |
|-----------|-----------------|--------------|
| novelty_known_max | 4% | 10% |
| novelty_novel_species_min | 4% | 10% |
| novelty_novel_species_max | 20% | 25% |
| novelty_novel_genus_min | 20% | 25% |
| novelty_novel_genus_max | 25% | 40% |
| uncertainty_known_max | 1.5% | 5% |
| uncertainty_novel_species_max | 1.5% | 5% |
| uncertainty_novel_genus_max | 1.5% | 5% |
| uncertainty_conserved_min | 5% | 10% |

In protein mode `ScoringConfig.get_effective_thresholds()` returns the fixed constants of `core/protein_constants.py`. Custom novelty/uncertainty thresholds (CLI flags, YAML `--config`, presets) and adaptive thresholds are therefore ignored in protein mode, both for `legacy_call` and for the Bayesian Gaussian parameters. Placement uncertainty is computed from the ANI matrix in both modes.

**AAI genus boundaries** (Riesco & Trujillo, 2024):
- AAI > 65%: Same genus
- AAI 58-65%: Genus boundary zone
- AAI < 58%: Different genus

These AAI rules are applied after the Bayesian call and Stage 2 refinement, to reads with N in the novel-genus range and U below the confident-placement maximum, by both `score classify --aai` (`VectorizedClassifier`) and the programmatic `ANIWeightedClassifier`. They override both `taxonomic_call` and `legacy_call`. The AAI value used is the maximum AAI between the top genome and any competing hit genome with AAI > 0; reads with no such value keep their ANI-based call.

**Multi-factor confidence score scaling** (constants used only by the unused `calculate_confidence_score()`, Section 3.1):

| Parameter | Nucleotide | Protein |
|-----------|------------|---------|
| margin_divisor_known | 5.0 | 10.0 |
| margin_divisor_novel_species | 7.5 | 15.0 |
| margin_divisor_novel_genus | 5.0 | 15.0 |
| identity_gap_thresholds | (5, 2, 1, 0.5) | (10, 5, 2, 1) |
| identity_score_base | 70.0 | 50.0 |
| identity_score_range | 25.0 | 40.0 |

---

## 7. Quality Control Metrics

metadarkmatter computes QC metrics at two stages: before classification (to identify problematic inputs) and after classification (to flag potentially unreliable results).

### 7.1 Pre-Classification QC

Pre-classification metrics are computed after alignment loading and quality filtering:

| Metric | Calculation | Warning Threshold |
|--------|-------------|-------------------|
| Filter rate | `removed_alignments / total_alignments` | > 50% |
| Genome coverage | `|genomes_in_alignments ∩ genomes_in_ani| / genomes_in_ani` | < 50% (only when > 3 ANI genomes) |
| Single-hit fraction | `reads_with_one_alignment / total_reads` | > 80% |
| Mean alignment length | Mean of alignment lengths (after filtering) | (informational) |
| Mean identity | Mean of percent identities (after filtering) | < 80% |

**Filter rate** indicates how many alignments are removed before classification: by the quality thresholds (e.g., minimum alignment fraction) and, when family validation is active, by removal of out-of-family hits. A high filter rate suggests either stringent filters or low-quality alignments.

**Genome coverage** measures the fraction of ANI matrix genomes that appear in the filtered alignments. Low coverage suggests that many reference genomes have no matching reads, which may affect classification accuracy.

**Single-hit fraction** counts reads with exactly one alignment row after filtering. For these reads, placement uncertainty cannot be directly measured from competing hits (see Section 3.2).

### 7.2 Post-Classification QC

Post-classification metrics assess the overall classification quality:

| Metric | Calculation | Warning Threshold |
|--------|-------------|-------------------|
| Ambiguous fraction | `reads with taxonomic_call == "Ambiguous" / total_reads` | > 50% |
| Low confidence fraction | `reads with low_confidence == True / total` | > 30% |
| Novel fraction | `(novel_species + novel_genus) / total` | (informational) |

The ambiguous fraction counts only the `Ambiguous` category (not Conserved Region, Species Boundary or Ambiguous Within Genus). `low_confidence` is `posterior_entropy > 1.5` (Off-target rows are always flagged); if the column is absent, `confidence_score < 50` on the 0-100 scale is used instead. Before PR #13 this metric compared the 0-100 score with 0.5.

QC metrics are saved as JSON when `--qc-output` is specified.

---

## 8. Threshold Sensitivity Analysis

Sensitivity analysis evaluates how classification counts change across a range of threshold values. This helps assess whether results are robust to the specific threshold choice or whether small threshold changes lead to large shifts in classification counts.

**Method:**

The analysis sweeps the `novelty_known_max` threshold (and `novelty_novel_species_min` with it) across a user-specified range (default: 2-8%) while `uncertainty_known_max`, `uncertainty_novel_species_max` and `uncertainty_novel_genus_max` are swept jointly over a second range (default: 0.5-3.0%). At each step, all reads are re-classified from the pre-computed metrics in a saved `score classify` output (no re-computation of alignment metrics). All other settings are taken from the base configuration (`--config`, or the defaults).

```
for novelty_thresh, uncertainty_thresh in zip(linspace(2.0, 8.0, 9), linspace(0.5, 3.0, 9)):
    counts[category] = legacy_cascade(reads, novelty_thresh, uncertainty_thresh)
```

The re-classification uses the hard-threshold cascade of Section 2.3 (`apply_classification_thresholds`), not the Bayesian model, so the counts describe the sensitivity of `legacy_call`. Off-target rows are kept unchanged.

**Output:** A TSV (one row per threshold point and category) or JSON with category counts as a function of threshold position. Stable results (flat lines) indicate robustness; large count changes near the default threshold indicate sensitivity.

**Interpretation:**
- Flat count profiles across thresholds suggest classification is robust
- Sharp transitions at or near the default threshold indicate threshold-dependent results
- Gradual transitions are expected near true biological boundaries

---

## 9. Adaptive Thresholds from ANI Distribution

The default species boundary (96% ANI) is appropriate for many bacterial lineages but may not fit all taxa equally. Rapidly evolving lineages may have narrower within-species ANI distributions, while conserved lineages may show broader distributions.

### 9.1 Gaussian Mixture Model

When `--adaptive-thresholds` is enabled, the species boundary is estimated from the actual ANI matrix using a 2-component Gaussian Mixture Model:

**Step 1:** Extract upper-triangle pairwise ANI values from the ANI matrix, excluding self-comparisons and pairs that were not computed by the ANI tool (stored as 0). The default value of 70% that classification substitutes for uncomputed pairs is not included; before PR #13 these substituted values entered the fit and could collapse one component onto 70%.

**Step 2:** Fit a 2-component GMM:
```
GMM(n_components=2) -> (mu_low, sigma_low), (mu_high, sigma_high)
```

The two components typically correspond to:
- **Within-species comparisons** (high ANI, narrow distribution)
- **Between-species comparisons** (lower ANI, broader distribution)

**Step 3:** Compute the species boundary as the weighted midpoint:
```
w_low = 1 / sigma_low
w_high = 1 / sigma_high
boundary = (mu_low * w_low + mu_high * w_high) / (w_low + w_high)
```

The boundary is clamped to the range 80-99% ANI. The fit uses `n_init=5`, `max_iter=200` and `random_state = METADARKMATTER_SEED` (default 42).

**Step 4:** `novelty_known_max` and `novelty_novel_species_min` are set to `100 - boundary`; all other settings of the base configuration are kept. Because the Known Species and Novel Species Gaussians are derived from these thresholds (Section 10.1), the adaptive boundary also shifts the Bayesian model, not only `legacy_call`. In protein mode the adjusted thresholds are ignored (Section 6).

### 9.2 Separation Quality

The GMM separation quality is assessed using the standardized gap between component means:

```
separation = |mu_high - mu_low| / (sigma_low + sigma_high)
confidence = min(1, separation / 3)          # reported, not used further
```

| Separation | Action |
|------------|--------|
| < 1.0 | Falls back to default threshold |
| ≥ 1.0 | Uses the GMM boundary |

A separate 3-component GMM (`detect_genus_boundary`, separation cut-off 0.8, boundary clamped to 70-90% ANI) is used only by the phylogenetic-neighborhood analysis in the report's novel-diversity tab; it does not change read classification.

### 9.3 Fallback Conditions

The adaptive method falls back to the default threshold (96% ANI) when:
- Fewer than 5 genomes in the ANI matrix
- Fewer than 3 computed (non-zero) pairwise ANI values
- GMM separation quality < 1.0
- GMM fails to converge
- scikit-learn is not installed

**Dependency:** Requires `scikit-learn>=1.3.0` (install with `pip install metadarkmatter[adaptive]`).

---

## 10. Bayesian Classification Framework

The Bayesian model produces the primary `taxonomic_call`. It assigns posterior probabilities P(category | N, U) over six categories and takes the maximum a posteriori (MAP) category, which is then refined in a discrete second stage. Near classification boundaries the posteriors are split across categories, making the classification uncertainty explicit. Implementation: `core/classification/bayesian.py` (`BayesianClassifier.classify_primary`); configuration: `BayesianConfig` in `models/config.py`.

### 10.1 Likelihood Model

Each of the six Stage 1 categories (Known Species, Novel Species, Novel Genus, Species Boundary, Ambiguous, Unclassified) is modeled as a 2D Gaussian with diagonal covariance, centered on its expected (novelty, uncertainty) values. The density includes its normalisation constant, so that broad categories are not favoured over narrow ones:

```
L(category | N, U) = 1 / (2π σ_N σ_U) × exp(-0.5 × (z_N² + z_U²))

where:
    z_N = (N - μ_N) / σ_N
    z_U = (U - μ_U) / σ_U
```

**Derivation from the thresholds** (`build_category_params`; `eff` = effective thresholds of Section 2.3 / Section 6):

| Category | μ_N | σ_N | μ_U | σ_U |
|----------|-----|-----|-----|-----|
| Known Species | known_max / 2 | known_max / 2 + 0.5 | u_known_max / 2 | u_known_max / 2 + 0.5 |
| Novel Species | (sp_min + sp_max) / 2 | (sp_max - sp_min) / 3 + 0.5 | u_sp_max / 2 | u_sp_max / 2 + 0.5 |
| Novel Genus | (gen_min + gen_max) / 2 | (gen_max - gen_min) / 3 + 0.5 | u_gen_max / 2 | u_gen_max / 2 + 0.5 |
| Species Boundary | (sp_min + sp_max) / 2 | (sp_max - sp_min) / 3 + 0.5 | (u_gen_max + u_conserved) / 2 | (u_conserved - u_gen_max) / 3 + 0.5 |
| Ambiguous | 15.0 (fixed) | 10.0 | u_conserved + 2 | 5.0 |
| Unclassified | gen_max + 5 | 5.0 | 5.0 | 5.0 |

**Resulting default parameters:**

| Category | μ_N | σ_N | μ_U | σ_U | protein μ_N | σ_N | μ_U | σ_U |
|----------|-----|-----|-----|-----|------|------|------|------|
| Known Species | 2.0 | 2.5 | 0.75 | 1.25 | 5.0 | 5.5 | 2.5 | 3.0 |
| Novel Species | 12.0 | 5.83 | 0.75 | 1.25 | 17.5 | 5.5 | 2.5 | 3.0 |
| Novel Genus | 22.5 | 2.17 | 0.75 | 1.25 | 32.5 | 5.5 | 2.5 | 3.0 |
| Species Boundary | 12.0 | 5.83 | 3.25 | 1.67 | 17.5 | 5.5 | 7.5 | 2.17 |
| Ambiguous | 15.0 | 10.0 | 7.0 | 5.0 | 15.0 | 10.0 | 12.0 | 5.0 |
| Unclassified | 30.0 | 5.0 | 5.0 | 5.0 | 45.0 | 5.0 | 5.0 | 5.0 |

When `bayesian.category_params` is set in a YAML config (produced by `scripts/calibrate_bayesian.py`), the fitted values override these defaults for the categories they specify. Fitted parameters are not shipped as defaults (`docs/CALIBRATION_RESULTS.md`).

Non-finite N or U values are set to 0 with a warning before the likelihood is evaluated.

### 10.2 Posterior Computation

Posteriors are computed with Bayes' rule:

```
P(c | N, U) = L(c | N, U) × π_c / Σ_k L(k | N, U) × π_k
```

The priors π default to uniform (1/6 each) and must sum to 1 (± 0.01) when set in `bayesian.priors`.

**Prior modulation for Ambiguous.** Before normalisation, the Ambiguous prior of a read is multiplied by:
- `identity_gap_boost` (default **2.0**) if `identity_gap` is not null, `identity_gap < identity_gap_ambiguous_max` (default 2%) and `num_ambiguous_hits > 1`;
- `single_hit_boost` (default **1.5**) if `num_ambiguous_hits == 1` and `N ≥ novelty_novel_species_min` (effective threshold; 4% nucleotide, 10% protein).

The two conditions are mutually exclusive (`num_ambiguous_hits > 1` versus `== 1`). (The scalar `compute_posteriors` API compares against the configured `novelty_novel_species_min` rather than the effective threshold, which differs in protein mode.) The boosted priors are not renormalised separately; the posterior normalisation absorbs them. If the sum of the weighted likelihoods underflows (< 1e-300), the read receives uniform posteriors (1/6).

With the default Gaussians these boosts are modest relative to the likelihood differences: they do not change the MAP call in the low-uncertainty region (Section 2.5).

**Output columns:** `p_known_species`, `p_novel_species`, `p_novel_genus`, `p_species_boundary`, `p_ambiguous`, `p_unclassified` (sum to 1.0; rounded to 4 decimals), and `posterior_entropy`.

### 10.3 Stage 2 Refinement

The Stage 1 MAP category is refined with phylogenetic context (`apply_stage2_refinement`):

| Stage 1 MAP | Condition | `taxonomic_call` |
|---|---|---|
| Ambiguous | `ambiguity_scope == "across_genera"` | Conserved Region |
| Novel Genus | `genus_uncertainty ≥ genus_uncertainty_ambiguous_min` (10%) AND at least one secondary genome | Ambiguous Within Genus |
| any other | - | unchanged |

`ambiguity_scope` and `genus_uncertainty` are defined in Section 2.2. Off-target reads (family validation) bypass the model and are reported with `taxonomic_call = legacy_call = "Off-target"`.

### 10.4 Shannon Entropy and Confidence

The posterior entropy quantifies classification confidence as a single scalar:

```
H = -Σ p_i × log2(p_i)   over categories with p_i > 0
confidence_score = clip((1 - H / log2(6)) × 100, 0, 100)
```

For six categories, entropy ranges from 0 (all probability mass on one category) to log2(6) ≈ 2.585 (uniform posterior). `low_confidence` is set when H > `bayesian.entropy_threshold` (default 1.5, equivalent to confidence_score < 42.0). An optional empirical mapping (`bayesian.entropy_calibration`, produced by `scripts/calibrate_entropy.py`) replaces the linear map by piecewise-linear interpolation, clipped to [0, 100].

| Entropy | Interpretation |
|---------|---------------|
| H < 0.5 | Confident assignment to one category |
| 0.5 ≤ H < 1.0 | Moderate confidence; one category dominant |
| 1.0 ≤ H ≤ 1.5 | Low confidence; multiple categories plausible |
| H > 1.5 | Flagged `low_confidence`; read lies near a classification boundary |

These bands are descriptive; the posteriors are not calibrated probabilities (`docs/STATISTICAL_ASSUMPTIONS.md`, sections 4 and 10).

### 10.5 Clade-Novelty Rule

Placement uncertainty U measures how uncertain the assignment to a **particular** reference is. In a clade of close relatives (several species of one genus at 85-95% ANI to one another), U is high even for a read that is clearly divergent from every one of them, and such reads fall into Ambiguous or Species Boundary instead of Novel Species. The read is then novel whichever reference it is closest to; only its placement within the clade is uncertain.

With `clade_novel_species = true` (the default; disable with `score classify --no-clade-novel-species` or `clade_novelty: {enabled: false}` in a YAML config), the following rule is applied after Stage 2 (nucleotide mode only; ignored with a warning in protein mode):

| Condition | Value |
|---|---|
| `taxonomic_call` | Ambiguous or Species Boundary |
| `ambiguity_scope` | `within_species` or `within_genus` (reads whose competitors span genera stay Conserved Region) |
| novelty to the best reference, N | `>= clade_novelty_min` (default 5.0) and `< clade_novelty_max` (default 18.6) |
| novelty to the best competing reference, `100 - second_hit_identity` | `>= clade_novelty_min` |

Selected reads get `taxonomic_call = "Novel Species"` (and therefore `diversity_status = Novel`, `is_novel = true`). The posteriors, entropy, `placement_uncertainty` and `ambiguity_scope` are unchanged, so the placement ambiguity remains visible, and an extra boolean column `clade_novel_species` records which reads the rule changed. The defaults lie between the legacy Known Species threshold (4%) and the measured Bayesian Known/Novel Species boundary (about 6.2%), and at the measured Bayesian Novel Species/Novel Genus boundary (18.6%), both at U = 0 (Section 2.5).

In a leave-one-out benchmark on three families (Francisellaceae, Lactobacillaceae, Enterobacteriaceae; `scripts/benchmark_leave_one_out.py`), the rule raised read-level novel-species recall from 0.47-0.52 to 0.59-0.70, left novel-versus-known precision at or above 0.989 and Novel Genus recall unchanged, and increased the fraction of genus-level reads called Novel Species. Details and limitations: [CLADE_NOVELTY_VALIDATION.md](CLADE_NOVELTY_VALIDATION.md). On the strength of that validation the rule is on by default; setting it off reproduces the previous classification exactly.

### 10.6 Interpretation

The MAP category and the legacy cascade often agree away from boundaries, but they are not equivalent (Section 2.5). Near boundaries the posteriors show uncertainty that the hard thresholds do not.

**Example** (default nucleotide configuration, multi-hit read, identity gap 5%): a read at (N = 3.8%, U = 1.4%) lies just inside the Known Species thresholds. The model gives p_known = 0.752, p_novel_species = 0.156, p_species_boundary = 0.072, p_ambiguous = 0.020, p_novel_genus ≈ 0, p_unclassified ≈ 0; MAP = Known Species, H = 1.11, confidence_score = 57.0. Both `taxonomic_call` and `legacy_call` are Known Species, but the entropy indicates a boundary case.

**Vectorized implementation:** `classify_primary` evaluates all reads with numpy broadcasting over an (n_reads × 6) array, without Python loops. The scalar `compute_posteriors` method (used by `ANIWeightedClassifier`) implements the same model.

---

## 11. References

1. **Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S.** (2018). High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. *Nature Communications* 9:5114. https://doi.org/10.1038/s41467-018-07641-9

2. **Goris J, Konstantinidis KT, Klappenbach JA, Coenye T, Vandamme P, Tiedje JM.** (2007). DNA-DNA hybridization values and their relationship to whole-genome sequence similarities. *International Journal of Systematic and Evolutionary Microbiology* 57:81-91. https://doi.org/10.1099/ijs.0.64483-0

3. **Konstantinidis KT, Tiedje JM.** (2005). Genomic insights that advance the species definition for prokaryotes. *Proceedings of the National Academy of Sciences* 102:2567-2572. https://doi.org/10.1073/pnas.0409727102

4. **Qin QL, Xie BB, Zhang XY, Chen XL, Zhou BC, Zhou J, Oren A, Zhang YZ.** (2014). A proposed genus boundary for the prokaryotes based on genomic insights. *Journal of Bacteriology* 196:2210-2215. https://doi.org/10.1128/JB.01688-14

5. **Parks DH, Chuvochina M, Waite DW, Rinke C, Skarshewski A, Chaumeil PA, Hugenholtz P.** (2018). A standardized bacterial taxonomy based on genome phylogeny substantially revises the tree of life. *Nature Biotechnology* 36:996-1004. https://doi.org/10.1038/nbt.4229

6. **Riesco R, Trujillo ME.** (2024). Update on the proposed minimal standards for the use of genome data for the taxonomy of prokaryotes. *International Journal of Systematic and Evolutionary Microbiology* 74:006300. https://doi.org/10.1099/ijsem.0.006300

7. **Rodriguez-R LM, Konstantinidis KT.** (2014). Bypassing cultivation to identify bacterial species. *Microbe Magazine* 9:111-118.

8. **Luo C, Rodriguez-R LM, Konstantinidis KT.** (2014). MyTaxa: an advanced taxonomic classifier for genomic and metagenomic sequences. *Nucleic Acids Research* 42:e73. https://doi.org/10.1093/nar/gku169

---

## Appendix: Constants Reference

Nucleotide threshold constants are defined in `src/metadarkmatter/core/constants.py` (protein constants in `core/protein_constants.py`). At run time the nucleotide values come from `ScoringConfig` (`models/config.py`), whose defaults equal these constants:

```python
# Novelty Index Thresholds
NOVELTY_KNOWN_MAX = 4.0          # 96% ANI species boundary (Jain et al. 2018)
NOVELTY_NOVEL_SPECIES_MIN = 4.0
NOVELTY_NOVEL_SPECIES_MAX = 20.0
NOVELTY_NOVEL_GENUS_MIN = 20.0
NOVELTY_NOVEL_GENUS_MAX = 25.0

# Placement Uncertainty Thresholds
UNCERTAINTY_CONFIDENT_MAX = 1.5   # Known / Novel Species / Novel Genus maximum
UNCERTAINTY_CONSERVED_MIN = 5.0

# ANI Boundaries
ANI_SPECIES_BOUNDARY_LOW = 95.0
ANI_SPECIES_BOUNDARY_HIGH = 96.0
ANI_GENUS_BOUNDARY_LOW = 75.0
ANI_GENUS_BOUNDARY_HIGH = 80.0
ANI_DEFAULT_UNRELATED = 70.0      # substituted for uncomputed pairs

# AAI Boundaries
AAI_GENUS_BOUNDARY_HIGH = 65.0
AAI_GENUS_BOUNDARY_LOW = 58.0
AAI_DEFAULT_UNRELATED = 50.0

# Inferred uncertainty (single-hit reads)
INFERRED_UNCERTAINTY_BASE = 5.0
INFERRED_UNCERTAINTY_NOVEL_SPECIES_SLOPE = 1.0
INFERRED_UNCERTAINTY_NOVEL_GENUS_SLOPE = 1.5
INFERRED_UNCERTAINTY_MAX = 35.0
```

Other `ScoringConfig` defaults referred to in this document: `bitscore_threshold_pct = 95`, `genus_bitscore_threshold_pct = 90` (programmatic API only), `uncertainty_mode = "second"`, `coverage_weight_mode = "linear"`, `coverage_weight_strength = 0.5`, `min_alignment_length = 100`, `min_alignment_fraction = 0.3`, `identity_gap_ambiguous_max = 2.0`, `genus_uncertainty_ambiguous_min = 10.0`, `single_hit_uncertainty_threshold = 10.0`, `same_species_ani_threshold = 95`, `same_genus_ani_threshold = 80`, `family_ratio_threshold = 0.8`, `merge_band_tolerance = 2.0`; `BayesianConfig`: uniform priors, `identity_gap_boost = 2.0`, `single_hit_boost = 1.5`, `entropy_threshold = 1.5`.

**Default configuration changes (v0.2):**

| Parameter | Previous Default | Current Default | Rationale |
|-----------|-----------------|-----------------|-----------|
| `novelty_known_max` | 5.0 | 4.0 | 96% ANI species boundary |
| `uncertainty_known_max` | 2.0 | 1.5 | Stricter confident placement |
| `coverage_weight_mode` | "none" | "linear" | Reduce conserved-gene bias |
| `min_alignment_fraction` | 0.0 | 0.3 | Filter low-coverage alignments |
| Enhanced scoring | opt-in flag | always on | Inferred uncertainty, quality metrics (legacy sub-scores written only with `--include-legacy-scores`) |
| Single-hit inference | opt-in flag | always on | `inferred_uncertainty` column; used by `legacy_call` only (Section 2.5) |
| Primary call | threshold cascade | Bayesian MAP + Stage 2 | Cascade retained as `legacy_call` |
