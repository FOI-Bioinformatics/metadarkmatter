# Computational Methods

This document describes the computational methods used by metadarkmatter for detecting novel bacterial diversity from metagenomic sequencing data. All calculations are presented with their mathematical formulations and biological interpretations.

**Version:** Corresponds to metadarkmatter v0.2
**Last updated:** 2026-02-10

---

## Table of Contents

1. [Overview](#1-overview)
2. [Read Classification Algorithm](#2-read-classification-algorithm)
   - 2.1 [Novelty Index](#21-novelty-index-n)
   - 2.2 [Placement Uncertainty](#22-placement-uncertainty-u)
   - 2.3 [Classification Decision Tree](#23-classification-decision-tree)
   - 2.4 [Threshold Justification](#24-threshold-justification)
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
10. [Bayesian Confidence Framework](#10-bayesian-confidence-framework)
    - 10.1 [Likelihood Model](#101-likelihood-model)
    - 10.2 [Posterior Computation](#102-posterior-computation)
    - 10.3 [Shannon Entropy](#103-shannon-entropy)
    - 10.4 [Interpretation](#104-interpretation)
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
- `pident_top` = percent identity of the best BLAST hit (highest bitscore)
- N is expressed as a percentage (0-100%)

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
threshold = 0.95 × bitscore_top
ambiguous_hits = {hit : hit.bitscore ≥ threshold}
competing_genomes = unique genomes from ambiguous_hits
```

**Step 2: Calculate uncertainty from ANI**

```
secondary_genomes = competing_genomes - {G_top}

if |secondary_genomes| = 0:
    U = 0  (no competition)
else:
    U = 100 - max(ANI(G_top, G_i)) for all G_i in secondary_genomes
```

Where:
- `G_top` = genome with highest bitscore hit
- `ANI(G_top, G_i)` = Average Nucleotide Identity between genome pairs (from pre-computed matrix)

**Interpretation:**

| Uncertainty | ANI Between Competitors | Biological Interpretation |
|-------------|-------------------------|---------------------------|
| U < 1.5% | ANI > 98.5% | Same species (confident placement) |
| 1.5% ≤ U < 5% | 95% < ANI ≤ 98.5% | Species boundary zone (ambiguous) |
| U ≥ 5% | ANI ≤ 95% | Different species (conserved gene region) |

**Biological rationale:** When competing genomes share high ANI (>98%), they represent the same species and placement ambiguity reflects strain-level variation. When competing genomes share low ANI (<95%), the read likely originates from a conserved gene present across multiple species.

### 2.3 Classification Decision Tree

Classification follows a hierarchical decision tree where uncertainty takes priority over novelty:

```
START: Read with BLAST hits
         |
         v
    Calculate N and U
         |
         v
   Is U ≥ 5%? ----YES----> AMBIGUOUS (conserved region)
         |
        NO
         v
   Is 1.5% ≤ U < 5%? --YES--> SPECIES BOUNDARY
         |
        NO (U < 1.5%)
         v
   Is N < 4%? ----YES----> KNOWN SPECIES
         |
        NO
         v
   Is 4% ≤ N < 20%? --YES--> NOVEL SPECIES
         |
        NO
         v
   Is 20% ≤ N ≤ 25%? --YES--> NOVEL GENUS
         |
        NO
         v
   UNCLASSIFIED (N > 25%)
```

**Classification categories:**

| Category | Criteria | Biological Interpretation |
|----------|----------|---------------------------|
| Known Species | N < 4% AND U < 1.5% | Matches characterized species |
| Novel Species | 4% ≤ N < 20% AND U < 1.5% | Significant divergence, confident placement |
| Novel Genus | 20% ≤ N ≤ 25% AND U < 1.5% | Genus-level divergence, confident placement |
| Ambiguous | U ≥ 5% OR 1.5% ≤ U < 5% | Conserved gene or species boundary |
| Species Boundary | 1.5% ≤ U < 5% | Matches multiple closely related species |
| Unclassified | N > 25% OR edge cases | Very high divergence, manual review needed |

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

---

## 3. Enhanced Scoring Metrics

### 3.1 Confidence Score

The confidence score (0-100) quantifies classification reliability without changing the classification category.

**Formula:**

```
confidence_score = margin_score + placement_score + alignment_score
```

**Component 1: Margin from threshold boundaries (0-40 points)**

```python
if taxonomic_call == "Known Species":
    margin = novelty_known_max - novelty_index  # Distance from 5% boundary
    margin_score = min(40, (margin / 5.0) × 40)

elif taxonomic_call == "Novel Species":
    margin_lower = novelty_index - 5.0   # Distance from lower boundary
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
if placement_uncertainty < 2.0:
    uncertainty_margin = 2.0 - placement_uncertainty
    margin_score += min(10, (uncertainty_margin / 2.0) × 10)
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

**Component 3: Alignment quality proxy (0-20 points)**

```python
alignment_score = min(20, max(0, (top_hit_identity - 70) / 25 × 20))
```

**Score interpretation:**

| Score Range | Interpretation | Recommended Action |
|-------------|----------------|-------------------|
| 80-100 | High confidence | Suitable for publication |
| 60-79 | Moderate confidence | Review recommended |
| 40-59 | Low confidence | Manual verification needed |
| < 40 | Very low confidence | Treat with caution |

### 3.2 Inferred Uncertainty for Single-Hit Reads

For reads with only one BLAST hit, placement uncertainty cannot be measured from ANI between competing genomes. Instead, uncertainty is inferred from the novelty level.

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

**Output field:** `uncertainty_type` = "inferred" or "measured"

### 3.3 Alignment Quality Score

The alignment quality score incorporates BLAST statistics beyond percent identity.

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

### 3.5 Discovery Score

The discovery score prioritizes novel findings for experimental validation. Only calculated for Novel Species and Novel Genus classifications.

**Formula:**

```
discovery_score = novelty_pts + quality_pts + confidence_pts
```

**Components:**

```python
# Novelty component (0-40 points)
if taxonomic_call == "Novel Species":
    # 5-20% novelty → 15-37.5 points
    novelty_pts = 15 + (novelty_index - 5) × 1.5

elif taxonomic_call == "Novel Genus":
    # 20-25% novelty → 30-40 points
    novelty_pts = 30 + (novelty_index - 20) × 2.0

novelty_pts = min(40, novelty_pts)

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
cluster_key = (best_match_genome, novelty_band, taxonomic_call)
```

Where:
- `best_match_genome` = closest reference genome accession
- `novelty_band` = floor(novelty_index / 5) × 5 (5% intervals: 5-10%, 10-15%, etc.)
- `taxonomic_call` = "Novel Species" or "Novel Genus"

**Algorithm steps:**

1. Filter to novel reads only (`taxonomic_call ∈ {"Novel Species", "Novel Genus"}`)
2. Compute novelty band for each read
3. Group by clustering key
4. Filter clusters with `read_count ≥ min_cluster_size` (default: 3)
5. Generate cluster IDs: `NSP_001`, `NSP_002`... for Novel Species; `NGN_001`... for Novel Genus
6. Assign confidence ratings
7. Compute cluster statistics

**Cluster statistics:**

| Metric | Calculation |
|--------|-------------|
| read_count | Number of reads in cluster |
| mean_novelty_index | Mean of read novelty values |
| novelty_min, novelty_max | Range of novelty in cluster |
| mean_placement_uncertainty | Mean uncertainty (using effective uncertainty) |
| mean_discovery_score | Mean discovery score (if enhanced scoring enabled) |

### 4.2 Cluster Confidence Ratings

Clusters are assigned confidence ratings based on multiple criteria.

**Rating criteria:**

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

For reads with `num_ambiguous_hits ≤ 1`, the effective uncertainty used for cluster statistics is:

```python
if placement_uncertainty < 0.1:  # Essentially zero (single hit)
    effective_uncertainty = min(novelty_index × 0.4, 15.0)
else:
    effective_uncertainty = placement_uncertainty
```

---

## 5. Coverage Weighting

Coverage weighting adjusts hit selection to prioritize longer alignments over short conserved domains. This is enabled by default (`coverage_weight_mode="linear"`) to reduce conserved-gene bias in classification.

**Weighted score calculation:**

```python
coverage = (qend - qstart + 1) / read_length
weight = calculate_weight(coverage, mode, strength)
weighted_score = bitscore × weight
```

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

**AAI genus boundaries** (Riesco & Trujillo, 2024):
- AAI > 65%: Same genus
- AAI 58-65%: Genus boundary zone
- AAI < 58%: Different genus

**Confidence score scaling (protein mode):**

| Parameter | Nucleotide | Protein |
|-----------|------------|---------|
| margin_divisor_known | 5.0 | 10.0 |
| margin_divisor_novel_species | 7.5 | 7.5 |
| margin_divisor_novel_genus | 5.0 | 7.5 |
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
| Filter rate | `filtered_alignments / total_alignments` | > 50% |
| Genome coverage | `genomes_in_alignments / genomes_in_ani` | < 50% |
| Single-hit fraction | `single_hit_reads / total_reads` | > 80% |
| Mean alignment length | Mean of BLAST alignment lengths | (informational) |
| Mean identity | Mean of BLAST percent identities | < 80% |

**Filter rate** indicates how many alignments are removed by quality thresholds (e.g., minimum alignment fraction). A high filter rate suggests either stringent filters or low-quality alignments.

**Genome coverage** measures the fraction of ANI matrix genomes that appear in the alignment file. Low coverage suggests that many reference genomes have no matching reads, which may affect classification accuracy.

**Single-hit fraction** indicates how many reads align to only one genome. For these reads, placement uncertainty cannot be directly measured from competing hits and must be inferred from the novelty level (see Section 3.2).

### 7.2 Post-Classification QC

Post-classification metrics assess the overall classification quality:

| Metric | Calculation | Warning Threshold |
|--------|-------------|-------------------|
| Ambiguous fraction | `ambiguous_reads / total_reads` | > 50% |
| Low confidence fraction | `reads with confidence < 0.5 / total` | > 30% |
| Novel fraction | `(novel_species + novel_genus) / total` | (informational) |

QC metrics are saved as JSON when `--qc-output` is specified.

---

## 8. Threshold Sensitivity Analysis

Sensitivity analysis evaluates how classification counts change across a range of threshold values. This helps assess whether results are robust to the specific threshold choice or whether small threshold changes lead to large shifts in classification counts.

**Method:**

The analysis sweeps the `novelty_known_max` threshold across a user-specified range (default: 2-8%) with proportional co-variation of uncertainty thresholds (default: 0.5-3.0%). At each step, all reads are re-classified using the threshold-extraction approach (no re-computation of alignment metrics).

```
for novelty_thresh in linspace(2.0, 8.0, steps=9):
    uncertainty_thresh = proportional(0.5, 3.0, steps=9)
    counts[category] = reclassify(reads, novelty_thresh, uncertainty_thresh)
```

**Output:** A JSON or chart showing category counts as a function of threshold position. Stable results (flat lines) indicate robustness; large count changes near the default threshold indicate sensitivity.

**Interpretation:**
- Flat count profiles across thresholds suggest classification is robust
- Sharp transitions at or near the default threshold indicate threshold-dependent results
- Gradual transitions are expected near true biological boundaries

---

## 9. Adaptive Thresholds from ANI Distribution

The default species boundary (96% ANI) is appropriate for many bacterial lineages but may not fit all taxa equally. Rapidly evolving lineages may have narrower within-species ANI distributions, while conserved lineages may show broader distributions.

### 9.1 Gaussian Mixture Model

When `--adaptive-thresholds` is enabled, the species boundary is estimated from the actual ANI matrix using a 2-component Gaussian Mixture Model:

**Step 1:** Extract upper-triangle pairwise ANI values from the ANI matrix (excluding self-comparisons and zero entries).

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

The boundary is clamped to the range 80-99% ANI.

### 9.2 Separation Quality

The GMM separation quality is assessed using the standardized gap between component means:

```
separation = |mu_high - mu_low| / (sigma_low + sigma_high)
```

| Separation | Assessment | Action |
|------------|------------|--------|
| < 1.0 | Components not distinct | Falls back to default threshold |
| 1.0 - 2.0 | Moderate separation | Uses GMM boundary with moderate confidence |
| > 2.0 | Well separated | Uses GMM boundary with high confidence |

### 9.3 Fallback Conditions

The adaptive method falls back to the default threshold (96% ANI) when:
- Fewer than 5 genomes in the ANI matrix
- Fewer than 3 non-zero pairwise ANI values
- GMM separation quality < 1.0
- GMM fails to converge
- scikit-learn is not installed

**Dependency:** Requires `scikit-learn>=1.3.0` (install with `pip install metadarkmatter[adaptive]`).

---

## 10. Bayesian Confidence Framework

The Bayesian framework replaces hard threshold boundaries with posterior probabilities P(category | novelty, uncertainty). Near classification boundaries, posteriors will be split across categories (e.g., 55% Novel Species, 40% Known Species), making the classification uncertainty explicit.

### 10.1 Likelihood Model

Each classification category is modeled as a 2D Gaussian distribution centered on its expected (novelty, uncertainty) values:

```
L(category | N, U) = exp(-0.5 * (z_N^2 + z_U^2))

where:
    z_N = (N - mu_N) / sigma_N
    z_U = (U - mu_U) / sigma_U
```

**Category parameters (nucleotide mode defaults):**

| Category | mu_N | sigma_N | mu_U | sigma_U |
|----------|------|---------|------|---------|
| Known Species | 2.0 | 2.5 | 0.75 | 1.25 |
| Novel Species | 12.0 | 5.83 | 0.75 | 1.25 |
| Novel Genus | 22.5 | 2.17 | 0.75 | 1.25 |
| Ambiguous | 15.0 | 10.0 | 5.0 | 5.0 |

Parameters are derived from the scoring configuration thresholds: means are placed at category centers, and sigmas are scaled from the category width.

### 10.2 Posterior Computation

Posteriors are computed using Bayes' rule with uniform priors:

```
P(category | N, U) = L(category | N, U) / sum(L(all categories | N, U))
```

**Ambiguous category boosting:** Two additional factors increase the Ambiguous likelihood:
- If `identity_gap < 2%` and `num_hits > 1`: likelihood multiplied by 3x
- If `placement_uncertainty > uncertainty_conserved_min`: likelihood multiplied by 2x

These factors ensure that reads with ambiguous alignment characteristics are appropriately assigned higher Ambiguous posteriors.

**Output columns:**
- `p_known_species`, `p_novel_species`, `p_novel_genus`, `p_ambiguous` (sum to 1.0)
- `bayesian_category`: Maximum A Posteriori (MAP) classification
- `posterior_entropy`: Shannon entropy of the posterior distribution

### 10.3 Shannon Entropy

The posterior entropy quantifies classification confidence as a single scalar:

```
H = -sum(p_i * log2(p_i)) for all categories where p_i > 0
```

For 4 categories, entropy ranges from 0 (all probability mass on one category) to 2.0 (uniform distribution across all categories).

| Entropy | Interpretation |
|---------|---------------|
| H < 0.5 | Confident assignment to one category |
| 0.5 ≤ H < 1.0 | Moderate confidence; one category dominant |
| 1.0 ≤ H < 1.5 | Low confidence; multiple categories plausible |
| H ≥ 1.5 | Near-uniform; read lies at a classification boundary |

### 10.4 Interpretation

The Bayesian framework is complementary to threshold-based classification. The MAP category often agrees with the threshold-based `taxonomic_call`, but near boundaries the posteriors reveal uncertainty that the threshold-based system hides.

**Example:** A read at (N=3.8%, U=1.4%) falls just within the Known Species boundary by threshold classification, but the Bayesian posteriors might be: p_known=0.52, p_novel=0.38, p_ambiguous=0.08, p_genus=0.02. The moderate entropy (H=1.3) flags this as a boundary case.

**Vectorized implementation:** The classify_dataframe method uses numpy broadcasting for efficient computation over large DataFrames without Python loops.

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

All classification thresholds are defined in `src/metadarkmatter/core/constants.py`:

```python
# Novelty Index Thresholds
NOVELTY_KNOWN_MAX = 4.0          # 96% ANI species boundary (Jain et al. 2018)
NOVELTY_NOVEL_SPECIES_MIN = 4.0
NOVELTY_NOVEL_SPECIES_MAX = 20.0
NOVELTY_NOVEL_GENUS_MIN = 20.0
NOVELTY_NOVEL_GENUS_MAX = 25.0

# Placement Uncertainty Thresholds
UNCERTAINTY_CONFIDENT_MAX = 1.5   # Tightened for fewer false-confident calls
UNCERTAINTY_CONSERVED_MIN = 5.0
UNCERTAINTY_AMBIGUOUS_MIN = 1.5
UNCERTAINTY_AMBIGUOUS_MAX = 5.0

# ANI Boundaries
ANI_SPECIES_BOUNDARY_LOW = 95.0
ANI_SPECIES_BOUNDARY_HIGH = 96.0
ANI_GENUS_BOUNDARY = 75.0

# AAI Boundaries
AAI_GENUS_BOUNDARY_HIGH = 65.0
AAI_GENUS_BOUNDARY_LOW = 58.0
```

**Default configuration changes (v0.2):**

| Parameter | Previous Default | Current Default | Rationale |
|-----------|-----------------|-----------------|-----------|
| `novelty_known_max` | 5.0 | 4.0 | 96% ANI species boundary |
| `uncertainty_known_max` | 2.0 | 1.5 | Stricter confident placement |
| `coverage_weight_mode` | "none" | "linear" | Reduce conserved-gene bias |
| `min_alignment_fraction` | 0.0 | 0.3 | Filter low-coverage alignments |
| Enhanced scoring | opt-in flag | always on | Inferred uncertainty, quality metrics |
| Single-hit inference | opt-in flag | always on | No more 0% uncertainty for single hits |
