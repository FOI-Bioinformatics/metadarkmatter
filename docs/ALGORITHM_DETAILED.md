# Classification Algorithm - Detailed Explanation

Complete technical explanation of how metadarkmatter classifies reads using ANI-weighted placement uncertainty.

**For comprehensive methods documentation with all calculations at scientific standard, see [METHODS.md](METHODS.md).**

**For details on BLAST/MMseqs2 output statistics (pident, bitscore, evalue, etc.), see [ALIGNMENT_OUTPUT_STATISTICS.md](ALIGNMENT_OUTPUT_STATISTICS.md).**

## Overview

The classification algorithm uses a two-metric decision tree:
1. **Novelty Index (N)**: Measures sequence divergence from closest reference
2. **Placement Uncertainty (U)**: Measures ambiguity when multiple genomes compete

Both metrics are calculated from BLAST alignment results combined with pre-computed ANI (or AAI) matrix lookups.

---

## Algorithm Steps

### Step 1: Parse BLAST Results and Group by Read

For each read with BLAST hits:

```
Input: BLAST tabular output with columns:
  qseqid, sseqid, pident, length, mismatch, gapopen,
  qstart, qend, sstart, send, evalue, bitscore

Process:
  1. Group all hits by read ID (qseqid)
  2. For each read, collect all hits with their:
     - Target genome (extracted from sseqid prefix)
     - Percent identity (pident)
     - Bitscore
```

### Step 2: Identify Top Hit and Ambiguous Competitors

For each read:

```python
# Find top hit (highest bitscore)
top_hit = max(hits, key=lambda h: h.bitscore)
G_top = top_hit.genome
pident_top = top_hit.pident

# Define bitscore threshold (95% of top bitscore)
threshold = 0.95 * top_hit.bitscore

# Find ambiguous hits (competitive recruitment)
ambiguous_hits = [h for h in hits if h.bitscore >= threshold]
competing_genomes = {h.genome for h in ambiguous_hits}  # Unique genomes
```

**Key concept**: Ambiguous hits represent genomes that matched the read nearly as well as the top hit. If multiple genomes compete, we need to assess whether they're similar (same species) or divergent (conserved gene).

#### Optional: Coverage-Weighted Hit Selection

When coverage weighting is enabled (`--coverage-weight-mode linear|log|sigmoid`), the top hit selection uses a weighted score instead of raw bitscore:

```python
def calculate_weighted_score(hit, mode, strength):
    """Calculate coverage-weighted bitscore."""
    coverage = (hit.qend - hit.qstart + 1) / hit.qlen
    weight = calculate_coverage_weight(coverage, mode, strength)
    return hit.bitscore * weight

def calculate_coverage_weight(coverage, mode, strength):
    """Calculate weight factor based on coverage."""
    min_weight = 1.0 - strength
    max_weight = 1.0 + strength
    weight_range = max_weight - min_weight

    if mode == "linear":
        normalized = coverage
    elif mode == "log":
        normalized = log(1 + 9 * coverage) / log(10)
    elif mode == "sigmoid":
        normalized = 1.0 / (1.0 + exp(-10.0 * (coverage - 0.6)))

    return min_weight + weight_range * normalized

# With coverage weighting enabled:
top_hit = max(hits, key=lambda h: calculate_weighted_score(h, mode, strength))
```

**Effect:** Short high-identity hits (conserved domains) are penalized relative to longer alignments that span more of the read. This helps distinguish true species matches from conserved gene fragments.

**When to use:** Enable coverage weighting when short conserved regions (e.g., 16S rRNA) are dominating classification despite lower overall alignment coverage.

### Step 3: Calculate Novelty Index (N)

```python
N = 100 - pident_top
```

**Interpretation:**
- N < 5%: Read shares >95% identity with reference (at species boundary)
- N = 5-20%: Read shares 80-95% identity (below species boundary)
- N = 20-25%: Read shares 75-80% identity (genus-level divergence)
- N > 25%: Very high divergence

### Step 4: Calculate Placement Uncertainty (U)

This is where the ANI (or AAI) matrix is used:

```python
# Remove top genome from competing set
secondary_genomes = competing_genomes - {G_top}

if len(secondary_genomes) == 0:
    # Only one genome matched - no uncertainty
    U = 0.0
else:
    # Look up ANI between top genome and each competitor
    ani_values = []
    for G_i in secondary_genomes:
        ani = ANI_matrix.lookup(G_top, G_i)  # Query pre-computed matrix
        ani_values.append(ani)

    # Placement uncertainty = 100 - highest ANI to any competitor
    max_ani = max(ani_values)
    U = 100 - max_ani
```

**Key concept**: If the top genome and competitors share high ANI (>98%), they're the same species and placement is confident. If they share low ANI (<95%), they're different species and the read matches a conserved region.

**ANI matrix lookup:**
- Matrix is symmetric: ANI(A, B) = ANI(B, A)
- Self-ANI is always 100%
- Missing pairs default to 0% (different species assumed)

### Step 5: Apply Decision Tree

Classification uses a hierarchical decision tree where uncertainty takes priority:

```
START
  |
  v
Is U >= 5%?  # Competitors are different species (ANI < 95%)
  |
  YES → CONSERVED REGION (conserved gene across taxa)
  |
  NO
  v
Is 2% <= U < 5%?  # Competitors in species boundary zone (95-98% ANI)
  |
  YES → AMBIGUOUS (species boundary, unclear placement)
  |
  NO (U < 2%)  # Competitors are same species (ANI > 98%) OR no competitors
  v
Is N < 5%?  # High identity (>95%)
  |
  YES → KNOWN SPECIES
  |
  NO
  v
Is 5% <= N < 20%?  # Moderate divergence (80-95% identity)
  |
  YES → NOVEL SPECIES (confident, divergent from references)
  |
  NO
  v
Is 20% <= N <= 25%?  # High divergence (75-80% identity)
  |
  YES → NOVEL GENUS (confident, genus-level novelty)
  |
  NO
  v
UNCLASSIFIED (N > 25% or edge case)
```

**Python-like pseudocode:**

```python
def classify_read(N, U):
    # Priority 1: Check uncertainty (conserved regions)
    if U >= 5.0:
        return "Conserved Region"

    # Priority 2: Check species boundary zone
    if 2.0 <= U < 5.0:
        return "Ambiguous"

    # Priority 3: Confident placement (U < 2%), check novelty
    if N < 5.0:
        return "Known Species"
    elif 5.0 <= N < 20.0:
        return "Novel Species"
    elif 20.0 <= N <= 25.0:
        return "Novel Genus"
    else:
        return "Unclassified"
```

---

## ANI vs AAI Usage

### Nucleotide Mode (Default)

**Uses:** BLASTN alignment + ANI matrix

```bash
metadarkmatter score classify \
    --alignment sample.blast.tsv.gz \
    --ani ani_matrix.csv \
    --output classifications.csv
```

**Algorithm:**
- ANI matrix computed with skani or fastANI (nucleotide genome-to-genome)
- U = 100 - max(ANI between competing genomes)
- Thresholds calibrated to nucleotide divergence rates

### Protein Mode

**Uses:** BLASTX alignment + AAI matrix (optionally with ANI as fallback)

```bash
metadarkmatter score classify \
    --alignment sample.blastx.tsv.gz \
    --ani ani_matrix.csv \
    --aai aai_matrix.csv \
    --alignment-mode protein \
    --output classifications.csv
```

**Algorithm:**
- AAI matrix computed with Diamond BLASTP (protein-to-protein)
- **If AAI available**: U = 100 - max(AAI between competing genomes)
- **If AAI missing for a pair, fallback to ANI**: U = 100 - max(ANI for that pair)
- Thresholds adjusted for slower protein divergence:
  - Known Species: N < 10%, U < 5% (vs 5%/2% for nucleotide)
  - Novel Species: 10% <= N < 25%, U < 5%
  - Novel Genus: 25% <= N <= 40%, U < 5%

**Why use AAI?**
- Proteins evolve more slowly than nucleotides
- For genus-level novelty, nucleotide alignment may degrade while protein alignment remains
- AAI provides more accurate similarity estimates for divergent taxa

### Dual Matrix Strategy

When both `--ani` and `--aai` are provided:

```python
def calculate_uncertainty(G_top, secondary_genomes, ani_matrix, aai_matrix, mode):
    uncertainties = []

    for G_i in secondary_genomes:
        if mode == "protein" and aai_matrix.has_pair(G_top, G_i):
            # Protein mode: prefer AAI if available
            similarity = aai_matrix.lookup(G_top, G_i)
        elif ani_matrix.has_pair(G_top, G_i):
            # Fallback to ANI
            similarity = ani_matrix.lookup(G_top, G_i)
        else:
            # No data: assume different species
            similarity = 0.0

        uncertainties.append(100 - similarity)

    return min(uncertainties) if uncertainties else 0.0
```

**Decision logic:**
1. **Protein mode** (`--alignment-mode protein`):
   - Try AAI matrix first
   - Fallback to ANI if AAI missing
   - Use protein thresholds

2. **Nucleotide mode** (default):
   - Use ANI matrix only
   - Ignore AAI even if provided
   - Use nucleotide thresholds

---

## Example Walkthrough

### Example 1: Known Species

**Input:**
```
read_001 hits:
  GCF_000195955.1 (98.5% identity, bitscore 450)
  GCF_000242755.1 (97.2% identity, bitscore 435)  # 435/450 = 0.967 (>95%)
```

**Calculation:**
```python
# Step 1: Top hit
G_top = "GCF_000195955.1"
pident_top = 98.5

# Step 2: Ambiguous hits (bitscore >= 0.95 * 450 = 427.5)
competing_genomes = {"GCF_000195955.1", "GCF_000242755.1"}
secondary = {"GCF_000242755.1"}

# Step 3: Novelty
N = 100 - 98.5 = 1.5%

# Step 4: Uncertainty
ANI(GCF_000195955.1, GCF_000242755.1) = 98.7%  # Look up in matrix
U = 100 - 98.7 = 1.3%

# Step 5: Classify
U < 2% and N < 5% → KNOWN SPECIES
```

**Interpretation:** Both genomes are the same species (98.7% ANI). Read matches with high identity.

### Example 2: Novel Species

**Input:**
```
read_002 hits:
  GCF_000195955.1 (89.5% identity, bitscore 380)
  GCF_000242755.1 (88.1% identity, bitscore 365)  # 365/380 = 0.961 (>95%)
```

**Calculation:**
```python
N = 100 - 89.5 = 10.5%
ANI(GCF_000195955.1, GCF_000242755.1) = 98.7%
U = 100 - 98.7 = 1.3%

# Classify
U < 2% and 5% <= N < 20% → NOVEL SPECIES
```

**Interpretation:** Competing genomes are same species, but read shows significant divergence (89.5%). Likely a novel species in this genus.

### Example 3: Conserved Region

**Input:**
```
read_003 hits:
  GCF_000195955.1 (95.0% identity, bitscore 400)
  GCF_003456789.1 (94.5% identity, bitscore 395)  # Different genus
```

**Calculation:**
```python
N = 100 - 95.0 = 5.0%
ANI(GCF_000195955.1, GCF_003456789.1) = 78.3%  # Different genera
U = 100 - 78.3 = 21.7%

# Classify
U >= 5% → CONSERVED REGION
```

**Interpretation:** Read matches genomes from different genera with similar scores. This is a conserved gene (e.g., 16S rRNA, housekeeping gene).

### Example 4: Novel Genus

**Input:**
```
read_004 hits:
  GCF_000195955.1 (82.0% identity, bitscore 320)
  (no other hits above threshold)
```

**Calculation:**
```python
N = 100 - 82.0 = 18.0%
secondary_genomes = {} (empty)
U = 0.0  # No competition

# Classify
U < 2% and 20% <= N <= 25% → NOVEL GENUS
```

**Interpretation:** Read shows genus-level divergence with confident placement to single genome lineage.

---

## Advanced Considerations

### Ambiguous Hit Threshold

The 95% bitscore threshold for competitive recruitment:
- **Too strict (e.g., 99%)**: May miss biologically relevant competitors
- **Too relaxed (e.g., 90%)**: May include spurious low-quality hits
- **95% (default)**: Balances sensitivity and specificity based on empirical testing

### Missing ANI Pairs

If ANI(G_top, G_i) is not in the matrix:
```python
# Conservative assumption: treat as different species
default_ani = 0.0
U = 100 - 0.0 = 100.0  # Maximum uncertainty
→ Classified as "Conserved Region"
```

**Best practice:** Validate ANI matrix covers all BLAST database genomes:
```bash
metadarkmatter ani validate --ani ani_matrix.csv --genomes genomes/
```

### Edge Cases

1. **No BLAST hits**: Read not classified (filtered out)
2. **All hits below threshold**: Only top hit considered, U = 0
3. **N > 25%**: Unclassified (very high divergence, manual review)
4. **Negative bitscore**: Invalid hit, filtered out

---

## Computational Complexity

**Per read:**
1. Parse hits: O(h) where h = number of hits for this read
2. Find top hit: O(h)
3. Filter ambiguous: O(h)
4. ANI lookups: O(k) where k = number of competing genomes (typically 1-10)
5. Classification: O(1)

**Total: O(h + k) ≈ O(h)** since k << h in practice

**For dataset:**
- Reads: R
- Average hits per read: H
- **Total complexity: O(R × H)**

**Bottleneck:** BLAST alignment (not classification)

---

## Implementation Notes

The actual implementation is in:
- `src/metadarkmatter/core/classification/classifiers/base.py`
- `src/metadarkmatter/core/classification/metrics.py`

Key optimization:
- ANI matrix loaded once into memory (dict lookup: O(1))
- BLAST results streamed (doesn't load entire file into RAM)
- Parallel processing available for large datasets

---

## Summary

**The algorithm in one sentence:**

For each read, find the closest reference genome and competing genomes, then use ANI (or AAI) between these genomes to determine if the read's divergence represents novel diversity or just a conserved gene region.

**Key insight:**

High sequence identity alone doesn't guarantee novel diversity - a read might match multiple divergent genomes equally well (conserved region). The ANI-weighted uncertainty metric distinguishes confident placements from ambiguous ones.
