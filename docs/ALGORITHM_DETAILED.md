# Classification Algorithm - Detailed Explanation

Complete technical explanation of how metadarkmatter classifies reads using ANI-weighted placement uncertainty.

**For comprehensive methods documentation with all calculations at scientific standard, see [METHODS.md](METHODS.md).**

**For details on BLAST/MMseqs2 output statistics (pident, bitscore, evalue, etc.), see [ALIGNMENT_OUTPUT_STATISTICS.md](ALIGNMENT_OUTPUT_STATISTICS.md).**

## Overview

The classification algorithm is based on two metrics:
1. **Novelty Index (N)**: Measures sequence divergence from closest reference
2. **Placement Uncertainty (U)**: Measures ambiguity when multiple genomes compete

Both metrics are calculated from BLAST/MMseqs2 alignment results combined with pre-computed ANI matrix lookups. The primary call (`taxonomic_call`) is the maximum a posteriori category of a six-category Bayesian model in (N, U) space, refined by a discrete second stage; the hard-threshold decision tree of Step 5 is reported alongside as `legacy_call`. The two do not have identical decision boundaries ([METHODS.md](METHODS.md) section 2.5).

The description below follows the CLI classifier (`VectorizedClassifier`, `core/classification/classifiers/vectorized.py`).

---

## Algorithm Steps

### Step 1: Parse BLAST Results and Group by Read

For each read with BLAST hits:

```
Input: BLAST tabular output (outfmt 6) with columns:
  qseqid, sseqid, pident, length, mismatch, gapopen,
  qstart, qend, sstart, send, evalue, bitscore [, qlen]

Process:
  1. Normalise reverse-strand coordinates (qstart > qend -> swap)
  2. Filter hits: length >= 100 bp, alignment fraction
     (qend - qstart + 1) / qlen >= 0.3 (qend used when qlen is absent)
  3. Group all hits by read ID (qseqid)
  4. For each read, collect all hits with their:
     - Target genome (extracted from the sseqid prefix, {accession}|{contig})
     - Species representative of that genome (from genome_metadata.tsv;
       the genome itself when no representative column is present)
     - Percent identity (pident)
     - Bitscore and coverage-weighted bitscore
```

### Step 2: Identify Top Hit and Ambiguous Competitors

For each read:

```python
# Find top hit (highest coverage-weighted bitscore; ties: higher pident,
# then genome accession ascending)
top_hit = min(hits, key=lambda h: (-h.weighted_bitscore, -h.pident, h.genome))
G_top = top_hit.genome
pident_top = top_hit.pident

# Define bitscore threshold (95% of top weighted bitscore)
threshold = 0.95 * top_hit.weighted_bitscore

# Find ambiguous hits (competitive recruitment)
ambiguous_hits = [h for h in hits if h.weighted_bitscore >= threshold]
num_ambiguous_hits = len(ambiguous_hits)          # hits, not genomes
competing_reps = {rep(h.genome) for h in ambiguous_hits}
```

**Key concept**: Ambiguous hits represent genomes that matched the read nearly as well as the top hit. If multiple genomes compete, we need to assess whether they're similar (same species) or divergent (conserved gene).

#### Coverage-Weighted Hit Selection (default: linear)

Coverage weighting is enabled by default (`--coverage-weight-mode linear`, strength 0.5); `--coverage-weight-mode none` uses the raw bitscore. The weighted score is used both for top-hit selection and for the 95% window:

```python
def calculate_weighted_score(hit, mode, strength):
    """Calculate coverage-weighted bitscore."""
    coverage = min(1.0, max(0.0, (hit.qend - hit.qstart + 1) / hit.qlen))  # qend if no qlen
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

```

**Effect:** Short high-identity hits (conserved domains) are penalized relative to longer alignments that span more of the read. This helps distinguish true species matches from conserved gene fragments.

**When to change:** Use `none` to reproduce raw-bitscore behaviour, or `sigmoid` (as in the `coverage-strict` preset) when short conserved regions dominate despite low alignment coverage.

### Step 3: Calculate Novelty Index (N)

```python
N = 100 - pident_top
```

**Interpretation:**
- N < 4%: Read shares >96% identity with reference (at species boundary)
- N = 4-20%: Read shares 80-96% identity (below species boundary)
- N = 20-25%: Read shares 75-80% identity (genus-level divergence)
- N > 25%: Very high divergence

### Step 4: Calculate Placement Uncertainty (U)

This is where the ANI (or AAI) matrix is used:

```python
# Remove the top genome's representative from the competing set
secondary = competing_reps - {rep(G_top)}

if num_ambiguous_hits <= 1 or len(secondary) == 0:
    # Only one genome (representative) matched - no measurable uncertainty
    U = 0.0
elif uncertainty_mode == "second":            # default
    # ANI to the best-ranked secondary representative, using the same
    # order as top-hit selection (weighted bitscore, pident, accession)
    G_2 = first(sorted(secondary_hits, key=lambda h: (-h.weighted_bitscore, -h.pident, rep(h.genome))))
    U = 100 - ANI_matrix.lookup(rep(G_top), rep(G_2))
else:                                         # uncertainty_mode == "max"
    U = 100 - max(ANI_matrix.lookup(rep(G_top), G_i) for G_i in secondary)
```

**Key concept**: If the top genome and competitors share high ANI (>98.5%), they're the same species and placement is confident. If they share low ANI (<95%), they're different species and the read matches a conserved region.

**ANI matrix lookup:**
- Values must be on a 0-100 scale; a 0-1 scaled matrix is rejected
- Pairs supplied in only one direction are mirrored; asymmetric pairs are reported (`--strict-ani` turns this into an error)
- Self-ANI is always 100%
- Pairs present in the matrix but not computed by the ANI tool (stored as 0) are looked up as 70% (U = 30%)
- A competing genome absent from the matrix gives U = 100%

### Step 5: Classify

**Primary call (`taxonomic_call`).** The pair (N, U) is scored under six 2D Gaussians (Known Species, Novel Species, Novel Genus, Species Boundary, Ambiguous, Unclassified) whose means and widths are derived from the thresholds below. The Ambiguous prior is multiplied by 2.0 when the identity gap is < 2% with multiple hits, and by 1.5 for single-hit reads with N >= 4%. The MAP category is then refined: Ambiguous with hits across genera becomes Conserved Region; Novel Genus with genus uncertainty >= 10% and at least one secondary genome becomes Ambiguous Within Genus. `confidence_score = (1 - H / log2(6)) * 100`, where H is the posterior entropy. Full details: [METHODS.md](METHODS.md) section 10.

**Legacy call (`legacy_call`).** The hard-threshold cascade, in rule order (nucleotide defaults):

```
START
  |
  v
identity_gap < 2% and more than one ambiguous hit?
  |
  YES -> AMBIGUOUS
  |
  NO
  v
U >= 5%?  # Competitors are different species (ANI <= 95%)
  |
  YES -> CONSERVED REGION if hits span genera, otherwise AMBIGUOUS
  |
  NO
  v
1.5% <= U < 5%?  # Competitors in species boundary zone (95-98.5% ANI)
  |
  YES -> SPECIES BOUNDARY
  |
  NO (U < 1.5%)  # Competitors are same species (ANI > 98.5%) OR no competitors
  v
N < 4%?  # High identity (>96%)
  |
  YES -> KNOWN SPECIES
  |
  NO
  v
Single ambiguous hit and inferred uncertainty >= 10% (i.e. N >= 7%)?
  |
  YES -> AMBIGUOUS
  |
  NO
  v
4% <= N < 20%?  # Moderate divergence (80-96% identity)
  |
  YES -> NOVEL SPECIES
  |
  NO
  v
20% <= N <= 25% and (genus uncertainty < 10% or no secondary genomes)?
  |
  YES -> NOVEL GENUS
  |
  NO
  v
UNCLASSIFIED (N > 25% or edge case)
```

**Python-like pseudocode (legacy cascade):**

```python
def legacy_call(N, U, num_hits, identity_gap, scope, genus_U, n_secondary):
    if identity_gap is not None and num_hits > 1 and identity_gap < 2.0:
        return "Ambiguous"
    if U >= 5.0:
        return "Conserved Region" if scope == "across_genera" else "Ambiguous"
    if U >= 1.5:
        return "Species Boundary"
    if N < 4.0:
        return "Known Species"
    if num_hits <= 1 and inferred_uncertainty(N) >= 10.0:   # N >= 7
        return "Ambiguous"
    if 4.0 <= N < 20.0:
        return "Novel Species"
    if 20.0 <= N <= 25.0 and (genus_U < 10.0 or n_secondary == 0):
        return "Novel Genus"
    return "Unclassified"
```

**Where the two calls differ.** With the default parameters the Bayesian boundaries at U = 0 lie at N ~ 6.2% (Known -> Novel Species), ~18.6% (Novel Species -> Novel Genus) and ~27.6% (Novel Genus -> Unclassified), rather than 4/20/25%. The Bayesian model does not use the inferred uncertainty of single-hit reads (it sees U = 0), so single-hit reads with N >= 7% are Novel Species / Novel Genus in `taxonomic_call` but Ambiguous in `legacy_call`. The prior boosts do not reproduce the identity-gap rule: a multi-hit read with a gap < 2% at low U keeps its Known/Novel MAP call. See [METHODS.md](METHODS.md) section 2.5.

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
- ANI matrix computed with fastANI or skani (nucleotide genome-to-genome)
- U = 100 - ANI(top genome, second genome) by default (`--uncertainty-mode second`); `--uncertainty-mode max` uses the maximum ANI to any competing genome
- Thresholds calibrated to nucleotide divergence rates

### Protein Mode

**Uses:** BLASTX alignment + ANI matrix (an AAI matrix can be supplied but is not used by `score classify`; see below)

```bash
metadarkmatter score classify \
    --alignment sample.blastx.tsv.gz \
    --ani ani_matrix.csv \
    --alignment-mode protein \
    --output classifications.csv
```

**Algorithm:**
- N = 100 - protein identity of the best BLASTX hit
- U is computed from the ANI matrix exactly as in nucleotide mode
- Thresholds adjusted for slower protein divergence (fixed constants from `core/protein_constants.py`; custom, YAML and adaptive thresholds are ignored in protein mode):
  - Known Species: N < 10%, U < 5% (vs 4%/1.5% for nucleotide)
  - Novel Species: 10% <= N < 25%, U < 5%
  - Novel Genus: 25% <= N <= 40%, U < 5%
  - Conserved Region / Ambiguous: U >= 10%

**Why use protein mode?**
- Proteins evolve more slowly than nucleotides
- For genus-level novelty, nucleotide alignment may degrade while protein alignment remains

### AAI Matrix

`score classify --aai aai_matrix.csv` loads the AAI matrix and passes it to the CLI classifier, which builds an AAI lookup table but does not currently use it: placement uncertainty and all classification rules use the ANI matrix only, in both modes.

AAI is used by the programmatic `ANIWeightedClassifier` (`core/classification/classifiers/base.py`) when `use_aai_for_genus` is true: for reads with N in the novel-genus range (20-25% nucleotide) and U below the confident-placement maximum, the maximum AAI from the top genome to the other hit genomes overrides the Bayesian call (AAI > 65% -> Novel Species; 58-65% -> Ambiguous Within Genus; < 58% -> Novel Genus; Riesco & Trujillo 2024). This is not reachable from `score classify`.

---

## Example Walkthrough

The values below were computed with the default configuration (`BayesianClassifier(ScoringConfig()).classify_primary()` and `apply_legacy_thresholds()`). All hits are assumed to cover the full read, so coverage weighting scales every bitscore by the same factor and the 95% window can be read off the raw bitscores.

### Example 1: Known Species

**Input:**
```
read_001 hits:
  GCF_000195955.1 (98.5% identity, bitscore 450)
  GCF_000242755.1 (96.2% identity, bitscore 435)  # 435/450 = 0.967 (>95%)
```

**Calculation:**
```python
# Step 1: Top hit
G_top = "GCF_000195955.1"
pident_top = 98.5

# Step 2: Ambiguous hits (bitscore >= 0.95 * 450 = 427.5)
num_ambiguous_hits = 2
secondary = {"GCF_000242755.1"}
identity_gap = 98.5 - 96.2 = 2.3%

# Step 3: Novelty
N = 100 - 98.5 = 1.5%

# Step 4: Uncertainty
ANI(GCF_000195955.1, GCF_000242755.1) = 98.7%  # Look up in matrix
U = 100 - 98.7 = 1.3%

# Step 5: Classify
taxonomic_call: Known Species (p_known = 0.879, p_novel_species = 0.076,
                H = 0.69, confidence_score = 73.5)
legacy_call:    U < 1.5% and N < 4% -> Known Species
```

**Interpretation:** Both genomes are the same species (98.7% ANI). Read matches with high identity. If the second hit had 97.2% identity (gap 1.3%), `legacy_call` would be Ambiguous (identity-gap rule), whereas `taxonomic_call` would remain Known Species (confidence_score 71.0).

### Example 2: Novel Species

**Input:**
```
read_002 hits:
  GCF_000195955.1 (89.5% identity, bitscore 380)
  GCF_000242755.1 (87.1% identity, bitscore 365)  # 365/380 = 0.961 (>95%)
```

**Calculation:**
```python
N = 100 - 89.5 = 10.5%
identity_gap = 2.4%
ANI(GCF_000195955.1, GCF_000242755.1) = 98.7%
U = 100 - 98.7 = 1.3%

# Classify
taxonomic_call: Novel Species (p_novel_species = 0.666, p_species_boundary = 0.277,
                H = 1.17, confidence_score = 54.9)
legacy_call:    U < 1.5% and 4% <= N < 20% -> Novel Species
```

**Interpretation:** Competing genomes are same species, but read shows significant divergence (89.5%). Likely a novel species in this genus. The posterior assigns about a quarter of the mass to Species Boundary because U = 1.3% is close to the Species Boundary Gaussian.

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
identity_gap = 0.5%
ANI(GCF_000195955.1, GCF_003456789.1) = 78.3%  # < 80%: ambiguity_scope = across_genera
U = 100 - 78.3 = 21.7%

# Classify
taxonomic_call: Stage 1 MAP Ambiguous (p_ambiguous = 1.00)
                + across_genera -> Conserved Region (confidence_score = 100.0)
legacy_call:    identity gap < 2% with 2 hits -> Ambiguous (rule evaluated
                before the Conserved Region rule)
```

**Interpretation:** Read matches genomes from different genera with similar scores. This is a conserved gene (e.g., 16S rRNA, housekeeping gene). With an identity gap of 2% or more, `legacy_call` would also be Conserved Region.

### Example 4: Novel Genus (single hit)

**Input:**
```
read_004 hits:
  GCF_000195955.1 (78.0% identity, bitscore 320)
  (no other hits above threshold)
```

**Calculation:**
```python
N = 100 - 78.0 = 22.0%
num_ambiguous_hits = 1
U = 0.0  # No competition
inferred_uncertainty = 23.0 + (22.0 - 20.0) * 1.5 = 26.0%

# Classify
taxonomic_call: Novel Genus (p_novel_genus = 0.869, H = 0.77, confidence_score = 70.1)
legacy_call:    single hit with inferred uncertainty >= 10% -> Ambiguous
```

**Interpretation:** Read shows genus-level divergence and only one reference genome lineage was hit. The two calls differ because the Bayesian model sees U = 0 and does not use the inferred uncertainty, whereas the legacy cascade treats single-hit reads with N >= 7% as Ambiguous ([METHODS.md](METHODS.md) section 2.5).

---

## Advanced Considerations

### Ambiguous Hit Threshold

The 95% bitscore threshold (`--bitscore-threshold`) for competitive recruitment, applied to the coverage-weighted bitscore:
- **Too strict (e.g., 99%)**: May miss biologically relevant competitors
- **Too relaxed (e.g., 90%)**: May include spurious low-quality hits
- **95% (default)**: A heuristic compromise between sensitivity and specificity; it has not been formally benchmarked

### Missing ANI Pairs

```python
# Pair in the matrix but not computed by fastANI/skani (stored as 0):
ani = 70.0            # ANI_DEFAULT_UNRELATED
U = 100 - 70.0 = 30.0

# Competing genome absent from the ANI matrix:
U = 100.0             # maximum uncertainty
```

Either value is far above the 5% conserved-region threshold, so such reads are called Ambiguous or Conserved Region (depending on the scope and the MAP).

**Best practice:** Validate ANI matrix covers all BLAST database genomes:
```bash
metadarkmatter ani validate --ani ani_matrix.csv --blast sample.blast.tsv.gz
```

### Edge Cases

1. **No hits after filtering**: Read not classified (absent from output)
2. **All other hits below threshold**: Only top hit considered, U = 0
3. **N > 25%**: Unclassified in `legacy_call`; in `taxonomic_call` Novel Genus remains the MAP up to N ~ 27.6%
4. **Reverse-strand hits** (`qstart > qend`): coordinates are normalised before coverage is computed
5. **Off-target reads** (with `--target-family` or a metadata family): Off-target, U = null, confidence_score = 10.0

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
- `src/metadarkmatter/core/classification/classifiers/vectorized.py` (`VectorizedClassifier`, used by `score classify`)
- `src/metadarkmatter/core/classification/bayesian.py` (posteriors, Stage 2 refinement, entropy)
- `src/metadarkmatter/core/classification/thresholds.py` (legacy cascade)
- `src/metadarkmatter/core/classification/classifiers/base.py` (`ANIWeightedClassifier`, per-read programmatic API)

Key optimization:
- ANI matrix loaded once into a float64 NumPy array; lookups are polars joins against a pre-built symmetric pair table
- The default in-memory mode loads the whole alignment file; `--streaming --chunk-size N` processes reads in batches with bounded memory (a read's hits must be contiguous, as BLAST/MMseqs2 emit them)
- All per-read operations are vectorized in polars (no Python loop over reads)

---

## Summary

**The algorithm in one sentence:**

For each read, find the closest reference genome and competing genomes, then use ANI (or AAI) between these genomes to determine if the read's divergence represents novel diversity or just a conserved gene region.

**Key insight:**

High sequence identity alone doesn't guarantee novel diversity - a read might match multiple divergent genomes equally well (conserved region). The ANI-weighted uncertainty metric distinguishes confident placements from ambiguous ones.
