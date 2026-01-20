# ANI-Weighted Placement Uncertainty: Statistical Framework

This document describes the statistical metrics and classification thresholds used by metadarkmatter to distinguish novel bacterial diversity from known taxa in environmental DNA samples.

**For step-by-step algorithm explanation with ANI/AAI matrix usage, see [ALGORITHM_DETAILED.md](ALGORITHM_DETAILED.md).**

## Core Metrics

### Novelty Index (N)

**Definition:** Measures divergence from the closest known reference genome.

```
N = 100 - pident
```

Where `pident` is the percent identity of the best BLAST hit (highest bitscore).

| Novelty Index | BLAST Identity | Interpretation |
|---------------|----------------|----------------|
| N < 5% | pident > 95% | At or above species boundary |
| 5% <= N < 15% | 85% < pident <= 95% | Below species boundary |
| 15% <= N <= 25% | 75% <= pident <= 85% | Genus-level divergence |
| N > 25% | pident < 75% | Very high divergence |

### Placement Uncertainty (U)

**Definition:** Measures placement ambiguity when a read matches multiple genomes with similar scores.

```
U = 100 - max(ANI(G_top, G_i)) for all G_i in ambiguous hit set
```

Where:
- `G_top` = genome with highest bitscore hit
- Ambiguous hit set = genomes with bitscore >= 95% of top bitscore
- ANI = Average Nucleotide Identity between genome pairs

| Uncertainty | ANI Between Competitors | Interpretation |
|-------------|------------------------|----------------|
| U < 2% | ANI > 98% | Same species, confident placement |
| 2% <= U < 5% | 95% < ANI <= 98% | Species boundary zone, ambiguous |
| U >= 5% | ANI <= 95% | Different species, conserved gene |

## Classification Space

The two metrics create a 2D classification space:

```
Uncertainty (U)
     ^
     |
 5%+ |  [Conserved Region]   [Conserved Region]   [Conserved Region]
     |
 2-5%|  [Ambiguous]          [Ambiguous]          [Ambiguous]
     |
 <2% |  [Known Species]      [Novel Species]      [Novel Genus]
     |
     +------------------------------------------------------> Novelty (N)
         0%    5%           15%           25%
```

## Classification Categories

### Known Species
- **Criteria:** N < 5% AND U < 2%
- **Meaning:** Read matches a known species (>95% identity) with confident placement (competitors share >98% ANI)
- **Biological interpretation:** The read originates from a characterized species in the reference database

### Novel Species
- **Criteria:** 5% <= N < 15% AND U < 2%
- **Meaning:** Moderate divergence from reference (85-95% identity) with confident placement
- **Biological interpretation:** Potential novel species - significant divergence from known taxa but placement is unambiguous

### Novel Genus
- **Criteria:** 15% <= N <= 25% AND U < 2%
- **Meaning:** High divergence (75-85% identity) with confident placement
- **Biological interpretation:** Potential novel genus - substantial divergence suggesting higher taxonomic novelty

### Ambiguous
- **Criteria:** 2% <= U < 5% (regardless of N)
- **Meaning:** Read matches multiple genomes that share 95-98% ANI (species boundary zone)
- **Biological interpretation:** Cannot confidently assign to a single species; may be shared gene or recent divergence

### Conserved Region
- **Criteria:** U >= 5% (regardless of N)
- **Meaning:** Read matches multiple genomes that share <95% ANI (different species)
- **Biological interpretation:** Likely a conserved gene (16S rRNA, housekeeping genes) present across distantly related genomes

### Unclassified
- **Criteria:** Does not fit any above category
- **Meaning:** Edge cases, typically very high novelty (N > 25%) or unusual metric combinations
- **Biological interpretation:** Requires manual review; may represent contamination, chimeric reads, or extremely novel diversity

## Decision Tree

```
START: Read with BLAST hits
         |
         v
    Calculate N and U
         |
         v
   Is U >= 5%? ----YES----> CONSERVED REGION
         |
        NO
         v
   Is 2% <= U < 5%? --YES--> AMBIGUOUS
         |
        NO (U < 2%)
         v
   Is N < 5%? ----YES----> KNOWN SPECIES
         |
        NO
         v
   Is 5% <= N < 15%? --YES--> NOVEL SPECIES
         |
        NO
         v
   Is 15% <= N <= 25%? --YES--> NOVEL GENUS
         |
        NO
         v
   UNCLASSIFIED
```

## Threshold Justification

### Species Boundary (95-96% ANI)

The 95-96% ANI threshold for prokaryotic species is well-established:

- **Jain et al. 2018** (Nature Communications): Analyzed 90,000+ genomes, found 95-96% ANI corresponds to 70% DNA-DNA hybridization (classical species standard)
- **Goris et al. 2007** (PNAS): Demonstrated ANI >95% correlates with >70% DDH for diverse bacterial groups
- **Konstantinidis & Tiedje 2005** (PNAS): Proposed ANI as gold standard for species delineation

### Confident Placement (>98% ANI)

When competing genomes share >98% ANI:
- They represent the same species (potentially different strains)
- Placement ambiguity reflects strain-level variation, not taxonomic uncertainty
- The read can be confidently assigned to the species level

### Genus Boundary (~70-80% ANI)

The genus boundary is less precisely defined than species:
- **Qin et al. 2014**: Found genus boundaries vary from 65-83% ANI
- We use 75% (N=25%) as a conservative genus boundary
- Combined with 85% (N=15%) for novel species maximum

## GTDB Method Comparison

metadarkmatter's methods are designed to align with the Genome Taxonomy Database (GTDB) while extending the approach to read-level classification.

### Alignment with GTDB Standards

| Aspect | GTDB | metadarkmatter | Compatibility |
|--------|------|----------------|---------------|
| Species boundary | 95% ANI | 95% pident (N < 5%) | **Aligned** |
| Method | Genome-to-genome (FastANI) | Read-to-genome (BLAST pident) | Different scope |
| Alignment fraction | Requires AF >= 50% | Configurable (default 0%) | **Optional filter** |
| Rank normalization | RED (phylogenetic) | ANI-based clustering | Similar concept |
| Higher ranks | Phylogenetic trees | ANI-based genus clustering | Simplified approach |

### Key Differences

1. **Scope of comparison**
   - GTDB: Whole-genome ANI between assembled genomes
   - metadarkmatter: Per-read BLAST identity to reference genomes
   - Implication: Read-level pident is a proxy for genome-level ANI

2. **Phylogenetic context**
   - GTDB: Uses RED (Relative Evolutionary Divergence) from phylogenetic trees
   - metadarkmatter: Uses ANI-based genus clustering without trees
   - Both approaches group genomes by evolutionary relatedness

3. **Novel taxa handling**
   - GTDB: Assigns alphanumeric placeholders (e.g., "Firmicutes_A")
   - metadarkmatter: Generates deterministic names (e.g., "Francisella sp. MDM-A7X")

### GTDB-Compatible Presets

Use `--preset gtdb-strict` for maximum GTDB compatibility:
```bash
metadarkmatter score classify --alignment results.tsv --ani ani.csv \
    --preset gtdb-strict --output classifications.csv
```

| Preset | Species Boundary | Alignment Fraction | Use Case |
|--------|------------------|-------------------|----------|
| `gtdb-strict` | 95% ANI | >= 50% | Publication-quality, GTDB-compatible |
| `gtdb-relaxed` | 97% ANI | >= 30% | Exploratory, allows close relatives |
| `conservative` | 96% ANI | None | Balanced, stricter uncertainty |
| `literature-strict` | 96% ANI | >= 50% | Literature-backed, high confidence required |
| `default` | 95% ANI | None | Standard metadarkmatter defaults |

### Literature-Strict Preset

The `literature-strict` preset implements thresholds derived from a comprehensive literature review:

```bash
metadarkmatter score classify --alignment results.tsv --ani ani.csv \
    --preset literature-strict --output classifications.csv
```

**Key parameters:**
- **96% ANI species boundary** (4% novelty): Stricter than GTDB's 95%
- **1.5% uncertainty threshold**: Requires 98.5% ANI between competing genomes
- **60% confidence threshold**: Higher confidence required for classification
- **50% alignment fraction**: GTDB-compatible alignment quality requirement
- **22% novelty max for genus**: ~78% identity, more conservative genus boundary

**Literature basis:**
- Jain et al. 2018: 95-96% ANI species boundary validated on 90K genomes
- Riesco & Trujillo 2024: Modern species delineation standards
- Parks et al. 2020: GTDB alignment fraction requirements

## Additional Quality Metrics

### Confidence Score (0-100)

The confidence score integrates multiple quality factors to quantify classification reliability without changing the classification category itself.

```
confidence_score = margin_score + placement_certainty + alignment_quality
```

**Components:**

1. **Margin from threshold boundaries (0-40 pts)**: The further a classification falls from decision thresholds, the more confident the call.
   - Known Species: Distance from 5% novelty boundary
   - Novel Species: Minimum distance from both 5% and 20% novelty boundaries
   - Novel Genus: Minimum distance from both 20% and 25% boundaries
   - Ambiguous/Unclassified: Base score of 10 pts

2. **Placement certainty (0-40 pts)**: Based on hit count and identity gap to secondary genomes
   - Hit count: 1 hit = 20 pts, 2-3 hits = 15 pts, 4-5 hits = 10 pts, 6-10 hits = 5 pts
   - Identity gap: >= 5% = 20 pts, >= 2% = 15 pts, >= 1% = 10 pts, >= 0.5% = 5 pts

3. **Alignment quality proxy (0-20 pts)**: Derived from top hit identity
   - Scaled linearly from 70% (0 pts) to 95%+ (20 pts)

| Score Range | Interpretation | Action |
|-------------|----------------|--------|
| 80-100 | High confidence | Suitable for publication |
| 60-79 | Moderate confidence | Review recommended |
| 40-59 | Low confidence | Manual verification needed |
| < 40 | Very low confidence | Treat with caution |

**Example confidence scenarios:**
- Read at 98% identity, single hit, 2% novelty margin: ~90 (high)
- Read at 88% identity, 3 hits, 5% from boundaries: ~75 (moderate)
- Read at 95.1% identity, 5 hits, 0.1% from boundary: ~35 (low, borderline)

### Phylogenetic Context Metrics

**Genus uncertainty**: Based on minimum ANI to secondary genomes
- Low value: All ambiguous hits within same genus
- High value: Hits span multiple genera (potential conserved gene)

**Ambiguity scope**: Categorizes the phylogenetic breadth of ambiguous hits
- `unambiguous`: Single genome hit
- `within_species`: All secondary hits share >= 95% ANI with best hit
- `within_genus`: All secondary hits share >= 80% ANI with best hit
- `across_genera`: Hits span multiple genera

**Number of competing genera**: Count of secondary genomes from different genera

### Example Output with New Metrics

```
read_id,best_match_genome,top_hit_identity,novelty_index,placement_uncertainty,
genus_uncertainty,ambiguity_scope,num_ambiguous_hits,num_competing_genera,
confidence_score,taxonomic_call,is_novel

read001,GCF_000195955.2,98.5,1.5,0.8,0.8,within_species,3,0,87.2,Known Species,false
read002,GCF_000008985.1,91.2,8.8,1.2,15.5,within_genus,4,0,72.1,Novel Species,true
read003,GCF_000123456.1,99.1,0.9,17.7,28.3,across_genera,6,2,45.8,Conserved Region,false
```

## Worked Examples

### Example 1: Known Bacterium
```
Read: ERR123456.1
Best hit: GCF_000123456 at 97.5% identity
Second hit: GCF_000789012 at 96.8% identity (bitscore within 95%)
ANI(GCF_000123456, GCF_000789012) = 99.2%

N = 100 - 97.5 = 2.5% (< 5%)
U = 100 - 99.2 = 0.8% (< 2%)

Classification: KNOWN SPECIES
Interpretation: Matches E. coli strains K-12 and BL21
```

### Example 2: Potential Novel Species
```
Read: ERR789012.1
Best hit: GCF_000234567 at 91.2% identity
Second hit: GCF_000345678 at 90.5% identity
ANI(GCF_000234567, GCF_000345678) = 98.5%

N = 100 - 91.2 = 8.8% (5-15%)
U = 100 - 98.5 = 1.5% (< 2%)

Classification: NOVEL SPECIES
Interpretation: Divergent from known Pseudomonas, but placement is confident
```

### Example 3: Conserved Gene
```
Read: ERR345678.1
Best hit: GCF_000456789 at 99.1% identity
Second hit: GCF_000567890 at 98.9% identity
ANI(GCF_000456789, GCF_000567890) = 82.3%

N = 100 - 99.1 = 0.9% (< 5%)
U = 100 - 82.3 = 17.7% (>= 5%)

Classification: CONSERVED REGION
Interpretation: Likely 16S rRNA or housekeeping gene shared across genera
```

### Example 4: Ambiguous Placement
```
Read: ERR456789.1
Best hit: GCF_000678901 at 94.2% identity
Second hit: GCF_000789012 at 93.8% identity
ANI(GCF_000678901, GCF_000789012) = 96.5%

N = 100 - 94.2 = 5.8% (5-15%)
U = 100 - 96.5 = 3.5% (2-5%)

Classification: AMBIGUOUS
Interpretation: Competing genomes are at species boundary; cannot confidently place
```

## Summary Table

| Category | Novelty (N) | Uncertainty (U) | BLAST Identity | ANI Range |
|----------|-------------|-----------------|----------------|-----------|
| Known Species | < 5% | < 2% | > 95% | > 98% |
| Novel Species | 5% to <15% | < 2% | 85-95% | > 98% |
| Novel Genus | 15-25% | < 2% | 75-85% | > 98% |
| Ambiguous | any | 2% to <5% | any | >95% to 98% |
| Conserved Region | any | >= 5% | any | <= 95% |
| Unclassified | > 25% or edge | < 2% | < 75% | > 98% |

## References

1. Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S. (2018). High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. *Nature Communications* 9:5114.

2. Goris J, Konstantinidis KT, Klappenbach JA, Coenye T, Vandamme P, Tiedje JM. (2007). DNA-DNA hybridization values and their relationship to whole-genome sequence similarities. *International Journal of Systematic and Evolutionary Microbiology* 57:81-91.

3. Konstantinidis KT, Tiedje JM. (2005). Genomic insights that advance the species definition for prokaryotes. *Proceedings of the National Academy of Sciences* 102:2567-2572.

4. Qin QL, Xie BB, Zhang XY, Chen XL, Zhou BC, Zhou J, Oren A, Zhang YZ. (2014). A proposed genus boundary for the prokaryotes based on genomic insights. *Journal of Bacteriology* 196:2210-2215.

5. Parks DH, Chuvochina M, Waite DW, Rinke C, Skarshewski A, Chaumeil PA, Hugenholtz P. (2018). A standardized bacterial taxonomy based on genome phylogeny substantially revises the tree of life. *Nature Biotechnology* 36:996-1004.

6. Riesco R, Trujillo ME. (2024). Update on the proposed minimal standards for the use of genome data for the taxonomy of prokaryotes. *International Journal of Systematic and Evolutionary Microbiology* 74:006300.
