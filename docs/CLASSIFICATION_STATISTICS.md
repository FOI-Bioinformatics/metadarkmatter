# ANI-Weighted Placement Uncertainty: Statistical Framework

This document describes the statistical metrics and classification thresholds used by metadarkmatter to distinguish novel bacterial diversity from known taxa in environmental DNA samples.

**For comprehensive methods documentation with all calculations at scientific standard, see [METHODS.md](METHODS.md).**

**For step-by-step algorithm explanation with ANI/AAI matrix usage, see [ALGORITHM_DETAILED.md](ALGORITHM_DETAILED.md).**

**Scope of the thresholds.** The thresholds in this document define the hard-threshold cascade that is reported as `legacy_call`. The primary `taxonomic_call` is the maximum a posteriori category of a six-category Bayesian model whose Gaussians are derived from these thresholds; its decision boundaries are close to, but not the same as, the thresholds (for example Known -> Novel Species at N ~ 6.2% rather than 4%, and no use of the inferred uncertainty for single-hit reads). See [METHODS.md](METHODS.md) sections 2.5 and 10.

## Core Metrics

### Novelty Index (N)

**Definition:** Measures divergence from the closest known reference genome.

```
N = 100 - pident
```

Where `pident` is the percent identity of the best hit (highest coverage-weighted bitscore; ties broken by higher `pident`, then accession).

| Novelty Index | BLAST Identity | Interpretation |
|---------------|----------------|----------------|
| N < 4% | pident > 96% | At or above species boundary |
| 4% <= N < 20% | 80% < pident <= 96% | Below species boundary |
| 20% <= N <= 25% | 75% <= pident <= 80% | Genus-level divergence |
| N > 25% | pident < 75% | Very high divergence |

### Placement Uncertainty (U)

**Definition:** Measures placement ambiguity when a read matches multiple genomes with similar scores.

```
U = 100 - ANI(G_top, G_2)                                   (default, --uncertainty-mode second)
U = 100 - max(ANI(G_top, G_i)) for all G_i in the set      (--uncertainty-mode max)
```

Where:
- `G_top` = genome of the best hit
- Ambiguous hit set = hits with coverage-weighted bitscore >= 95% of the top weighted bitscore; genomes are collapsed to their species representatives
- `G_2` = the highest-ranked competing representative (same ordering as best-hit selection)
- ANI = Average Nucleotide Identity between genome pairs
- U = 0 when there is no competing representative (single-hit reads)

| Uncertainty | ANI Between Competitors | Interpretation |
|-------------|------------------------|----------------|
| U < 1.5% | ANI > 98.5% | Same species, confident placement |
| 1.5% <= U < 5% | 95% < ANI <= 98.5% | Species boundary zone, ambiguous |
| U >= 5% | ANI <= 95% | Different species, conserved gene |

## Classification Space

The two metrics create a 2D classification space:

```
Uncertainty (U)
     ^
     |
 5%+ |  [Ambiguous, or Conserved Region when hits span genera]
     |
1.5-5|  [Species Boundary]   [Species Boundary]   [Species Boundary]
     |
 <1.5%|  [Known Species]      [Novel Species]      [Novel Genus]      [Unclassified]
     |
     +------------------------------------------------------------------> Novelty (N)
         0%    4%           20%           25%
```

This is the layout of `legacy_call`. In addition, reads with an identity gap < 2% and more than one competing hit, and single-hit reads with N >= 7%, are Ambiguous in `legacy_call` regardless of their position in the diagram.

## Classification Categories

### Known Species
- **Criteria:** N < 4% AND U < 1.5%
- **Meaning:** Read matches a known species (>96% identity) with confident placement (competitors share >98.5% ANI)
- **Biological interpretation:** The read originates from a characterized species in the reference database

### Novel Species
- **Criteria:** 4% <= N < 20% AND U < 1.5%
- **Meaning:** Moderate divergence from reference (80-96% identity) with confident placement
- **Biological interpretation:** Potential novel species - significant divergence from known taxa but placement is unambiguous

### Novel Genus
- **Criteria:** 20% <= N <= 25% AND U < 1.5%
- **Meaning:** High divergence (75-80% identity) with confident placement
- **Biological interpretation:** Potential novel genus - substantial divergence suggesting higher taxonomic novelty

### Species Boundary
- **Criteria:** 1.5% <= U < 5% (regardless of N)
- **Meaning:** Read matches multiple genomes that share 95-98.5% ANI (species boundary zone)
- **Biological interpretation:** Cannot confidently assign to a single species; may be shared gene or recent divergence

### Ambiguous
- **Criteria:** U >= 5% with competitors in the same genus; or identity gap < 2% with more than one competing hit; or (legacy only) a single hit with N >= 7%
- **Meaning:** Read matches multiple genomes that share <95% ANI, or competing hits cannot be separated by identity
- **Biological interpretation:** Placement cannot be resolved from the available references

### Conserved Region
- **Criteria:** U >= 5% AND the competing hits span genera (`ambiguity_scope = across_genera`)
- **Meaning:** Read matches genomes from different genera with similar scores
- **Biological interpretation:** Likely a conserved gene (16S rRNA, housekeeping genes) present across distantly related genomes

### Ambiguous Within Genus
- **Criteria:** (Bayesian `taxonomic_call` only) MAP category Novel Genus with genus uncertainty >= 10% and at least one competing genome
- **Meaning:** Genus-level divergence, but the competing hits are to species with low mutual ANI

### Unclassified
- **Criteria:** Does not fit any above category
- **Meaning:** Edge cases, typically very high novelty (N > 25%) or unusual metric combinations
- **Biological interpretation:** Requires manual review; may represent contamination, chimeric reads, or extremely novel diversity

## Decision Tree

Legacy cascade (`legacy_call`), rules in evaluation order:

```
START: Read with BLAST hits
         |
         v
    Calculate N, U, identity gap
         |
         v
   Identity gap < 2% and > 1 hit? --YES--> AMBIGUOUS
         |
        NO
         v
   Is U >= 5%? ----YES----> CONSERVED REGION (hits span genera) or AMBIGUOUS
         |
        NO
         v
   Is 1.5% <= U < 5%? --YES--> SPECIES BOUNDARY
         |
        NO (U < 1.5%)
         v
   Is N < 4%? ----YES----> KNOWN SPECIES
         |
        NO
         v
   Single hit and N >= 7%? --YES--> AMBIGUOUS
         |
        NO
         v
   Is 4% <= N < 20%? --YES--> NOVEL SPECIES
         |
        NO
         v
   Is 20% <= N <= 25% (and genus uncertainty < 10%)? --YES--> NOVEL GENUS
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

When competing genomes share >98.5% ANI:
- They represent the same species (potentially different strains)
- Placement ambiguity reflects strain-level variation, not taxonomic uncertainty
- The read can be confidently assigned to the species level

### Genus Boundary (~70-80% ANI)

The genus boundary is less precisely defined than species:
- **Qin et al. 2014**: Found genus boundaries vary from 65-83% ANI
- We use 75% (N=25%) as a conservative genus boundary
- Combined with 80% (N=20%) for novel species maximum

## GTDB Method Comparison

metadarkmatter's methods are designed to align with the Genome Taxonomy Database (GTDB) while extending the approach to read-level classification.

### Alignment with GTDB Standards

| Aspect | GTDB | metadarkmatter | Compatibility |
|--------|------|----------------|---------------|
| Species boundary | 95% ANI | 96% pident (N < 4%) | **Aligned** |
| Method | Genome-to-genome (FastANI) | Read-to-genome (BLAST pident) | Different scope |
| Alignment fraction | Requires AF >= 50% | Configurable (default 30%) | **Optional filter** |
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

| Preset | Known Species (N <) | Novel Species N range | Novel Genus N range | Alignment Fraction | Other |
|--------|------|------|------|------|------|
| `default` | 4% (96% identity) | 4-20% | 20-25% | >= 30% | Package defaults |
| `gtdb-strict` | 5% (95%) | 5-15% | 15-22% | >= 50% | min length 100 bp |
| `gtdb-relaxed` | 3% (97%) | 3-15% | 15-25% | >= 30% | min length 50 bp |
| `conservative` | 4% (96%) | 4-15% | 15-25% | >= 30% | min length 100 bp |
| `literature-strict` | 4% (96%) | 4-15% | 15-22% | >= 50% | `confidence_threshold = 60` |
| `coverage-strict` | 4% (96%) | 4-15% | 15-25% | >= 30% | sigmoid coverage weighting, strength 0.7 |
| `adaptive` | 4% (96%) | 4-20% | 20-25% | >= 30% | intended for use with `--adaptive-thresholds` |

Presets change the thresholds and therefore both `legacy_call` and the Bayesian Gaussian parameters. Uncertainty thresholds are the defaults (1.5% / 5%) in every preset. Presets do not apply in protein mode, which uses fixed protein thresholds (`--alignment-mode protein`). Flags passed explicitly on the command line take precedence over the preset.

### Literature-Strict Preset

The `literature-strict` preset implements thresholds derived from a comprehensive literature review:

```bash
metadarkmatter score classify --alignment results.tsv --ani ani.csv \
    --preset literature-strict --output classifications.csv
```

**Key parameters:**
- **96% ANI species boundary** (4% novelty): Stricter than GTDB's 95%
- **1.5% uncertainty threshold**: Requires 98.5% ANI between competing genomes
- **60% confidence threshold**: `ScoringConfig.confidence_threshold`; this field is currently not used by the CLI classifier (`low_confidence` is based on posterior entropy), so it has no effect on the output
- **50% alignment fraction**: GTDB-compatible alignment quality requirement
- **15% novelty max for novel species** (85% identity)
- **22% novelty max for genus**: ~78% identity, more conservative genus boundary

**Literature basis:**
- Jain et al. 2018: 95-96% ANI species boundary validated on 90K genomes
- Riesco & Trujillo 2024: Modern species delineation standards
- Parks et al. 2020: GTDB alignment fraction requirements

## Additional Quality Metrics

### Confidence Score (0-100)

`confidence_score` is derived from the posterior entropy H of the Bayesian model:

```
confidence_score = (1 - H / log2(6)) * 100
```

H = 0 (all posterior mass on one category) gives 100; a uniform posterior over the six categories (H = log2(6) ~ 2.585) gives 0. Reads with H > 1.5 (confidence_score below ~42) are flagged `low_confidence`. Off-target reads have a fixed confidence_score of 10.0.

The score measures how clearly one category dominates the posterior, not the probability that the call is correct; the posteriors are not calibrated against labelled data ([STATISTICAL_ASSUMPTIONS.md](STATISTICAL_ASSUMPTIONS.md), sections 4 and 10). A read deep inside the Ambiguous Gaussian can have confidence_score = 100 while its taxonomic placement is unresolved.

**Example values** (default configuration, recomputed with the current code):
- N = 1.5%, U = 0.8%, 3 hits, identity gap 2.5%: Known Species, confidence_score 76.3
- N = 8.8%, U = 1.2%, 4 hits, identity gap 3.0%: Novel Species, confidence_score 49.3
- N = 3.8%, U = 1.4%, identity gap 5%: Known Species, confidence_score 57.0 (near the Known/Species Boundary/Novel Species overlap)

An earlier multi-factor score (margin from thresholds + hit count/identity gap + identity) is still defined in `core/constants.py` (`calculate_confidence_score`) but is not used to produce output ([METHODS.md](METHODS.md) section 3.1).

### Phylogenetic Context Metrics

**Genus uncertainty**: Based on minimum ANI to secondary genomes
- Low value: All ambiguous hits within same genus
- High value: Hits span multiple genera (potential conserved gene)

**Ambiguity scope**: Categorizes the phylogenetic breadth of ambiguous hits
- `unambiguous`: Single genome hit
- `within_species`: All secondary hits share >= 95% ANI with best hit
- `within_genus`: All secondary hits share >= 80% ANI with best hit
- `across_genera`: Hits span multiple genera

**Number of competing genera**: Count of secondary genomes from different genera (computed internally; not written to the output file)

**Known issue:** `ambiguity_scope` compares hit counts with the number of distinct secondary representatives, so a competing genome with two hits in the bitscore window is labelled `across_genera` even at high ANI ([METHODS.md](METHODS.md) section 2.2).

### Example Output with New Metrics

Selected columns (the full output also contains `second_hit_identity`, `identity_gap`, `diversity_status`, `low_confidence`, the six posterior columns, `posterior_entropy`, `inferred_uncertainty` and `uncertainty_type`); values recomputed with the default configuration:

```
read_id,best_match_genome,top_hit_identity,novelty_index,placement_uncertainty,
genus_uncertainty,ambiguity_scope,num_ambiguous_hits,
confidence_score,taxonomic_call,legacy_call,is_novel

read001,GCF_000195955.2,98.5,1.5,0.8,0.8,within_species,3,76.3,Known Species,Known Species,false
read002,GCF_000008985.1,91.2,8.8,1.2,15.5,within_genus,4,49.3,Novel Species,Novel Species,true
read003,GCF_000123456.1,99.1,0.9,22.3,28.3,across_genera,6,100.0,Conserved Region,Ambiguous,false
```

(read001 and read002 assume identity gaps of 2.5% and 3.0%; read003 has an identity gap of 0.2%, which makes `legacy_call` Ambiguous before the Conserved Region rule is reached.)

## Worked Examples

Calls were recomputed with the default configuration (`taxonomic_call` from the Bayesian model, `legacy_call` from the cascade). "Second hit" is the best hit to a different genome, within 95% of the top bitscore.

### Example 1: Known Bacterium
```
Read: ERR123456.1
Best hit: GCF_000123456 at 97.5% identity
Second hit: GCF_000789012 at 95.3% identity (bitscore within 95%)
ANI(GCF_000123456, GCF_000789012) = 99.2%

N = 100 - 97.5 = 2.5% (< 4%)
U = 100 - 99.2 = 0.8% (< 1.5%)
Identity gap = 2.2%

taxonomic_call: KNOWN SPECIES (p_known = 0.863, confidence_score = 71.9)
legacy_call:    KNOWN SPECIES
Interpretation: Matches two strains of the same species
```

With a second hit at 96.8% (gap 0.7%), `legacy_call` would be Ambiguous (identity-gap rule) while `taxonomic_call` stays Known Species.

### Example 2: Potential Novel Species
```
Read: ERR789012.1
Best hit: GCF_000234567 at 91.2% identity
Second hit: GCF_000345678 at 89.0% identity
ANI(GCF_000234567, GCF_000345678) = 99.0%

N = 100 - 91.2 = 8.8% (4-20%)
U = 100 - 99.0 = 1.0% (< 1.5%)
Identity gap = 2.2%

taxonomic_call: NOVEL SPECIES (p_novel_species = 0.692, p_species_boundary = 0.213,
                confidence_score = 51.3)
legacy_call:    NOVEL SPECIES
Interpretation: Divergent from the known references, but placement is confident
```

### Example 3: Conserved Gene
```
Read: ERR345678.1
Best hit: GCF_000456789 at 99.1% identity
Second hit: GCF_000567890 at 98.9% identity
ANI(GCF_000456789, GCF_000567890) = 77.7%  (< 80%: different genera)

N = 100 - 99.1 = 0.9% (< 4%)
U = 100 - 77.7 = 22.3% (>= 5%)
Identity gap = 0.2%

taxonomic_call: CONSERVED REGION (Stage 1 MAP Ambiguous, p_ambiguous = 1.00,
                hits across genera)
legacy_call:    AMBIGUOUS (identity-gap rule is evaluated first)
Interpretation: Likely 16S rRNA or housekeeping gene shared across genera
```

If the competing genomes shared 80-95% ANI (same genus by the 80% convention), both calls would be Ambiguous rather than Conserved Region.

### Example 4: Species Boundary
```
Read: ERR456789.1
Best hit: GCF_000678901 at 94.2% identity
Second hit: GCF_000789012 at 91.8% identity
ANI(GCF_000678901, GCF_000789012) = 96.5%

N = 100 - 94.2 = 5.8% (4-20%)
U = 100 - 96.5 = 3.5% (1.5-5%)
Identity gap = 2.4%

taxonomic_call: SPECIES BOUNDARY (p_species_boundary = 0.689, confidence_score = 46.5)
legacy_call:    SPECIES BOUNDARY
Interpretation: Competing genomes are at the species boundary; cannot confidently place
```

With a second hit at 93.8% identity (gap 0.4%), `legacy_call` would be Ambiguous, `taxonomic_call` would remain Species Boundary, and posterior_entropy (1.51) would exceed 1.5, setting `low_confidence`.

## Summary Table

Hard thresholds (`legacy_call`, nucleotide defaults):

| Category | Novelty (N) | Uncertainty (U) | BLAST Identity | ANI Range |
|----------|-------------|-----------------|----------------|-----------|
| Known Species | < 4% | < 1.5% | > 96% | > 98.5% |
| Novel Species | 4% to <20% | < 1.5% | 80-96% | > 98.5% |
| Novel Genus | 20-25% | < 1.5% | 75-80% | > 98.5% |
| Species Boundary | any | 1.5% to <5% | any | >95% to 98.5% |
| Ambiguous | any | >= 5% (same genus) | any | <= 95% |
| Conserved Region | any | >= 5% (across genera) | any | <= 95% |
| Unclassified | > 25% or edge | < 1.5% | < 75% | > 98.5% |

Additional legacy rules: identity gap < 2% with several hits -> Ambiguous; single hit with N >= 7% -> Ambiguous. The corresponding `taxonomic_call` boundaries are given in [METHODS.md](METHODS.md) section 2.5.

## References

1. Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S. (2018). High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. *Nature Communications* 9:5114.

2. Goris J, Konstantinidis KT, Klappenbach JA, Coenye T, Vandamme P, Tiedje JM. (2007). DNA-DNA hybridization values and their relationship to whole-genome sequence similarities. *International Journal of Systematic and Evolutionary Microbiology* 57:81-91.

3. Konstantinidis KT, Tiedje JM. (2005). Genomic insights that advance the species definition for prokaryotes. *Proceedings of the National Academy of Sciences* 102:2567-2572.

4. Qin QL, Xie BB, Zhang XY, Chen XL, Zhou BC, Zhou J, Oren A, Zhang YZ. (2014). A proposed genus boundary for the prokaryotes based on genomic insights. *Journal of Bacteriology* 196:2210-2215.

5. Parks DH, Chuvochina M, Waite DW, Rinke C, Skarshewski A, Chaumeil PA, Hugenholtz P. (2018). A standardized bacterial taxonomy based on genome phylogeny substantially revises the tree of life. *Nature Biotechnology* 36:996-1004.

6. Riesco R, Trujillo ME. (2024). Update on the proposed minimal standards for the use of genome data for the taxonomy of prokaryotes. *International Journal of Systematic and Evolutionary Microbiology* 74:006300.
