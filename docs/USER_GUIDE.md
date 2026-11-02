# Metadarkmatter User Guide

Practical guide for using metadarkmatter to detect novel microbial diversity.

## Prerequisites

- Familiarity with command-line interfaces
- Basic understanding of BLAST output formats
- Environmental DNA metagenomic sequencing data

## Documentation Map

For comprehensive coverage, see specialized documentation:

- **Installation & Setup**: [README.md](../README.md)
- **Hands-on Tutorial**: [TUTORIAL_ENVIRONMENTAL_SPECIES.md](TUTORIAL_ENVIRONMENTAL_SPECIES.md)
- **Complete Command Reference**: [CLI_REFERENCE.md](CLI_REFERENCE.md)
- **Performance Optimization**: [PERFORMANCE.md](PERFORMANCE.md)
- **Troubleshooting**: [TROUBLESHOOTING.md](TROUBLESHOOTING.md)
- **Workflow Strategy**: [WORKFLOW.md](WORKFLOW.md)
- **Algorithm Details**: [REFERENCE.md](REFERENCE.md)

---

## Quick Start

For the complete end-to-end workflow (genome download through report
generation), see the [Tutorial](TUTORIAL_ENVIRONMENTAL_SPECIES.md).

If you already have alignment results and an ANI matrix, the core
classification step is:

    metadarkmatter score classify \
        --alignment results.blast.tsv.gz \
        --ani ani_matrix.csv \
        --metadata genome_metadata.tsv \
        --output classifications.csv

For importing external BLAST/MMseqs2 results, see the
[Workflow Guide](WORKFLOW.md#importing-external-alignment-results).

---

## Input File Formats

### BLAST Tabular Format

Metadarkmatter accepts BLAST tabular output (`-outfmt 6`) with 12, 13 or 14 tab-separated columns:

```
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore [qlen] [qcov]
```

`blast align` and `mmseqs2 search` write the 13-column form (with `qlen`),
which lets the alignment-fraction filter use the true read length.
Reverse-strand coordinates (`qstart > qend`, as emitted by MMseqs2 and
DIAMOND) are normalized on input.

**Column requirements:**
- `qseqid`: Read identifier
- `sseqid`: Reference genome identifier (format: `{accession}|{contig_id}`)
- `pident`: Percent identity (0-100)
- `bitscore`: BLAST bitscore (used for top hit selection)

#### Generating BLAST Results

Use `metadarkmatter blast align` which handles:
- FASTQ/FASTA input (automatic conversion)
- Competitive recruitment parameters
- Proper output format (13 columns; gzipped when `--output` ends in `.gz`)

```bash
metadarkmatter blast align \
    --query reads.fastq.gz \
    --database blastdb/pangenome \
    --output sample.blast.tsv.gz
```

Alternatively, use BLAST directly:
```bash
blastn -query reads.fasta \
       -db blastdb/pangenome \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" \
       -out sample.blast.tsv \
       -max_target_seqs 500
```

If the database was not built with `metadarkmatter blast makedb`, subject IDs
will not be in `{accession}|{contig_id}` form; classify such results with
`--id-mapping` (see the Workflow Guide link above).

#### Compressed Files

BLAST output can be gzip-compressed (`.tsv.gz`). This is recommended for large files:
- 10M alignments: ~500 MB uncompressed → ~50 MB compressed
- Decompression is automatic for files ending in `.gz`

### ANI Matrix Format

Square CSV (or TSV) matrix of pairwise ANI values, as written by `ani compute`:

```csv
genome,GCF_000195955.1,GCF_000242755.1
GCF_000195955.1,100.00,92.15
GCF_000242755.1,92.15,100.00
```

**Requirements:**
- First column holds row labels (`genome`); header row holds the same accessions
- Labels must match the genome accessions in the alignment (`{accession}|...` prefix)
- ANI values: 0-100 (percent). Matrices on a 0-1 scale are rejected with an error
- Missing pairs (value 0) are treated as 70% ANI (unrelated) for lookups
- Pairs present in only one direction are mirrored; forward/reverse values
  differing by more than 0.5 produce a warning (an error with `--strict-ani`)
- Diagonal values should be 100.0

Generate with:
```bash
metadarkmatter ani compute \
    --genomes genomes/ \
    --output ani_matrix.csv \
    --threads 16
```

### Genome Metadata Format

TSV file with genome taxonomic information (created automatically by `download genomes`):

```tsv
accession	species	genus	family	representative	gtdb_taxonomy
GCF_000195955.1	Francisella tularensis	Francisella	Francisellaceae	GCF_000195955.1	d__Bacteria;p__...
```

**Required columns:**
- `accession`: Genome accession (must match BLAST sseqid prefix)
- `species`: Species name
- `genus`: Genus name

**Optional columns:**
- `family`: Family name (used to infer `--target-family` for family validation)
- `representative`: Species representative accession (used to map hits to
  representatives for ANI lookups; if absent, every genome maps to itself)
- `gtdb_taxonomy`: Full GTDB taxonomy string

This file enables species-level tracking in classifications and reports.

---

## Output Interpretation

### Classification Results File

Main output contains one row per read. The most commonly used columns are:

| Column | Type | Description |
|--------|------|-------------|
| read_id | string | Original read identifier from the alignment |
| best_match_genome | string | Genome of the best hit (coverage-weighted bitscore) |
| top_hit_identity | float | Percent identity of the best hit (0-100) |
| novelty_index | float | 100 - top_hit_identity; measures divergence |
| placement_uncertainty | float | 100 - ANI(best genome, second-best competing genome) with the default `--uncertainty-mode second`; 0 for single-hit reads; null for Off-target reads |
| num_ambiguous_hits | int | Hits within 95% of the top bitscore (`--bitscore-threshold`) |
| taxonomic_call | string | Primary call: Bayesian MAP category (see below) |
| legacy_call | string | Call from the hard-threshold cascade (for comparison) |
| confidence_score | float | `(1 - posterior_entropy / log2(6)) * 100` (0-100); fixed at 10.0 for Off-target reads |
| is_novel | bool | True if `taxonomic_call` is Novel Species or Novel Genus |
| diversity_status | string | Known / Novel / Uncertain / Off-target |
| low_confidence | bool | True when posterior entropy exceeds 1.5 bits |

Further columns: `genus_uncertainty`, `ambiguity_scope`, `second_hit_identity`,
`identity_gap`, six posterior columns (`p_known_species`, `p_novel_species`,
`p_novel_genus`, `p_species_boundary`, `p_ambiguous`, `p_unclassified`),
`posterior_entropy`, `inferred_uncertainty` and `uncertainty_type`
(single-hit reads), and the family-validation columns
(`family_bitscore_ratio`, `family_identity_gap`, `in_family_hit_fraction`,
`external_best_genome`, `external_best_identity`, `in_family_identity`,
`in_family_novelty_index`).

**Optional columns** (when `--metadata` is provided, not in `--streaming` mode):
- `species`, `genus`, `family` of the best match genome

### Taxonomic Classification Categories

`taxonomic_call` is the maximum a posteriori category of a Bayesian
classifier with six 2D Gaussian categories in (novelty, uncertainty)
space: Known Species, Novel Species, Novel Genus, Species Boundary,
Ambiguous and Unclassified. A second stage refines Ambiguous reads whose
hits span genera to **Conserved Region**, and Novel Genus reads with high
genus uncertainty to **Ambiguous Within Genus**. With family validation
(inferred from `--metadata` or set with `--target-family`), reads whose
best in-family hit is clearly weaker than the best overall hit are
**Off-target**.

The criteria listed below are the hard thresholds of the `legacy_call`
cascade (nucleotide mode). They describe the intended regions, but the
Bayesian boundaries used for `taxonomic_call` are not identical: with the
default configuration and U = 0, the Known -> Novel Species transition
lies near N ~ 6.2%, Novel Species -> Novel Genus near N ~ 18.6%, and
Novel Genus -> Unclassified near N ~ 27.6%. For single-hit reads, the
inferred uncertainty affects only `legacy_call`.

#### 1. Known Species

**Legacy criteria:**
- Novelty Index < 4%
- Placement Uncertainty < 1.5%

**Interpretation:** Read maps with high identity and low ambiguity to a known reference genome. Likely originates from the matched species or closely related strain.

**Biological significance:** Known diversity, not a discovery target.

#### 2. Novel Species

**Legacy criteria:**
- Novelty Index: 4-20%
- Placement Uncertainty < 1.5%

**Interpretation:** Moderate divergence from closest reference but unambiguous placement to a single lineage. Suggests a novel species within a known genus.

**Biological significance:** Putative novel species. Recommend validation through:
- Whole-genome recovery (binning, assembly)
- 16S rRNA analysis
- Targeted cultivation

#### 3. Novel Genus

**Legacy criteria:**
- Novelty Index: 20-25%
- Placement Uncertainty < 1.5%

**Interpretation:** High divergence from all references but moderate-confidence placement. Suggests a novel genus.

**Biological significance:** Putative novel genus. Priority targets for:
- Single-cell genomics
- Cultivation efforts
- Phylogenetic characterization

#### 4. Species Boundary, Ambiguous, Conserved Region

**Legacy criteria:**
- Species Boundary: Placement Uncertainty 1.5-5%
- Ambiguous / Conserved Region: Placement Uncertainty >= 5% (Conserved Region when the competing hits span genera)

**Interpretation:** Read maps with similar scores to multiple divergent genomes. Often indicates a conserved genomic region (housekeeping genes, rRNA, mobile elements). Taxonomic placement is ambiguous.

**Biological significance:** Exclude from novel diversity estimates. Cannot be confidently assigned to specific lineage.

#### 5. Unclassified and Off-target

- **Unclassified**: novelty beyond the Novel Genus range (legacy: N > 25%).
- **Off-target**: best hit lies outside the target family (only with family validation).

### Threshold Interpretation

| Metric | Legacy threshold | Biological Meaning |
|--------|-----------|-------------------|
| Novelty < 4% | Known Species | >96% identity to reference |
| Novelty 4-20% | Novel Species | 80-96% identity (within genus) |
| Novelty 20-25% | Novel Genus | 75-80% identity (genus-level divergence) |
| Uncertainty < 1.5% | Confident | ANI >98.5% between competing genomes |
| Uncertainty >= 1.5% | Species Boundary / Ambiguous | ANI <98.5% between competing genomes |

**Note:** Thresholds are based on 95-96% ANI species boundary (literature-backed).

**For complete algorithm explanation including ANI/AAI matrix usage and decision tree logic, see [ALGORITHM_DETAILED.md](ALGORITHM_DETAILED.md).**

---

## Report Generation

### Single-Sample Report

Generate interactive HTML report:

```bash
metadarkmatter report generate \
    --classifications classifications.csv \
    --metadata genome_metadata.tsv \
    --output report.html
```

**Report tabs:**
- **Summary**: Metric cards, diversity bar, key findings, top species (with metadata)
- **Classification**: Novelty/uncertainty scatter and distributions, confidence analysis
- **Novel Diversity**: Novel clusters and their phylogenetic context (shown only when novel reads exist)
- **Reference**: Top-20 species bar chart (with metadata), ANI/AAI heatmaps (hierarchically clustered, with `--ani`/`--aai`), phylogenetic tree panel, recruitment plots
- **Data**: Interactive classification table

The default `--report-mode offline` produces a self-contained HTML file;
`--report-mode cdn` is smaller but needs network access to render.

### Multi-Sample Comparison

Compare multiple samples:

```bash
# Sample names are taken from the classification filenames in --input-dir
metadarkmatter report multi \
    --input-dir classifications/ \
    --pattern "*.csv" \
    --output comparison.html
```

**Comparison visualizations:**
- Classification distribution across samples (stacked, grouped and heatmap views)
- Sample diversity overview
- Novelty index distribution per sample
- Novel diversity trend across samples

---

## Common Workflows

### Extract Novel Candidates

After classification, extract reads representing novel diversity:

```bash
metadarkmatter score extract-novel \
    --classifications classifications.csv \
    --output novel_candidates.csv \
    --read-ids novel_reads.txt
```

**Outputs:**
- `novel_candidates.csv`: By default one row per closest reference genome
  (`read_count`, `mean_novelty`, `min_novelty`, `max_novelty`,
  `mean_identity`, `mean_uncertainty`, `dominant_category`); use
  `--no-group` for one row per read with all metrics
- `novel_reads.txt`: Read IDs only (for downstream extraction)

Reads are selected by `taxonomic_call` (Novel Species and/or Novel Genus,
`--category`), or by `--min-novelty` / `--max-uncertainty` when given.

Use read IDs to extract sequences:
```bash
seqtk subseq reads.fastq novel_reads.txt > novel_reads.fastq
```

### Batch Processing

Process multiple samples with shared ANI matrix:

```bash
metadarkmatter score batch \
    --alignment-dir blast_results/ \
    --ani ani_matrix.csv \
    --output-dir classifications/
```

By default processes all files matching `*.blast.tsv.gz` in the directory (change with `--pattern`).

---

## Performance Modes

For large datasets, use performance optimization:

```bash
# Standard mode (suitable for most datasets)
metadarkmatter score classify \
    --alignment huge_sample.blast.tsv.gz \
    --ani ani_matrix.csv \
    --output classifications.csv

# Streaming mode (very large files, bounded memory)
metadarkmatter score classify \
    --alignment huge_sample.blast.tsv.gz \
    --ani ani_matrix.csv \
    --output classifications.csv \
    --streaming
```

Streaming mode requires each read's hits to be contiguous in the file (as
BLAST and MMseqs2 write them) and stops with an error otherwise. It does
not support `--id-mapping`/`--genomes`; `--metadata` is used for the
representative mapping and family inference, but species/genus/family
columns are not added. `--chunk-size` (default 1,000,000 rows) trades
memory for fewer batches.

For details, see [PERFORMANCE.md](PERFORMANCE.md).

---

## Tips and Best Practices

### ANI Matrix Coverage

Validate ANI matrix covers all genomes in BLAST database:

```bash
metadarkmatter ani validate \
    --ani ani_matrix.csv \
    --blast sample.blast.tsv.gz
```

### Species-Level Tracking

Always use the `--metadata` flag for species-level insights (it adds
species/genus/family columns and infers the target family):

```bash
metadarkmatter score classify \
    --alignment sample.blast.tsv.gz \
    --ani ani_matrix.csv \
    --metadata genome_metadata.tsv \
    --output classifications.csv
```

### Protein Mode for Divergent Taxa

For genus-level or higher novelty, use protein mode:

```bash
# Use BLASTX results; placement uncertainty still comes from the ANI matrix
metadarkmatter score classify \
    --alignment sample.blastx.tsv.gz \
    --ani ani_matrix.csv \
    --alignment-mode protein \
    --output classifications.csv
```

Protein mode legacy thresholds (fixed; custom, YAML or adaptive novelty
thresholds are not applied in protein mode):
- Known Species: N < 10%, U < 5%
- Novel Species: 10% <= N < 25%, U < 5%
- Novel Genus: 25% <= N <= 40%, U < 5%

---

## Next Steps

- **Performance optimization**: [PERFORMANCE.md](PERFORMANCE.md)
- **Troubleshooting issues**: [TROUBLESHOOTING.md](TROUBLESHOOTING.md)
- **Complete workflow**: [TUTORIAL_ENVIRONMENTAL_SPECIES.md](TUTORIAL_ENVIRONMENTAL_SPECIES.md)
- **Algorithm details**: [REFERENCE.md](REFERENCE.md)
