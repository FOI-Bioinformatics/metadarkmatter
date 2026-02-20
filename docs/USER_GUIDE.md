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

Metadarkmatter accepts standard BLAST tabular output (`-outfmt 6`) with exactly 12 columns:

```
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
```

**Column requirements:**
- `qseqid`: Read identifier
- `sseqid`: Reference genome identifier (format: `{accession}|{contig_id}`)
- `pident`: Percent identity (0-100)
- `bitscore`: BLAST bitscore (used for top hit selection)

#### Generating BLAST Results

Use `metadarkmatter blast align` which handles:
- FASTQ/FASTA input (automatic conversion)
- Competitive recruitment parameters
- Proper output format

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
       -outfmt 6 \
       -out sample.blast.tsv \
       -max_target_seqs 500
```

#### Compressed Files

BLAST output can be gzip-compressed (`.tsv.gz`). This is recommended for large files:
- 10M alignments: ~500 MB uncompressed → ~50 MB compressed
- Automatic detection and decompression during parsing

### ANI Matrix Format

CSV file with pairwise ANI values between reference genomes:

```csv
genome1,genome2,ani
GCF_000195955.1,GCF_000195955.1,100.00
GCF_000195955.1,GCF_000242755.1,92.15
GCF_000242755.1,GCF_000195955.1,92.15
GCF_000242755.1,GCF_000242755.1,100.00
```

**Requirements:**
- Column names: `genome1`, `genome2`, `ani`
- ANI values: 0-100 (percent)
- Symmetric matrix (optional but recommended)
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
accession	species	genus	family	gtdb_taxonomy
GCF_000195955.1	Francisella_A tularensis	Francisella_A	Francisellaceae	d__Bacteria;p__...
```

**Required columns:**
- `accession`: Genome accession (must match BLAST sseqid prefix)
- `species`: Species name
- `genus`: Genus name
- `family`: Family name
- `gtdb_taxonomy`: Full GTDB taxonomy string

This file enables species-level tracking in classifications and reports.

---

## Output Interpretation

### Classification Results File

Main output contains per-read classifications:

| Column | Type | Description |
|--------|------|-------------|
| read_id | string | Original read identifier from BLAST |
| best_match_genome | string | Genome with highest BLAST bitscore |
| top_hit_identity | float | Percent identity of best BLAST hit (0-100) |
| novelty_index | float | 100 - top_hit_identity; measures divergence |
| placement_uncertainty | float | 100 - ANI(top, secondary); measures ambiguity |
| num_ambiguous_hits | int | Hits within 95% of top bitscore |
| taxonomic_call | string | Classification category |
| is_novel | bool | True if novel species or novel genus |

**Optional columns** (when `--metadata` provided):
- `species`: Species of best match genome
- `genus`: Genus of best match genome

### Taxonomic Classification Categories

Reads are assigned to four categories based on novelty and placement uncertainty:

#### 1. Known Species

**Criteria:**
- Novelty Index < 4%
- Placement Uncertainty < 1.5%

**Interpretation:** Read maps with high identity and low ambiguity to a known reference genome. Likely originates from the matched species or closely related strain.

**Example:**
```csv
read_001,GCF_000195955.1,98.75,1.25,0.18,2,Known Species,false
```

**Biological significance:** Known diversity, not a discovery target.

#### 2. Novel Species

**Criteria:**
- Novelty Index: 4-20%
- Placement Uncertainty < 1.5%

**Interpretation:** Moderate divergence from closest reference but unambiguous placement to a single lineage. Suggests a novel species within a known genus.

**Example:**
```csv
read_002,GCF_000195955.1,91.20,8.80,0.35,5,Novel Species,true
```

**Biological significance:** Putative novel species. Recommend validation through:
- Whole-genome recovery (binning, assembly)
- 16S rRNA analysis
- Targeted cultivation

#### 3. Novel Genus

**Criteria:**
- Novelty Index: 20-25%
- Placement Uncertainty < 1.5%

**Interpretation:** High divergence from all references but moderate-confidence placement. Suggests a novel genus.

**Example:**
```csv
read_003,GCF_000242755.1,82.15,17.85,1.42,3,Novel Genus,true
```

**Biological significance:** Putative novel genus. Priority targets for:
- Single-cell genomics
- Cultivation efforts
- Phylogenetic characterization

#### 4. Conserved Region

**Criteria:**
- Placement Uncertainty >= 5%

**Interpretation:** Read maps with similar scores to multiple divergent genomes. Indicates conserved genomic region (housekeeping genes, rRNA, mobile elements). Taxonomic placement is ambiguous.

**Example:**
```csv
read_004,GCF_000195955.1,95.30,4.70,6.25,15,Conserved Region,false
```

**Biological significance:** Exclude from novel diversity estimates. Cannot be confidently assigned to specific lineage.

### Threshold Interpretation

| Metric | Threshold | Biological Meaning |
|--------|-----------|-------------------|
| Novelty < 4% | Known Species | >96% identity to reference |
| Novelty 4-20% | Novel Species | 80-96% identity (within genus) |
| Novelty 20-25% | Novel Genus | 75-80% identity (genus-level divergence) |
| Uncertainty < 1.5% | Confident | ANI >98.5% between competing genomes |
| Uncertainty >= 1.5% | Ambiguous | ANI <98.5% between competing genomes |

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

**Report sections:**
- **Overview**: Summary statistics and classification distribution
- **ANI Heatmap**: Genome similarity matrix (hierarchically clustered)
- **Species Breakdown**: Composition by species (when metadata provided)
- **Novel Candidates**: Detailed table of novel reads

### Multi-Sample Comparison

Compare multiple samples:

```bash
metadarkmatter report multi \
    --classifications sample1.csv sample2.csv sample3.csv \
    --names "Site A" "Site B" "Site C" \
    --output comparison.html
```

**Comparison visualizations:**
- Classification distribution across samples
- Novelty distribution patterns
- Species composition (with metadata)

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
- `novel_candidates.csv`: Novel reads with full metrics
- `novel_reads.txt`: Read IDs only (for downstream extraction)

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
    --metadata genome_metadata.tsv \
    --output-dir classifications/
```

Automatically processes all `.blast.tsv` or `.blast.tsv.gz` files in the directory.

---

## Performance Modes

For large datasets, use performance optimization:

```bash
# Standard mode (suitable for most datasets)
metadarkmatter score classify \
    --alignment huge_sample.blast.tsv.gz \
    --ani ani_matrix.csv \
    --output classifications.csv

# Streaming mode (100M+ alignments, memory-efficient)
metadarkmatter score classify \
    --alignment huge_sample.blast.tsv.gz \
    --ani ani_matrix.csv \
    --output classifications.csv \
    --streaming
```

For details, see [PERFORMANCE.md](PERFORMANCE.md).

---

## Tips and Best Practices

### ANI Matrix Coverage

Validate ANI matrix covers all genomes in BLAST database:

```bash
metadarkmatter ani validate \
    --ani ani_matrix.csv \
    --genomes genomes/
```

### Species-Level Tracking

Always use `--metadata` flag for species-level insights:

```bash
metadarkmatter score classify \
    --alignment sample.blast.tsv.gz \
    --ani ani_matrix.csv \
    --metadata genome_metadata.tsv \  # Enables species tracking
    --output classifications.csv
```

### Protein Mode for Divergent Taxa

For genus-level or higher novelty, use protein mode:

```bash
# Use BLASTX results and AAI matrix
metadarkmatter score classify \
    --alignment sample.blastx.tsv.gz \
    --ani aai_matrix.csv \
    --alignment-mode protein \  # Different thresholds
    --output classifications.csv
```

Protein mode thresholds:
- Known Species: N < 10%, U < 5%
- Novel Species: 10% <= N < 25%, U < 5%
- Novel Genus: 25% <= N <= 40%, U < 5%

---

## Next Steps

- **Performance optimization**: [PERFORMANCE.md](PERFORMANCE.md)
- **Troubleshooting issues**: [TROUBLESHOOTING.md](TROUBLESHOOTING.md)
- **Complete workflow**: [TUTORIAL_ENVIRONMENTAL_SPECIES.md](TUTORIAL_ENVIRONMENTAL_SPECIES.md)
- **Algorithm details**: [REFERENCE.md](REFERENCE.md)
