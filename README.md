# Metadarkmatter

ANI-weighted placement uncertainty for detecting novel microbial diversity in environmental eDNA samples.

## Overview

Metadarkmatter analyzes environmental DNA (eDNA) metagenomic data to detect novel bacterial taxa ("microbial dark matter") using competitive read recruitment. The tool calculates:

- **Novelty Index**: Sequence divergence from reference genomes
- **Placement Uncertainty**: Ambiguity when reads match multiple genomes
- **Species-Level Tracking**: Automatic metadata propagation from GTDB taxonomy

Reads are classified as: Known Species, Novel Species, Novel Genus, or Conserved Region.

## Installation

### Python Package

```bash
pip install -e .
```

**Requirements:** Python >= 3.11

### External Software Dependencies

Metadarkmatter requires several external bioinformatics tools (non-Python software):

#### Required Tools

Install via conda/mamba:

```bash
conda create -n metadarkmatter -c conda-forge  -c bioconda kraken2 krakentools blast skani seqtk mmseqs2 diamond pyani
```

| Tool | Purpose | When Needed |
|------|---------|--------------|
| **Kraken2** | Taxonomic classification of reads | Classifying reads to target bacterial family |
| **KrakenTools** | Extract reads by taxid | Filtering reads to target bacterial family |
| **BLAST+** | Nucleotide alignment (blastn) |  |
| **skani** | Faster ANI computation | |
| **pyani** | Fast ANI computation | |
| **MMseqs2** | Fast sequence search (BLAST alternative) | Datasets with >100,000 reads |
| **seqtk** | Sequence manipulation | Assembly workflows, read extraction |
| **diamond** | Fast AAI computation  |  |

**Note:** BLAST accepts FASTQ files directly (automatic conversion), so seqtk is only needed for specialized workflows like assembly.

## Quick Start

```bash
# 1. Download reference genomes (auto-creates genome_metadata.tsv)
metadarkmatter download genomes list "f__Francisellaceae" --output genomes.tsv
metadarkmatter download genomes fetch --accessions genomes.tsv --output-dir genomes/

# 2. Classify reads with Kraken2
metadarkmatter kraken2 classify \
  --reads-1 sample_R1.fastq.gz \
  --kraken-db /path/to/kraken_db \
  --output kraken_output/

# 3. Extract target family reads
metadarkmatter kraken2 extract \
  --kraken-output kraken_output/sample.kraken \
  --reads-1 sample_R1.fastq.gz \
  --taxid FAMILY_TAXID \
  --output extraction/

# 4. Build BLAST database (standardizes contig headers)
metadarkmatter blast makedb --genomes genomes/ --output blastdb/pangenome

# 5. Run competitive BLAST alignment (accepts FASTQ/FASTA, gzipped or plain)
metadarkmatter blast align --query extraction/sample_taxid262_R1.fastq.gz \
  --database blastdb/pangenome --output sample.blast.tsv.gz

# 6. Compute ANI matrix
metadarkmatter ani compute --genomes genomes/ --output ani_matrix.csv

# 7. Classify reads (with species-level tracking)
metadarkmatter score classify \
  --blast sample.blast.tsv.gz \
  --ani ani_matrix.csv \
  --metadata genome_metadata.tsv \
  --output classifications.csv \
  --parallel

# 8. Generate report (with species breakdown)
metadarkmatter report generate \
  --classifications classifications.csv \
  --metadata genome_metadata.tsv \
  --output report.html
```

## Commands

| Command | Description |
|---------|-------------|
| `download genomes list` | Query GTDB for family genomes (creates metadata) |
| `download genomes fetch` | Download genomes from NCBI |
| `kraken2 classify` | Run Kraken2 classification |
| `kraken2 extract` | Extract reads for target taxid |
| `blast makedb` | Build BLAST database (standardizes headers) |
| `blast align` | Run competitive BLAST alignment |
| `ani compute` | Compute ANI matrix (skani/fastANI) |
| `ani validate` | Validate ANI matrix coverage |
| `score classify` | ANI-weighted classification (supports `--metadata`) |
| `score batch` | Batch process multiple samples |
| `score extract-novel` | Extract candidate novel species/genera |
| `report generate` | Create HTML report (supports `--metadata` for species tab) |
| `report multi` | Multi-sample comparison |

## Key Features

### Species-Level Tracking

Metadarkmatter automatically tracks species metadata throughout the pipeline:

- **Automatic metadata**: `genome_metadata.tsv` created during genome download with species, genus, family, and full GTDB taxonomy
- **Standardized headers**: FASTA headers rewritten to `{accession}|{contig_id}` format for reliable genome identification from multi-contig draft genomes
- **Species aggregation**: `--metadata` option in `score classify` adds species/genus columns to output
- **Species breakdown**: HTML reports include a dedicated "Species Breakdown" tab with composition charts

### Novel Diversity Detection

- Identify **Novel Species** candidates (5-15% novelty, <2% uncertainty)
- Identify **Novel Genus** candidates (15-25% novelty, <2% uncertainty)
- Extract candidate reads with `score extract-novel` for targeted assembly
- Literature-backed thresholds based on 95-96% ANI species boundary

### Performance

**Classification Performance** (BLAST results → classifications):

| Dataset | Mode | Runtime | RAM |
|---------|------|---------|-----|
| < 1M alignments | default | 2-5 min | 4 GB |
| 1-10M | `--fast` | 5-10 min | 8 GB |
| 10-100M | `--parallel` | 15-45 min | 16 GB |
| 100M+ | `--streaming` | 1-2 hr | 16 GB |

**Alignment Performance** (reads → BLAST results):

BLAST is efficient for typical workflows. For very large datasets (>100K reads), consider MMseqs2:

| Reads | BLAST Time | MMseqs2 Time | When to Use |
|-------|-----------|--------------|-------------|
| <10K | 10s-10min | Slower | Use BLAST |
| 100K | 30-60 min | 5-10 min | Either tool works |
| 1M+ | 3-6 hours | 15-30 min | MMseqs2 recommended |

**Note**: BLAST accepts FASTQ directly (automatic conversion to FASTA). See [Tutorial](docs/TUTORIAL_ENVIRONMENTAL_SPECIES.md) for workflow details.

## Documentation

- **[Workflow Guide](docs/WORKFLOW.md)** - Step-by-step analysis tutorial
- **[CLI Reference](docs/CLI_REFERENCE.md)** - Complete command documentation
- **[User Guide](docs/USER_GUIDE.md)** - Detailed usage examples
- **[Performance Guide](docs/PERFORMANCE.md)** - Optimization tips
- **[Troubleshooting](docs/TROUBLESHOOTING.md)** - Common issues
- **[API Reference](docs/API_REFERENCE.md)** - Python API documentation
- **[Tutorial](docs/TUTORIAL_ENVIRONMENTAL_SPECIES.md)** - Finding novel diversity in environmental samples

## License

MIT License
