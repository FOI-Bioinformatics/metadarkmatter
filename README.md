# Metadarkmatter

[![Python](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)

ANI-weighted placement uncertainty for detecting novel microbial diversity in environmental eDNA samples.

## Overview

Metadarkmatter identifies novel bacterial taxa ("microbial dark matter") in metagenomic data using competitive read recruitment with ANI-weighted placement uncertainty. The tool distinguishes between true novelty and taxonomic ambiguity by analyzing both sequence divergence and phylogenetic placement confidence.

**Core Innovation:** Combines BLAST alignment with Average Nucleotide Identity (ANI) matrices to separate reads from genuinely novel taxa from reads that simply match multiple reference genomes with similar identity.

### Key Metrics

- **Novelty Index (N)**: Sequence divergence from closest reference genome
- **Placement Uncertainty (U)**: Phylogenetic ambiguity based on ANI between competing genomes
- **Classification**: Known Species, Novel Species, Novel Genus, or Conserved Region

### How It Works

1. **Competitive alignment**: Reads aligned against all reference genomes simultaneously
2. **ANI-weighted scoring**: Top hits evaluated using genome-genome ANI matrix
3. **Uncertainty quantification**: Placement confidence computed from ANI between competing matches
4. **Threshold-based classification**: Literature-backed boundaries for species (95% ANI) and genus (75% ANI)

## Installation

### 1. Install Python Package

```bash
# Clone repository
git clone https://github.com/FOI-Bioinformatics/metadarkmatter.git
cd metadarkmatter

# Install with pip
pip install -e .
```

**Requirements:** Python >= 3.11

### 2. Install External Tools

Metadarkmatter requires external bioinformatics tools. Install via conda/mamba:

```bash
conda create -n metadarkmatter -c conda-forge -c bioconda \
  kraken2 krakentools blast skani mmseqs2 diamond pyani
conda activate metadarkmatter
```

**Core Tools:**
- **Kraken2 + KrakenTools**: Taxonomic read classification and extraction
- **BLAST+**: Nucleotide alignment (accepts FASTQ directly)
- **skani**: Fast ANI computation
- **MMseqs2**: Fast sequence search for large datasets (>100K reads)
- **Diamond**: Protein alignment and AAI computation (for protein mode)

**Optional:** `seqtk` (assembly workflows), `pyani` (alternative ANI tool)

## Quick Start

Complete workflow for detecting novel diversity in Francisellaceae:

```bash
# 1. Download reference genomes from GTDB
mdm download genomes list "f__Francisellaceae" --output genomes.tsv
mdm download genomes fetch --accessions genomes.tsv --output-dir genomes/

# 2. Extract target family reads (requires Kraken2 database)
mdm kraken2 classify --reads-1 sample_R1.fastq.gz --kraken-db db/ --output kraken_out/
mdm kraken2 extract --kraken-output kraken_out/sample.kraken \
  --reads-1 sample_R1.fastq.gz --taxid 119060 --output extraction/

# 3. Build BLAST database and align
mdm blast makedb --genomes genomes/ --output blastdb/pangenome
mdm blast align --query extraction/reads_R1.fastq.gz \
  --database blastdb/pangenome --output sample.blast.tsv.gz --threads 16

# 4. Compute ANI matrix
mdm ani compute --genomes genomes/ --output ani_matrix.csv --threads 16

# 5. Classify reads and generate report
mdm score classify --alignment sample.blast.tsv.gz --ani ani_matrix.csv \
  --metadata genome_metadata.tsv --output classifications.csv
mdm report generate --classifications classifications.csv \
  --metadata genome_metadata.tsv --output report.html
```

**Note:** Both `metadarkmatter` and `mdm` commands are available. See [Tutorial](docs/TUTORIAL_ENVIRONMENTAL_SPECIES.md) for detailed walkthrough.

## Commands

| Command | Description |
|---------|-------------|
| `download genomes list` | Query GTDB for family genomes (creates metadata) |
| `download genomes fetch` | Download genomes from NCBI |
| `kraken2 classify` | Run Kraken2 classification |
| `kraken2 extract` | Extract reads for target taxid |
| `blast makedb` | Build BLAST database (standardizes headers) |
| `blast align` | Run competitive BLAST alignment |
| `mmseqs2 makedb` | Build MMseqs2 database (for large datasets) |
| `mmseqs2 search` | MMseqs2 search (supports paired-end) |
| `blastx align` | Protein-level alignment with Diamond |
| `ani compute` | Compute ANI matrix (skani/fastANI) |
| `aai compute` | Compute AAI matrix (Diamond) |
| `score classify` | ANI-weighted classification (core algorithm) |
| `score batch` | Batch process multiple samples |
| `score extract-novel` | Extract candidate novel species/genera |
| `report generate` | Create HTML report with interactive tabs |
| `report multi` | Multi-sample comparison |
| `util generate-mapping` | Generate contig-to-genome ID mapping |

## Key Features

### Species-Level Tracking

Metadarkmatter automatically tracks species metadata throughout the pipeline:

- **Automatic metadata**: `genome_metadata.tsv` created during genome download with species, genus, family, and full GTDB taxonomy
- **Standardized headers**: FASTA headers rewritten to `{accession}|{contig_id}` format for reliable genome identification from multi-contig draft genomes
- **Species aggregation**: `--metadata` option in `score classify` adds species/genus columns to output
- **Species breakdown**: HTML reports include a dedicated "Species Breakdown" tab with composition charts

### Novel Diversity Detection

- Identify **Novel Species** candidates (4-20% novelty, <1.5% uncertainty)
- Identify **Novel Genus** candidates (20-25% novelty, <1.5% uncertainty)
- Extract candidate reads with `score extract-novel` for targeted assembly
- Literature-backed thresholds based on 95-96% ANI species boundary

### Coverage-Weighted Hit Selection

Optional feature to prioritize alignments spanning larger portions of reads:

- **Problem solved**: Short conserved domains (e.g., 16S rRNA) with high identity can dominate over longer, more informative alignments
- **Solution**: Weight bitscore by alignment coverage to favor hits that explain more of the read
- **Modes**: Linear, logarithmic, or sigmoid weighting functions
- **Usage**: `--coverage-weight-mode linear` in `score classify`

### Performance

**Classification Performance** (BLAST results → classifications):

All classification uses the Polars-based vectorized engine with automatic parallelization.

| Dataset | Mode | Runtime | RAM |
|---------|------|---------|-----|
| < 10M alignments | default | 2-10 min | 4-8 GB |
| 10-100M | default | 15-45 min | 16 GB |
| 100M+ | `--streaming` | 1-2 hr | 16 GB (constant) |

**Alignment Performance** (reads → BLAST results):

BLAST is efficient for typical workflows. For very large datasets (>100K reads), consider MMseqs2:

| Reads | BLAST Time | MMseqs2 Time | When to Use |
|-------|-----------|--------------|-------------|
| <10K | 10s-10min | Slower | Use BLAST |
| 100K | 30-60 min | 5-10 min | Either tool works |
| 1M+ | 3-6 hours | 15-30 min | MMseqs2 recommended |

**Note**: BLAST accepts FASTQ directly (automatic conversion to FASTA). See [Tutorial](docs/TUTORIAL_ENVIRONMENTAL_SPECIES.md) for workflow details.

## Use Cases

### Environmental Monitoring
- Track novel pathogen emergence in environmental samples
- Monitor microbial diversity shifts in ecosystems
- Identify candidate novel species for targeted isolation

### Biosurveillance
- Detect divergent bacterial lineages in clinical or environmental eDNA
- Quantify taxonomic coverage gaps relative to reference databases
- Prioritize samples for further characterization

### Research Applications
- Characterize uncultured microbial diversity
- Validate reference genome coverage for target taxa
- Generate candidate lists for genome-resolved metagenomics

## Documentation

### Getting Started
- **[Tutorial](docs/TUTORIAL_ENVIRONMENTAL_SPECIES.md)** - Complete walkthrough for environmental samples
- **[Workflow Guide](docs/WORKFLOW.md)** - Step-by-step analysis patterns
- **[User Guide](docs/USER_GUIDE.md)** - Detailed usage examples

### Reference
- **[CLI Reference](docs/CLI_REFERENCE.md)** - Complete command documentation
- **[API Reference](docs/API_REFERENCE.md)** - Python API documentation
- **[Algorithm Details](docs/CLASSIFICATION_STATISTICS.md)** - Statistical framework and literature references

### Advanced
- **[Performance Guide](docs/PERFORMANCE.md)** - Optimization strategies
- **[Troubleshooting](docs/TROUBLESHOOTING.md)** - Common issues and solutions
- **[Architecture](docs/ARCHITECTURE.md)** - System design and internals

## Citation

If you use metadarkmatter in your research, please cite:

```bibtex
@software{metadarkmatter2026,
  author = {Metadarkmatter Team},
  title = {Metadarkmatter: ANI-weighted placement uncertainty for detecting novel microbial diversity},
  year = {2026},
  url = {https://github.com/FOI-Bioinformatics/metadarkmatter}
}
```

### Related Methods

This tool implements concepts from:
- **ANI species boundary**: Jain et al. (2018). "High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries." *Nature Communications* 9:5114.
- **Competitive recruitment**: Rodriguez-R et al. (2018). "The Microbial Genomes Atlas (MiGA) webserver." *Nucleic Acids Research*.

## Contributing

Contributions are welcome! Please:
- Open an issue for bug reports or feature requests
- Submit pull requests for code contributions
- Follow existing code style (ruff formatting, type hints)
- Add tests for new features

## Support

- **Issues**: [GitHub Issues](https://github.com/FOI-Bioinformatics/metadarkmatter/issues)
- **Documentation**: [docs/](docs/)
- **Tutorial**: [docs/TUTORIAL_ENVIRONMENTAL_SPECIES.md](docs/TUTORIAL_ENVIRONMENTAL_SPECIES.md)

## License

MIT License - see [LICENSE](LICENSE) file for details.
