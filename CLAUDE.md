# CLAUDE.md

Quick reference for Claude Code when working with this repository.

## Project Overview

**metadarkmatter** is a microbial ecology tool for detecting novel bacterial taxa ("microbial dark matter") from metagenomic sequencing data.

**Core Innovation:** ANI-weighted placement uncertainty - combines BLAST alignment with ANI/AAI matrices to distinguish ambiguous placement within known taxa from confident identification of novel species or genera.

**Supported Workflows:**
- Nucleotide-level classification (BLASTN + ANI)
- Protein-level classification (BLASTX + AAI) for highly divergent taxa

## Quick Start

```bash
# 1. Download reference genomes
metadarkmatter download genomes list "f__Francisellaceae" --output genomes.tsv
metadarkmatter download genomes fetch --accessions genomes.tsv --output-dir genomes/

# 2. Extract family reads with Kraken2
metadarkmatter kraken2 extract --kraken-output sample.kraken --reads-1 sample_R1.fastq.gz \
  --reads-2 sample_R2.fastq.gz --taxid 119060 --output extraction/

# 3. Build BLAST database and align (accepts FASTQ directly)
metadarkmatter blast makedb --genomes genomes/ --output blastdb/pangenome
metadarkmatter blast align --query extraction/reads_R1.fastq.gz --database blastdb/pangenome \
  --output sample.blast.tsv.gz --threads 16

# 4. Compute ANI matrix
metadarkmatter ani compute --genomes genomes/ --output ani_matrix.csv --threads 16

# 5. Classify reads (core algorithm)
metadarkmatter score classify --blast sample.blast.tsv.gz --ani ani_matrix.csv \
  --metadata genome_metadata.tsv --output classifications.csv --parallel

# 6. Generate HTML report
metadarkmatter report generate --classifications classifications.csv \
  --metadata genome_metadata.tsv --output report.html
```

## Workflow Notes

**BLAST Input Formats:**
- BLAST accepts FASTQ, FASTA, or gzipped versions (.gz)
- FASTQ files are automatically converted to FASTA (transparent to user)
- No manual conversion with seqtk required

**MMseqs2 for Large Datasets:**
- Available for datasets with >100,000 reads
- For <10K reads: MMseqs2 is slower than BLAST (use BLAST)
- For 100K-1M+ reads: MMseqs2 provides 5-15x speedup
- See [Tutorial](docs/TUTORIAL_ENVIRONMENTAL_SPECIES.md) Step 7 for decision guide

## Documentation

- **[docs/TUTORIAL_ENVIRONMENTAL_SPECIES.md](docs/TUTORIAL_ENVIRONMENTAL_SPECIES.md)** - Complete tutorial for finding novel diversity
- **[docs/REFERENCE.md](docs/REFERENCE.md)** - CLI reference, algorithm details, advanced features
- **[docs/WORKFLOW.md](docs/WORKFLOW.md)** - Detailed workflow guidance and database strategies
- **[docs/CLASSIFICATION_STATISTICS.md](docs/CLASSIFICATION_STATISTICS.md)** - Statistical framework and literature references

## Important Conventions

- Use modest scientific language in documentation and code comments
- Avoid Unicode in any Nextflow workflow files (if workflows are added)
- FASTA headers are standardized to `{accession}|{contig_id}` format in BLAST databases
- Classification thresholds differ between nucleotide (default) and protein mode (`--alignment-mode protein`)

## Critical Known Issues

### ANI/AAI Heatmap Clustering (FIXED)

**Status:** ✅ FIXED - scipy added to dependencies, clustering now works reliably.

**What Was Fixed:**
1. **Added scipy>=1.11.0 to dependencies** (`pyproject.toml:33`)
2. **Extracted duplicate clustering logic** to shared `_perform_hierarchical_clustering()` method
3. **Made titles conditional** - only shows "(Hierarchically Clustered)" when clustering succeeds
4. **Added logging** - warns users if scipy unavailable (graceful degradation)

**Location:** `src/metadarkmatter/visualization/report/generator.py`
- Shared clustering method: `_perform_hierarchical_clustering()` (lines 1039-1103)
- ANI heatmap: `_build_ani_section()` uses shared method
- AAI heatmap: `_build_aai_section()` uses shared method

**What This Means:**
- Clustering now works by default for all users
- Heatmaps show genomes grouped by similarity (species/genus blocks visible)
- Titles accurately reflect whether clustering was applied
- Logs warning if scipy somehow unavailable (rare edge case)

**Impact:**
- Better scientific interpretation of ANI/AAI matrices
- Visual clustering makes taxonomic structure immediately apparent
- Same-species genomes now grouped together in heatmap

## Key Architecture

```
src/metadarkmatter/
├── cli/              # Typer CLI commands (entry points)
├── core/             # Core algorithms (ANI classification, parsers, metadata)
│   ├── constants.py          # Nucleotide thresholds (default mode)
│   └── protein_constants.py  # Protein thresholds (--alignment-mode protein)
├── external/         # External tool wrappers (BLAST, MMseqs2, Kraken2, etc.)
├── models/           # Pydantic data models
├── clients/          # API clients (GTDB with retry logic)
└── visualization/    # Plotly charts and HTML report generation
    └── report/
        └── generator.py  # Single-sample reports with heatmap clustering
```

## Classification Thresholds

### Nucleotide Mode (BLASTN, default)
- Known Species: N < 5%, U < 2%
- Novel Species: 5% <= N < 20%, U < 2%
- Novel Genus: 20% <= N <= 25%, U < 2%

### Protein Mode (BLASTX, `--alignment-mode protein`)
- Known Species: N < 10%, U < 5%
- Novel Species: 10% <= N < 25%, U < 5%
- Novel Genus: 25% <= N <= 40%, U < 5%

Where:
- **N** = Novelty Index = `100 - TopHitIdentity`
- **U** = Placement Uncertainty = `100 - max(ANI between competing genomes)`

See [docs/REFERENCE.md](docs/REFERENCE.md) for complete threshold tables and interpretation.
