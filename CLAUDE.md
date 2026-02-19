# CLAUDE.md

Quick reference for Claude Code when working with this repository.

## Project Overview

**metadarkmatter** is a microbial ecology tool for detecting novel bacterial taxa ("microbial dark matter") from metagenomic sequencing data.

**Core Innovation:** ANI-weighted placement uncertainty - combines BLAST alignment with ANI/AAI matrices to distinguish ambiguous placement within known taxa from confident identification of novel species or genera.

**Supported Workflows:**
- Nucleotide-level classification (BLASTN + ANI)
- Protein-level classification (BLASTX + AAI) for highly divergent taxa
- External alignment import (BLAST/MMseqs2 results run outside metadarkmatter)

## Quick Start

```bash
# 1. Download reference genomes
metadarkmatter download genomes list "f__Francisellaceae" --output genomes.tsv
metadarkmatter download genomes fetch --accessions genomes.tsv --output-dir genomes/

# 2. Extract family reads with Kraken2
metadarkmatter kraken2 extract --kraken-output sample.kraken --reads-1 sample_R1.fastq.gz \
  --reads-2 sample_R2.fastq.gz --taxid 119060 --output extraction/

# 3. Compute ANI matrix
metadarkmatter ani compute --genomes genomes/ --output ani_matrix.csv --threads 16

# 4. Sequence alignment (choose ONE option)

# Option A: BLAST (single-end only, best for <100K reads)
metadarkmatter blast makedb --genomes genomes/ --output blastdb/pangenome
metadarkmatter blast align --query extraction/reads_R1.fastq.gz --database blastdb/pangenome \
  --output sample.blast.tsv.gz --threads 16

# Option B: MMseqs2 (supports paired-end, best for >100K reads, 5-100x faster)
metadarkmatter mmseqs2 makedb --genomes genomes/ --output mmseqs_db/pangenome
metadarkmatter mmseqs2 search \
  --query-1 extraction/reads_R1.fastq.gz --query-2 extraction/reads_R2.fastq.gz \
  --database mmseqs_db/pangenome --output sample.mmseqs2.tsv.gz --threads 16

# 5. Classify reads (core algorithm - works with either BLAST or MMseqs2 output)
metadarkmatter score classify --alignment sample.blast.tsv.gz --ani ani_matrix.csv \
  --metadata genome_metadata.tsv --output classifications.csv --parallel

# 6. Generate HTML report
metadarkmatter report generate --classifications classifications.csv \
  --metadata genome_metadata.tsv --output report.html
```

## Importing External Alignment Results

External BLAST/MMseqs2 results can be classified using ID mapping to translate contig IDs to genome accessions:

```bash
# Generate mapping from genome directory
metadarkmatter util generate-mapping --genomes genomes/ --output id_mapping.tsv

# Validate mapping against BLAST file
metadarkmatter util validate-mapping id_mapping.tsv --blast external_results.tsv.gz

# Classify with mapping (requires --parallel)
metadarkmatter score classify --alignment external_results.tsv.gz --ani ani_matrix.csv \
  --id-mapping id_mapping.tsv --output classifications.csv --parallel
```

**Key Files:**
- `cli/util.py` - `generate-mapping`, `validate-mapping` commands
- `core/id_mapping.py` - `ContigIdMapping` class (from_genome_dir, from_tsv, to_tsv, transform_column)
- `core/genome_utils.py` - `extract_accession_from_filename()`

**Requirements:**
- Alignment must be BLAST tabular format (outfmt 6), 12 or 13 columns, tab-separated
- ID mapping requires `--parallel` mode (not supported with `--fast` or `--streaming`)
- Genome accessions in mapping must match ANI matrix row/column labels

## Workflow Notes

**Alignment Tool Selection:**

| Tool | Input Support | Best For | Speed |
|------|---------------|----------|-------|
| BLAST | Single-end only (`--query`) | <100K reads | Baseline |
| MMseqs2 | Single-end or paired-end (`--query-1`, `--query-2`) | >100K reads | 5-100x faster |

**BLAST:**
- Accepts FASTQ, FASTA, or gzipped versions (.gz)
- FASTQ files are automatically converted to FASTA (transparent to user)
- Single-end input only (`--query`)
- Recommended for smaller datasets where setup time matters

**MMseqs2:**
- Accepts FASTQ, FASTA, or gzipped versions (.gz)
- Supports both single-end (`--query`) and paired-end (`--query-1`, `--query-2`)
- For paired-end: reads are concatenated internally (standard MMseqs2 approach)
- 5-100x faster than BLAST for large datasets
- Recommended for >100K reads

See [Tutorial](docs/TUTORIAL_ENVIRONMENTAL_SPECIES.md) Step 7 for detailed decision guide

## Documentation

- **[docs/METHODS.md](docs/METHODS.md)** - Comprehensive methods documentation with all calculations (scientific standard)
- **[docs/TUTORIAL_ENVIRONMENTAL_SPECIES.md](docs/TUTORIAL_ENVIRONMENTAL_SPECIES.md)** - Complete tutorial for finding novel diversity
- **[docs/REFERENCE.md](docs/REFERENCE.md)** - CLI reference, algorithm details, advanced features
- **[docs/WORKFLOW.md](docs/WORKFLOW.md)** - Detailed workflow guidance and database strategies
- **[docs/CLASSIFICATION_STATISTICS.md](docs/CLASSIFICATION_STATISTICS.md)** - Statistical framework and literature references

## Important Conventions

- Use modest scientific language in documentation and code comments
- Avoid Unicode in any Nextflow workflow files (if workflows are added)
- FASTA headers are standardized to `{accession}|{contig_id}` format in BLAST databases
- Classification thresholds differ between nucleotide (default) and protein mode (`--alignment-mode protein`)

## Phylogeny Tab

The HTML report includes an interactive phylogenetic tree tab showing novel clusters in evolutionary context.

**CLI Options:**
```bash
# Default: auto-generate NJ tree from ANI matrix
metadarkmatter report generate --classifications results.csv --ani ani.csv --output report.html

# Use custom tree
metadarkmatter report generate --classifications results.csv --ani ani.csv --tree custom.nwk --output report.html

# Skip phylogeny tab (faster for large datasets)
metadarkmatter report generate --classifications results.csv --ani ani.csv --no-phylogeny --output report.html
```

**Key Files:**
- `core/phylogeny/tree_builder.py` - `ani_to_newick()`, `load_user_tree()`
- `core/phylogeny/placement.py` - `NovelCluster`, `extract_novel_clusters()`, `place_novel_clusters()`
- `visualization/report/templates.py` - `PHYLOGENY_SECTION_TEMPLATE`, `PHYLOTREE_JS_TEMPLATE`
- `visualization/report/generator.py` - `_build_phylogeny_section()`

**Requirements:**
- ANI matrix with >= 3 genomes
- BioPython (dependency in pyproject.toml)
- D3.js loaded from CDN in report

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
│   ├── protein_constants.py  # Protein thresholds (--alignment-mode protein)
│   ├── classification/       # Classification pipeline
│   │   ├── thresholds.py     # Extracted threshold application logic
│   │   ├── qc.py             # Pre/post-classification QC metrics
│   │   ├── sensitivity.py    # Threshold sensitivity analysis
│   │   ├── adaptive.py       # GMM-based adaptive threshold detection
│   │   ├── bayesian.py       # Bayesian posterior probabilities
│   │   ├── ani_matrix.py     # ANIMatrix class
│   │   └── classifiers/      # Classifier implementations
│   │       ├── base.py       # Base classifier with coverage weighting
│   │       ├── vectorized.py # Polars-based vectorized classifier
│   │       └── parallel.py   # Multiprocessing classifier
│   └── phylogeny/            # Phylogenetic tree building and novel cluster placement
│       ├── tree_builder.py   # ANI-to-Newick conversion, user tree loading
│       └── placement.py      # Novel cluster extraction and tree placement
├── external/         # External tool wrappers (BLAST, MMseqs2, Kraken2, etc.)
├── models/           # Pydantic data models
├── clients/          # API clients (GTDB with retry logic)
└── visualization/    # Plotly charts and HTML report generation
    └── report/
        ├── generator.py      # Single-sample reports (with Bayesian tab)
        └── components/       # Report building blocks
            ├── clustering.py # Hierarchical clustering
            └── heatmap_builder.py # ANI/AAI heatmaps
```

## Classification Thresholds

### Nucleotide Mode (BLASTN, default)
- Known Species: N < 4%, U < 1.5%
- Novel Species: 4% <= N < 20%, U < 1.5%
- Novel Genus: 20% <= N <= 25%, U < 1.5%

### Protein Mode (BLASTX, `--alignment-mode protein`)
- Known Species: N < 10%, U < 5%
- Novel Species: 10% <= N < 25%, U < 5%
- Novel Genus: 25% <= N <= 40%, U < 5%

Where:
- **N** = Novelty Index = `100 - TopHitIdentity`
- **U** = Placement Uncertainty = `100 - max(ANI between competing genomes)`

See [docs/REFERENCE.md](docs/REFERENCE.md) for complete threshold tables and interpretation.

## Coverage-Weighted Hit Selection

Linear coverage weighting is enabled by default to reduce bias from short conserved domains.

**Weighting Modes:**
- `none`: Disables coverage weighting (raw bitscore)
- `linear` (default): `weight = min + (max - min) * coverage`
- `log`: `weight = min + (max - min) * log(1 + 9*coverage) / log(10)`
- `sigmoid`: `weight = min + (max - min) / (1 + exp(-10*(coverage - 0.6)))`

Where `min = 1 - strength`, `max = 1 + strength`, `strength` defaults to 0.5.

**Key Files:**
- `models/blast.py` - `BlastHit.calculate_coverage()`, `calculate_weighted_score()`
- `models/config.py` - `ScoringConfig.coverage_weight_mode`, `coverage_weight_strength`
- `core/classification/classifiers/vectorized.py` - Polars-based vectorized classifier

## Advanced Classification Features

### QC Metrics (`--qc-output`)
- Pre-classification: filter rate, genome coverage, single-hit fraction, mean identity
- Post-classification: ambiguous fraction, low-confidence fraction
- Automated warnings for problematic inputs
- Key file: `core/classification/qc.py`

### Adaptive Thresholds (`--adaptive-thresholds`)
- Detects species boundary from ANI matrix via 2-component GMM
- Falls back to default 96% ANI if GMM does not converge
- Requires: `pip install metadarkmatter[adaptive]` (scikit-learn)
- Key file: `core/classification/adaptive.py`

### Bayesian Confidence (`--bayesian`)
- Posterior probabilities P(category | novelty, uncertainty) via 2D Gaussian likelihoods
- Shannon entropy as confidence metric (0 = confident, 2.0 = uniform)
- Adds 6 columns: p_known_species, p_novel_species, p_novel_genus, p_ambiguous, bayesian_category, posterior_entropy
- Vectorized numpy implementation (no Python loops)
- Key file: `core/classification/bayesian.py`

### Family Validation (`--target-family`)
- Detects off-target reads from broad-database alignments
- Partitions hits by ANI matrix membership (in-family vs external)
- Reads with best_in_family / best_overall bitscore < 0.8 are Off-target
- New output columns: family_bitscore_ratio, family_identity_gap, in_family_hit_fraction
- Key file: `core/classification/classifiers/vectorized.py`

### Sensitivity Analysis (`score sensitivity` subcommand)
- Sweeps novelty/uncertainty thresholds across configurable range
- Re-classifies reads at each threshold point
- Key file: `core/classification/sensitivity.py`

### Enhanced Scoring (always on)
- Inferred uncertainty for single-hit reads
- Alignment quality, identity confidence, placement confidence
- Discovery score for novel taxa prioritization
- All computed automatically (no flags needed)
