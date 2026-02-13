# metadarkmatter Reference Guide

Complete technical reference for the metadarkmatter package.

## Table of Contents

- [CLI Workflows](#cli-workflows)
  - [Standard Nucleotide Workflow](#standard-nucleotide-workflow)
  - [Protein-Level Workflow](#protein-level-workflow)
- [Command Reference](#command-reference)
- [Core Algorithm](#core-algorithm)
  - [Nucleotide Mode Thresholds](#nucleotide-mode-thresholds)
  - [Protein Mode Thresholds](#protein-mode-thresholds)
- [GTDB Compatibility](#gtdb-compatibility)
- [Database Strategy](#database-strategy)
- [Project Structure](#project-structure)
- [Advanced Features](#advanced-features)
  - [Performance Modes](#performance-modes)
  - [ANI Matrix Computation](#ani-matrix-computation)
  - [Novel Species Extraction](#novel-species-extraction)
  - [Species-Level Tracking](#species-level-tracking)
- [Development Guide](#development-guide)

---

## CLI Workflows

### Standard Nucleotide Workflow

Complete pipeline for nucleotide-level (BLASTN) read classification:

```bash
# Step 1: Download reference genomes for your bacterial family
# Note: This automatically creates genome_metadata.tsv with species/genus info
metadarkmatter download genomes list "f__YourFamily" --output genomes.tsv
metadarkmatter download genomes fetch --accessions genomes.tsv --output-dir genomes/ --include-protein

# Step 1b: Generate missing protein files with Prodigal (if any)
metadarkmatter proteins predict --genomes genomes/ --missing-only --threads 16

# Step 2a: Run Kraken2 taxonomic classification
metadarkmatter kraken2 classify \
  --reads-1 sample_R1.fastq.gz \
  --reads-2 sample_R2.fastq.gz \
  --kraken-db /path/to/kraken2_db \
  --output kraken_output/

# Step 2b: Extract reads for target family (outputs compressed by default)
metadarkmatter kraken2 extract \
  --kraken-output kraken_output/sample.kraken \
  --reads-1 sample_R1.fastq.gz \
  --reads-2 sample_R2.fastq.gz \
  --taxid FAMILY_TAXID \
  --output extraction/

# Step 3: Build BLAST database from reference genomes
# Note: Headers are standardized to {accession}|{contig_id} format
# Creates contig_mapping.tsv for multi-contig draft genome tracking
metadarkmatter blast makedb \
  --genomes genomes/ \
  --output blastdb/pangenome

# Step 4: Run competitive BLAST alignment (accepts compressed query)
metadarkmatter blast align \
  --query extraction/sample_taxidFAMILY_TAXID_R1.fastq.gz \
  --database blastdb/pangenome \
  --output sample.blast.tsv.gz \
  --threads 16

# Step 5: Compute ANI matrix from reference genomes
metadarkmatter ani compute \
  --genomes genomes/ \
  --output ani_matrix.csv \
  --threads 16

# Step 6: ANI-weighted classification (core algorithm)
# Use --metadata to add species/genus columns and species-level aggregation
metadarkmatter score classify \
  --alignment sample.blast.tsv.gz \
  --ani ani_matrix.csv \
  --metadata genome_metadata.tsv \
  --output classifications.csv \
  --summary summary.json \
  --parallel

# Step 7: Generate HTML report with species breakdown
metadarkmatter report generate \
  --classifications classifications.csv \
  --metadata genome_metadata.tsv \
  --output report.html

# Optional: Multi-sample comparison
metadarkmatter report multi \
  --input-dir results/ \
  --output comparison.html
```

### Protein-Level Workflow

For highly divergent taxa where nucleotide-based alignment may miss homology, use protein-level (BLASTX) classification:

```bash
# Step 3b (Alternative): Build Diamond protein database from reference proteomes
metadarkmatter blastx makedb \
  --proteins genomes/proteins/ \
  --output blastdb/panproteome \
  --threads 16

# Step 4b (Alternative): Run BLASTX alignment (DNA reads vs protein database)
metadarkmatter blastx align \
  --query extraction/sample_taxidFAMILY_TAXID_R1.fastq.gz \
  --database blastdb/panproteome \
  --output sample.blastx.tsv.gz \
  --threads 16

# Step 6b (Alternative): Protein-level classification (uses protein thresholds)
metadarkmatter score classify \
  --alignment sample.blastx.tsv.gz \
  --ani ani_matrix.csv \
  --aai aai_matrix.csv \
  --alignment-mode protein \
  --output classifications.csv \
  --summary summary.json \
  --parallel
```

**When to use protein-level classification:**
- Genus-level novelty where nucleotide identity degrades faster than protein identity
- Samples with highly divergent sequences (>25% nucleotide divergence)
- When BLASTN fails to detect homology that BLASTX can capture

---

## Command Reference

| Command | Purpose |
|---------|---------|
| `download genomes list` | Query GTDB API for family genomes |
| `download genomes fetch` | Download genomes from NCBI (use `--include-protein` for AAI) |
| `proteins predict` | Predict proteins with Prodigal (use `--missing-only` for gaps) |
| `kraken2 classify` | Run Kraken2 classification |
| `kraken2 extract` | Extract reads for target taxid (compressed output) |
| `blast makedb` | Build BLAST nucleotide database |
| `blast align` | BLASTN competitive alignment (accepts .gz query) |
| `blastx makedb` | Build Diamond protein database |
| `blastx align` | BLASTX alignment (DNA query vs protein database) |
| `ani compute` | Compute ANI matrix (auto-detects skani/fastANI) |
| `ani validate` | Validate ANI matrix coverage against BLAST results |
| `aai compute` | Compute AAI matrix using Diamond BLASTP |
| `aai validate` | Validate AAI matrix coverage against ANI matrix |
| `map reads` | Bowtie2 competitive mapping |
| `score classify` | ANI-weighted classification (supports `--alignment-mode protein`) |
| `score batch` | Batch process multiple samples |
| `score extract-novel` | Extract candidate novel species/genera |
| `visualize recruitment` | Recruitment plots |
| `visualize summary` | Summary charts |
| `report generate` | Single-sample HTML report |
| `report multi` | Multi-sample comparison |

---

## Core Algorithm

**ANI-Weighted Placement Uncertainty:**

1. For each read, find top BLAST hit (max bitscore)
2. Find ambiguous hits (within 95% of top bitscore)
3. Calculate **Novelty Index**: `N = 100 - TopHitIdentity`
4. Calculate **Placement Uncertainty**: `U = 100 - max(ANI between competing genomes)`
5. Classify based on thresholds (mode-dependent, literature-backed)

**For comprehensive methods documentation with all calculations at scientific standard, see [METHODS.md](METHODS.md).**

**For step-by-step algorithm explanation with examples, see [ALGORITHM_DETAILED.md](ALGORITHM_DETAILED.md).**

### Nucleotide Mode Thresholds

Default mode for BLASTN alignment:

| Category | Novelty (N) | Uncertainty (U) | Interpretation |
|----------|-------------|-----------------|----------------|
| Known Species | N < 4% | U < 1.5% | >96% identity, confident placement |
| **Novel Species** | 4% <= N < 20% | U < 1.5% | 80-96% identity, potential new species |
| **Novel Genus** | 20% <= N <= 25% | U < 1.5% | 75-80% identity, potential new genus |
| Ambiguous | any | 1.5% <= U < 5% | Species boundary zone (95-98.5% ANI) |
| Conserved Region | any | U >= 5% | Conserved gene across genera (<95% ANI) |
| Unclassified | N > 25% or edge | U < 1.5% | Requires manual review |

**Threshold basis:** 96% ANI = prokaryotic species boundary (Jain et al. 2018, Nature Communications). See [METHODS.md](METHODS.md) Section 2.4 for detailed justification.

### Protein Mode Thresholds

For BLASTX alignment (`--alignment-mode protein`):

| Category | Novelty (N) | Uncertainty (U) | Interpretation |
|----------|-------------|-----------------|----------------|
| Known Species | N < 10% | U < 5% | >90% protein identity, confident placement |
| **Novel Species** | 10% <= N < 25% | U < 5% | 75-90% identity, potential new species |
| **Novel Genus** | 25% <= N <= 40% | U < 5% | 60-75% identity, potential new genus |
| Ambiguous | any | 5% <= U < 10% | Species boundary zone |
| Conserved Region | any | U >= 10% | Highly conserved protein region |
| Unclassified | N > 40% or edge | U < 5% | Requires manual review |

Protein thresholds are calibrated based on the observation that amino acid sequences diverge more slowly than nucleotide sequences.

See `docs/CLASSIFICATION_STATISTICS.md` for detailed statistical framework and literature references.

---

## GTDB Compatibility

metadarkmatter aligns with GTDB (Genome Taxonomy Database) standards for species delineation.

### Preset Configurations

```bash
# GTDB-strict: 95% ANI, requires 50% alignment fraction
metadarkmatter score classify --alignment results.tsv --ani ani.csv \
    --preset gtdb-strict --output classifications.csv

# Literature-strict: 96% ANI, 1.5% uncertainty, 60% confidence threshold
# Based on comprehensive literature review (Jain 2018, Riesco 2024)
metadarkmatter score classify --alignment results.tsv --ani ani.csv \
    --preset literature-strict --output classifications.csv

# Available presets: gtdb-strict, gtdb-relaxed, conservative, literature-strict, coverage-strict, gtdb-coverage
```

### Alignment Quality Filters

```bash
# Custom alignment filters (GTDB requires AF >= 50%)
metadarkmatter score classify --alignment results.tsv --ani ani.csv \
    --min-alignment-length 100 \
    --min-alignment-fraction 0.5 \
    --output classifications.csv
```

### Output Metrics

The classification output includes additional metrics for quality assessment:

| Metric | Description |
|--------|-------------|
| `confidence_score` | 0-100 score integrating margin from thresholds, placement certainty, and alignment quality. Higher values (80+) indicate classifications far from decision boundaries. |
| `genus_uncertainty` | Based on minimum ANI to secondary genomes (low = same genus, high = multiple genera) |
| `ambiguity_scope` | `unambiguous`, `within_species`, `within_genus`, or `across_genera` |
| `num_competing_genera` | Count of secondary genomes from different genera |

**Confidence Score Components:**
- Margin from threshold boundaries (0-40 pts): Distance from classification cutoffs
- Placement certainty (0-40 pts): Hit count and identity gap to secondary genomes
- Alignment quality proxy (0-20 pts): Derived from top hit identity

### Novel Taxon Naming

When using metadata, novel taxa receive GTDB-style placeholder names:

```python
from metadarkmatter.core.metadata import GenomeMetadata

metadata = GenomeMetadata.from_file(Path("genome_metadata.tsv"))
enriched = metadata.join_classifications_with_novel_names(classifications)
# Adds 'suggested_name' column:
# - Novel Species: "Francisella sp. MDM-A7X"
# - Novel Genus: "Francisellaceae gen. nov. MDM-B3K"
```

---

## Database Strategy

The BLAST database choice affects results interpretation.

| Strategy | Database | Use Case |
|----------|----------|----------|
| **Family DB** | Target family genomes only | Fast screening, most analyses |
| **GTDB DB** | All representative bacteria | Validation, discovery |
| **Hybrid** | Family DB + GTDB validation | Publication-quality results |

**Key insight:** Reads misclassified by Kraken2 will appear as "highly novel" when using a family-only database. For publication, validate novel hits against GTDB to rule out misclassification.

See `docs/WORKFLOW.md` for detailed guidance on database strategy.

---

## Project Structure

```
metadarkmatter/
├── src/metadarkmatter/
│   ├── cli/              # CLI commands (typer)
│   │   ├── main.py       # Entry point, app definition
│   │   ├── score.py      # ANI-weighted classification (supports --alignment-mode)
│   │   ├── kraken2.py    # Kraken2 classify/extract
│   │   ├── blast.py      # BLAST makedb/align
│   │   ├── blastx.py     # Diamond BLASTX makedb/align (protein-level)
│   │   ├── ani.py        # ANI compute/validate
│   │   ├── map.py        # Bowtie2 mapping
│   │   ├── download.py   # Genome acquisition (saves metadata)
│   │   ├── visualize.py  # Recruitment/summary plots
│   │   ├── report.py     # HTML reports (with species tab)
│   │   └── utils.py      # Shared utilities (progress, sample names)
│   ├── core/             # Core algorithms
│   │   ├── ani_placement.py    # Main classifier (ANIMatrix, classifiers)
│   │   ├── ani_matrix_builder.py  # ANI matrix construction
│   │   ├── parsers.py          # Streaming BLAST/ANI parsers
│   │   ├── recruitment.py      # BAM recruitment data extraction
│   │   ├── genome_utils.py     # Genome concatenation, header rewriting
│   │   ├── metadata.py         # Species/genus metadata handling
│   │   ├── constants.py        # Centralized nucleotide constants and thresholds
│   │   ├── protein_constants.py  # Protein-specific thresholds for BLASTX
│   │   ├── exceptions.py       # Custom exception classes
│   │   └── io_utils.py         # DataFrame I/O utilities
│   ├── external/         # External tool wrappers
│   │   ├── base.py       # ExternalTool base class, path validation
│   │   ├── blast.py      # BLAST+ (makeblastdb, blastn)
│   │   ├── diamond.py    # Diamond (makedb, blastp, blastx)
│   │   ├── bowtie2.py    # Bowtie2 aligner
│   │   ├── kraken.py     # Kraken2 classifier
│   │   ├── fastani.py    # FastANI wrapper
│   │   ├── skani.py      # skani wrapper
│   │   ├── samtools.py   # Samtools wrapper
│   │   └── ncbi_datasets.py  # NCBI datasets CLI
│   ├── models/           # Pydantic data models
│   │   ├── blast.py      # BlastHit, BlastResult
│   │   ├── classification.py   # ReadClassification, TaxonomicSummary
│   │   ├── config.py     # ScoringConfig
│   │   └── genomes.py    # AccessionList, GTDBGenome
│   ├── clients/          # API clients
│   │   └── gtdb.py       # GTDB API client (with retry logic)
│   ├── visualization/    # Plotly charts and reports
│   │   ├── plots/        # Chart components
│   │   │   ├── base.py   # BasePlot class
│   │   │   ├── distributions.py    # Histograms
│   │   │   ├── scatter_2d.py       # Scatter plots
│   │   │   ├── classification_charts.py  # Donut/bar charts
│   │   │   └── multi_sample.py     # Multi-sample comparisons
│   │   ├── recruitment_plots.py    # Recruitment visualizations
│   │   └── report/       # HTML report generation
│   │       ├── generator.py      # Single-sample reports
│   │       ├── multi_generator.py # Multi-sample reports
│   │       ├── styles.py         # CSS themes
│   │       └── templates.py      # HTML templates
│   ├── io/               # I/O utilities (placeholder)
│   └── utils/            # General utilities (placeholder)
├── tests/                # pytest test suite
│   ├── unit/             # Unit tests
│   ├── e2e/              # End-to-end tests
│   ├── integration/      # Integration tests (require external tools)
│   ├── conftest.py       # pytest fixtures
│   └── factories.py      # Test data factories
├── docs/                 # Documentation
└── config/               # Default configuration
```

---

## Advanced Features

### Coverage-Weighted Hit Selection

By default, coverage-weighted scoring prioritizes hits that cover more of the query sequence, reducing bias from short conserved domains (e.g., 16S rRNA fragments) that may align with high identity but low coverage. Linear coverage weighting is enabled by default (`--coverage-weight-mode linear`).

**When to Use:**
- Reads contain conserved domains that align with high identity but low coverage
- Short high-identity hits (50bp) are outranking longer moderate-identity hits (280bp)
- You want to weight alignment quality by how much of the read was aligned

**Basic Usage:**
```bash
# Enable coverage weighting with linear mode
metadarkmatter score classify --alignment sample.blast.tsv.gz --ani ani.csv \
    --coverage-weight-mode linear --output classifications.csv

# Use a preset with coverage weighting
metadarkmatter score classify --alignment sample.blast.tsv.gz --ani ani.csv \
    --preset coverage-linear --output classifications.csv
```

**Coverage Weighting Modes:**

| Mode | Formula | Behavior |
|------|---------|----------|
| `none` | `weighted_score = bitscore` | Disables coverage weighting |
| `linear` | `weight = min + (max - min) × coverage` | Default; gradual penalty for low coverage |
| `log` | `weight = min + (max - min) × log(1 + 9×coverage) / log(10)` | Rewards any coverage, diminishing returns |
| `sigmoid` | `weight = min + (max - min) / (1 + exp(-10×(coverage - 0.6)))` | Sharp threshold around 60% coverage |

Where:
- `coverage = (qend - qstart + 1) / qlen`
- `min = 1.0 - strength`, `max = 1.0 + strength`
- `strength` controls magnitude (default 0.5, range 0.0-1.0)

**Example:**

```
Read length: 300bp
Hit A: 98% identity, 50bp aligned, bitscore=200 (coverage=0.17)
Hit B: 92% identity, 280bp aligned, bitscore=180 (coverage=0.93)

With mode="none" (default): Hit A wins (bitscore 200 > 180)
With mode="linear", strength=0.5:
  - Hit A weighted: 200 × 0.58 = 116
  - Hit B weighted: 180 × 0.97 = 174
  Hit B wins (weighted 174 > 116)
```

**Coverage Weighting Presets:**

| Preset | Mode | Strength | Description |
|--------|------|----------|-------------|
| `coverage-linear` | linear | 0.5 | Balanced coverage weighting |
| `coverage-strict` | sigmoid | 0.7 | Enforces >60% coverage threshold |
| `coverage-gentle` | log | 0.3 | Mild preference for higher coverage |
| `gtdb-coverage` | linear | 0.5 | GTDB-style with 50% alignment fraction |

**CLI Parameters:**

```bash
# Explicit parameters
metadarkmatter score classify --alignment sample.blast.tsv.gz --ani ani.csv \
    --coverage-weight-mode linear \
    --coverage-weight-strength 0.5 \
    --output classifications.csv
```

**Note:** The `--alignment` parameter accepts both BLAST and MMseqs2 tabular output (identical 13-column format with qlen).

### Quality Control

Pre- and post-classification QC metrics help identify problematic inputs and flag unreliable results:

```bash
# Save QC metrics as JSON
metadarkmatter score classify --alignment sample.blast.tsv.gz --ani ani.csv \
    --qc-output qc_metrics.json --output classifications.csv --parallel
```

QC warnings are generated when:
- Filter rate exceeds 50%
- Genome coverage below 50%
- Single-hit fraction above 80%
- Mean identity below 80%
- Ambiguous fraction above 50%

### Adaptive Thresholds

Detect the species boundary from the ANI matrix distribution instead of using the fixed 96% default:

```bash
# Enable adaptive threshold detection
metadarkmatter score classify --alignment sample.blast.tsv.gz --ani ani.csv \
    --adaptive-thresholds --output classifications.csv --parallel
```

This fits a 2-component Gaussian Mixture Model to pairwise ANI values, detecting natural within-species vs. between-species clusters. Falls back to the default threshold if the GMM does not converge or components are not well separated.

Requires: `pip install metadarkmatter[adaptive]` (installs scikit-learn).

See [METHODS.md](METHODS.md) Section 9 for the mathematical details.

### Bayesian Confidence

Add posterior probabilities to classification output, providing continuous confidence scores near threshold boundaries:

```bash
# Add Bayesian posterior columns
metadarkmatter score classify --alignment sample.blast.tsv.gz --ani ani.csv \
    --bayesian --output classifications.csv --parallel
```

This adds 6 columns to the output:
- `p_known_species`, `p_novel_species`, `p_novel_genus`, `p_ambiguous` (sum to 1.0)
- `bayesian_category` (MAP classification)
- `posterior_entropy` (0 = confident, 2.0 = uniform across all categories)

The HTML report will include a Bayesian Confidence tab with entropy distribution, posterior bar charts, and a confidence landscape plot.

See [METHODS.md](METHODS.md) Section 10 for the likelihood model and interpretation.

### Threshold Sensitivity Analysis

Assess whether classification results are robust to threshold choice:

```bash
# Run sensitivity analysis
metadarkmatter score sensitivity --alignment sample.blast.tsv.gz --ani ani.csv \
    --output sensitivity.json --parallel
```

Sweeps the novelty/uncertainty thresholds across a range and reports how classification counts change at each point.

### Combined Advanced Options

All new features can be combined:

```bash
metadarkmatter score classify --alignment sample.blast.tsv.gz --ani ani.csv \
    --adaptive-thresholds \
    --bayesian \
    --qc-output qc_metrics.json \
    --output classifications.csv \
    --summary summary.json \
    --parallel --verbose
```

### Performance Modes

For large datasets, use appropriate processing mode:

```bash
# Standard mode (< 1M reads)
metadarkmatter score classify --alignment sample.blast.tsv.gz --ani ani.csv --output out.csv

# Fast mode (1-10M reads)
metadarkmatter score classify --alignment sample.blast.tsv.gz --ani ani.csv --output out.csv --fast

# Parallel mode (10-100M reads) - recommended
metadarkmatter score classify --alignment sample.blast.tsv.gz --ani ani.csv --output out.csv --parallel

# Streaming mode (100M+ reads, memory-constrained)
metadarkmatter score classify --alignment sample.blast.tsv.gz --ani ani.csv --output out.csv --streaming
```

### ANI Matrix Computation

The `ani compute` command supports multiple backends with auto-detection:

```bash
# Auto-detect best available backend (prefers skani > fastANI)
metadarkmatter ani compute --genomes genomes/ --output ani_matrix.csv --threads 16

# Explicitly use skani (faster)
metadarkmatter ani compute --genomes genomes/ --output ani_matrix.csv --backend skani

# Explicitly use fastANI
metadarkmatter ani compute --genomes genomes/ --output ani_matrix.csv --backend fastani

# Validate ANI matrix coverage against BLAST results
metadarkmatter ani validate --ani ani_matrix.csv --alignment sample.blast.tsv.gz
```

Backend comparison:
- **skani**: Faster, memory-efficient, good for large genome sets (preferred)
- **fastANI**: Original tool, well-validated, widely used in publications

### Novel Species Extraction

After classification, use the `score extract-novel` command to identify reads representing novel diversity:

```bash
# Extract all novel candidates (species + genus) grouped by reference genome
metadarkmatter score extract-novel \
  --classifications classifications.csv \
  --output novel_candidates.csv \
  --read-ids novel_read_ids.txt

# Extract only novel species candidates
metadarkmatter score extract-novel \
  --classifications classifications.csv \
  --output novel_species.csv \
  --category species

# Extract highly novel reads (custom threshold)
metadarkmatter score extract-novel \
  --classifications classifications.csv \
  --output highly_novel.csv \
  --min-novelty 20

# Use read IDs with seqtk for targeted assembly
seqtk subseq extraction/sample_R1.fastq.gz novel_read_ids.txt > novel_reads.fastq
```

Alternatively, use Python for custom filtering:

```python
import polars as pl

df = pl.read_csv("classifications.csv")

# Filter and group by genome
candidates = df.filter(
    pl.col("taxonomic_call").is_in(["Novel Species", "Novel Genus"])
).group_by("best_match_genome").agg([
    pl.len().alias("read_count"),
    pl.col("novelty_index").mean().alias("mean_novelty"),
]).sort("read_count", descending=True)

print(candidates)
```

### Species-Level Tracking

metadarkmatter tracks species-level information throughout the pipeline.

#### Metadata Files

| File | Created By | Contents |
|------|------------|----------|
| `genome_metadata.tsv` | `download genomes list` | accession, species, genus, family, gtdb_taxonomy |
| `contig_mapping.tsv` | `blast makedb` | contig_id, genome_accession, original_header |

#### Standardized FASTA Headers

When building databases, FASTA headers are rewritten to a standardized format:

```
>{accession}|{original_contig_id}
```

This enables reliable genome identification from BLAST hits, especially for multi-contig draft genomes.

#### Using Metadata in Classification

```bash
# Classification with species-level aggregation
metadarkmatter score classify \
  --alignment sample.blast.tsv.gz \
  --ani ani_matrix.csv \
  --metadata genome_metadata.tsv \
  --output classifications.csv \
  --summary summary.json

# Output CSV includes species and genus columns
# Summary JSON includes species_hit_counts
```

#### Using Metadata in Reports

```bash
# Report with species breakdown tab
metadarkmatter report generate \
  --classifications classifications.csv \
  --metadata genome_metadata.tsv \
  --output report.html
```

The report includes a "Species Breakdown" tab with:
- Species read count bar chart (top 20)
- Species composition pie chart (top 10)
- Unique species count vs reference species count

#### Python API for Metadata

```python
from metadarkmatter.core.metadata import GenomeMetadata
import polars as pl

# Load metadata
metadata = GenomeMetadata.from_file(Path("genome_metadata.tsv"))

# Join with classifications
df = pl.read_csv("classifications.csv")
enriched = metadata.join_classifications(df)
# enriched now has 'species' and 'genus' columns

# Aggregate by species
species_stats = metadata.aggregate_by_species(enriched)
print(species_stats)
# Columns: species, read_count, mean_novelty, mean_identity, mean_uncertainty, genome_count
```

---

## Development Guide

### Running Tests

```bash
# Run all unit tests
pytest tests/ --ignore=tests/integration/ -v

# Run with coverage
pytest tests/ --ignore=tests/integration/ --cov=metadarkmatter

# Run specific test file
pytest tests/unit/test_ani_placement.py -v
```

### Code Quality

```bash
# Lint check
ruff check src/metadarkmatter/

# Format check
ruff format --check src/metadarkmatter/

# Type checking (if mypy configured)
mypy src/metadarkmatter/
```

### Code Quality and Security

#### Security Hardening

The codebase includes several security measures:

- **Path validation**: All file paths passed to subprocess calls are validated via `validate_path_safe()` in `external/base.py`. This checks for null bytes, path traversal patterns, and unusual characters.
- **Zip Slip protection**: NCBI dataset extraction validates that extracted files stay within the output directory.
- **Input validation**: CLI parameters include bounds checking (e.g., `--taxid` must be >= 1, chunk sizes have min/max bounds).
- **Accession validation**: NCBI accessions are validated against expected patterns before use.

#### API Client Resilience

- **GTDB client**: Implements exponential backoff retry logic for transient failures (5xx errors, rate limiting, connection errors). Configurable via `max_retries`, `retry_delay`, `retry_backoff` parameters.
- **NCBI downloads**: Rate limiting between batch requests (default 0.5s delay) to avoid overwhelming the API.

#### Constants Module

Centralized constants are defined in `core/constants.py` (nucleotide) and `core/protein_constants.py` (protein):

```python
# Nucleotide thresholds (default mode)
from metadarkmatter.core.constants import (
    UNKNOWN_GENOME,           # Default placeholder for unknown genomes
    BLAST_OUTFMT_12COL,       # Standard BLAST output format string
    ANI_SPECIES_BOUNDARY_LOW, # 95% ANI threshold
    CATEGORY_NOVEL_SPECIES,   # Classification category strings
    NOVELTY_KNOWN_MAX,        # 5% novelty threshold
    UNCERTAINTY_CONFIDENT_MAX,# 2% uncertainty threshold
    calculate_confidence_score,  # Margin-based confidence calculation
)

# Protein thresholds (for --alignment-mode protein)
from metadarkmatter.core.protein_constants import (
    PROTEIN_NOVELTY_KNOWN_MAX,           # 10% novelty threshold
    PROTEIN_NOVELTY_NOVEL_SPECIES_MIN,   # 10% novelty threshold
    PROTEIN_NOVELTY_NOVEL_SPECIES_MAX,   # 25% novelty threshold
    PROTEIN_NOVELTY_NOVEL_GENUS_MIN,     # 25% novelty threshold
    PROTEIN_NOVELTY_NOVEL_GENUS_MAX,     # 40% novelty threshold
    PROTEIN_UNCERTAINTY_CONFIDENT_MAX,   # 5% uncertainty threshold
    PROTEIN_UNCERTAINTY_CONSERVED_MIN,   # 10% uncertainty threshold
    calculate_protein_confidence_score,  # Protein-specific confidence
)
```

#### Error Handling

CLI commands use specific exception handlers:
- `PolarsError`: Data processing issues (malformed files)
- `FileNotFoundError`: Missing input files
- `PermissionError`: File access issues
- `MemoryError`: Suggests using `--streaming` mode

#### Type Hints

Complex data structures use type aliases defined in `core/ani_placement.py`:

```python
HitData: TypeAlias = tuple[str, str, float, float, str]  # BLAST hit tuple
ReadChunk: TypeAlias = tuple[str, list[HitData]]         # Read with hits
ClassificationDict: TypeAlias = dict[str, Any]           # Classification result
```

### Dependency Injection for Testing

External tool wrappers support dependency injection for testing:

```python
from metadarkmatter.external.base import ExternalTool

# Mock executable resolver for tests
def mock_resolver(name: str) -> str | None:
    return f"/usr/bin/{name}"

ExternalTool._executable_resolver = mock_resolver
```
