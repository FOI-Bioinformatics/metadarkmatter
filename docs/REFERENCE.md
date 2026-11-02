# metadarkmatter Reference Guide

Complete technical reference for the metadarkmatter package.

## Table of Contents

- [CLI Workflows](#cli-workflows)
  - [Command Summary](#command-summary)
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

For complete step-by-step workflows, see the
[Tutorial](TUTORIAL_ENVIRONMENTAL_SPECIES.md) or the
[Workflow Guide](WORKFLOW.md).

### Command Summary

| Step | Nucleotide Mode | Protein Mode |
|------|----------------|--------------|
| 1. Download genomes | `download genomes list` + `fetch` | Same |
| 2. Extract reads | `kraken2 classify` + `extract` | Same |
| 3. Build database | `blast makedb` or `mmseqs2 makedb` | `blastx makedb` |
| 4. Align reads | `blast align` or `mmseqs2 search` | `blastx align` |
| 5. Compute matrix | `ani compute` | `ani compute` + `aai compute` |
| 6. Classify | `score classify` | `score classify --alignment-mode protein` |
| 7. Report | `report generate` | Same |

### Protein-Level Workflow

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
| `kraken2 extract` | Extract reads for target taxid (gzip output by default) |
| `blast makedb` | Build BLAST nucleotide database |
| `blast align` | BLASTN competitive alignment (accepts .gz query) |
| `mmseqs2 makedb` | Build MMseqs2 nucleotide database |
| `mmseqs2 search` | MMseqs2 search with BLAST-compatible output (single- or paired-end) |
| `blastx makedb` | Build Diamond protein database |
| `blastx align` | BLASTX alignment (DNA query vs protein database) |
| `ani compute` | Compute ANI matrix (auto prefers skani, falls back to fastANI) |
| `ani validate` | Validate ANI matrix coverage against BLAST results |
| `aai compute` | Compute AAI matrix using Diamond BLASTP |
| `aai validate` | Validate AAI matrix coverage against ANI matrix |
| `map reads` | Bowtie2 competitive mapping |
| `score classify` | ANI-weighted classification (supports `--alignment-mode protein`) |
| `score batch` | Batch process multiple samples |
| `score extract-novel` | Extract candidate novel species/genera |
| `score export-config` | Write an editable YAML configuration |
| `score sensitivity` | Sweep hard thresholds on a saved classification file |
| `score evaluate` | Accuracy, per-category precision/recall, ECE against a truth column |
| `util generate-mapping` | Generate contig-to-genome ID mapping from FASTA directory |
| `util validate-mapping` | Validate ID mapping file (optionally against BLAST results) |
| `util generate-metadata` | Build genome_metadata.tsv from a genome directory |
| `util select-representatives` | Assign species representatives in metadata |
| `tree build` | Build a Newick tree (nj, upgma, mashtree) |
| `validate ani/blast/classifications/fasta/fastq` | Input/output format checks |
| `visualize recruitment` | Recruitment plots |
| `visualize summary` | Summary charts |
| `report generate` | Single-sample HTML report |
| `report multi` | Multi-sample comparison |
| `report summarize` | Summarize novel diversity into putative taxon clusters |
| `report vendor-d3` / `vendor-cdn-sri` | Cache D3 locally / compute SRI hashes for CDN tags |
| `run pipeline` | Run the full pipeline end to end |
| `doctor` | Report environment, dependencies and tool versions |

Option-level details for every command are in [CLI_REFERENCE.md](CLI_REFERENCE.md).

---

## Core Algorithm

**ANI-Weighted Placement Uncertainty:**

1. For each read, find top BLAST hit (max coverage-weighted bitscore; linear weighting by default)
2. Find ambiguous hits (within 95% of top bitscore)
3. Calculate **Novelty Index**: `N = 100 - TopHitIdentity`
4. Calculate **Placement Uncertainty**: `U = 100 - ANI(best genome, second-best genome)` with the default `--uncertainty-mode second`; `--uncertainty-mode max` uses `U = 100 - max(ANI between the best genome and competing genomes)`
5. Classify: the primary `taxonomic_call` is the Bayesian MAP category over (N, U) with a discrete refinement step; the mode-dependent hard thresholds below are reported as `legacy_call`

**For comprehensive methods documentation with all calculations at scientific standard, see [METHODS.md](METHODS.md).**

**For step-by-step algorithm explanation with examples, see [ALGORITHM_DETAILED.md](ALGORITHM_DETAILED.md).**

### Nucleotide Mode Thresholds

Default mode for BLASTN/MMseqs2 alignment. These are the hard thresholds applied in order to produce `legacy_call` (`core/classification/thresholds.py`):

| Category | Novelty (N) | Uncertainty (U) | Interpretation |
|----------|-------------|-----------------|----------------|
| Ambiguous | any | identity gap < 2% with >1 ambiguous hit | Near-tied hits to different genomes |
| Conserved Region | any | U >= 5% and hits span genera | Conserved gene across genera (<95% ANI) |
| Ambiguous | any | U >= 5% (within genus) | High placement ambiguity |
| Species Boundary | any | 1.5% <= U < 5% | Species boundary zone (95-98.5% ANI) |
| Known Species | N < 4% | U < 1.5% | >96% identity, confident placement |
| Ambiguous | N >= 4% | single hit, inferred U >= 10% | Single-hit read, placement not measurable |
| **Novel Species** | 4% <= N < 20% | U < 1.5% | 80-96% identity, potential new species |
| **Novel Genus** | 20% <= N <= 25% | U < 1.5% | 75-80% identity, potential new genus |
| Unclassified | N > 25% or edge | U < 1.5% | Requires manual review |

**Relation to `taxonomic_call`:** the primary call is the Bayesian MAP category (Gaussian likelihoods centred on these regions), so its decision boundaries do not coincide exactly with the table. With default parameters at U = 0, the measured transitions are approximately N = 6.2 (Known -> Novel Species), N = 18.6 (Novel Species -> Novel Genus) and N = 27.6 (Novel Genus -> Unclassified); at N = 1, Species Boundary becomes the MAP call from about U = 3.8 and Ambiguous from about U = 5.5. Single-hit reads use inferred uncertainty only in `legacy_call`. Stage 2 refinement relabels Ambiguous reads with cross-genus hits as Conserved Region and some Novel Genus reads as Ambiguous Within Genus.

**Threshold basis:** 96% ANI = prokaryotic species boundary (Jain et al. 2018, Nature Communications). See [METHODS.md](METHODS.md) Section 2.4 for detailed justification.

### Protein Mode Thresholds

For BLASTX alignment (`--alignment-mode protein`):

| Category | Novelty (N) | Uncertainty (U) | Interpretation |
|----------|-------------|-----------------|----------------|
| Known Species | N < 10% | U < 5% | >90% protein identity, confident placement |
| **Novel Species** | 10% <= N < 25% | U < 5% | 75-90% identity, potential new species |
| **Novel Genus** | 25% <= N <= 40% | U < 5% | 60-75% identity, potential new genus |
| Species Boundary | any | 5% <= U < 10% | Species boundary zone |
| Conserved Region / Ambiguous | any | U >= 10% | Conserved Region when hits span genera, otherwise Ambiguous |
| Unclassified | N > 40% or edge | U < 5% | Requires manual review |

As in nucleotide mode, these thresholds drive `legacy_call`; the same identity-gap and single-hit rules apply. In protein mode the fixed values above are used regardless of `--preset`, `--config` or `--adaptive-thresholds`. The Bayesian `taxonomic_call` novelty transitions at U = 0 are measured at about N = 11.3, 25.1 and 40.4.

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

# Available presets: default, gtdb-strict, gtdb-relaxed, conservative, literature-strict, coverage-strict, adaptive
# Explicitly passed flags override preset values.
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
| `confidence_score` | 0-100, derived from posterior entropy: `(1 - posterior_entropy / log2(6)) * 100`. 10.0 for Off-target reads. |
| `genus_uncertainty` | Based on minimum ANI to secondary genomes (low = same genus, high = multiple genera) |
| `ambiguity_scope` | `unambiguous`, `within_species`, `within_genus`, or `across_genera` |

The older margin-based score (`calculate_confidence_score()` in `core/constants.py`: threshold margin, placement certainty, alignment quality) is no longer the `confidence_score` column; its sub-scores are available with `--include-legacy-scores`.

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
│   │   ├── main.py       # Entry point, app definition, global options
│   │   ├── score.py      # classify, batch, extract-novel, export-config, sensitivity, evaluate
│   │   ├── kraken2.py    # Kraken2 classify/extract
│   │   ├── blast.py      # BLAST makedb/align
│   │   ├── blastx.py     # Diamond BLASTX makedb/align (protein-level)
│   │   ├── mmseqs2.py    # MMseqs2 makedb/search
│   │   ├── ani.py        # ANI compute/validate
│   │   ├── aai.py        # AAI compute/validate
│   │   ├── map.py        # Bowtie2 mapping
│   │   ├── proteins.py   # Prodigal protein prediction
│   │   ├── download.py   # Genome acquisition (saves metadata)
│   │   ├── visualize.py  # Recruitment/summary plots
│   │   ├── report.py     # HTML reports (generate, multi, summarize, vendor-*)
│   │   ├── mapping.py    # 'util' commands (generate-mapping, validate-mapping, generate-metadata, select-representatives)
│   │   ├── tree.py       # tree build
│   │   ├── validate.py   # validate ani/blast/classifications/fasta/fastq
│   │   ├── doctor.py     # Environment report
│   │   ├── run.py        # run pipeline
│   │   └── utils.py      # Shared utilities (progress, sample names)
│   ├── core/             # Core algorithms
│   │   ├── classification/     # ANIMatrix, classifiers, bayesian.py, thresholds.py, qc.py, adaptive.py, sensitivity.py, evaluation.py
│   │   ├── ani_placement.py    # Backward-compatible re-exports of the classification package
│   │   ├── ani_matrix_builder.py  # ANI matrix construction
│   │   ├── aai_matrix_builder.py  # AAI matrix construction
│   │   ├── parsers.py          # Streaming BLAST/ANI parsers
│   │   ├── recruitment.py      # BAM recruitment data extraction
│   │   ├── genome_utils.py     # Genome concatenation, header rewriting
│   │   ├── id_mapping.py       # Contig-to-genome ID mapping for external results
│   │   ├── metadata.py         # Species/genus metadata handling
│   │   ├── novel_diversity/    # Novel read clustering
│   │   ├── phylogeny/          # Tree building and novel cluster placement
│   │   ├── pipeline/           # 'run pipeline' orchestrator
│   │   ├── constants.py        # Centralized nucleotide constants and thresholds
│   │   ├── protein_constants.py  # Protein-specific thresholds for BLASTX
│   │   ├── exceptions.py       # Custom exception classes
│   │   └── io_utils.py         # DataFrame I/O utilities
│   ├── external/         # External tool wrappers
│   │   ├── base.py       # ExternalTool base class, path validation
│   │   ├── blast.py      # BLAST+ (makeblastdb, blastn)
│   │   ├── diamond.py    # Diamond (makedb, blastp, blastx)
│   │   ├── mmseqs2.py    # MMseqs2
│   │   ├── bowtie2.py    # Bowtie2 aligner
│   │   ├── kraken.py     # Kraken2 classifier
│   │   ├── fastani.py    # FastANI wrapper
│   │   ├── skani.py      # skani wrapper
│   │   ├── prodigal.py   # Prodigal wrapper
│   │   ├── mashtree.py   # Mashtree wrapper
│   │   ├── samtools.py   # Samtools wrapper
│   │   └── ncbi_datasets.py  # NCBI datasets CLI
│   ├── models/           # Pydantic data models
│   │   ├── blast.py      # BlastHit, BlastResult
│   │   ├── classification.py   # ReadClassification, TaxonomicSummary
│   │   ├── config.py     # ScoringConfig
│   │   ├── pipeline_config.py  # 'run pipeline' config
│   │   └── genomes.py    # AccessionList, GTDBGenome
│   ├── clients/          # API clients
│   │   └── gtdb.py       # GTDB API client (with retry logic)
│   └── visualization/    # Plotly charts and reports
│       ├── plots/        # Chart components
│       ├── recruitment_plots.py    # Recruitment visualizations
│       └── report/       # HTML report generation (generator.py facade + sections/ mixins)
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

# Use a preset with stronger (sigmoid) coverage weighting
metadarkmatter score classify --alignment sample.blast.tsv.gz --ani ani.csv \
    --preset coverage-strict --output classifications.csv
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

With mode="none": Hit A wins (bitscore 200 > 180)
With mode="linear" (default), strength=0.5:
  - Hit A weighted: 200 × 0.67 = 133
  - Hit B weighted: 180 × 1.43 = 258
  Hit B wins (weighted 258 > 133)
```

**Coverage Weighting Presets:**

| Preset | Mode | Strength | Description |
|--------|------|----------|-------------|
| `coverage-strict` | sigmoid | 0.7 | Sharp preference for >60% coverage (also narrows the novel species range to N < 15%) |

All other presets use the default linear mode with strength 0.5; other modes are selected with `--coverage-weight-mode`.

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
    --qc-output qc_metrics.json --output classifications.csv
```

QC warnings are generated when:
- Filter rate exceeds 50%
- Genome coverage below 50%
- Single-hit fraction above 80%
- Mean identity below 80%
- Ambiguous fraction above 50%
- Low-confidence fraction above 30%

### Adaptive Thresholds

Detect the species boundary from the ANI matrix distribution instead of using the fixed 96% default:

```bash
# Enable adaptive threshold detection
metadarkmatter score classify --alignment sample.blast.tsv.gz --ani ani.csv \
    --adaptive-thresholds --output classifications.csv
```

This fits a 2-component Gaussian Mixture Model to pairwise ANI values, detecting natural within-species vs. between-species clusters. Falls back to the default threshold if the GMM does not converge or components are not well separated.

Requires: `pip install metadarkmatter[adaptive]` (installs scikit-learn).

See [METHODS.md](METHODS.md) Section 9 for the mathematical details.

### Bayesian Classification (always on)

`score classify` always computes posterior probabilities; the `--bayesian` flag is deprecated and only prints a warning. `taxonomic_call` is the maximum a posteriori (MAP) category from 2D Gaussian likelihoods over (novelty, uncertainty), followed by a discrete Stage 2 refinement (Conserved Region, Ambiguous Within Genus). The hard-threshold result is kept in `legacy_call`.

Posterior columns in the output:
- `p_known_species`, `p_novel_species`, `p_novel_genus`, `p_species_boundary`, `p_ambiguous`, `p_unclassified` (sum to 1.0)
- `posterior_entropy` (0 = confident, log2(6) = 2.585 = uniform across the six categories)
- `confidence_score = (1 - posterior_entropy / log2(6)) * 100`

Parameters can be edited in a YAML file (`score export-config`, then `--config`). Calibration scripts exist for diagnostics, but fitted parameters reduced accuracy in the three-family validation; see [CALIBRATION_RESULTS.md](CALIBRATION_RESULTS.md).

The HTML report shows Bayesian confidence panels (entropy distribution, posterior charts) in the Classification tab.

See [METHODS.md](METHODS.md) Section 10 for the likelihood model and interpretation.

### Family Validation (auto-inferred from `--metadata`, or `--target-family`)

When running BLAST or MMseqs2 against a broad database (e.g., all bacteria rather than a single family), reads may produce hits outside the target family. Family validation detects these off-target reads.

**How it works:**
1. Each BLAST hit is classified as in-family (genome exists in ANI matrix) or external
2. Per-read family metrics are computed (bitscore ratio, identity gap, hit fraction)
3. Reads with `best_in_family_bitscore / best_overall_bitscore < 0.8` are classified as "Off-target"
4. Remaining reads are classified using only in-family hits

**Usage:**

The recommended invocation passes `--metadata` and lets the classifier infer the target family from the most common entry in the `family` column. Pass `--target-family` only when you need to override that choice.

```bash
# Common case: auto-infer from metadata.
metadarkmatter score classify \
    --alignment broad_results.tsv.gz \
    --ani family_ani_matrix.csv \
    --metadata genome_metadata.tsv \
    --output classifications.csv

# Override: pin the target family explicitly.
metadarkmatter score classify \
    --alignment broad_results.tsv.gz \
    --ani family_ani_matrix.csv \
    --target-family "f__Francisellaceae" \
    --output classifications.csv
```

**Output columns** (when family validation is active):
- `family_bitscore_ratio`: Best in-family / best overall bitscore (0.0-1.0)
- `family_identity_gap`: Best in-family identity - best external identity
- `in_family_hit_fraction`: Fraction of hits from in-family genomes
- `external_best_genome`: Best-matching genome outside the family
- `external_best_identity`: Percent identity of best external hit
- `in_family_identity`, `in_family_novelty_index`: Identity and novelty of the best in-family hit

Off-target reads carry `placement_uncertainty = null` (the matched genome is outside the ANI matrix so a uncertainty value would be misleading) and a fixed `confidence_score = 10.0` with `low_confidence = True`. Filter downstream on `taxonomic_call == "Off-target"` rather than on uncertainty.

### Threshold Sensitivity Analysis

Available both as a CLI command (`metadarkmatter score sensitivity --classifications <csv> --output <tsv>`) and as a library API:

```python
from metadarkmatter.core.classification.sensitivity import run_sensitivity_analysis

results = run_sensitivity_analysis(classifications_df, config)
```

Sweeps the novelty/uncertainty thresholds across a configurable range and reports how classification counts change at each point. The sweep re-applies the hard-threshold cascade (the rules behind `legacy_call`) to the saved metrics, keeping the other settings of the base config (`--config`).

### Combined Advanced Options

All new features can be combined:

```bash
metadarkmatter score classify --alignment sample.blast.tsv.gz --ani ani.csv \
    --adaptive-thresholds \
    --qc-output qc_metrics.json \
    --output classifications.csv \
    --summary summary.json \
    --verbose
```

### Performance Modes

Classification uses the VectorizedClassifier by default, which handles most dataset sizes efficiently. For large (more than about 5M alignments) or memory-constrained runs, streaming mode is available:

```bash
# Default mode (handles most datasets efficiently)
metadarkmatter score classify --alignment sample.blast.tsv.gz --ani ani.csv --output out.csv

# Streaming mode (bounded memory)
metadarkmatter score classify --alignment sample.blast.tsv.gz --ani ani.csv --output out.csv \
    --streaming --chunk-size 500000
```

In streaming mode `--id-mapping`/`--genomes` are rejected with an error, `--metadata` still provides the representative mapping and family inference (but no species/genus/family columns), and each read's hits must be contiguous in the input.

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
metadarkmatter ani validate --ani ani_matrix.csv --blast sample.blast.tsv.gz
```

Backend comparison:
- **skani**: Faster, memory-efficient, good for large genome sets (chosen first by `--backend auto`)
- **fastANI**: Original tool, well-validated, widely used in publications (`--fastani-batches N` runs N parallel batches)

`--metadata genome_metadata.tsv --representatives-only` restricts the matrix to species representatives with either backend. The classifier expects values on a 0-100 scale (0-1 matrices are rejected), mirrors pairs present in only one direction, and uses 70% ANI for missing pairs.

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

### Importing External Alignment Results

If alignment was performed outside metadarkmatter (e.g., on a cluster or via a separate pipeline), the results can be imported using ID mapping to translate contig IDs to genome accessions:

```bash
# Generate mapping from genome directory
metadarkmatter util generate-mapping --genomes genomes/ --output id_mapping.tsv

# Validate mapping covers the alignment subject IDs
metadarkmatter util validate-mapping id_mapping.tsv --blast external_results.tsv.gz

# Classify with mapping (not supported in --streaming mode)
metadarkmatter score classify \
    --alignment external_results.tsv.gz \
    --ani ani_matrix.csv \
    --id-mapping id_mapping.tsv \
    --output classifications.csv
```

The `--genomes` flag on `score classify` can also auto-generate the mapping at runtime, but pre-generating with `util generate-mapping` is more efficient when processing multiple samples.

See [WORKFLOW.md](WORKFLOW.md#importing-external-alignment-results) for the complete external import guide.

### Species-Level Tracking

metadarkmatter tracks species-level information throughout the pipeline.

#### Metadata Files

| File | Created By | Contents |
|------|------------|----------|
| `genome_metadata.tsv` | `download genomes list` (written next to `--output`) | accession, species, genus, family, representative, gtdb_taxonomy |
| `{db_prefix}_contig_mapping.tsv` | `blast makedb` | contig_id, genome_accession, original_header |

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

# Output CSV includes species, genus and family columns (not in --streaming mode)
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

With metadata, the report adds a top-species table to the Summary tab and a species read-count bar chart to the Reference tab.

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
    ANI_SPECIES_BOUNDARY_LOW, # 95% ANI threshold
    CATEGORY_NOVEL_SPECIES,   # Classification category strings
    NOVELTY_KNOWN_MAX,        # 4% novelty threshold
    UNCERTAINTY_CONFIDENT_MAX,# 1.5% uncertainty threshold
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
)
```

#### Error Handling

CLI commands use specific exception handlers:
- `PolarsError`: Data processing issues (malformed files)
- `FileNotFoundError`: Missing input files
- `PermissionError`: File access issues
- `MemoryError`: Suggests using `--streaming` mode

### Dependency Injection for Testing

External tool wrappers support dependency injection for testing:

```python
from metadarkmatter.external.base import ExternalTool

# Mock executable resolver for tests
def mock_resolver(name: str) -> str | None:
    return f"/usr/bin/{name}"

ExternalTool.set_executable_resolver(mock_resolver)
# ... run tests ...
ExternalTool.reset_executable_resolver()
```
