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

The `mdm` and `metadarkmatter` CLIs are interchangeable. Run
`mdm doctor` first to verify Python deps, external tools, and the
project environment variables (`METADARKMATTER_SEED`,
`METADARKMATTER_GTDB_CACHE_DIR`).

```bash
# 1. Download reference genomes (writes BOTH genomes.tsv and
#    genome_metadata.tsv; the metadata file is consumed by steps 5/6).
mdm download genomes list "f__Francisellaceae" --output genomes.tsv
mdm download genomes fetch --accessions genomes.tsv --output-dir genomes/

# 2. Extract family reads with Kraken2 (requires a built kraken2 DB).
mdm kraken2 extract --kraken-output sample.kraken --reads-1 sample_R1.fastq.gz \
  --reads-2 sample_R2.fastq.gz --taxid 119060 --output extraction/

# 3. Compute ANI matrix (auto-selects fastANI, falls back to skani).
mdm ani compute --genomes genomes/ --output ani_matrix.csv --threads 16

# 4. Sequence alignment (choose ONE option)

# Option A: BLAST (single-end only). Default --task is blastn (sensitive,
# slow); for closely-related references add '--task dc-megablast
# --word-size 11' for ~10x speedup with similar 80-100% ANI sensitivity.
mdm blast makedb --genomes genomes/ --output blastdb/pangenome
mdm blast align --query extraction/reads_R1.fastq.gz --database blastdb/pangenome \
  --output sample.blast.tsv.gz --threads 16

# Option B: MMseqs2 (supports paired-end, best for >100K reads, 5-100x faster)
mdm mmseqs2 makedb --genomes genomes/ --output mmseqs_db/pangenome
mdm mmseqs2 search \
  --query-1 extraction/reads_R1.fastq.gz --query-2 extraction/reads_R2.fastq.gz \
  --database mmseqs_db/pangenome --output sample.mmseqs2.tsv.gz --threads 16

# 5. Classify reads (target family auto-inferred from --metadata).
#    For >5M alignments add '--streaming --chunk-size 500000' to keep
#    memory bounded (vectorized mode estimates ~18 GB at 10M hits;
#    streaming finishes the same workload in seconds).
mdm score classify --alignment sample.blast.tsv.gz --ani ani_matrix.csv \
  --metadata genome_metadata.tsv --output classifications.csv

# 6. Generate HTML report (offline self-contained by default; pass
#    '--report-mode cdn' for the smaller-but-network-dependent variant).
mdm report generate --classifications classifications.csv \
  --metadata genome_metadata.tsv --output report.html
```

End-to-end real-data verification: `docs/CASE_STUDIES/SRR25038281_chesapeake_cb52.md`.

## Importing External Alignment Results

External BLAST/MMseqs2 results can be classified using ID mapping to translate contig IDs to genome accessions:

```bash
# Generate mapping from genome directory
metadarkmatter util generate-mapping --genomes genomes/ --output id_mapping.tsv

# Validate mapping against BLAST file
metadarkmatter util validate-mapping id_mapping.tsv --blast external_results.tsv.gz

# Classify with mapping
metadarkmatter score classify --alignment external_results.tsv.gz --ani ani_matrix.csv \
  --id-mapping id_mapping.tsv --output classifications.csv
```

**Key Files:**
- `cli/mapping.py` - `generate-mapping`, `validate-mapping`, `select-representatives` commands
- `core/id_mapping.py` - `ContigIdMapping` class (from_genome_dir, from_tsv, to_tsv, transform_column)
- `core/genome_utils.py` - `extract_accession_from_filename()`

**Requirements:**
- Alignment must be BLAST tabular format (outfmt 6), 12 or 13 columns, tab-separated
- ID mapping is not supported with `--streaming` mode
- Genome accessions in mapping must match ANI matrix row/column labels

## Representative Genome Architecture

For large families (>2000 genomes), computing a full ANI matrix is impractical (O(N^2)). The representative genome architecture decouples alignment from classification:

- **Align** against ALL genomes (full sensitivity)
- **Compute ANI** for species representatives only (one per species from GTDB)
- **Classify** by mapping each BLAST hit genome to its representative for ANI lookup

```bash
# Download all genomes with auto-identified representatives
metadarkmatter download genomes list "f__Pseudomonadaceae" --all-genomes --output genomes.tsv

# Compute ANI for representatives only (~85 instead of ~2000)
metadarkmatter ani compute --genomes genomes/ --output ani_matrix.csv \
    --metadata genome_metadata.tsv --representatives-only --threads 16

# Build alignment database from ALL genomes
metadarkmatter blast makedb --genomes genomes/ --output blastdb/pangenome

# Classify (metadata provides representative mapping automatically)
metadarkmatter score classify --alignment sample.blast.tsv.gz --ani ani_matrix.csv \
    --metadata genome_metadata.tsv --output classifications.csv
```

**Key Files:**
- `models/genomes.py` - `GenomeAccession.is_representative`, `AccessionList.representative_map`
- `core/metadata.py` - `GenomeMetadata.get_representative()`, `build_representative_mapping()`, `has_representatives`, `representative_count`
- `cli/mapping.py` - `select-representatives` command
- `cli/ani.py` - `--metadata` + `--representatives-only` flags on compute; `--metadata` on validate
- `cli/score.py` - Wires representative mapping from metadata to VectorizedClassifier
- `core/classification/classifiers/vectorized.py` - `representative_mapping` parameter, `representative_genome` column for ANI lookups

**Metadata format** (`genome_metadata.tsv`):
```
accession	species	genus	family	representative	gtdb_taxonomy
GCF_000195955.2	F. tularensis	Francisella	Francisellaceae	GCF_000195955.2	d__...
GCF_001234567.1	F. tularensis	Francisella	Francisellaceae	GCF_000195955.2	d__...
```

**Backwards compatibility:** If `representative` column is absent, every genome maps to itself (current behavior preserved).

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
- **[docs/STATISTICAL_ASSUMPTIONS.md](docs/STATISTICAL_ASSUMPTIONS.md)** - Load-bearing modelling assumptions, sources, and empirically observed failure modes
- **[docs/STATISTICS_AUDIT.md](docs/STATISTICS_AUDIT.md)** - 2026-05 audit of classifier math, ANI symmetry fix, calibration guardrails, and the per-read calibration experiment orchestrator (`scripts/run_per_read_calibration.sh`)
- **[docs/CALIBRATION.md](docs/CALIBRATION.md)** - Bayesian calibration pipeline (build benchmark -> fit Gaussians -> fit entropy -> evaluate)
- **[docs/CALIBRATION_RESULTS.md](docs/CALIBRATION_RESULTS.md)** - Three-family validation. **Hand-tuned defaults beat every fitted YAML; do not ship calibrated configs as the default.** Calibration scripts remain as diagnostic tools.
- **[docs/EXTERNAL_TOOLS.md](docs/EXTERNAL_TOOLS.md)** - Tested external tool versions and reproducibility checklist
- **[docs/PRODUCTION_READINESS.md](docs/PRODUCTION_READINESS.md)** - Production-readiness audit (8.4/10 internal-lab composite)
- **[docs/CASE_STUDIES/](docs/CASE_STUDIES/)** - End-to-end runs against real published samples (start with `SRR25038281_chesapeake_cb52.md`)

## v0.2.0 production-readiness additions

The audit pass that shipped in v0.2.0 (see `CHANGELOG.md`) introduced
several cross-cutting features that aren't obvious from the per-command
documentation:

- **`mdm doctor`** - reports Python version, dependency versions, external
  tool paths/versions, and project env vars. Run before any long
  pipeline invocation.
- **`mdm validate fasta` / `mdm validate fastq`** - streaming format
  validators (no BioPython dependency); useful pre-flight checks.
- **Global `--dry-run`** at the top-level CLI and the `MDM_DRY_RUN=1`
  env var. Every subcommand that supports per-command `--dry-run`
  also honours the global flag.
- **`METADARKMATTER_SEED`** env var controls every randomness-using
  site (plot subsampling, adaptive GMM fits). Default 42.
- **`METADARKMATTER_GTDB_CACHE_DIR`** env var enables an on-disk
  cache for GTDB API responses with a 7-day TTL. Empty value disables.
- **Offline self-contained reports** are the default; SRI integrity
  attributes on CDN-mode tags via `mdm report vendor-cdn-sri`.
- **Configurable thresholds** that were previously hardcoded:
  `ScoringConfig.same_genus_ani_threshold`,
  `same_species_ani_threshold`, `merge_band_tolerance` (see
  `models/config.py`).
- **Container build**: `Containerfile` at the repo root pins all
  external tools via bioconda; `uv.lock` pins every transitive Python
  dependency for `uv sync --frozen`.
- **Determinism guard**: `tests/integration/test_pipeline_determinism.py`
  runs `score classify` twice and asserts SHA-256 equality of the
  output CSVs.

## Calibration: diagnostic only, not for shipping defaults

`scripts/calibrate_bayesian.py` and `scripts/calibrate_entropy.py`
exist and run end-to-end; the three-family validation (`docs/CALIBRATION_RESULTS.md`)
empirically shows that fitted Gaussian parameters **regress
classification accuracy** on every test corpus, sometimes
catastrophically (Lactobacillus calibration drops Francisella
accuracy from 22% to 3.6%). The mechanism is documented in
`STATISTICAL_ASSUMPTIONS.md` section 13: per-genome ANI labels do
not map cleanly to per-read (novelty, uncertainty) clusters.

Treat the scripts as diagnostic tools (per-category mean, sigma, and
empirical r(N, U) are informative) and do NOT commit fitted
calibrated YAMLs to `configs/` as new defaults without a fresh
validation against a held-out family.

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
- BioPython (runtime dependency in pyproject.toml)
- D3.js for the rendered tree. Reports built with the default
  `--report-mode offline` inline a vendored copy of D3 if present
  (`mdm report vendor-d3` downloads and caches it next to the package
  assets); otherwise they emit a CDN `<script src>` tag. CDN tags
  carry an `integrity=sha384-...` attribute when
  `mdm report vendor-cdn-sri` has been run.

## ANI/AAI Heatmap Clustering

scipy is a dependency; hierarchical clustering is enabled by default for all heatmaps.

**Key Files:**
- `visualization/report/components/clustering.py` - `perform_hierarchical_clustering()` (shared by ANI and AAI)
- `visualization/report/components/heatmap_builder.py` - `HeatmapConfig` dataclass, `_build_similarity_heatmap()` shared renderer, thin wrappers `build_ani_heatmap()` / `build_aai_heatmap()`

Graceful degradation: logs a warning if scipy is unavailable and returns the unclustered matrix.

## Key Architecture

```
src/metadarkmatter/
├── cli/              # Typer CLI commands (entry points; doctor, score sensitivity/evaluate, validate fasta/fastq added in v0.2.0)
├── core/             # Core algorithms (ANI classification, parsers, metadata)
│   ├── constants.py          # Nucleotide thresholds (default mode)
│   ├── protein_constants.py  # Protein thresholds (--alignment-mode protein)
│   ├── random.py             # METADARKMATTER_SEED resolution (v0.2.0)
│   ├── runtime.py            # Process-wide --dry-run via MDM_DRY_RUN (v0.2.0)
│   ├── sequence_validation.py # FASTA/FASTQ streaming validators (v0.2.0)
│   ├── classification/       # Classification pipeline
│   │   ├── bayesian.py       # Bayesian-primary classifier (6-cat 2D Gaussians + Stage 2)
│   │   ├── thresholds.py     # Legacy threshold cascade + apply_legacy_thresholds()
│   │   ├── qc.py             # Pre/post-classification QC metrics
│   │   ├── sensitivity.py    # Threshold sensitivity analysis (also exposed via 'mdm score sensitivity')
│   │   ├── evaluation.py     # Accuracy / ECE / confusion matrix (v0.2.0; exposed via 'mdm score evaluate')
│   │   ├── adaptive.py       # GMM-based adaptive threshold detection
│   │   ├── ani_matrix.py     # ANIMatrix class
│   │   └── classifiers/      # Classifier implementations
│   │       ├── base.py       # ANIWeightedClassifier (programmatic API, Bayesian-primary)
│   │       └── vectorized.py # VectorizedClassifier (sole CLI classifier, Bayesian-primary)
│   ├── novel_diversity/      # Novel taxa clustering and models
│   │   ├── clustering.py     # Cluster novel reads by genome proximity
│   │   └── models.py         # NovelClusterResult, NovelReadInfo dataclasses
│   └── phylogeny/            # Phylogenetic tree building and novel cluster placement
│       ├── tree_builder.py   # ANI-to-Newick conversion, user tree loading
│       └── placement.py      # Novel cluster extraction and tree placement
├── external/         # External tool wrappers (BLAST, MMseqs2, Kraken2, fastANI with --fastani-batches, etc.)
├── models/           # Pydantic data models (ScoringConfig with YAML from_yaml/to_yaml; BayesianConfig.category_params + entropy_calibration in v0.2.0)
├── clients/          # API clients (GTDB with retry logic, on-disk TTL cache via METADARKMATTER_GTDB_CACHE_DIR)
└── visualization/    # Plotly charts and HTML report generation
    └── report/
        ├── generator.py      # Single-sample reports; ReportConfig.report_mode='offline' default in v0.2.0
        ├── multi_generator.py # Multi-sample comparative reports
        ├── assets/           # Bundled JS (Plotly inlined; D3 vendored on demand) + SRI cache (v0.2.0)
        └── components/       # Report building blocks
            ├── clustering.py          # perform_hierarchical_clustering()
            ├── heatmap_builder.py     # HeatmapConfig, unified ANI/AAI heatmaps
            └── extended_matrix_builder.py # Novel-to-reference extended matrices
```

Companion artifacts at the repo root:

- `Containerfile` + `.dockerignore` - reproducible internal-lab image.
- `uv.lock` - pinned Python deps; `uv sync --frozen` recreates the env.
- `scripts/` - calibration pipeline (build_corpus, build_metrics_tsv,
  calibrate_bayesian, calibrate_entropy, run_calibration_evaluation,
  build_synthetic_benchmark) + the canonical end-to-end
  `run_pipeline.sh`.
- `.github/workflows/test.yml` - ruff + mypy + pytest under
  `uv sync --frozen` on Python 3.11 and 3.12.

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
- `models/config.py` - `ScoringConfig` (sole config class, includes `get_effective_thresholds()`)
- `core/classification/classifiers/vectorized.py` - VectorizedClassifier (sole CLI classifier)

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

### Bayesian-Primary Classification (always on)
- Primary classifier: 6-category 2D Gaussian posteriors in (novelty, uncertainty) space
- Categories: Known Species, Novel Species, Novel Genus, Species Boundary, Ambiguous, Unclassified
- Two-stage: Stage 1 (Bayesian MAP) + Stage 2 (discrete refinement for Conserved Region, Ambiguous Within Genus)
- Prior modulation: identity_gap_boost (2.0x) and single_hit_boost (1.5x) for Ambiguous prior
- Shannon entropy as confidence metric (0 = confident, 2.585 = uniform over 6 categories)
- `confidence_score = (1 - entropy / log2(6)) * 100` (0-100 scale)
- Output columns: taxonomic_call (Bayesian MAP), legacy_call (hard thresholds), 6 posterior columns, posterior_entropy, confidence_score
- `--config` loads YAML config; `score export-config` generates editable YAML
- `--include-legacy-scores` adds old sub-scores (alignment_quality, etc.)
- `--bayesian` flag deprecated (always-on); numeric flags deprecated when `--config` used
- Key files: `core/classification/bayesian.py`, `models/config.py` (YAML methods), `cli/score.py` (export-config)

### Family Validation (`--target-family`, or inferred from `--metadata`)
- Detects off-target reads from broad-database alignments
- Pass `--target-family` explicitly, or let it auto-infer the most
  common family from the `family` column in `genome_metadata.tsv`
- Partitions hits by ANI matrix membership (in-family vs external)
- Reads with best_in_family / best_overall bitscore < 0.8 are Off-target
- Off-target rows now report `placement_uncertainty=null` and a low
  fixed `confidence_score=10.0` (previously 0.0, which was misleading)
- Output columns: family_bitscore_ratio, family_identity_gap, in_family_hit_fraction
- Key file: `core/classification/classifiers/vectorized.py`

### Sensitivity Analysis (`score sensitivity`)
- Sweeps novelty/uncertainty thresholds across a configurable range
- Re-classifies reads at each threshold point and emits a long-form TSV
  (one row per threshold-point / category) or JSON
- Library API is also available via `run_sensitivity_analysis()` in
  `core/classification/sensitivity.py`
- Key file: `cli/score.py` (the `sensitivity` subcommand)

### Evaluation (`score evaluate`)
- Scores predictions against a ground-truth column; reports top-1
  accuracy, per-category precision/recall, confusion matrix, and
  expected calibration error (ECE)
- Used to close the calibration loop (build benchmark -> calibrate ->
  evaluate)
- Key files: `cli/score.py` (the `evaluate` subcommand),
  `core/classification/evaluation.py`

### Enhanced Scoring (always on)
- Inferred uncertainty for single-hit reads
- Alignment quality, identity confidence, placement confidence
- Discovery score for novel taxa prioritization
- All computed automatically (no flags needed)
