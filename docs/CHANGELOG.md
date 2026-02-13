# Changelog

All notable changes to metadarkmatter are documented in this file.

## [0.2.0] - 2026-02-10

### Breaking Changes

**Default Threshold Updates**
- Species boundary tightened: `novelty_known_max` changed from 5.0% to 4.0% (96% ANI, Jain et al. 2018)
- Uncertainty threshold tightened: `uncertainty_known_max` changed from 2.0% to 1.5%
- Coverage weighting enabled by default: `coverage_weight_mode` changed from "none" to "linear"
- Minimum alignment fraction: `min_alignment_fraction` changed from 0.0 to 0.3

**Removed Feature Flags**
- Removed `--enhanced-scoring` CLI flag (now always on)
- Removed `--infer-single-hit-uncertainty` CLI flag (now always on)
- Removed `--use-inferred-for-single-hits` CLI flag (now always on)
- Removed presets: `default`, `balanced-conservative`, `coverage-linear`, `coverage-gentle`
- Removed `enhanced_scoring`, `infer_single_hit_uncertainty`, `use_inferred_for_single_hits` config fields

### Added

**Quality Control Metrics (Phase 1)**
- New `core/classification/qc.py` module with pre- and post-classification QC metrics
- Pre-classification checks: filter rate, genome coverage, single-hit fraction, mean identity
- Post-classification checks: ambiguous fraction, low-confidence fraction
- Automated warnings for problematic inputs (high filter rate, low genome coverage, etc.)
- New `--qc-output` CLI flag to save QC metrics as JSON

**Threshold Sensitivity Analysis (Phase 2)**
- New `core/classification/sensitivity.py` module for robustness assessment
- Sweeps novelty/uncertainty thresholds across configurable ranges
- Re-classifies reads at each threshold without re-computing alignment metrics
- New `core/classification/thresholds.py` with extracted threshold application logic
- New `metadarkmatter score sensitivity` subcommand

**Adaptive Threshold Detection (Phase 3)**
- New `core/classification/adaptive.py` module using Gaussian Mixture Models
- Detects natural species boundary from ANI matrix distribution
- Falls back to default threshold when GMM does not converge or separation is poor
- Species boundary clamped to 80-99% ANI range
- New `--adaptive-thresholds` CLI flag
- Optional dependency: `scikit-learn>=1.3.0` (install with `pip install metadarkmatter[adaptive]`)

**Bayesian Confidence Framework (Phase 4)**
- New `core/classification/bayesian.py` with vectorized numpy implementation
- Posterior probabilities P(category | novelty, uncertainty) using 2D Gaussian likelihoods
- Shannon entropy as single-scalar confidence metric (0 = confident, 2.0 = uniform)
- MAP (Maximum A Posteriori) classification column
- Ambiguous boosting for high-uncertainty and low-gap reads
- New `--bayesian` CLI flag adds 6 columns to output: p_known_species, p_novel_species, p_novel_genus, p_ambiguous, bayesian_category, posterior_entropy
- Bayesian Confidence tab in HTML report with 3 interactive visualizations

### Changed

- Enhanced scoring metrics (inferred uncertainty, alignment quality, identity/placement confidence, discovery score) are now always computed
- Coverage weighting (linear mode) is now the default for all classifications
- Minimum alignment fraction defaults to 30%, filtering low-coverage BLAST hits
- Constants updated in `core/constants.py` to reflect new threshold values

---

## [Unreleased] - 2026-01-21

### Added

**Novel Diversity Analysis**
- New `core/novel_diversity/` module for analyzing reads classified as Novel Species or Novel Genus
- Clusters novel reads based on nearest reference genome and divergence level
- Assigns confidence ratings (High/Medium/Low) based on read count, uncertainty, and discovery score
- Generates suggested names for novel clusters (e.g., "Francisella sp. nov. NSP_001")
- Provides phylogenetic context placement within family/genus hierarchy

**Novel Diversity Tab in Reports**
- New "Novel Diversity" tab appears when novel reads are detected
- Summary metrics: total clusters, novel species/genus counts, high-confidence clusters
- Cluster quality scatter plot showing novelty vs discovery score
- Phylogenetic context sunburst chart for taxonomic hierarchy visualization
- Phylogenetic Context Heatmap showing novel clusters alongside reference genomes
- Unified "Novel Clusters" table with integrated phylogenetic placement information

**Extended Matrix Builder**
- New `extended_matrix_builder.py` for building similarity matrices that include novel clusters
- Estimates ANI/AAI between novel clusters and reference genomes
- Supports hierarchical clustering for visual grouping

**Enhanced Scoring System**
- Discovery scores for prioritizing novel taxa validation
- Identity confidence and placement confidence metrics
- Inferred uncertainty for single-hit reads
- Alignment quality scoring

**Scientific Methods Documentation**
- New `docs/METHODS.md` with comprehensive documentation of all calculations
- Publication-ready methods section covering classification algorithm, enhanced scoring, and novel diversity clustering
- Complete formulas for Novelty Index, Placement Uncertainty, confidence scores, and discovery scores
- Documentation of coverage weighting modes and protein mode thresholds
- Updated all existing documentation to reference METHODS.md as authoritative source

### Changed

**Report Tab Reorganization**
- Reordered tabs for logical interpretation flow: Overview, Distributions, Species, Genomes, Novel Diversity, Reference ANI, Reference AAI, Discovery Scores, Recruitment, Data
- Renamed "ANI Matrix" to "Reference ANI" for clarity
- Renamed "AAI Matrix" to "Reference AAI" for clarity
- Renamed "Enhanced Scoring" to "Discovery Scores" for clarity

**Heatmap Improvements**
- Phylogenetic Context Heatmap now uses dynamic width based on matrix size (1000/1200/1400px for <30/30-50/>50 entities)
- Reference AAI heatmap scale changed from 40-75% to 40-100% for consistency with Novel Diversity AAI heatmap
- Improved colorscale consistency across all AAI visualizations

**Novel Diversity Table Consolidation**
- Merged "Cluster Details" and "Phylogenetic Context" tables into single "Novel Clusters" table
- New "Phylogenetic Placement" column combines suggested name, context narrative, and nearest reference
- Reduced redundancy while preserving all information

### Fixed

- AAI heatmap scale inconsistency between Reference AAI and Novel Diversity tabs

---

## [Previous] - 2026-01-20

### Added

**BLAST FASTQ Support**
- BLAST now accepts FASTQ files directly (`.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`)
- Automatic conversion to FASTA using `seqtk` (primary) or Python fallback
- Transparent temporary file management with automatic cleanup
- Users no longer need to manually convert FASTQ to FASTA

**MMseqs2 Integration**
- MMseqs2 support as BLAST alternative for very large datasets (>100K reads)
- Multi-step workflow with query database caching for performance
- BLAST-compatible output format for seamless pipeline integration
- CLI commands: `metadarkmatter mmseqs2 makedb` and `metadarkmatter mmseqs2 search`

**Project Organization**
- External software dependencies clearly documented in all main files
- Organized test data into `test_data/` directory
- Added `.gitignore` for Python artifacts and test outputs

### Changed

**Documentation Updates**
- Streamlined MMseqs2 documentation with realistic performance expectations
- MMseqs2 is 20x slower than BLAST for small datasets (<10K reads)
- MMseqs2 provides 5-15x speedup only for datasets >100K reads
- Updated all documentation to reflect BLAST as primary method
- Made `seqtk` optional (only needed for assembly workflows)
- Reorganized documentation with tutorial as primary entry point

**Performance Guidance**
- Added realistic performance tables based on actual testing
- Clear dataset size thresholds for tool selection
- BLAST recommended for <100K reads, MMseqs2 for larger datasets

### Fixed

**MMseqs2 Issues**
- Fixed duplicate `-q` parameter warning in CLI
- Implemented query database caching to prevent 30+ min overhead per run
- Added missing `--search-type` parameter to search command
- Fixed `get_executable_path()` method name error

**Documentation Accuracy**
- Corrected misleading "100-1000x speedup" claims for all datasets
- Removed overstated MMseqs2 performance claims
- Added warnings about dataset size requirements

### Developer Notes

**File Reorganization**
- Moved documentation summaries from root to `docs/`
- Consolidated test files into `test_data/demo_outputs/` and `test_data/e2e_tests/`
- Created CHANGELOG.md (this file) from CLEANUP_SUMMARY.md, DOCUMENTATION_UPDATE_SUMMARY.md

**External Dependencies**
- Required: kraken2, krakentools, blast, skani
- Optional: mmseqs2 (>100K reads), seqtk (assembly workflows)

## Implementation Details

For detailed implementation notes, see:
- `docs/MMSEQS2_IMPLEMENTATION_SUMMARY.md` - Complete MMseqs2 integration details
- `docs/DOCUMENTATION_UPDATE_SUMMARY.md` - Documentation update details
