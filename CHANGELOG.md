# Changelog

All notable changes to metadarkmatter are recorded here. The format
loosely follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versions use [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] — 2026-05-25

Production-readiness release for internal lab use. Every numbered
phase of the production-readiness plan has either shipped or has a
working scaffold; the calibration validation experiment ran end-to-end
on three families and produced a defensible negative result.

### Added

- `metadarkmatter doctor` — environment, dependency, and external-tool
  diagnostic command. Reports project env vars
  (`METADARKMATTER_SEED`, `METADARKMATTER_GTDB_CACHE_DIR`).
- `metadarkmatter score sensitivity` — sweep classification thresholds
  and emit category counts as TSV or JSON.
- `metadarkmatter score evaluate` — score classifications against a
  ground-truth column; emits per-category precision/recall, confusion
  matrix, and expected calibration error (ECE).
- `metadarkmatter validate fasta` and `validate fastq` — streaming
  format validators (no BioPython dependency).
- `metadarkmatter report vendor-d3` and `report vendor-cdn-sri` —
  cache Plotly/D3 assets and their SRI hashes for offline reports.
- Global `--dry-run` flag and `MDM_DRY_RUN` env var honoured by all
  subcommands that previously supported per-command `--dry-run`.
- `--chunk-size` flag on `score classify` for memory-constrained
  streaming runs.
- `--fastani-batches` flag on `ani compute` for parallel pair-batch
  fastANI runs.
- `--allow-duplicates` flag on `util generate-mapping`; default
  behaviour now raises on duplicate contig IDs.
- `METADARKMATTER_SEED` env var for process-wide random-state control.
- `scripts/build_corpus.py`, `build_synthetic_benchmark.py`,
  `build_metrics_tsv.py`, `calibrate_bayesian.py`,
  `calibrate_entropy.py`, `run_calibration_evaluation.py`,
  `run_pipeline.sh` — full calibration validation toolchain.
- `Containerfile` and `.dockerignore` for the internal-lab image.
- `uv.lock` pinning all 68 transitive Python dependencies.
- GitHub Actions workflow with ruff, mypy, and pytest under `uv sync
  --frozen` on Python 3.11 and 3.12.
- Empirical entropy → confidence isotonic calibration in
  `BayesianConfig.entropy_calibration`; `entropy_to_confidence`
  interpolates between knots when set.
- `BayesianConfig.category_params` for per-category Gaussian overrides
  fitted by `calibrate_bayesian.py`.
- `ScoringConfig.same_genus_ani_threshold`, `same_species_ani_threshold`,
  `merge_band_tolerance` — previously hardcoded literals are now
  configurable.
- New `core/sequence_validation.py`, `core/runtime.py`,
  `core/random.py`, `core/classification/evaluation.py` modules.
- New docs: `STATISTICAL_ASSUMPTIONS.md`, `CALIBRATION.md`,
  `CALIBRATION_RESULTS.md`, `EXTERNAL_TOOLS.md`, `PRODUCTION_READINESS.md`.
- Tests: 1838 passing (up from 1719), 1 deliberate skip for the
  calibration regression guard.

### Changed

- Offline self-contained HTML reports are the new default
  (`--report-mode offline`). Plotly is inlined from the Python package;
  D3 is inlined when vendored via `report vendor-d3`, otherwise CDN.
- BLAST `ToolExecutionError` stderr capture: 500 → 4000 chars.
- GTDB API client now caches successful responses on disk with a 7-day
  TTL and honours `Retry-After` on 429s.
- MMseqs2 multistep workflow now translates subprocess failures to
  `ToolExecutionError` / `ToolTimeoutError` instead of raw
  `CalledProcessError`.
- BAM extraction in `report generate` now degrades to a warning on
  `RuntimeError` (pysam/htslib parse failures) instead of exiting 1.
- `score classify` off-target reads no longer report misleading
  `placement_uncertainty = 0.0`; they now use `null` and a separate
  low-confidence flag.
- `calibrate_bayesian.py` default `--min-samples` raised from 100 to
  5000 after empirical evidence of overfitting at the lower threshold.
- `Containerfile` defaults to dc-megablast-compatible BLAST CLI flags.

### Fixed

- Polars 1.37+ compatibility: replaced
  `col("weighted_bitscore").get(1)` (which now raises on single-hit
  groups) with `col("weighted_bitscore").slice(1, 1).first()`.
- Five long-standing pre-existing test failures resolved (four MMseqs2
  multistep + one BAM corruption-handling).
- `pandas` declared as a runtime dependency (was imported by four
  modules but missing from `pyproject.toml`).
- Off-target classification confidence corrected (see Changed above).
- Duplicate contig IDs in `ContigIdMapping.from_genome_dir` raise
  by default instead of silently overwriting.

### Documented but not shipped

- Calibrated Bayesian YAML defaults. Three-family validation
  (Francisellaceae, Pseudomonas, Lactobacillus) shows fitted
  Gaussians regress accuracy versus the hand-tuned defaults on
  every corpus, in some cases catastrophically. The hand-tuned
  defaults remain the package default; the calibration scripts
  stay as diagnostic tools. See `docs/STATISTICAL_ASSUMPTIONS.md`
  section 13 and `docs/CALIBRATION_RESULTS.md`.

### Production-readiness score

- Composite (internal-lab target): **5.7/10 → 8.4/10**
- Scientific dimension stays at 7/10 (calibration validation
  intentionally produced a negative result, documented as a
  contribution).

## [0.1.0]

Initial release; superseded by the work in 0.2.0.
