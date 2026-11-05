# Changelog

All notable changes to metadarkmatter are recorded here. The format
loosely follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versions use [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Production-readiness and UX hardening pass.

### Added

- `metadarkmatter run pipeline` — single end-to-end orchestration command
  that chains download -> kraken2 extract -> align -> ANI -> classify ->
  report. Reuses the existing subcommands as steps, with a structured output
  directory, per-step checkpointing (skips steps whose outputs already exist),
  `--from` resume, `--dry-run` planning, and a `run_manifest.json` recording
  provenance and per-step status/timing. Driven by CLI flags or a `--config`
  YAML (`PipelineConfig`).
- Centralized CLI error handling: `MetadarkmatterError` subclasses now print a
  clean `message` + `suggestion`; common library/runtime errors are translated
  to friendly messages; full tracebacks are shown only under the new global
  `--debug` flag (`MDM_DEBUG=1`).
- Shared CLI option types (`cli/options.py`) for consistent
  `--quiet/--verbose/--dry-run/--threads`.
- CI: macOS matrix leg, experimental Python 3.13 leg, a `pip-audit`
  dependency-vulnerability job, and a tool-dependent job that runs the
  determinism guard plus integration/e2e tests against bioconda tools.
- `.pre-commit-config.yaml` (ruff + mypy mirror of CI) and an internal,
  tag-triggered `release.yml` (CHANGELOG/version check, container build+smoke,
  GitHub release; no public publish).
- Unit tests for the previously untested `core/ani_matrix_builder.py`, the
  pipeline orchestrator, and the centralized error handler.

- Opt-in clade-novelty rule (`score classify --clade-novel-species`,
  `ScoringConfig.clade_novel_species`, YAML `clade_novelty`): reads called
  Ambiguous or Species Boundary because their competitors are close
  relatives within one genus, but that are at least 5% divergent from both
  the best and the second-best reference, are called Novel Species. Off by
  default. In a leave-one-out benchmark on Francisellaceae,
  Lactobacillaceae and Enterobacteriaceae it raised read-level
  novel-species recall from 0.47-0.52 to 0.59-0.70 with novel-versus-known
  precision unchanged (>= 0.989); see docs/CLADE_NOVELTY_VALIDATION.md.
- `scripts/benchmark_leave_one_out.py`: leave-one-out novelty benchmark that
  compares ScoringConfigs on any genome set.

### Changed

- mypy is now **blocking** in CI. `strict` stays on with a documented set of
  deferred error codes (mostly Polars-expression operator typing and untyped
  third-party libraries); the enforced core catches signature/return/overload
  bugs. See `docs/PRODUCTION_READINESS.md` for the deferred list.
- Coverage gate is now enforced in CI (`--cov-fail-under=78`).
- Removed redundant inline `try/except` blocks in `score classify` / `aai
  validate` so typed errors flow to the centralized handler.
- Behaviour changes from the 2026-11 pipeline audit (PR #13):
  - `blast align`, `mmseqs2 search` and `blastx align` choose compression from
    the output extension (`.gz` is gzipped, `.tsv` stays plain); `--compress`
    forces gzip. Previously output was always gzipped.
  - `mmseqs2 search --query-1/--query-2` gives mates distinct `/1` and `/2`
    IDs, so the two mates are classified as separate reads instead of being
    pooled under one shared ID.
  - `blast align --word-size` defaults per `--task` (blastn 7, dc-megablast
    11, megablast 28).
  - External tools no longer have a default one-hour timeout, and a timeout is
    now reported as a tool error.
  - An ANI matrix on a 0-1 scale is rejected; a pair given in only one
    direction is mirrored. ANI values are stored as float64.
  - `score classify --streaming` refuses `--id-mapping`/`--genomes` (the
    mapping used to be ignored without warning), and stops with an error when a
    read's hits are not contiguous in the file.
  - `mdm run pipeline`: values from a `--config` YAML are no longer overridden
    by CLI defaults; `--from STEP` re-runs that step and all later ones; new
    `--extracted-reads-2`; the run manifest records start and finish times.

### Fixed

- Fixes from the 2026-11 pipeline audit (PR #13), checked against a real
  Francisellaceae BLAST/MMseqs2 run:
  - MMseqs2/DIAMOND reverse-strand hits (`qstart > qend`) had negative
    coverage, so the correct hit was down-weighted or filtered out. Reads from
    reference genomes were called Novel Species (689 of about 3,000 in the test
    run) and single-end MMseqs2 input lost reads. Coordinates are now
    normalized.
  - With tied hits, the "second" competing genome (default
    `uncertainty_mode`) came from an unstable sort, so
    `placement_uncertainty` depended on batch size and streaming disagreed
    with in-memory mode.
  - Streaming: crashes on integer-looking bitscores, numeric read IDs and `#`
    headers; Parquet schema changing between batches; misaligned CSV columns
    for batches containing only Off-target reads.
  - `--adaptive-thresholds` dropped the user's configuration, and uncomputed
    ANI pairs pulled the detected boundary down to 80%. `score sensitivity`
    crashed on `score classify` output. The QC low-confidence fraction
    compared the 0-100 confidence score against 0.5.
  - `report generate` crashed on Off-target rows and on empty CSVs. Reports
    were not byte-identical across processes. Metadata names were inserted
    into the HTML unescaped. Novel clustering dropped reads before adjacent
    novelty bands were merged.
  - `mdm run pipeline` could not finish from raw reads: extracted-read
    filenames did not match what `kraken2 extract` writes, kraken2 output was
    passed as a file instead of a directory, `--metadata` was missing on first
    runs, and MMseqs2 single-end input failed. `scripts/run_pipeline.sh`
    failed at the extract step because no `.kreport` was written.
  - `ani compute --representatives-only` was ignored by the fastANI backend.
  - `--preset` silently dropped explicitly passed flags. The global
    `--dry-run` was ignored by `score classify` and `blast makedb/align`.
  - `blast align`: gzipped FASTA queries were parsed as FASTQ, dropping every
    other read; multi-volume and dotted-prefix databases were rejected.
  - MMseqs2 translated search (`--search-type 2`) built an amino-acid query
    database and found nothing.
  - ID mapping now handles NCBI-style `ref|ACC|` seqids and numeric contig
    IDs. Genome concatenation handles files without a trailing newline. GTDB
    placeholder taxa such as `f__CAG-138` are accepted.
  - `anyio` and `certifi` bumped in `uv.lock` for published advisories.
- Fixes found during the 2026-11 documentation audit:
  - `score classify --aai` loaded the AAI matrix but never used it. The CLI
    classifier now applies the same AAI genus rule as `ANIWeightedClassifier`
    (novel-genus range, low U: AAI > 65% -> Novel Species, 58-65% ->
    Ambiguous Within Genus, < 58% -> Novel Genus) to both `taxonomic_call` and
    `legacy_call`. Output with `--aai` can therefore differ from before for
    novel-genus candidates.
  - `ambiguity_scope` counted hit rows instead of distinct genomes, so a
    competitor with several contig hits made a read `across_genera`, which can
    turn Ambiguous into Conserved Region.
  - `score classify --summary` wrote 0 for `diversity_known/novel/uncertain`,
    and its species counts were empty because the summary was built from the
    table before the metadata join.
  - `BlastHit.genome_name` returned the part after `|` for non-accession IDs
    (e.g. MAG bins `bin_01|contig_5`), unlike the CLI classifier and the
    `{genome}|{contig}` header standard.
  - `ANIWeightedClassifier` reported U = 100 (Ambiguous) when all ambiguous
    hits were to the best genome; it now reports U = 0, as the CLI classifier
    does. Its single-hit prior boost used the nucleotide novelty threshold in
    protein mode.
  - Help text: `--bayesian` is described as deprecated, `--single-hit-threshold`
    no longer refers to a non-existent flag, and `score batch --workers` is
    marked as unused. The `pyproject.toml` repository URL now points at
    FOI-Bioinformatics/metadarkmatter.
- Fixes from the second 2026-11 audit round:
  - Protein mode (`--alignment-mode protein`) classified 0 reads with the
    default filters: `min_alignment_length` (100) was compared with BLASTX
    alignment length in amino acids. It now uses the query span in
    nucleotides. A run that classifies 0 reads writes a header-only file and
    warns (stale outputs were left in place).
  - Subject IDs without an `{accession}|` prefix made every read Known
    Species with no warning; this is now an error pointing to `--id-mapping`,
    as is 0% overlap between hit genomes and the ANI matrix. Hit genomes
    missing from the ANI matrix are reported with a count.
  - Blank metadata representatives, duplicate accessions, NaN/NA ANI cells
    and single-genome matrices are handled; non-finite N/U give Unclassified.
  - `validate blast` and `validate classifications` rejected the pipeline's
    own outputs. `score batch` exits 1 when a file fails; `score sensitivity`
    labels its output `legacy_call`; `score extract-novel` option fixes;
    `score evaluate` tolerates unlabelled reads.
  - `run pipeline` re-ran nothing when inputs changed and treated
    interrupted steps as done; steps now carry a fingerprint and a
    completion marker. BLAST with R2 reads warns.
  - Report: the Diversity Classification sunburst was blank whenever there
    were no off-target reads; the Methods panel now describes the Bayesian
    call; ANI pair counts and missing-ANI display corrected; the phylogeny
    tree and the Novel Diversity tab now share one set of clusters (they
    used the same IDs for different clusters); cluster uncertainty is the
    measured value (was a fabricated estimate for single-hit reads);
    multi-sample report shows all categories offline; mobile layout.
  - `ani compute` warns when many genome pairs are undetected (skani left
    49% of a 10-genome set undetected) and adds `--skani-mode`/`--skani-c`;
    fastANI asymmetric pairs are averaged; recruitment identity uses the NM
    tag; Bowtie2, KrakenTools and skani argument fixes; `map reads` reuses
    an existing index. New flags: `download genomes fetch --require-all`,
    `doctor --strict`, `util validate-mapping --min-coverage`.

- Streaming classification (`score classify --streaming` / `stream_to_file`)
  silently bypassed the alignment-quality filters (`min_alignment_length`,
  `min_alignment_fraction`, `max_evalue`, `min_percent_identity`,
  `min_bitscore`, `min_read_length`, `min_query_coverage`) that the in-memory
  path applied, so it classified more reads than the default path for the same
  input. The filters are now shared by both paths; a regression test asserts
  identical read sets.
- Streaming parquet output is now written incrementally (one row group per
  partition) instead of re-reading and rewriting the whole file per partition,
  cutting streaming time from O(partitions^2) to O(rows).
- `metadarkmatter map reads` against a genome directory crashed with a
  `TypeError` (passed a non-existent `prefix_headers` argument to
  `concatenate_genomes`). Surfaced by enabling blocking mypy.
- Brought the lint baseline to green under ruff 0.15.13 (the lockfile's pinned
  version), which had regressed to 469 reported issues.

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
