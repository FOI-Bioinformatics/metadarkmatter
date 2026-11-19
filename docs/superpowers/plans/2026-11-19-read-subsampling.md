# Read Subsampling Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Cap extracted reads with a reproducible, mate-preserving random subsample in `kraken2 extract`, `run pipeline` and `run cohort`.

**Architecture:** A pure-Python two-pass sampler in `core/subsample.py` (count records, choose indices with a seeded RNG, stream the chosen records). `kraken2 extract` applies it in place after extraction and writes a JSON sidecar. The pipeline forwards the option to the extract step argv, so the existing argv-based fingerprint re-runs extraction and downstream steps when it changes.

**Tech Stack:** Python 3.11, typer, pydantic, pytest.

**Spec:** `docs/superpowers/specs/2026-11-19-read-subsampling-design.md`

## Global Constraints

- Option name: `--max-read-pairs` (int, min 1, default none = no subsampling); `--seed` (int; default `core.random.get_seed()`, i.e. `METADARKMATTER_SEED`, default 42).
- Same input, N and seed give identical output files; input record order preserved.
- R1 and R2 keep the same records (same indices); R1/R2 record-count mismatch is an error.
- When total <= N, files are not rewritten and kept == total.
- Output file names unchanged: `{sample}_taxid{taxid}_R1/_R2.fastq[.gz]`. Sidecar: `{sample}_taxid{taxid}_subsample.json` with keys `total`, `kept`, `max_read_pairs`, `seed`.
- Console line (exact): `Subsampled {kept:,} of {total:,} read pairs (seed {seed})` (use "reads" instead of "read pairs" for single-end).
- No new external dependency. Modest scientific language. Run tests with `uv run --frozen pytest <path> -v`.

## Review Focus

1. FASTQ records whose quality line starts with "@" (valid FASTQ): parsing must be by 4-line records, not by "@" detection. Test in Task 1.
2. Interrupted or failing subsample must not leave a truncated file under the final name: write to a temp file in the same directory and `os.replace`. Test in Task 1 (simulate failure on R2 write; R1 and R2 originals intact).
3. `--skip-if-exists` with existing extracted files: no re-subsampling, sidecar untouched. Test in Task 2.
4. Gzip and plain FASTQ both supported, output compression follows the file extension. Test in Task 1.
5. `run pipeline --extracted-reads-1 ... --max-read-pairs N`: warning that the option only acts on extraction. Test in Task 2.

---

### Task 1: Core sampler

**Files:**
- Create: `src/metadarkmatter/core/subsample.py`
- Test: `tests/unit/test_subsample.py`

**Interfaces:**
- Produces:
  - `@dataclass(frozen=True) class SubsampleResult: total: int; kept: int; seed: int`
  - `subsample_reads(reads_1: Path, reads_2: Path | None, max_pairs: int, seed: int) -> SubsampleResult` — subsamples in place (temp file + `os.replace` per mate, both temp files written before either replace). Reads with `core.paired_reads.open_reads`; writes gzip when the path ends with `.gz`. Index choice: `sorted(random.Random(seed).sample(range(total), max_pairs))`. Raises `ValueError` on R1/R2 count mismatch or a truncated record (line count not a multiple of 4).

- [ ] **Step 1: Write failing tests**: `test_subsample_keeps_n_pairs_and_matches_mates` (100 pairs, N=10: 10 records each, read IDs equal between R1 and R2 after stripping /1 /2 or comment), `test_subsample_is_reproducible` (two copies, same seed -> byte-identical; different seed -> different IDs), `test_subsample_preserves_input_order`, `test_subsample_no_op_when_fewer_reads` (N >= total: files unchanged by mtime/content, kept == total), `test_subsample_single_end`, `test_subsample_gzip_and_plain`, `test_subsample_quality_line_starting_with_at` (Review Focus 1), `test_subsample_count_mismatch_raises`, `test_subsample_failure_leaves_originals` (patch the writer to raise on the R2 temp file; both originals byte-identical, no temp files left) (Review Focus 2).
- [ ] **Step 2: Run** `uv run --frozen pytest tests/unit/test_subsample.py -v` - expect FAIL (module missing).
- [ ] **Step 3: Implement** `core/subsample.py`.
- [ ] **Step 4: Run** - expect PASS.
- [ ] **Step 5: Commit** `feat(subsample): reproducible mate-preserving read subsampling`

### Task 2: CLI, pipeline and docs

**Files:**
- Modify: `src/metadarkmatter/cli/kraken2.py` (`extract_reads`)
- Modify: `src/metadarkmatter/models/pipeline_config.py` (`max_read_pairs: int | None = None` with `ge=1`, `seed: int | None = None`)
- Modify: `src/metadarkmatter/cli/run.py` (`--max-read-pairs`, `--seed` on `run pipeline` and `run cohort`; warning when `--extracted-reads-1` is given with `--max-read-pairs`, using the existing warnings pattern)
- Modify: `src/metadarkmatter/core/pipeline/orchestrator.py` (append `--max-read-pairs N` and, if set, `--seed S` to `ke_argv`)
- Modify: `CLAUDE.md` (Quick Start step 2 note), `docs/REFERENCE.md`, `CHANGELOG.md` (Unreleased/Added)
- Test: `tests/unit/test_cli_kraken2_subsample.py` (new), `tests/unit/test_pipeline_design.py`, `tests/unit/test_cli_pipeline_regressions.py`

**Interfaces:**
- Consumes: `subsample_reads`, `SubsampleResult` (Task 1); `core.random.get_seed`.

- [ ] **Step 1: Write failing tests**
  - `test_extract_max_read_pairs_subsamples_and_writes_sidecar` (patch the KrakenTools run to write small FASTQ outputs; assert N records, sidecar keys/values, console line text).
  - `test_extract_without_max_read_pairs_unchanged` (no sidecar, files untouched).
  - `test_extract_skip_if_exists_does_not_resubsample` (Review Focus 3).
  - `test_extract_seed_default_from_env` (`METADARKMATTER_SEED=7` -> sidecar seed 7).
  - `test_extract_argv_includes_max_read_pairs` (pipeline design: config with `max_read_pairs=250000, seed=42` -> `ke_argv` contains `["--max-read-pairs", "250000", "--seed", "42"]`; absent when None).
  - `test_pipeline_warns_max_read_pairs_with_extracted_reads` (Review Focus 5).
- [ ] **Step 2: Run** the new/changed test files - expect FAIL.
- [ ] **Step 3: Implement** CLI options, config fields, orchestrator argv, warning, docs.
- [ ] **Step 4: Run** `uv run --frozen pytest tests/unit tests/integration -q`, `uv run --frozen ruff check src tests`, `uv run --frozen mypy src` - expect all PASS.
- [ ] **Step 5: Commit** `feat(pipeline): --max-read-pairs subsampling after extraction`
