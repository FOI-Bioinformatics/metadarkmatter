# Read subsampling after Kraken2 extraction

Date: 2026-11-19
Status: draft for approval

## Intent

Samples dominated by the target family (e.g. B_iv_01_M_S1: 10.3M Brucellaceae
read pairs) produce alignments too large for a laptop, because every read hits
most of the closely related references. The user wants to cap the number of
reads inside metadarkmatter instead of running seqtk by hand.

## Goals

- `mdm kraken2 extract --max-read-pairs N` keeps a uniform random subsample of
  at most N read pairs (or N reads for single-end input) of the extracted reads.
- The same option on `mdm run pipeline` and `mdm run cohort`, so the
  subsample replaces the extracted reads and every downstream step uses it.
- Reproducible: the same input, N and seed give the same subsample.
- Paired-end mates stay matched.
- No new external dependency (no seqtk).

## Non-goals

- Subsampling options on the aligners (approach B, rejected).
- A separate `util subsample` command (approach C, rejected).
- Showing the subsample in the HTML report (the report does not read the run
  manifest today; adding that plumbing is out of scope).

## Approaches

- A (chosen): option on extraction and pipeline; the subsample replaces the
  extracted files.
- B: option on `mmseqs2 search` / `blast align`; rejected because the cap would
  be split across two commands and the pipeline would still write all
  extracted reads first.
- C: separate subsample command plus pipeline flag; rejected as an extra step
  for manual users.

## Design

- New module `core/subsample.py`:
  `subsample_reads(reads_1, reads_2, output_1, output_2, max_pairs, seed) -> SubsampleResult(total, kept, seed)`.
  Two passes: count records in R1 (and check R2 has the same count), choose
  `max_pairs` record indices with `random.Random(seed).sample`, then stream
  both files and write the chosen records. Input order is preserved. When
  `total <= max_pairs`, files are left unchanged and `kept == total`.
  Reads/writes gzip by extension; R1/R2 count mismatch raises a clear error.
- `kraken2 extract`: options `--max-read-pairs INT` (min 1, default none) and
  `--seed INT` (default from `METADARKMATTER_SEED` via `core/random.get_seed`,
  i.e. 42). After extraction (and compression), subsample in place via a
  temporary file and atomic replace, keeping the existing output names so the
  pipeline needs no path changes. Write
  `{sample}_taxid{taxid}_subsample.json` next to the reads with total, kept,
  max_read_pairs and seed. Print one line: "Subsampled N of M read pairs
  (seed S)".
- `--skip-if-exists` keeps working: an existing extraction is not
  re-subsampled.
- Pipeline: `PipelineConfig.max_read_pairs: int | None` and `seed` are passed
  to the kraken-extract step argv; both are part of the step's fingerprint so
  changing them re-runs extraction and everything downstream. The manifest
  records them through the existing config dump. CLI options on
  `run pipeline` and `run cohort`. When `--extracted-reads-1` is supplied
  (extraction skipped), `--max-read-pairs` triggers a warning that it only
  acts on the extraction step.
- Docs: CLAUDE.md quick start note, docs/REFERENCE.md options, CHANGELOG.

## Left to the builder

- Exact record parsing (reuse helpers in `core/paired_reads.py` where they
  fit).
- Test fixtures (small FASTQ files in tmp_path).
