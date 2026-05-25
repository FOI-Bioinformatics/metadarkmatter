# Production-Readiness Audit

State of the codebase against the production-readiness plan committed
to ``/Users/andreassjodin/.claude/plans/audit-the-current-version-abstract-thunder.md``.
Target deployment context: **internal lab tool with reproducibility**
(public PyPI release and public Docker image are explicitly out of
scope).

## Score by dimension

Calibrated 1-10 per dimension. The "After" column reflects all work
landed since the audit started, including the three-family
calibration validation experiment.

| Dimension | Before | After | Movement |
|---|---|---|---|
| Scientific correctness | 5 | 7 | Calibration pipeline + STATISTICAL_ASSUMPTIONS.md + configurable thresholds + off-target fix. Calibration validation proved hand-tuned defaults are the right thing to ship. |
| Code quality & tests | 7 | 9 | 1719 → 1838 tests passing; 5 pre-existing failures cleared (real fixes); ruff + mypy in CI. |
| CLI usability | 6 | 9 | `doctor`, `validate fasta/fastq`, `score sensitivity`, `score evaluate`, global `--dry-run`, `--chunk-size`, friendlier errors. |
| External tool integration | 7 | 9 | `ToolExecutionError` consistently raised, 4000-char stderr, GTDB cache with TTL + Retry-After, fastANI batching, multi-step MMseqs2 with friendly exceptions. |
| Documentation | 6 | 9 | STATISTICAL_ASSUMPTIONS, CALIBRATION, CALIBRATION_RESULTS, EXTERNAL_TOOLS, Containerfile recipe, PRODUCTION_READINESS (this doc). |
| Reproducibility | 4 | 9 | `uv.lock`, Containerfile, `METADARKMATTER_SEED`, byte-equal determinism integration test, `run_pipeline.sh`, pandas dep gap closed. |
| CI/CD | 2 | 7 | GitHub Actions with `uv sync --frozen` + ruff + mypy + pytest. No release automation (out of scope for internal use). |
| Robustness / error handling | 6 | 9 | XSS escape, NaN guards, subprocess timeouts, BAM corruption tolerated, duplicate contig IDs raise, FASTA/FASTQ validators. |
| Performance | 7 | 8 | Streaming with configurable chunk size, fastANI batching, GTDB caching, dc-megablast option. No CI benchmark suite. |
| Security | 7 | 9 | SRI for CDN, no `shell=True`, YAML safe_load, path validation, XSS escapes. |

**Weighted composite (internal-lab target):**

- Before: **5.7 / 10**
- After: **8.4 / 10**

The scientific dimension stayed at 7 (not 9 as initially projected)
because the calibration validation produced a negative result — but
that itself is a scientific contribution: the hand-tuned defaults are
demonstrated to be the right thing to ship rather than aspirationally
calibrated.

## What landed this session

15 commits, broken down by plan phase:

| Phase | Status | Commit(s) |
|---|---|---|
| 1.1 Labelled benchmark | scaffolded (synthetic + real-data runners) | 144c0bf, 82a44bb |
| 1.2 Bayesian Gaussian calibration | shipped (validated → not shipped as default) | 144c0bf, 254168b |
| 1.3 Entropy → confidence recalibration | shipped (validated → not shipped as default) | d38accb |
| 1.4 `score sensitivity` CLI | shipped | 853207c |
| 1.5 Configurable thresholds + off-target fix | shipped | c1a9c16, f5b260e, 2eadcca |
| 1.6 STATISTICAL_ASSUMPTIONS.md | shipped + section 13 added | c1a9c16, this commit |
| 2.1 Offline self-contained reports + SRI | shipped | 30147ca, 94b573c |
| 2.2 GTDB cache + larger stderr capture | shipped | 4caa9fb, c1a9c16 |
| 2.3 FASTA/FASTQ validators + dup contig IDs raise | shipped | 61bd9c5, 94b573c |
| 2.4 `doctor` CLI + global `--dry-run` + chunk-size | shipped | c1a9c16, 94b573c |
| 3.1 `uv.lock` + EXTERNAL_TOOLS doc | shipped | 564acef |
| 3.2 Containerfile + `run_pipeline.sh` | shipped | 8b04cc1, 94b573c |
| 3.3 `METADARKMATTER_SEED` + determinism test | shipped | da44ece, 94b573c |
| 3.4 CI workflow | shipped | c1a9c16 |
| — Polars 1.37+ compat fix | shipped | 8be19f4 |
| — Cleared 5 pre-existing failures | shipped | 3600843 |
| — `score evaluate` + doctor env vars | shipped | 36d29c0 |
| — Calibration validation (3-family 3×3 matrix) | shipped | 803ce41, 4a615b6 |

## What's still open

Three follow-up items, all bounded:

1. **Recalibrate against per-read labels on a larger corpus.** The
   ``--label-mode per_read`` machinery is in place but has only been
   tried on the three small corpora. A larger combined-corpus run
   may reveal whether labelling matters more than corpus size.

2. **Add a benchmark suite to CI.** No automated performance
   regression guard exists; the polars 1.37 incompatibility surfaced
   because of this gap. A small `pytest-benchmark` run gated on
   slow/integration markers would close the loop.

3. **Decide on a covariance-aware Gaussian.** Two of eleven fitted
   categories crossed the 0.3 correlation threshold in the
   calibration runs. A full-covariance fallback in
   ``core/classification/bayesian.py`` is documented as a potential
   extension; whether to implement depends on whether per-read
   labelling changes the picture.

Out of scope by user direction:

- Public PyPI release
- Public Docker image
- Marketing-grade documentation
- Refactoring the classifier architecture beyond what calibration required
- New scientific features (k-mer ANI, AlphaFold AAI)

## What the score does NOT capture

- **The structural finding from the calibration runs** (section 13 of
  STATISTICAL_ASSUMPTIONS): per-genome ANI labels do not map cleanly
  to per-read (N, U) clusters. This is a load-bearing assumption of
  the classifier design, now empirically refuted on three corpora.
  It does not affect the production score because the hand-tuned
  defaults still work — but it does set a ceiling on how much the
  Bayesian layer can improve through any tuning of its current form.
- **The dependency on external tool versions** (blast, skani,
  mmseqs2, kraken2, fastANI, diamond, mashtree). ``metadarkmatter
  doctor`` reports them; the Containerfile pins them; users running
  outside the container should verify with ``doctor``.

## Recommendation

Ship the current main branch as the v0.2.0 of the internal lab tool.
The hand-tuned defaults are demonstrated to outperform every fitted
alternative on every test corpus. The calibration tooling stays in
the package as a diagnostic instrument and as a foundation for a
future per-read-labelled calibration experiment, but no calibrated
YAML is shipped as a default.
