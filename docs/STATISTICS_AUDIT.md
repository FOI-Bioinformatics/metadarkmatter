# Statistics and Calculations Audit

Status: 2026-05-26. Audit of the classifier math, the calibration
pipeline, and the proposed read-simulation experiment. Pair this with
`docs/STATISTICAL_ASSUMPTIONS.md` for the load-bearing modelling
assumptions and `docs/CALIBRATION_RESULTS.md` for the three-family
empirical validation results.

Line references and code descriptions below were re-checked against
the code on 2026-11-02, after PR #13; changes made by that PR are
summarised in "Changes after this audit" at the end. The audit's
verdicts and experimental conclusions are unchanged.

## Summary

| Area | Verdict |
|---|---|
| Bayesian 6-category math | Correct; normalised 2-D Gaussian densities, priors validated, boosts pre-normalisation |
| Vectorised N / U computation | Correct, including the polars 1.37+ `slice(1,1).first()` fix |
| Legacy threshold cascade | Correct; thresholds sourced from literature |
| Coverage weighting | Mathematically sound; "empirical strength" parameter is unvalidated |
| ANI matrix | One latent failure mode fixed in this pass (symmetry now checked) |
| Calibration on per-genome labels | Demonstrably regresses accuracy on every test family |
| Calibration on per-read labels | Plausible path forward; gated by cross-family validation |

No critical bugs were found in the core classifier. One silent
failure mode in `ANIMatrix` (asymmetric input dicts overwriting one
direction) has been fixed in this audit pass. See "Fixes shipped"
below.

## Findings by area

### Bayesian 6-category classifier

Files: `src/metadarkmatter/core/classification/bayesian.py:288-644`,
`src/metadarkmatter/models/config.py:57-126`.

- Per-category 2-D Gaussian log-likelihoods are formed as
  `-log(2 * pi * sigma_n * sigma_u) - 0.5 * (z_n**2 + z_u**2)` (the
  normalisation constant is pre-computed) and exponentiated before
  normalisation; there is no log-sum-exp step. If the weighted sum
  underflows (< 1e-300) the read receives uniform posteriors. With the
  default parameters and N, U in [0, 100] this guard is not reached.
- Priors default to uniform (1/6) and are validated to sum to 1.0
  within 0.01 (`config.py:116-122`).
- Prior boost mechanics (`bayesian.py:611-624`) multiply
  `identity_gap_boost = 2.0` and `single_hit_boost = 1.5` into the
  Ambiguous prior *before* the per-category posterior is renormalised
  by the sum. This is the documented intent.
- Diagonal covariance assumption is defensible for most categories;
  see `STATISTICAL_ASSUMPTIONS.md` section 4 for the two categories
  that exceed `|r| > 0.3`.
- Entropy to confidence is a linear map
  `(1 - H / log2(6)) * 100` (`bayesian.py:235`), with optional
  piecewise-linear calibration via `BayesianConfig.entropy_calibration`.

### Vectorised classifier

File: `src/metadarkmatter/core/classification/classifiers/vectorized.py:421-1078`.

- `novelty_index = 100 - top_hit_identity` (line 896). Correct.
- Second-best hit gap uses
  `pl.col(...).slice(1, 1).first()` (commit 8be19f4) after the
  polars 1.37+ regression that broke `.get(1)`. The replacement is
  null-safe and verified on both polars 1.36 and 1.40.
- `placement_uncertainty = 100 - ANI(best, second genome)` under the
  default `uncertainty_mode = "second"` (`100 - max_competing_ani`
  with `uncertainty_mode = "max"`) is correct;
  off-target reads receive a fixed `confidence_score = 10.0` and
  null uncertainty rather than the prior 0.0 (less misleading).
- Coverage weighting (linear, log, sigmoid; `vectorized.py:220-253`)
  is mathematically straightforward. The `strength = 0.5` default is
  marked "empirical" in `STATISTICAL_ASSUMPTIONS.md` section 6 with no
  cited validation. Treat as a tunable, not as physics.

### Legacy threshold cascade

Files: `src/metadarkmatter/core/classification/thresholds.py`,
`core/constants.py`, `core/protein_constants.py`.

- ANI species boundary 95-96%: Jain et al. 2018, Goris et al. 2007.
- Genus boundary ~80%: Konstantinidis & Tiedje 2005, Riesco & Trujillo 2024.
- AAI genus 58-65%: Riesco & Trujillo 2024.
- `apply_legacy_thresholds` is a pure rule application; thresholds are
  validated for ordering in `config.py:511-560`.
- The cascade is used for `legacy_call` and `score sensitivity` only;
  the primary `taxonomic_call` is the Bayesian MAP, whose decision
  boundaries do not coincide with these thresholds (see
  `docs/METHODS.md` section 2.5).

### ANI matrix

File: `src/metadarkmatter/core/classification/ani_matrix.py`.

- Lookups are correct. Symmetric lookups assume the input nested
  dict is symmetric; the underlying NumPy array stores one value per
  `(i, j)` cell as written. (At the time of the audit the vectorized
  classifier read only the upper triangle, so a pair supplied only in
  the lower triangle was replaced by the default ANI; since PR #13
  one-directional pairs are mirrored.)
- **Prior latent bug**: an asymmetric input dict would overwrite one
  direction silently, corrupting placement uncertainty for any reads
  whose competing genomes hit the asymmetric pair. **Fixed in this
  pass** - see "Fixes shipped".
- Diagonal is set to 100.0 explicitly. Missing pairs return
  `default_ani` (70%), documented in `STATISTICAL_ASSUMPTIONS.md`
  section 11.

### Calibration scripts

Files: `scripts/calibrate_bayesian.py`, `scripts/calibrate_entropy.py`,
`src/metadarkmatter/core/classification/evaluation.py`,
`docs/CALIBRATION_RESULTS.md`.

- Gaussian fit math is correct (`np.mean`, `np.std(ddof=1)`,
  `np.corrcoef`).
- The empirical failure is structural, not numerical. Per
  `STATISTICAL_ASSUMPTIONS.md` section 13, per-genome category labels
  do not map cleanly to per-read `(N, U)` clusters: a single source
  genome produces reads from conserved (high-identity) and divergent
  (low-identity) regions, all carrying the same per-genome label.
  Fitting Gaussians on this collapses the means toward overlapping
  territory and breaks discrimination. The three-family table in
  `CALIBRATION_RESULTS.md` shows every fitted YAML regresses accuracy
  on every test corpus.
- The synthetic benchmark generator
  (`scripts/build_synthetic_benchmark.py:55-72`) samples from the
  hand-tuned defaults, which makes it useful for shape testing but
  **circular** for validating calibration.

### Read simulation as a calibration path

`scripts/build_corpus.py` already emits a per-target InSilicoSeq
script, and `scripts/build_metrics_tsv.py --label-mode per_read`
already attaches per-read truth labels via ANI from the source genome
to the best-match reference. The remaining gap was a single
orchestrator; that orchestrator now exists at
`scripts/run_per_read_calibration.sh`.

Per-read labelling is a legitimate attempt to escape the section-13
problem because the label is finally assigned at the unit of
inference. It does **not** automatically solve the problem - even
with per-read labels, fitted Gaussians can still pull means together
on small or skewed corpora. The cross-family validation gate in
`docs/CALIBRATION_RESULTS.md` (calibrate on one family, evaluate on
two others; no family worse than -2pp vs the hand-tuned defaults) is
the explicit defence against overfitting and against overstating
novelty.

## Fixes shipped in this audit pass

1. **ANI matrix symmetry check.**
   `ANIMatrix.__init__` and `ANIMatrix.from_file` now accept
   `symmetry_check={"off", "warn", "strict"}` with default `"warn"`.
   Tolerance defaults to 0.5 ANI units. Behaviour:
   - `"warn"` logs a single message naming up to five disagreeing
     pairs.
   - `"strict"` raises `ValueError`.
   - `"off"` preserves prior behaviour for callers that need it.
   Unit tests in `tests/unit/test_ani_matrix_class.py::TestSymmetryCheck`.
2. **`--strict-ani` flag on `mdm score classify`**, threaded into
   both `ANIMatrix.from_file` call sites (dry-run preview and main
   classify path).
3. **Guardrails on `scripts/calibrate_bayesian.py`.**
   `--label-mode {per_genome, per_read}` is now required so the
   operator must explicitly declare the source of the truth labels.
   `--i-have-validated` is required when `--output` writes inside
   `configs/`, refusing to silently install an unvalidated fitted
   YAML as a package default. `per_genome` mode prints a warning
   (to stdout) pointing at `CALIBRATION_RESULTS.md`.
4. **`scripts/run_per_read_calibration.sh`** - end-to-end
   orchestrator that wires `build_corpus` -> `simulate.sh` ->
   reference BLAST DB -> classify -> `build_metrics_tsv --label-mode
   per_read` -> `calibrate_bayesian --label-mode per_read` ->
   re-classify -> `score evaluate` and prints a default-vs-calibrated
   accuracy delta. Writes the fitted YAML outside `configs/` by
   design so the guardrail is not tripped.

## How to use the experiment

```bash
METADARKMATTER_SEED=42 scripts/run_per_read_calibration.sh \
    f__Francisellaceae corpora/francisella 16
cat corpora/francisella/summary.txt
```

Repeat on `f__Pseudomonadaceae` and `f__Lactobacillaceae` to fill out
the cross-family validation matrix in
`docs/CALIBRATION_RESULTS.md`. Only after that matrix shows the
fitted config beating the hand-tuned defaults on at least two of
three families, with no family worse than -2pp, should the YAML be
proposed for `configs/` (and even then `--i-have-validated` must be
explicitly passed).

## Better options beyond per-read calibration

These are out of scope for this audit pass but are the natural next
investigations if the per-read fit also fails to clear the gate:

- **Data-driven boost factors.** Replace the heuristic
  `identity_gap_boost = 2.0` and `single_hit_boost = 1.5` with logistic
  regression coefficients fitted on the per-read corpus. Caveat: same
  cross-family generalisation question applies.
- **Full-covariance Gaussians for the two categories with `|r| > 0.3`.**
  `STATISTICAL_ASSUMPTIONS.md` section 4 already flags the candidates;
  the model code change is small.
- **Conformal prediction over the novelty score.** Distribution-free
  coverage guarantees would give a principled "this read is novel
  with confidence at least 1 - alpha" statement and explicitly cap
  the rate of overstated novelty.

None of these should ship without a held-out validation family.

## Changes after this audit (PR #13)

PR #13 changed classifier behaviour in ways that bear on the findings
above (see `CHANGELOG.md` for the full list):

- **Reverse-strand coordinates.** Hits with `qstart > qend`
  (MMseqs2/DIAMOND minus strand) are normalised at parse time and in
  the vectorized filters. Previously the true best hit received a
  negative coverage and was down-weighted or filtered.
- **Second-genome tie-breaking.** The genome used for
  `uncertainty_mode = "second"` is chosen with the same total order as
  best-hit selection (weighted bitscore, `pident`, accession). The
  earlier unstable sort made `placement_uncertainty` depend on batch
  size.
- **ANI matrix.** Stored as float64; one-directional pairs are
  mirrored; 0-1 scaled matrices are rejected.
- **Adaptive thresholds** ignore uncomputed ANI pairs and preserve the
  base configuration.
- **`min_alignment_fraction`** uses `qlen` when available.
- **QC low-confidence fraction** uses the `low_confidence` flag (or
  the 0-100 `confidence_score` scale) instead of comparing the
  0-100 score with 0.5.

A post-PR re-check also documented, without changing code or defaults,
that the Bayesian MAP boundaries differ from the hard thresholds (e.g.
Known -> Novel Species at N ~ 6.2% rather than 4%) and that the
Bayesian call does not use the inferred uncertainty of single-hit reads
(`docs/METHODS.md` section 2.5).
