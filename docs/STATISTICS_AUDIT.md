# Statistics and Calculations Audit

Status: 2026-05-26. Audit of the classifier math, the calibration
pipeline, and the proposed read-simulation experiment. Pair this with
`docs/STATISTICAL_ASSUMPTIONS.md` for the load-bearing modelling
assumptions and `docs/CALIBRATION_RESULTS.md` for the three-family
empirical validation results.

## Summary

| Area | Verdict |
|---|---|
| Bayesian 6-category math | Correct; log-space stable, priors normalised, boosts pre-normalisation |
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

Files: `src/metadarkmatter/core/classification/bayesian.py:314-643`,
`src/metadarkmatter/models/config.py:20-97`.

- Per-category 2-D Gaussian likelihoods are computed in log-space
  with pre-computed `log(2 * pi * sigma_n * sigma_u)` normalisation
  constants; standardised `z_n`, `z_u` are exponentiated as
  `exp(-0.5 * (z_n**2 + z_u**2))`. Numerically stable for the
  observed parameter ranges.
- Priors default to uniform (1/6) and are validated to sum to 1.0
  (`config.py:90-97`).
- Prior boost mechanics (`bayesian.py:616-624`) multiply
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

File: `src/metadarkmatter/core/classification/classifiers/vectorized.py:600-884`.

- `novelty_index = 100 - top_hit_identity` (line 795). Correct.
- Second-best hit gap uses
  `pl.col(...).slice(1, 1).first()` (commit 8be19f4) after the
  polars 1.37+ regression that broke `.get(1)`. The replacement is
  null-safe and verified on both polars 1.36 and 1.40.
- `placement_uncertainty = 100 - max_competing_ani` is correct;
  off-target reads receive a fixed `confidence_score = 10.0` and
  null uncertainty rather than the prior 0.0 (less misleading).
- Coverage weighting (linear, log, sigmoid; `vectorized.py:578-596`)
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
  validated for monotonic ordering in `config.py:244-256`.

### ANI matrix

File: `src/metadarkmatter/core/classification/ani_matrix.py`.

- Lookups are correct. Symmetric lookups assume the input nested
  dict is symmetric; the underlying NumPy array stores only one
  value per `(i, j)` cell as written.
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
  (`scripts/build_synthetic_benchmark.py:59-72`) samples from the
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
   YAML as a package default. `per_genome` mode emits a stderr
   warning pointing at `CALIBRATION_RESULTS.md`.
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
