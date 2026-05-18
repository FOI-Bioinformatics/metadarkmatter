# Bayesian Calibration Pipeline

The Bayesian classifier uses six 2D Gaussians in (novelty,
placement_uncertainty) space. By default the per-category means and
sigmas are hand-tuned from the threshold table in
``ScoringConfig``. This pipeline replaces those hand-tuned values with
empirical fits to a labelled corpus and writes them as a config YAML.

This page is a walk-through; the assumptions the calibration relies on
(diagonal-covariance Gaussian, fixed prior modulation, linear
entropy→confidence map) are documented in
``docs/STATISTICAL_ASSUMPTIONS.md``.

## Inputs

A labelled TSV with three columns:

| Column | Type | Description |
|--------|------|-------------|
| ``novelty_index`` | float (0-100) | ``100 - top_hit_identity`` per read |
| ``placement_uncertainty`` | float (0-100) | ``100 - max(ANI between competing genomes)`` |
| ``true_category`` | string | one of the six Bayesian categories |

Categories are: ``Known Species``, ``Novel Species``, ``Novel Genus``,
``Species Boundary``, ``Ambiguous``, ``Unclassified``.

The categories must be assigned by an external ground truth (e.g.
simulated reads from genomes at known ANI distance, or
expert-labelled environmental reads). The calibration script
intentionally trusts the labels.

## Synthetic benchmark (smoke test)

For quick verification of the pipeline before plugging in real
labelled data, generate a synthetic corpus that samples each category
from the package's hand-tuned default Gaussians:

```bash
python scripts/build_synthetic_benchmark.py \
    --output benchmarks/synthetic.tsv \
    --per-category 500 \
    --seed 42
```

The output is a labelled TSV with 6×500 rows. It is *not* a substitute
for a real corpus; fitting to it will recover the defaults plus
sampling noise.

## Fitting

```bash
python scripts/calibrate_bayesian.py \
    --benchmark benchmarks/synthetic.tsv \
    --output configs/bayesian_calibrated.yaml \
    --min-samples 20 \
    --sigma-floor 0.5
```

The script prints a one-row-per-category summary table:

```
Category                  N     mu_N    mu_U   sg_N   sg_U   r(N,U)
------------------------------------------------------------------------
Known Species           500     2.01    0.81   1.81   0.87   -0.049
Novel Species           500    11.87    0.95   5.63   0.97   +0.005
...
```

Rows where ``|r(N,U)|`` exceeds ``--corr-threshold`` (default 0.3) are
flagged with a trailing ``!``. Those categories have meaningful
correlation between novelty and uncertainty, and the diagonal-
covariance Gaussian will under-fit them. They are good candidates for
a follow-up full-covariance treatment (not currently implemented in
the production classifier).

Categories with fewer than ``--min-samples`` rows are skipped, not
filled with defaults. The unspecified categories in the output YAML
will continue to use the hand-tuned defaults at classification time.

## Using the calibrated config

```bash
metadarkmatter score classify \
    --alignment sample.blast.tsv.gz \
    --ani ani_matrix.csv \
    --config configs/bayesian_calibrated.yaml \
    --output classifications.csv
```

``build_category_params`` in
``src/metadarkmatter/core/classification/bayesian.py`` consults
``config.bayesian.category_params`` first and falls back to the
hand-tuned defaults for any category the fitted dictionary does not
specify. This means a partially-fitted config is always safe to use.

## Entropy → confidence recalibration

The Bayesian classifier reports a 0-100 confidence score derived from
posterior entropy. By default the mapping is linear:
``(1 - entropy / log2(6)) * 100``. ``scripts/calibrate_entropy.py``
replaces this with an empirical curve fitted to the same labelled
corpus:

```bash
python scripts/calibrate_entropy.py \
    --benchmark benchmarks/synthetic.tsv \
    --base-config configs/bayesian_calibrated.yaml \
    --output configs/bayesian_calibrated.yaml \
    --bins 20
```

The script computes posterior entropy for every labelled row (using
the current Gaussian parameters), then fits an isotonic regression of
accuracy on entropy (or a per-bin monotone-enforced fallback when
scikit-learn is not installed). The output is a list of
``[entropy, confidence]`` knots written into
``bayesian.entropy_calibration``; the classifier interpolates between
them at lookup time and clips to ``[0, 100]``.

When ``entropy_calibration`` is ``None`` (the default), the linear
mapping is used unchanged. Partially-calibrated configs are safe:
``category_params`` and ``entropy_calibration`` are independent fields
and either may be set without the other.

## Known limitations

- The calibration treats the labels as ground truth; if the corpus is
  itself biased (e.g. ANI labels approximated from BLAST identity),
  the fit will inherit that bias.
- Only the means and standard deviations are fitted; the priors and
  the prior-modulation boosts are taken from ``BayesianConfig``
  defaults unless explicitly overridden in the base config.
- The diagonal-covariance Gaussian assumption is unchanged; categories
  flagged with high ``|r(N,U)|`` would benefit from a full-covariance
  treatment that is not currently implemented in the runtime
  classifier.
