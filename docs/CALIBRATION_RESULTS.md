# Calibration Validation Results

This document records the empirical results of running the calibration
pipeline (``scripts/build_corpus.py`` → ``scripts/calibrate_bayesian.py``
→ ``scripts/calibrate_entropy.py`` → ``metadarkmatter score evaluate``)
against two small synthetic corpora. The findings shift the default
shipment strategy and inform the next experiments.

## Corpora

| Corpus | GTDB taxon | Genomes | Targets | Categories populated (per-genome) |
|---|---|---|---|---|
| Francisellaceae | f__Francisellaceae | 48 | 11 | Ambiguous (1), Novel Species (7), Unclassified (3) |
| Pseudomonas | g__Pseudomonas | 34 | 10 | Known Species (2), Novel Species (7), Ambiguous (1) |

Both corpora use the InSilicoSeq HiSeq simulation profile, ~10x depth,
75/25 stratified split, BLASTN alignment against the reference subset,
skani-derived reference ANI matrix. Francisellaceae uses the default
96% Known Species threshold; Pseudomonas uses ``--ani-known-species
95.0`` (GTDB).

## Empirical Gaussian fits

### Francisellaceae

| Category | N | mu_N | mu_U | sigma_N | sigma_U | r(N, U) |
|---|---|---|---|---|---|---|
| Novel Species | 27,736 | 9.01 | 4.51 | 7.02 | 7.37 | +0.081 |
| Ambiguous | 2,034 | 14.93 | 6.57 | 8.83 | 10.03 | +0.239 |
| Unclassified | 6,971 | 19.36 | 4.41 | 6.41 | 9.29 | -0.054 |

### Pseudomonas (genus)

| Category | N | mu_N | mu_U | sigma_N | sigma_U | r(N, U) |
|---|---|---|---|---|---|---|
| Known Species | 4,926 | 3.06 | 5.68 | 3.57 | 4.36 | -0.250 |
| Novel Species | 27,051 | 8.09 | 5.63 | 6.87 | 5.08 | -0.296 |
| Ambiguous | 4,514 | 4.24 | 7.07 | 4.77 | 3.62 | **-0.338 !** |

The Pseudomonas Ambiguous row is the first empirical correlation
|r| > 0.3 we have seen, supporting the documented limitation of the
diagonal-covariance assumption.

## Cross-corpus evaluation matrix

| Evaluated on → | Pseudomonas | Francisellaceae |
|---|---|---|
| Baseline (no calibration) | 26.99% | 22.24% |
| Self-calibrated | 20.33% | 14.23% |
| **Cross-calibrated** | **38.20%** | 18.58% |

Cross-calibration outperforms self-calibration in both directions
and outperforms the baseline on Pseudomonas. The Francisella cross
configuration applied to Pseudomonas gives +11 pp over baseline and
+18 pp over Pseudomonas's own self-fit.

## Interpretation

Three pieces of evidence point to the same explanation:

1. **Self-calibration on small corpora overfits.** Both self-fit
   entropy isotonic curves collapsed to near-zero confidence across
   the entropy range. The fitted Gaussians for under-sampled
   categories are too narrow.

2. **Cross-calibration regularises by retaining hand-tuned defaults
   for missing categories.** Francisella's fitted YAML leaves Known
   Species at the package default because the corpus has no Known
   Species samples; that default then matches Pseudomonas's 4,926
   Known Species reads better than Pseudomonas's own
   too-narrow fitted Gaussian.

3. **Asymmetric benefit.** Cross helped Pseudomonas (+11 pp vs
   baseline) more than Francisella (-3.7 pp vs baseline) because
   Pseudomonas's baseline was further from the truth to start with.

## Decisions

1. **Raise ``--min-samples``.** The current default of 100 admits
   fits on cohorts where overfitting dominates. Raise to 5000 or
   higher; under-supported categories should fall through to
   hand-tuned defaults.

2. **Pursue a single global default.** Calibration on a *combined*
   corpus from multiple families is the empirically motivated path
   for shipping ``configs/bayesian_default.yaml``.

3. **Revise the shipment criteria** in ``docs/CALIBRATION.md``. The
   strict "no top-1 accuracy regression" rule rejects exactly the
   configs that improve cross-family transfer. Replace with a
   cross-validation accuracy criterion measured on a held-out
   family.

## v2 follow-up: raising ``--min-samples`` did not preserve the cross-cal gain

After observing the v1 cross-cal wins, ``--min-samples`` was raised
from 100 to a default of 5000 and the cross matrix re-evaluated. The
expectation was that under-supported categories would fall back to
hand-tuned defaults and stabilise the result.

| v2 calibration → eval | Pseudomonas | Francisellaceae |
|---|---|---|
| Baseline | 26.99% | 22.24% |
| v2 self | 24.30% | 13.40% |
| v2 cross | 17.46% | 13.03% |

Compared to v1:

| Metric | v1 | v2 | Δ |
|---|---|---|---|
| Cross on Pseudomonas | 38.20% | 17.46% | **-20.74 pp** |
| Self on Pseudomonas | 20.33% | 24.30% | +3.97 pp |
| Self on Francisellaceae | 14.23% | 13.40% | -0.83 pp |

The v1 "cross beats self by 18 pp on Pseudomonas" result did **not**
survive raising ``--min-samples``. v2 self-cal on Pseudomonas
improved (closer to baseline, the more regularised behaviour we
wanted), but v2 cross-cal regressed sharply. The most plausible
explanation: the v1 cross win was driven by the spurious-helpful
fitted Ambiguous Gaussian from Francisella (a 2,034-sample fit
that happened to fit Pseudomonas's Ambiguous reads well). With v2
defaulting Ambiguous in both YAMLs, that accident is gone.

Revised conclusion: on these two small corpora, no calibration
configuration consistently beats the baseline. **The v1 cross
result was overfitting noise, not generalisation.** A third
family with a more balanced category distribution is now in
flight to test whether a calibration fit on a richer corpus
behaves differently.

## Caveats

- Both corpora are small (10-11 holdout genomes, ~36k reads each)
  and heavily skewed toward Novel Species.
- ECE comparisons across self/cross are confounded by the
  entropy-curve collapse in both self-runs (``ECE ≈ 1 - accuracy``
  when confidence is uniformly near zero). The accuracy numbers
  remain meaningful.
- The cross-family test was a 2x2; preliminary results from a
  100-genome g__Lactobacillus corpus (5 of 6 categories populated)
  will be added in a follow-up update.
