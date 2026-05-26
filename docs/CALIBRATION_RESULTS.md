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

## v3: third family (g__Lactobacillus) and the full 3x3 matrix

A third corpus was built from g__Lactobacillus (100 GTDB
representatives, 25 holdout targets, **5 of 6 categories populated**:
Known Species, Novel Species, Species Boundary, Ambiguous,
Unclassified — every category except Novel Genus). dc-megablast was
used for the alignment step because default blastn was too slow for
the larger reference set; 2.7M alignments, 89k reads classified.

Empirical Gaussian fits (--min-samples 5000):

| Category | N | mu_N | mu_U | sigma_N | sigma_U | r(N, U) |
|---|---|---|---|---|---|---|
| Known Species | 12,286 | 4.82 | 2.99 | 6.61 | 5.49 | +0.118 |
| Novel Species | 43,236 | 10.37 | 5.12 | 8.19 | 9.10 | -0.012 |
| Species Boundary | 13,760 | 2.64 | 5.09 | 3.89 | 5.69 | +0.076 |
| Ambiguous | 13,951 | 13.80 | 7.89 | 6.87 | 10.29 | -0.064 |
| Unclassified | 5,873 | 17.57 | 9.72 | 7.05 | 12.22 | **-0.384 !** |

Second category to be flagged |r| > 0.3 (Unclassified, -0.384).

### The 3x3 cross matrix (top-1 accuracy)

| Calibration ↓ / Eval → | F | P | L |
|---|---|---|---|
| **Baseline** | **22.24%** | **26.99%** | **27.10%** |
| F-cal (v2) | 13.40% | 17.46% | 19.14% |
| P-cal (v2) | 13.03% | 24.30% | 18.82% |
| L-cal | 3.58% | 2.22% | 21.99% |

**Every off-baseline cell regresses.** No calibration configuration —
self or cross — improves accuracy over the baseline on any corpus.
L-cal is catastrophically bad on F and P (~3% and ~2% accuracy).

### Why L-cal destroys F and P accuracy

The fitted L Species Boundary Gaussian has ``mu_N = 2.64``: empirically
Species Boundary reads in Lactobacillus have low novelty (they come
from genomes with a close 95-96% ANI neighbour, so reads mostly
match well). The package default Species Boundary Gaussian has
``mu_N = 12.0`` — positioned mid-novelty as a catch-all.

When L-cal is applied to F or P alignments, every low-novelty read
(in F/P these are mostly Known Species truth) now matches the
relocated Species Boundary Gaussian and gets misclassified. The
overall accuracy collapses.

Similar story for L's Unclassified fit at ``mu_N = 17.57`` vs default
``30.0`` — pulled toward Novel Species territory.

## Revised conclusion

The current calibration approach — refitting per-category 2D Gaussian
means and sigmas from a labelled corpus — does **not** improve
classification accuracy on any of the three small corpora we tried.

The mechanism is that the per-genome ANI label does not map cleanly
to a per-read ``(novelty_index, placement_uncertainty)`` location.
Different categories' reads occupy overlapping regions of (N, U)
space because conserved/divergent gene content varies across the
reads of a single source genome. Fitting Gaussians to empirical (N, U)
distributions of differently-labelled reads pulls category means
into one another's threshold regions, breaking the classifier's
ability to discriminate.

The diagonal-covariance assumption is only part of the problem;
Lactobacillus's correlations were all |r| < 0.4 yet the fitted
config still regressed. The deeper issue is that the categories
themselves are defined by genome-level ANI thresholds rather than
by clusters in read-level (N, U) space.

### Concrete implications

1. **Do not ship a calibrated default.** None of the three calibrated
   YAMLs improve on the hand-tuned defaults.

2. **The shipment criteria in ``docs/CALIBRATION.md`` correctly
   reject all three candidates.** Keep them as is.

3. **The calibration scripts are valuable as diagnostic tools** —
   they reveal the (N, U) distribution per category, which is
   informative even when the fitted parameters cannot be shipped.
   Empirical r(N, U) flags categories where the diagonal Gaussian
   underfits (Pseudomonas Ambiguous -0.338, Lactobacillus
   Unclassified -0.384 in this run).

4. **A genuinely improved Bayesian classifier likely requires**
   either (a) labels defined directly in (N, U) space (clustering
   per-read metrics rather than per-genome ANI), or (b) a richer
   feature set than just (N, U) per read, or (c) acknowledging the
   current six categories are a per-genome taxonomy that cannot be
   recovered exactly from per-read evidence.

## Caveats

- All three corpora are small (10-25 holdout genomes). A much larger
  benchmark might support fitting categories empirically. But the
  3.6% accuracy from L-cal on F suggests the issue is structural
  rather than statistical.
- ECE comparisons across self/cross are confounded by entropy-curve
  collapse in the self-runs. Accuracy numbers are the reliable signal.
- Three families are enough to invalidate the "calibration improves
  accuracy" hypothesis. The "calibration improves ECE in exchange
  for accuracy" finding from v1 is still valid for self-runs;
  cross-runs offer no consistent benefit.

## Cross-family validation gate (operational)

Any future fitted YAML proposed for `configs/` must clear the gate
below. The gate is asymmetric by design: it is built to reject
overstated novelty, not to certify universal correctness. Headline
accuracy alone is not a sufficient pass criterion because calibration
can raise accuracy while simultaneously inflating the rate at which
true Known Species reads are called Novel Species or Novel Genus.

### Procedure

1. Pick three families spanning the family-size and divergence range.
   `f__Francisellaceae`, `f__Pseudomonadaceae`, `f__Lactobacillaceae`
   are the working defaults because they cover the failure cases
   already documented above.
2. For each family `F`, run
   `scripts/run_per_read_calibration.sh F corpora/F` to produce
   `corpora/F/bayesian_per_read.yaml` and the per-read metrics TSV.
3. For each test family `X != F`, run
   `metadarkmatter score classify --config corpora/F/bayesian_per_read.yaml`
   over `corpora/X` and run `metadarkmatter score evaluate
   --predictions corpora/X/metrics_with_F_config.tsv
   --truth-column true_category`.
4. Compare against the hand-tuned-defaults baseline for the same
   test family.

This yields a 3 x 3 matrix of (train_family, test_family) results;
ignore the diagonal (trivially good self-fit).

### Pass criteria

Let `baseline_X` be the hand-tuned defaults evaluated on family `X`,
and `cal_F_X` be the fit-from-F config evaluated on `X`. Define:

- `acc` = top-1 accuracy.
- `fnr_ks` = `P(predicted in {Novel Species, Novel Genus} | true = Known Species)`.
  This is the **overstated-novelty metric**.
- `fnr_ng` = `P(predicted = Novel Genus | true = Novel Species)`.

For every off-diagonal cell (F, X):

1. `acc(cal_F_X) >= acc(baseline_X) - 0.02` (no worse than 2pp).
2. `fnr_ks(cal_F_X) <= fnr_ks(baseline_X) * 1.10` (false-novel rate
   on Known Species rises at most 10% relative).
3. `fnr_ng(cal_F_X) <= fnr_ng(baseline_X) * 1.10`.
4. `acc(cal_F_X) >= acc(baseline_X)` on at least 2 of 3 off-diagonal
   cells.

Ship the YAML only if all four conditions hold. Pasting the 3 x 3
matrix into this document is the evidence requirement.

### Robustness checks before reading the matrix

- Seed sweep on one family: re-run the orchestrator with
  `METADARKMATTER_SEED` set to three different values. If fitted
  category means jitter by more than the fitted sigmas, the corpus
  is too small; increase per-target depth before doing the
  cross-family run.
- Per-category sigma vs centroid spacing: within each true category,
  the empirical (N, U) sigma should be at least 2x the gap between
  adjacent category centroids. Overconfident Gaussians are the
  mechanism by which calibration overstates novelty.

### Why the bar is set here

Three families is the minimum that lets you distinguish "this
calibration genuinely generalises" from "this calibration is just
less bad on the family it was fit on." Five would be better; the
information per added family drops fast and simulation cost grows
linearly. The 10% relative tolerance on `fnr_ks` / `fnr_ng` matches
the noise level seen across the existing self-runs in this document,
so a passing config is one whose novelty inflation is indistinguishable
from natural per-family variance.
