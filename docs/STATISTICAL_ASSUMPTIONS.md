# Statistical Assumptions

This document enumerates the statistical and modelling assumptions baked
into the metadarkmatter classifier, the sources for the default
thresholds, and the failure modes that arise when those assumptions
break. It is intended as a companion to `docs/METHODS.md` and is the
canonical reference for users who need to defend or audit a
classification result.

## 1. ANI species boundary at 95-96%

**Assumption.** Two genomes belong to the same prokaryotic species when
their average nucleotide identity (ANI) is at least 95-96%.

**Source.** Jain et al. 2018 (*Nature Communications*) demonstrated a
bimodal ANI distribution with a strong gap around 95% for >90,000 prokaryotic
genome pairs. Goris et al. 2007 (*PNAS*) had earlier proposed the same
threshold as a proxy for 70% DNA-DNA hybridisation. GTDB uses 95%.

**Defaults.** `ScoringConfig.same_species_ani_threshold = 95.0`;
`novelty_known_max = 4.0` (i.e. `pident > 96%`).

**Failure modes.** The 95% boundary is an empirical mode of the genome-pair
distribution, not a universal law. Several taxa break it: some
*Mycobacterium*, *Streptomyces*, and *Pseudomonas* clades show within-species
ANI well below 95%; some *Bordetella* and *Yersinia* species are barely
above. For these clades, set `same_species_ani_threshold` and the
corresponding novelty thresholds explicitly via `--config`, do not rely
on the defaults.

## 2. ANI genus boundary at ~80%

**Assumption.** Two genomes belong to the same genus when their ANI is at
least ~80%.

**Source.** This is a convention rather than a sharp empirical mode. The
ANI signal becomes unreliable below ~80% because of saturation in shared
genomic content; AAI is preferred at this divergence level
(Konstantinidis & Tiedje 2005, Riesco & Trujillo 2024).

**Defaults.** `ScoringConfig.same_genus_ani_threshold = 80.0`.

**Failure modes.** Genus circumscriptions vary across phyla. Reads
falling between 75-85% ANI to the closest genome should be treated as
uncertain at the genus level; use AAI (`--alignment-mode protein`) for
firm genus calls.

## 3. AAI genus boundary at 58-65%

**Assumption.** Two genomes belong to the same genus when their AAI
falls in 58-65%.

**Source.** Riesco & Trujillo 2024 reviewed AAI distributions across
>10,000 genera and recommended 58-65% as the genus boundary range, with
the lower end appropriate for fast-evolving clades.

**Defaults.** `aai_genus_boundary_low = 58.0`,
`aai_genus_boundary_high = 65.0`.

**Failure modes.** Like ANI, this range is a population summary, not a
per-clade guarantee. For lineages where the genus has been recently
redefined (e.g. the split of *Clostridium*), set explicit thresholds.

## 4. Independent 2D Gaussians per Bayesian category

**Assumption.** For each of the six categories
(`known_species`, `novel_species`, `novel_genus`, `species_boundary`,
`ambiguous`, `unclassified`), the joint likelihood over
(novelty, placement_uncertainty) is a 2D Gaussian with diagonal covariance.

**Source.** Convenience; the categories are reasonably localised in
(N, U) space, and a Gaussian gives a closed-form log-likelihood. There is
no peer-reviewed precedent for this exact form in microbial classification.

**Defaults.** Means and sigmas are hand-set in
`core/classification/bayesian.py:114-163` based on the centre and half-width
of each category's threshold region.

**Failure modes.**
- *Correlation*: novelty and placement uncertainty are not strictly
  independent in real data; reads with high N often have high U because few
  close competitors exist. The diagonal-covariance Gaussian therefore
  underestimates the joint likelihood for off-diagonal mass.
- *Tails*: the Gaussian has light tails, so reads in the long tail of the
  novelty distribution (e.g. spurious high-N hits from sequencing errors)
  get suppressed posteriors that overstate the "Ambiguous" category.
- *Calibration*: the means and sigmas are not fitted to labelled data, so
  posterior probabilities are not calibrated — they preserve the ranking
  of categories but should not be used as probabilities in formal
  inference without first running the calibration pipeline.

**Empirical evidence (added after calibration validation).** Three
synthetic-but-labelled corpora (f__Francisellaceae, g__Pseudomonas,
g__Lactobacillus) ran the full calibration pipeline; the per-category
empirical Pearson correlations are summarised below. Values |r| > 0.3
(flagged with `!`) violate the diagonal-covariance assumption.

| Category | F r(N,U) | P r(N,U) | L r(N,U) |
|---|---|---|---|
| Known Species | — | -0.250 | +0.118 |
| Novel Species | +0.081 | -0.296 | -0.012 |
| Species Boundary | — | — | +0.076 |
| Ambiguous | +0.239 | **-0.338 !** | -0.064 |
| Unclassified | -0.054 | — | **-0.384 !** |

Two of the eleven fitted categories cross the 0.3 threshold. The
diagonal-Gaussian assumption is empirically defensible for most
categories but materially underfits at least two; a full-covariance
extension would be warranted if calibration on a larger corpus were
revived. See ``docs/CALIBRATION_RESULTS.md`` for the full matrix and
the deeper structural issue documented in section 13 below.

## 5. Prior modulation for Ambiguous

**Assumption.** The prior on the Ambiguous category is multiplied by
`identity_gap_boost` (default 2.0) when `identity_gap < threshold` and
`num_hits > 1`, and by `single_hit_boost` (default 1.5) for single-hit
reads with novelty in the novel-species range or higher.

**Source.** Heuristic. These boosts capture the intuition that small
identity gaps and single-hit evidence are weak signals.

**Failure modes.** The boost factors are not data-driven. On datasets
dominated by single-hit reads (low-coverage environmental samples), the
single-hit boost can dominate and inflate the Ambiguous fraction.

## 6. Linear coverage weighting

**Assumption.** A hit's effective bitscore is the raw bitscore multiplied
by `min + (max - min) * coverage`, where `min = 1 - strength`,
`max = 1 + strength`, and `strength` defaults to 0.5.

**Source.** Empirical: short conserved-domain hits inflate bitscores
relative to long, divergent hits. Linear weighting is the simplest mode
that produces a monotone correction without distorting the score
distribution at low strength values.

**Defaults.** `coverage_weight_mode = "linear"`,
`coverage_weight_strength = 0.5`. Weight ranges from 0.5x at zero coverage
to 1.5x at full coverage.

**Failure modes.** For target databases where most hits are full-coverage
(e.g. closely related references), the weighting has little effect; for
fragmented contigs or partial assemblies, the weight magnitude depends on
how coverage is computed (qlen-based vs alignment-based). When in doubt,
benchmark with `coverage_weight_mode = "none"` and compare.

## 7. Inferred uncertainty for single-hit reads

**Assumption.** A read with exactly one ANI-matrix hit can have its
placement uncertainty inferred from its novelty index using a piecewise-
linear function (slopes 0.5, 1.0, 1.5 across novelty bands).

**Source.** Construction: reads at low novelty have a small uncertainty
band, reads near the species boundary have a wider band, and reads at the
genus boundary have the widest band. The slopes are not fitted.

**Failure modes.** Inferred uncertainty is a smooth proxy for the missing
empirical uncertainty; it cannot recover real ANI competition that was
not measured. Treat single-hit Bayesian calls as lower-confidence than
multi-hit calls.

## 8. Adjacent novelty-band merging

**Assumption.** Two clusters that share `(representative_genome,
taxonomic_call)` and occupy adjacent novelty bands are merged when their
novelty ranges overlap within `merge_band_tolerance` percentage points.

**Source.** Heuristic to avoid hard band edges inflating cluster counts.

**Defaults.** `merge_band_tolerance = 2.0`.

**Failure modes.** A larger tolerance collapses biologically distinct
sub-populations into a single cluster; a smaller tolerance fragments a
single population into adjacent-band singletons.

## 9. Off-target classification

**Assumption.** A read whose best in-family bitscore is less than
`family_ratio_threshold * best_overall_bitscore` is labelled
`Off-target` and given a low fixed confidence; its placement uncertainty
is reported as `null` because the matching genome is outside the ANI
matrix.

**Source.** Empirical heuristic, intended for broad-database alignments
where a target family has been declared.

**Defaults.** `family_ratio_threshold = 0.8`,
off-target `confidence_score = 10.0`, `low_confidence = True`.

**Failure modes.** The threshold is not calibrated against false-positive
rates; reads near the threshold should be treated as ambiguous rather
than confidently off-target.

## 10. Entropy-to-confidence linear scaling

**Assumption.** Posterior entropy maps linearly to a 0-100 confidence
score via `(1 - entropy / log2(6)) * 100`.

**Source.** Construction: maximum entropy under six categories is
`log2(6)`, so the formula maps a uniform posterior to 0 and a delta
posterior to 100.

**Failure modes.** Linear scaling assumes uniform information content
across the entropy range. Empirically, classifier accuracy as a function
of entropy is non-linear; the linear mapping over-rewards mid-range
posteriors. The Phase 1 calibration pipeline replaces this with an
isotonic or logistic fit against a labelled holdout set.

## 11. ANI matrix completeness

**Assumption.** Pairs missing from the ANI matrix are treated as 0.0% ANI
(i.e. unrelated) rather than as missing data.

**Source.** Construction: missing pairs are typically those below the
ANI tool's reporting threshold, so 0.0 is the conservative substitute.

**Failure modes.** When the ANI matrix is computed only over
representatives, missing pairs may reflect computational shortcuts rather
than true non-similarity. For representative-only matrices, ensure the
representative mapping in the metadata file is comprehensive.

## 12. Reproducibility

**Assumption.** Given the same input alignment, ANI matrix, metadata,
and config, the classifier produces deterministic output.

**Verification.** Sorting in the vectorized classifier uses
`(weighted_bitscore, pident, genome_name)` for stable tie-breaking.
Random-state-dependent components (adaptive GMM, visualisation
sampling) accept a `random_state` seed.

**Failure modes.** External tools (FastANI, MMseqs2, Diamond) may be
non-deterministic across runs due to multi-threading; pin tool versions
and thread counts when bit-identical reruns are required.

## 13. Per-genome labels map cleanly to per-read (N, U) clusters

**Assumption (the load-bearing one for calibration).** A read drawn
from a source genome labelled with a per-genome category C tends to
produce ``(novelty_index, placement_uncertainty)`` values that fall
near the same locus in (N, U) space as the package's hand-tuned
Gaussian for C. Equivalently: the six categories form well-separated
clusters in (N, U) space when reads are grouped by their source
genome's per-genome ANI category.

**Source.** Implicit in the design of the Bayesian classifier — if
this assumption held, calibrating the per-category Gaussian means
from labelled data would be an unambiguously useful exercise.

**Empirical verdict: REFUTED on three synthetic corpora.** Running
the calibration pipeline end-to-end on f__Francisellaceae,
g__Pseudomonas, and g__Lactobacillus and evaluating the resulting
configs in a 3x3 cross-family matrix shows every fitted YAML
regresses accuracy versus the hand-tuned defaults, in some cases
catastrophically (Lactobacillus's calibrated config drops
Francisellaceae accuracy from 22% to 3.6%).

**Mechanism.** A single source genome produces reads from both
conserved and divergent regions. Conserved-region reads from a
"Novel Species" genome match a reference at high identity (low N);
divergent-region reads from the same genome match at low identity
(high N). The per-genome label is "Novel Species" for all of them,
but their (N, U) values span a wide range that overlaps with Known
Species, Species Boundary, and even Conserved Region territory.

Concrete example from Lactobacillus: "Species Boundary" target
genomes (best reference ~95-96% ANI) produce reads whose mean
novelty is ``mu_N = 2.64``, very close to "Known Species" territory.
The package default Gaussian places Species Boundary at
``mu_N = 12.0`` (mid-range, catch-all). Replacing the default with
the empirical 2.64 reassigns every low-novelty read in another
corpus (which is mostly Known Species truth) to Species Boundary,
collapsing overall accuracy.

**Failure modes.**
- *Calibration is not a tuning lever*. Fitting per-category Gaussian
  means from per-genome-labelled data does not improve accuracy on
  any corpus tested.
- *Cross-corpus transfer is worse than self-transfer in the L→F and
  L→P directions*. The richer Lactobacillus fit happens to relocate
  Gaussians furthest from the defaults, which produces the largest
  damage when applied to a different corpus.

**Practical consequence.** Ship the hand-tuned defaults. Use the
calibration scripts (``scripts/calibrate_bayesian.py``,
``scripts/calibrate_entropy.py``) as diagnostic tools — the per-
category empirical mean, sigma, and correlation r(N, U) are
informative for understanding category overlap — but do not commit
the fitted YAMLs as a default config.

**Open question.** A more faithful labelling scheme might use
per-read truth derived from ``ANI(source_target, best_match_genome)``
for each read individually rather than the source genome's
per-genome category. ``scripts/build_metrics_tsv.py --label-mode
per_read`` implements this; preliminary results on small corpora
(see ``docs/CALIBRATION_RESULTS.md``) show it produces more
internally consistent metrics but does not change the structural
finding. A targeted experiment with per-read labels on a larger
corpus would be the natural next step if revisiting calibration.
