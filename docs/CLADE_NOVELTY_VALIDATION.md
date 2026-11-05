# Clade-novelty rule: leave-one-out validation

This document records the validation of the clade-novelty Stage 2 rule
(`ScoringConfig.clade_novel_species`, on by default; `score classify --no-clade-novel-species` disables it).
The rule is described in [METHODS.md](METHODS.md#105-clade-novelty-rule).

## Problem

Placement uncertainty U measures how unsure the classifier is about **which**
reference a read belongs to. In a clade of close relatives (several species
of one genus at 85-95% ANI to each other), U is high even for reads that are
clearly divergent from **every** reference. Such reads were called Ambiguous
or Species Boundary instead of Novel Species. In the Francisellaceae data from
the 2026-11 audit, only 12-14% of reads from a held-out species at 91% ANI to
its nearest reference were called Novel Species; the missed reads had median
novelty 8-9% to both the best and the second-best reference.

## Rule

After Stage 1/2 classification (nucleotide mode only), a read is called Novel
Species when:

- `taxonomic_call` is Ambiguous or Species Boundary,
- `ambiguity_scope` is `within_species` or `within_genus`,
- `novelty_index >= 5` and `100 - second_hit_identity >= 5`,
- `novelty_index < 18.6`.

`clade_novelty_min` (5.0) lies between the legacy Known Species threshold (4%)
and the measured Bayesian Known/Novel Species boundary at U = 0 (about 6.2%);
`clade_novelty_max` (18.6) is the measured Bayesian Novel Species / Novel
Genus boundary at U = 0. Posteriors, entropy, `placement_uncertainty` and
`ambiguity_scope` are not changed; the `clade_novel_species` column records
which reads the rule changed.

## Method

`scripts/benchmark_leave_one_out.py` (seed 42):

1. Simulate 400 single-end 150 bp reads per genome (uniform substitution
   errors, rate 0.005).
2. Compute the ANI matrix with fastANI (`mdm ani compute --backend fastani`)
   and align all reads against all genomes with MMseqs2.
3. For each genome g: drop g's hits and g's ANI row/column, then classify g's
   reads. Truth from the full ANI matrix: nearest ANI >= 96 Known Species,
   >= 80 Novel Species, otherwise Novel Genus.
4. Control: nothing removed; every read's genome is in the database, so any
   Novel Species / Novel Genus call is a false-novel call.

Metrics are read-level. "Novel-species recall" uses genomes with nearest ANI
in [80, 94); genomes in [94, 96) are reported separately because Species
Boundary is a defensible call there. "Any-novel" counts Novel Species or Novel
Genus as novel and uses all genomes whose truth is not Known Species.

## Genome sets

All genomes are GTDB species representatives.

- **Francisellaceae** (design set, 12 genomes): GCF_000156715.1,
  GCF_000815225.1, GCF_000833475.1, GCF_001042525.2, GCF_001653955.1,
  GCF_001880225.1, GCF_001885235.1, GCF_003290445.1, GCF_007923265.1,
  GCF_009823375.1, GCF_012222825.1, GCF_019852205.1.
- **Lactobacillaceae** (held out, 16 genomes; 10 novel-species and 6
  novel-genus truths): GCF_000014425.1, GCF_000159355.1, GCF_000160835.1,
  GCF_000160855.1, GCF_000193635.1, GCF_000372485.1, GCF_000829035.1,
  GCF_000829055.1, GCF_001434365.1, GCF_001436025.1, GCF_003096575.1,
  GCF_007991875.1, GCF_013394085.1, GCF_018987235.1, GCF_034298135.1,
  GCF_900636965.1.
- **Enterobacteriaceae** (held out, 16 genomes; 13 novel-species and 3
  novel-genus truths; dense Escherichia, Salmonella, Klebsiella,
  Enterobacter and Citrobacter clades): GCF_000006945.2, GCF_000007365.1,
  GCF_000025565.1, GCF_000026225.1, GCF_000215745.1, GCF_000252995.1,
  GCF_000742135.1, GCF_000759775.1, GCF_000828055.2, GCF_001521715.1,
  GCF_001598695.1, GCF_002075345.1, GCF_003697165.2, GCF_011064845.1,
  GCF_017299535.1, GCF_901472495.1.

## Results

Default configuration -> with the clade-novelty rule.

| Metric | Francisellaceae | Lactobacillaceae | Enterobacteriaceae |
|---|---|---|---|
| Novel-species recall, nearest ANI 80-94% | 0.483 -> 0.695 | 0.521 -> 0.588 | 0.467 -> 0.632 |
| Novel Species precision | 0.930 -> 0.917 | 0.922 -> 0.926 | 0.955 -> 0.936 |
| Any-novel recall | 0.350 -> 0.517 | 0.794 -> 0.852 | 0.507 -> 0.656 |
| Any-novel precision | 0.988 -> 0.989 | 0.995 -> 0.994 | 0.992 -> 0.994 |
| Novel Genus recall | 0.318 -> 0.318 | 0.549 -> 0.549 | 0.357 -> 0.357 |
| Novel-genus reads called Novel Species | 0.161 -> 0.323 | 0.205 -> 0.216 | 0.191 -> 0.398 |
| Control false-novel reads | 17/4787 -> 22/4787 | 12/6315 -> 14/6315 | 19/6393 -> 20/6393 |

Observations:

- Novel-species recall improves in every family (+7 to +21 points). The gain
  is largest in dense clades (Enterobacteriaceae: Klebsiella oxytoca
  0.27 -> 0.70, K. aerogenes 0.33 -> 0.69, Escherichia albertii
  0.42 -> 0.68) and smallest where the novel-species genomes have few close
  relatives in the reference set (Lactobacillaceae, nearest ANI 81-86%).
- Novel-versus-known discrimination is unchanged: any-novel precision stays
  at or above 0.989, and the control false-novel count rises by 1-5 reads out
  of about 5,000-6,000.
- Novel Genus recall is unchanged: the rule never removes a Novel Genus
  call.
- The cost is rank confusion for genus-level genomes: more of their reads
  are called Novel Species instead of Ambiguous. In Enterobacteriaceae most
  of the increase comes from Serratia marcescens, whose nearest ANI (79.88%)
  is 0.12 points below the 80% genus line used for the benchmark truth and
  whose reads have median novelty 16.9%; the benchmark's genus truth is
  itself uncertain at that distance.

## Limitations

- Simulated reads with uniform substitution errors and no indels; one read
  length (150 bp); MMseqs2 only (BLAST not benchmarked).
- Reference sets are GTDB species representatives only; strain-rich
  reference sets (many genomes per species) were not tested.
- Truth is assigned from genome-level ANI thresholds (96/80), while
  individual reads vary around the genome-level divergence.
- Three families, 12-16 genomes each.

## Status

The rule improved novel-species recall in all three families, including the
two held-out families, without a measurable loss in novel-versus-known
precision, which meets the held-out-family requirement in CLAUDE.md
("Calibration: diagnostic only"). It is therefore **on by default**.
`--no-clade-novel-species` (or `clade_novelty: {enabled: false}` in a YAML
config) reproduces the previous classification exactly. Re-run the benchmark
with `--compare-clade-rule` when changing the rule's thresholds.

## Reproducing

```bash
export PATH=/path/to/bioconda/bin:$PATH   # fastANI, mmseqs
python scripts/benchmark_leave_one_out.py \
    --genomes genomes/ --out loo_bench/ --threads 4 --compare-clade-rule
```

Genome files must be named `{accession}.fna`; download them with
`mdm download genomes fetch`.
