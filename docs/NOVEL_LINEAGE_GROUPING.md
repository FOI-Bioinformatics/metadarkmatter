# Novel lineage grouping: multi-holdout validation

This document records the validation of the genus grouping of novel reads
(`NovelDiversityAnalyzer(grouping=...)`, `report generate --cluster-grouping`,
`report summarize --grouping`). The method is described in
[METHODS.md](METHODS.md#41-clustering-algorithm). Grouping only changes how
the report and `report summarize` group reads that are already called Novel
Species or Novel Genus; it does not change any read's classification.

## Problem

The Novel Diversity tab is meant to answer "how many novel organisms are in
this sample, and which reads belong to each?". The previous (legacy) method
grouped novel reads by (ANI neighbourhood of the nearest reference, 5%
novelty band, call), with neighbourhoods formed by single linkage at 80% ANI
over the whole reference set. Two failure modes were observed:

- **Chaining in dense families.** With five Enterobacteriaceae genomes held
  out at once, 80% single linkage joined the Escherichia, Salmonella,
  Klebsiella, Enterobacter and Citrobacter references into one
  neighbourhood. The five novel organisms fell into four groups, and on
  average 26% of a group's reads came from its majority organism.
- **Splitting.** Reads of one novel organism spread over its roughly
  equidistant relatives. Reads of a held-out *Klebsiella aerogenes*
  (86% ANI to its nearest remaining reference) were split about evenly over
  three *Klebsiella* references, and the band and call split separated them
  further.

## Method compared

| Method | Group key |
|---|---|
| legacy | (80% ANI neighbourhood over all references, 5% novelty band, call) |
| genus (new default) | GTDB genus of the nearest reference, split within the genus by single linkage at ANI >= 80% |
| reference (fallback without metadata) | nearest reference (species representative) |

Genus and reference groups are not split by band or call; a group's call is
the majority call of its reads.

An ANI-only genus boundary (single linkage at a fixed ANI between
references, without metadata) was also tested and rejected: agreement with
GTDB genera (adjusted Rand index) in Enterobacteriaceae was 0.10, 0.49 and
0.73 at 80, 82 and 85% ANI, while in Francisellaceae 80% gave 1.0. No single
threshold works across families.

## Benchmark

`scripts/benchmark_lineage_grouping.py` uses the output directory of
`scripts/benchmark_leave_one_out.py` (simulated 150 bp reads, uniform 0.5%
substitution errors, fastANI, MMseqs2; see
[CLADE_NOVELTY_VALIDATION.md](CLADE_NOVELTY_VALIDATION.md)). Several genomes
are held out **at once**, so that the sample contains several novel
organisms. Their hits and ANI rows/columns are removed, all reads are
classified with the default `ScoringConfig`, and novel reads are grouped with
each method. The true source of each read comes from its read ID.

Metrics over reads in groups of at least 3 reads:

- **purity**: fraction of grouped reads that belong to their group's majority source;
- **pair precision / recall**: of read pairs in one group, the fraction from one source / of read pairs from one source, the fraction in one group;
- **completeness**: per held-out genome, the fraction of its novel reads that are in its largest group.

Held-out sets:

- **Enterobacteriaceae** (design set): *E. albertii* GCF_000759775.1 (nearest
  remaining ANI 90.4%), *K. aerogenes* GCF_000215745.1 (86.4%), *S. bongori*
  GCF_000252995.1 (90.1%), *E. asburiae* GCF_001521715.1 (89.1%), *Serratia*
  GCF_017299535.1 (79.9%).
- **Lactobacillaceae** (held out from the design): *Lacticaseibacillus
  paracasei* GCF_000829035.1 (81.2%), *L. acidophilus* GCF_034298135.1
  (81.9%), *L. helveticus* GCF_000160855.1 (82.4%), *L. gasseri*
  GCF_000014425.1 (85.5%), *Limosilactobacillus fermentum* GCF_013394085.1
  (novel genus).
- **Francisellaceae, full pipeline**: a simulated paired-end sample from
  Francisellaceae genomes run through `mdm mmseqs2 search`, `mdm ani compute
  --backend fastani` and `mdm score classify` with metadata, where
  GCF_007923265.1 (*F. salimarina*) and GCF_003290445.1 are absent from the
  reference database.

## Results

| Set | Method | Groups | Purity | Pair precision | Pair recall | Completeness |
|---|---|---|---|---|---|---|
| Enterobacteriaceae | legacy | 4 | 0.264 | 0.215 | 0.825 | 0.94 / 0.87 / 0.93 / 0.92 / 0.61 |
| | **genus** | 6 | **0.809** | 0.704 | 0.750 | 0.95 / 0.87 / 0.85 / 0.85 / 0.54 |
| | reference | 10 | 0.809 | 0.757 | 0.553 | 0.63 / 0.30 / 0.85 / 0.85 / 0.23 |
| Lactobacillaceae | legacy | 13 | 0.760 | 0.687 | 0.538 | 0.53 / 0.64 / 0.55 / 0.82 / 0.31 |
| | **genus** | 6 | 0.746 | 0.637 | **0.844** | 1.00 / 0.90 / 0.88 / 0.94 / 0.49 |
| | reference | 7 | 0.746 | 0.615 | 0.759 | 0.51 / 0.90 / 0.88 / 0.94 / 0.49 |
| Francisellaceae (pipeline) | legacy | 6 | 0.833 | 0.716 | 0.775 | 0.94 / 0.37 |
| | **genus** | 3 | 0.771 | 0.622 | **0.901** | 0.99 / 0.68 |
| | reference | 10 | 0.822 | 0.724 | 0.224 | 0.39 / 0.15 |

Observations:

- In the dense family, genus grouping raises purity from 0.26 to 0.81 and
  pair precision from 0.22 to 0.70; the legacy method's high recall there
  comes from putting nearly everything into one group.
- In the two other sets, genus grouping trades a small loss of purity
  (0.01-0.06) for a large gain in recall (+0.13 to +0.31): reads of one
  organism are kept together (completeness of the main group up to 0.99-1.00
  where the legacy method had 0.37-0.64).
- The within-genus 80% split matters for the held-out family: without it,
  the three held-out *Lactobacillus* species and their relatives formed one
  group (purity 0.569, pair precision 0.390). With it, *L. gasseri*
  (nearest reference 85.5% ANI) is separated from *L. acidophilus* and
  *L. helveticus*.
- Nearest-reference grouping (the fallback without metadata) has similar
  purity but splits organisms (*K. aerogenes* completeness 0.30, Francisella
  0.39/0.15), so `report` warns when it is used.

## Limitations

- Novel species of one genus whose reads hit the same reference are pooled:
  *L. acidophilus* and *L. helveticus* reads both hit the same
  *Lactobacillus* reference and share one group.
- A novel genus whose reads hit an existing genus is pooled with it
  (*Serratia* reads that hit *Klebsiella*; completeness 0.54).
- Requires GTDB genus names in the metadata (`mdm download genomes list`
  always writes them); otherwise the nearest-reference fallback applies.
- Simulated reads only (5-genome held-out sets in two families, and one
  2-genome set run through the full pipeline). Groups are candidates for follow-up, not counts of species.

## Reproducing

```bash
python scripts/benchmark_leave_one_out.py --genomes genomes/ --out loo_bench/ --threads 4
python scripts/benchmark_lineage_grouping.py --bench loo_bench/ \
    --metadata genome_metadata.tsv \
    --held GCF_000759775.1,GCF_000215745.1,GCF_000252995.1,GCF_001521715.1,GCF_017299535.1
```

`--cluster-grouping legacy` (report) and `--grouping legacy` (summarize)
reproduce the previous groups exactly; this was checked against the code
before this change on the Francisellaceae example.
