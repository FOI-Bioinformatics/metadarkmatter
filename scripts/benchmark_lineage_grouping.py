#!/usr/bin/env python
"""
Multi-holdout benchmark for grouping novel reads into novel lineage groups.

Several genomes are held out at once from a leave-one-out benchmark
directory (``scripts/benchmark_leave_one_out.py --out DIR``, which holds
``alignment.tsv.gz`` and ``ani_matrix.csv`` with read IDs ``r{N}_{accession}``).
Their hits and ANI rows/columns are removed, all reads are classified with
the default ``ScoringConfig``, and the novel reads are grouped with each
``NovelDiversityAnalyzer`` grouping method. Each read's true source genome
is taken from its read ID.

Metrics (over reads in groups of >= 3 reads):

    purity          fraction of grouped reads that belong to their group's
                    majority source genome
    pair precision  of read pairs in one group, the fraction from one source
    pair recall     of read pairs from one source, the fraction in one group
    completeness    per held-out genome, the fraction of its novel reads that
                    are in its largest group

Usage:

    python scripts/benchmark_lineage_grouping.py --bench loo_bench/ \
        --metadata genome_metadata.tsv \
        --held GCF_000759775.1,GCF_000215745.1,GCF_000252995.1
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import polars as pl

from metadarkmatter.core.classification.ani_matrix import ANIMatrix
from metadarkmatter.core.classification.classifiers.vectorized import VectorizedClassifier
from metadarkmatter.core.metadata import GenomeMetadata
from metadarkmatter.core.novel_diversity.clustering import NovelDiversityAnalyzer
from metadarkmatter.core.parsers import StreamingBlastParser
from metadarkmatter.models.config import ScoringConfig

NOVEL = ["Novel Species", "Novel Genus"]


def score_groups(
    groups: dict[str, list[str]], src_of: dict[str, str], novel_by_src: dict[str, list[str]],
    held: list[str],
) -> dict[str, object]:
    """Purity, pair precision/recall and per-genome completeness."""
    label = {r: g for g, rs in groups.items() for r in rs}

    def pairs(n: int) -> float:
        return n * (n - 1) / 2

    cont = Counter((label[r], src_of[r]) for r in label)
    tp = sum(pairs(n) for n in cont.values())
    same_group = sum(pairs(n) for n in Counter(label.values()).values())
    same_src = sum(pairs(n) for n in Counter(src_of[r] for r in label).values())
    purity = sum(max(Counter(src_of[r] for r in rs).values()) for rs in groups.values())
    completeness = {}
    for h in held:
        counts = Counter(label[r] for r in novel_by_src.get(h, []) if r in label)
        n = len(novel_by_src.get(h, []))
        completeness[h] = round(max(counts.values()) / n, 3) if counts and n else 0.0
    return {
        "groups": len(groups),
        "reads_grouped": len(label),
        "purity": round(purity / max(1, len(label)), 3),
        "pair_precision": round(tp / same_group, 3) if same_group else 1.0,
        "pair_recall": round(tp / same_src, 3) if same_src else 1.0,
        "completeness": completeness,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bench", type=Path, required=True, help="benchmark_leave_one_out.py output directory")
    ap.add_argument("--metadata", type=Path, required=True, help="genome_metadata.tsv (genus column)")
    ap.add_argument("--held", required=True, help="comma-separated accessions to hold out together")
    ap.add_argument("--json", type=Path, help="write metrics as JSON")
    args = ap.parse_args()
    held = [h.strip() for h in args.held.split(",") if h.strip()]

    ani_df = pl.read_csv(args.bench / "ani_matrix.csv")
    names = [c for c in ani_df.columns if c != "genome"]
    full = {r["genome"]: {g: float(r[g]) for g in names} for r in ani_df.iter_rows(named=True)}
    missing = [h for h in held if h not in full]
    if missing:
        ap.error(f"held-out accessions not in the ANI matrix: {missing}")
    ani = ANIMatrix({a: {b: v for b, v in row.items() if b not in held}
                     for a, row in full.items() if a not in held})
    aln = (
        StreamingBlastParser(args.bench / "alignment.tsv.gz").parse_eager()
        .with_columns(pl.col("sseqid").str.extract(r"^([^|]+)\|").alias("genome_name"))
        .filter(~pl.col("genome_name").is_in(held))
    )
    cls = VectorizedClassifier(ani, config=ScoringConfig()).classify_file(aln).with_columns(
        pl.col("read_id").str.extract(r"_(GC[AF]_\d+\.\d+)").alias("src")
    )
    md = GenomeMetadata.from_file(args.metadata)
    novel = cls.filter(pl.col("taxonomic_call").is_in(NOVEL))
    src_of = dict(zip(cls["read_id"], cls["src"], strict=True))
    novel_by_src: dict[str, list[str]] = {}
    for rid, src in zip(novel["read_id"], novel["src"], strict=True):
        novel_by_src.setdefault(src, []).append(rid)

    nearest = {h: round(max(v for g, v in full[h].items() if g not in held and g != h), 2) for h in held}
    print(f"held out (nearest remaining ANI): {nearest}")
    print(f"novel reads: {novel.height}; from genomes still in the database: "
          f"{novel.filter(~pl.col('src').is_in(held)).height}")
    results: dict[str, object] = {"held": nearest}
    for name, kwargs in (
        ("legacy", {"metadata": md, "grouping": "legacy"}),
        ("genus", {"metadata": md, "grouping": "genus"}),
        ("reference", {"metadata": md, "grouping": "reference"}),
    ):
        an = NovelDiversityAnalyzer(cls, ani_matrix=ani, **kwargs)
        groups = {c.cluster_id: list(c.read_ids) for c in an.cluster_novel_reads()}
        m = score_groups(groups, src_of, novel_by_src, held)
        results[name] = m
        print(f"  {name:10s} groups={m['groups']:3d} purity={m['purity']:.3f} "
              f"pairP={m['pair_precision']:.3f} pairR={m['pair_recall']:.3f} "
              f"completeness={list(m['completeness'].values())}")  # type: ignore[union-attr]
    if args.json:
        args.json.write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
