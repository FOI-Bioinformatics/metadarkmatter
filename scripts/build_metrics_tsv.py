#!/usr/bin/env python
"""
Join ground-truth labels into a classifications file.

Phase 3 helper. Takes:

  * a classifications CSV/Parquet produced by ``score classify`` whose
    ``read_id`` column has been tagged at simulation time with the
    source genome accession (format: ``{accession}__{original_id}``),
  * the ``target_to_label.tsv`` written by ``build_corpus.py``,

and writes a TSV containing the three columns the calibration scripts
require: ``novelty_index``, ``placement_uncertainty``, ``true_category``.

The default ``--include-confidence`` flag also keeps
``confidence_score`` and the predicted ``taxonomic_call`` so the same
file can be fed to ``score evaluate`` for the baseline-vs-calibrated
comparison.

Usage:

    python scripts/build_metrics_tsv.py \\
        --classifications corpora/francisella/baseline.csv \\
        --labels corpora/francisella/target_to_label.tsv \\
        --output corpora/francisella/metrics.tsv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import polars as pl

REQUIRED_INPUT_COLUMNS = ("read_id", "novelty_index", "placement_uncertainty")


def split_accession(read_id: str) -> str:
    """Recover the source genome accession from a tagged read header.

    Tagged headers look like ``GCF_000123456.1__read_42``. Reads that
    do not carry the tag fall back to ``"unknown"`` so the calibration
    script can skip them via ``--min-samples``.
    """
    if "__" not in read_id:
        return "unknown"
    return read_id.split("__", 1)[0]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--classifications",
        required=True,
        type=Path,
        help="score classify output CSV/Parquet with tagged read_ids.",
    )
    parser.add_argument(
        "--labels",
        required=True,
        type=Path,
        help="target_to_label.tsv from build_corpus.py.",
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--include-confidence",
        action="store_true",
        default=True,
        help=(
            "Also carry through confidence_score and taxonomic_call "
            "(needed for score evaluate). On by default."
        ),
    )
    parser.add_argument(
        "--no-confidence",
        dest="include_confidence",
        action="store_false",
        help="Drop confidence_score and taxonomic_call (metrics-only TSV).",
    )
    args = parser.parse_args()

    # Auto-detect format from the extension.
    suffix = args.classifications.suffix.lower()
    if suffix == ".parquet":
        df = pl.read_parquet(args.classifications)
    elif suffix in (".tsv",):
        df = pl.read_csv(args.classifications, separator="\t")
    else:
        df = pl.read_csv(args.classifications)

    missing = set(REQUIRED_INPUT_COLUMNS) - set(df.columns)
    if missing:
        raise SystemExit(
            f"Classification file is missing required columns: {sorted(missing)}"
        )

    labels = pl.read_csv(args.labels, separator="\t")
    if {"target_accession", "true_category"} - set(labels.columns):
        raise SystemExit(
            "Labels file must contain target_accession and true_category columns."
        )

    df = df.with_columns(
        pl.col("read_id")
        .map_elements(split_accession, return_dtype=pl.Utf8)
        .alias("source_accession"),
    )
    joined = df.join(
        labels.select(
            pl.col("target_accession").alias("source_accession"),
            pl.col("true_category"),
        ),
        on="source_accession",
        how="inner",
    )

    if joined.is_empty():
        raise SystemExit(
            "No read_ids could be matched to labelled targets. Verify that "
            "the simulation step tagged read headers with the source "
            "accession (format: {accession}__{original_id})."
        )

    keep_cols = [
        "novelty_index",
        "placement_uncertainty",
        "true_category",
    ]
    if args.include_confidence:
        for extra in ("confidence_score", "taxonomic_call"):
            if extra in joined.columns:
                keep_cols.append(extra)

    out = joined.select(keep_cols)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.write_csv(args.output, separator="\t")
    print(
        f"Wrote {len(out)} rows to {args.output}. "
        f"Categories: {out.group_by('true_category').len().sort('true_category')}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
