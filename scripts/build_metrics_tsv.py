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

# Per-read ANI label thresholds. These mirror build_corpus.py so the
# two scripts agree on what counts as Known/Novel Species/Genus.
ANI_KNOWN_SPECIES = 96.0
ANI_NOVEL_SPECIES = 80.0
ANI_NOVEL_GENUS = 75.0


def split_accession(read_id: str) -> str:
    """Recover the source genome accession from a tagged read header.

    Tagged headers look like ``GCF_000123456.1__read_42``. Reads that
    do not carry the tag fall back to ``"unknown"`` so the calibration
    script can skip them via ``--min-samples``.
    """
    if "__" not in read_id:
        return "unknown"
    return read_id.split("__", 1)[0]


def _category_from_ani(ani: float) -> str:
    """Map a single ANI value to a Bayesian category.

    Used by ``per_read`` label mode where each read's truth is the
    ANI between the source target and the reference the classifier
    picked for that read. Species Boundary and Ambiguous are not
    produced from a single ANI value alone; they require the
    classifier's own Stage 2 reasoning to identify and would need
    second-best context that this mode does not have.
    """
    if ani >= ANI_KNOWN_SPECIES:
        return "Known Species"
    if ani >= ANI_NOVEL_SPECIES:
        return "Novel Species"
    if ani >= ANI_NOVEL_GENUS:
        return "Novel Genus"
    return "Unclassified"


def _load_ani_matrix(path: Path) -> dict[str, dict[str, float]]:
    """Load a metadarkmatter ANI CSV into a nested dict for O(1) lookup."""
    df = pl.read_csv(path)
    rows = df.to_dicts()
    ref_genomes = df.columns[1:]
    matrix: dict[str, dict[str, float]] = {}
    for row in rows:
        # The first column is conventionally named 'genome' but is just
        # the row label; fall back to the first column name in case the
        # convention changes.
        src = row[df.columns[0]]
        matrix[src] = {g: float(row[g]) for g in ref_genomes}
    return matrix


def _apply_per_read_labels(
    df: pl.DataFrame,
    ani_matrix: dict[str, dict[str, float]],
) -> pl.DataFrame:
    """Compute a per-read true_category column.

    For every row we look up ``ANI(source_accession, best_match_genome)``
    and map the value through ``_category_from_ani``. Reads whose
    source-target pair is missing from the matrix (e.g. an off-target
    hit outside the family) get the conservative ``Unclassified``
    label.
    """

    def label(source: str | None, predicted: str | None) -> str:
        if source is None or predicted is None:
            return "Unclassified"
        ani = ani_matrix.get(source, {}).get(predicted, 0.0)
        return _category_from_ani(ani)

    return df.with_columns(
        pl.struct(["source_accession", "best_match_genome"])
        .map_elements(
            lambda row: label(row.get("source_accession"), row.get("best_match_genome")),
            return_dtype=pl.Utf8,
        )
        .alias("true_category")
    )


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
    parser.add_argument(
        "--label-mode",
        choices=("per_genome", "per_read"),
        default="per_genome",
        help=(
            "How to assign truth labels. 'per_genome' (default) uses the "
            "label from --labels for every read produced by a given source "
            "target. 'per_read' looks up the actual ANI between the source "
            "target and the genome the classifier picked for each read, "
            "producing a finer-grained truth that accounts for reads "
            "hitting conserved regions across species. Requires --ani-matrix."
        ),
    )
    parser.add_argument(
        "--ani-matrix",
        type=Path,
        default=None,
        help=(
            "Path to the full-family ANI matrix (e.g. ani_full.csv from "
            "build_corpus.py). Required when --label-mode=per_read."
        ),
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

    # Restrict to reads from labelled targets so the eval is well-defined.
    target_set = set(labels["target_accession"].to_list())
    df = df.filter(pl.col("source_accession").is_in(target_set))

    if args.label_mode == "per_read":
        if args.ani_matrix is None:
            raise SystemExit(
                "--label-mode=per_read requires --ani-matrix pointing at "
                "the full-family ANI CSV from build_corpus.py."
            )
        if "best_match_genome" not in df.columns:
            raise SystemExit(
                "per_read mode requires a 'best_match_genome' column in "
                "the classifications file."
            )
        ani_matrix = _load_ani_matrix(args.ani_matrix)
        joined = _apply_per_read_labels(df, ani_matrix)
    else:
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
