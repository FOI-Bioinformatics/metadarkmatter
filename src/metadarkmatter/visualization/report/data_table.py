"""Compact payload builder for the client-side Data tab table.

Instead of emitting one ``<tr>`` per read (roughly 840 bytes each), the
report embeds a column-oriented JSON document: numbers are rounded, and
repeated strings (genome names, categories, ...) are dictionary-encoded as
integer indices. The browser renders only the visible page of rows.

Sampling (``--max-table-rows``) is deterministic: see
:func:`select_table_rows`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import polars as pl

from metadarkmatter.visualization.report.templates import get_cell_class

if TYPE_CHECKING:
    from metadarkmatter.core.metadata import GenomeMetadata

SAMPLING_DESCRIPTION = (
    "balanced across classification categories (rare categories are kept "
    "in full), highest novelty index first within a category"
)


@dataclass(frozen=True)
class ColumnSpec:
    """Definition of one data-table column."""

    key: str
    label: str
    kind: str  # str | genome | species | genus | dict | bool | num
    source: str | None  # DataFrame column (None for derived columns)
    default_visible: bool = False
    optional: bool = False  # only shown when the source holds data
    digits: int = 2
    title: str = ""


# Order is the display order (and the legacy column-index order).
COLUMNS: tuple[ColumnSpec, ...] = (
    ColumnSpec("read_id", "Read ID", "str", "read_id", True),
    ColumnSpec("best_match", "Best Match", "genome", "best_match_genome", True),
    ColumnSpec("species", "Species", "species", None),
    ColumnSpec("genus", "Genus", "genus", None),
    ColumnSpec("identity", "Identity %", "num", "top_hit_identity", True),
    ColumnSpec("novelty", "Novelty", "num", "novelty_index", True),
    ColumnSpec("uncertainty", "Uncertainty", "num", "placement_uncertainty", True),
    ColumnSpec("ambiguous_hits", "Ambig. Hits", "num", "num_ambiguous_hits", digits=0),
    ColumnSpec("identity_gap", "Identity Gap", "num", "identity_gap"),
    ColumnSpec("is_novel", "Is Novel", "bool", "is_novel"),
    ColumnSpec("classification", "Classification", "dict", "taxonomic_call", True),
    ColumnSpec("uncertainty_type", "Unc. Type", "dict", "uncertainty_type", optional=True),
    ColumnSpec("inferred_uncertainty", "Inferred Unc.", "num", "inferred_uncertainty",
               optional=True, digits=1),
    ColumnSpec("discovery_score", "Discovery", "num", "discovery_score",
               optional=True, digits=1),
    ColumnSpec("identity_confidence", "Id. Conf.", "num", "identity_confidence",
               optional=True, digits=1),
    ColumnSpec("placement_confidence", "Pl. Conf.", "num", "placement_confidence",
               optional=True, digits=1),
    ColumnSpec("legacy_call", "Legacy Call", "dict", "legacy_call", optional=True),
    ColumnSpec("confidence_score", "Confidence", "num", "confidence_score",
               optional=True, digits=1),
)


def _has_data(df: pl.DataFrame, column: str) -> bool:
    return column in df.columns and df[column].drop_nulls().len() > 0


def active_columns(df: pl.DataFrame) -> list[ColumnSpec]:
    """Columns to show for ``df`` (optional columns only when they hold data)."""
    cols: list[ColumnSpec] = []
    for spec in COLUMNS:
        if spec.optional and not (spec.source and _has_data(df, spec.source)):
            continue
        cols.append(spec)
    return cols


def select_table_rows(df: pl.DataFrame, max_rows: int | None) -> pl.DataFrame:
    """Return the rows shown in the Data tab table, in original order.

    With no limit (or when ``df`` already fits) every row is returned.
    Otherwise a deterministic *balanced stratified* sample is taken:

    * strata are the ``taxonomic_call`` categories;
    * quotas are allocated by water-filling: categories smaller than the
      fair share are kept in full, and the unused budget is shared among
      the larger categories;
    * inside a category the rows with the highest ``novelty_index`` are
      kept (ties broken by ``read_id`` then original position);
    * the selected rows are returned in their original file order.

    The result depends only on the data, never on hash seeds or threads.
    """
    n = len(df)
    if max_rows is None or n <= max_rows:
        return df
    if max_rows <= 0:
        return df.head(0)

    indexed = df.with_row_index("__row")
    if "taxonomic_call" not in df.columns:
        strata_key = pl.lit("all").alias("__stratum")
    else:
        strata_key = pl.col("taxonomic_call").cast(pl.String).fill_null("").alias("__stratum")
    indexed = indexed.with_columns(strata_key)

    sizes = indexed.group_by("__stratum").len().to_dicts()
    ordered = sorted(sizes, key=lambda r: (r["len"], r["__stratum"]))
    remaining = max_rows
    quotas: dict[str, int] = {}
    for i, rec in enumerate(ordered):
        share = remaining // (len(ordered) - i)
        take = min(rec["len"], share)
        quotas[rec["__stratum"]] = take
        remaining -= take

    sort_cols = ["__novelty_sort", "__rid", "__row"]
    novelty = (
        pl.col("novelty_index").cast(pl.Float64, strict=False)
        if "novelty_index" in df.columns
        else pl.lit(0.0)
    )
    rid = (
        pl.col("read_id").cast(pl.String).fill_null("")
        if "read_id" in df.columns
        else pl.lit("")
    )
    keyed = indexed.with_columns(
        novelty.fill_nan(None).fill_null(float("-inf")).alias("__novelty_sort"),
        rid.alias("__rid"),
    )
    picked: list[pl.DataFrame] = []
    for stratum in sorted(quotas):
        q = quotas[stratum]
        if q <= 0:
            continue
        sub = keyed.filter(pl.col("__stratum") == stratum).sort(
            sort_cols, descending=[True, False, False]
        )
        picked.append(sub.head(q).select("__row"))
    rows = pl.concat(picked).sort("__row")["__row"].to_list() if picked else []
    return df[rows]


def _round_list(series: pl.Series, digits: int) -> list[float | int | None]:
    values = series.cast(pl.Float64, strict=False).fill_nan(None).to_list()
    out: list[float | int | None] = []
    for v in values:
        if v is None:
            out.append(None)
        elif digits == 0:
            out.append(round(v))
        else:
            out.append(round(v, digits) + 0.0)  # +0.0 turns -0.0 into 0.0
    return out


def _dict_encode(values: list[Any]) -> tuple[list[str], list[int]]:
    """Dictionary-encode strings; dictionary is sorted for determinism."""
    cleaned = ["" if v is None else str(v) for v in values]
    dictionary = sorted(set(cleaned))
    lookup = {v: i for i, v in enumerate(dictionary)}
    return dictionary, [lookup[v] for v in cleaned]


def _genome_taxa(
    genomes: list[str],
    df: pl.DataFrame,
    genome_metadata: GenomeMetadata | None,
) -> tuple[list[str], list[str]]:
    """Species and genus for each genome in ``genomes`` (same order)."""
    species: list[str] = []
    genus: list[str] = []
    df_map: dict[str, tuple[str, str]] = {}
    # The table's own species/genus columns are the fallback for genomes the
    # metadata does not cover (or has no value for).
    if "species" in df.columns and "best_match_genome" in df.columns:
        cols = ["best_match_genome", "species"] + (["genus"] if "genus" in df.columns else [])
        pairs = df.select(cols).unique().sort(cols, nulls_last=True)
        for row in pairs.iter_rows(named=True):
            g = row["best_match_genome"]
            if g not in df_map and row["species"] is not None:
                df_map[g] = (str(row["species"]), str(row.get("genus") or ""))
    for g in genomes:
        sp, ge = df_map.get(g, ("", ""))
        if genome_metadata is not None:
            sp = genome_metadata.get_species(g) or sp
            ge = genome_metadata.get_genus(g) or ge
        species.append(sp)
        genus.append(ge)
    return species, genus


def build_table_payload(
    df: pl.DataFrame,
    genome_metadata: GenomeMetadata | None,
    *,
    page_size: int,
    total_reads: int | None = None,
) -> dict[str, Any]:
    """Build the JSON-serialisable payload for the client-side table.

    Args:
        df: The (possibly sampled) rows to embed.
        genome_metadata: Optional metadata for species/genus lookups.
        page_size: Rows per page.
        total_reads: Rows in the full classification file (for the notice).
    """
    columns = active_columns(df)
    n = len(df)

    cols: dict[str, list[Any]] = {}
    dicts: dict[str, list[str]] = {}
    col_defs: list[dict[str, Any]] = []

    genomes_dict: list[str] = []
    for spec in columns:
        entry: dict[str, Any] = {
            "key": spec.key,
            "label": spec.label,
            "kind": spec.kind,
            "visible": spec.default_visible,
        }
        if spec.kind == "num":
            entry["digits"] = spec.digits
        col_defs.append(entry)

        if spec.kind == "str":
            values = df[spec.source].cast(pl.String).fill_null("").to_list() if spec.source in df.columns else [""] * n
            cols[spec.key] = values
        elif spec.kind == "genome":
            values = df[spec.source].cast(pl.String).fill_null("").to_list() if spec.source in df.columns else [""] * n
            genomes_dict, cols["genome"] = _dict_encode(values)
        elif spec.kind == "dict":
            if spec.source in df.columns:
                dicts[spec.key], cols[spec.key] = _dict_encode(df[spec.source].to_list())
            else:
                dicts[spec.key], cols[spec.key] = [""], [0] * n
        elif spec.kind == "bool":
            if spec.source in df.columns and df[spec.source].null_count() < n:
                raw = df[spec.source].to_list()
            else:
                raw = [None] * n
            calls = df["taxonomic_call"].to_list() if "taxonomic_call" in df.columns else [None] * n
            cols[spec.key] = [
                1 if (r if r is not None else c in ("Novel Species", "Novel Genus")) else 0
                for r, c in zip(raw, calls, strict=True)
            ]
        elif spec.kind == "num":
            if spec.source in df.columns:
                series = df[spec.source]
                if spec.key == "identity_gap":
                    series = series.cast(pl.Float64, strict=False).abs()
                cols[spec.key] = _round_list(series, spec.digits)
            else:
                cols[spec.key] = [None] * n

    species_dict, genus_dict = _genome_taxa(genomes_dict, df, genome_metadata)

    classes = {c: get_cell_class(c) for c in dicts.get("classification", [])}
    return {
        "n": n,
        "total": total_reads if total_reads is not None else n,
        "pageSize": page_size,
        "columns": col_defs,
        "genomes": genomes_dict,
        "species": species_dict,
        "genus": genus_dict,
        "dicts": dicts,
        "cellClass": {k: v for k, v in classes.items() if v},
        "cols": cols,
    }
