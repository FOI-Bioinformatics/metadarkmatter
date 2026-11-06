"""Pure builders for the call-comparison and species-breakdown tables.

These functions take data frames and return plain data / HTML strings so
they can be unit-tested without constructing a full report.
"""

from __future__ import annotations

import html
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import polars as pl

from metadarkmatter.visualization.plots.base import TAXONOMY_COLORS

if TYPE_CHECKING:
    from metadarkmatter.core.metadata import GenomeMetadata

NO_CALL = "(none)"
OFF_TARGET = "Off-target"
NOVEL_CALLS: tuple[str, ...] = ("Novel Species", "Novel Genus")

# Maximum species/genome rows rendered in the HTML table (CSV has all rows).
MAX_BREAKDOWN_TABLE_ROWS = 500
# Maximum genomes listed per species in the "best-match genomes" cell.
MAX_LISTED_GENOMES = 3
# Largest similarity matrix embedded for CSV download (n x n cells).
MAX_EXPORT_MATRIX_GENOMES = 500


def order_categories(categories: set[str] | list[str]) -> list[str]:
    """Canonical category order (palette order), then any extras sorted."""
    present = set(categories)
    ordered = [c for c in TAXONOMY_COLORS if c in present]
    ordered += sorted(present - set(TAXONOMY_COLORS))
    return ordered


def _esc(value: object) -> str:
    return html.escape("" if value is None else str(value))


def _swatch(category: str) -> str:
    color = TAXONOMY_COLORS.get(category, "#999999")
    return (
        f'<span class="legend-dot" style="background:{color};display:inline-block;'
        f'width:0.7em;height:0.7em;border-radius:50%;margin-right:0.35em;'
        f'box-shadow:inset 0 0 0 1px rgba(0,0,0,0.3)" aria-hidden="true"></span>'
    )


# ---------------------------------------------------------------------------
# taxonomic_call x legacy_call
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CallComparison:
    """Confusion table of the Bayesian call against the legacy cascade."""

    rows: list[str]  # taxonomic_call categories (row labels)
    cols: list[str]  # legacy_call categories (column labels)
    counts: dict[tuple[str, str], int]
    total: int
    agree: int

    @property
    def agreement_pct(self) -> float:
        return self.agree / self.total * 100 if self.total else 0.0

    def row_total(self, row: str) -> int:
        return sum(self.counts.get((row, c), 0) for c in self.cols)

    def col_total(self, col: str) -> int:
        return sum(self.counts.get((r, col), 0) for r in self.rows)


def compute_call_comparison(df: pl.DataFrame) -> CallComparison | None:
    """Cross-tabulate ``taxonomic_call`` (rows) against ``legacy_call`` (columns).

    Returns ``None`` when either column is missing or the frame is empty.
    Null calls are shown as ``(none)``.
    """
    if "taxonomic_call" not in df.columns or "legacy_call" not in df.columns or len(df) == 0:
        return None
    pairs = (
        df.select(
            pl.col("taxonomic_call").cast(pl.String).fill_null(NO_CALL).alias("t"),
            pl.col("legacy_call").cast(pl.String).fill_null(NO_CALL).alias("l"),
        )
        .group_by(["t", "l"])
        .len()
        .sort(["t", "l"])
    )
    counts: dict[tuple[str, str], int] = {}
    for rec in pairs.iter_rows(named=True):
        counts[(rec["t"], rec["l"])] = int(rec["len"])
    rows = order_categories({t for t, _ in counts})
    cols = order_categories({lcl for _, lcl in counts})
    total = sum(counts.values())
    agree = sum(n for (t, lcl), n in counts.items() if t == lcl)
    return CallComparison(rows, cols, counts, total, agree)


CALL_COMPARISON_INTRO = (
    "<p class=\"mdm-note\"><strong>taxonomic_call</strong> is the primary call: the most "
    "probable category under a Bayesian model of each read's novelty and placement "
    "uncertainty, followed by two discrete refinement rules. <strong>legacy_call</strong> "
    "is the earlier cascade of fixed thresholds on the same two values (for example "
    "Known Species below 4% novelty and 1.5% uncertainty). The Bayesian decision "
    "boundaries do not coincide with the legacy thresholds (with default settings a "
    "multi-hit read changes from Known Species to Novel Species at roughly 6% novelty "
    "rather than 4%), so disagreement is expected near category boundaries and is not "
    "by itself evidence of an error. Cells on the diagonal are agreements.</p>"
)


def build_call_comparison_html(cmp: CallComparison) -> str:
    """HTML confusion table (counts with row/column totals and agreement)."""
    head = "".join(f'<th scope="col">{_swatch(c)}{_esc(c)}</th>' for c in cmp.cols)
    body: list[str] = []
    for r in cmp.rows:
        cells = []
        for c in cmp.cols:
            n = cmp.counts.get((r, c), 0)
            cls = []
            if r == c and n:
                cls.append("diag")
            if n == 0:
                cls.append("zero")
            class_attr = f' class="{" ".join(cls)}"' if cls else ""
            cells.append(f"<td{class_attr}>{n:,}</td>")
        body.append(
            f'<tr><th scope="row">{_swatch(r)}{_esc(r)}</th>{"".join(cells)}'
            f'<td class="total">{cmp.row_total(r):,}</td></tr>'
        )
    foot_cells = "".join(f'<td class="total">{cmp.col_total(c):,}</td>' for c in cmp.cols)
    return (
        CALL_COMPARISON_INTRO
        + f'<p class="mdm-note" id="call-agreement"><strong>Agreement:</strong> '
        f"{cmp.agree:,} of {cmp.total:,} reads ({cmp.agreement_pct:.1f}%) have the same "
        "taxonomic_call and legacy_call.</p>"
        + '<div class="mdm-table-wrap"><table class="mdm-table" id="call-comparison-table">'
        "<caption>Reads by taxonomic_call (rows) and legacy_call (columns)</caption>"
        f'<thead><tr><th scope="col">taxonomic_call \\ legacy_call</th>{head}'
        '<th scope="col" class="total">Total</th></tr></thead>'
        f'<tbody>{"".join(body)}</tbody>'
        f'<tfoot><tr><th scope="row">Total</th>{foot_cells}'
        f'<td class="total">{cmp.total:,}</td></tr></tfoot></table></div>'
    )


def call_comparison_csv(cmp: CallComparison) -> tuple[list[str], list[list[Any]]]:
    """Header and rows for the call-comparison CSV download."""
    header = ["taxonomic_call", *cmp.cols, "total"]
    rows: list[list[Any]] = [
        [r, *[cmp.counts.get((r, c), 0) for c in cmp.cols], cmp.row_total(r)]
        for r in cmp.rows
    ]
    rows.append(["total", *[cmp.col_total(c) for c in cmp.cols], cmp.total])
    return header, rows


# ---------------------------------------------------------------------------
# Per-species / per-genome breakdown
# ---------------------------------------------------------------------------


def _genome_labels(
    df: pl.DataFrame, genome_metadata: GenomeMetadata | None
) -> dict[str, str]:
    """Genome accession -> species label (accession when no species known)."""
    genomes = sorted(set(df["best_match_genome"].drop_nulls().to_list()))
    labels: dict[str, str] = {}
    if "species" in df.columns:
        pairs = (
            df.select("best_match_genome", "species")
            .drop_nulls()
            .unique()
            .sort(["best_match_genome", "species"])
        )
        for rec in pairs.iter_rows(named=True):
            labels.setdefault(rec["best_match_genome"], str(rec["species"]))
    for g in genomes:
        if g in labels and labels[g] not in ("", "Unknown"):
            continue
        sp = genome_metadata.get_species(g) if genome_metadata is not None else None
        labels[g] = sp if sp and sp != "Unknown" else g
    return labels


def compute_species_breakdown(
    df: pl.DataFrame, genome_metadata: GenomeMetadata | None
) -> tuple[pl.DataFrame, list[str], bool]:
    """Aggregate on-target reads per species (or per genome without metadata).

    Returns:
        ``(table, categories, uses_species)``. ``table`` has columns
        ``label, genomes, n_genomes, reads, <one per category>, novel_reads,
        novel_pct, mean_identity`` sorted by reads (descending) then label.
        ``uses_species`` is False when labels are genome accessions.
    """
    needed = {"best_match_genome", "taxonomic_call"}
    if not needed <= set(df.columns) or len(df) == 0:
        return pl.DataFrame(), [], False
    base = df.filter(pl.col("taxonomic_call") != OFF_TARGET)
    if len(base) == 0:
        return pl.DataFrame(), [], False

    labels = _genome_labels(base, genome_metadata)
    uses_species = any(lbl != g for g, lbl in labels.items())
    base = base.with_columns(
        pl.col("best_match_genome")
        .replace_strict(labels, default=None, return_dtype=pl.String)
        .fill_null(pl.col("best_match_genome").cast(pl.String).fill_null("(unknown)"))
        .alias("__label")
    )
    categories = order_categories(set(base["taxonomic_call"].drop_nulls().to_list()))
    ident = (
        pl.col("top_hit_identity").cast(pl.Float64, strict=False).fill_nan(None).mean()
        if "top_hit_identity" in base.columns
        else pl.lit(None, dtype=pl.Float64)
    )
    aggs = [pl.len().alias("reads")]
    aggs += [(pl.col("taxonomic_call") == c).sum().alias(c) for c in categories]
    aggs.append(ident.alias("mean_identity"))
    table = base.group_by("__label").agg(aggs).rename({"__label": "label"})

    # best-match genomes per label (most reads first, ties by accession)
    per_genome = (
        base.group_by(["__label", "best_match_genome"])
        .len()
        .sort(["__label", "len", "best_match_genome"], descending=[False, True, False])
    )
    listing: dict[str, list[str]] = {}
    counts: dict[str, int] = {}
    for rec in per_genome.iter_rows(named=True):
        lab = rec["__label"]
        counts[lab] = counts.get(lab, 0) + 1
        bucket = listing.setdefault(lab, [])
        if len(bucket) < MAX_LISTED_GENOMES:
            bucket.append(f"{rec['best_match_genome']} ({rec['len']:,})")
    novel_expr = sum(
        (pl.col(c) for c in NOVEL_CALLS if c in categories), start=pl.lit(0)
    )
    table = (
        table.with_columns(
            novel_expr.alias("novel_reads"),
            pl.col("label").replace_strict(counts, default=0, return_dtype=pl.Int64).alias("n_genomes"),
            pl.col("label")
            .replace_strict(
                {k: ", ".join(v) for k, v in listing.items()}, default="", return_dtype=pl.String
            )
            .alias("genomes"),
        )
        .with_columns((pl.col("novel_reads") / pl.col("reads") * 100).alias("novel_pct"))
        .sort(["reads", "label"], descending=[True, False])
    )
    return table, categories, uses_species


def build_species_table_html(
    table: pl.DataFrame, categories: list[str], *, uses_species: bool
) -> str:
    """Sortable HTML table (top ``MAX_BREAKDOWN_TABLE_ROWS`` rows)."""
    first = "Species" if uses_species else "Genome"
    shown = table.head(MAX_BREAKDOWN_TABLE_ROWS)

    def th(label: str, kind: str, extra: str = "", sort: str = "none") -> str:
        return (
            f'<th scope="col" class="sortable {extra}" data-type="{kind}" aria-sort="{sort}">'
            f'<button type="button" class="sort-btn">{label} '
            '<span class="sort-icon" aria-hidden="true"></span></button></th>'
        )

    head = (
        th(first, "str")
        + th("Best-match genome(s)", "str")
        + th("Reads", "num", sort="descending")
        + "".join(th(f"{_swatch(c)}{_esc(c)}", "num") for c in categories)
        + th("Novel %", "num")
        + th("Mean identity %", "num")
    )
    body: list[str] = []
    for rec in shown.iter_rows(named=True):
        ident = rec["mean_identity"]
        cells = [
            f'<td data-v="{_esc(rec["label"]).lower()}">{_esc(rec["label"])}</td>',
            f'<td data-v="{_esc(rec["genomes"]).lower()}">{_esc(rec["genomes"])}'
            + (
                f" +{rec['n_genomes'] - MAX_LISTED_GENOMES} more"
                if rec["n_genomes"] > MAX_LISTED_GENOMES
                else ""
            )
            + "</td>",
            f'<td data-v="{rec["reads"]}">{rec["reads"]:,}</td>',
        ]
        for c in categories:
            n = rec[c]
            cells.append(f'<td data-v="{n}"' + (' class="zero"' if not n else "") + f">{n:,}</td>")
        cells.append(f'<td data-v="{rec["novel_pct"]:.4f}">{rec["novel_pct"]:.1f}</td>')
        cells.append(
            f'<td data-v="{ident:.4f}">{ident:.2f}</td>' if ident is not None
            else '<td data-v="">-</td>'
        )
        body.append("<tr>" + "".join(cells) + "</tr>")
    note = ""
    if len(table) > len(shown):
        note = (
            f'<p class="mdm-note">Showing the {len(shown):,} rows with most reads of '
            f"{len(table):,}; the CSV download contains all rows.</p>"
        )
    return (
        '<div class="mdm-table-wrap"><table class="mdm-table" id="species-breakdown-table" '
        'data-sortable="1">'
        f"<caption>On-target reads per {first.lower()} by classification</caption>"
        f'<thead><tr>{head}</tr></thead><tbody>{"".join(body)}</tbody></table></div>{note}'
    )


def species_breakdown_csv(
    table: pl.DataFrame, categories: list[str], *, uses_species: bool
) -> tuple[list[str], list[list[Any]]]:
    """Header and rows for the species-breakdown CSV (all rows)."""
    header = [
        "species" if uses_species else "genome", "best_match_genomes", "n_genomes", "reads",
        *categories, "novel_reads", "novel_pct", "mean_identity",
    ]
    rows: list[list[Any]] = []
    for rec in table.iter_rows(named=True):
        ident = rec["mean_identity"]
        rows.append([
            rec["label"], rec["genomes"], rec["n_genomes"], rec["reads"],
            *[rec[c] for c in categories], rec["novel_reads"],
            round(rec["novel_pct"], 2), None if ident is None else round(ident, 3),
        ])
    return header, rows


# ---------------------------------------------------------------------------
# Matrices and cluster tables (CSV exports)
# ---------------------------------------------------------------------------


def matrix_csv(matrix: pl.DataFrame) -> tuple[list[str], list[list[Any]]] | None:
    """Header and rows for a similarity-matrix CSV download.

    Returns ``None`` when the matrix is empty or larger than
    ``MAX_EXPORT_MATRIX_GENOMES`` (to keep the report small).
    """
    if len(matrix) == 0 or len(matrix) > MAX_EXPORT_MATRIX_GENOMES:
        return None
    if "genome" in matrix.columns:
        names = [str(g) for g in matrix["genome"].to_list()]
        value_cols = [c for c in matrix.columns if c != "genome"]
        values = matrix.select(value_cols)
    else:
        value_cols = list(matrix.columns)
        names = list(value_cols)
        values = matrix
    rows: list[list[Any]] = []
    for name, rec in zip(names, values.iter_rows(), strict=True):
        rows.append([
            name,
            *[None if v is None or v != v else round(float(v), 4) for v in rec],
        ])
    return ["genome", *value_cols], rows


CLUSTER_CSV_HEADER: list[str] = [
    "cluster_id", "taxonomic_call", "read_count", "suggested_name", "nearest_genome",
    "nearest_species", "nearest_genus", "nearest_family", "closest_known_taxon",
    "estimated_ani", "mean_novelty_index", "mean_placement_uncertainty", "confidence",
    "phylogenetic_context", "neighborhood_context", "placement_support",
    "novel_species_reads", "novel_genus_reads", "contributing_genomes", "grouping",
]


def clusters_csv(clusters: list[Any]) -> tuple[list[str], list[list[Any]]]:
    """Header and rows for the novel-cluster CSV download."""
    rows: list[list[Any]] = []
    for c in clusters:
        nbr = c.neighborhood
        rows.append([
            c.cluster_id, c.taxonomic_call, c.read_count, c.suggested_name, c.nearest_genome,
            c.nearest_species, c.nearest_genus, c.nearest_family, c.closest_known_taxon,
            c.estimated_ani, round(c.mean_novelty_index, 3), round(c.mean_placement_uncertainty, 3),
            c.confidence, c.phylogenetic_context,
            nbr.phylogenetic_context if nbr is not None else None,
            round(nbr.placement_support, 1) if nbr is not None else None,
            getattr(c, "novel_species_reads", None), getattr(c, "novel_genus_reads", None),
            ";".join(getattr(c, "contributing_genomes", None) or []),
            getattr(c, "grouping", None),
        ])
    return CLUSTER_CSV_HEADER, rows


SORTABLE_TABLE_JS: str = r"""
// Client-side sorting for tables marked data-sortable (values in data-v)
document.addEventListener('DOMContentLoaded', function() {
    document.querySelectorAll('table[data-sortable]').forEach(function(table) {
        var heads = Array.prototype.slice.call(table.querySelectorAll('thead th'));
        heads.forEach(function(th, col) {
            var btn = th.querySelector('button.sort-btn');
            if (!btn) return;
            var icon = th.querySelector('.sort-icon');
            if (icon && th.getAttribute('aria-sort') === 'descending') icon.textContent = '▼';
            btn.addEventListener('click', function() {
                var asc = th.getAttribute('aria-sort') !== 'ascending';
                var numeric = th.getAttribute('data-type') === 'num';
                var tbody = table.tBodies[0];
                var rows = Array.prototype.slice.call(tbody.rows);
                var val = function(r) {
                    var v = r.cells[col].getAttribute('data-v');
                    if (numeric) return v === '' || v === null ? null : parseFloat(v);
                    return v === null ? '' : v;
                };
                rows.sort(function(a, b) {
                    var x = val(a), y = val(b);
                    if (x === null && y === null) return 0;
                    if (x === null) return 1;
                    if (y === null) return -1;
                    var c = numeric ? x - y : (x < y ? -1 : (x > y ? 1 : 0));
                    return asc ? c : -c;
                });
                rows.forEach(function(r) { tbody.appendChild(r); });
                heads.forEach(function(h) {
                    h.setAttribute('aria-sort', 'none');
                    var i = h.querySelector('.sort-icon');
                    if (i) i.textContent = '';
                });
                th.setAttribute('aria-sort', asc ? 'ascending' : 'descending');
                if (icon) icon.textContent = asc ? '▲' : '▼';
            });
        });
    });
});
"""
