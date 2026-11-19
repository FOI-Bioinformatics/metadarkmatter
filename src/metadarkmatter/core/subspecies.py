"""NCBI species labels within a GTDB species ("subspecies").

This module derives an NCBI species name from an assembly organism name and
selects a small panel of representative genomes per NCBI species within one
GTDB species.
"""

from __future__ import annotations

import csv
import dataclasses
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Final

import polars as pl

from metadarkmatter.external.ncbi_datasets import AssemblySummary

if TYPE_CHECKING:
    from metadarkmatter.core.metadata import GenomeMetadata

UNRESOLVED: Final = "unresolved"
UNLABELLED: Final = "unlabelled"
_KEY_SEP: Final = "\t"
_UNKNOWN_SPECIES: Final = "Unknown species"

_LEVEL_ORDER: Final = {
    "complete genome": 0,
    "chromosome": 1,
    "scaffold": 2,
    "contig": 3,
}
_OTHER_LEVEL: Final = 4


@dataclass(frozen=True)
class PanelEntry:
    """One genome selected for the subspecies panel."""

    accession: str
    gtdb_species: str
    subspecies: str
    ncbi_organism_name: str
    assembly_level: str
    refseq_category: str
    contig_n50: int
    rank: int
    n_genomes_in_subspecies: int


def ncbi_species_name(organism_name: str) -> str | None:
    """Return the first two name tokens as the NCBI species, or None.

    A leading "Candidatus " is removed. None is returned when fewer than two
    tokens remain, when the second token is "sp.", or when the first token is
    not a capitalised alphabetic genus name (for example "[Brucella]").
    """
    tokens = organism_name.split()
    if tokens and tokens[0] == "Candidatus":
        tokens = tokens[1:]
    if len(tokens) < 2:
        return None
    genus, species = tokens[0], tokens[1]
    if species == "sp.":
        return None
    if not (genus.isalpha() and genus[0].isupper()):
        return None
    return f"{genus} {species}"


def _rank_key(s: AssemblySummary) -> tuple[int, int, int, int, str]:
    return (
        0 if s.refseq_category.strip().lower() == "reference genome" else 1,
        _LEVEL_ORDER.get(s.assembly_level.strip().lower(), _OTHER_LEVEL),
        0 if s.accession.startswith("GCF_") else 1,
        -s.contig_n50,
        s.accession,
    )


def _entry(
    s: AssemblySummary, gtdb_species: str, subspecies: str, rank: int, n: int
) -> PanelEntry:
    return PanelEntry(
        accession=s.accession,
        gtdb_species=gtdb_species,
        subspecies=subspecies,
        ncbi_organism_name=s.organism_name,
        assembly_level=s.assembly_level,
        refseq_category=s.refseq_category,
        contig_n50=s.contig_n50,
        rank=rank,
        n_genomes_in_subspecies=n,
    )


def select_subspecies_panel(
    summaries: list[AssemblySummary],
    gtdb_species: str,
    gtdb_representative: str,
    *,
    per_subspecies: int = 2,
    per_subspecies_large: int = 3,
    large_threshold: int = 50,
) -> list[PanelEntry]:
    """Select top-ranked genomes per NCBI species within one GTDB species.

    Genomes are grouped by NCBI species name and ranked by: RefSeq reference
    genome first, assembly level (Complete Genome, Chromosome, Scaffold,
    Contig, other), GCF before GCA, contig N50 descending, accession
    ascending. The top ``per_subspecies`` are kept, or ``per_subspecies_large``
    when the group has more than ``large_threshold`` genomes.

    The GTDB representative is always included. If its organism name has no
    parsable NCBI species it is kept with an empty subspecies label. If no
    summary exists for it, it is not added. Output is sorted by
    (subspecies, rank).
    """
    groups: dict[str, list[AssemblySummary]] = {}
    for s in summaries:
        name = ncbi_species_name(s.organism_name)
        if name is not None:
            groups.setdefault(name, []).append(s)

    entries: list[PanelEntry] = []
    for name, members in groups.items():
        ranked = sorted(members, key=_rank_key)
        keep = per_subspecies_large if len(ranked) > large_threshold else per_subspecies
        for i, s in enumerate(ranked, start=1):
            if i <= keep or s.accession == gtdb_representative:
                entries.append(_entry(s, gtdb_species, name, i, len(ranked)))

    if not any(e.accession == gtdb_representative for e in entries):
        rep = next((s for s in summaries if s.accession == gtdb_representative), None)
        if rep is not None:
            entries.append(_entry(rep, gtdb_species, "", 1, 1))

    entries.sort(key=lambda e: (e.subspecies, e.rank, e.accession))
    return entries


def write_panel_tsv(entries: list[PanelEntry], path: Path) -> None:
    """Write panel entries as a tab-separated file in dataclass field order."""
    fields = [f.name for f in dataclasses.fields(PanelEntry)]
    with Path(path).open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for e in entries:
            writer.writerow([getattr(e, f) for f in fields])


def assign_subspecies(df: pl.DataFrame, subspecies_map: dict[str, str]) -> pl.DataFrame:
    """Add per-read NCBI species labels from the tied top-hit genomes.

    The input needs ``taxonomic_call`` and ``top_hit_genomes`` (list of
    accessions). Labels are the sorted unique ``subspecies_map`` values of the
    top-hit genomes that have a label. Reads whose call is not "Known Species",
    or with no labelled top-hit genome (including a null list), get null in
    both new columns. One label gives that label in ``closest_subspecies``;
    several give ``UNRESOLVED``. ``subspecies_candidates`` holds all labels
    joined by ";". The ``top_hit_genomes`` column is dropped.

    Labels are computed once per distinct top-hit genome list, so the cost
    scales with the number of distinct lists rather than the number of reads.
    """
    out = df.drop("top_hit_genomes")
    is_known = pl.col("taxonomic_call") == "Known Species"
    if not subspecies_map or df.is_empty() or not df.select(is_known.any()).item():
        return out.with_columns(
            pl.lit(None, dtype=pl.Utf8).alias("closest_subspecies"),
            pl.lit(None, dtype=pl.Utf8).alias("subspecies_candidates"),
        )

    # One string key per row (null for non-Known-Species rows and null lists).
    # Accessions come from tab-separated alignment fields, so they cannot
    # contain the tab separator.
    keys = df.select(
        pl.when(is_known)
        .then(pl.col("top_hit_genomes").list.join(_KEY_SEP))
        .otherwise(None)
        .alias("_key")
    )["_key"]

    closest: dict[str, str] = {}
    candidates: dict[str, str] = {}
    for key in keys.drop_nulls().unique().sort().to_list():
        if not key:
            continue
        labels = sorted(
            {subspecies_map[g] for g in key.split(_KEY_SEP) if g in subspecies_map}
        )
        if labels:
            closest[key] = labels[0] if len(labels) == 1 else UNRESOLVED
            candidates[key] = ";".join(labels)

    null_utf8 = pl.lit(None, dtype=pl.Utf8)
    return out.with_columns(
        keys.replace_strict(closest, default=null_utf8, return_dtype=pl.Utf8).alias(
            "closest_subspecies"
        ),
        keys.replace_strict(candidates, default=null_utf8, return_dtype=pl.Utf8).alias(
            "subspecies_candidates"
        ),
    )


_SUMMARY_SCHEMA: Final = {
    "gtdb_species": pl.Utf8,
    "subspecies": pl.Utf8,
    "reads": pl.UInt32,
    "pct": pl.Float64,
}


def summarize_subspecies(
    df: pl.DataFrame, metadata: GenomeMetadata | None
) -> pl.DataFrame:
    """Count Known Species reads per subspecies within each GTDB species.

    Columns: gtdb_species, subspecies, reads, pct. The percentage is relative
    to all Known Species reads of the GTDB species, including reads without a
    subspecies label, which are listed as "unlabelled". The result is sorted by
    gtdb_species, reads (descending) and subspecies, and is empty when no subspecies column is present.
    """
    if "closest_subspecies" not in df.columns:
        return pl.DataFrame(schema=_SUMMARY_SCHEMA)

    known = df
    if "taxonomic_call" in df.columns:
        known = df.filter(pl.col("taxonomic_call") == "Known Species")

    if "species" in known.columns:
        species = known["species"].cast(pl.Utf8).fill_null(_UNKNOWN_SPECIES)
    else:
        species = pl.Series(
            "species",
            [
                (metadata.get_species(g) if metadata is not None else None)
                or _UNKNOWN_SPECIES
                for g in known["best_match_genome"].to_list()
            ],
            dtype=pl.Utf8,
        )
    known = known.select(
        species.alias("gtdb_species"), pl.col("closest_subspecies").alias("subspecies")
    )
    if known.is_empty():
        return pl.DataFrame(schema=_SUMMARY_SCHEMA)

    totals = known.group_by("gtdb_species").len().rename({"len": "_total"})
    return (
        known.with_columns(pl.col("subspecies").fill_null(UNLABELLED))
        .group_by("gtdb_species", "subspecies")
        .len()
        .rename({"len": "reads"})
        .join(totals, on="gtdb_species")
        .with_columns((pl.col("reads") / pl.col("_total") * 100).alias("pct"))
        .drop("_total")
        .sort(["gtdb_species", "reads", "subspecies"], descending=[False, True, False])
    )
