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
from typing import Final

import polars as pl

from metadarkmatter.external.ncbi_datasets import AssemblySummary

UNRESOLVED: Final = "unresolved"

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
    """
    labels: list[list[str]] = []
    for call, genomes in zip(
        df["taxonomic_call"].to_list(), df["top_hit_genomes"].to_list(), strict=True
    ):
        if call != "Known Species" or not genomes:
            labels.append([])
            continue
        labels.append(sorted({subspecies_map[g] for g in genomes if g in subspecies_map}))

    closest = [
        None if not lab else (lab[0] if len(lab) == 1 else UNRESOLVED) for lab in labels
    ]
    candidates = [";".join(lab) if lab else None for lab in labels]
    return df.drop("top_hit_genomes").with_columns(
        pl.Series("closest_subspecies", closest, dtype=pl.Utf8),
        pl.Series("subspecies_candidates", candidates, dtype=pl.Utf8),
    )
