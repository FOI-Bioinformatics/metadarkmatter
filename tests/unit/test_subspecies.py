"""Tests for NCBI species parsing and subspecies panel selection."""

from __future__ import annotations

import dataclasses
from pathlib import Path

from metadarkmatter.core.subspecies import (
    PanelEntry,
    ncbi_species_name,
    select_subspecies_panel,
    write_panel_tsv,
)
from metadarkmatter.external.ncbi_datasets import AssemblySummary


def _s(acc, name="Brucella abortus 1", level="Complete Genome", cat="na", n50=1):
    return AssemblySummary(acc, name, level, cat, n50)


def test_ncbi_species_name_cases():
    assert ncbi_species_name("Brucella abortus 2308") == "Brucella abortus"
    assert ncbi_species_name("Brucella sp. 09RB8471") is None
    assert (
        ncbi_species_name("Candidatus Liberibacter asiaticus str. psy62")
        == "Liberibacter asiaticus"
    )
    assert ncbi_species_name("[Brucella] x") is None
    assert ncbi_species_name("Brucella") is None


def test_select_panel_ranking():
    sums = [
        _s("GCF_000000001.1", level="Scaffold", cat="reference genome", n50=1),
        _s("GCA_000000002.1", n50=9_000_000),
        _s("GCF_000000003.1", n50=2_000_000),
        _s("GCF_000000004.1", n50=3_000_000),
    ]
    panel = select_subspecies_panel(
        sums, "Brucella abortus", "GCF_000000001.1", per_subspecies=2
    )
    assert [e.accession for e in panel] == ["GCF_000000001.1", "GCF_000000004.1"]
    assert [e.rank for e in panel] == [1, 2]
    assert panel[0].n_genomes_in_subspecies == 4


def test_select_panel_large_group_keeps_three():
    def many(n):
        return [_s(f"GCF_{i:09d}.1", n50=i) for i in range(1, n + 1)]

    assert len(select_subspecies_panel(many(51), "g", "GCF_000000051.1")) == 3
    assert len(select_subspecies_panel(many(50), "g", "GCF_000000050.1")) == 2


def test_select_panel_always_keeps_representative():
    sums = [_s(f"GCF_{i:09d}.1", n50=100 - i) for i in range(1, 6)]
    rep = "GCF_000000005.1"
    panel = select_subspecies_panel(sums, "g", rep)
    accs = [e.accession for e in panel]
    assert rep in accs
    assert len(panel) == 3
    assert [e.rank for e in panel if e.accession == rep] == [5]


def test_select_panel_skips_sp():
    sums = [
        _s("GCF_000000001.1", name="Brucella sp."),
        _s("GCF_000000002.1", name="Brucella sp. X"),
    ]
    assert select_subspecies_panel(sums, "g", "GCF_999999999.1") == []


def test_select_panel_unparsable_representative_kept_unlabelled():
    sums = [
        _s("GCF_000000001.1", name="Brucella sp. X"),
        _s("GCF_000000002.1"),
    ]
    panel = select_subspecies_panel(sums, "g", "GCF_000000001.1")
    by_acc = {e.accession: e for e in panel}
    assert by_acc["GCF_000000001.1"].subspecies == ""
    assert "GCF_000000002.1" in by_acc
    assert [e.subspecies for e in panel] == ["", "Brucella abortus"]


def test_write_panel_tsv_header(tmp_path: Path):
    entry = PanelEntry("GCF_1", "g", "A b", "A b 1", "Complete Genome", "na", 5, 1, 1)
    out = tmp_path / "p.tsv"
    write_panel_tsv([entry], out)
    lines = out.read_text().splitlines()
    assert lines[0] == "\t".join(f.name for f in dataclasses.fields(PanelEntry))
    assert lines[1].split("\t")[0] == "GCF_1"
