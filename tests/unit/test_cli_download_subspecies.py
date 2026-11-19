"""Tests for 'download genomes list --ncbi-subspecies'."""

from __future__ import annotations

from unittest.mock import patch

import polars as pl

from metadarkmatter.cli.main import app
from metadarkmatter.clients.gtdb import GTDBAPIError, GTDBGenome, GTDBQueryResult
from metadarkmatter.external.base import ToolExecutionError, ToolNotFoundError
from metadarkmatter.external.ncbi_datasets import AssemblySummary

TAX = "d__Bacteria;p__P;c__C;o__O;f__F;g__Genus;s__Genus alpha"
TAX_B = "d__Bacteria;p__P;c__C;o__O;f__F;g__Genus;s__Genus beta"

REP = "GCF_000000001.1"
OTHER_A = "GCF_000000002.1"
OTHER_B = "GCF_000000003.1"
OTHER_C = "GCA_000000004.1"
BETA = "GCF_000000009.1"


def _result(genomes: list[GTDBGenome]) -> GTDBQueryResult:
    species_counts: dict[str, int] = {}
    genus_counts: dict[str, int] = {}
    for g in genomes:
        species_counts[g.species] = species_counts.get(g.species, 0) + 1
        genus_counts["Genus"] = genus_counts.get("Genus", 0) + 1
    return GTDBQueryResult(
        taxon="g__Genus",
        genomes=tuple(genomes),
        total_count=len(genomes),
        genus_counts=genus_counts,
        species_counts=species_counts,
    )


def _g(acc: str, species: str = "Genus alpha", tax: str = TAX) -> GTDBGenome:
    return GTDBGenome(accession=acc, gtdb_taxonomy=tax, species=species, genome_size=1000)


def _s(acc: str, org: str, level: str = "Complete Genome", n50: int = 1000) -> AssemblySummary:
    return AssemblySummary(
        accession=acc,
        organism_name=org,
        assembly_level=level,
        refseq_category="",
        contig_n50=n50,
    )


REPS_ONLY = _result([_g(REP), _g(BETA, "Genus beta", TAX_B)])
ALL = _result(
    [_g(REP), _g(OTHER_A), _g(OTHER_B), _g(OTHER_C), _g(BETA, "Genus beta", TAX_B)]
)
SUMMARIES = [
    _s(REP, "Genus alpha strain1"),
    _s(OTHER_A, "Genus alpha strain2", n50=900),
    _s(OTHER_B, "Genus gamma strain3"),
    _s(OTHER_C, "Genus gamma strain4", level="Contig", n50=10),
]


def _query(taxon: str, *, representatives_only: bool = True) -> GTDBQueryResult:
    if taxon == "g__Genus":
        return REPS_ONLY if representatives_only else ALL
    assert taxon == "s__Genus alpha"
    assert representatives_only is False
    return _result([_g(REP), _g(OTHER_A), _g(OTHER_B), _g(OTHER_C)])


def _run(cli_runner, tmp_path, extra, *, summaries=None, side_effect=None, base=None):
    out = tmp_path / "list.tsv"
    args = ["download", "genomes", "list", "g__Genus", "-o", str(out), *(base or []), *extra]
    with (
        patch("metadarkmatter.cli.download.GTDBClient.query_genomes", side_effect=_query),
        patch(
            "metadarkmatter.cli.download.NCBIDatasets.fetch_assembly_summaries",
            return_value=SUMMARIES if summaries is None else summaries,
            side_effect=side_effect,
        ),
    ):
        return cli_runner.invoke(app, args)


def _meta(tmp_path) -> pl.DataFrame:
    return pl.read_csv(tmp_path / "genome_metadata.tsv", separator="\t")


def test_list_ncbi_subspecies_reps_only_appends_panel(cli_runner, tmp_path):
    result = _run(cli_runner, tmp_path, ["--ncbi-subspecies", "s__Genus alpha"])
    assert result.exit_code == 0, result.output
    df = _meta(tmp_path)
    alpha = df.filter(pl.col("species") == "Genus alpha")
    assert set(alpha["accession"]) == {REP, OTHER_A, OTHER_B, OTHER_C}
    assert set(alpha["representative"]) == {REP}
    assert "subspecies" in df.columns
    labels = dict(zip(alpha["accession"], alpha["subspecies"], strict=True))
    assert labels[REP] == "Genus alpha"
    assert labels[OTHER_A] == "Genus alpha"
    assert labels[OTHER_B] == "Genus gamma"
    beta = df.filter(pl.col("species") == "Genus beta")
    assert beta["representative"].to_list() == [BETA]
    assert (tmp_path / "subspecies_panel.tsv").exists()
    # Appended genomes follow the original ones in sorted order.
    appended = df["accession"].to_list()[2:]
    assert appended == sorted(appended)


def test_list_ncbi_subspecies_all_genomes_labels_all(cli_runner, tmp_path):
    result = _run(
        cli_runner, tmp_path, ["--ncbi-subspecies", "s__Genus alpha"], base=["--all-genomes"]
    )
    assert result.exit_code == 0, result.output
    df = _meta(tmp_path)
    alpha = df.filter(pl.col("species") == "Genus alpha")
    assert len(alpha) == 4
    labels = dict(zip(alpha["accession"], alpha["subspecies"], strict=True))
    assert labels[OTHER_C] == "Genus gamma"
    assert labels[REP] == "Genus alpha"
    assert (tmp_path / "subspecies_panel.tsv").exists()


def test_list_ncbi_subspecies_unknown_species_exits(cli_runner, tmp_path):
    result = _run(cli_runner, tmp_path, ["--ncbi-subspecies", "s__Genus missing"])
    assert result.exit_code == 1
    assert "s__Genus missing" in result.output
    assert not (tmp_path / "genome_metadata.tsv").exists()


def test_list_ncbi_subspecies_requires_s_prefix(cli_runner, tmp_path):
    result = _run(cli_runner, tmp_path, ["--ncbi-subspecies", "Genus alpha"])
    assert result.exit_code == 1
    assert not (tmp_path / "genome_metadata.tsv").exists()


def test_list_ncbi_subspecies_datasets_missing_exits(cli_runner, tmp_path):
    result = _run(
        cli_runner,
        tmp_path,
        ["--ncbi-subspecies", "s__Genus alpha"],
        side_effect=ToolNotFoundError("datasets", "conda install -c conda-forge ncbi-datasets-cli"),
    )
    assert result.exit_code == 1
    assert "ncbi-datasets-cli" in result.output
    assert not (tmp_path / "genome_metadata.tsv").exists()
    assert not (tmp_path / "subspecies_panel.tsv").exists()


def test_list_ncbi_subspecies_datasets_failure_exits(cli_runner, tmp_path):
    result = _run(
        cli_runner,
        tmp_path,
        ["--ncbi-subspecies", "s__Genus alpha"],
        side_effect=ToolExecutionError("datasets", ["datasets"], 1, "boom"),
    )
    assert result.exit_code == 1
    assert not (tmp_path / "genome_metadata.tsv").exists()
    assert not (tmp_path / "list.tsv").exists()


def test_list_without_subspecies_unchanged(cli_runner, tmp_path):
    result = _run(cli_runner, tmp_path, [])
    assert result.exit_code == 0, result.output
    assert "subspecies" not in _meta(tmp_path).columns
    assert not (tmp_path / "subspecies_panel.tsv").exists()


def test_list_ncbi_subspecies_gtdb_species_query_failure_exits(cli_runner, tmp_path):
    def query(taxon, *, representatives_only=True):
        if taxon.startswith("s__"):
            raise GTDBAPIError("GTDB unavailable")
        return _query(taxon, representatives_only=representatives_only)

    out = tmp_path / "list.tsv"
    with (
        patch("metadarkmatter.cli.download.GTDBClient.query_genomes", side_effect=query),
        patch(
            "metadarkmatter.cli.download.NCBIDatasets.fetch_assembly_summaries",
            return_value=SUMMARIES,
        ),
    ):
        result = cli_runner.invoke(
            app,
            ["download", "genomes", "list", "g__Genus", "-o", str(out),
             "--ncbi-subspecies", "s__Genus alpha"],
        )
    assert result.exit_code == 1
    assert not (tmp_path / "genome_metadata.tsv").exists()
    assert not out.exists()
