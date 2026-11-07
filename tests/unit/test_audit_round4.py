"""Regression tests for the round-4 audit fixes."""

from __future__ import annotations

import json
from pathlib import Path

import httpx
import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.main import app
from metadarkmatter.clients.gtdb import GTDBAPIError, GTDBClient
from metadarkmatter.core.exceptions import ConfigurationError

runner = CliRunner()


# ------------------------------------------------------------- GTDB client


def _client(tmp_path: Path, handler) -> GTDBClient:
    c = GTDBClient(cache_dir=tmp_path / "cache", max_retries=1, retry_delay=0.0)
    c._client = httpx.Client(transport=httpx.MockTransport(handler), base_url="https://x")
    return c


@pytest.mark.parametrize("content", [b"\xff\xfe", b"[]", b"null", b'"x"',
                                     b'{"timestamp": "abc", "data": {"rows": [1]}}'])
def test_bad_cache_entry_is_a_miss(tmp_path, content):
    calls = []

    def handler(request):
        calls.append(request)
        return httpx.Response(200, json={"rows": [{"gid": "G"}]})

    c = _client(tmp_path, handler)
    assert c._get("/t", {"a": "1"}) == {"rows": [{"gid": "G"}]}
    cache_file = c._cache_path("/t", {"a": "1"})
    cache_file.write_bytes(content)
    assert c._get("/t", {"a": "1"}) == {"rows": [{"gid": "G"}]}
    assert len(calls) == 2


def test_non_json_response_is_api_error(tmp_path):
    c = _client(tmp_path, lambda r: httpx.Response(200, text="<html>proxy</html>"))
    with pytest.raises(GTDBAPIError, match="not JSON"):
        c._get("/t", {})


def test_retry_after_http_date_does_not_crash(tmp_path):
    responses = iter([
        httpx.Response(429, headers={"Retry-After": "Wed, 21 Oct 2026 07:28:00 GMT"}),
        httpx.Response(200, json={"rows": [{"gid": "G"}]}),
    ])
    c = _client(tmp_path, lambda r: next(responses))
    assert c._get("/t", {}) == {"rows": [{"gid": "G"}]}


def test_empty_result_not_cached(tmp_path):
    c = _client(tmp_path, lambda r: httpx.Response(200, json={"rows": []}))
    c._get("/t", {})
    assert not c._cache_path("/t", {}).exists()


# ---------------------------------------------------- metadata and mapping


def test_select_representatives_without_genus_or_family(tmp_path):
    m = tmp_path / "m.tsv"
    m.write_text("accession\tspecies\nGCF_000000001.1\tA b\nGCF_000000002.1\tA b\n")
    out = tmp_path / "o.tsv"
    res = runner.invoke(app, ["util", "select-representatives", "-m", str(m), "-o", str(out)])
    assert res.exit_code == 0, res.output
    assert out.read_text().splitlines()[0].split("\t")[:3] == ["accession", "species", "representative"]


def test_family_override_drops_family_token():
    from metadarkmatter.cli.mapping import _parse_structured_filename

    r = _parse_structured_filename(
        "Francisellaceae_Francisella_tularensis_GCF_000195955.2.fasta", "Francisellaceae")
    assert (r["family"], r["genus"], r["species"]) == (
        "Francisellaceae", "Francisella", "Francisella tularensis")
    r = _parse_structured_filename("Francisella_tularensis_GCF_000195955.2.fasta", "Fam")
    assert (r["genus"], r["species"]) == ("Francisella", "Francisella tularensis")


def test_accession_list_with_empty_species(tmp_path):
    from metadarkmatter.models.genomes import AccessionList

    p = tmp_path / "g.tsv"
    p.write_text("accession\tgtdb_taxonomy\tspecies\n"
                 "GCF_000000001.1\td__B;p__P;c__C;o__O;f__F;g__G;s__\t\n")
    acc = AccessionList.from_tsv(p).accessions[0]
    assert acc.species == ""


def test_download_list_refuses_output_named_genome_metadata(tmp_path, monkeypatch):
    from metadarkmatter.clients import gtdb

    def no_network(*a, **k):  # pragma: no cover - must not be reached
        raise AssertionError("queried GTDB before validating --output")

    monkeypatch.setattr(gtdb.GTDBClient, "_get", no_network)
    out = tmp_path / "genome_metadata.tsv"
    res = runner.invoke(app, ["download", "genomes", "list", "f__Fam", "-o", str(out)])
    assert res.exit_code == 1
    assert "must not be named genome_metadata.tsv" in res.output
    assert not out.exists()


def test_empty_fasta_header_is_an_error(tmp_path):
    from metadarkmatter.core.id_mapping import ContigIdMapping

    g = tmp_path / "gen"
    g.mkdir()
    (g / "GCF_000000001.1.fna").write_text(">ctg3\nAC\n>\nAA\n>ctg4\nCC\n")
    with pytest.raises(ConfigurationError, match="Empty FASTA header"):
        ContigIdMapping.from_genome_dir(g)
    res = runner.invoke(app, ["util", "generate-mapping", "-g", str(g),
                              "-o", str(tmp_path / "m.tsv")])
    assert res.exit_code != 0
