"""Regression tests for the round-4 audit fixes."""

from __future__ import annotations

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


# ------------------------------------------------- validators, seed, kraken2


def _validate(tmp_path, name, content: bytes, kind="fastq"):
    p = tmp_path / name
    p.write_bytes(content)
    return runner.invoke(app, ["validate", kind, "-i", str(p)])


def test_fastq_trailing_blank_line_is_valid(tmp_path):
    res = _validate(tmp_path, "t.fq", b"@r1\nACGT\n+\nIIII\n\n")
    assert res.exit_code == 0, res.output


@pytest.mark.parametrize("content", [b"@r1\nACGT\n+\nIIII\n@r2\n\n+\n\n",
                                     b"@r1\nAC>GT\n+\nIIIII\n"])
def test_fastq_empty_or_invalid_sequence_rejected(tmp_path, content):
    assert _validate(tmp_path, "e.fq", content).exit_code == 1


@pytest.mark.parametrize("kind,record", [("fasta", b">a\nACGT\n" * 200),
                                         ("fastq", b"@a\nACGT\n+\nIIII\n" * 200)])
def test_truncated_gzip_reported_as_issue(tmp_path, kind, record):
    import gzip

    data = gzip.compress(record)
    res = _validate(tmp_path, f"x.{kind}.gz", data[: len(data) // 2], kind)
    assert res.exit_code == 1
    assert "truncated or corrupt" in res.output
    assert "Unexpected error" not in res.output


@pytest.mark.parametrize("raw", ["-5", "99999999999999999999"])
def test_out_of_range_seed_falls_back(monkeypatch, raw):
    from metadarkmatter.core.random import DEFAULT_SEED, get_seed

    monkeypatch.setenv("METADARKMATTER_SEED", raw)
    assert get_seed() == DEFAULT_SEED


def _kraken_inputs(tmp_path):
    kr = tmp_path / "s.kraken"
    kr.write_text("C\tr1\t263\t150\t263:116\nC\tr2\t262\t150\t262:116\n"
                  "C\tr3\t262\t150\t262:116\n")
    (tmp_path / "s.kreport").write_text(
        " 100.00\t3\t0\tF\t34064\tFrancisellaceae\n"
        "100.00\t3\t2\tG\t262\t  Francisella\n"
        " 33.33\t1\t1\tS\t263\t    Francisella tularensis\n")
    reads = tmp_path / "s_R1.fastq"
    reads.write_text("".join(f"@r{i}\nACGT\n+\nIIII\n" for i in (1, 2, 3)))
    return kr, reads


def test_kraken_extract_dry_run_creates_nothing(tmp_path):
    from metadarkmatter.external.base import ExternalTool

    kr, reads = _kraken_inputs(tmp_path)
    ExternalTool.set_executable_resolver(lambda name: f"/usr/bin/{name}")
    ExternalTool.clear_cache()
    try:
        res = runner.invoke(app, ["kraken2", "extract", "-k", str(kr), "-1", str(reads),
                                  "-t", "262", "-o", str(tmp_path / "out"), "--dry-run"])
    finally:
        ExternalTool.reset_executable_resolver()
    assert res.exit_code == 0, res.output
    assert not (tmp_path / "out").exists()


def test_kraken_exact_match_expected_reads(tmp_path):
    from metadarkmatter.external.kraken import KrakenReport

    _kraken_inputs(tmp_path)
    report = KrakenReport(tmp_path / "s.kreport")
    assert (report.get_taxon_reads(262), report.get_clade_reads(262)) == (2, 3)


# ---------------------------------------------------------- analysis tools


def test_adaptive_boundary_outside_species_range_falls_back():
    pytest.importorskip("sklearn")
    from metadarkmatter.core.classification.adaptive import detect_species_boundary
    from metadarkmatter.core.classification.ani_matrix import ANIMatrix

    # Two groups of genomes at ~78% and ~88% ANI: a genus-level split.
    names = [f"G{i}" for i in range(8)]
    m = {a: {b: (100.0 if a == b else (88.0 if (i < 4) == (j < 4) else 78.0) + 0.1 * ((i + j) % 3))
             for j, b in enumerate(names)} for i, a in enumerate(names)}
    res = detect_species_boundary(ANIMatrix(m))
    assert res.method == "fallback" and res.species_boundary == pytest.approx(96.0)


def test_adaptive_config_keeps_preset_when_boundary_too_high():
    from metadarkmatter.core.classification.adaptive import (
        AdaptiveThresholds,
        build_adaptive_config,
    )
    from metadarkmatter.models.config import ScoringConfig

    base = ScoringConfig(novelty_novel_species_max=15.0, novelty_novel_genus_min=15.0)
    ad = AdaptiveThresholds(species_boundary=83.0, novelty_known_max=17.0, confidence=0.7,
                            method="gmm", ani_values_used=45)
    assert build_adaptive_config(base, ad) == base


def _metrics_df():
    import polars as pl

    return pl.DataFrame({
        "read_id": ["a", "b"], "novelty_index": [3.0, 12.0], "placement_uncertainty": [0.5, 0.2],
        "num_ambiguous_hits": [1, 1], "identity_gap": [None, None], "taxonomic_call": ["x", "y"],
    })


@pytest.mark.parametrize("kw", [{"novelty_range": (2.0, 25.0)}, {"uncertainty_range": (0.5, 8.0)}])
def test_sensitivity_range_outside_config_is_clear_error(kw):
    from metadarkmatter.core.classification.sensitivity import run_sensitivity_analysis
    from metadarkmatter.models.config import ScoringConfig

    with pytest.raises(ConfigurationError, match="Sweep range outside the configuration"):
        run_sensitivity_analysis(_metrics_df(), ScoringConfig(), **kw)


def test_sensitivity_protein_mode_is_error():
    from metadarkmatter.core.classification.sensitivity import run_sensitivity_analysis
    from metadarkmatter.models.config import ScoringConfig

    with pytest.raises(ConfigurationError, match="protein mode"):
        run_sensitivity_analysis(_metrics_df(), ScoringConfig(alignment_mode="protein"))


def test_ece_ignores_rows_without_confidence():
    import polars as pl

    from metadarkmatter.core.classification.evaluation import evaluate_classifications

    df = pl.DataFrame({"true_category": ["A", "A"], "taxonomic_call": ["A", "B"],
                       "confidence_score": [float("nan"), 100.0]})
    assert evaluate_classifications(df).expected_calibration_error == pytest.approx(1.0)


# ---------------------------------------------------------- ANI / AAI / protein


def test_skani_full_matrix_paths_with_spaces(tmp_path):
    from metadarkmatter.core.ani_matrix_builder import parse_skani_output

    out = tmp_path / "skani.txt"
    out.write_text("2\n/g dir/GCF_000000001.1.fna\t100.0\t91.5\n"
                   "/g dir/GCF_000000002.1.fna\t91.5\t100.0\n")
    m = parse_skani_output(out, full_matrix=True)
    assert m["GCF_000000001.1"]["GCF_000000002.1"] == 91.5


def test_fastani_keeps_genomes_without_rows(tmp_path):
    from metadarkmatter.core.ani_matrix_builder import parse_fastani_output

    out = tmp_path / "fa.txt"
    out.write_text("/g/A.fna\t/g/B.fna\t95.0\t10\t12\n/g/B.fna\t/g/A.fna\t95.2\t10\t12\n")
    m = parse_fastani_output(out, genomes=["A", "B", "C"])
    assert set(m) == {"A", "B", "C"} and m["A"]["C"] == 0.0 and m["C"]["C"] == 100.0


def test_ani_matrix_csv_named_tsv_is_read(tmp_path):
    from metadarkmatter.core.classification.ani_matrix import ANIMatrix

    p = tmp_path / "ani.tsv"
    p.write_text("genome,A,B\nA,100,92\nB,92,100\n")
    assert ANIMatrix.from_file(p).get_ani("A", "B") == 92.0


def test_aai_rejects_fraction_scale_and_nan(tmp_path):
    from metadarkmatter.core.aai_matrix_builder import parse_aai_csv

    p = tmp_path / "aai.csv"
    p.write_text("genome,A,B\nA,1.0,0.75\nB,0.75,1.0\n")
    with pytest.raises(ValueError, match="fractions"):
        parse_aai_csv(p)
    p.write_text("genome,A,B\nA,100,NaN\nB,NaN,100\n")
    assert parse_aai_csv(p)["A"]["B"] == 0.0


def test_dmnd_path_keeps_dotted_prefix():
    from metadarkmatter.cli.blastx import _dmnd_path

    assert _dmnd_path(Path("bx/pan.v2")).name == "pan.v2.dmnd"
    assert _dmnd_path(Path("bx/pan.dmnd")).name == "pan.dmnd"


@pytest.mark.parametrize("cmd", [["blast", "align"], ["blastx", "align"]])
def test_short_q_means_quiet_not_query(cmd):
    from typer.main import get_command

    group = get_command(app)
    params = group.commands[cmd[0]].commands[cmd[1]].params
    query = next(p for p in params if p.name == "query")
    assert "-q" not in query.opts


def test_genome_labels_agree_across_ani_and_blastdb():
    from metadarkmatter.core.aai_matrix_builder import extract_genome_name_from_filename
    from metadarkmatter.core.ani_matrix_builder import extract_genome_name_from_path
    from metadarkmatter.core.genome_utils import extract_accession_from_filename

    for name in ["strainX_genomic.fna", "Fam_Gen_sp_GCA_000710735.1.fasta", "custom.fa"]:
        assert extract_accession_from_filename(name) == extract_genome_name_from_path(name)
    assert extract_genome_name_from_filename(
        Path("Fam_Gen_sp_GCA_000710735.1.faa")) == "GCA_000710735.1"


def test_aai_compute_honours_protein_pattern(tmp_path, monkeypatch):
    from metadarkmatter.core import aai_matrix_builder as amb

    prot = tmp_path / "prot"
    prot.mkdir()
    for name in ("GCF_000000001.1.faa", "GCF_000000002.1.faa", "GCF_000000003.1.pep",
                 "GCF_000000004.1.fasta"):
        (prot / name).write_text(">p\nMKV\n")
    out = tmp_path / "concat.faa"
    gm = amb.concatenate_protein_files(prot, out, patterns=["*.pep"])
    assert set(gm) == {"GCF_000000003.1"}
    gm = amb.concatenate_protein_files(prot, out)
    # The default no longer includes *.fasta (nucleotide genomes live there too).
    assert set(gm) == {"GCF_000000001.1", "GCF_000000002.1"}
