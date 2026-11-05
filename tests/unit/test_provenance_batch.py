"""Tests for the provenance sidecar and for 'score batch' parity with classify."""

from __future__ import annotations

import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.main import app
from metadarkmatter.cli.score import (
    DEFAULT_BATCH_PATTERNS,
    sample_name_from_alignment,
)
from metadarkmatter.core.provenance import (
    PROVENANCE_SCHEMA,
    build_provenance,
    provenance_path_for,
    sha256_file,
    summarize_classifications,
    write_provenance,
)
from metadarkmatter.models.config import ScoringConfig

runner = CliRunner()

GENOMES = ["GCF_000000001.1", "GCF_000000002.1", "GCF_000000003.1"]
FILTERS = ["--min-alignment-length", "0", "--min-alignment-fraction", "0"]


def _write_ani(path: Path) -> Path:
    vals = [[100.0, 96.0, 80.0], [96.0, 100.0, 79.0], [80.0, 79.0, 100.0]]
    rows = {"genome": GENOMES}
    for j, g in enumerate(GENOMES):
        rows[g] = [vals[i][j] for i in range(3)]
    pl.DataFrame(rows).write_csv(path)
    return path


def _write_metadata(path: Path) -> Path:
    path.write_text(
        "accession\tspecies\tgenus\tfamily\trepresentative\n"
        + "".join(
            f"{g}\tSpecies {i}\tGenus {i}\tf__Testaceae\t{g}\n"
            for i, g in enumerate(GENOMES)
        )
    )
    return path


def _write_alignment(path: Path, n_reads: int = 6) -> Path:
    """Reads with a top hit at varying identity and a competing second hit."""
    lines = []
    for i in range(n_reads):
        ident = 99.0 - 4.0 * i
        for g_idx, genome in enumerate(GENOMES[:2]):
            pident = ident - g_idx * 1.0
            bits = 250.0 - 4.0 * i - g_idx
            lines.append(
                f"read_{i}/1\t{genome}|contig1\t{pident:.3f}\t150\t2\t0\t1\t150\t"
                f"100\t250\t1e-50\t{bits:.1f}\t150"
            )
    text = "\n".join(lines) + "\n"
    if path.name.endswith(".gz"):
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return path


@pytest.fixture
def data(tmp_path: Path) -> dict[str, Path]:
    return {
        "ani": _write_ani(tmp_path / "ani.csv"),
        "metadata": _write_metadata(tmp_path / "meta.tsv"),
        "aln": _write_alignment(tmp_path / "s1.blast.tsv.gz"),
    }


def _classify(data, out: Path, *extra: str, summary: Path | None = None):
    args = [
        "score", "classify", "--alignment", str(data["aln"]), "--ani", str(data["ani"]),
        "--output", str(out), *FILTERS, *extra,
    ]
    if summary is not None:
        args += ["--summary", str(summary)]
    res = runner.invoke(app, args)
    assert res.exit_code == 0, res.output
    return res


# ---------------------------------------------------------------- module unit


def test_sha256_file_streams_and_matches_hashlib(tmp_path):
    p = tmp_path / "big.bin"
    payload = b"abc" * (1024 * 1024 + 17)  # spans several 1 MiB chunks
    p.write_bytes(payload)
    assert sha256_file(p) == hashlib.sha256(payload).hexdigest()


def test_build_provenance_schema_and_null_slots(tmp_path):
    f = tmp_path / "a.tsv"
    f.write_text("x\n")
    prov = build_provenance(
        ["mdm", "score", "classify"],
        {"alignment": f, "ani": f},
        ScoringConfig(),
        {"outputs": {"classifications": "/x", "summary": None}, "results": {}},
    )
    assert prov["schema"] == PROVENANCE_SCHEMA
    assert prov["command"] == ["mdm", "score", "classify"]
    assert isinstance(prov["seed"], int)
    assert set(prov["inputs"]) == {"alignment", "ani", "metadata", "aai", "id_mapping"}
    assert prov["inputs"]["metadata"] is None
    rec = prov["inputs"]["alignment"]
    assert rec == {
        "path": str(f.resolve()),
        "size": 2,
        "sha256": hashlib.sha256(b"x\n").hexdigest(),
    }
    assert prov["config"] == ScoringConfig().model_dump(mode="json")
    assert prov["created_utc"].endswith("Z")
    # JSON-serialisable
    json.dumps(prov)


def test_build_provenance_uses_env_seed(monkeypatch):
    monkeypatch.setenv("METADARKMATTER_SEED", "7")
    assert build_provenance([], {}, None)["seed"] == 7
    assert build_provenance([], {}, None)["config"] is None


def test_summarize_classifications_sorted():
    df = pl.DataFrame({"taxonomic_call": ["Novel Species", "Ambiguous", "Novel Species"]})
    out = summarize_classifications(df)
    assert out["total_reads"] == 3
    assert list(out["taxonomic_call_counts"]) == ["Ambiguous", "Novel Species"]
    assert out["taxonomic_call_counts"]["Novel Species"] == 2
    assert summarize_classifications(df.lazy()) == out


def test_write_provenance_roundtrip(tmp_path):
    p = tmp_path / "sub" / "x.provenance.json"
    write_provenance(p, {"a": 1})
    assert json.loads(p.read_text()) == {"a": 1}
    assert provenance_path_for(tmp_path / "c.csv").name == "c.csv.provenance.json"


# ------------------------------------------------------------ classify sidecar


def test_classify_writes_provenance(tmp_path, data):
    out = tmp_path / "out" / "classifications.csv"
    summ = tmp_path / "out" / "summary.json"
    _classify(data, out, "--metadata", str(data["metadata"]), summary=summ)

    side = tmp_path / "out" / "classifications.csv.provenance.json"
    assert side.exists()
    prov = json.loads(side.read_text())
    assert prov["schema"] == "metadarkmatter.provenance/1"
    assert prov["inputs"]["alignment"]["sha256"] == sha256_file(data["aln"])
    assert prov["inputs"]["alignment"]["path"] == str(data["aln"].resolve())
    assert prov["inputs"]["ani"]["size"] == data["ani"].stat().st_size
    assert prov["inputs"]["metadata"]["sha256"] == sha256_file(data["metadata"])
    assert prov["inputs"]["aai"] is None
    assert prov["inputs"]["id_mapping"] is None
    assert prov["outputs"] == {
        "classifications": str(out.resolve()),
        "summary": str(summ.resolve()),
    }

    df = pl.read_csv(out)
    counts = Counter(df["taxonomic_call"].to_list())
    assert prov["results"]["total_reads"] == df.height
    assert prov["results"]["taxonomic_call_counts"] == dict(sorted(counts.items()))
    # Target family was inferred from the metadata and is recorded as applied.
    assert prov["config"]["target_family"] == "f__Testaceae"

    summary = json.loads(summ.read_text())
    assert summary["total_reads"] == df.height
    assert summary["provenance"] == prov


def test_no_provenance_flag_and_output_unchanged(tmp_path, data):
    with_p = tmp_path / "a" / "c.csv"
    without = tmp_path / "b" / "c.csv"
    _classify(data, with_p)
    summ = tmp_path / "b" / "s.json"
    _classify(data, without, "--no-provenance", summary=summ)
    assert not provenance_path_for(without).exists()
    assert "provenance" not in json.loads(summ.read_text())
    # No extra columns / identical table either way.
    assert with_p.read_bytes() == without.read_bytes()


def test_provenance_deterministic_except_timestamp(tmp_path, data):
    out = tmp_path / "c.csv"
    _classify(data, out)
    first = json.loads(provenance_path_for(out).read_text())
    _classify(data, out)
    second = json.loads(provenance_path_for(out).read_text())
    first.pop("created_utc")
    second.pop("created_utc")
    assert first == second


def test_streaming_classify_writes_provenance(tmp_path, data):
    out = tmp_path / "s.csv"
    _classify(data, out, "--streaming", "--chunk-size", "100")
    prov = json.loads(provenance_path_for(out).read_text())
    df = pl.read_csv(out)
    assert prov["results"]["total_reads"] == df.height > 0


# ----------------------------------------------------------------- batch parity


def test_sample_name_from_alignment():
    assert sample_name_from_alignment(Path("s1.blast.tsv.gz")) == "s1"
    assert sample_name_from_alignment(Path("s1.mmseqs2.tsv.gz")) == "s1"
    assert sample_name_from_alignment(Path("s1.mmseqs2.tsv")) == "s1"
    assert sample_name_from_alignment(Path("a.b.tsv")) == "a.b"
    assert sample_name_from_alignment(Path("weird.out")) == "weird"
    assert "*.mmseqs2.tsv.gz" in DEFAULT_BATCH_PATTERNS


def _batch(data, in_dir: Path, out_dir: Path, *extra: str):
    return runner.invoke(
        app,
        [
            "score", "batch", "--alignment-dir", str(in_dir), "--ani", str(data["ani"]),
            "--output-dir", str(out_dir), *FILTERS, *extra,
        ],
    )


def test_batch_matches_classify_columns_and_values(tmp_path, data):
    """Same input + options => identical classification table from batch."""
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    _write_alignment(in_dir / "blastsamp.blast.tsv.gz")
    _write_alignment(in_dir / "mmseqssamp.mmseqs2.tsv.gz")

    opts = ["--metadata", str(data["metadata"]), "--no-clade-novel-species"]
    res = _batch(data, in_dir, tmp_path / "out", *opts)
    assert res.exit_code == 0, res.output
    assert "Found 2 alignment files" in res.output

    ref = tmp_path / "ref.csv"
    data = {**data, "aln": in_dir / "blastsamp.blast.tsv.gz"}
    _classify(data, ref, *opts)

    ref_df = pl.read_csv(ref)
    assert "species" in ref_df.columns  # metadata was honoured
    for sample in ("blastsamp", "mmseqssamp"):
        got = pl.read_csv(tmp_path / "out" / f"{sample}_classifications.csv")
        assert got.columns == ref_df.columns
        assert got.equals(ref_df)
        prov = json.loads(
            (tmp_path / "out" / f"{sample}_classifications.csv.provenance.json").read_text()
        )
        assert prov["inputs"]["metadata"] is not None
        assert prov["results"]["total_reads"] == ref_df.height
        summary = json.loads((tmp_path / "out" / f"{sample}_summary.json").read_text())
        assert summary["provenance"]["outputs"]["summary"].endswith(f"{sample}_summary.json")


def test_batch_config_flags_reach_classifier(tmp_path, data):
    """--config, --target-family and --clade flags are honoured (classify parity)."""
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    _write_alignment(in_dir / "s.blast.tsv.gz")
    cfg = tmp_path / "cfg.yaml"
    ScoringConfig(include_legacy_scores=True).to_yaml(cfg)
    res = _batch(
        data, in_dir, tmp_path / "o", "--config", str(cfg), "--target-family", "f__X",
        "--no-clade-novel-species",
    )
    assert res.exit_code == 0, res.output
    prov = json.loads((tmp_path / "o" / "s_classifications.csv.provenance.json").read_text())
    assert prov["config"]["include_legacy_scores"] is True
    assert prov["config"]["target_family"] == "f__X"
    assert prov["config"]["clade_novel_species"] is False


def test_batch_streaming_matches_streaming_classify(tmp_path, data):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    _write_alignment(in_dir / "s.blast.tsv.gz")
    res = _batch(data, in_dir, tmp_path / "o", "--streaming", "--chunk-size", "100")
    assert res.exit_code == 0, res.output
    ref = tmp_path / "ref.csv"
    _classify({**data, "aln": in_dir / "s.blast.tsv.gz"}, ref, "--streaming")
    assert pl.read_csv(tmp_path / "o" / "s_classifications.csv").equals(pl.read_csv(ref))


def test_batch_no_provenance(tmp_path, data):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    _write_alignment(in_dir / "s.blast.tsv.gz")
    res = _batch(data, in_dir, tmp_path / "o", "--no-provenance")
    assert res.exit_code == 0
    assert not list((tmp_path / "o").glob("*.provenance.json"))
    summary = json.loads((tmp_path / "o" / "s_summary.json").read_text())
    assert "provenance" not in summary


def test_batch_sample_name_collision_is_an_error(tmp_path, data):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    _write_alignment(in_dir / "s1.blast.tsv.gz")
    _write_alignment(in_dir / "s1.mmseqs2.tsv.gz")
    res = _batch(data, in_dir, tmp_path / "o")
    assert res.exit_code == 1
    assert "same sample name" in res.output
    assert "s1.blast.tsv.gz" in res.output and "s1.mmseqs2.tsv.gz" in res.output
    assert not (tmp_path / "o").exists() or not list((tmp_path / "o").iterdir())
    # Narrowing the pattern resolves it.
    res = _batch(data, in_dir, tmp_path / "o", "--pattern", "*.mmseqs2.tsv.gz")
    assert res.exit_code == 0, res.output
    assert (tmp_path / "o" / "s1_classifications.csv").exists()


def test_batch_workers_is_deprecated_noop(tmp_path, data):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    _write_alignment(in_dir / "s.blast.tsv.gz")
    res = _batch(data, in_dir, tmp_path / "o", "--workers", "4")
    assert res.exit_code == 0
    assert "--workers is deprecated" in res.output


def test_batch_rejects_streaming_with_id_mapping(tmp_path, data):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    _write_alignment(in_dir / "s.blast.tsv.gz")
    mapping = tmp_path / "map.tsv"
    mapping.write_text("contig_id\taccession\nc\tGCF_000000001.1\n")
    res = _batch(data, in_dir, tmp_path / "o", "--streaming", "--id-mapping", str(mapping))
    assert res.exit_code != 0
