"""Tests for streaming read export in 'score extract-novel'."""

from __future__ import annotations

import gzip
from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.main import app
from metadarkmatter.core.exceptions import ConfigurationError
from metadarkmatter.core.read_export import (
    extract_reads,
    iter_records,
    mate_output_path,
    normalize_read_id,
)

runner = CliRunner()


def _fastq(path: Path, names: list[str]) -> Path:
    text = "".join(f"@{n} 1:N:0\nACGTACGT\n+\nIIIIIIII\n" for n in names)
    if path.name.endswith(".gz"):
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return path


def _read_names(path: Path) -> list[str]:
    return [r.header.decode().split()[0] for r in iter_records(path)]


def _classifications(path: Path, rows: list[tuple[str, str, str]]) -> Path:
    """rows: (read_id, genome, call)."""
    pl.DataFrame({
        "read_id": [r[0] for r in rows],
        "best_match_genome": [r[1] for r in rows],
        "taxonomic_call": [r[2] for r in rows],
        "novelty_index": [8.0] * len(rows),
        "placement_uncertainty": [0.5] * len(rows),
        "top_hit_identity": [92.0] * len(rows),
        "num_ambiguous_hits": [1] * len(rows),
        "identity_gap": [None] * len(rows),
    }, schema_overrides={"identity_gap": pl.Float64}).write_csv(path)
    return path


# ------------------------------------------------------------------ unit


def test_normalize_read_id():
    assert normalize_read_id("r1/1") == "r1"
    assert normalize_read_id("r1/2") == "r1"
    assert normalize_read_id("@r1/1 1:N:0:ACGT") == "r1"
    assert normalize_read_id(">r1 desc") == "r1"
    assert normalize_read_id("r1") == "r1"
    assert normalize_read_id("a/b/1") == "a/b"


def test_mate_output_path():
    assert mate_output_path(Path("d/novel.fastq.gz"), 1) == Path("d/novel_R1.fastq.gz")
    assert mate_output_path(Path("novel.fa"), 2) == Path("novel_R2.fa")
    assert mate_output_path(Path("novel.fastq"), None) == Path("novel.fastq")
    assert mate_output_path(Path("odd.out.gz"), 1) == Path("odd.out_R1.gz")


def test_iter_records_fasta_multiline_and_gzip(tmp_path):
    p = tmp_path / "x.fa.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(">a desc\nACGT\nTTTT\n>b\nGG\n")
    recs = list(iter_records(p))
    assert [(r.read_id, r.seq) for r in recs] == [("a", b"ACGTTTTT"), ("b", b"GG")]


def test_extract_reads_single_end_suffix_normalisation(tmp_path):
    r1 = _fastq(tmp_path / "r.fastq.gz", ["a", "b", "c", "d"])
    rep = extract_reads(
        reads_1=r1, ids={"a", "c", "zzz"}, fastq_out=tmp_path / "o" / "n.fastq.gz"
    )
    assert _read_names(tmp_path / "o" / "n.fastq.gz") == ["a", "c"]
    assert (rep.requested, rep.found_1, rep.written_1) == (3, 2, 2)
    assert rep.missing == ["zzz"]


def test_extract_reads_paired_mate_suffix_in_one_side_only(tmp_path):
    # Classification IDs carry /1 (normalised to bare); R1 FASTQ has no suffix,
    # R2 FASTQ has /2. Either mate being selected pulls the whole pair.
    r1 = _fastq(tmp_path / "s_R1.fastq", ["p1", "p2", "p3"])
    r2 = _fastq(tmp_path / "s_R2.fastq", ["p1/2", "p2/2", "p3/2"])
    rep = extract_reads(
        reads_1=r1, reads_2=r2, ids={"p1", "p3"}, fastq_out=tmp_path / "novel.fastq"
    )
    assert _read_names(tmp_path / "novel_R1.fastq") == ["p1", "p3"]
    assert _read_names(tmp_path / "novel_R2.fastq") == ["p1/2", "p3/2"]
    assert (rep.found_1, rep.found_2) == (2, 2)
    assert not (tmp_path / "novel.fastq").exists()


def test_fastq_to_fasta_conversion_and_gz(tmp_path):
    r1 = _fastq(tmp_path / "r.fastq", ["a", "b"])
    extract_reads(reads_1=r1, ids={"b"}, fasta_out=tmp_path / "n.fasta.gz")
    with gzip.open(tmp_path / "n.fasta.gz", "rt") as fh:
        assert fh.read() == ">b 1:N:0\nACGTACGT\n"


def test_fastq_out_from_fasta_input_is_an_error(tmp_path):
    fa = tmp_path / "r.fa"
    fa.write_text(">a\nACGT\n")
    with pytest.raises(ConfigurationError):
        extract_reads(reads_1=fa, ids={"a"}, fastq_out=tmp_path / "o.fastq")


def test_per_lineage_many_lineages_exceed_open_handle_cap(tmp_path):
    n = 150  # more than the 64 simultaneously open handles
    names = [f"r{i}" for i in range(n)]
    # interleave lineages so evicted handles are re-opened in append mode
    lineage_of = {name: f"G{i % 100}" for i, name in enumerate(names)}
    r1 = _fastq(tmp_path / "r.fastq.gz", names)
    rep = extract_reads(
        reads_1=r1, ids=set(names), per_lineage_dir=tmp_path / "lin", lineage_of=lineage_of
    )
    assert rep.found_1 == n
    total = 0
    for f in (tmp_path / "lin").glob("*.fastq.gz"):
        total += len(_read_names(f))
    assert total == n
    assert len(list((tmp_path / "lin").glob("*.fastq.gz"))) == 100


# ------------------------------------------------------------------ CLI


@pytest.fixture
def cls_and_reads(tmp_path):
    rows = [
        ("p1/1", "GCF_A", "Novel Species"),
        ("p1/2", "GCF_A", "Novel Species"),
        ("p2/1", "GCF_B", "Novel Genus"),
        ("p3/1", "GCF_A", "Known Species"),  # not novel -> not exported
        ("ghost/1", "GCF_B", "Novel Species"),  # absent from the FASTQ
    ]
    c = _classifications(tmp_path / "c.csv", rows)
    r1 = _fastq(tmp_path / "s_R1.fastq.gz", ["p1/1", "p2/1", "p3/1", "p4/1"])
    r2 = _fastq(tmp_path / "s_R2.fastq.gz", ["p1/2", "p2/2", "p3/2", "p4/2"])
    return c, r1, r2


def test_cli_extract_novel_paired_fastq_per_lineage(tmp_path, cls_and_reads):
    c, r1, r2 = cls_and_reads
    res = runner.invoke(app, [
        "score", "extract-novel", "-c", str(c), "-o", str(tmp_path / "n.csv"),
        "--reads-1", str(r1), "--reads-2", str(r2),
        "--fastq-out", str(tmp_path / "out" / "novel.fastq.gz"),
        "--per-lineage-dir", str(tmp_path / "lin"),
    ])
    assert res.exit_code == 0, res.output
    o1 = tmp_path / "out" / "novel_R1.fastq.gz"
    o2 = tmp_path / "out" / "novel_R2.fastq.gz"
    assert _read_names(o1) == ["p1/1", "p2/1"]
    assert _read_names(o2) == ["p1/2", "p2/2"]
    # requested = p1, p2, ghost; found 2
    assert "Requested read IDs:  3" in res.output
    assert "Found in reads-1:    2" in res.output
    assert "not found in reads-1" in res.output and "ghost" in res.output
    lin = {f.name: f for f in (tmp_path / "lin").iterdir()}
    assert sorted(lin) == [
        "GCF_A_R1.fastq.gz", "GCF_A_R2.fastq.gz", "GCF_B_R1.fastq.gz", "GCF_B_R2.fastq.gz",
    ]
    assert _read_names(lin["GCF_A_R1.fastq.gz"]) == ["p1/1"]
    assert _read_names(lin["GCF_B_R2.fastq.gz"]) == ["p2/2"]


def test_cli_extract_novel_single_end_fasta_out_plain(tmp_path, cls_and_reads):
    c, r1, _ = cls_and_reads
    out = tmp_path / "n.fasta"
    res = runner.invoke(app, [
        "score", "extract-novel", "-c", str(c), "-o", str(tmp_path / "n.csv"),
        "--reads-1", str(r1), "--fasta-out", str(out), "--category", "species",
    ])
    assert res.exit_code == 0, res.output
    assert out.read_text().count(">") == 1  # only p1 (p2 is Novel Genus)


def test_cli_extract_novel_respects_min_reads_lineage_set(tmp_path, cls_and_reads):
    c, r1, _ = cls_and_reads
    out = tmp_path / "n.fastq"
    res = runner.invoke(app, [
        "score", "extract-novel", "-c", str(c), "-o", str(tmp_path / "n.csv"),
        "--reads-1", str(r1), "--fastq-out", str(out), "--min-reads", "2",
    ])
    assert res.exit_code == 0, res.output
    # GCF_A has 2 mate rows (p1/1, p1/2), GCF_B has 2 rows (p2/1, ghost/1); both pass.
    assert _read_names(out) == ["p1/1", "p2/1"]


@pytest.mark.parametrize(
    "args",
    [
        ["--reads-2", "R2"],
        ["--fastq-out", "OUT"],
        ["--reads-1", "R1"],
    ],
)
def test_cli_extract_novel_option_validation(tmp_path, cls_and_reads, args):
    c, r1, r2 = cls_and_reads
    args = [{"R1": str(r1), "R2": str(r2), "OUT": str(tmp_path / "o.fastq")}.get(a, a)
            for a in args]
    res = runner.invoke(app, [
        "score", "extract-novel", "-c", str(c), "-o", str(tmp_path / "n.csv"), *args,
    ])
    assert res.exit_code == 1
    assert "Error" in res.output


def test_cli_extract_novel_without_reads_is_unchanged(tmp_path, cls_and_reads):
    c, _, _ = cls_and_reads
    res = runner.invoke(app, [
        "score", "extract-novel", "-c", str(c), "-o", str(tmp_path / "n.csv"),
        "-r", str(tmp_path / "ids.txt"),
    ])
    assert res.exit_code == 0, res.output
    assert (tmp_path / "ids.txt").exists()
