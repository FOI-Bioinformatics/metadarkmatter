"""Paired-end ``blast align``: argument handling, mate IDs, and a real-blastn check."""

from __future__ import annotations

import gzip
import shutil
import subprocess
from pathlib import Path

import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.main import app
from metadarkmatter.core.paired_reads import combine_mates
from metadarkmatter.external.base import ExternalTool

runner = CliRunner()


def _fastq(path: Path, records: list[tuple[str, str]], gz: bool = False) -> Path:
    text = "".join(f"@{rid}\n{seq}\n+\n{'I' * len(seq)}\n" for rid, seq in records)
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return path


def test_combine_mates_suffixes_and_orders(tmp_path: Path):
    r1 = _fastq(tmp_path / "r1.fastq", [("a", "ACGT"), ("b/1", "ACGA")])
    r2 = _fastq(tmp_path / "r2.fastq.gz", [("a", "TTTT"), ("b/2", "GGGG")], gz=True)
    out = tmp_path / "all.fastq"
    assert combine_mates(r1, r2, out) == 4
    ids = [ln[1:].strip() for ln in out.read_text().splitlines()[::4]]
    assert ids == ["a/1", "b/1", "a/2", "b/2"]


def _align(tmp_path: Path, *args: str):
    db = tmp_path / "db"
    return runner.invoke(
        app,
        ["blast", "align", "--database", str(db), "--output", str(tmp_path / "o.tsv.gz"), *args],
    )


def test_align_requires_an_input_mode(tmp_path: Path):
    res = _align(tmp_path, "--dry-run")
    assert res.exit_code == 1
    assert "--query" in res.output


def test_align_rejects_mixed_and_partial_pairs(tmp_path: Path):
    r1 = _fastq(tmp_path / "r1.fastq", [("a", "ACGT")])
    r2 = _fastq(tmp_path / "r2.fastq", [("a", "ACGT")])
    mixed = _align(tmp_path, "--query", str(r1), "--query-1", str(r1), "--query-2", str(r2))
    assert mixed.exit_code == 1 and "Choose one input mode" in mixed.output
    partial = _align(tmp_path, "--query-1", str(r1))
    assert partial.exit_code == 1 and "requires both" in partial.output


def test_align_rejects_mixed_formats(tmp_path: Path):
    r1 = _fastq(tmp_path / "r1.fastq", [("a", "ACGT")])
    r2 = tmp_path / "r2.fasta"
    r2.write_text(">a\nACGT\n")
    res = _align(tmp_path, "--query-1", str(r1), "--query-2", str(r2))
    assert res.exit_code == 1 and "same format" in res.output


def test_align_paired_dry_run_describes_combination(tmp_path: Path):
    r1 = _fastq(tmp_path / "r1.fastq", [("a", "ACGT")])
    r2 = _fastq(tmp_path / "r2.fastq", [("a", "ACGT")])
    res = _align(tmp_path, "--query-1", str(r1), "--query-2", str(r2), "--dry-run")
    assert res.exit_code == 0, res.output
    text = " ".join(res.output.split())
    assert "Dry run complete" in text
    assert "combined into one temp query" in text or "not found; required" in text


@pytest.mark.skipif(
    shutil.which("blastn") is None or shutil.which("makeblastdb") is None,
    reason="BLAST+ not installed",
)
def test_paired_align_with_real_blastn(tmp_path: Path):
    import random

    # Other tests in this worker may have mocked the resolver or cached a
    # "not found" result; this test needs the real executables.
    ExternalTool.reset_executable_resolver()
    ExternalTool.clear_cache()

    rng = random.Random(7)
    genome = "".join(rng.choice("ACGT") for _ in range(3000))
    fna = tmp_path / "GCF_000000001.1.fna"
    fna.write_text(f">c1\n{genome}\n")
    db = tmp_path / "db" / "pan"
    db.parent.mkdir()
    res = runner.invoke(app, ["blast", "makedb", "--genomes", str(tmp_path), "--output", str(db), "-q"])
    assert res.exit_code == 0, res.output
    assert subprocess.run(["blastn", "-version"], capture_output=True).returncode == 0

    # Mates share one ID, as in SRA-style data.
    r1_records, r2_records = [], []
    for i in range(5):
        start = 100 + i * 300
        r1_records.append((f"read{i}", genome[start : start + 100]))
        r2_records.append((f"read{i}", genome[start + 150 : start + 250]))
    r1 = _fastq(tmp_path / "r1.fastq", r1_records)
    r2 = _fastq(tmp_path / "r2.fastq.gz", r2_records, gz=True)
    out = tmp_path / "paired.tsv.gz"
    res = runner.invoke(
        app,
        ["blast", "align", "--query-1", str(r1), "--query-2", str(r2), "--database", str(db),
         "--output", str(out), "--task", "megablast", "--threads", "1", "-q"],
    )
    assert res.exit_code == 0, res.output
    with gzip.open(out, "rt") as fh:
        qids = [ln.split("\t")[0] for ln in fh]
    distinct = list(dict.fromkeys(qids))
    assert sorted(distinct) == sorted(
        [f"read{i}/1" for i in range(5)] + [f"read{i}/2" for i in range(5)]
    )
    # Hits of each read are contiguous (required by streaming classify).
    assert len(distinct) == len(set(qids))
    assert [q for i, q in enumerate(qids) if i == 0 or q != qids[i - 1]] == distinct
