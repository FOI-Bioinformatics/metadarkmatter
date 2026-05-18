"""Tests for FASTA/FASTQ validators and their CLI surface."""

from __future__ import annotations

import gzip
from pathlib import Path

from typer.testing import CliRunner

from metadarkmatter.cli.main import app
from metadarkmatter.core.sequence_validation import (
    validate_fasta,
    validate_fastq,
)


def test_validate_fasta_accepts_well_formed(tmp_path: Path) -> None:
    p = tmp_path / "ok.fasta"
    p.write_text(">read_1 first\nACGTACGT\nACGTAC\n>read_2\nNNNGT\n")
    r = validate_fasta(p)
    assert r.ok, r.issues
    assert r.record_count == 2
    assert r.sample_record_id == "read_1"


def test_validate_fasta_flags_invalid_characters(tmp_path: Path) -> None:
    p = tmp_path / "bad.fasta"
    p.write_text(">r1\nACGTAC123\n")
    r = validate_fasta(p)
    assert not r.ok
    assert any("invalid nucleotide" in m for m in r.issues)


def test_validate_fasta_flags_empty_record(tmp_path: Path) -> None:
    p = tmp_path / "empty.fasta"
    p.write_text(">r1\n>r2\nACGT\n")
    r = validate_fasta(p)
    assert not r.ok
    assert any("no sequence body" in m for m in r.issues)


def test_validate_fasta_reads_gzip(tmp_path: Path) -> None:
    p = tmp_path / "ok.fasta.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(">r1\nACGT\n")
    r = validate_fasta(p)
    assert r.ok
    assert r.record_count == 1


def test_validate_fasta_protein_mode(tmp_path: Path) -> None:
    p = tmp_path / "p.faa"
    p.write_text(">r1\nMKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWE\n")
    r = validate_fasta(p, sequence_type="protein")
    assert r.ok, r.issues


def test_validate_fastq_accepts_well_formed(tmp_path: Path) -> None:
    p = tmp_path / "ok.fastq"
    p.write_text(
        "@r1\nACGT\n+\nIIII\n"
        "@r2 extra\nAACC\n+\nJJJJ\n"
    )
    r = validate_fastq(p)
    assert r.ok, r.issues
    assert r.record_count == 2
    assert r.sample_record_id == "r1"


def test_validate_fastq_flags_length_mismatch(tmp_path: Path) -> None:
    p = tmp_path / "bad.fastq"
    p.write_text("@r1\nACGT\n+\nII\n")
    r = validate_fastq(p)
    assert not r.ok
    assert any("sequence length" in m for m in r.issues)


def test_validate_fastq_flags_truncated_record(tmp_path: Path) -> None:
    p = tmp_path / "trunc.fastq"
    p.write_text("@r1\nACGT\n+\nIIII\n@r2\nAA\n")
    r = validate_fastq(p)
    assert not r.ok
    assert any("Trailing partial record" in m for m in r.issues)


def test_validate_fastq_flags_missing_header(tmp_path: Path) -> None:
    p = tmp_path / "bad.fastq"
    p.write_text("r1\nACGT\n+\nIIII\n")
    r = validate_fastq(p)
    assert not r.ok
    assert any("must start with '@'" in m for m in r.issues)


def test_cli_validate_fasta_happy(tmp_path: Path) -> None:
    runner = CliRunner()
    p = tmp_path / "ok.fa"
    p.write_text(">r1\nACGT\n")
    res = runner.invoke(app, ["validate", "fasta", "--input", str(p)])
    assert res.exit_code == 0, res.output
    assert "valid" in res.output


def test_cli_validate_fastq_reports_errors(tmp_path: Path) -> None:
    runner = CliRunner()
    p = tmp_path / "bad.fq"
    p.write_text("@r1\nACGT\n+\nII\n")
    res = runner.invoke(app, ["validate", "fastq", "--input", str(p)])
    assert res.exit_code == 1
    assert "sequence length" in res.output
