"""Tests for --max-read-pairs / --seed on 'kraken2 extract'."""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from unittest.mock import patch

from typer.testing import CliRunner

from metadarkmatter.cli.main import app
from metadarkmatter.external.base import ToolResult

N_RECORDS = 20


def _fastq(prefix: str, n: int = N_RECORDS) -> str:
    return "".join(f"@{prefix}{i}\nACGT\n+\nIIII\n" for i in range(n))


def _fake_run(**kwargs):
    for key, prefix in (("output_1", "a"), ("output_2", "b")):
        p = kwargs.get(key)
        if p is not None:
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(_fastq(prefix))
    return ToolResult(
        command=("mock",), return_code=0, stdout="", stderr="", elapsed_seconds=0.1
    )


def _invoke(tmp_path: Path, extra: list[str], env: dict[str, str] | None = None):
    kraken = tmp_path / "s.kraken"
    kraken.write_text("C\tr\t262\t150\t262:150")
    (tmp_path / "s.kreport").write_text("100.00\t1000\t0\tR\t1\troot\n")
    reads = tmp_path / "r1.fastq"
    reads.write_text("@r\nA\n+\nI\n")
    reads2 = tmp_path / "r2.fastq"
    reads2.write_text("@r\nA\n+\nI\n")
    out = tmp_path / "ext"
    with (
        patch("metadarkmatter.cli.kraken2.ExtractKrakenReads.check_available",
              return_value=True),
        patch("metadarkmatter.cli.kraken2.ExtractKrakenReads.run_or_raise",
              side_effect=_fake_run),
        patch("metadarkmatter.cli.kraken2.ExtractKrakenReads.run",
              return_value=ToolResult(command=("mock",), return_code=0, stdout="",
                                      stderr="", elapsed_seconds=0.0)),
        patch("metadarkmatter.cli.kraken2.KrakenReport"),
    ):
        result = CliRunner().invoke(
            app,
            ["kraken2", "extract", "--kraken-output", str(kraken),
             "--reads-1", str(reads), "--reads-2", str(reads2),
             "--taxid", "262", "--output", str(out), *extra],
            env=env,
        )
    return result, out


def _count(path: Path) -> int:
    with gzip.open(path, "rt") as fh:
        return sum(1 for _ in fh) // 4


def test_extract_max_read_pairs_subsamples_and_writes_sidecar(tmp_path):
    result, out = _invoke(tmp_path, ["--max-read-pairs", "5", "--seed", "3"])
    assert result.exit_code == 0, result.output
    assert _count(out / "s_taxid262_R1.fastq.gz") == 5
    assert _count(out / "s_taxid262_R2.fastq.gz") == 5
    side = json.loads((out / "s_taxid262_subsample.json").read_text())
    assert side == {"total": N_RECORDS, "kept": 5, "max_read_pairs": 5, "seed": 3}
    assert "Subsampled 5 of 20 read pairs (seed 3)" in " ".join(result.output.split())


def test_extract_quiet_suppresses_subsample_line(tmp_path):
    result, out = _invoke(tmp_path, ["--max-read-pairs", "5", "--quiet"])
    assert "Subsampled" not in result.output
    assert (out / "s_taxid262_subsample.json").exists()


def test_extract_without_max_read_pairs_unchanged(tmp_path):
    result, out = _invoke(tmp_path, [])
    assert result.exit_code == 0, result.output
    assert _count(out / "s_taxid262_R1.fastq.gz") == N_RECORDS
    assert not (out / "s_taxid262_subsample.json").exists()


def test_extract_skip_if_exists_does_not_resubsample(tmp_path):
    out = tmp_path / "ext"
    out.mkdir()
    r1 = out / "s_taxid262_R1.fastq.gz"
    r2 = out / "s_taxid262_R2.fastq.gz"
    for p in (r1, r2):
        with gzip.open(p, "wt") as fh:
            fh.write(_fastq("x"))
    side = out / "s_taxid262_subsample.json"
    side.write_text('{"sentinel": 1}')
    before = r1.read_bytes()
    result, _ = _invoke(tmp_path, ["--max-read-pairs", "5", "--skip-if-exists"])
    assert result.exit_code == 0, result.output
    assert r1.read_bytes() == before
    assert side.read_text() == '{"sentinel": 1}'


def test_extract_seed_default_from_env(tmp_path):
    result, out = _invoke(tmp_path, ["--max-read-pairs", "5"],
                          env={"METADARKMATTER_SEED": "7"})
    assert result.exit_code == 0, result.output
    side = json.loads((out / "s_taxid262_subsample.json").read_text())
    assert side["seed"] == 7


def test_extract_dry_run_does_no_subsampling(tmp_path):
    result, out = _invoke(tmp_path, ["--max-read-pairs", "5", "--dry-run"])
    assert result.exit_code == 0, result.output
    assert not (out / "s_taxid262_subsample.json").exists()


def test_extract_max_read_pairs_rejects_zero(tmp_path):
    result, _ = _invoke(tmp_path, ["--max-read-pairs", "0"])
    assert result.exit_code != 0
