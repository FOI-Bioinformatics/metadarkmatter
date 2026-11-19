"""Tests for core.subsample."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from metadarkmatter.core import subsample as sub
from metadarkmatter.core.paired_reads import open_reads
from metadarkmatter.core.subsample import SubsampleResult, subsample_reads


def _write(path: Path, n: int, mate: int, quality_at: bool = False) -> None:
    lines = []
    for i in range(n):
        qual = "@IIII" if quality_at else "IIII"
        lines += [f"@read{i}/{mate} comment", "ACGT", "+", qual]
    data = "\n".join(lines) + "\n"
    if path.suffix == ".gz":
        with gzip.open(path, "wt") as fh:
            fh.write(data)
    else:
        path.write_text(data)


def _records(path: Path) -> list[list[str]]:
    with open_reads(path) as fh:
        lines = fh.read().splitlines()
    return [lines[i : i + 4] for i in range(0, len(lines), 4)]


def _ids(path: Path) -> list[str]:
    return [r[0].split()[0][1:].rsplit("/", 1)[0] for r in _records(path)]


def _pair(tmp_path: Path, n: int = 100, suffix: str = ".fastq") -> tuple[Path, Path]:
    r1, r2 = tmp_path / f"a_1{suffix}", tmp_path / f"a_2{suffix}"
    _write(r1, n, 1)
    _write(r2, n, 2)
    return r1, r2


def test_subsample_keeps_n_pairs_and_matches_mates(tmp_path):
    r1, r2 = _pair(tmp_path)
    res = subsample_reads(r1, r2, 10, 42)
    assert res == SubsampleResult(total=100, kept=10, seed=42)
    assert len(_records(r1)) == 10 and len(_records(r2)) == 10
    assert _ids(r1) == _ids(r2)


def test_subsample_is_reproducible(tmp_path):
    dirs = [tmp_path / n for n in "abc"]
    for d in dirs:
        d.mkdir()
    pairs = [_pair(d) for d in dirs]
    subsample_reads(*pairs[0], 10, 7)
    subsample_reads(*pairs[1], 10, 7)
    subsample_reads(*pairs[2], 10, 8)
    assert pairs[0][0].read_bytes() == pairs[1][0].read_bytes()
    assert pairs[0][1].read_bytes() == pairs[1][1].read_bytes()
    assert _ids(pairs[0][0]) != _ids(pairs[2][0])


def test_subsample_preserves_input_order(tmp_path):
    r1, r2 = _pair(tmp_path)
    subsample_reads(r1, r2, 20, 1)
    nums = [int(i[4:]) for i in _ids(r1)]
    assert nums == sorted(nums)


def test_subsample_no_op_when_fewer_reads(tmp_path):
    r1, r2 = _pair(tmp_path, n=5)
    before = (r1.read_bytes(), r2.read_bytes())
    mtime = r1.stat().st_mtime_ns
    res = subsample_reads(r1, r2, 5, 1)
    assert res.total == 5 and res.kept == 5
    assert (r1.read_bytes(), r2.read_bytes()) == before
    assert r1.stat().st_mtime_ns == mtime
    res = subsample_reads(r1, r2, 50, 1)
    assert res.kept == 5


def test_subsample_single_end(tmp_path):
    r1 = tmp_path / "s.fastq"
    _write(r1, 30, 1)
    res = subsample_reads(r1, None, 4, 3)
    assert res.kept == 4
    assert len(_records(r1)) == 4


def test_subsample_gzip_and_plain(tmp_path):
    gz1, gz2 = _pair(tmp_path, suffix=".fastq.gz")
    subsample_reads(gz1, gz2, 10, 5)
    assert gz1.read_bytes()[:2] == b"\x1f\x8b"
    assert len(_records(gz1)) == 10
    assert _ids(gz1) == _ids(gz2)
    plain_dir = tmp_path / "plain"
    plain_dir.mkdir()
    p1, p2 = _pair(plain_dir)
    subsample_reads(p1, p2, 10, 5)
    assert p1.read_bytes()[:1] == b"@"
    assert _ids(p1) == _ids(gz1)


def test_subsample_quality_line_starting_with_at(tmp_path):
    r1 = tmp_path / "q.fastq"
    _write(r1, 20, 1, quality_at=True)
    res = subsample_reads(r1, None, 7, 2)
    assert res.total == 20 and res.kept == 7
    for rec in _records(r1):
        assert rec[0].startswith("@read") and rec[3] == "@IIII"


def test_subsample_count_mismatch_raises(tmp_path):
    r1 = tmp_path / "m_1.fastq"
    r2 = tmp_path / "m_2.fastq"
    _write(r1, 10, 1)
    _write(r2, 9, 2)
    before = (r1.read_bytes(), r2.read_bytes())
    with pytest.raises(ValueError, match="mismatch"):
        subsample_reads(r1, r2, 3, 1)
    assert (r1.read_bytes(), r2.read_bytes()) == before


def test_subsample_truncated_record_raises(tmp_path):
    r1 = tmp_path / "t.fastq"
    r1.write_text("@a\nACGT\n+\nIIII\n@b\nAC\n")
    with pytest.raises(ValueError, match="truncated"):
        subsample_reads(r1, None, 1, 1)


def test_subsample_failure_leaves_originals(tmp_path, monkeypatch):
    r1, r2 = _pair(tmp_path)
    before = (r1.read_bytes(), r2.read_bytes())
    real = sub._write_subset
    calls = {"n": 0}

    def flaky(*args, **kwargs):
        calls["n"] += 1
        if calls["n"] == 2:
            raise OSError("disk full")
        return real(*args, **kwargs)

    monkeypatch.setattr(sub, "_write_subset", flaky)
    with pytest.raises(OSError):
        subsample_reads(r1, r2, 10, 1)
    assert (r1.read_bytes(), r2.read_bytes()) == before
    assert sorted(p.name for p in tmp_path.iterdir()) == sorted([r1.name, r2.name])
