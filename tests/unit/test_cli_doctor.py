"""Tests for `mdm doctor` database/cache checks and install hints."""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from metadarkmatter.cli.doctor import _check_kraken_db
from metadarkmatter.cli.main import app

runner = CliRunner()


def test_check_kraken_db_missing_path(tmp_path: Path):
    ok, detail = _check_kraken_db(tmp_path / "nope")
    assert not ok
    assert "does not exist" in detail


def test_check_kraken_db_incomplete(tmp_path: Path):
    (tmp_path / "hash.k2d").write_text("x")  # only one of three index files
    ok, detail = _check_kraken_db(tmp_path)
    assert not ok
    assert "missing index files" in detail
    assert "opts.k2d" in detail and "taxo.k2d" in detail


def test_check_kraken_db_valid(tmp_path: Path):
    for f in ("hash.k2d", "opts.k2d", "taxo.k2d"):
        (tmp_path / f).write_text("x")
    ok, detail = _check_kraken_db(tmp_path)
    assert ok
    assert "valid Kraken2 DB" in detail


def test_doctor_runs_and_shows_db_section():
    result = runner.invoke(app, ["doctor"])
    assert result.exit_code == 0
    assert "Databases & caches" in result.output
    assert "GTDB cache" in result.output
    # Without --kraken-db, the Kraken2 DB row is informational.
    assert "not checked" in result.output


def test_doctor_validates_kraken_db(tmp_path: Path):
    for f in ("hash.k2d", "opts.k2d", "taxo.k2d"):
        (tmp_path / f).write_text("x")
    result = runner.invoke(app, ["doctor", "--kraken-db", str(tmp_path)])
    assert result.exit_code == 0
    assert "valid Kraken2 DB" in result.output


def test_doctor_gtdb_cache_disabled(monkeypatch):
    monkeypatch.setenv("METADARKMATTER_GTDB_CACHE_DIR", "")
    result = runner.invoke(app, ["doctor"])
    assert result.exit_code == 0
    assert "disabled" in result.output
