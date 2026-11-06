"""Regression tests for the round-3 audit fixes (extract-novel, read export, provenance)."""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.main import app
from metadarkmatter.core.exceptions import ConfigurationError
from metadarkmatter.core.read_export import _LineageWriters, extract_reads

from .test_provenance_batch import FILTERS, _write_alignment, _write_ani

runner = CliRunner()


def _novel_csv(path: Path, rows: list[tuple[str, str, str]]) -> Path:
    pl.DataFrame({
        "read_id": [r for r, _g, _c in rows],
        "best_match_genome": [g for _r, g, _c in rows],
        "top_hit_identity": [90.0] * len(rows),
        "novelty_index": [10.0] * len(rows),
        "placement_uncertainty": [0.5] * len(rows),
        "taxonomic_call": [c for _r, _g, c in rows],
    }).write_csv(path)
    return path


class TestExtractNovel:
    def test_header_only_input(self, tmp_path):
        cls = tmp_path / "c.csv"
        cls.write_text("read_id,best_match_genome,top_hit_identity,novelty_index,"
                       "placement_uncertainty,taxonomic_call,clade_novel_species\n")
        out, ids = tmp_path / "n.csv", tmp_path / "ids.txt"
        res = runner.invoke(app, ["score", "extract-novel", "-c", str(cls), "-o", str(out),
                                  "-r", str(ids), "--min-novelty", "5"])
        assert res.exit_code == 0, res.output
        assert out.exists() and ids.exists()

    def test_grouped_output_deterministic_on_ties(self, tmp_path):
        rows = []
        for g in range(30):
            rows += [(f"g{g}a", f"G{g:02d}", "Novel Species"), (f"g{g}b", f"G{g:02d}", "Novel Genus")]
        cls = _novel_csv(tmp_path / "tie.csv", rows)
        outs = []
        for seed in range(3):
            out = tmp_path / f"o{seed}.csv"
            res = runner.invoke(app, ["score", "extract-novel", "-c", str(cls), "-o", str(out),
                                      "-t", "all", "-q"], env={"PYTHONHASHSEED": str(seed)})
            assert res.exit_code == 0, res.output
            outs.append(out.read_bytes())
        assert outs[0] == outs[1] == outs[2]
        df = pl.read_csv(tmp_path / "o0.csv")
        assert df["best_match_genome"].to_list() == sorted(df["best_match_genome"].to_list())
        assert set(df["dominant_category"]) == {"Novel Genus"}  # alphabetical tie-break

    def test_min_reads_applies_without_grouping(self, tmp_path):
        cls = _novel_csv(tmp_path / "c.csv", [
            ("a", "G1", "Novel Species"), ("b", "G1", "Novel Species"), ("c", "G2", "Novel Species"),
        ])
        out = tmp_path / "n.csv"
        res = runner.invoke(app, ["score", "extract-novel", "-c", str(cls), "-o", str(out),
                                  "--no-group", "--min-reads", "2", "-q"])
        assert res.exit_code == 0, res.output
        assert sorted(pl.read_csv(out)["read_id"]) == ["a", "b"]


class TestReadExport:
    def test_refuses_to_overwrite_input(self, tmp_path):
        reads = tmp_path / "r.fq"
        reads.write_text("@a\nACGT\n+\nIIII\n")
        with pytest.raises(ConfigurationError, match="input read file"):
            extract_reads(reads_1=reads, ids={"a"}, fastq_out=reads)
        assert reads.read_text().startswith("@a")

    def test_refuses_paired_mate_output_equal_to_input(self, tmp_path):
        r1, r2 = tmp_path / "x_R1.fq", tmp_path / "x_R2.fq"
        for p in (r1, r2):
            p.write_text("@a\nACGT\n+\nIIII\n")
        with pytest.raises(ConfigurationError):
            extract_reads(reads_1=r1, reads_2=r2, ids={"a"}, fastq_out=tmp_path / "x.fq")

    def test_lineage_names_unique(self, tmp_path):
        w = _LineageWriters(tmp_path, ["a b", "a_b", "a_b_1"], "fastq", paired=False)
        assert len(set(w.paths.values())) == 3


class TestClassifyCompanions:
    def _run(self, tmp_path, *extra):
        aln = _write_alignment(tmp_path / "s.blast.tsv")
        ani = _write_ani(tmp_path / "ani.csv")
        return runner.invoke(app, ["score", "classify", "--alignment", str(aln), "--ani", str(ani),
                                   "--output", str(tmp_path / "c.csv"), *FILTERS, *extra])

    def test_stale_summary_removed_when_nothing_classified(self, tmp_path):
        summary = tmp_path / "s.json"
        assert self._run(tmp_path, "--summary", str(summary)).exit_code == 0
        assert summary.exists()
        res = self._run(tmp_path, "--summary", str(summary), "--min-alignment-length", "1000")
        assert res.exit_code == 0, res.output
        assert not summary.exists()
        assert "Summary written" not in res.output

    def test_stale_provenance_removed_with_no_provenance(self, tmp_path):
        assert self._run(tmp_path).exit_code == 0
        sidecar = tmp_path / "c.csv.provenance.json"
        assert sidecar.exists()
        assert self._run(tmp_path, "--no-provenance").exit_code == 0
        assert not sidecar.exists()

    def test_genomes_directory_recorded(self, tmp_path):
        gdir = tmp_path / "genomes"
        gdir.mkdir()
        for acc in ("GCF_000000001.1", "GCF_000000002.1"):
            (gdir / f"{acc}.fna").write_text(f">{acc}_contig1\nACGT\n")
        res = self._run(tmp_path, "--genomes", str(gdir))
        assert res.exit_code == 0, res.output
        prov = json.loads((tmp_path / "c.csv.provenance.json").read_text())
        rec = prov["inputs"]["genomes"]
        assert rec["files"] == 2 and len(rec["listing_sha256"]) == 64
