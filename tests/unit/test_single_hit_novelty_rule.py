"""Tests for the single-hit novelty Stage 2 rule (``single_hit_novel_species``)."""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.main import app
from metadarkmatter.core.classification.ani_matrix import ANIMatrix
from metadarkmatter.core.classification.bayesian import (
    apply_single_hit_novelty_rule,
    single_hit_novel_species_mask,
)
from metadarkmatter.core.classification.classifiers.vectorized import (
    VectorizedClassifier,
)
from metadarkmatter.models.config import ScoringConfig

# Single-hit reads at several novelty levels plus one multi-hit read.
ROWS = [
    ("n2", "A|c", 98.0, 260),   # N = 2: Known Species either way
    ("n45", "A|c", 95.5, 240),  # N = 4.5: Bayesian Known, rule -> Novel Species
    ("n55", "A|c", 94.5, 235),  # N = 5.5: same
    ("n8", "A|c", 92.0, 220),   # N = 8: already Novel Species
    ("multi", "A|c", 95.5, 240),
    ("multi", "B|c", 95.0, 238),  # competitor -> not unambiguous
]


def _ani() -> ANIMatrix:
    return ANIMatrix({"A": {"A": 100.0, "B": 99.0}, "B": {"A": 99.0, "B": 100.0}})


def _write(path: Path) -> Path:
    path.write_text(
        "".join(
            f"{q}\t{s}\t{p}\t150\t0\t0\t1\t150\t1\t150\t1e-40\t{b}\t150\n"
            for q, s, p, b in ROWS
        )
    )
    return path


def _calls(df: pl.DataFrame) -> dict[str, str]:
    return dict(zip(df["read_id"], df["taxonomic_call"], strict=True))


class TestMask:
    def _sel(self, call="Known Species", scope="unambiguous", n=4.5, tmin=4.0) -> bool:
        df = pl.DataFrame({"taxonomic_call": [call], "ambiguity_scope": [scope], "novelty_index": [n]})
        return bool(df.select(single_hit_novel_species_mask(tmin))[0, 0])

    def test_selected(self):
        assert self._sel()

    @pytest.mark.parametrize(
        "kw",
        [{"n": 3.9}, {"scope": "within_genus"}, {"call": "Ambiguous"}, {"call": "Novel Genus"}],
    )
    def test_not_selected(self, kw):
        assert not self._sel(**kw)

    def test_apply_keeps_posteriors(self):
        df = pl.DataFrame({
            "taxonomic_call": ["Known Species"], "ambiguity_scope": ["unambiguous"],
            "novelty_index": [5.0], "p_known_species": [0.6],
        })
        row = apply_single_hit_novelty_rule(df, 4.0).row(0, named=True)
        assert row["taxonomic_call"] == "Novel Species"
        assert row["single_hit_novel_species"] is True
        assert row["p_known_species"] == 0.6


class TestClassifier:
    def test_on_by_default(self, tmp_path):
        out = VectorizedClassifier(_ani(), config=ScoringConfig()).classify_file(
            _write(tmp_path / "a.tsv")
        )
        calls = _calls(out)
        assert calls["n2"] == "Known Species"
        assert calls["n45"] == "Novel Species"
        assert calls["n55"] == "Novel Species"
        assert calls["n8"] == "Novel Species"
        flags = dict(zip(out["read_id"], out["single_hit_novel_species"], strict=True))
        assert flags["n45"] and flags["n55"]
        assert not flags["n2"] and not flags["n8"]
        assert not flags["multi"]  # has a competing reference: rule never applies
        row = out.filter(pl.col("read_id") == "n45").row(0, named=True)
        assert row["diversity_status"] == "Novel" and row["is_novel"] is True

    def test_disabled_restores_previous(self, tmp_path):
        out = VectorizedClassifier(
            _ani(), config=ScoringConfig(single_hit_novel_species=False)
        ).classify_file(_write(tmp_path / "a.tsv"))
        assert "single_hit_novel_species" not in out.columns
        assert _calls(out)["n45"] == "Known Species"

    def test_threshold_follows_novelty_known_max(self, tmp_path):
        cfg = ScoringConfig(novelty_known_max=5.0, novelty_novel_species_min=5.0)
        out = VectorizedClassifier(_ani(), config=cfg).classify_file(_write(tmp_path / "a.tsv"))
        calls = _calls(out)
        assert calls["n45"] == "Known Species"
        assert calls["n55"] == "Novel Species"

    def test_streaming_matches_in_memory(self, tmp_path):
        path = _write(tmp_path / "a.tsv")
        clf = VectorizedClassifier(_ani(), config=ScoringConfig())
        mem = clf.classify_file(path).sort("read_id")
        out = tmp_path / "s.parquet"
        clf.stream_to_file(path, out, output_format="parquet", partition_size=2)
        assert pl.read_parquet(out).sort("read_id").equals(mem)

    def test_not_applied_in_protein_mode(self, tmp_path):
        out = VectorizedClassifier(
            _ani(), config=ScoringConfig(alignment_mode="protein")
        ).classify_file(_write(tmp_path / "a.tsv"))
        assert "single_hit_novel_species" not in out.columns


class TestConfigAndCli:
    def test_yaml_round_trip(self, tmp_path):
        path = tmp_path / "c.yaml"
        ScoringConfig(single_hit_novel_species=False).to_yaml(path)
        assert "single_hit_novelty" in path.read_text()
        assert ScoringConfig.from_yaml(path).single_hit_novel_species is False

    @pytest.mark.parametrize(
        ("flag", "expected"),
        [([], "Novel Species"), (["--no-single-hit-novel-species"], "Known Species")],
    )
    def test_cli_flag(self, tmp_path, flag, expected):
        aln = _write(tmp_path / "a.tsv")
        ani = tmp_path / "ani.csv"
        ani.write_text("genome,A,B\nA,100,99\nB,99,100\n")
        out = tmp_path / "o.csv"
        result = CliRunner().invoke(
            app,
            ["score", "classify", "--alignment", str(aln), "--ani", str(ani),
             "--output", str(out), "-q", *flag],
        )
        assert result.exit_code == 0, result.output
        assert _calls(pl.read_csv(out))["n45"] == expected
