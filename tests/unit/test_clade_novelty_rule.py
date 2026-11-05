"""Tests for the opt-in clade-novelty Stage 2 rule (``clade_novel_species``)."""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.main import app
from metadarkmatter.core.classification.ani_matrix import ANIMatrix
from metadarkmatter.core.classification.bayesian import (
    apply_clade_novelty_rule,
    clade_novel_species_mask,
)
from metadarkmatter.core.classification.classifiers.vectorized import (
    VectorizedClassifier,
)
from metadarkmatter.models.config import ScoringConfig

# A/B are two references of one genus (92% ANI). r_clade is ~9% divergent
# from both (novel species, placement ambiguous between A and B); r_known is
# a conserved-region read nearly identical to both (genuinely ambiguous);
# r_single is a single-hit novel read (unaffected by the rule).
ROWS = [
    ("r_clade", "A|c", 91.0, 220),
    ("r_clade", "B|c", 90.6, 218),
    ("r_known", "A|c", 99.5, 280),
    ("r_known", "B|c", 99.0, 278),
    ("r_single", "A|c", 90.0, 200),
]


def _ani() -> ANIMatrix:
    return ANIMatrix({"A": {"A": 100.0, "B": 92.0}, "B": {"A": 92.0, "B": 100.0}})


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
    def _frame(self, **overrides):
        base = {
            "taxonomic_call": ["Ambiguous"],
            "ambiguity_scope": ["within_genus"],
            "novelty_index": [9.0],
            "second_hit_identity": [90.6],
        }
        base.update({k: [v] for k, v in overrides.items()})
        return pl.DataFrame(base)

    def _selected(self, **overrides) -> bool:
        df = self._frame(**overrides)
        return bool(df.select(clade_novel_species_mask(5.0, 18.6))[0, 0])

    def test_selects_divergent_within_genus(self):
        assert self._selected()

    def test_species_boundary_and_within_species_selected(self):
        assert self._selected(
            taxonomic_call="Species Boundary", ambiguity_scope="within_species"
        )

    @pytest.mark.parametrize(
        "overrides",
        [
            {"ambiguity_scope": "across_genera"},  # stays Conserved Region
            {"taxonomic_call": "Known Species"},
            {"taxonomic_call": "Conserved Region"},
            {"novelty_index": 4.9},  # not divergent enough from best hit
            {"second_hit_identity": 96.0},  # close to a competitor
            {"novelty_index": 18.6},  # genus-level: left to genus categories
            {"second_hit_identity": None},  # no competitor
        ],
    )
    def test_not_selected(self, overrides):
        assert not self._selected(**overrides)

    def test_apply_keeps_posteriors_and_scope(self):
        df = self._frame().with_columns(
            pl.lit(0.7).alias("p_ambiguous"), pl.lit(8.0).alias("placement_uncertainty")
        )
        out = apply_clade_novelty_rule(df, 5.0, 18.6)
        row = out.row(0, named=True)
        assert row["taxonomic_call"] == "Novel Species"
        assert row["clade_novel_species"] is True
        assert row["p_ambiguous"] == 0.7
        assert row["placement_uncertainty"] == 8.0
        assert row["ambiguity_scope"] == "within_genus"


class TestClassifier:
    def test_default_off_unchanged_schema(self, tmp_path):
        out = VectorizedClassifier(_ani(), config=ScoringConfig()).classify_file(
            _write(tmp_path / "a.tsv")
        )
        assert "clade_novel_species" not in out.columns
        assert _calls(out)["r_clade"] == "Ambiguous"

    def test_enabled(self, tmp_path):
        out = VectorizedClassifier(
            _ani(), config=ScoringConfig(clade_novel_species=True)
        ).classify_file(_write(tmp_path / "a.tsv"))
        calls = _calls(out)
        assert calls["r_clade"] == "Novel Species"
        assert calls["r_known"] == "Ambiguous"
        assert calls["r_single"] == "Novel Species"
        flags = dict(zip(out["read_id"], out["clade_novel_species"], strict=True))
        assert flags == {"r_clade": True, "r_known": False, "r_single": False}
        row = out.filter(pl.col("read_id") == "r_clade").row(0, named=True)
        assert row["diversity_status"] == "Novel"
        assert row["is_novel"] is True
        assert row["placement_uncertainty"] == pytest.approx(8.0)

    def test_streaming_matches_in_memory(self, tmp_path):
        path = _write(tmp_path / "a.tsv")
        clf = VectorizedClassifier(_ani(), config=ScoringConfig(clade_novel_species=True))
        mem = clf.classify_file(path).sort("read_id")
        out = tmp_path / "s.parquet"
        clf.stream_to_file(path, out, output_format="parquet", partition_size=2)
        streamed = pl.read_parquet(out).sort("read_id")
        assert streamed.columns == mem.columns
        assert streamed.equals(mem)

    def test_not_applied_in_protein_mode(self, tmp_path, caplog):
        cfg = ScoringConfig(clade_novel_species=True, alignment_mode="protein")
        with caplog.at_level("WARNING"):
            clf = VectorizedClassifier(_ani(), config=cfg)
        assert "only validated for nucleotide" in caplog.text
        out = clf.classify_file(_write(tmp_path / "a.tsv"))
        assert "clade_novel_species" not in out.columns


class TestConfig:
    def test_bounds_validated(self):
        with pytest.raises(ValueError, match="clade_novelty_min"):
            ScoringConfig(clade_novelty_min=20.0, clade_novelty_max=18.6)

    def test_yaml_round_trip(self, tmp_path):
        cfg = ScoringConfig(
            clade_novel_species=True, clade_novelty_min=6.0, clade_novelty_max=15.0
        )
        path = tmp_path / "c.yaml"
        cfg.to_yaml(path)
        assert "clade_novelty" in path.read_text()
        back = ScoringConfig.from_yaml(path)
        assert back.clade_novel_species is True
        assert back.clade_novelty_min == 6.0
        assert back.clade_novelty_max == 15.0

    def test_cli_flag(self, tmp_path):
        aln = _write(tmp_path / "a.tsv")
        ani = tmp_path / "ani.csv"
        ani.write_text("genome,A,B\nA,100,92\nB,92,100\n")
        out = tmp_path / "o.csv"
        result = CliRunner().invoke(
            app,
            ["score", "classify", "--alignment", str(aln), "--ani", str(ani),
             "--output", str(out), "--clade-novel-species", "-q"],
        )
        assert result.exit_code == 0, result.output
        df = pl.read_csv(out)
        assert "clade_novel_species" in df.columns
        assert _calls(df)["r_clade"] == "Novel Species"
