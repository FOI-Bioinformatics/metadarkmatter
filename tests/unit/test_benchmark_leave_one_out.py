"""Unit tests for scripts/benchmark_leave_one_out.py (pure parts only)."""

from __future__ import annotations

import importlib.util
import random
from pathlib import Path

import polars as pl
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "benchmark_leave_one_out.py"
_spec = importlib.util.spec_from_file_location("benchmark_leave_one_out", SCRIPT)
bench = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bench)


def _genome(n: int = 2000, seed: int = 1) -> str:
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


class TestSimulator:
    def test_deterministic_given_seed(self):
        g = _genome()
        a = bench.simulate_reads(g, "GCF_1.1", 20, 100, 0.01, random.Random(7))
        b = bench.simulate_reads(g, "GCF_1.1", 20, 100, 0.01, random.Random(7))
        c = bench.simulate_reads(g, "GCF_1.1", 20, 100, 0.01, random.Random(8))
        assert a == b
        assert a != c

    def test_ids_lengths_and_no_global_state(self):
        g = _genome()
        random.seed(0)
        before = random.random()
        random.seed(0)
        reads = bench.simulate_reads(g, "GCF_1.1", 10, 80, 0.0, random.Random(1))
        assert random.random() == before
        assert [r[0] for r in reads] == [f"r{i}_GCF_1.1" for i in range(1, 11)]
        assert all(len(s) == 80 and s in g for _, s in reads)

    def test_error_rate_introduces_substitutions(self):
        g = _genome()
        reads = bench.simulate_reads(g, "A", 30, 100, 0.2, random.Random(3))
        assert any(s not in g for _, s in reads)

    def test_simulate_all_deterministic_and_unique_ids(self, tmp_path):
        files = {}
        for i, acc in enumerate(["GCF_2.1", "GCF_1.1"]):
            p = tmp_path / f"{acc}.fna"
            p.write_text(f">c{i}\n{_genome(1500, i)}\n")
            files[acc] = p
        r1 = bench.simulate_all(files, 5, 100, 0.005, 42)
        r2 = bench.simulate_all(files, 5, 100, 0.005, 42)
        assert r1 == r2
        ids = [r[0] for r in r1]
        assert len(set(ids)) == len(ids) == 10

    def test_accession_from_path(self):
        assert bench.accession_from_path(Path("x/GCF_1.1.fna.gz")) == "GCF_1.1"
        assert bench.accession_from_path(Path("GCF_1.1.fna")) == "GCF_1.1"


class TestTruth:
    @pytest.mark.parametrize(
        "ani,expected",
        [(99.0, "Known Species"), (96.0, "Known Species"), (95.9, "Novel Species"),
         (80.0, "Novel Species"), (79.9, "Novel Genus"), (70.0, "Novel Genus")],
    )
    def test_assign_truth(self, ani, expected):
        assert bench.assign_truth(ani) == expected


class TestLeaveOneOut:
    FULL = {
        "A": {"A": 100.0, "B": 90.0, "C": 75.0},
        "B": {"A": 90.0, "B": 100.0, "C": 76.0},
        "C": {"A": 75.0, "B": 76.0, "C": 100.0},
    }

    def test_ani_without_and_nearest(self):
        sub = bench.ani_without(self.FULL, "A")
        assert set(sub) == {"B", "C"}
        assert all("A" not in row for row in sub.values())
        assert bench.nearest_ani(self.FULL, "A") == 90.0
        assert bench.nearest_ani(self.FULL, "C") == 76.0

    def test_load_ani_dict(self, tmp_path):
        p = tmp_path / "ani.csv"
        pl.DataFrame({"genome": ["A", "B"], "A": [100.0, 90.0], "B": [90.0, 100.0]}).write_csv(p)
        assert bench.load_ani_dict(p)["B"]["A"] == 90.0

    def test_hit_filtering(self):
        aln = pl.DataFrame({
            "qseqid": ["r1_GCF_1.1", "r1_GCF_1.1", "r2_GCF_2.1"],
            "sseqid": ["GCF_1.1|c1", "GCF_2.1|c9", "GCF_1.1|c1"],
            "pident": [100.0, 85.0, 84.0],
        })
        prep = bench.prepare_alignment(aln)
        assert prep["src"].to_list() == ["GCF_1.1", "GCF_1.1", "GCF_2.1"]
        out = bench.leave_one_out_hits(prep, "GCF_1.1")
        assert out["genome_name"].to_list() == ["GCF_2.1"]
        assert "src" not in out.columns


def _results() -> pl.DataFrame:
    rows = []

    def add(cfg, setting, truth, ani, calls, src="g"):
        for c in calls:
            rows.append({"config": cfg, "setting": setting, "truth": truth,
                         "true_nearest_ani": ani, "taxonomic_call": c, "src": src})

    for cfg, ns_calls in (("base", ["Novel Species", "Known Species"]),
                          ("alt", ["Novel Species", "Novel Species"])):
        add(cfg, "loo", "Novel Species", 85.0, ns_calls, "a")           # NS group
        add(cfg, "loo", "Novel Species", 95.0, ["Novel Species", "Known Species"], "b")  # boundary
        add(cfg, "loo", "Novel Genus", 70.0, ["Novel Genus", "Novel Species", "Novel Genus", "Unclassified"], "c")
        add(cfg, "control", "Known Species", 100.0, ["Known Species"] * 3 + ["Novel Species"], "a")
    return pl.DataFrame(rows)


class TestMetrics:
    def test_values(self):
        m = bench.compute_metrics(_results())
        b, a = m["base"], m["alt"]
        assert b["novel_species_recall"] == 0.5 and a["novel_species_recall"] == 1.0
        assert b["boundary_called_novel_species"] == 0.5
        assert b["novel_genus_recall"] == 0.5
        assert b["novel_genus_called_novel_species"] == 0.25
        # base: called NS = 1(a)+1(b)+1(c)+1(control) = 4; truth NS among them = 2
        assert b["novel_species_precision"] == 0.5
        # a: NS calls 2(a)+1(b)+1(c)+1(ctl)=5; truth NS 3
        assert a["novel_species_precision"] == pytest.approx(0.6)
        # base: truth novel rows 2+2+4 = 8, called novel: 1+1+3 = 5
        assert b["any_novel_recall"] == pytest.approx(5 / 8)
        # base: called novel = 5 loo + 1 control = 6, truth novel 5
        assert b["any_novel_precision"] == pytest.approx(5 / 6)
        assert b["control_false_novel"] == 1 and b["n_reads_control"] == 4

    def test_empty_denominator_is_none(self):
        df = _results().filter(pl.col("true_nearest_ani") != 70.0)
        assert bench.compute_metrics(df)["base"]["novel_genus_recall"] is None

    def test_per_genome_and_summary(self):
        res = _results()
        t = bench.per_genome_table(res)
        assert set(t.columns) >= {"accession", "nearest_ani", "truth", "n_reads",
                                  "recall_base", "recall_alt"}
        row = t.filter(pl.col("accession") == "a").row(0, named=True)
        assert row["recall_base"] == 0.5 and row["recall_alt"] == 1.0
        assert "NS precision" in bench.format_summary(bench.compute_metrics(res))


class TestConfigs:
    def test_default_and_yaml_parsing(self, tmp_path):
        assert list(bench.build_configs([], False)) == ["default"]
        p = tmp_path / "c.yaml"
        from metadarkmatter.models.config import ScoringConfig
        ScoringConfig().to_yaml(p)
        cfgs = bench.build_configs([f"mine={p}"], False)
        assert list(cfgs) == ["mine"]
        with pytest.raises(SystemExit):
            bench.build_configs(["nopath"], False)

    def test_clade_rule(self):
        from metadarkmatter.models.config import ScoringConfig
        if "clade_novel_species" not in ScoringConfig.model_fields:
            pytest.skip("clade_novel_species not yet in ScoringConfig")
        cfgs = bench.build_configs([], True)
        assert list(cfgs) == ["default", "clade_rule"]
        assert cfgs["clade_rule"].clade_novel_species is True
