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


# ------------------------------------------------------------ classifier core


class TestConfigPrecedence:
    def test_explicit_default_flag_overrides_yaml(self, tmp_path):
        cfg = tmp_path / "c.yaml"
        cfg.write_text("alignment_mode: protein\nfilters:\n  min_alignment_length: 50\n"
                       "uncertainty:\n  mode: max\n")
        aln = _write_alignment(tmp_path / "s.blast.tsv")
        ani = _write_ani(tmp_path / "ani.csv")
        out = tmp_path / "c.csv"
        base = ["score", "classify", "--alignment", str(aln), "--ani", str(ani),
                "--output", str(out), "--config", str(cfg)]
        res = runner.invoke(app, [*base, "--min-alignment-length", "100",
                                  "--uncertainty-mode", "second", "--alignment-mode", "nucleotide"])
        assert res.exit_code == 0, res.output
        c = json.loads((tmp_path / "c.csv.provenance.json").read_text())["config"]
        assert (c["alignment_mode"], c["min_alignment_length"], c["uncertainty_mode"]) == (
            "nucleotide", 100, "second")
        # Without flags the YAML values stay.
        assert runner.invoke(app, base).exit_code == 0
        c = json.loads((tmp_path / "c.csv.provenance.json").read_text())["config"]
        assert (c["alignment_mode"], c["min_alignment_length"], c["uncertainty_mode"]) == (
            "protein", 50, "max")


class TestConfigValidation:
    @pytest.mark.parametrize("kw", [{"novelty_novel_genus_max": 15.0},
                                    {"uncertainty_conserved_min": 0.0}])
    def test_reversed_ranges_rejected(self, kw):
        from metadarkmatter.models.config import ScoringConfig
        with pytest.raises(ValueError, match="must be <="):
            ScoringConfig(**kw)

    def test_taxonomy_yaml_round_trip(self, tmp_path):
        from metadarkmatter.models.config import ScoringConfig
        c = ScoringConfig(same_genus_ani_threshold=82.0, use_aai_for_genus=False,
                          aai_genus_boundary_low=60.0)
        p = tmp_path / "c.yaml"
        c.to_yaml(p)
        b = ScoringConfig.from_yaml(p)
        assert (b.same_genus_ani_threshold, b.use_aai_for_genus, b.aai_genus_boundary_low) == (
            82.0, False, 60.0)


class TestClassifyInputs:
    def _base(self, tmp_path, aln):
        ani = _write_ani(tmp_path / "ani.csv")
        return ["score", "classify", "--alignment", str(aln), "--ani", str(ani),
                "--output", str(tmp_path / "c.csv"), *FILTERS]

    def test_streaming_with_qc_output_is_an_error(self, tmp_path):
        aln = _write_alignment(tmp_path / "s.blast.tsv")
        res = runner.invoke(app, [*self._base(tmp_path, aln), "--streaming",
                                  "--qc-output", str(tmp_path / "qc.json")])
        assert res.exit_code == 1
        assert "--qc-output" in res.output

    def test_gzip_without_gz_suffix(self, tmp_path):
        gz = _write_alignment(tmp_path / "s.blast.tsv.gz")
        plain_name = tmp_path / "s.blast.tsv"
        gz.rename(plain_name)
        res = runner.invoke(app, self._base(tmp_path, plain_name))
        assert res.exit_code == 0, res.output
        assert pl.read_csv(tmp_path / "c.csv").height > 0

    def test_comment_line_keeps_coverage_check(self, tmp_path, caplog):
        aln = _write_alignment(tmp_path / "s.blast.tsv")
        aln.write_text("# BLASTN 2.15\n" + aln.read_text())
        with caplog.at_level("WARNING"):
            res = runner.invoke(app, self._base(tmp_path, aln))
        assert res.exit_code == 0, res.output
        assert "Could not validate" not in caplog.text


def test_clade_rule_stays_inside_preset_species_range():
    from metadarkmatter.core.classification.ani_matrix import ANIMatrix
    from metadarkmatter.core.classification.classifiers.vectorized import VectorizedClassifier
    from metadarkmatter.models.config import ScoringConfig

    ani = ANIMatrix({"A": {"A": 100.0, "B": 92.0}, "B": {"A": 92.0, "B": 100.0}})
    rows = [("r16", "A|c", 84.0, 200), ("r16", "B|c", 83.7, 199),   # N = 16
            ("r9", "A|c", 91.0, 220), ("r9", "B|c", 90.6, 218)]     # N = 9
    df = pl.DataFrame({
        "qseqid": [r[0] for r in rows], "sseqid": [r[1] for r in rows],
        "pident": [r[2] for r in rows], "length": [150] * 4, "mismatch": [0] * 4,
        "gapopen": [0] * 4, "qstart": [1] * 4, "qend": [150] * 4, "sstart": [1] * 4,
        "send": [150] * 4, "evalue": [1e-40] * 4, "bitscore": [float(r[3]) for r in rows],
        "qlen": [150] * 4,
    }).with_columns(pl.col("sseqid").str.extract(r"^([^|]+)\|").alias("genome_name"))
    cfg = ScoringConfig(novelty_novel_species_max=15.0, novelty_novel_genus_min=15.0)
    out = VectorizedClassifier(ani, config=cfg).classify_file(df)
    flags = dict(zip(out["read_id"], out["clade_novel_species"], strict=True))
    assert flags == {"r9": True, "r16": False}


# ------------------------------------------------------- alignment / pipeline


def _fq(path: Path, ids: list[str]) -> Path:
    path.write_text("".join(f"@{i}\nACGTACGT\n+\nIIIIIIII\n" for i in ids))
    return path


@pytest.fixture
def fake_tools():
    from metadarkmatter.external.base import ExternalTool

    ExternalTool.set_executable_resolver(lambda name: f"/usr/bin/{name}")
    ExternalTool.clear_cache()
    yield
    ExternalTool.reset_executable_resolver()


def test_mmseqs2_dry_run_touches_nothing(tmp_path, fake_tools):
    r1, r2 = _fq(tmp_path / "a_R1.fq", ["x"]), _fq(tmp_path / "a_R2.fq", ["x"])
    db = tmp_path / "db"
    db.write_text("")
    out = tmp_path / "new" / "deep" / "o.tsv.gz"
    res = runner.invoke(app, ["mmseqs2", "search", "--query-1", str(r1), "--query-2", str(r2),
                              "--database", str(db), "--output", str(out), "--dry-run"])
    assert res.exit_code == 0, res.output
    assert not (tmp_path / "new").exists()


@pytest.mark.parametrize("tool", ["blast", "mmseqs2"])
def test_compressed_output_temp_is_derived_from_output(tmp_path, fake_tools, tool):
    """--output res.blast.gz used to use (then delete) an unrelated res.tsv."""
    q = _fq(tmp_path / "q.fq", ["r1"])
    db = tmp_path / "db" / "pan"
    db.parent.mkdir()
    for ext in (".nhr", ".nin", ".nsq", ".dbtype", ""):
        (tmp_path / "db" / f"pan{ext}").write_text("")
    args = (["blast", "align", "--query", str(q)] if tool == "blast"
            else ["mmseqs2", "search", "--query", str(q)])
    res = runner.invoke(app, [*args, "--database", str(db),
                              "--output", str(tmp_path / "res.blast.gz"), "--dry-run"],
                        env={"COLUMNS": "300"})
    assert res.exit_code == 0, res.output
    assert "res.blast.gz.partial" in res.output
    assert "res.tsv" not in res.output


def test_sample_sheet_with_utf8_bom(tmp_path):
    from metadarkmatter.core.pipeline.cohort import read_samples_tsv

    reads = _fq(tmp_path / "b.fq", ["x"])
    sheet = tmp_path / "s.tsv"
    sheet.write_bytes(b"\xef\xbb\xbfsample\treads_1\ns1\t" + reads.name.encode() + b"\n")
    assert [s.name for s in read_samples_tsv(sheet)] == ["s1"]


def test_run_pipeline_dry_run_with_markup_in_path(tmp_path):
    g = tmp_path / "g"
    g.mkdir()
    (g / "GCF_000000001.1.fna").write_text(">c\nACGT\n")
    r1 = _fq(tmp_path / "x_R1.fq", ["x"])
    res = runner.invoke(app, ["run", "pipeline", "--genomes", str(g), "--extracted-reads-1",
                              str(r1), "--output-dir", str(tmp_path / "out[/bold]"),
                              "--no-report", "--dry-run"])
    assert res.exit_code == 0, res.output
    assert "[dry-run]" in res.output


# ---------------------------------------------------------- report / grouping


def _cls_df(rows: list[tuple[str, str, float, str]]) -> pl.DataFrame:
    return pl.DataFrame({
        "read_id": [r[0] for r in rows], "best_match_genome": [r[1] for r in rows],
        "top_hit_identity": [100.0 - r[2] for r in rows], "novelty_index": [r[2] for r in rows],
        "placement_uncertainty": [0.5] * len(rows), "num_ambiguous_hits": [1] * len(rows),
        "taxonomic_call": [r[3] for r in rows], "confidence_score": [80.0] * len(rows),
    })


def _md(tmp_path: Path, genomes: dict[str, tuple[str, str]]) -> Path:
    p = tmp_path / "md.tsv"
    p.write_text("accession\tspecies\tgenus\tfamily\n" + "".join(
        f"{a}\t{sp}\t{ge}\tFam\n" for a, (ge, sp) in genomes.items()))
    return p


def test_genus_grouping_reference_missing_from_ani_does_not_chain(tmp_path):
    from metadarkmatter.core.classification.ani_matrix import ANIMatrix
    from metadarkmatter.core.metadata import GenomeMetadata
    from metadarkmatter.core.novel_diversity.clustering import NovelDiversityAnalyzer

    md = GenomeMetadata.from_file(_md(tmp_path, {"A": ("G", "G a"), "B": ("G", "G b"),
                                                 "C": ("G", "G c")}))
    ani = ANIMatrix({"A": {"A": 100.0, "B": 70.0}, "B": {"A": 70.0, "B": 100.0}})
    df = _cls_df([(f"a{i}", "A", 10.0, "Novel Species") for i in range(3)]
                 + [(f"b{i}", "B", 10.0, "Novel Species") for i in range(3)]
                 + [("c0", "C", 10.0, "Novel Species")])
    clusters = NovelDiversityAnalyzer(df, metadata=md, ani_matrix=ani,
                                      min_cluster_size=1).cluster_novel_reads()
    assert sorted(c.read_count for c in clusters) == [1, 3, 3]


def test_summary_counts_reads_in_groups(tmp_path):
    from metadarkmatter.core.metadata import GenomeMetadata
    from metadarkmatter.core.novel_diversity.clustering import NovelDiversityAnalyzer

    md = GenomeMetadata.from_file(_md(tmp_path, {"A": ("G", "G a"), "Z": ("H", "H z")}))
    df = _cls_df([(f"a{i}", "A", 10.0, "Novel Species") for i in range(4)]
                 + [("a9", "A", 22.0, "Novel Genus"), ("z0", "Z", 22.0, "Novel Genus")])
    s = NovelDiversityAnalyzer(df, metadata=md, min_cluster_size=2).get_summary()
    assert (s.novel_species_clusters, s.novel_species_cluster_reads) == (1, 5)
    assert (s.novel_genus_clusters, s.novel_genus_cluster_reads) == (0, 0)
    assert (s.grouped_reads, s.total_novel_reads) == (5, 6)


def test_summary_without_groups_keeps_read_totals():
    from metadarkmatter.core.novel_diversity.clustering import NovelDiversityAnalyzer

    df = _cls_df([("a", "A", 10.0, "Novel Species")])
    s = NovelDiversityAnalyzer(df, min_cluster_size=3).get_summary()
    assert s.total_clusters == 0 and s.total_novel_reads == 1


def test_band_size_must_be_whole():
    from metadarkmatter.core.novel_diversity.clustering import NovelDiversityAnalyzer

    with pytest.raises(ValueError, match="whole number"):
        NovelDiversityAnalyzer(_cls_df([("a", "A", 10.0, "Novel Species")]),
                               novelty_band_size=2.5)


def test_species_breakdown_ignores_nan_identity():
    from metadarkmatter.visualization.report.breakdowns import compute_species_breakdown

    df = _cls_df([("a", "A", 10.0, "Novel Species"), ("b", "A", 12.0, "Novel Species")])
    df = df.with_columns(pl.Series("top_hit_identity", [float("nan"), 88.0]))
    table, _cats, _uses_species = compute_species_breakdown(df, None)
    assert table["mean_identity"].to_list() == [88.0]


def test_data_table_taxa_fall_back_to_table_columns(tmp_path):
    from metadarkmatter.core.metadata import GenomeMetadata
    from metadarkmatter.visualization.report.data_table import _genome_taxa

    md = GenomeMetadata.from_file(_md(tmp_path, {"A": ("G", "G a")}))
    df = pl.DataFrame({"best_match_genome": ["A", "B"], "species": ["G a", "K b"],
                       "genus": ["G", "K"]})
    assert _genome_taxa(["A", "B"], df, md) == (["G a", "K b"], ["G", "K"])


def test_report_multi_rejects_sample_name_clash(tmp_path):
    d = tmp_path / "multi"
    d.mkdir()
    for name in ("s1.csv", "s1_classifications.csv"):
        _cls_df([("a", "A", 10.0, "Novel Species")]).write_csv(d / name)
    res = runner.invoke(app, ["report", "multi", "-i", str(d), "-p", "*",
                              "-o", str(tmp_path / "m.html")])
    assert res.exit_code == 1
    assert "same sample name" in res.output
