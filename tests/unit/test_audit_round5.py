"""Regression tests for the round-5 audit fixes."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.main import app
from metadarkmatter.core.exceptions import ConfigurationError
from metadarkmatter.core.genome_utils import concatenate_genomes_with_mapping, find_genome_files

runner = CliRunner()


def _genome(path: Path, text: str) -> Path:
    if path.name.endswith(".gz"):
        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text)
    return path


def test_mixed_plain_and_gzipped_genomes_all_used(tmp_path):
    g = tmp_path / "g"
    g.mkdir()
    _genome(g / "GCF_000000001.1.fna", ">a\nACGT\n")
    _genome(g / "GCF_000000002.1.fna.gz", ">b\nACGT\n")
    _genome(g / "GCF_000000003.1.fasta", ">c\nACGT\n")
    (g / "notes.txt").write_text("x")
    assert [p.name for p in find_genome_files(g)] == [
        "GCF_000000001.1.fna", "GCF_000000002.1.fna.gz", "GCF_000000003.1.fasta"]
    n, contigs = concatenate_genomes_with_mapping(g, tmp_path / "o.fa", tmp_path / "m.tsv")
    assert (n, contigs) == (3, 3)


def test_explicit_pattern_is_respected_with_fallback(tmp_path):
    g = tmp_path / "g"
    g.mkdir()
    _genome(g / "GCF_000000001.1.fna", ">a\nACGT\n")
    _genome(g / "GCF_000000002.1.fa", ">b\nACGT\n")
    assert [p.name for p in find_genome_files(g, "*.fa")] == ["GCF_000000002.1.fa"]
    assert len(find_genome_files(g, "*.xyz")) == 2


@pytest.mark.parametrize("files,match", [
    ({"GCF_000000001.1.fna": ">a\nAC\n", "GCF_000000001.1_ASM1_genomic.fna": ">b\nAC\n"},
     "map to accession"),
    ({"GCF_000000009.1.fna": "<html>error</html>\nACGT\n"}, "before the first FASTA header"),
    ({"GCF_000000005.1.fna": ""}, "No FASTA records"),
])
def test_concatenation_errors(tmp_path, files, match):
    g = tmp_path / "g"
    g.mkdir()
    for name, text in files.items():
        _genome(g / name, text)
    with pytest.raises(ConfigurationError, match=match):
        concatenate_genomes_with_mapping(g, tmp_path / "o.fa", tmp_path / "m.tsv")
    assert not (tmp_path / "o.fa").exists() or match != "No FASTA records"


def test_empty_genome_file_not_counted(tmp_path):
    g = tmp_path / "g"
    g.mkdir()
    _genome(g / "GCF_000000001.1.fna", ">a\nACGT\n")
    _genome(g / "GCF_000000002.1.fna", "")
    n, contigs = concatenate_genomes_with_mapping(g, tmp_path / "o.fa", tmp_path / "m.tsv")
    assert (n, contigs) == (1, 1)


def test_corrupt_gzip_genome_is_configuration_error(tmp_path):
    g = tmp_path / "g"
    g.mkdir()
    import random

    rng = random.Random(1)
    seq = "".join(rng.choice("ACGT") for _ in range(20000))
    data = bytearray(gzip.compress((">a\n" + seq + "\n").encode()))
    assert len(data) > 2000
    mid = len(data) // 2
    data[mid:mid + 64] = bytes(64)
    (g / "GCF_000000001.1.fna.gz").write_bytes(bytes(data))
    with pytest.raises(ConfigurationError, match="Could not read genome file"):
        concatenate_genomes_with_mapping(g, tmp_path / "o.fa", tmp_path / "m.tsv")


def test_aai_default_pattern_keeps_gz_proteins(tmp_path):
    from metadarkmatter.cli.aai import _find_protein_files

    (tmp_path / "GCF_000000001.1.faa").write_text(">p\nMKV\n")
    (tmp_path / "GCF_000000002.1.faa.gz").write_bytes(gzip.compress(b">p\nMKV\n"))
    files, patterns = _find_protein_files(tmp_path, "*.faa")
    assert len(files) == 2 and "*.faa.gz" in patterns


def test_proteins_accession_from_filename():
    from metadarkmatter.core.genome_utils import extract_accession_from_filename

    assert extract_accession_from_filename("GCF_000833475.1.fasta.gz") == "GCF_000833475.1"
    assert extract_accession_from_filename(
        "GCF_000815225.1_ASM81522v1_genomic.fna") == "GCF_000815225.1"


@pytest.mark.parametrize("header,expected", [
    ('genome\t"A, sp. 1"\t"B, sp. 2"\n', "\t"),
    ("genome\tA,1\tB,2\n", "\t"),
    ("genome,A,B\n", ","),
])
def test_sniff_separator_with_commas_in_names(tmp_path, header, expected):
    from metadarkmatter.core.parsers import sniff_matrix_separator

    p = tmp_path / "m.txt"
    p.write_text(header)
    assert sniff_matrix_separator(p) == expected


def test_blast_makedb_rejects_space_in_output(tmp_path):
    g = tmp_path / "g"
    g.mkdir()
    _genome(g / "GCF_000000001.1.fna", ">a\nACGT\n")
    res = runner.invoke(app, ["blast", "makedb", "-g", str(g), "-o", str(tmp_path / "my db")])
    assert res.exit_code == 1
    assert "must not contain spaces" in res.output
    assert not list(tmp_path.glob("my db*"))


def test_mmseqs2_makedb_dry_run(tmp_path):
    from metadarkmatter.external.base import ExternalTool

    g = tmp_path / "g"
    g.mkdir()
    _genome(g / "GCF_000000001.1.fna", ">a\nACGT\n")
    ExternalTool.set_executable_resolver(lambda name: f"/usr/bin/{name}")
    ExternalTool.clear_cache()
    try:
        res = runner.invoke(app, ["mmseqs2", "makedb", "-g", str(g),
                                  "-o", str(tmp_path / "out" / "db"), "--dry-run"])
    finally:
        ExternalTool.reset_executable_resolver()
    assert res.exit_code == 0, res.output
    assert not (tmp_path / "out").exists()


def test_capture_output_false_keeps_stderr(tmp_path):
    import sys

    from metadarkmatter.external.base import ExternalTool

    class Py(ExternalTool):
        TOOL_NAME = "python"

        def build_command(self, **kwargs):
            return [sys.executable, "-c", "import sys; sys.stderr.write('boom ' * 2000); sys.exit(3)"]

    ExternalTool.set_executable_resolver(lambda name: sys.executable)
    ExternalTool.clear_cache()
    try:
        result = Py().run(capture_output=False)
    finally:
        ExternalTool.reset_executable_resolver()
    assert result.return_code == 3
    assert "boom" in result.stderr and len(result.stderr) <= 4000


def test_cli_command_normalises_module_argv(monkeypatch):
    from metadarkmatter.cli.score import _cli_command

    monkeypatch.setattr("sys.argv", ["/x/metadarkmatter/cli/main.py", "score", "classify"])
    assert _cli_command() == ["mdm", "score", "classify"]


# ------------------------------------------------------------ plots / report


def test_report_thresholds_follow_alignment_mode():
    from metadarkmatter.models.config import ScoringConfig
    from metadarkmatter.visualization.plots.base import ThresholdConfig

    t = ThresholdConfig.from_effective(
        ScoringConfig(alignment_mode="protein").get_effective_thresholds())
    assert t.novelty_known_max == 10.0 and t.uncertainty_known_max == 5.0
    n = ThresholdConfig.from_effective(ScoringConfig().get_effective_thresholds())
    assert (n.novelty_known_max, n.identity_known_min) == (4.0, 96.0)


def test_scatter_legend_counts_use_all_reads():
    import polars as pl

    from metadarkmatter.visualization.plots.scatter_2d import NoveltyUncertaintyScatter

    df = pl.DataFrame({
        "read_id": [f"r{i}" for i in range(200)],
        "novelty_index": [1.0] * 190 + [22.0] * 10,
        "placement_uncertainty": [0.2] * 200,
        "num_ambiguous_hits": [2] * 200,
        "taxonomic_call": ["Known Species"] * 190 + ["Novel Genus"] * 10,
    })
    fig = NoveltyUncertaintyScatter(df, max_points=50).create_figure()
    names = [t.name for t in fig.data if t.name]
    assert any("(190)" in n for n in names) and any("(10)" in n for n in names)
    assert "points shown: 50 of 200" in (fig.layout.legend.title.text or "")


def test_extended_matrix_aligns_columns_and_fills_triangle():
    import polars as pl

    from metadarkmatter.visualization.report.components.extended_matrix_builder import (
        build_extended_similarity_matrix,
    )

    m = pl.DataFrame({"genome": ["A", "B", "C"], "C": [90.0, 85.0, 100.0],
                      "B": [95.0, 100.0, 0.0], "A": [100.0, 0.0, 0.0]})
    mat, labels, _ = build_extended_similarity_matrix(m, [], {}, max_references=3)
    idx = {lab: i for i, lab in enumerate(labels)}
    assert mat[idx["A"], idx["B"]] == mat[idx["B"], idx["A"]] == 95.0
    assert mat[idx["B"], idx["C"]] == mat[idx["C"], idx["B"]] == 85.0


def test_genome_positions_offset_contigs():
    import polars as pl

    from metadarkmatter.core.recruitment import add_genome_positions

    df = pl.DataFrame({"genome_name": ["G", "G", "H"], "contig": ["G|c1", "G|c2", "H|c1"],
                       "position": [10, 10, 5]})
    out = add_genome_positions(df, {"G|c1": 1000, "G|c2": 500, "H|c1": 300})
    assert out["genome_position"].to_list() == [10, 1010, 5]


def test_multi_report_without_placement_uncertainty(tmp_path):
    import polars as pl

    d = tmp_path / "multi"
    d.mkdir()
    for s in ("s1", "s2"):
        pl.DataFrame({"read_id": ["a", "b"], "best_match_genome": ["G", "G"],
                      "novelty_index": [1.0, 10.0],
                      "taxonomic_call": ["Known Species", "Novel Species"]}).write_csv(
            d / f"{s}_classifications.csv")
    res = runner.invoke(app, ["report", "multi", "-i", str(d), "-o", str(tmp_path / "m.html")])
    assert res.exit_code == 0, res.output


def test_report_with_inf_and_nan_values(tmp_path):
    import polars as pl

    cls = tmp_path / "c.csv"
    pl.DataFrame({"read_id": ["a", "b", "c"], "best_match_genome": ["G", "G", "G"],
                  "top_hit_identity": [99.0, 90.0, 85.0],
                  "novelty_index": [1.0, float("inf"), float("nan")],
                  "placement_uncertainty": [0.1, 0.2, 0.3],
                  "taxonomic_call": ["Known Species", "Novel Species", "Novel Species"]}
                 ).write_csv(cls)
    res = runner.invoke(app, ["report", "generate", "-c", str(cls),
                              "-o", str(tmp_path / "r.html")])
    assert res.exit_code == 0, res.output
