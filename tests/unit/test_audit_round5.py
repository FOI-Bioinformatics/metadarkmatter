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
