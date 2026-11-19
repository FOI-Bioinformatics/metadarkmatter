"""Hidden files (macOS AppleDouble ``._*`` sidecars) are not genomes."""

from __future__ import annotations

from pathlib import Path

import pytest

from metadarkmatter.cli.aai import _find_protein_files
from metadarkmatter.core.aai_matrix_builder import concatenate_protein_files
from metadarkmatter.core.genome_utils import (
    concatenate_genomes_with_mapping,
    find_genome_files,
    is_hidden_file,
)
from metadarkmatter.external.bowtie2 import concatenate_genomes
from metadarkmatter.external.fastani import create_genome_list_file

ACC = "GCA_001898995.1"


def _make(directory: Path, _suffix: str, names: list[str], text: str) -> None:
    for name in names:
        # Sidecars hold binary metadata, not FASTA.
        (directory / name).write_text(text if not name.startswith(".") else "\x00\x05\x16")


def test_is_hidden_file() -> None:
    assert is_hidden_file(Path("._A.fna"))
    assert is_hidden_file(Path("/x/.hidden.fna"))
    assert not is_hidden_file(Path("/x/A.fna"))


def test_find_genome_files_skips_hidden(tmp_path: Path) -> None:
    _make(tmp_path, ".fna", ["A.fna", "._A.fna", ".hidden.fna"], ">c\nACGT\n")
    assert [p.name for p in find_genome_files(tmp_path)] == ["A.fna"]
    assert [p.name for p in find_genome_files(tmp_path, "*.fna")] == ["A.fna"]


def test_concatenate_genomes_no_duplicate_accession_for_sidecar(tmp_path: Path) -> None:
    _make(tmp_path, ".fna", [f"{ACC}.fna", f"._{ACC}.fna"], ">c1\nACGT\n")
    count, contigs = concatenate_genomes_with_mapping(
        tmp_path, tmp_path / "out.fasta", tmp_path / "map.tsv"
    )
    assert (count, contigs) == (1, 1)


def test_bowtie2_concatenate_genomes_skips_sidecar(tmp_path: Path) -> None:
    _make(tmp_path, ".fna", [f"{ACC}.fna", f"._{ACC}.fna"], ">c1\nACGT\n")
    accessions = concatenate_genomes(tmp_path, tmp_path / "pan.fasta")
    assert accessions == [ACC]


def test_fastani_genome_list_skips_sidecar(tmp_path: Path) -> None:
    _make(tmp_path, ".fna", ["A.fna", "._A.fna"], ">c\nACGT\n")
    paths = create_genome_list_file(tmp_path, tmp_path / "list.txt")
    assert [p.name for p in paths] == ["A.fna"]


def test_concatenate_protein_files_skips_sidecar(tmp_path: Path) -> None:
    _make(tmp_path, ".faa", [f"{ACC}.faa", f"._{ACC}.faa"], ">p1\nMKV\n")
    mapping = concatenate_protein_files(tmp_path, tmp_path / "all.faa")
    assert len(mapping) == 1
    assert "._" not in next(iter(mapping.values()))


@pytest.mark.parametrize("pattern", ["*.faa", "*.faa.gz"])
def test_cli_aai_find_protein_files_skips_sidecar(tmp_path: Path, pattern: str) -> None:
    _make(tmp_path, ".faa", ["A.faa", "._A.faa"], ">p\nMKV\n")
    files, _ = _find_protein_files(tmp_path, pattern)
    assert [p.name for p in files] == ["A.faa"]
