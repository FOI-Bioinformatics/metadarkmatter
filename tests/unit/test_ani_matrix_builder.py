"""Tests for core.ani_matrix_builder (fastANI/skani parsing and CSV conversion)."""

from __future__ import annotations

import gzip
from pathlib import Path

import polars as pl
import pytest

from metadarkmatter.core.ani_matrix_builder import (
    ani_dict_to_csv,
    extract_genome_name_from_path,
    parse_fastani_output,
    parse_skani_output,
    validate_ani_coverage,
)


class TestExtractGenomeNameFromPath:
    def test_ncbi_genomic_filename(self):
        assert (
            extract_genome_name_from_path("/path/to/GCF_000123456.1_ASM123v1_genomic.fna")
            == "GCF_000123456.1"
        )

    def test_prefixed_filename(self):
        assert (
            extract_genome_name_from_path("Francisellaceae_Genus_species_GCA_000710735.1.fasta")
            == "GCA_000710735.1"
        )

    def test_pipe_format_sseqid(self):
        assert extract_genome_name_from_path("GCF_000123456.1|NZ_CP007557.1") == "GCF_000123456.1"

    def test_custom_pipe_name_without_accession(self):
        assert extract_genome_name_from_path("custom_genome|contig1") == "custom_genome"

    def test_simple_filename_fallback(self):
        assert extract_genome_name_from_path("my_genome.fasta") == "my_genome"

    def test_compound_gz_extension(self):
        assert extract_genome_name_from_path("GCF_000123456.1.fna.gz") == "GCF_000123456.1"


class TestParseFastaniOutput:
    def test_symmetric_with_diagonal_and_missing(self, tmp_path: Path):
        out = tmp_path / "fastani.tsv"
        # Only one direction provided; parser must symmetrize and fill diagonal.
        out.write_text(
            "GCF_000000001.1_genomic.fna\tGCF_000000002.1_genomic.fna\t98.5\t100\t120\n"
        )
        result = parse_fastani_output(out)

        assert result["GCF_000000001.1"]["GCF_000000002.1"] == 98.5
        # symmetric
        assert result["GCF_000000002.1"]["GCF_000000001.1"] == 98.5
        # diagonal
        assert result["GCF_000000001.1"]["GCF_000000001.1"] == 100.0
        assert result["GCF_000000002.1"]["GCF_000000002.1"] == 100.0

    def test_skips_comments_and_short_lines(self, tmp_path: Path):
        out = tmp_path / "fastani.tsv"
        out.write_text(
            "# header comment\n"
            "\n"
            "GCF_000000001.1.fna\tGCF_000000002.1.fna\t97.0\n"
            "truncated_line_no_ani\n"
        )
        result = parse_fastani_output(out)
        assert result["GCF_000000001.1"]["GCF_000000002.1"] == 97.0


class TestParseSkaniOutput:
    def test_full_matrix(self, tmp_path: Path):
        out = tmp_path / "skani_full.txt"
        out.write_text(
            "3\n"
            "GCF_000000001.1.fna 100.0 95.0 80.0\n"
            "GCF_000000002.1.fna 95.0 100.0 82.0\n"
            "GCF_000000003.1.fna 80.0 82.0 100.0\n"
        )
        result = parse_skani_output(out, full_matrix=True)
        assert result["GCF_000000001.1"]["GCF_000000002.1"] == 95.0
        assert result["GCF_000000003.1"]["GCF_000000002.1"] == 82.0
        assert result["GCF_000000003.1"]["GCF_000000003.1"] == 100.0

    def test_full_matrix_empty_file(self, tmp_path: Path):
        out = tmp_path / "empty.txt"
        out.write_text("")
        assert parse_skani_output(out, full_matrix=True) == {}

    def test_pairwise(self, tmp_path: Path):
        out = tmp_path / "skani_pair.tsv"
        out.write_text(
            "Ref_file\tQuery_file\tANI\tAF_ref\tAF_query\n"
            "GCF_000000001.1.fna\tGCF_000000002.1.fna\t96.4\t90\t91\n"
        )
        result = parse_skani_output(out, full_matrix=False)
        assert result["GCF_000000001.1"]["GCF_000000002.1"] == 96.4
        assert result["GCF_000000002.1"]["GCF_000000001.1"] == 96.4
        assert result["GCF_000000001.1"]["GCF_000000001.1"] == 100.0


class TestAniDictToCsv:
    def test_writes_square_matrix(self, tmp_path: Path):
        ani = {
            "B": {"B": 100.0, "A": 95.0},
            "A": {"A": 100.0, "B": 95.0},
        }
        out = tmp_path / "ani.csv"
        n = ani_dict_to_csv(ani, out)
        assert n == 2

        df = pl.read_csv(out)
        # genomes sorted -> A, B
        assert df["genome"].to_list() == ["A", "B"]
        assert "A" in df.columns and "B" in df.columns

    def test_compressed_output(self, tmp_path: Path):
        ani = {"A": {"A": 100.0}}
        out = tmp_path / "ani.csv.gz"
        ani_dict_to_csv(ani, out)
        with gzip.open(out, "rt") as fh:
            content = fh.read()
        assert "genome" in content
        assert "A" in content

    def test_empty_dict_raises(self, tmp_path: Path):
        with pytest.raises(ValueError, match="empty"):
            ani_dict_to_csv({}, tmp_path / "x.csv")


class TestValidateAniCoverage:
    def test_full_coverage(self):
        matched, total, pct, missing = validate_ani_coverage(
            {"A", "B", "C"}, {"A", "B"}
        )
        assert matched == 2
        assert total == 2
        assert pct == 100.0
        assert missing == set()

    def test_partial_coverage(self):
        matched, total, pct, missing = validate_ani_coverage({"A"}, {"A", "B", "C"})
        assert matched == 1
        assert total == 3
        assert pct == pytest.approx(33.333, abs=0.01)
        assert missing == {"B", "C"}

    def test_empty_blast_genomes(self):
        matched, total, pct, missing = validate_ani_coverage({"A"}, set())
        assert matched == 0
        assert total == 0
        assert pct == 0.0
        assert missing == set()
