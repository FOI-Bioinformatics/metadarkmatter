"""
Unit tests for ContigIdMapping class.

Tests ID mapping generation, transformation, and serialization for
external tool results processing.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import polars as pl
import pytest

from metadarkmatter.core.id_mapping import ContigIdMapping


class TestContigIdMappingFromGenomeDir:
    """Tests for ContigIdMapping.from_genome_dir()."""

    def test_generates_mapping_from_single_fasta(self, temp_dir: Path):
        """Should extract accession and map contigs from a single FASTA."""
        # Create a test genome file
        genome_file = temp_dir / "GCF_000195955.2_genomic.fna"
        genome_file.write_text(
            ">NZ_CP007557.1 Francisella tularensis subsp. tularensis\n"
            "ATGCATGCATGC\n"
            ">NZ_CP007558.1 Francisella tularensis plasmid\n"
            "GCATGCATGCAT\n"
        )

        mapping = ContigIdMapping.from_genome_dir(temp_dir)

        assert len(mapping) == 2
        assert "NZ_CP007557.1" in mapping
        assert "NZ_CP007558.1" in mapping
        assert mapping.contig_to_accession["NZ_CP007557.1"] == "GCF_000195955.2"
        assert mapping.contig_to_accession["NZ_CP007558.1"] == "GCF_000195955.2"

    def test_generates_mapping_from_multiple_fastas(self, temp_dir: Path):
        """Should handle multiple genome files."""
        # Create two genome files
        (temp_dir / "GCF_000123456.1_genomic.fna").write_text(
            ">NC_001234.1 Genome 1 contig\nATGC\n"
        )
        (temp_dir / "GCA_000789012.2_genomic.fna").write_text(
            ">scaffold_001 Genome 2 contig\nGCAT\n"
        )

        mapping = ContigIdMapping.from_genome_dir(temp_dir)

        assert len(mapping) == 2
        assert mapping.contig_to_accession["NC_001234.1"] == "GCF_000123456.1"
        assert mapping.contig_to_accession["scaffold_001"] == "GCA_000789012.2"

    def test_extracts_contig_id_from_complex_header(self, temp_dir: Path):
        """Should extract only the first word from FASTA headers."""
        genome_file = temp_dir / "GCF_000111222.1_genomic.fna"
        genome_file.write_text(
            ">NZ_ABCD01000001.1 description text species info\nATGC\n"
        )

        mapping = ContigIdMapping.from_genome_dir(temp_dir)

        assert len(mapping) == 1
        assert "NZ_ABCD01000001.1" in mapping
        # Should NOT include the description
        assert "NZ_ABCD01000001.1 description" not in mapping

    def test_handles_fa_extension(self, temp_dir: Path):
        """Should handle .fa files with alternative patterns."""
        genome_file = temp_dir / "GCF_000333444.1_genomic.fa"
        genome_file.write_text(">contig1\nATGC\n")

        mapping = ContigIdMapping.from_genome_dir(temp_dir, pattern="*.fa")

        assert len(mapping) == 1
        assert mapping.contig_to_accession["contig1"] == "GCF_000333444.1"

    def test_raises_for_nonexistent_directory(self):
        """Should raise FileNotFoundError for missing directory."""
        with pytest.raises(FileNotFoundError, match="not found"):
            ContigIdMapping.from_genome_dir(Path("/nonexistent/path"))

    def test_raises_for_empty_directory(self, temp_dir: Path):
        """Should raise ValueError when no genome files found."""
        with pytest.raises(ValueError, match="No genome files found"):
            ContigIdMapping.from_genome_dir(temp_dir)


class TestContigIdMappingTransform:
    """Tests for ID transformation methods."""

    @pytest.fixture
    def sample_mapping(self) -> ContigIdMapping:
        """Create a sample mapping for testing."""
        return ContigIdMapping(
            contig_to_accession={
                "NZ_CP007557.1": "GCF_000195955.2",
                "NC_006570.2": "GCF_000008985.1",
                "scaffold_001": "GCA_000789012.1",
            }
        )

    def test_transform_sseqid_with_known_contig(self, sample_mapping: ContigIdMapping):
        """Should transform known contig ID to standardized format."""
        result = sample_mapping.transform_sseqid("NZ_CP007557.1")
        assert result == "GCF_000195955.2|NZ_CP007557.1"

    def test_transform_sseqid_preserves_unknown(self, sample_mapping: ContigIdMapping):
        """Should pass through unknown contig IDs unchanged."""
        result = sample_mapping.transform_sseqid("unknown_contig_123")
        assert result == "unknown_contig_123"

    def test_transform_sseqid_already_standardized(self, sample_mapping: ContigIdMapping):
        """Should not double-transform already standardized IDs."""
        result = sample_mapping.transform_sseqid("GCF_000195955.2|NZ_CP007557.1")
        assert result == "GCF_000195955.2|NZ_CP007557.1"

    def test_transform_column_vectorized(self, sample_mapping: ContigIdMapping):
        """Should transform entire DataFrame column efficiently."""
        df = pl.DataFrame({
            "sseqid": [
                "NZ_CP007557.1",
                "NC_006570.2",
                "unknown_contig",
            ],
            "pident": [98.5, 95.0, 90.0],
        })

        result = sample_mapping.transform_column(df, "sseqid")

        expected = [
            "GCF_000195955.2|NZ_CP007557.1",
            "GCF_000008985.1|NC_006570.2",
            "unknown_contig",
        ]
        assert result["sseqid"].to_list() == expected
        # Verify other columns preserved
        assert result["pident"].to_list() == [98.5, 95.0, 90.0]

    def test_transform_column_preserves_already_standardized(
        self, sample_mapping: ContigIdMapping
    ):
        """Should not double-transform already standardized IDs in column."""
        df = pl.DataFrame({
            "sseqid": [
                "GCF_000195955.2|NZ_CP007557.1",
                "NZ_CP007557.1",
            ]
        })

        result = sample_mapping.transform_column(df, "sseqid")

        # First one unchanged (already has pipe), second one transformed
        assert result["sseqid"][0] == "GCF_000195955.2|NZ_CP007557.1"
        assert result["sseqid"][1] == "GCF_000195955.2|NZ_CP007557.1"

    def test_transform_column_raises_for_missing_column(
        self, sample_mapping: ContigIdMapping
    ):
        """Should raise ValueError when column doesn't exist."""
        df = pl.DataFrame({"other_col": ["a", "b"]})

        with pytest.raises(ValueError, match="not found"):
            sample_mapping.transform_column(df, "sseqid")


class TestContigIdMappingTsv:
    """Tests for TSV serialization."""

    @pytest.fixture
    def sample_mapping(self) -> ContigIdMapping:
        """Create a sample mapping for testing."""
        return ContigIdMapping(
            contig_to_accession={
                "NZ_CP007557.1": "GCF_000195955.2",
                "NC_006570.2": "GCF_000008985.1",
            }
        )

    def test_to_tsv_and_from_tsv_roundtrip(
        self, sample_mapping: ContigIdMapping, temp_dir: Path
    ):
        """Should roundtrip through TSV file."""
        tsv_path = temp_dir / "mapping.tsv"

        sample_mapping.to_tsv(tsv_path)
        loaded = ContigIdMapping.from_tsv(tsv_path)

        assert len(loaded) == len(sample_mapping)
        assert loaded.contig_to_accession == sample_mapping.contig_to_accession

    def test_from_tsv_validates_columns(self, temp_dir: Path):
        """Should raise error for invalid TSV format."""
        tsv_path = temp_dir / "invalid.tsv"
        tsv_path.write_text("wrong_col1\twrong_col2\nvalue1\tvalue2\n")

        with pytest.raises(ValueError, match="must have columns"):
            ContigIdMapping.from_tsv(tsv_path)

    def test_to_tsv_raises_for_empty_mapping(self, temp_dir: Path):
        """Should raise error when saving empty mapping."""
        empty = ContigIdMapping()

        with pytest.raises(ValueError, match="empty mapping"):
            empty.to_tsv(temp_dir / "empty.tsv")


class TestContigIdMappingContains:
    """Tests for __contains__ and __len__ methods."""

    def test_len(self):
        """Should return correct length."""
        mapping = ContigIdMapping(
            contig_to_accession={"a": "X", "b": "Y", "c": "Z"}
        )
        assert len(mapping) == 3

    def test_contains_present_key(self):
        """Should return True for present contig ID."""
        mapping = ContigIdMapping(contig_to_accession={"contig1": "GCF_123"})
        assert "contig1" in mapping

    def test_contains_absent_key(self):
        """Should return False for absent contig ID."""
        mapping = ContigIdMapping(contig_to_accession={"contig1": "GCF_123"})
        assert "contig2" not in mapping


class TestContigIdMappingEdgeCases:
    """Tests for edge cases and error handling."""

    def test_empty_mapping_transform(self):
        """Should handle empty mapping gracefully."""
        empty = ContigIdMapping()

        # transform_sseqid should pass through
        assert empty.transform_sseqid("any_id") == "any_id"

        # transform_column should return unchanged
        df = pl.DataFrame({"sseqid": ["id1", "id2"]})
        result = empty.transform_column(df, "sseqid")
        assert result["sseqid"].to_list() == ["id1", "id2"]

    def test_duplicate_contigs_last_wins(self, temp_dir: Path):
        """When duplicate contig IDs exist, last genome file should win."""
        # Create two genome files with overlapping contig ID
        (temp_dir / "GCF_000111111.1_genomic.fna").write_text(">contig1\nATGC\n")
        (temp_dir / "GCF_000222222.1_genomic.fna").write_text(">contig1\nGCAT\n")

        # This should log a warning but not fail
        mapping = ContigIdMapping.from_genome_dir(temp_dir)

        # Last one wins (alphabetical order means 222222 comes after 111111)
        assert len(mapping) == 1
        assert mapping.contig_to_accession["contig1"] == "GCF_000222222.1"
