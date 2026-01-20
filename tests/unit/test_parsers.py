"""
Unit tests for streaming parsers.

Tests StreamingBlastParser and ANIMatrixParser including:
- BLAST tabular file parsing
- Streaming and chunked reading
- ANI matrix validation
- Edge cases (empty files, malformed input)
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from metadarkmatter.core.parsers import (
    ANIMatrixParser,
    BlastHitFast,
    BlastResultFast,
    StreamingBlastParser,
    extract_genome_name_expr,
)


class TestExtractGenomeNameExpr:
    """Tests for vectorized genome name extraction.

    The extract_genome_name_expr() function expects standardized headers
    in the format {accession}|{contig_id} and extracts the accession part.
    """

    def test_standardized_pipe_format(self):
        """Should extract accession from standardized pipe format."""
        df = pl.DataFrame({"sseqid": ["GCF_000123456.1|NZ_CP012345.1"]})
        result = df.with_columns([extract_genome_name_expr()])

        assert result["genome_name"][0] == "GCF_000123456.1"

    def test_genbank_pipe_format(self):
        """Should extract GCA accession from pipe format."""
        df = pl.DataFrame({"sseqid": ["GCA_000789012.1|contig_001"]})
        result = df.with_columns([extract_genome_name_expr()])

        assert result["genome_name"][0] == "GCA_000789012.1"

    def test_complex_contig_id(self):
        """Should handle complex contig IDs after pipe."""
        df = pl.DataFrame({"sseqid": ["GCF_000195955.2|NZ_CP000439.1"]})
        result = df.with_columns([extract_genome_name_expr()])

        assert result["genome_name"][0] == "GCF_000195955.2"

    def test_no_pipe_returns_unknown(self):
        """Should return unknown for non-standardized formats."""
        df = pl.DataFrame({"sseqid": ["custom_contig_1"]})
        result = df.with_columns([extract_genome_name_expr()])

        assert result["genome_name"][0] == "unknown"

    def test_null_handling(self):
        """Should handle null values."""
        df = pl.DataFrame({"sseqid": [None, "GCF_000123456.1|contig1"]})
        result = df.with_columns([extract_genome_name_expr()])

        assert result["genome_name"][0] == "unknown"
        assert result["genome_name"][1] == "GCF_000123456.1"

    def test_multiple_standardized_formats(self):
        """Should handle multiple standardized headers."""
        df = pl.DataFrame({
            "sseqid": [
                "GCF_000123456.1|NZ_CP012345.1",
                "GCA_000789012.2|scaffold_1",
                "GCF_000195955.2|NZ_ABCD01000001.1",
                "GCA_001234567.1|contig_42",
            ]
        })
        result = df.with_columns([extract_genome_name_expr()])

        assert result["genome_name"].to_list() == [
            "GCF_000123456.1",
            "GCA_000789012.2",
            "GCF_000195955.2",
            "GCA_001234567.1",
        ]

    def test_legacy_format_returns_unknown(self):
        """Legacy formats without pipe return unknown."""
        df = pl.DataFrame({
            "sseqid": [
                "GCF_000123456.1_ASM123v1",
                "NZ_CP012345.1",
            ]
        })
        result = df.with_columns([extract_genome_name_expr()])

        # Legacy formats should return unknown since they don't have pipe
        assert result["genome_name"].to_list() == ["unknown", "unknown"]


class TestBlastHitFast:
    """Tests for lightweight BlastHitFast NamedTuple."""

    def test_create_hit_fast(self):
        """Should create BlastHitFast with all fields."""
        hit = BlastHitFast(
            qseqid="read_001",
            sseqid="GCF_000123456.1",
            pident=98.5,
            bitscore=250.0,
            genome_name="GCF_000123456.1",
        )

        assert hit.qseqid == "read_001"
        assert hit.pident == 98.5
        assert hit.bitscore == 250.0
        assert hit.genome_name == "GCF_000123456.1"

    def test_hit_fast_is_immutable(self):
        """BlastHitFast should be immutable (NamedTuple)."""
        hit = BlastHitFast("read_001", "GCF_000123456.1", 98.5, 250.0, "GCF_000123456.1")

        with pytest.raises(AttributeError):
            hit.pident = 50.0


class TestBlastResultFast:
    """Tests for lightweight BlastResultFast NamedTuple."""

    def test_create_result_fast(self):
        """Should create BlastResultFast with hits tuple."""
        hits = (
            BlastHitFast("read_001", "GCF_000123456.1", 98.5, 250.0, "GCF_000123456.1"),
            BlastHitFast("read_001", "GCF_000789012.1", 95.0, 240.0, "GCF_000789012.1"),
        )
        result = BlastResultFast(read_id="read_001", hits=hits)

        assert result.read_id == "read_001"
        assert len(result.hits) == 2

    def test_result_fast_empty_hits(self):
        """Should handle empty hits tuple."""
        result = BlastResultFast(read_id="read_001", hits=())

        assert result.read_id == "read_001"
        assert len(result.hits) == 0


class TestStreamingBlastParser:
    """Tests for StreamingBlastParser."""

    def test_parser_file_not_found(self, temp_dir):
        """Should raise FileNotFoundError for missing file."""
        with pytest.raises(FileNotFoundError):
            StreamingBlastParser(temp_dir / "nonexistent.blast.tsv")

    def test_parse_lazy_returns_lazyframe(self, temp_blast_file):
        """parse_lazy should return Polars LazyFrame."""
        parser = StreamingBlastParser(temp_blast_file)
        lazy_df = parser.parse_lazy()

        assert isinstance(lazy_df, pl.LazyFrame)

    def test_parse_eager_returns_dataframe(self, temp_blast_file):
        """parse_eager should return materialized DataFrame."""
        parser = StreamingBlastParser(temp_blast_file)
        df = parser.parse_eager()

        assert isinstance(df, pl.DataFrame)
        assert len(df) > 0

    def test_parse_eager_correct_columns(self, temp_blast_file):
        """Parsed DataFrame should have correct column names."""
        parser = StreamingBlastParser(temp_blast_file)
        df = parser.parse_eager()

        expected_columns = [
            "qseqid", "sseqid", "pident", "length", "mismatch",
            "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"
        ]
        assert df.columns == expected_columns

    def test_parse_eager_correct_types(self, temp_blast_file):
        """Parsed DataFrame should have correct column types."""
        parser = StreamingBlastParser(temp_blast_file)
        df = parser.parse_eager()

        assert df["qseqid"].dtype == pl.Utf8
        assert df["sseqid"].dtype == pl.Utf8
        assert df["pident"].dtype == pl.Float64
        assert df["length"].dtype == pl.Int64
        assert df["bitscore"].dtype == pl.Float64

    def test_iter_reads_yields_blast_results(self, temp_blast_file):
        """iter_reads should yield BlastResult objects."""
        from metadarkmatter.models.blast import BlastResult

        parser = StreamingBlastParser(temp_blast_file)
        results = list(parser.iter_reads())

        assert len(results) > 0
        assert all(isinstance(r, BlastResult) for r in results)

    def test_iter_reads_groups_by_read_id(self, temp_blast_file):
        """iter_reads should group hits by read_id."""
        parser = StreamingBlastParser(temp_blast_file)
        results = list(parser.iter_reads())

        # Each result should have a unique read_id
        read_ids = [r.read_id for r in results]
        assert len(read_ids) == len(set(read_ids))

    def test_iter_reads_hits_sorted_by_bitscore(self, temp_blast_file):
        """Hits within each result should be sorted by bitscore descending."""
        parser = StreamingBlastParser(temp_blast_file)

        for result in parser.iter_reads():
            if len(result.hits) > 1:
                bitscores = [h.bitscore for h in result.hits]
                assert bitscores == sorted(bitscores, reverse=True)

    def test_iter_reads_fast_yields_blast_result_fast(self, temp_blast_file):
        """iter_reads_fast should yield BlastResultFast objects."""
        parser = StreamingBlastParser(temp_blast_file)
        results = list(parser.iter_reads_fast())

        assert len(results) > 0
        assert all(isinstance(r, BlastResultFast) for r in results)

    def test_iter_reads_fast_has_genome_name(self, temp_blast_file):
        """iter_reads_fast should have pre-extracted genome names."""
        parser = StreamingBlastParser(temp_blast_file)

        for result in parser.iter_reads_fast():
            for hit in result.hits:
                assert hit.genome_name is not None
                assert hit.genome_name != ""

    def test_get_top_hits_per_read(self, temp_blast_file):
        """get_top_hits_per_read should return one row per read."""
        parser = StreamingBlastParser(temp_blast_file)
        top_hits = parser.get_top_hits_per_read()

        # Should have unique qseqid values
        assert len(top_hits) == top_hits["qseqid"].n_unique()

    def test_get_ambiguous_hits_default_threshold(self, temp_blast_file):
        """get_ambiguous_hits should filter by 95% of top bitscore."""
        parser = StreamingBlastParser(temp_blast_file)
        ambiguous = parser.get_ambiguous_hits(bitscore_threshold_pct=95.0)

        assert "max_bitscore" in ambiguous.columns
        # All rows should have bitscore >= 95% of max
        for row in ambiguous.iter_rows(named=True):
            assert row["bitscore"] >= row["max_bitscore"] * 0.95

    def test_count_hits_per_read(self, temp_blast_file):
        """count_hits_per_read should return hit counts."""
        parser = StreamingBlastParser(temp_blast_file)
        counts = parser.count_hits_per_read()

        assert "qseqid" in counts.columns
        assert "hit_count" in counts.columns
        assert all(counts["hit_count"] > 0)

    def test_get_summary_stats(self, temp_blast_file):
        """get_summary_stats should return summary dictionary."""
        parser = StreamingBlastParser(temp_blast_file)
        stats = parser.get_summary_stats()

        assert "total_hits" in stats
        assert "unique_reads" in stats
        assert "mean_pident" in stats
        assert "mean_bitscore" in stats
        assert stats["total_hits"] > 0

    def test_empty_file_handling(self, empty_blast_file):
        """Should raise NoDataError for empty BLAST file."""
        import polars

        parser = StreamingBlastParser(empty_blast_file)

        # Polars raises NoDataError for empty CSV files
        with pytest.raises(polars.exceptions.NoDataError):
            list(parser.iter_reads())

    def test_custom_chunk_size(self, temp_blast_file):
        """Should respect custom chunk_size parameter."""
        # Use minimum valid chunk_size (1000)
        parser = StreamingBlastParser(
            temp_blast_file, chunk_size=StreamingBlastParser.MIN_CHUNK_SIZE
        )

        # Should still parse correctly with minimum chunks
        results = list(parser.iter_reads())
        assert len(results) > 0

    def test_chunk_size_too_small_raises(self, temp_blast_file):
        """Should raise ValueError for chunk_size below minimum."""
        with pytest.raises(ValueError, match="chunk_size must be between"):
            StreamingBlastParser(
                temp_blast_file,
                chunk_size=StreamingBlastParser.MIN_CHUNK_SIZE - 1,
            )

    def test_chunk_size_too_large_raises(self, temp_blast_file):
        """Should raise ValueError for chunk_size above maximum."""
        with pytest.raises(ValueError, match="chunk_size must be between"):
            StreamingBlastParser(
                temp_blast_file,
                chunk_size=StreamingBlastParser.MAX_CHUNK_SIZE + 1,
            )


class TestANIMatrixParser:
    """Tests for ANIMatrixParser."""

    def test_parser_file_not_found(self, temp_dir):
        """Should raise FileNotFoundError for missing file."""
        with pytest.raises(FileNotFoundError):
            ANIMatrixParser(temp_dir / "nonexistent.ani.csv")

    def test_parse_returns_dataframe(self, temp_ani_file):
        """parse should return Polars DataFrame."""
        parser = ANIMatrixParser(temp_ani_file)
        df = parser.parse()

        assert isinstance(df, pl.DataFrame)

    def test_parse_correct_shape(self, temp_ani_file):
        """Parsed matrix should be square."""
        parser = ANIMatrixParser(temp_ani_file)
        df = parser.parse()

        # First column is genome names, rest are values
        num_genomes = len(df)
        num_value_cols = len(df.columns) - 1
        assert num_genomes == num_value_cols

    def test_to_dict_returns_nested_dict(self, temp_ani_file):
        """to_dict should return nested dictionary."""
        parser = ANIMatrixParser(temp_ani_file)
        ani_dict = parser.to_dict()

        assert isinstance(ani_dict, dict)
        # Each value should be a dict
        for genome, inner in ani_dict.items():
            assert isinstance(inner, dict)

    def test_to_dict_diagonal_is_100(self, temp_ani_file):
        """Diagonal values (self-ANI) should be 100."""
        parser = ANIMatrixParser(temp_ani_file)
        ani_dict = parser.to_dict()

        for genome, inner in ani_dict.items():
            assert inner[genome] == 100.0

    def test_to_dict_symmetric(self, temp_ani_file):
        """Matrix should be symmetric."""
        parser = ANIMatrixParser(temp_ani_file)
        ani_dict = parser.to_dict()

        genomes = list(ani_dict.keys())
        for i, g1 in enumerate(genomes):
            for g2 in genomes[i + 1:]:
                assert ani_dict[g1][g2] == ani_dict[g2][g1]

    def test_validate_non_square_matrix(self, temp_dir):
        """Should reject non-square matrix."""
        # Create non-square matrix
        ani_path = temp_dir / "non_square.csv"
        df = pl.DataFrame({
            "genome": ["A", "B"],
            "A": [100.0, 95.0],
            "B": [95.0, 100.0],
            "C": [80.0, 85.0],  # Extra column
        })
        df.write_csv(ani_path)

        parser = ANIMatrixParser(ani_path)
        with pytest.raises(ValueError) as exc_info:
            parser.parse()

        assert "not square" in str(exc_info.value)

    def test_validate_row_column_mismatch(self, temp_dir):
        """Should reject matrix with mismatched row/column names."""
        ani_path = temp_dir / "mismatch.csv"
        df = pl.DataFrame({
            "genome": ["A", "B"],
            "X": [100.0, 95.0],  # Should be "A"
            "Y": [95.0, 100.0],  # Should be "B"
        })
        df.write_csv(ani_path)

        parser = ANIMatrixParser(ani_path)
        with pytest.raises(ValueError) as exc_info:
            parser.parse()

        assert "mismatch" in str(exc_info.value).lower()

    def test_validate_out_of_range_values(self, temp_dir):
        """Should reject ANI values outside 0-100 range."""
        ani_path = temp_dir / "out_of_range.csv"
        df = pl.DataFrame({
            "genome": ["A", "B"],
            "A": [100.0, 150.0],  # Invalid: 150 > 100
            "B": [150.0, 100.0],
        })
        df.write_csv(ani_path)

        parser = ANIMatrixParser(ani_path)
        with pytest.raises(ValueError) as exc_info:
            parser.parse()

        assert "out of range" in str(exc_info.value).lower()

    def test_validate_negative_values(self, temp_dir):
        """Should reject negative ANI values."""
        ani_path = temp_dir / "negative.csv"
        df = pl.DataFrame({
            "genome": ["A", "B"],
            "A": [100.0, -5.0],  # Invalid: negative
            "B": [-5.0, 100.0],
        })
        df.write_csv(ani_path)

        parser = ANIMatrixParser(ani_path)
        with pytest.raises(ValueError) as exc_info:
            parser.parse()

        assert "out of range" in str(exc_info.value).lower()

    def test_validate_empty_matrix(self, temp_dir):
        """Should reject empty matrix."""
        ani_path = temp_dir / "empty.csv"
        ani_path.write_text("genome\n")  # Header only

        parser = ANIMatrixParser(ani_path)
        with pytest.raises(ValueError) as exc_info:
            parser.parse()

        assert "empty" in str(exc_info.value).lower()

    def test_tsv_format_detection(self, temp_dir):
        """Should correctly parse TSV format."""
        ani_path = temp_dir / "matrix.tsv"
        df = pl.DataFrame({
            "genome": ["A", "B"],
            "A": [100.0, 95.0],
            "B": [95.0, 100.0],
        })
        df.write_csv(ani_path, separator="\t")

        parser = ANIMatrixParser(ani_path)
        ani_dict = parser.to_dict()

        assert ani_dict["A"]["B"] == 95.0


class TestStreamingBlastParserChunked:
    """Tests for chunked streaming behavior."""

    def test_iter_reads_chunked_same_as_iter_reads(self, temp_blast_file):
        """iter_reads_chunked should produce same results as iter_reads."""
        # Use minimum valid chunk_size
        parser = StreamingBlastParser(
            temp_blast_file, chunk_size=StreamingBlastParser.MIN_CHUNK_SIZE
        )

        chunked_results = list(parser.iter_reads_chunked())
        parser2 = StreamingBlastParser(temp_blast_file)
        regular_results = list(parser2.iter_reads())

        # Same number of results
        assert len(chunked_results) == len(regular_results)

        # Same read IDs
        chunked_ids = {r.read_id for r in chunked_results}
        regular_ids = {r.read_id for r in regular_results}
        assert chunked_ids == regular_ids

    def test_iter_reads_fast_with_small_chunks(self, temp_blast_file):
        """iter_reads_fast should work with small chunk sizes."""
        # Use minimum valid chunk_size
        parser = StreamingBlastParser(
            temp_blast_file, chunk_size=StreamingBlastParser.MIN_CHUNK_SIZE
        )

        results = list(parser.iter_reads_fast())
        assert len(results) > 0

        # Each result should have valid genome names
        for result in results:
            for hit in result.hits:
                assert hit.genome_name is not None

    def test_boundary_handling_across_chunks(self, temp_dir):
        """Reads spanning chunk boundaries should be handled correctly."""
        # Create file with many reads to test chunking with MIN_CHUNK_SIZE
        blast_path = temp_dir / "boundary_test.blast.tsv"
        rows = []
        num_reads = 50  # Enough reads for chunking to matter
        hits_per_read = 5

        for read_idx in range(num_reads):
            for hit_idx in range(hits_per_read):
                rows.append({
                    "qseqid": f"read_{read_idx:03d}",
                    "sseqid": f"GCF_{hit_idx:09d}.1",
                    "pident": 95.0,
                    "length": 150,
                    "mismatch": 7,
                    "gapopen": 0,
                    "qstart": 1,
                    "qend": 150,
                    "sstart": 1000,
                    "send": 1150,
                    "evalue": 1e-50,
                    "bitscore": 250.0 - hit_idx,
                })

        df = pl.DataFrame(rows)
        df.write_csv(blast_path, separator="\t", include_header=False)

        # Parse with minimum chunk size
        parser = StreamingBlastParser(
            blast_path, chunk_size=StreamingBlastParser.MIN_CHUNK_SIZE
        )
        results = list(parser.iter_reads())

        assert len(results) == num_reads
        for result in results:
            assert result.num_hits == hits_per_read
