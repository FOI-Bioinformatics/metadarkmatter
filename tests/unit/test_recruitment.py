"""
Unit tests for recruitment data extraction module.

Tests genome name extraction, CIGAR-based identity calculation,
SAM line parsing, BAM streaming, DataFrame loading, and aggregation.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import polars as pl
import pytest

from metadarkmatter.core.recruitment import (
    AlignmentRecord,
    aggregate_by_genome,
    calculate_percent_identity,
    extract_genome_name,
    get_contig_lengths,
    load_recruitment_data,
    parse_sam_line,
    stream_bam_alignments,
)


# =============================================================================
# AlignmentRecord dataclass tests
# =============================================================================


class TestAlignmentRecord:
    """Tests for the AlignmentRecord frozen dataclass."""

    def test_create_record(self):
        """Should create a valid AlignmentRecord with all fields."""
        record = AlignmentRecord(
            read_id="read_001",
            genome_name="GCF_000123456.1",
            contig="GCF_000123456.1_scaffold_1",
            position=1000,
            mapq=60,
            cigar="100M",
            percent_identity=98.5,
            alignment_length=100,
        )
        assert record.read_id == "read_001"
        assert record.genome_name == "GCF_000123456.1"
        assert record.contig == "GCF_000123456.1_scaffold_1"
        assert record.position == 1000
        assert record.mapq == 60
        assert record.cigar == "100M"
        assert record.percent_identity == 98.5
        assert record.alignment_length == 100

    def test_record_is_frozen(self):
        """AlignmentRecord should be immutable."""
        record = AlignmentRecord(
            read_id="read_001",
            genome_name="genome",
            contig="contig_1",
            position=1,
            mapq=30,
            cigar="50M",
            percent_identity=95.0,
            alignment_length=50,
        )
        with pytest.raises(AttributeError):
            record.read_id = "different_read"


# =============================================================================
# extract_genome_name tests
# =============================================================================


class TestExtractGenomeName:
    """Tests for genome name extraction from contig headers."""

    def test_ncbi_gcf_accession(self):
        """Should extract accession from GCF_ prefixed contig."""
        result = extract_genome_name("GCF_000123456.1_ASM123v1_scaffold_1")
        assert result == "GCF_000123456.1"

    def test_ncbi_gca_accession(self):
        """Should extract accession from GCA_ prefixed contig."""
        result = extract_genome_name("GCA_000789012.2_ASM789v2_contig_5")
        assert result == "GCA_000789012.2"

    def test_gcf_accession_only(self):
        """Should handle bare GCF accession."""
        result = extract_genome_name("GCF_000111222.3")
        assert result == "GCF_000111222.3"

    def test_generic_format_with_digit_second_part(self):
        """Should return first part when second part is a digit."""
        result = extract_genome_name("genome_1")
        assert result == "genome"

    def test_generic_format_with_scaffold_keyword(self):
        """Should return first part when second part starts with scaffold."""
        result = extract_genome_name("mygenome_scaffold_3")
        assert result == "mygenome"

    def test_generic_format_with_contig_keyword(self):
        """Should return first part when second part starts with contig."""
        result = extract_genome_name("mygenome_contig_7")
        assert result == "mygenome"

    def test_generic_format_with_version_prefix(self):
        """Should return first part when second part starts with v."""
        result = extract_genome_name("assembly_v2_contig_1")
        assert result == "assembly"

    def test_generic_format_species_name(self):
        """Should join first two parts for species-like names."""
        result = extract_genome_name("Francisella_tularensis_contig1")
        assert result == "Francisella_tularensis"

    def test_single_part_contig(self):
        """Should return the contig as-is when no underscore present."""
        result = extract_genome_name("simplecontig")
        assert result == "simplecontig"

    def test_gcf_no_match_falls_through(self):
        """GCF_ prefix that does not match accession pattern should fall through."""
        # GCF_notanumber does not match the regex for accession
        result = extract_genome_name("GCF_notanumber_stuff")
        assert result == "GCF_notanumber"


# =============================================================================
# calculate_percent_identity tests
# =============================================================================


class TestCalculatePercentIdentity:
    """Tests for CIGAR-based percent identity calculation."""

    def test_perfect_match(self):
        """All M operations should yield 100% identity."""
        assert calculate_percent_identity("100M") == 100.0

    def test_exact_match_operator(self):
        """Sequence match (=) should count as matches."""
        assert calculate_percent_identity("80=20X") == 80.0

    def test_with_insertions(self):
        """Insertions should reduce identity."""
        # 90M + 10I => matches=90, total=100
        result = calculate_percent_identity("90M10I")
        assert result == pytest.approx(90.0)

    def test_with_deletions(self):
        """Deletions should reduce identity."""
        # 80M + 20D => matches=80, total=100
        result = calculate_percent_identity("80M20D")
        assert result == pytest.approx(80.0)

    def test_with_mismatches(self):
        """X operations should not count as matches."""
        # 75M + 25X => matches=75, total=100
        result = calculate_percent_identity("75M25X")
        assert result == pytest.approx(75.0)

    def test_complex_cigar(self):
        """Should handle a realistic multi-operation CIGAR."""
        # 50M2I30M5D15M => matches=50+30+15=95, total=50+2+30+5+15=102
        result = calculate_percent_identity("50M2I30M5D15M")
        assert result == pytest.approx(100.0 * 95 / 102)

    def test_soft_clipping_ignored(self):
        """Soft clipped bases (S) should not affect identity."""
        # 10S80M10S => matches=80, total=80
        result = calculate_percent_identity("10S80M10S")
        assert result == 100.0

    def test_hard_clipping_ignored(self):
        """Hard clipped bases (H) should not affect identity."""
        # 5H90M5H => matches=90, total=90
        result = calculate_percent_identity("5H90M5H")
        assert result == 100.0

    def test_empty_cigar(self):
        """Empty CIGAR string should return 0.0."""
        assert calculate_percent_identity("") == 0.0

    def test_only_soft_clip(self):
        """CIGAR with only soft clipping should return 0.0."""
        assert calculate_percent_identity("100S") == 0.0

    def test_padding_ignored(self):
        """Padding (P) and skipped regions (N) should not affect identity."""
        # 50M100N50M => matches=100, total=100
        result = calculate_percent_identity("50M100N50M")
        assert result == 100.0


# =============================================================================
# parse_sam_line tests
# =============================================================================


class TestParseSamLine:
    """Tests for SAM line parsing."""

    def _make_sam_line(
        self,
        read_id: str = "read_001",
        flag: int = 0,
        contig: str = "GCF_000123456.1_scaffold_1",
        pos: int = 1000,
        mapq: int = 60,
        cigar: str = "100M",
        rnext: str = "*",
        pnext: int = 0,
        tlen: int = 0,
        seq: str = "A" * 100,
        qual: str = "I" * 100,
    ) -> str:
        """Build a SAM-formatted line from individual fields."""
        return "\t".join([
            read_id,
            str(flag),
            contig,
            str(pos),
            str(mapq),
            cigar,
            rnext,
            str(pnext),
            str(tlen),
            seq,
            qual,
        ])

    def test_parse_valid_line(self):
        """Should parse a valid mapped SAM line."""
        line = self._make_sam_line()
        record = parse_sam_line(line)

        assert record is not None
        assert record.read_id == "read_001"
        assert record.contig == "GCF_000123456.1_scaffold_1"
        assert record.genome_name == "GCF_000123456.1"
        assert record.position == 1000
        assert record.mapq == 60
        assert record.cigar == "100M"
        assert record.percent_identity == 100.0
        assert record.alignment_length == 100

    def test_skip_header_line(self):
        """Should return None for SAM header lines."""
        assert parse_sam_line("@HD\tVN:1.6\tSO:coordinate") is None

    def test_skip_sq_header(self):
        """Should return None for @SQ header lines."""
        assert parse_sam_line("@SQ\tSN:chr1\tLN:248956422") is None

    def test_skip_short_line(self):
        """Should return None for lines with fewer than 11 fields."""
        assert parse_sam_line("read_001\t0\tchr1") is None

    def test_skip_unmapped_read(self):
        """Should return None for unmapped reads (flag & 4)."""
        line = self._make_sam_line(flag=4, contig="*", cigar="*")
        assert parse_sam_line(line) is None

    def test_skip_star_cigar(self):
        """Should return None when CIGAR is *."""
        line = self._make_sam_line(cigar="*")
        assert parse_sam_line(line) is None

    def test_alignment_length_includes_indels(self):
        """Alignment length should include M, =, X, I, D operations."""
        # 50M5I30M10D10M => 50+5+30+10+10 = 105
        line = self._make_sam_line(cigar="50M5I30M10D10M", seq="A" * 95)
        record = parse_sam_line(line)
        assert record is not None
        assert record.alignment_length == 105

    def test_secondary_alignment_parsed(self):
        """Secondary alignments (flag 256) should be parsed if mapped."""
        line = self._make_sam_line(flag=256)
        record = parse_sam_line(line)
        assert record is not None

    def test_supplementary_and_unmapped_skipped(self):
        """Supplementary + unmapped (flag 4 set) should be skipped."""
        line = self._make_sam_line(flag=2052)  # 4 + 2048
        assert parse_sam_line(line) is None

    def test_trailing_whitespace_handled(self):
        """Trailing whitespace on a SAM line should be handled."""
        line = self._make_sam_line() + "  \n"
        record = parse_sam_line(line)
        assert record is not None
        assert record.read_id == "read_001"


# =============================================================================
# stream_bam_alignments tests
# =============================================================================


class TestStreamBamAlignments:
    """Tests for BAM streaming via samtools."""

    def test_stream_yields_records(self, tmp_path: Path):
        """Should yield AlignmentRecord objects from samtools output."""
        bam_path = tmp_path / "test.bam"
        bam_path.touch()

        sam_output = (
            "read_001\t0\tGCF_000123456.1_scaffold_1\t100\t60\t100M\t*\t0\t0\t"
            + "A" * 100
            + "\t"
            + "I" * 100
            + "\n"
            "read_002\t0\tGCF_000123456.1_scaffold_2\t200\t30\t50M\t*\t0\t0\t"
            + "A" * 50
            + "\t"
            + "I" * 50
            + "\n"
        )

        mock_process = MagicMock()
        mock_process.stdout = iter(sam_output.splitlines(keepends=True))
        mock_process.wait.return_value = 0

        with patch("metadarkmatter.core.recruitment.subprocess.Popen", return_value=mock_process):
            records = list(stream_bam_alignments(bam_path))

        assert len(records) == 2
        assert records[0].read_id == "read_001"
        assert records[1].read_id == "read_002"

    def test_stream_filters_by_identity(self, tmp_path: Path):
        """Should filter records below minimum identity."""
        bam_path = tmp_path / "test.bam"
        bam_path.touch()

        # 80M20X => 80% identity
        sam_output = (
            "read_low\t0\tcontig_1\t100\t60\t80M20X\t*\t0\t0\t"
            + "A" * 100
            + "\t"
            + "I" * 100
            + "\n"
            "read_high\t0\tcontig_2\t200\t60\t100M\t*\t0\t0\t"
            + "A" * 100
            + "\t"
            + "I" * 100
            + "\n"
        )

        mock_process = MagicMock()
        mock_process.stdout = iter(sam_output.splitlines(keepends=True))
        mock_process.wait.return_value = 0

        with patch("metadarkmatter.core.recruitment.subprocess.Popen", return_value=mock_process):
            records = list(stream_bam_alignments(bam_path, min_identity=90.0))

        assert len(records) == 1
        assert records[0].read_id == "read_high"

    def test_stream_passes_mapq_filter(self, tmp_path: Path):
        """Should pass min_mapq to samtools via -q flag."""
        bam_path = tmp_path / "test.bam"
        bam_path.touch()

        mock_process = MagicMock()
        mock_process.stdout = iter([])
        mock_process.wait.return_value = 0

        with patch("metadarkmatter.core.recruitment.subprocess.Popen", return_value=mock_process) as mock_popen:
            list(stream_bam_alignments(bam_path, min_mapq=20))

        called_cmd = mock_popen.call_args[0][0]
        assert "-q" in called_cmd
        assert "20" in called_cmd

    def test_stream_no_mapq_flag_when_zero(self, tmp_path: Path):
        """Should omit -q flag when min_mapq is 0."""
        bam_path = tmp_path / "test.bam"
        bam_path.touch()

        mock_process = MagicMock()
        mock_process.stdout = iter([])
        mock_process.wait.return_value = 0

        with patch("metadarkmatter.core.recruitment.subprocess.Popen", return_value=mock_process) as mock_popen:
            list(stream_bam_alignments(bam_path, min_mapq=0))

        called_cmd = mock_popen.call_args[0][0]
        assert "-q" not in called_cmd

    def test_stream_returns_empty_when_stdout_none(self, tmp_path: Path):
        """Should return no records if process.stdout is None."""
        bam_path = tmp_path / "test.bam"
        bam_path.touch()

        mock_process = MagicMock()
        mock_process.stdout = None

        with patch("metadarkmatter.core.recruitment.subprocess.Popen", return_value=mock_process):
            records = list(stream_bam_alignments(bam_path))

        assert len(records) == 0

    def test_stream_skips_header_and_unmapped(self, tmp_path: Path):
        """Should skip header lines and unmapped reads in stream."""
        bam_path = tmp_path / "test.bam"
        bam_path.touch()

        sam_output = (
            "@HD\tVN:1.6\n"
            "@SQ\tSN:contig_1\tLN:5000\n"
            "unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"
            "mapped\t0\tcontig_1\t100\t60\t50M\t*\t0\t0\t"
            + "A" * 50
            + "\t"
            + "I" * 50
            + "\n"
        )

        mock_process = MagicMock()
        mock_process.stdout = iter(sam_output.splitlines(keepends=True))
        mock_process.wait.return_value = 0

        with patch("metadarkmatter.core.recruitment.subprocess.Popen", return_value=mock_process):
            records = list(stream_bam_alignments(bam_path))

        assert len(records) == 1
        assert records[0].read_id == "mapped"

    def test_stream_file_not_found(self, tmp_path: Path):
        """Should raise FileNotFoundError for non-existent BAM."""
        bam_path = tmp_path / "nonexistent.bam"

        with pytest.raises(FileNotFoundError):
            list(stream_bam_alignments(bam_path))


# =============================================================================
# load_recruitment_data tests
# =============================================================================


class TestLoadRecruitmentData:
    """Tests for loading recruitment data into a Polars DataFrame."""

    def test_load_returns_dataframe(self, tmp_path: Path):
        """Should return a Polars DataFrame with expected columns."""
        bam_path = tmp_path / "test.bam"
        bam_path.touch()

        sam_output = (
            "read_001\t0\tGCF_000123456.1_scaffold_1\t100\t60\t100M\t*\t0\t0\t"
            + "A" * 100
            + "\t"
            + "I" * 100
            + "\n"
        )

        mock_process = MagicMock()
        mock_process.stdout = iter(sam_output.splitlines(keepends=True))
        mock_process.wait.return_value = 0

        with patch("metadarkmatter.core.recruitment.subprocess.Popen", return_value=mock_process):
            df = load_recruitment_data(bam_path)

        assert isinstance(df, pl.DataFrame)
        expected_cols = {
            "read_id",
            "genome_name",
            "contig",
            "position",
            "mapq",
            "percent_identity",
            "alignment_length",
        }
        assert set(df.columns) == expected_cols
        assert len(df) == 1

    def test_load_respects_max_records(self, tmp_path: Path):
        """Should stop loading after max_records."""
        bam_path = tmp_path / "test.bam"
        bam_path.touch()

        lines = []
        for i in range(10):
            lines.append(
                f"read_{i:03d}\t0\tcontig_1\t{100 + i}\t60\t50M\t*\t0\t0\t"
                + "A" * 50
                + "\t"
                + "I" * 50
                + "\n"
            )

        mock_process = MagicMock()
        mock_process.stdout = iter(lines)
        mock_process.wait.return_value = 0

        with patch("metadarkmatter.core.recruitment.subprocess.Popen", return_value=mock_process):
            df = load_recruitment_data(bam_path, max_records=3)

        assert len(df) == 3

    def test_load_empty_bam(self, tmp_path: Path):
        """Should return empty DataFrame for BAM with no alignments."""
        bam_path = tmp_path / "test.bam"
        bam_path.touch()

        mock_process = MagicMock()
        mock_process.stdout = iter([])
        mock_process.wait.return_value = 0

        with patch("metadarkmatter.core.recruitment.subprocess.Popen", return_value=mock_process):
            df = load_recruitment_data(bam_path)

        assert isinstance(df, pl.DataFrame)
        assert len(df) == 0

    def test_load_passes_identity_filter(self, tmp_path: Path):
        """Should pass min_identity to stream_bam_alignments."""
        bam_path = tmp_path / "test.bam"
        bam_path.touch()

        # 80M20X => 80% identity; 100M => 100% identity
        sam_output = (
            "read_low\t0\tcontig_1\t100\t60\t80M20X\t*\t0\t0\t"
            + "A" * 100
            + "\t"
            + "I" * 100
            + "\n"
            "read_high\t0\tcontig_2\t200\t60\t100M\t*\t0\t0\t"
            + "A" * 100
            + "\t"
            + "I" * 100
            + "\n"
        )

        mock_process = MagicMock()
        mock_process.stdout = iter(sam_output.splitlines(keepends=True))
        mock_process.wait.return_value = 0

        with patch("metadarkmatter.core.recruitment.subprocess.Popen", return_value=mock_process):
            df = load_recruitment_data(bam_path, min_identity=90.0)

        assert len(df) == 1
        assert df["read_id"][0] == "read_high"

    def test_load_dataframe_values(self, tmp_path: Path):
        """Should populate DataFrame values correctly."""
        bam_path = tmp_path / "test.bam"
        bam_path.touch()

        sam_output = (
            "read_001\t0\tGCF_000123456.1_scaffold_1\t500\t42\t80M20I\t*\t0\t0\t"
            + "A" * 100
            + "\t"
            + "I" * 100
            + "\n"
        )

        mock_process = MagicMock()
        mock_process.stdout = iter(sam_output.splitlines(keepends=True))
        mock_process.wait.return_value = 0

        with patch("metadarkmatter.core.recruitment.subprocess.Popen", return_value=mock_process):
            df = load_recruitment_data(bam_path, min_identity=0.0)

        assert df["read_id"][0] == "read_001"
        assert df["genome_name"][0] == "GCF_000123456.1"
        assert df["contig"][0] == "GCF_000123456.1_scaffold_1"
        assert df["position"][0] == 500
        assert df["mapq"][0] == 42
        assert df["alignment_length"][0] == 100  # 80 + 20
        assert df["percent_identity"][0] == pytest.approx(80.0)


# =============================================================================
# get_contig_lengths tests
# =============================================================================


class TestGetContigLengths:
    """Tests for extracting contig lengths from BAM header."""

    def test_parse_sq_lines(self, tmp_path: Path):
        """Should extract contig names and lengths from @SQ lines."""
        bam_path = tmp_path / "test.bam"
        bam_path.touch()

        header = (
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:contig_1\tLN:5000\n"
            "@SQ\tSN:contig_2\tLN:3200\n"
            "@SQ\tSN:contig_3\tLN:10500\n"
        )

        mock_result = MagicMock()
        mock_result.stdout = header

        with patch("metadarkmatter.core.recruitment.subprocess.run", return_value=mock_result):
            lengths = get_contig_lengths(bam_path)

        assert lengths == {
            "contig_1": 5000,
            "contig_2": 3200,
            "contig_3": 10500,
        }

    def test_parse_empty_header(self, tmp_path: Path):
        """Should return empty dict for a header with no @SQ lines."""
        bam_path = tmp_path / "test.bam"
        bam_path.touch()

        header = "@HD\tVN:1.6\tSO:coordinate\n"

        mock_result = MagicMock()
        mock_result.stdout = header

        with patch("metadarkmatter.core.recruitment.subprocess.run", return_value=mock_result):
            lengths = get_contig_lengths(bam_path)

        assert lengths == {}

    def test_skips_non_sq_lines(self, tmp_path: Path):
        """Should skip non-@SQ header lines and blank lines."""
        bam_path = tmp_path / "test.bam"
        bam_path.touch()

        header = (
            "@HD\tVN:1.6\n"
            "@RG\tID:sample1\tSM:sample1\n"
            "@SQ\tSN:chr1\tLN:248956422\n"
            "@PG\tID:bwa\tPN:bwa\n"
            "\n"
        )

        mock_result = MagicMock()
        mock_result.stdout = header

        with patch("metadarkmatter.core.recruitment.subprocess.run", return_value=mock_result):
            lengths = get_contig_lengths(bam_path)

        assert lengths == {"chr1": 248956422}

    def test_file_not_found(self, tmp_path: Path):
        """Should raise FileNotFoundError for non-existent BAM."""
        bam_path = tmp_path / "nonexistent.bam"

        with pytest.raises(FileNotFoundError):
            get_contig_lengths(bam_path)


# =============================================================================
# aggregate_by_genome tests
# =============================================================================


class TestAggregateByGenome:
    """Tests for per-genome aggregation of recruitment data."""

    @pytest.fixture
    def recruitment_df(self) -> pl.DataFrame:
        """Sample recruitment DataFrame with multiple genomes."""
        return pl.DataFrame({
            "read_id": [
                "r1", "r2", "r3", "r4", "r5",
                "r6", "r7", "r8",
            ],
            "genome_name": [
                "genome_A", "genome_A", "genome_A", "genome_A", "genome_A",
                "genome_B", "genome_B", "genome_B",
            ],
            "contig": [
                "c1", "c1", "c2", "c2", "c3",
                "c4", "c4", "c5",
            ],
            "position": [100, 200, 300, 400, 500, 600, 700, 800],
            "mapq": [60, 50, 40, 30, 20, 60, 55, 45],
            "percent_identity": [
                99.0, 97.0, 95.0, 93.0, 91.0,
                88.0, 86.0, 84.0,
            ],
            "alignment_length": [100, 100, 100, 100, 100, 100, 100, 100],
        })

    def test_aggregate_columns(self, recruitment_df: pl.DataFrame):
        """Should produce expected output columns."""
        result = aggregate_by_genome(recruitment_df)
        expected_cols = {
            "genome_name",
            "num_reads",
            "mean_identity",
            "median_identity",
            "min_identity",
            "max_identity",
        }
        assert set(result.columns) == expected_cols

    def test_aggregate_num_reads(self, recruitment_df: pl.DataFrame):
        """Should count reads per genome correctly."""
        result = aggregate_by_genome(recruitment_df)

        genome_a = result.filter(pl.col("genome_name") == "genome_A")
        genome_b = result.filter(pl.col("genome_name") == "genome_B")

        assert genome_a["num_reads"][0] == 5
        assert genome_b["num_reads"][0] == 3

    def test_aggregate_mean_identity(self, recruitment_df: pl.DataFrame):
        """Should calculate mean identity per genome."""
        result = aggregate_by_genome(recruitment_df)

        genome_a = result.filter(pl.col("genome_name") == "genome_A")
        expected_mean_a = (99.0 + 97.0 + 95.0 + 93.0 + 91.0) / 5.0
        assert genome_a["mean_identity"][0] == pytest.approx(expected_mean_a)

    def test_aggregate_min_max_identity(self, recruitment_df: pl.DataFrame):
        """Should track min and max identity per genome."""
        result = aggregate_by_genome(recruitment_df)

        genome_b = result.filter(pl.col("genome_name") == "genome_B")
        assert genome_b["min_identity"][0] == pytest.approx(84.0)
        assert genome_b["max_identity"][0] == pytest.approx(88.0)

    def test_aggregate_sorted_by_reads_descending(self, recruitment_df: pl.DataFrame):
        """Should sort output by num_reads descending."""
        result = aggregate_by_genome(recruitment_df)
        assert result["genome_name"][0] == "genome_A"
        assert result["genome_name"][1] == "genome_B"

    def test_aggregate_single_genome(self):
        """Should handle a single genome correctly."""
        df = pl.DataFrame({
            "read_id": ["r1", "r2"],
            "genome_name": ["only_genome", "only_genome"],
            "contig": ["c1", "c1"],
            "position": [100, 200],
            "mapq": [60, 60],
            "percent_identity": [95.0, 85.0],
            "alignment_length": [100, 100],
        })

        result = aggregate_by_genome(df)
        assert len(result) == 1
        assert result["num_reads"][0] == 2
        assert result["mean_identity"][0] == pytest.approx(90.0)
        assert result["median_identity"][0] == pytest.approx(90.0)
        assert result["min_identity"][0] == pytest.approx(85.0)
        assert result["max_identity"][0] == pytest.approx(95.0)

    def test_aggregate_empty_dataframe(self):
        """Should handle empty DataFrame gracefully."""
        df = pl.DataFrame({
            "read_id": [],
            "genome_name": [],
            "contig": [],
            "position": [],
            "mapq": [],
            "percent_identity": [],
            "alignment_length": [],
        }).cast({
            "read_id": pl.Utf8,
            "genome_name": pl.Utf8,
            "contig": pl.Utf8,
            "position": pl.Int64,
            "mapq": pl.Int64,
            "percent_identity": pl.Float64,
            "alignment_length": pl.Int64,
        })

        result = aggregate_by_genome(df)
        assert len(result) == 0
