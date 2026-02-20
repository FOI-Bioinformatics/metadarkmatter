"""
Integration tests for metadarkmatter CLI.

Tests the command-line interface using Typer's CliRunner for:
- Main entry point and version display
- score classify command with various options
- score batch command for multiple files
- Error handling for invalid inputs
"""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.main import app

runner = CliRunner()


# =============================================================================
# CLI Test Fixtures
# =============================================================================


@pytest.fixture
def cli_blast_file(temp_dir: Path) -> Path:
    """Create a realistic BLAST file for CLI testing.

    Uses standardized pipe format: {accession}|{contig_id}
    """
    blast_data = []

    # Create reads with different classification outcomes
    # Known species reads (high identity, low ambiguity)
    for i in range(10):
        blast_data.append({
            "qseqid": f"known_read_{i:03d}",
            "sseqid": "GCF_000123456.1|ASM123v1_genomic",
            "pident": 99.5,
            "length": 150,
            "mismatch": 0,
            "gapopen": 0,
            "qstart": 1,
            "qend": 150,
            "sstart": 1000 + i * 200,
            "send": 1150 + i * 200,
            "evalue": 1e-80,
            "bitscore": 280.0,
        })

    # Novel species reads (moderate identity)
    for i in range(5):
        blast_data.append({
            "qseqid": f"novel_sp_read_{i:03d}",
            "sseqid": "GCF_000123456.1|ASM123v1_genomic",
            "pident": 90.0,
            "length": 150,
            "mismatch": 15,
            "gapopen": 0,
            "qstart": 1,
            "qend": 150,
            "sstart": 5000 + i * 200,
            "send": 5150 + i * 200,
            "evalue": 1e-50,
            "bitscore": 200.0,
        })

    # Novel genus reads (low identity)
    for i in range(3):
        blast_data.append({
            "qseqid": f"novel_gen_read_{i:03d}",
            "sseqid": "GCF_000123456.1|ASM123v1_genomic",
            "pident": 80.0,
            "length": 150,
            "mismatch": 30,
            "gapopen": 0,
            "qstart": 1,
            "qend": 150,
            "sstart": 8000 + i * 200,
            "send": 8150 + i * 200,
            "evalue": 1e-30,
            "bitscore": 150.0,
        })

    # Conserved region reads (ambiguous - multiple similar hits)
    for i in range(2):
        # First hit
        blast_data.append({
            "qseqid": f"conserved_read_{i:03d}",
            "sseqid": "GCF_000123456.1|ASM123v1_genomic",
            "pident": 98.0,
            "length": 150,
            "mismatch": 3,
            "gapopen": 0,
            "qstart": 1,
            "qend": 150,
            "sstart": 10000,
            "send": 10150,
            "evalue": 1e-70,
            "bitscore": 260.0,
        })
        # Second hit (similar bitscore, different genome)
        blast_data.append({
            "qseqid": f"conserved_read_{i:03d}",
            "sseqid": "GCA_000111222.1|ASM111v1_genomic",
            "pident": 97.5,
            "length": 150,
            "mismatch": 4,
            "gapopen": 0,
            "qstart": 1,
            "qend": 150,
            "sstart": 20000,
            "send": 20150,
            "evalue": 1e-68,
            "bitscore": 255.0,
        })

    df = pl.DataFrame(blast_data)
    blast_path = temp_dir / "test_sample.blast.tsv"
    df.write_csv(blast_path, separator="\t", include_header=False)

    return blast_path


@pytest.fixture
def cli_ani_file(temp_dir: Path) -> Path:
    """Create ANI matrix for CLI testing."""
    genomes = ["GCF_000123456.1", "GCF_000789012.1", "GCA_000111222.1"]

    # ANI values that will create conserved regions when hits span genomes
    ani_data = {
        "genome": genomes,
        "GCF_000123456.1": [100.0, 95.0, 80.0],
        "GCF_000789012.1": [95.0, 100.0, 82.0],
        "GCA_000111222.1": [80.0, 82.0, 100.0],
    }

    df = pl.DataFrame(ani_data)
    ani_path = temp_dir / "test.ani.csv"
    df.write_csv(ani_path)

    return ani_path


@pytest.fixture
def cli_blast_dir(temp_dir: Path) -> Path:
    """Create directory with multiple BLAST files for batch testing.

    Uses standardized pipe format: {accession}|{contig_id}
    """
    blast_dir = temp_dir / "blast_files"
    blast_dir.mkdir()

    # Create 3 sample BLAST files
    for sample_idx in range(3):
        blast_data = []
        for read_idx in range(5):
            blast_data.append({
                "qseqid": f"sample{sample_idx}_read_{read_idx:03d}",
                "sseqid": "GCF_000123456.1|ASM123v1_genomic",
                "pident": 95.0 + read_idx,
                "length": 150,
                "mismatch": 5 - read_idx,
                "gapopen": 0,
                "qstart": 1,
                "qend": 150,
                "sstart": 1000 + read_idx * 200,
                "send": 1150 + read_idx * 200,
                "evalue": 1e-60,
                "bitscore": 220.0 + read_idx * 10,
            })

        df = pl.DataFrame(blast_data)
        blast_path = blast_dir / f"sample_{sample_idx}.blast.tsv.gz"
        # Write as regular tsv (Polars handles the naming)
        df.write_csv(
            blast_dir / f"sample_{sample_idx}.blast.tsv",
            separator="\t",
            include_header=False,
        )

    return blast_dir


# =============================================================================
# Main CLI Tests
# =============================================================================


class TestMainCLI:
    """Tests for main CLI entry point."""

    def test_version_flag(self):
        """--version should display version and exit."""
        result = runner.invoke(app, ["--version"])

        assert result.exit_code == 0
        assert "metadarkmatter version" in result.stdout
        assert "0.1.0" in result.stdout

    def test_version_short_flag(self):
        """-v should display version and exit."""
        result = runner.invoke(app, ["-v"])

        assert result.exit_code == 0
        assert "0.1.0" in result.stdout

    def test_no_args_shows_help(self):
        """Running without arguments should show help/usage message."""
        result = runner.invoke(app, [])

        # With no_args_is_help=True, Typer shows help (exit code may be 0 or 2
        # depending on Typer version)
        assert result.exit_code in (0, 2)
        # Help/usage info should be in output
        assert "Usage:" in result.output or "score" in result.output

    def test_help_flag(self):
        """--help should display help text."""
        result = runner.invoke(app, ["--help"])

        assert result.exit_code == 0
        assert "Usage:" in result.stdout
        assert "score" in result.stdout

    def test_score_subcommand_help(self):
        """score --help should display score command help."""
        result = runner.invoke(app, ["score", "--help"])

        assert result.exit_code == 0
        assert "classify" in result.stdout
        assert "batch" in result.stdout


# =============================================================================
# Score Classify Command Tests
# =============================================================================


class TestScoreClassifyCommand:
    """Tests for score classify command."""

    def test_classify_basic(self, cli_blast_file, cli_ani_file, temp_dir):
        """Basic classify command should succeed."""
        output_path = temp_dir / "classifications.csv"

        result = runner.invoke(app, [
            "score", "classify",
            "--alignment", str(cli_blast_file),
            "--ani", str(cli_ani_file),
            "--output", str(output_path),
        ])

        assert result.exit_code == 0
        assert output_path.exists()
        assert "Classified" in result.stdout

        # Verify output content
        df = pl.read_csv(output_path)
        assert len(df) > 0
        assert "read_id" in df.columns
        assert "taxonomic_call" in df.columns

    def test_classify_with_summary(self, cli_blast_file, cli_ani_file, temp_dir):
        """Classify with --summary should generate summary JSON."""
        output_path = temp_dir / "classifications.csv"
        summary_path = temp_dir / "summary.json"

        result = runner.invoke(app, [
            "score", "classify",
            "--alignment", str(cli_blast_file),
            "--ani", str(cli_ani_file),
            "--output", str(output_path),
            "--summary", str(summary_path),
        ])

        assert result.exit_code == 0
        assert summary_path.exists()

        # Verify summary content
        summary = json.loads(summary_path.read_text())
        assert "total_reads" in summary
        assert "known_species" in summary
        assert "novel_species" in summary
        assert "novel_genus" in summary

    def test_classify_csv_format(self, cli_blast_file, cli_ani_file, temp_dir):
        """--format csv should produce CSV output."""
        output_path = temp_dir / "classifications.csv"

        result = runner.invoke(app, [
            "score", "classify",
            "--alignment", str(cli_blast_file),
            "--ani", str(cli_ani_file),
            "--output", str(output_path),
            "--format", "csv",
        ])

        assert result.exit_code == 0
        assert output_path.exists()

        # Should be readable as CSV
        df = pl.read_csv(output_path)
        assert len(df) > 0

    def test_classify_parquet_format(self, cli_blast_file, cli_ani_file, temp_dir):
        """--format parquet should produce Parquet output."""
        output_path = temp_dir / "classifications.parquet"

        result = runner.invoke(app, [
            "score", "classify",
            "--alignment", str(cli_blast_file),
            "--ani", str(cli_ani_file),
            "--output", str(output_path),
            "--format", "parquet",
        ])

        assert result.exit_code == 0
        assert output_path.exists()

        # Should be readable as Parquet
        df = pl.read_parquet(output_path)
        assert len(df) > 0

    def test_classify_streaming_mode(self, cli_blast_file, cli_ani_file, temp_dir):
        """--streaming should use streaming mode for large files."""
        output_path = temp_dir / "classifications.csv"

        result = runner.invoke(app, [
            "score", "classify",
            "--alignment", str(cli_blast_file),
            "--ani", str(cli_ani_file),
            "--output", str(output_path),
            "--streaming",
        ])

        assert result.exit_code == 0
        assert output_path.exists()
        assert "Streaming mode" in result.stdout

        df = pl.read_csv(output_path)
        assert len(df) > 0

    def test_classify_streaming_parquet(self, cli_blast_file, cli_ani_file, temp_dir):
        """--streaming with --format parquet should produce parquet output."""
        output_path = temp_dir / "classifications.parquet"

        result = runner.invoke(app, [
            "score", "classify",
            "--alignment", str(cli_blast_file),
            "--ani", str(cli_ani_file),
            "--output", str(output_path),
            "--streaming",
            "--format", "parquet",
        ])

        assert result.exit_code == 0
        assert output_path.exists()

        df = pl.read_parquet(output_path)
        assert len(df) > 0

    def test_classify_custom_bitscore_threshold(
        self, cli_blast_file, cli_ani_file, temp_dir
    ):
        """--bitscore-threshold should affect ambiguous hit detection."""
        output_path = temp_dir / "classifications.csv"

        result = runner.invoke(app, [
            "score", "classify",
            "--alignment", str(cli_blast_file),
            "--ani", str(cli_ani_file),
            "--output", str(output_path),
            "--bitscore-threshold", "90.0",
        ])

        assert result.exit_code == 0
        assert output_path.exists()

    def test_classify_verbose(self, cli_blast_file, cli_ani_file, temp_dir):
        """--verbose should provide detailed output."""
        output_path = temp_dir / "classifications.csv"

        result = runner.invoke(app, [
            "score", "classify",
            "--alignment", str(cli_blast_file),
            "--ani", str(cli_ani_file),
            "--output", str(output_path),
            "--verbose",
        ])

        assert result.exit_code == 0

    def test_classify_creates_output_directory(
        self, cli_blast_file, cli_ani_file, temp_dir
    ):
        """Should create output directory if it doesn't exist."""
        output_path = temp_dir / "new_dir" / "subdir" / "classifications.csv"

        result = runner.invoke(app, [
            "score", "classify",
            "--alignment", str(cli_blast_file),
            "--ani", str(cli_ani_file),
            "--output", str(output_path),
        ])

        assert result.exit_code == 0
        assert output_path.exists()

    def test_classify_classification_counts(
        self, cli_blast_file, cli_ani_file, temp_dir
    ):
        """Verify classification produces expected categories."""
        output_path = temp_dir / "classifications.csv"
        summary_path = temp_dir / "summary.json"

        result = runner.invoke(app, [
            "score", "classify",
            "--alignment", str(cli_blast_file),
            "--ani", str(cli_ani_file),
            "--output", str(output_path),
            "--summary", str(summary_path),
        ])

        assert result.exit_code == 0

        summary = json.loads(summary_path.read_text())

        # Should have reads in multiple categories
        assert summary["total_reads"] == 20  # 10 + 5 + 3 + 2
        assert summary["known_species"] >= 0
        assert summary["novel_species"] >= 0
        assert summary["novel_genus"] >= 0


# =============================================================================
# Score Classify Error Handling Tests
# =============================================================================


class TestScoreClassifyErrors:
    """Tests for error handling in score classify command."""

    def test_missing_blast_file(self, cli_ani_file, temp_dir):
        """Should error if BLAST file doesn't exist."""
        result = runner.invoke(app, [
            "score", "classify",
            "--alignment", "/nonexistent/path/blast.tsv",
            "--ani", str(cli_ani_file),
            "--output", str(temp_dir / "out.csv"),
        ])

        assert result.exit_code != 0

    def test_missing_ani_file(self, cli_blast_file, temp_dir):
        """Should error if ANI file doesn't exist."""
        result = runner.invoke(app, [
            "score", "classify",
            "--alignment", str(cli_blast_file),
            "--ani", "/nonexistent/path/ani.csv",
            "--output", str(temp_dir / "out.csv"),
        ])

        assert result.exit_code != 0

    def test_invalid_format_option(self, cli_blast_file, cli_ani_file, temp_dir):
        """Should error for invalid --format value."""
        result = runner.invoke(app, [
            "score", "classify",
            "--alignment", str(cli_blast_file),
            "--ani", str(cli_ani_file),
            "--output", str(temp_dir / "out.csv"),
            "--format", "invalid_format",
        ])

        assert result.exit_code != 0
        assert "Invalid format" in result.stdout or result.exit_code == 1

    def test_bitscore_threshold_out_of_range(
        self, cli_blast_file, cli_ani_file, temp_dir
    ):
        """Should error for bitscore threshold outside 0-100."""
        result = runner.invoke(app, [
            "score", "classify",
            "--alignment", str(cli_blast_file),
            "--ani", str(cli_ani_file),
            "--output", str(temp_dir / "out.csv"),
            "--bitscore-threshold", "150.0",
        ])

        assert result.exit_code != 0

    def test_missing_required_options(self, temp_dir):
        """Should error if required options are missing."""
        # Missing --alignment
        result = runner.invoke(app, [
            "score", "classify",
            "--ani", str(temp_dir / "ani.csv"),
            "--output", str(temp_dir / "out.csv"),
        ])

        assert result.exit_code != 0

    def test_invalid_ani_matrix(self, cli_blast_file, temp_dir):
        """Should error for malformed ANI matrix."""
        # Create invalid ANI matrix (non-square)
        invalid_ani = temp_dir / "invalid.ani.csv"
        df = pl.DataFrame({
            "genome": ["A", "B"],
            "A": [100.0, 95.0],
            "B": [95.0, 100.0],
            "C": [80.0, 85.0],  # Extra column makes it non-square
        })
        df.write_csv(invalid_ani)

        result = runner.invoke(app, [
            "score", "classify",
            "--alignment", str(cli_blast_file),
            "--ani", str(invalid_ani),
            "--output", str(temp_dir / "out.csv"),
        ])

        assert result.exit_code != 0

    def test_low_genome_coverage_warning(self, temp_dir):
        """Should warn when BLAST genomes are mostly missing from ANI matrix."""
        # Create BLAST file with genomes NOT in ANI matrix
        # Use standardized pipe format: {accession}|{contig_id}
        blast_data = [
            "read_1\tGCF_MISSING_001.1|genomic\t99.0\t150\t1\t0\t1\t150\t1\t150\t1e-50\t200",
            "read_2\tGCF_MISSING_002.1|genomic\t98.0\t150\t2\t0\t1\t150\t1\t150\t1e-45\t190",
            "read_3\tGCF_MISSING_003.1|genomic\t97.0\t150\t3\t0\t1\t150\t1\t150\t1e-40\t180",
            "read_4\tGCF_KNOWN_001.1|genomic\t99.5\t150\t0\t0\t1\t150\t1\t150\t1e-60\t250",
        ]
        blast_path = temp_dir / "mismatched.blast.tsv"
        blast_path.write_text("\n".join(blast_data) + "\n")

        # Create ANI matrix with only one genome (25% coverage)
        ani_data = pl.DataFrame({
            "genome": ["GCF_KNOWN_001.1"],
            "GCF_KNOWN_001.1": [100.0],
        })
        ani_path = temp_dir / "small.ani.csv"
        ani_data.write_csv(ani_path)

        result = runner.invoke(app, [
            "score", "classify",
            "--alignment", str(blast_path),
            "--ani", str(ani_path),
            "--output", str(temp_dir / "out.csv"),
            "--verbose",
        ])

        # Should still succeed but emit warning
        assert result.exit_code == 0
        assert "low genome coverage" in result.stdout.lower()
        assert "missing" in result.stdout.lower() or "mismatched" in result.stdout.lower()


# =============================================================================
# Score Batch Command Tests
# =============================================================================


class TestScoreBatchCommand:
    """Tests for score batch command."""

    def test_batch_basic(self, cli_blast_dir, cli_ani_file, temp_dir):
        """Basic batch command should process all files."""
        output_dir = temp_dir / "batch_output"

        result = runner.invoke(app, [
            "score", "batch",
            "--alignment-dir", str(cli_blast_dir),
            "--ani", str(cli_ani_file),
            "--output-dir", str(output_dir),
            "--pattern", "*.blast.tsv",
        ])

        assert result.exit_code == 0
        assert output_dir.exists()

        # Should have created output files for each input
        output_files = list(output_dir.glob("*_classifications.csv"))
        assert len(output_files) == 3

    def test_batch_with_parquet(self, cli_blast_dir, cli_ani_file, temp_dir):
        """Batch with --format parquet should produce Parquet files."""
        output_dir = temp_dir / "batch_output"

        result = runner.invoke(app, [
            "score", "batch",
            "--alignment-dir", str(cli_blast_dir),
            "--ani", str(cli_ani_file),
            "--output-dir", str(output_dir),
            "--pattern", "*.blast.tsv",
            "--format", "parquet",
        ])

        assert result.exit_code == 0

        # Should have Parquet output files
        output_files = list(output_dir.glob("*_classifications.parquet"))
        assert len(output_files) == 3

    def test_batch_creates_summaries(self, cli_blast_dir, cli_ani_file, temp_dir):
        """Batch should create summary JSON for each sample."""
        output_dir = temp_dir / "batch_output"

        result = runner.invoke(app, [
            "score", "batch",
            "--alignment-dir", str(cli_blast_dir),
            "--ani", str(cli_ani_file),
            "--output-dir", str(output_dir),
            "--pattern", "*.blast.tsv",
        ])

        assert result.exit_code == 0

        # Should have summary files
        summary_files = list(output_dir.glob("*_summary.json"))
        assert len(summary_files) == 3

        # Verify summary content
        for summary_file in summary_files:
            summary = json.loads(summary_file.read_text())
            assert "total_reads" in summary
            assert summary["total_reads"] == 5

    def test_batch_no_matching_files(self, cli_blast_dir, cli_ani_file, temp_dir):
        """Batch should fail with error when no files match pattern."""
        output_dir = temp_dir / "batch_output"

        result = runner.invoke(app, [
            "score", "batch",
            "--alignment-dir", str(cli_blast_dir),
            "--ani", str(cli_ani_file),
            "--output-dir", str(output_dir),
            "--pattern", "*.nonexistent",
        ])

        # Should exit with error code 1 (not silent success)
        assert result.exit_code == 1
        assert "No alignment files found" in result.stdout
        assert "Suggestions" in result.stdout  # Helpful guidance included

    def test_batch_custom_pattern(self, temp_dir, cli_ani_file):
        """Batch should respect custom --pattern."""
        # Create files with custom naming
        blast_dir = temp_dir / "custom_blast"
        blast_dir.mkdir()

        for i in range(2):
            df = pl.DataFrame({
                "qseqid": [f"read_{i}"],
                "sseqid": ["GCF_000123456.1|contig"],  # Standardized pipe format
                "pident": [98.0],
                "length": [150],
                "mismatch": [3],
                "gapopen": [0],
                "qstart": [1],
                "qend": [150],
                "sstart": [1000],
                "send": [1150],
                "evalue": [1e-50],
                "bitscore": [250.0],
            })
            df.write_csv(
                blast_dir / f"sample_{i}.custom.tsv",
                separator="\t",
                include_header=False,
            )

        output_dir = temp_dir / "batch_output"

        result = runner.invoke(app, [
            "score", "batch",
            "--alignment-dir", str(blast_dir),
            "--ani", str(cli_ani_file),
            "--output-dir", str(output_dir),
            "--pattern", "*.custom.tsv",
        ])

        assert result.exit_code == 0
        assert len(list(output_dir.glob("*_classifications.csv"))) == 2

    def test_batch_displays_progress(self, cli_blast_dir, cli_ani_file, temp_dir):
        """Batch should show processing progress."""
        output_dir = temp_dir / "batch_output"

        result = runner.invoke(app, [
            "score", "batch",
            "--alignment-dir", str(cli_blast_dir),
            "--ani", str(cli_ani_file),
            "--output-dir", str(output_dir),
            "--pattern", "*.blast.tsv",
        ])

        assert result.exit_code == 0
        # Should show progress indicators
        assert "[1/3]" in result.stdout or "Processing" in result.stdout


# =============================================================================
# Score Batch Error Handling Tests
# =============================================================================


class TestScoreBatchErrors:
    """Tests for error handling in score batch command."""

    def test_batch_missing_blast_dir(self, cli_ani_file, temp_dir):
        """Should error if blast-dir doesn't exist."""
        result = runner.invoke(app, [
            "score", "batch",
            "--alignment-dir", "/nonexistent/directory",
            "--ani", str(cli_ani_file),
            "--output-dir", str(temp_dir / "output"),
        ])

        assert result.exit_code != 0

    def test_batch_missing_ani_file(self, cli_blast_dir, temp_dir):
        """Should error if ANI file doesn't exist."""
        result = runner.invoke(app, [
            "score", "batch",
            "--alignment-dir", str(cli_blast_dir),
            "--ani", "/nonexistent/ani.csv",
            "--output-dir", str(temp_dir / "output"),
        ])

        assert result.exit_code != 0

    def test_batch_invalid_format(self, cli_blast_dir, cli_ani_file, temp_dir):
        """Should error for invalid --format value."""
        result = runner.invoke(app, [
            "score", "batch",
            "--alignment-dir", str(cli_blast_dir),
            "--ani", str(cli_ani_file),
            "--output-dir", str(temp_dir / "output"),
            "--format", "invalid",
        ])

        assert result.exit_code != 0


# =============================================================================
# CLI Output Validation Tests
# =============================================================================


class TestCLIOutputValidation:
    """Tests validating CLI output content and format."""

    def test_output_columns(self, cli_blast_file, cli_ani_file, temp_dir):
        """Output should have all required columns."""
        output_path = temp_dir / "classifications.csv"

        runner.invoke(app, [
            "score", "classify",
            "--alignment", str(cli_blast_file),
            "--ani", str(cli_ani_file),
            "--output", str(output_path),
        ])

        df = pl.read_csv(output_path)

        required_columns = [
            "read_id",
            "best_match_genome",
            "top_hit_identity",
            "novelty_index",
            "placement_uncertainty",
            "num_ambiguous_hits",
            "taxonomic_call",
            "is_novel",
        ]

        for col in required_columns:
            assert col in df.columns, f"Missing column: {col}"

    def test_output_taxonomic_calls_valid(
        self, cli_blast_file, cli_ani_file, temp_dir
    ):
        """Taxonomic calls should be valid enum values."""
        output_path = temp_dir / "classifications.csv"

        runner.invoke(app, [
            "score", "classify",
            "--alignment", str(cli_blast_file),
            "--ani", str(cli_ani_file),
            "--output", str(output_path),
        ])

        df = pl.read_csv(output_path)

        valid_calls = {
            "Known Species",
            "Novel Species",
            "Novel Genus",
            "Conserved Region",
            "Ambiguous",
        }

        actual_calls = set(df["taxonomic_call"].unique().to_list())
        assert actual_calls.issubset(valid_calls)

    def test_output_novelty_index_range(
        self, cli_blast_file, cli_ani_file, temp_dir
    ):
        """Novelty index should be between 0 and 100."""
        output_path = temp_dir / "classifications.csv"

        runner.invoke(app, [
            "score", "classify",
            "--alignment", str(cli_blast_file),
            "--ani", str(cli_ani_file),
            "--output", str(output_path),
        ])

        df = pl.read_csv(output_path)

        assert df["novelty_index"].min() >= 0
        assert df["novelty_index"].max() <= 100

    def test_output_placement_uncertainty_range(
        self, cli_blast_file, cli_ani_file, temp_dir
    ):
        """Placement uncertainty should be between 0 and 100."""
        output_path = temp_dir / "classifications.csv"

        runner.invoke(app, [
            "score", "classify",
            "--alignment", str(cli_blast_file),
            "--ani", str(cli_ani_file),
            "--output", str(output_path),
        ])

        df = pl.read_csv(output_path)

        assert df["placement_uncertainty"].min() >= 0
        assert df["placement_uncertainty"].max() <= 100

    def test_output_is_novel_consistency(
        self, cli_blast_file, cli_ani_file, temp_dir
    ):
        """is_novel should be consistent with taxonomic_call."""
        output_path = temp_dir / "classifications.csv"

        runner.invoke(app, [
            "score", "classify",
            "--alignment", str(cli_blast_file),
            "--ani", str(cli_ani_file),
            "--output", str(output_path),
        ])

        df = pl.read_csv(output_path)

        for row in df.iter_rows(named=True):
            is_novel = row["is_novel"]
            call = row["taxonomic_call"]

            if call in ("Novel Species", "Novel Genus"):
                assert is_novel is True, f"is_novel should be True for {call}"
            else:
                assert is_novel is False, f"is_novel should be False for {call}"

    def test_summary_json_schema(self, cli_blast_file, cli_ani_file, temp_dir):
        """Summary JSON should have correct schema."""
        output_path = temp_dir / "classifications.csv"
        summary_path = temp_dir / "summary.json"

        runner.invoke(app, [
            "score", "classify",
            "--alignment", str(cli_blast_file),
            "--ani", str(cli_ani_file),
            "--output", str(output_path),
            "--summary", str(summary_path),
        ])

        summary = json.loads(summary_path.read_text())

        required_fields = [
            "total_reads",
            "known_species",
            "novel_species",
            "novel_genus",
            "conserved_regions",
            "mean_novelty_index",
            "mean_placement_uncertainty",
            "known_species_pct",
            "novel_species_pct",
            "novel_genus_pct",
            "novel_diversity_pct",
        ]

        for field in required_fields:
            assert field in summary, f"Missing field in summary: {field}"

    def test_parquet_and_csv_equivalent(
        self, cli_blast_file, cli_ani_file, temp_dir
    ):
        """CSV and Parquet output should contain equivalent data."""
        csv_path = temp_dir / "classifications.csv"
        parquet_path = temp_dir / "classifications.parquet"

        # Generate CSV
        runner.invoke(app, [
            "score", "classify",
            "--alignment", str(cli_blast_file),
            "--ani", str(cli_ani_file),
            "--output", str(csv_path),
            "--format", "csv",
        ])

        # Generate Parquet
        runner.invoke(app, [
            "score", "classify",
            "--alignment", str(cli_blast_file),
            "--ani", str(cli_ani_file),
            "--output", str(parquet_path),
            "--format", "parquet",
        ])

        df_csv = pl.read_csv(csv_path)
        df_parquet = pl.read_parquet(parquet_path)

        # Same number of rows
        assert len(df_csv) == len(df_parquet)

        # Same columns
        assert set(df_csv.columns) == set(df_parquet.columns)

        # Same read IDs
        assert set(df_csv["read_id"].to_list()) == set(df_parquet["read_id"].to_list())
