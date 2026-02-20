"""
Unit tests for visualize and mmseqs2 CLI commands.

Tests argument validation, file handling, tool availability checks,
and output handling for the visualize and mmseqs2 subcommands.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from unittest.mock import MagicMock, patch

import polars as pl
import pytest

from metadarkmatter.cli.main import app
from metadarkmatter.external.base import ToolExecutionError, ToolResult


# =============================================================================
# Fixtures specific to these tests
# =============================================================================


@pytest.fixture
def cli_runner():
    """Typer CLI runner for testing commands."""
    from typer.testing import CliRunner

    return CliRunner()


@pytest.fixture
def temp_dir(tmp_path: Path) -> Path:
    """Temporary directory for test artifacts."""
    return tmp_path


@pytest.fixture
def temp_bam_file(temp_dir: Path) -> Path:
    """Create a dummy BAM file (empty, for path validation only)."""
    bam_path = temp_dir / "sample.bam"
    bam_path.write_bytes(b"\x00" * 64)
    return bam_path


@pytest.fixture
def temp_classification_file(temp_dir: Path) -> Path:
    """Temporary classification CSV file with test data."""
    data = {
        "read_id": [f"read_{i:03d}" for i in range(10)],
        "best_match_genome": ["GCF_000123456.1"] * 5 + ["GCF_000789012.1"] * 5,
        "taxonomic_call": (
            ["Known Species"] * 3
            + ["Novel Species"] * 2
            + ["Novel Genus"] * 2
            + ["Conserved Region"] * 3
        ),
        "novelty_index": [1.0, 2.0, 3.0, 8.0, 12.0, 18.0, 22.0, 5.0, 6.0, 7.0],
        "placement_uncertainty": [
            0.2, 0.3, 0.4, 0.5, 0.6, 1.0, 1.5, 8.0, 9.0, 10.0,
        ],
        "top_hit_identity": [
            99.0, 98.0, 97.0, 92.0, 88.0, 82.0, 78.0, 95.0, 94.0, 93.0,
        ],
    }
    df = pl.DataFrame(data)
    output_path = temp_dir / "classifications.csv"
    df.write_csv(output_path)
    return output_path


@pytest.fixture
def temp_genome_dir(temp_dir: Path) -> Path:
    """Create a temporary directory with a genome FASTA file."""
    genome_dir = temp_dir / "genomes"
    genome_dir.mkdir()
    fasta_path = genome_dir / "GCF_000123456.1_ASM123v1_genomic.fna"
    fasta_path.write_text(">contig_1\nACGTACGTACGT\n>contig_2\nGGGGCCCCAAAA\n")
    return genome_dir


@pytest.fixture
def temp_single_fasta(temp_dir: Path) -> Path:
    """Create a single temporary FASTA file."""
    fasta_path = temp_dir / "pangenome.fasta"
    fasta_path.write_text(">seq1\nACGTACGT\n>seq2\nGGGGCCCC\n")
    return fasta_path


@pytest.fixture
def temp_query_fasta(temp_dir: Path) -> Path:
    """Create a temporary query FASTA file."""
    query_path = temp_dir / "reads.fasta"
    query_path.write_text(">read_1\nACGTACGT\n>read_2\nGGGGCCCC\n")
    return query_path


@pytest.fixture
def temp_query_fastq_r1(temp_dir: Path) -> Path:
    """Create a temporary forward FASTQ file."""
    fq_path = temp_dir / "reads_R1.fastq"
    fq_path.write_text(
        "@read_1\nACGTACGT\n+\nIIIIIIII\n"
        "@read_2\nGGGGCCCC\n+\nIIIIIIII\n"
    )
    return fq_path


@pytest.fixture
def temp_query_fastq_r2(temp_dir: Path) -> Path:
    """Create a temporary reverse FASTQ file."""
    fq_path = temp_dir / "reads_R2.fastq"
    fq_path.write_text(
        "@read_1\nTGCATGCA\n+\nIIIIIIII\n"
        "@read_2\nCCCCGGGG\n+\nIIIIIIII\n"
    )
    return fq_path


@pytest.fixture
def temp_query_fastq_gz_r1(temp_dir: Path) -> Path:
    """Create a temporary gzipped forward FASTQ file."""
    fq_path = temp_dir / "reads_R1.fastq.gz"
    with gzip.open(fq_path, "wt") as f:
        f.write("@read_1\nACGTACGT\n+\nIIIIIIII\n")
        f.write("@read_2\nGGGGCCCC\n+\nIIIIIIII\n")
    return fq_path


@pytest.fixture
def temp_query_fastq_gz_r2(temp_dir: Path) -> Path:
    """Create a temporary gzipped reverse FASTQ file."""
    fq_path = temp_dir / "reads_R2.fastq.gz"
    with gzip.open(fq_path, "wt") as f:
        f.write("@read_1\nTGCATGCA\n+\nIIIIIIII\n")
        f.write("@read_2\nCCCCGGGG\n+\nIIIIIIII\n")
    return fq_path


@pytest.fixture
def temp_mmseqs_db(temp_dir: Path) -> Path:
    """Create dummy MMseqs2 database files."""
    db_dir = temp_dir / "mmseqs_db"
    db_dir.mkdir()
    db_path = db_dir / "pangenome"
    # MMseqs2 databases consist of several files
    db_path.write_bytes(b"\x00" * 32)
    db_path.with_suffix(".dbtype").write_bytes(b"\x00")
    db_path.with_suffix(".index").write_bytes(b"\x00")
    return db_path


def _make_tool_result(
    command: tuple[str, ...] = ("mmseqs",),
    return_code: int = 0,
    elapsed: float = 1.0,
) -> ToolResult:
    """Create a mock ToolResult for testing."""
    return ToolResult(
        command=command,
        return_code=return_code,
        stdout="",
        stderr="",
        elapsed_seconds=elapsed,
    )


# =============================================================================
# Visualize: recruitment command tests
# =============================================================================


class TestVisualizeRecruitmentValidation:
    """Tests for visualize recruitment command argument validation."""

    def test_recruitment_requires_bam(self, cli_runner, temp_dir):
        """Test that --bam argument is required."""
        output = temp_dir / "plot.html"
        result = cli_runner.invoke(
            app,
            ["visualize", "recruitment", "--output", str(output)],
        )
        assert result.exit_code != 0
        assert "--bam" in result.output or "bam" in result.output.lower()

    def test_recruitment_requires_output(self, cli_runner, temp_bam_file):
        """Test that --output argument is required."""
        result = cli_runner.invoke(
            app,
            ["visualize", "recruitment", "--bam", str(temp_bam_file)],
        )
        assert result.exit_code != 0
        assert "--output" in result.output or "output" in result.output.lower()

    def test_recruitment_missing_bam_file(self, cli_runner, temp_dir):
        """Test error handling for non-existent BAM file."""
        output = temp_dir / "plot.html"
        result = cli_runner.invoke(
            app,
            [
                "visualize", "recruitment",
                "--bam", str(temp_dir / "nonexistent.bam"),
                "--output", str(output),
            ],
        )
        assert result.exit_code != 0
        assert (
            "not found" in result.output.lower()
            or "does not exist" in result.output.lower()
            or "invalid" in result.output.lower()
        )

    def test_recruitment_invalid_format(self, cli_runner, temp_bam_file, temp_dir):
        """Test error on invalid output format."""
        output = temp_dir / "plot.xyz"
        with patch(
            "metadarkmatter.cli.visualize.Samtools.check_available",
            return_value=True,
        ):
            result = cli_runner.invoke(
                app,
                [
                    "visualize", "recruitment",
                    "--bam", str(temp_bam_file),
                    "--output", str(output),
                    "--format", "xyz",
                ],
            )
        assert result.exit_code != 0
        assert "invalid" in result.output.lower() or "format" in result.output.lower()

    def test_recruitment_samtools_not_available(
        self, cli_runner, temp_bam_file, temp_dir
    ):
        """Test error when samtools is not installed."""
        output = temp_dir / "plot.html"
        with patch(
            "metadarkmatter.cli.visualize.Samtools.check_available",
            return_value=False,
        ):
            result = cli_runner.invoke(
                app,
                [
                    "visualize", "recruitment",
                    "--bam", str(temp_bam_file),
                    "--output", str(output),
                ],
            )
        assert result.exit_code != 0
        assert "samtools" in result.output.lower()

    @patch("metadarkmatter.cli.visualize.Samtools.check_available", return_value=True)
    @patch("metadarkmatter.cli.visualize.load_recruitment_data")
    def test_recruitment_bam_read_error(
        self, mock_load, mock_samtools, cli_runner, temp_bam_file, temp_dir
    ):
        """Test error handling when BAM file cannot be read."""
        mock_load.side_effect = RuntimeError("Invalid BAM format")
        output = temp_dir / "plot.html"
        result = cli_runner.invoke(
            app,
            [
                "visualize", "recruitment",
                "--bam", str(temp_bam_file),
                "--output", str(output),
            ],
        )
        assert result.exit_code != 0
        assert "error" in result.output.lower()

    @patch("metadarkmatter.cli.visualize.Samtools.check_available", return_value=True)
    @patch("metadarkmatter.cli.visualize.load_recruitment_data")
    def test_recruitment_empty_data(
        self, mock_load, mock_samtools, cli_runner, temp_bam_file, temp_dir
    ):
        """Test handling when no alignments pass filters."""
        mock_load.return_value = pl.DataFrame(
            schema={
                "read_id": pl.Utf8,
                "genome_name": pl.Utf8,
                "position": pl.Int64,
                "percent_identity": pl.Float64,
            }
        )
        output = temp_dir / "plot.html"
        result = cli_runner.invoke(
            app,
            [
                "visualize", "recruitment",
                "--bam", str(temp_bam_file),
                "--output", str(output),
            ],
        )
        # Exit code 0 because the command prints a warning and exits gracefully
        assert result.exit_code == 0
        assert "no alignments" in result.output.lower() or "warning" in result.output.lower()

    @patch("metadarkmatter.cli.visualize.Samtools.check_available", return_value=True)
    @patch("metadarkmatter.cli.visualize.load_recruitment_data")
    def test_recruitment_happy_path_html(
        self, mock_load, mock_samtools, cli_runner, temp_bam_file, temp_dir
    ):
        """Test successful recruitment plot generation (HTML format)."""
        # Create mock recruitment data
        mock_data = pl.DataFrame({
            "read_id": [f"read_{i}" for i in range(20)],
            "genome_name": ["GCF_000123456.1"] * 20,
            "position": list(range(100, 2100, 100)),
            "percent_identity": [95.0 + (i * 0.2) for i in range(20)],
        })
        mock_load.return_value = mock_data

        output = temp_dir / "plot.html"

        # Mock the visualization classes at their source module
        mock_fig = MagicMock()
        mock_generator_instance = MagicMock()
        mock_generator_instance.create_figure.return_value = mock_fig
        mock_generator_cls = MagicMock(return_value=mock_generator_instance)
        mock_band_cls = MagicMock()

        with patch(
            "metadarkmatter.visualization.RecruitmentPlotGenerator",
            mock_generator_cls,
        ), patch(
            "metadarkmatter.visualization.IdentityBand",
            mock_band_cls,
        ):
            result = cli_runner.invoke(
                app,
                [
                    "visualize", "recruitment",
                    "--bam", str(temp_bam_file),
                    "--output", str(output),
                    "--format", "html",
                    "--quiet",
                ],
            )

        assert result.exit_code == 0
        mock_fig.write_html.assert_called_once()

    @patch("metadarkmatter.cli.visualize.Samtools.check_available", return_value=True)
    @patch("metadarkmatter.cli.visualize.load_recruitment_data")
    def test_recruitment_happy_path_json(
        self, mock_load, mock_samtools, cli_runner, temp_bam_file, temp_dir
    ):
        """Test successful recruitment plot generation (JSON format)."""
        mock_data = pl.DataFrame({
            "read_id": [f"read_{i}" for i in range(5)],
            "genome_name": ["GCF_000123456.1"] * 5,
            "position": list(range(100, 600, 100)),
            "percent_identity": [95.0, 96.0, 97.0, 98.0, 99.0],
        })
        mock_load.return_value = mock_data

        output = temp_dir / "plot.json"
        mock_fig = MagicMock()
        mock_generator_instance = MagicMock()
        mock_generator_instance.create_figure.return_value = mock_fig
        mock_generator_cls = MagicMock(return_value=mock_generator_instance)

        with patch(
            "metadarkmatter.visualization.RecruitmentPlotGenerator",
            mock_generator_cls,
        ), patch("metadarkmatter.visualization.IdentityBand"):
            result = cli_runner.invoke(
                app,
                [
                    "visualize", "recruitment",
                    "--bam", str(temp_bam_file),
                    "--output", str(output),
                    "--format", "json",
                    "--quiet",
                ],
            )

        assert result.exit_code == 0
        mock_fig.write_json.assert_called_once()

    @patch("metadarkmatter.cli.visualize.Samtools.check_available", return_value=True)
    @patch("metadarkmatter.cli.visualize.load_recruitment_data")
    def test_recruitment_happy_path_png(
        self, mock_load, mock_samtools, cli_runner, temp_bam_file, temp_dir
    ):
        """Test successful recruitment plot generation (PNG format)."""
        mock_data = pl.DataFrame({
            "read_id": [f"read_{i}" for i in range(5)],
            "genome_name": ["GCF_000123456.1"] * 5,
            "position": list(range(100, 600, 100)),
            "percent_identity": [95.0, 96.0, 97.0, 98.0, 99.0],
        })
        mock_load.return_value = mock_data

        output = temp_dir / "plot.png"
        mock_fig = MagicMock()
        mock_generator_instance = MagicMock()
        mock_generator_instance.create_figure.return_value = mock_fig
        mock_generator_cls = MagicMock(return_value=mock_generator_instance)

        with patch(
            "metadarkmatter.visualization.RecruitmentPlotGenerator",
            mock_generator_cls,
        ), patch("metadarkmatter.visualization.IdentityBand"):
            result = cli_runner.invoke(
                app,
                [
                    "visualize", "recruitment",
                    "--bam", str(temp_bam_file),
                    "--output", str(output),
                    "--format", "png",
                    "--quiet",
                ],
            )

        assert result.exit_code == 0
        mock_fig.write_image.assert_called_once()

    @patch("metadarkmatter.cli.visualize.Samtools.check_available", return_value=True)
    @patch("metadarkmatter.cli.visualize.load_recruitment_data")
    def test_recruitment_multi_panel(
        self, mock_load, mock_samtools, cli_runner, temp_bam_file, temp_dir
    ):
        """Test multi-panel mode creates multi-genome figure."""
        mock_data = pl.DataFrame({
            "read_id": [f"read_{i}" for i in range(5)],
            "genome_name": ["GCF_000123456.1"] * 5,
            "position": list(range(100, 600, 100)),
            "percent_identity": [95.0, 96.0, 97.0, 98.0, 99.0],
        })
        mock_load.return_value = mock_data

        output = temp_dir / "plot.html"
        mock_fig = MagicMock()
        mock_generator_instance = MagicMock()
        mock_generator_instance.create_multi_genome_figure.return_value = mock_fig
        mock_generator_cls = MagicMock(return_value=mock_generator_instance)

        with patch(
            "metadarkmatter.visualization.RecruitmentPlotGenerator",
            mock_generator_cls,
        ), patch("metadarkmatter.visualization.IdentityBand"):
            result = cli_runner.invoke(
                app,
                [
                    "visualize", "recruitment",
                    "--bam", str(temp_bam_file),
                    "--output", str(output),
                    "--multi-panel",
                    "--quiet",
                ],
            )

        assert result.exit_code == 0
        mock_generator_instance.create_multi_genome_figure.assert_called_once()

    @patch("metadarkmatter.cli.visualize.Samtools.check_available", return_value=True)
    @patch("metadarkmatter.cli.visualize.load_recruitment_data")
    def test_recruitment_with_genome_filter(
        self, mock_load, mock_samtools, cli_runner, temp_bam_file, temp_dir
    ):
        """Test filtering to specific genomes."""
        mock_data = pl.DataFrame({
            "read_id": ["r1", "r2", "r3"],
            "genome_name": ["GCF_000123456.1", "GCF_000789012.1", "GCF_000123456.1"],
            "position": [100, 200, 300],
            "percent_identity": [95.0, 96.0, 97.0],
        })
        mock_load.return_value = mock_data

        output = temp_dir / "plot.html"
        mock_fig = MagicMock()
        mock_generator_instance = MagicMock()
        mock_generator_instance.create_figure.return_value = mock_fig
        mock_generator_cls = MagicMock(return_value=mock_generator_instance)

        with patch(
            "metadarkmatter.visualization.RecruitmentPlotGenerator",
            mock_generator_cls,
        ), patch("metadarkmatter.visualization.IdentityBand"):
            result = cli_runner.invoke(
                app,
                [
                    "visualize", "recruitment",
                    "--bam", str(temp_bam_file),
                    "--output", str(output),
                    "--genome", "GCF_000123456.1",
                    "--quiet",
                ],
            )

        assert result.exit_code == 0

    @patch("metadarkmatter.cli.visualize.Samtools.check_available", return_value=True)
    @patch("metadarkmatter.cli.visualize.load_recruitment_data")
    def test_recruitment_with_export_tsv(
        self, mock_load, mock_samtools, cli_runner, temp_bam_file, temp_dir
    ):
        """Test TSV export for anvi'o compatibility."""
        mock_data = pl.DataFrame({
            "read_id": [f"read_{i}" for i in range(5)],
            "genome_name": ["GCF_000123456.1"] * 5,
            "position": list(range(100, 600, 100)),
            "percent_identity": [95.0, 96.0, 97.0, 98.0, 99.0],
        })
        mock_load.return_value = mock_data

        output = temp_dir / "plot.html"
        export_tsv = temp_dir / "export_data.tsv"
        mock_fig = MagicMock()
        mock_generator_instance = MagicMock()
        mock_generator_instance.create_figure.return_value = mock_fig
        mock_generator_cls = MagicMock(return_value=mock_generator_instance)

        with patch(
            "metadarkmatter.visualization.RecruitmentPlotGenerator",
            mock_generator_cls,
        ), patch("metadarkmatter.visualization.IdentityBand"):
            result = cli_runner.invoke(
                app,
                [
                    "visualize", "recruitment",
                    "--bam", str(temp_bam_file),
                    "--output", str(output),
                    "--export-tsv", str(export_tsv),
                    "--quiet",
                ],
            )

        assert result.exit_code == 0
        mock_generator_instance.export_for_anvio.assert_called_once()

    @patch("metadarkmatter.cli.visualize.Samtools.check_available", return_value=True)
    @patch("metadarkmatter.cli.visualize.load_recruitment_data")
    def test_recruitment_plot_creation_error(
        self, mock_load, mock_samtools, cli_runner, temp_bam_file, temp_dir
    ):
        """Test error handling when plot creation fails."""
        mock_data = pl.DataFrame({
            "read_id": ["r1"],
            "genome_name": ["GCF_000123456.1"],
            "position": [100],
            "percent_identity": [95.0],
        })
        mock_load.return_value = mock_data

        output = temp_dir / "plot.html"
        mock_generator_instance = MagicMock()
        mock_generator_instance.create_figure.side_effect = ValueError("Plot error")
        mock_generator_cls = MagicMock(return_value=mock_generator_instance)

        with patch(
            "metadarkmatter.visualization.RecruitmentPlotGenerator",
            mock_generator_cls,
        ), patch("metadarkmatter.visualization.IdentityBand"):
            result = cli_runner.invoke(
                app,
                [
                    "visualize", "recruitment",
                    "--bam", str(temp_bam_file),
                    "--output", str(output),
                    "--quiet",
                ],
            )

        assert result.exit_code != 0
        assert "error" in result.output.lower()

    @patch("metadarkmatter.cli.visualize.Samtools.check_available", return_value=True)
    @patch("metadarkmatter.cli.visualize.load_recruitment_data")
    @patch("metadarkmatter.cli.visualize.aggregate_by_genome")
    def test_recruitment_verbose_shows_genome_stats(
        self, mock_agg, mock_load, mock_samtools, cli_runner, temp_bam_file, temp_dir
    ):
        """Test verbose mode displays genome statistics."""
        mock_data = pl.DataFrame({
            "read_id": [f"read_{i}" for i in range(5)],
            "genome_name": ["GCF_000123456.1"] * 5,
            "position": list(range(100, 600, 100)),
            "percent_identity": [95.0, 96.0, 97.0, 98.0, 99.0],
        })
        mock_load.return_value = mock_data

        mock_stats = pl.DataFrame({
            "genome_name": ["GCF_000123456.1"],
            "num_reads": [5],
            "mean_identity": [97.0],
        })
        mock_agg.return_value = mock_stats

        output = temp_dir / "plot.html"
        mock_fig = MagicMock()
        mock_generator_instance = MagicMock()
        mock_generator_instance.create_figure.return_value = mock_fig
        mock_generator_cls = MagicMock(return_value=mock_generator_instance)

        with patch(
            "metadarkmatter.visualization.RecruitmentPlotGenerator",
            mock_generator_cls,
        ), patch("metadarkmatter.visualization.IdentityBand"):
            result = cli_runner.invoke(
                app,
                [
                    "visualize", "recruitment",
                    "--bam", str(temp_bam_file),
                    "--output", str(output),
                    "--verbose",
                ],
            )

        assert result.exit_code == 0
        mock_agg.assert_called_once()


# =============================================================================
# Visualize: summary command tests
# =============================================================================


class TestVisualizeSummaryValidation:
    """Tests for visualize summary command argument validation."""

    def test_summary_requires_classifications(self, cli_runner, temp_dir):
        """Test that --classifications argument is required."""
        output = temp_dir / "summary.html"
        result = cli_runner.invoke(
            app,
            ["visualize", "summary", "--output", str(output)],
        )
        assert result.exit_code != 0
        assert (
            "--classifications" in result.output
            or "classifications" in result.output.lower()
        )

    def test_summary_requires_output(self, cli_runner, temp_classification_file):
        """Test that --output argument is required."""
        result = cli_runner.invoke(
            app,
            [
                "visualize", "summary",
                "--classifications", str(temp_classification_file),
            ],
        )
        assert result.exit_code != 0
        assert "--output" in result.output or "output" in result.output.lower()

    def test_summary_missing_input_file(self, cli_runner, temp_dir):
        """Test error handling for non-existent classification file."""
        output = temp_dir / "summary.html"
        result = cli_runner.invoke(
            app,
            [
                "visualize", "summary",
                "--classifications", str(temp_dir / "nonexistent.csv"),
                "--output", str(output),
            ],
        )
        assert result.exit_code != 0

    def test_summary_happy_path_html(
        self, cli_runner, temp_classification_file, temp_dir
    ):
        """Test successful summary plot generation in HTML format."""
        output = temp_dir / "summary.html"
        result = cli_runner.invoke(
            app,
            [
                "visualize", "summary",
                "--classifications", str(temp_classification_file),
                "--output", str(output),
                "--quiet",
            ],
        )
        assert result.exit_code == 0
        assert output.exists()

    def test_summary_invalid_format(
        self, cli_runner, temp_classification_file, temp_dir
    ):
        """Test error on invalid output format for summary."""
        output = temp_dir / "summary.xyz"
        result = cli_runner.invoke(
            app,
            [
                "visualize", "summary",
                "--classifications", str(temp_classification_file),
                "--output", str(output),
                "--format", "xyz",
                "--quiet",
            ],
        )
        assert result.exit_code != 0
        assert "format" in result.output.lower() or "invalid" in result.output.lower()

    def test_summary_custom_dimensions(
        self, cli_runner, temp_classification_file, temp_dir
    ):
        """Test summary plot with custom width and height."""
        output = temp_dir / "summary.html"
        result = cli_runner.invoke(
            app,
            [
                "visualize", "summary",
                "--classifications", str(temp_classification_file),
                "--output", str(output),
                "--width", "800",
                "--height", "400",
                "--quiet",
            ],
        )
        assert result.exit_code == 0
        assert output.exists()

    def test_summary_output_directory_created(
        self, cli_runner, temp_classification_file, temp_dir
    ):
        """Test that output directory is created when it does not exist."""
        output = temp_dir / "subdir" / "summary.html"
        result = cli_runner.invoke(
            app,
            [
                "visualize", "summary",
                "--classifications", str(temp_classification_file),
                "--output", str(output),
                "--quiet",
            ],
        )
        assert result.exit_code == 0
        assert output.exists()

    def test_summary_read_error(self, cli_runner, temp_dir):
        """Test error handling when classification file is unreadable."""
        bad_file = temp_dir / "bad.csv"
        bad_file.write_text("not,valid,csv\nwith\nmismatched,columns,a,b")
        output = temp_dir / "summary.html"
        result = cli_runner.invoke(
            app,
            [
                "visualize", "summary",
                "--classifications", str(bad_file),
                "--output", str(output),
                "--quiet",
            ],
        )
        # Should fail because 'taxonomic_call' column is missing
        assert result.exit_code != 0


# =============================================================================
# Visualize: help text tests
# =============================================================================


class TestVisualizeHelpText:
    """Tests for visualize CLI help text."""

    def test_visualize_help(self, cli_runner):
        """Test that visualize top-level help is shown."""
        result = cli_runner.invoke(app, ["visualize", "--help"])
        assert result.exit_code == 0
        assert "recruitment" in result.output.lower()
        assert "summary" in result.output.lower()

    def test_recruitment_help(self, cli_runner):
        """Test that recruitment help shows expected options."""
        result = cli_runner.invoke(app, ["visualize", "recruitment", "--help"])
        assert result.exit_code == 0
        assert "--bam" in result.output
        assert "--output" in result.output
        assert "--format" in result.output
        assert "--genome" in result.output
        assert "--export-tsv" in result.output
        assert "--multi-panel" in result.output

    def test_summary_help(self, cli_runner):
        """Test that summary help shows expected options."""
        result = cli_runner.invoke(app, ["visualize", "summary", "--help"])
        assert result.exit_code == 0
        assert "--classifications" in result.output
        assert "--output" in result.output
        assert "--format" in result.output


# =============================================================================
# MMseqs2: makedb command tests
# =============================================================================


class TestMMseqs2MakedbValidation:
    """Tests for mmseqs2 makedb command argument validation."""

    def test_makedb_requires_genomes(self, cli_runner, temp_dir):
        """Test that --genomes argument is required."""
        output = temp_dir / "mmseqs_db" / "pangenome"
        result = cli_runner.invoke(
            app,
            ["mmseqs2", "makedb", "--output", str(output)],
        )
        assert result.exit_code != 0
        assert "--genomes" in result.output or "genomes" in result.output.lower()

    def test_makedb_requires_output(self, cli_runner, temp_genome_dir):
        """Test that --output argument is required."""
        result = cli_runner.invoke(
            app,
            ["mmseqs2", "makedb", "--genomes", str(temp_genome_dir)],
        )
        assert result.exit_code != 0
        assert "--output" in result.output or "output" in result.output.lower()

    def test_makedb_missing_genomes_path(self, cli_runner, temp_dir):
        """Test error handling for non-existent genomes path."""
        output = temp_dir / "mmseqs_db" / "pangenome"
        result = cli_runner.invoke(
            app,
            [
                "mmseqs2", "makedb",
                "--genomes", str(temp_dir / "nonexistent"),
                "--output", str(output),
            ],
        )
        assert result.exit_code != 0

    def test_makedb_mmseqs2_not_available(
        self, cli_runner, temp_genome_dir, temp_dir
    ):
        """Test error when mmseqs2 is not installed."""
        output = temp_dir / "mmseqs_db" / "pangenome"
        with patch(
            "metadarkmatter.cli.mmseqs2.MMseqs2.check_available",
            return_value=False,
        ):
            result = cli_runner.invoke(
                app,
                [
                    "mmseqs2", "makedb",
                    "--genomes", str(temp_genome_dir),
                    "--output", str(output),
                ],
            )
        assert result.exit_code != 0
        assert "mmseqs2" in result.output.lower()

    @patch("metadarkmatter.cli.mmseqs2.MMseqs2.check_available", return_value=True)
    def test_makedb_skip_if_exists(
        self, mock_avail, cli_runner, temp_dir
    ):
        """Test --skip-if-exists when database already exists."""
        db_dir = temp_dir / "mmseqs_db"
        db_dir.mkdir()
        db_path = db_dir / "pangenome"
        # Create all expected database files
        db_path.write_bytes(b"\x00")
        db_path.with_suffix(".dbtype").write_bytes(b"\x00")
        db_path.with_suffix(".index").write_bytes(b"\x00")

        genome_dir = temp_dir / "genomes"
        genome_dir.mkdir()
        (genome_dir / "test.fna").write_text(">seq\nACGT\n")

        result = cli_runner.invoke(
            app,
            [
                "mmseqs2", "makedb",
                "--genomes", str(genome_dir),
                "--output", str(db_path),
                "--skip-if-exists",
            ],
        )
        assert result.exit_code == 0
        assert "skipping" in result.output.lower()

    @patch("metadarkmatter.cli.mmseqs2.MMseqs2.check_available", return_value=True)
    @patch("metadarkmatter.cli.mmseqs2._concatenate_genomes")
    @patch("metadarkmatter.cli.mmseqs2.MMseqs2.run_or_raise")
    def test_makedb_from_directory_happy_path(
        self,
        mock_run,
        mock_concat,
        mock_avail,
        cli_runner,
        temp_genome_dir,
        temp_dir,
    ):
        """Test successful database creation from a genome directory."""
        mock_concat.return_value = (3, 15)
        mock_run.return_value = _make_tool_result(elapsed=2.5)

        output = temp_dir / "mmseqs_db" / "pangenome"
        result = cli_runner.invoke(
            app,
            [
                "mmseqs2", "makedb",
                "--genomes", str(temp_genome_dir),
                "--output", str(output),
                "--quiet",
            ],
        )
        assert result.exit_code == 0
        mock_concat.assert_called_once()
        mock_run.assert_called_once()

    @patch("metadarkmatter.cli.mmseqs2.MMseqs2.check_available", return_value=True)
    @patch("metadarkmatter.cli.mmseqs2.MMseqs2.run_or_raise")
    def test_makedb_from_single_fasta(
        self, mock_run, mock_avail, cli_runner, temp_single_fasta, temp_dir
    ):
        """Test database creation from a single FASTA file."""
        mock_run.return_value = _make_tool_result(elapsed=1.0)

        output = temp_dir / "mmseqs_db" / "pangenome"
        result = cli_runner.invoke(
            app,
            [
                "mmseqs2", "makedb",
                "--genomes", str(temp_single_fasta),
                "--output", str(output),
                "--quiet",
            ],
        )
        assert result.exit_code == 0
        mock_run.assert_called_once()

    @patch("metadarkmatter.cli.mmseqs2.MMseqs2.check_available", return_value=True)
    @patch("metadarkmatter.cli.mmseqs2._concatenate_genomes")
    def test_makedb_no_genomes_found(
        self, mock_concat, mock_avail, cli_runner, temp_genome_dir, temp_dir
    ):
        """Test error when no genome files match the pattern."""
        mock_concat.return_value = (0, 0)

        output = temp_dir / "mmseqs_db" / "pangenome"
        result = cli_runner.invoke(
            app,
            [
                "mmseqs2", "makedb",
                "--genomes", str(temp_genome_dir),
                "--output", str(output),
            ],
        )
        assert result.exit_code != 0
        assert "no genome" in result.output.lower()

    @patch("metadarkmatter.cli.mmseqs2.MMseqs2.check_available", return_value=True)
    @patch("metadarkmatter.cli.mmseqs2._concatenate_genomes")
    @patch("metadarkmatter.cli.mmseqs2.MMseqs2.run_or_raise")
    def test_makedb_tool_execution_error(
        self,
        mock_run,
        mock_concat,
        mock_avail,
        cli_runner,
        temp_genome_dir,
        temp_dir,
    ):
        """Test error handling when MMseqs2 createdb fails."""
        mock_concat.return_value = (3, 15)
        mock_run.side_effect = ToolExecutionError(
            "mmseqs", ["mmseqs", "createdb"], 1, "createdb failed"
        )

        output = temp_dir / "mmseqs_db" / "pangenome"
        result = cli_runner.invoke(
            app,
            [
                "mmseqs2", "makedb",
                "--genomes", str(temp_genome_dir),
                "--output", str(output),
            ],
        )
        assert result.exit_code != 0
        assert "failed" in result.output.lower()

    @patch("metadarkmatter.cli.mmseqs2.MMseqs2.check_available", return_value=True)
    @patch("metadarkmatter.cli.mmseqs2.MMseqs2.run")
    def test_makedb_dry_run_from_directory(
        self, mock_run, mock_avail, cli_runner, temp_genome_dir, temp_dir
    ):
        """Test dry-run mode shows commands without executing."""
        mock_run.return_value = _make_tool_result(
            command=("mmseqs", "createdb", "input.fasta", "db")
        )
        output = temp_dir / "mmseqs_db" / "pangenome"
        result = cli_runner.invoke(
            app,
            [
                "mmseqs2", "makedb",
                "--genomes", str(temp_genome_dir),
                "--output", str(output),
                "--dry-run",
            ],
        )
        assert result.exit_code == 0
        assert "dry run" in result.output.lower()

    @patch("metadarkmatter.cli.mmseqs2.MMseqs2.check_available", return_value=True)
    @patch("metadarkmatter.cli.mmseqs2.MMseqs2.run")
    def test_makedb_dry_run_from_fasta(
        self, mock_run, mock_avail, cli_runner, temp_single_fasta, temp_dir
    ):
        """Test dry-run mode from a single FASTA file."""
        mock_run.return_value = _make_tool_result(
            command=("mmseqs", "createdb", "in.fasta", "db")
        )
        output = temp_dir / "mmseqs_db" / "pangenome"
        result = cli_runner.invoke(
            app,
            [
                "mmseqs2", "makedb",
                "--genomes", str(temp_single_fasta),
                "--output", str(output),
                "--dry-run",
            ],
        )
        assert result.exit_code == 0
        assert "dry run" in result.output.lower()


# =============================================================================
# MMseqs2: search command tests
# =============================================================================


class TestMMseqs2SearchValidation:
    """Tests for mmseqs2 search command argument validation."""

    def test_search_requires_query_or_paired(self, cli_runner, temp_dir, temp_mmseqs_db):
        """Test that at least one query input is required."""
        output = temp_dir / "results.tsv.gz"
        with patch(
            "metadarkmatter.cli.mmseqs2.MMseqs2.check_available",
            return_value=True,
        ):
            result = cli_runner.invoke(
                app,
                [
                    "mmseqs2", "search",
                    "--database", str(temp_mmseqs_db),
                    "--output", str(output),
                ],
            )
        assert result.exit_code != 0
        assert "must provide" in result.output.lower() or "error" in result.output.lower()

    def test_search_mutual_exclusivity(
        self, cli_runner, temp_query_fasta, temp_query_fastq_r1,
        temp_query_fastq_r2, temp_mmseqs_db, temp_dir,
    ):
        """Test that --query and --query-1/--query-2 are mutually exclusive."""
        output = temp_dir / "results.tsv.gz"
        with patch(
            "metadarkmatter.cli.mmseqs2.MMseqs2.check_available",
            return_value=True,
        ):
            result = cli_runner.invoke(
                app,
                [
                    "mmseqs2", "search",
                    "--query", str(temp_query_fasta),
                    "--query-1", str(temp_query_fastq_r1),
                    "--query-2", str(temp_query_fastq_r2),
                    "--database", str(temp_mmseqs_db),
                    "--output", str(output),
                ],
            )
        assert result.exit_code != 0
        assert "cannot" in result.output.lower()

    def test_search_paired_requires_both(
        self, cli_runner, temp_query_fastq_r1, temp_mmseqs_db, temp_dir
    ):
        """Test that paired-end mode requires both --query-1 and --query-2."""
        output = temp_dir / "results.tsv.gz"
        with patch(
            "metadarkmatter.cli.mmseqs2.MMseqs2.check_available",
            return_value=True,
        ):
            result = cli_runner.invoke(
                app,
                [
                    "mmseqs2", "search",
                    "--query-1", str(temp_query_fastq_r1),
                    "--database", str(temp_mmseqs_db),
                    "--output", str(output),
                ],
            )
        assert result.exit_code != 0
        assert "both" in result.output.lower() or "--query-2" in result.output

    def test_search_mmseqs2_not_available(
        self, cli_runner, temp_query_fasta, temp_mmseqs_db, temp_dir
    ):
        """Test error when mmseqs2 is not installed."""
        output = temp_dir / "results.tsv.gz"
        with patch(
            "metadarkmatter.cli.mmseqs2.MMseqs2.check_available",
            return_value=False,
        ):
            result = cli_runner.invoke(
                app,
                [
                    "mmseqs2", "search",
                    "--query", str(temp_query_fasta),
                    "--database", str(temp_mmseqs_db),
                    "--output", str(output),
                ],
            )
        assert result.exit_code != 0
        assert "mmseqs2" in result.output.lower()

    def test_search_database_not_found(
        self, cli_runner, temp_query_fasta, temp_dir
    ):
        """Test error when database path does not exist."""
        output = temp_dir / "results.tsv.gz"
        with patch(
            "metadarkmatter.cli.mmseqs2.MMseqs2.check_available",
            return_value=True,
        ):
            result = cli_runner.invoke(
                app,
                [
                    "mmseqs2", "search",
                    "--query", str(temp_query_fasta),
                    "--database", str(temp_dir / "nonexistent_db"),
                    "--output", str(output),
                ],
            )
        assert result.exit_code != 0
        assert "not found" in result.output.lower() or "database" in result.output.lower()

    @patch("metadarkmatter.cli.mmseqs2.MMseqs2.check_available", return_value=True)
    def test_search_skip_if_exists(
        self, mock_avail, cli_runner, temp_query_fasta, temp_mmseqs_db, temp_dir
    ):
        """Test --skip-if-exists when output file already exists."""
        output = temp_dir / "results.tsv.gz"
        output.write_bytes(b"\x00")

        result = cli_runner.invoke(
            app,
            [
                "mmseqs2", "search",
                "--query", str(temp_query_fasta),
                "--database", str(temp_mmseqs_db),
                "--output", str(output),
                "--skip-if-exists",
            ],
        )
        assert result.exit_code == 0
        assert "skipping" in result.output.lower()

    @patch("metadarkmatter.cli.mmseqs2.MMseqs2.check_available", return_value=True)
    def test_search_dry_run_single_end(
        self, mock_avail, cli_runner, temp_query_fasta, temp_mmseqs_db, temp_dir
    ):
        """Test dry-run mode for single-end search."""
        output = temp_dir / "results.tsv.gz"
        result = cli_runner.invoke(
            app,
            [
                "mmseqs2", "search",
                "--query", str(temp_query_fasta),
                "--database", str(temp_mmseqs_db),
                "--output", str(output),
                "--dry-run",
            ],
        )
        assert result.exit_code == 0
        assert "dry run" in result.output.lower()

    @patch("metadarkmatter.cli.mmseqs2.MMseqs2.check_available", return_value=True)
    def test_search_dry_run_paired_end(
        self,
        mock_avail,
        cli_runner,
        temp_query_fastq_r1,
        temp_query_fastq_r2,
        temp_mmseqs_db,
        temp_dir,
    ):
        """Test dry-run mode for paired-end search."""
        output = temp_dir / "results.tsv.gz"
        result = cli_runner.invoke(
            app,
            [
                "mmseqs2", "search",
                "--query-1", str(temp_query_fastq_r1),
                "--query-2", str(temp_query_fastq_r2),
                "--database", str(temp_mmseqs_db),
                "--output", str(output),
                "--dry-run",
            ],
        )
        assert result.exit_code == 0
        assert "dry run" in result.output.lower()

    @patch("metadarkmatter.cli.mmseqs2.MMseqs2.check_available", return_value=True)
    @patch("metadarkmatter.cli.mmseqs2.MMseqs2.search_multistep")
    def test_search_happy_path_single_end(
        self,
        mock_multistep,
        mock_avail,
        cli_runner,
        temp_query_fasta,
        temp_mmseqs_db,
        temp_dir,
    ):
        """Test successful single-end search."""
        # Prepare output file so counts and size work
        output = temp_dir / "results.tsv.gz"
        temp_tsv = output.with_suffix("").with_suffix(".tsv")

        def fake_search(**kwargs):
            # Write dummy TSV output that the command expects
            tsv_path = kwargs.get("output", temp_tsv)
            Path(tsv_path).write_text(
                "read_1\tgenome_1\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-10\t200\t150\n"
                "read_2\tgenome_1\t90.0\t100\t10\t0\t1\t100\t200\t300\t1e-8\t180\t150\n"
            )
            return {
                "query_db": str(temp_dir / "queryDB"),
                "result_db": str(temp_dir / "resultDB"),
                "createdb_time": 0.5,
                "search_time": 2.0,
                "convertalis_time": 0.3,
                "total_time": 2.8,
            }

        mock_multistep.side_effect = fake_search

        result = cli_runner.invoke(
            app,
            [
                "mmseqs2", "search",
                "--query", str(temp_query_fasta),
                "--database", str(temp_mmseqs_db),
                "--output", str(output),
                "--quiet",
            ],
        )
        assert result.exit_code == 0
        mock_multistep.assert_called_once()

    @patch("metadarkmatter.cli.mmseqs2.MMseqs2.check_available", return_value=True)
    @patch("metadarkmatter.cli.mmseqs2.MMseqs2.search_multistep")
    @patch("metadarkmatter.cli.mmseqs2._concatenate_paired_reads")
    def test_search_happy_path_paired_end(
        self,
        mock_concat_reads,
        mock_multistep,
        mock_avail,
        cli_runner,
        temp_query_fastq_r1,
        temp_query_fastq_r2,
        temp_mmseqs_db,
        temp_dir,
    ):
        """Test successful paired-end search."""
        mock_concat_reads.return_value = 4
        output = temp_dir / "results.tsv.gz"
        temp_tsv = output.with_suffix("").with_suffix(".tsv")

        def fake_search(**kwargs):
            tsv_path = kwargs.get("output", temp_tsv)
            Path(tsv_path).write_text(
                "read_1\tgenome_1\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-10\t200\t150\n"
            )
            return {
                "query_db": str(temp_dir / "queryDB"),
                "result_db": str(temp_dir / "resultDB"),
                "createdb_time": 0.5,
                "search_time": 2.0,
                "convertalis_time": 0.3,
                "total_time": 2.8,
            }

        mock_multistep.side_effect = fake_search

        result = cli_runner.invoke(
            app,
            [
                "mmseqs2", "search",
                "--query-1", str(temp_query_fastq_r1),
                "--query-2", str(temp_query_fastq_r2),
                "--database", str(temp_mmseqs_db),
                "--output", str(output),
                "--quiet",
            ],
        )
        assert result.exit_code == 0

    @patch("metadarkmatter.cli.mmseqs2.MMseqs2.check_available", return_value=True)
    @patch("metadarkmatter.cli.mmseqs2.MMseqs2.search_multistep")
    def test_search_execution_error(
        self,
        mock_multistep,
        mock_avail,
        cli_runner,
        temp_query_fasta,
        temp_mmseqs_db,
        temp_dir,
    ):
        """Test error handling when MMseqs2 search fails."""
        mock_multistep.side_effect = RuntimeError("Search failed")
        output = temp_dir / "results.tsv.gz"

        result = cli_runner.invoke(
            app,
            [
                "mmseqs2", "search",
                "--query", str(temp_query_fasta),
                "--database", str(temp_mmseqs_db),
                "--output", str(output),
            ],
        )
        assert result.exit_code != 0
        assert "failed" in result.output.lower()

    @patch("metadarkmatter.cli.mmseqs2.MMseqs2.check_available", return_value=True)
    @patch("metadarkmatter.cli.mmseqs2.MMseqs2.search_multistep")
    def test_search_no_compress(
        self,
        mock_multistep,
        mock_avail,
        cli_runner,
        temp_query_fasta,
        temp_mmseqs_db,
        temp_dir,
    ):
        """Test search with --no-compress flag."""
        output = temp_dir / "results.tsv"

        def fake_search(**kwargs):
            tsv_path = kwargs.get("output", output)
            Path(tsv_path).write_text(
                "read_1\tgenome_1\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-10\t200\t150\n"
            )
            return {
                "query_db": str(temp_dir / "queryDB"),
                "result_db": str(temp_dir / "resultDB"),
                "createdb_time": 0.5,
                "search_time": 2.0,
                "convertalis_time": 0.3,
                "total_time": 2.8,
            }

        mock_multistep.side_effect = fake_search

        result = cli_runner.invoke(
            app,
            [
                "mmseqs2", "search",
                "--query", str(temp_query_fasta),
                "--database", str(temp_mmseqs_db),
                "--output", str(output),
                "--no-compress",
                "--quiet",
            ],
        )
        assert result.exit_code == 0


# =============================================================================
# MMseqs2: help text tests
# =============================================================================


class TestMMseqs2HelpText:
    """Tests for mmseqs2 CLI help text."""

    def test_mmseqs2_help(self, cli_runner):
        """Test that mmseqs2 top-level help is shown."""
        result = cli_runner.invoke(app, ["mmseqs2", "--help"])
        assert result.exit_code == 0
        assert "makedb" in result.output.lower()
        assert "search" in result.output.lower()

    def test_makedb_help(self, cli_runner):
        """Test that makedb help shows expected options."""
        result = cli_runner.invoke(app, ["mmseqs2", "makedb", "--help"])
        assert result.exit_code == 0
        assert "--genomes" in result.output
        assert "--output" in result.output
        assert "--skip-if-exists" in result.output
        assert "--dry-run" in result.output
        assert "--genome-pattern" in result.output

    def test_search_help(self, cli_runner):
        """Test that search help shows expected options."""
        result = cli_runner.invoke(app, ["mmseqs2", "search", "--help"])
        assert result.exit_code == 0
        assert "--query" in result.output
        assert "--query-1" in result.output
        assert "--query-2" in result.output
        assert "--database" in result.output
        assert "--output" in result.output
        assert "--sensitivity" in result.output
        assert "--evalue" in result.output
        assert "--threads" in result.output
        assert "--dry-run" in result.output
        assert "--skip-if-exists" in result.output
        assert "--compress" in result.output


# =============================================================================
# MMseqs2: _concatenate_paired_reads unit tests
# =============================================================================


class TestConcatenatePairedReads:
    """Tests for the _concatenate_paired_reads helper function."""

    def test_concatenate_plain_fastq(
        self, temp_query_fastq_r1, temp_query_fastq_r2, temp_dir
    ):
        """Test concatenation of plain FASTQ files."""
        from rich.console import Console

        from metadarkmatter.cli.mmseqs2 import _concatenate_paired_reads

        output_path = temp_dir / "combined.fastq"
        console = Console()
        count = _concatenate_paired_reads(
            temp_query_fastq_r1,
            temp_query_fastq_r2,
            output_path,
            console,
            quiet=True,
        )
        assert output_path.exists()
        # 2 reads from R1 + 2 reads from R2 = 4 reads
        assert count == 4

    def test_concatenate_gzipped_fastq(
        self, temp_query_fastq_gz_r1, temp_query_fastq_gz_r2, temp_dir
    ):
        """Test concatenation of gzipped FASTQ files."""
        from rich.console import Console

        from metadarkmatter.cli.mmseqs2 import _concatenate_paired_reads

        output_path = temp_dir / "combined.fastq"
        console = Console()
        count = _concatenate_paired_reads(
            temp_query_fastq_gz_r1,
            temp_query_fastq_gz_r2,
            output_path,
            console,
            quiet=True,
        )
        assert output_path.exists()
        assert count == 4

    def test_concatenate_fasta(self, temp_dir):
        """Test concatenation of plain FASTA files."""
        from rich.console import Console

        from metadarkmatter.cli.mmseqs2 import _concatenate_paired_reads

        r1 = temp_dir / "reads_R1.fasta"
        r2 = temp_dir / "reads_R2.fasta"
        r1.write_text(">read_1\nACGT\n>read_2\nGGGG\n")
        r2.write_text(">read_3\nTTTT\n>read_4\nCCCC\n")

        output_path = temp_dir / "combined.fasta"
        console = Console()
        count = _concatenate_paired_reads(
            r1, r2, output_path, console, quiet=True,
        )
        assert output_path.exists()
        # FASTA: 2 lines per read, 4 reads total
        assert count == 4


# =============================================================================
# MMseqs2: _check_mmseqs2_available unit tests
# =============================================================================


class TestCheckMMseqs2Available:
    """Tests for the _check_mmseqs2_available helper."""

    def test_available(self):
        """Test when mmseqs2 is available."""
        from metadarkmatter.cli.mmseqs2 import _check_mmseqs2_available

        with patch(
            "metadarkmatter.cli.mmseqs2.MMseqs2.check_available",
            return_value=True,
        ):
            ok, missing = _check_mmseqs2_available()
        assert ok is True
        assert missing == ""

    def test_not_available(self):
        """Test when mmseqs2 is not available."""
        from metadarkmatter.cli.mmseqs2 import _check_mmseqs2_available

        with patch(
            "metadarkmatter.cli.mmseqs2.MMseqs2.check_available",
            return_value=False,
        ):
            ok, missing = _check_mmseqs2_available()
        assert ok is False
        assert missing == "mmseqs2"
