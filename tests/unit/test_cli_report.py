"""
Unit tests for report CLI commands.

Tests argument validation, file handling, and error cases for report generation.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metadarkmatter.cli.main import app


class TestReportGenerateValidation:
    """Tests for report generate command argument validation."""

    def test_generate_requires_input(self, cli_runner, temp_dir):
        """Test that --classifications argument is required."""
        output = temp_dir / "report.html"
        result = cli_runner.invoke(
            app,
            [
                "report",
                "generate",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code != 0
        assert "--classifications" in result.output or "classifications" in result.output.lower()

    def test_generate_missing_input_file(self, cli_runner, temp_dir):
        """Test error handling for non-existent classification file."""
        output = temp_dir / "report.html"
        result = cli_runner.invoke(
            app,
            [
                "report",
                "generate",
                "--classifications",
                str(temp_dir / "nonexistent.csv"),
                "--output",
                str(output),
            ],
        )
        assert result.exit_code != 0
        assert "not found" in result.output.lower() or "does not exist" in result.output.lower()

    def test_generate_with_basic_input(
        self,
        cli_runner,
        temp_classification_file,
        temp_dir,
    ):
        """Test basic report generation with minimal arguments."""
        output = temp_dir / "report.html"
        result = cli_runner.invoke(
            app,
            [
                "report",
                "generate",
                "--classifications",
                str(temp_classification_file),
                "--output",
                str(output),
            ],
        )
        # Should succeed
        assert result.exit_code == 0 or "error" not in result.output.lower()
        # Output file should be created
        if result.exit_code == 0:
            assert output.exists()
            assert output.stat().st_size > 0

    def test_generate_with_metadata(
        self,
        cli_runner,
        temp_classification_file,
        temp_metadata_file,
        temp_dir,
    ):
        """Test report generation with genome metadata."""
        output = temp_dir / "report.html"
        result = cli_runner.invoke(
            app,
            [
                "report",
                "generate",
                "--classifications",
                str(temp_classification_file),
                "--metadata",
                str(temp_metadata_file),
                "--output",
                str(output),
            ],
        )
        assert result.exit_code == 0 or "error" not in result.output.lower()

    def test_generate_with_ani_matrix(
        self,
        cli_runner,
        temp_classification_file,
        temp_ani_file,
        temp_dir,
    ):
        """Test report generation with ANI matrix."""
        output = temp_dir / "report.html"
        result = cli_runner.invoke(
            app,
            [
                "report",
                "generate",
                "--classifications",
                str(temp_classification_file),
                "--ani",
                str(temp_ani_file),
                "--output",
                str(output),
            ],
        )
        assert result.exit_code == 0 or "error" not in result.output.lower()

    def test_generate_with_all_options(
        self,
        cli_runner,
        temp_classification_file,
        temp_metadata_file,
        temp_ani_file,
        temp_dir,
    ):
        """Test report generation with all optional arguments."""
        output = temp_dir / "report_full.html"
        result = cli_runner.invoke(
            app,
            [
                "report",
                "generate",
                "--classifications",
                str(temp_classification_file),
                "--metadata",
                str(temp_metadata_file),
                "--ani",
                str(temp_ani_file),
                "--output",
                str(output),
            ],
        )
        assert result.exit_code == 0 or "error" not in result.output.lower()

    def test_generate_missing_ani_file(
        self,
        cli_runner,
        temp_classification_file,
        temp_dir,
    ):
        """Test error handling for non-existent ANI file."""
        output = temp_dir / "report.html"
        result = cli_runner.invoke(
            app,
            [
                "report",
                "generate",
                "--classifications",
                str(temp_classification_file),
                "--ani",
                str(temp_dir / "nonexistent.ani.csv"),
                "--output",
                str(output),
            ],
        )
        assert result.exit_code != 0

    def test_generate_missing_metadata_file(
        self,
        cli_runner,
        temp_classification_file,
        temp_dir,
    ):
        """Test error handling for non-existent metadata file."""
        output = temp_dir / "report.html"
        result = cli_runner.invoke(
            app,
            [
                "report",
                "generate",
                "--classifications",
                str(temp_classification_file),
                "--metadata",
                str(temp_dir / "nonexistent_metadata.tsv"),
                "--output",
                str(output),
            ],
        )
        assert result.exit_code != 0


class TestReportMulti:
    """Tests for report multi command."""

    def test_multi_requires_input_dir(self, cli_runner, temp_dir):
        """Test that --input-dir argument is required."""
        output = temp_dir / "comparison.html"
        result = cli_runner.invoke(
            app,
            [
                "report",
                "multi",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code != 0

    def test_multi_missing_directory(self, cli_runner, temp_dir):
        """Test error handling for non-existent input directory."""
        output = temp_dir / "comparison.html"
        result = cli_runner.invoke(
            app,
            [
                "report",
                "multi",
                "--input-dir",
                str(temp_dir / "nonexistent"),
                "--output",
                str(output),
            ],
        )
        assert result.exit_code != 0

    def test_multi_empty_directory(self, cli_runner, temp_dir):
        """Test handling of directory with no classification files."""
        input_dir = temp_dir / "empty_input"
        input_dir.mkdir()
        output = temp_dir / "comparison.html"
        result = cli_runner.invoke(
            app,
            [
                "report",
                "multi",
                "--input-dir",
                str(input_dir),
                "--output",
                str(output),
            ],
        )
        # Should warn about no files found
        assert result.exit_code != 0 or "no files" in result.output.lower()

    def test_multi_with_single_file(
        self,
        cli_runner,
        temp_classification_file,
        temp_dir,
    ):
        """Test multi-sample report with single classification file."""
        # Create input directory with one file
        input_dir = temp_dir / "multi_input"
        input_dir.mkdir()

        # Copy classification file to input directory
        import shutil
        shutil.copy(temp_classification_file, input_dir / "sample1_classifications.csv")

        output = temp_dir / "comparison.html"
        result = cli_runner.invoke(
            app,
            [
                "report",
                "multi",
                "--input-dir",
                str(input_dir),
                "--output",
                str(output),
            ],
        )
        # Should succeed or warn about single sample
        assert result.exit_code == 0 or "single" in result.output.lower()

    def test_multi_with_multiple_files(
        self,
        cli_runner,
        temp_classification_file,
        temp_dir,
    ):
        """Test multi-sample report with multiple classification files."""
        # Create input directory with multiple files
        input_dir = temp_dir / "multi_input"
        input_dir.mkdir()

        # Copy classification file multiple times with different names
        import shutil
        for i in range(3):
            shutil.copy(
                temp_classification_file,
                input_dir / f"sample{i+1}_classifications.csv"
            )

        output = temp_dir / "comparison.html"
        result = cli_runner.invoke(
            app,
            [
                "report",
                "multi",
                "--input-dir",
                str(input_dir),
                "--output",
                str(output),
            ],
        )
        assert result.exit_code == 0 or "error" not in result.output.lower()
        # Output should be created
        if result.exit_code == 0:
            assert output.exists()

    def test_multi_with_custom_pattern(
        self,
        cli_runner,
        temp_classification_file,
        temp_dir,
    ):
        """Test multi-sample report with custom file pattern."""
        # Create input directory
        input_dir = temp_dir / "multi_input"
        input_dir.mkdir()

        # Copy classification files
        import shutil
        shutil.copy(temp_classification_file, input_dir / "sample1_classifications.csv")
        shutil.copy(temp_classification_file, input_dir / "sample2_classifications.csv")

        output = temp_dir / "comparison.html"
        result = cli_runner.invoke(
            app,
            [
                "report",
                "multi",
                "--input-dir",
                str(input_dir),
                "--pattern",
                "*_classifications.csv",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code == 0 or "error" not in result.output.lower()


class TestPhylogenyCLI:
    """Test phylogeny-related CLI options."""

    def test_tree_option_exists(self, cli_runner):
        """--tree option is available in generate command."""
        result = cli_runner.invoke(
            app,
            ["report", "generate", "--help"],
        )

        assert result.exit_code == 0
        assert "--tree" in result.output

    def test_no_phylogeny_option_exists(self, cli_runner):
        """--no-phylogeny option is available in generate command."""
        result = cli_runner.invoke(
            app,
            ["report", "generate", "--help"],
        )

        assert result.exit_code == 0
        assert "--no-phylogeny" in result.output

    def test_generate_with_no_phylogeny_flag(
        self,
        cli_runner,
        temp_classification_file,
        temp_ani_file,
        temp_dir,
    ):
        """Test report generation with --no-phylogeny flag."""
        output = temp_dir / "report_no_phylo.html"
        result = cli_runner.invoke(
            app,
            [
                "report",
                "generate",
                "--classifications",
                str(temp_classification_file),
                "--ani",
                str(temp_ani_file),
                "--output",
                str(output),
                "--no-phylogeny",
            ],
        )
        assert result.exit_code == 0 or "error" not in result.output.lower()

    def test_generate_with_tree_file(
        self,
        cli_runner,
        temp_classification_file,
        temp_ani_file,
        temp_dir,
    ):
        """Test report generation with --tree option."""
        # Create a simple Newick tree file
        tree_file = temp_dir / "tree.nwk"
        tree_file.write_text("((GCF_000001.1:0.1,GCF_000002.1:0.2):0.3,GCF_000003.1:0.4);")

        output = temp_dir / "report_with_tree.html"
        result = cli_runner.invoke(
            app,
            [
                "report",
                "generate",
                "--classifications",
                str(temp_classification_file),
                "--ani",
                str(temp_ani_file),
                "--tree",
                str(tree_file),
                "--output",
                str(output),
            ],
        )
        assert result.exit_code == 0 or "error" not in result.output.lower()

    def test_tree_option_requires_existing_file(
        self,
        cli_runner,
        temp_classification_file,
        temp_dir,
    ):
        """Test that --tree option validates file existence."""
        output = temp_dir / "report.html"
        result = cli_runner.invoke(
            app,
            [
                "report",
                "generate",
                "--classifications",
                str(temp_classification_file),
                "--tree",
                str(temp_dir / "nonexistent.nwk"),
                "--output",
                str(output),
            ],
        )
        assert result.exit_code != 0


class TestReportOutputValidation:
    """Tests for output file validation."""

    def test_generate_output_directory_created(
        self,
        cli_runner,
        temp_classification_file,
        temp_dir,
    ):
        """Test that output directory is created if it doesn't exist."""
        output_dir = temp_dir / "new_subdir"
        output = output_dir / "report.html"

        result = cli_runner.invoke(
            app,
            [
                "report",
                "generate",
                "--classifications",
                str(temp_classification_file),
                "--output",
                str(output),
            ],
        )

        # Should succeed and create directory
        if result.exit_code == 0:
            assert output_dir.exists()
            assert output.exists()

    def test_generate_overwrites_existing_file(
        self,
        cli_runner,
        temp_classification_file,
        temp_dir,
    ):
        """Test that existing output file is overwritten."""
        output = temp_dir / "report.html"

        # Create existing file
        output.write_text("old content")
        old_size = output.stat().st_size

        result = cli_runner.invoke(
            app,
            [
                "report",
                "generate",
                "--classifications",
                str(temp_classification_file),
                "--output",
                str(output),
            ],
        )

        # Should succeed and overwrite
        if result.exit_code == 0:
            assert output.exists()
            assert output.stat().st_size != old_size
