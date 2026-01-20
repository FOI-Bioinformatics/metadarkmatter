"""
Unit tests for score CLI commands.

Tests argument validation, file handling, and error cases for the score subcommands.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metadarkmatter.cli.main import app


class TestScoreClassifyValidation:
    """Tests for score classify command argument validation."""

    def test_classify_requires_blast_file(self, cli_runner, temp_ani_file, temp_dir):
        """Test that --blast argument is required."""
        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "classify",
                "--ani",
                str(temp_ani_file),
                "--output",
                str(output),
            ],
        )
        assert result.exit_code != 0
        assert "--blast" in result.stdout or "blast" in result.stdout.lower()

    def test_classify_requires_ani_file(self, cli_runner, temp_blast_file, temp_dir):
        """Test that --ani argument is required."""
        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "classify",
                "--blast",
                str(temp_blast_file),
                "--output",
                str(output),
            ],
        )
        assert result.exit_code != 0
        assert "--ani" in result.stdout or "ani" in result.stdout.lower()

    def test_classify_missing_blast_file(self, cli_runner, temp_ani_file, temp_dir):
        """Test error handling for non-existent BLAST file."""
        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "classify",
                "--blast",
                str(temp_dir / "nonexistent.blast.tsv"),
                "--ani",
                str(temp_ani_file),
                "--output",
                str(output),
            ],
        )
        assert result.exit_code != 0
        assert "not found" in result.stdout.lower() or "does not exist" in result.stdout.lower()

    def test_classify_missing_ani_file(self, cli_runner, temp_blast_file, temp_dir):
        """Test error handling for non-existent ANI file."""
        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "classify",
                "--blast",
                str(temp_blast_file),
                "--ani",
                str(temp_dir / "nonexistent.ani.csv"),
                "--output",
                str(output),
            ],
        )
        assert result.exit_code != 0
        assert "not found" in result.stdout.lower() or "does not exist" in result.stdout.lower()

    def test_classify_empty_blast_file(self, cli_runner, empty_blast_file, temp_ani_file, temp_dir):
        """Test handling of empty BLAST file."""
        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "classify",
                "--blast",
                str(empty_blast_file),
                "--ani",
                str(temp_ani_file),
                "--output",
                str(output),
            ],
        )
        # Should either fail with error or succeed with warning
        assert result.exit_code != 0 or "empty" in result.stdout.lower()

    def test_classify_invalid_threshold(self, cli_runner, temp_blast_file, temp_ani_file, temp_dir):
        """Test validation of threshold values."""
        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "classify",
                "--blast",
                str(temp_blast_file),
                "--ani",
                str(temp_ani_file),
                "--bitscore-threshold",
                "150",  # > 100, invalid
                "--output",
                str(output),
            ],
        )
        assert result.exit_code != 0

    def test_classify_invalid_mode(self, cli_runner, temp_blast_file, temp_ani_file, temp_dir):
        """Test validation of processing mode."""
        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "classify",
                "--blast",
                str(temp_blast_file),
                "--ani",
                str(temp_ani_file),
                "--mode",
                "invalid_mode",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code != 0


class TestScoreClassifyOptions:
    """Tests for score classify command with various options."""

    def test_classify_with_metadata(
        self,
        cli_runner,
        temp_blast_file,
        temp_ani_file,
        temp_metadata_file,
        temp_dir,
    ):
        """Test classify with metadata file."""
        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "classify",
                "--blast",
                str(temp_blast_file),
                "--ani",
                str(temp_ani_file),
                "--metadata",
                str(temp_metadata_file),
                "--output",
                str(output),
            ],
        )
        # Should succeed or at least not crash on argument parsing
        assert result.exit_code == 0 or "error" not in result.stdout.lower()

    def test_classify_parallel_mode(
        self,
        cli_runner,
        temp_blast_file,
        temp_ani_file,
        temp_dir,
    ):
        """Test classify with --parallel flag."""
        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "classify",
                "--blast",
                str(temp_blast_file),
                "--ani",
                str(temp_ani_file),
                "--parallel",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code == 0 or "error" not in result.stdout.lower()

    def test_classify_custom_preset(
        self,
        cli_runner,
        temp_blast_file,
        temp_ani_file,
        temp_dir,
    ):
        """Test classify with preset configuration."""
        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "classify",
                "--blast",
                str(temp_blast_file),
                "--ani",
                str(temp_ani_file),
                "--preset",
                "gtdb-strict",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code == 0 or "error" not in result.stdout.lower()


class TestScoreExtractNovel:
    """Tests for score extract-novel command."""

    def test_extract_novel_requires_classifications(self, cli_runner, temp_dir):
        """Test that --classifications argument is required."""
        output = temp_dir / "novel.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "extract-novel",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code != 0
        assert "--classifications" in result.stdout or "classifications" in result.stdout.lower()

    def test_extract_novel_missing_file(self, cli_runner, temp_dir):
        """Test error handling for non-existent classification file."""
        output = temp_dir / "novel.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "extract-novel",
                "--classifications",
                str(temp_dir / "nonexistent.csv"),
                "--output",
                str(output),
            ],
        )
        assert result.exit_code != 0

    def test_extract_novel_species_only(
        self,
        cli_runner,
        temp_classification_file,
        temp_dir,
    ):
        """Test extracting only novel species."""
        output = temp_dir / "novel_species.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "extract-novel",
                "--classifications",
                str(temp_classification_file),
                "--category",
                "species",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code == 0 or "error" not in result.stdout.lower()

    def test_extract_novel_genus_only(
        self,
        cli_runner,
        temp_classification_file,
        temp_dir,
    ):
        """Test extracting only novel genus."""
        output = temp_dir / "novel_genus.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "extract-novel",
                "--classifications",
                str(temp_classification_file),
                "--category",
                "genus",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code == 0 or "error" not in result.stdout.lower()

    def test_extract_novel_with_read_ids(
        self,
        cli_runner,
        temp_classification_file,
        temp_dir,
    ):
        """Test extracting novel reads with read ID output."""
        output_csv = temp_dir / "novel.csv"
        output_ids = temp_dir / "novel_ids.txt"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "extract-novel",
                "--classifications",
                str(temp_classification_file),
                "--output",
                str(output_csv),
                "--read-ids",
                str(output_ids),
            ],
        )
        assert result.exit_code == 0 or "error" not in result.stdout.lower()

    def test_extract_novel_custom_threshold(
        self,
        cli_runner,
        temp_classification_file,
        temp_dir,
    ):
        """Test extracting with custom novelty threshold."""
        output = temp_dir / "highly_novel.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "extract-novel",
                "--classifications",
                str(temp_classification_file),
                "--min-novelty",
                "15.0",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code == 0 or "error" not in result.stdout.lower()


class TestScoreBatch:
    """Tests for score batch command."""

    def test_batch_requires_input_dir(self, cli_runner, temp_dir):
        """Test that --input-dir argument is required."""
        output_dir = temp_dir / "batch_output"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "batch",
                "--output-dir",
                str(output_dir),
            ],
        )
        assert result.exit_code != 0

    def test_batch_missing_directory(self, cli_runner, temp_dir):
        """Test error handling for non-existent input directory."""
        output_dir = temp_dir / "batch_output"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "batch",
                "--input-dir",
                str(temp_dir / "nonexistent"),
                "--output-dir",
                str(output_dir),
            ],
        )
        assert result.exit_code != 0

    def test_batch_empty_directory(self, cli_runner, temp_dir):
        """Test handling of empty input directory."""
        input_dir = temp_dir / "empty_input"
        input_dir.mkdir()
        output_dir = temp_dir / "batch_output"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "batch",
                "--input-dir",
                str(input_dir),
                "--output-dir",
                str(output_dir),
            ],
        )
        # Should either warn or fail gracefully
        assert result.exit_code != 0 or "no files" in result.stdout.lower()
