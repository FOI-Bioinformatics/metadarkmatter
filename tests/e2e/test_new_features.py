"""
E2E tests for new CLI features.

Tests for:
- Format/extension mismatch correction
- Dry-run mode
- Quiet mode
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable

import pytest

from tests.utils.assertions import CLIAssertions


pytestmark = pytest.mark.e2e


class TestFormatExtensionMismatch:
    """Tests for format/extension mismatch correction."""

    def test_parquet_format_with_csv_extension_warns(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """Should warn and correct when --format parquet with .csv extension."""
        blast_path, ani_path = standard_dataset

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="wrong_ext.csv",  # Wrong extension for parquet
            output_format="parquet",
        )

        CLIAssertions.assert_success(result)
        CLIAssertions.assert_output_contains(
            result, "warning", case_sensitive=False
        )
        CLIAssertions.assert_output_contains(
            result, ".parquet", case_sensitive=False
        )

        # Should have created .parquet file, not .csv
        assert (e2e_temp_dir / "wrong_ext.parquet").exists()
        assert not (e2e_temp_dir / "wrong_ext.csv").exists()

    def test_csv_format_with_parquet_extension_warns(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """Should warn and correct when --format csv with .parquet extension."""
        blast_path, ani_path = standard_dataset

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="wrong_ext.parquet",  # Wrong extension for csv
            output_format="csv",
        )

        CLIAssertions.assert_success(result)
        CLIAssertions.assert_output_contains(
            result, "warning", case_sensitive=False
        )

        # Should have created .csv file, not .parquet
        assert (e2e_temp_dir / "wrong_ext.csv").exists()
        assert not (e2e_temp_dir / "wrong_ext.parquet").exists()

    def test_matching_extension_no_warning(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """Should not warn when extension matches format."""
        blast_path, ani_path = standard_dataset

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="correct.parquet",
            output_format="parquet",
        )

        CLIAssertions.assert_success(result)
        # Should not have extension mismatch warning
        CLIAssertions.assert_output_not_contains(
            result, "doesn't match", case_sensitive=False
        )


class TestDryRunMode:
    """Tests for dry-run mode."""

    def test_dry_run_validates_inputs(
        self,
        standard_dataset: tuple[Path, Path],
        e2e_runner,
        e2e_temp_dir: Path,
    ):
        """Dry-run should validate inputs and show summary."""
        from metadarkmatter.cli.main import app

        blast_path, ani_path = standard_dataset

        result = e2e_runner.invoke(app, [
            "score", "classify",
            "--blast", str(blast_path),
            "--ani", str(ani_path),
            "--output", str(e2e_temp_dir / "dry_run.csv"),
            "--dry-run",
        ])

        CLIAssertions.assert_success(result)
        CLIAssertions.assert_output_contains(result, "DRY RUN")
        CLIAssertions.assert_output_contains(result, "BLAST File")
        CLIAssertions.assert_output_contains(result, "ANI Matrix")
        CLIAssertions.assert_output_contains(result, "Genome Coverage")
        CLIAssertions.assert_output_contains(result, "Validation complete")

    def test_dry_run_does_not_create_output(
        self,
        standard_dataset: tuple[Path, Path],
        e2e_runner,
        e2e_temp_dir: Path,
    ):
        """Dry-run should not create output files."""
        from metadarkmatter.cli.main import app

        blast_path, ani_path = standard_dataset
        output_path = e2e_temp_dir / "should_not_exist.csv"

        result = e2e_runner.invoke(app, [
            "score", "classify",
            "--blast", str(blast_path),
            "--ani", str(ani_path),
            "--output", str(output_path),
            "--dry-run",
        ])

        CLIAssertions.assert_success(result)
        assert not output_path.exists(), "Dry-run should not create output file"

    def test_dry_run_shows_file_size(
        self,
        standard_dataset: tuple[Path, Path],
        e2e_runner,
        e2e_temp_dir: Path,
    ):
        """Dry-run should show BLAST file size."""
        from metadarkmatter.cli.main import app

        blast_path, ani_path = standard_dataset

        result = e2e_runner.invoke(app, [
            "score", "classify",
            "--blast", str(blast_path),
            "--ani", str(ani_path),
            "--output", str(e2e_temp_dir / "dry_run.csv"),
            "--dry-run",
        ])

        CLIAssertions.assert_success(result)
        CLIAssertions.assert_output_contains(result, "Size:")
        CLIAssertions.assert_output_contains(result, "MB")

    def test_dry_run_shows_processing_mode(
        self,
        standard_dataset: tuple[Path, Path],
        e2e_runner,
        e2e_temp_dir: Path,
    ):
        """Dry-run should show configured processing mode."""
        from metadarkmatter.cli.main import app

        blast_path, ani_path = standard_dataset

        result = e2e_runner.invoke(app, [
            "score", "classify",
            "--blast", str(blast_path),
            "--ani", str(ani_path),
            "--output", str(e2e_temp_dir / "dry_run.csv"),
            "--parallel",
            "--dry-run",
        ])

        CLIAssertions.assert_success(result)
        CLIAssertions.assert_output_contains(result, "Mode:")
        CLIAssertions.assert_output_contains(result, "parallel")


class TestQuietMode:
    """Tests for quiet mode."""

    def test_quiet_suppresses_progress_output(
        self,
        standard_dataset: tuple[Path, Path],
        e2e_runner,
        e2e_temp_dir: Path,
    ):
        """Quiet mode should suppress progress messages."""
        from metadarkmatter.cli.main import app

        blast_path, ani_path = standard_dataset

        result = e2e_runner.invoke(app, [
            "score", "classify",
            "--blast", str(blast_path),
            "--ani", str(ani_path),
            "--output", str(e2e_temp_dir / "quiet.csv"),
            "--quiet",
        ])

        CLIAssertions.assert_success(result)
        # Should not have banner or progress messages
        CLIAssertions.assert_output_not_contains(
            result, "Metadarkmatter ANI-Weighted Classification"
        )

    def test_quiet_still_creates_output(
        self,
        standard_dataset: tuple[Path, Path],
        e2e_runner,
        e2e_temp_dir: Path,
    ):
        """Quiet mode should still create output files."""
        from metadarkmatter.cli.main import app

        blast_path, ani_path = standard_dataset
        output_path = e2e_temp_dir / "quiet_output.csv"

        result = e2e_runner.invoke(app, [
            "score", "classify",
            "--blast", str(blast_path),
            "--ani", str(ani_path),
            "--output", str(output_path),
            "--quiet",
        ])

        CLIAssertions.assert_success(result)
        assert output_path.exists(), "Quiet mode should still create output"

    def test_quiet_errors_still_shown(
        self,
        e2e_temp_dir: Path,
        e2e_runner,
    ):
        """Quiet mode should still show error messages."""
        from metadarkmatter.cli.main import app

        result = e2e_runner.invoke(app, [
            "score", "classify",
            "--blast", str(e2e_temp_dir / "nonexistent.blast.tsv"),
            "--ani", str(e2e_temp_dir / "nonexistent.ani.csv"),
            "--output", str(e2e_temp_dir / "quiet.csv"),
            "--quiet",
        ])

        # Errors should still be shown
        CLIAssertions.assert_failure(result, expected_code=2)

    def test_quiet_with_summary(
        self,
        standard_dataset: tuple[Path, Path],
        e2e_runner,
        e2e_temp_dir: Path,
    ):
        """Quiet mode should still create summary file."""
        from metadarkmatter.cli.main import app

        blast_path, ani_path = standard_dataset
        output_path = e2e_temp_dir / "quiet_with_summary.csv"
        summary_path = e2e_temp_dir / "quiet_summary.json"

        result = e2e_runner.invoke(app, [
            "score", "classify",
            "--blast", str(blast_path),
            "--ani", str(ani_path),
            "--output", str(output_path),
            "--summary", str(summary_path),
            "--quiet",
        ])

        CLIAssertions.assert_success(result)
        assert output_path.exists()
        assert summary_path.exists()

        # Should not show summary table
        CLIAssertions.assert_output_not_contains(result, "Classification Summary")
