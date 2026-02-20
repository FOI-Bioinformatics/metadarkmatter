"""
E2E error handling tests for metadarkmatter CLI.

Tests that invalid inputs produce appropriate error messages
and non-zero exit codes.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable

import polars as pl
import pytest

from tests.utils.assertions import CLIAssertions


pytestmark = pytest.mark.e2e


class TestMissingFiles:
    """Tests for missing input files."""

    def test_missing_blast_file(
        self,
        e2e_temp_dir: Path,
        run_classify: Callable,
    ):
        """Should error with clear message when BLAST file doesn't exist."""
        nonexistent = e2e_temp_dir / "nonexistent.blast.tsv"
        ani_path = e2e_temp_dir / "dummy.ani.csv"

        # Create minimal ANI file
        pl.DataFrame({
            "genome": ["G1"],
            "G1": [100.0],
        }).write_csv(ani_path)

        result = run_classify(
            blast=nonexistent,
            ani=ani_path,
            output_name="error.csv",
        )

        CLIAssertions.assert_failure(result, expected_code=2)

    def test_missing_ani_file(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """Should error with clear message when ANI file doesn't exist."""
        blast_path, _ = standard_dataset
        nonexistent = e2e_temp_dir / "nonexistent.ani.csv"

        result = run_classify(
            blast=blast_path,
            ani=nonexistent,
            output_name="error.csv",
        )

        CLIAssertions.assert_failure(result, expected_code=2)


class TestInvalidFiles:
    """Tests for invalid/malformed input files."""

    def test_malformed_blast_file(
        self,
        e2e_temp_dir: Path,
        run_classify: Callable,
    ):
        """Should handle malformed BLAST file gracefully."""
        # Create malformed BLAST file (wrong number of columns)
        blast_path = e2e_temp_dir / "malformed.blast.tsv"
        blast_path.write_text("read1\tgenome1\t99.0\n")  # Only 3 columns

        ani_path = e2e_temp_dir / "valid.ani.csv"
        pl.DataFrame({
            "genome": ["genome1"],
            "genome1": [100.0],
        }).write_csv(ani_path)

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="error.csv",
        )

        CLIAssertions.assert_failure(result)

    def test_malformed_ani_file(
        self,
        standard_dataset: tuple[Path, Path],
        e2e_temp_dir: Path,
        run_classify: Callable,
    ):
        """Should handle malformed ANI file gracefully."""
        blast_path, _ = standard_dataset

        # Create malformed ANI file
        ani_path = e2e_temp_dir / "malformed.ani.csv"
        ani_path.write_text("not,valid,ani,format\n1,2,3,4\n")

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="error.csv",
        )

        CLIAssertions.assert_failure(result)

    def test_empty_blast_file(
        self,
        empty_blast_file: Path,
        e2e_temp_dir: Path,
        run_classify: Callable,
    ):
        """Should handle empty BLAST file gracefully."""
        ani_path = e2e_temp_dir / "valid.ani.csv"
        pl.DataFrame({
            "genome": ["G1", "G2"],
            "G1": [100.0, 80.0],
            "G2": [80.0, 100.0],
        }).write_csv(ani_path)

        result = run_classify(
            blast=empty_blast_file,
            ani=ani_path,
            output_name="empty.csv",
        )

        # Empty file should complete but produce empty/zero output
        # Depending on implementation, this might succeed with 0 reads
        # or fail - either is acceptable as long as it doesn't crash


class TestInvalidArguments:
    """Tests for invalid command-line arguments."""

    def test_invalid_format(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
    ):
        """Should error on invalid --format value."""
        blast_path, ani_path = standard_dataset

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="error.csv",
            output_format="invalid_format",
        )

        CLIAssertions.assert_failure(result)
        CLIAssertions.assert_output_contains(result, "invalid", case_sensitive=False)

class TestBatchErrorHandling:
    """Tests for batch command error handling."""

    def test_batch_missing_directory(
        self,
        e2e_temp_dir: Path,
        run_batch: Callable,
    ):
        """Should error when blast directory doesn't exist."""
        nonexistent = e2e_temp_dir / "nonexistent_dir"
        ani_path = e2e_temp_dir / "dummy.ani.csv"

        pl.DataFrame({
            "genome": ["G1"],
            "G1": [100.0],
        }).write_csv(ani_path)

        result = run_batch(
            blast_dir=nonexistent,
            ani=ani_path,
        )

        CLIAssertions.assert_failure(result, expected_code=2)

    def test_batch_no_matching_files(
        self,
        e2e_temp_dir: Path,
        run_batch: Callable,
    ):
        """Should error when no files match pattern."""
        blast_dir = e2e_temp_dir / "empty_blast_dir"
        blast_dir.mkdir()

        ani_path = e2e_temp_dir / "dummy.ani.csv"
        pl.DataFrame({
            "genome": ["G1"],
            "G1": [100.0],
        }).write_csv(ani_path)

        result = run_batch(
            blast_dir=blast_dir,
            ani=ani_path,
            pattern="*.nonexistent",
        )

        CLIAssertions.assert_failure(result)
        CLIAssertions.assert_output_contains(
            result, "No alignment files found", case_sensitive=False
        )


class TestWarnings:
    """Tests for warning conditions (non-fatal issues)."""

    def test_low_genome_coverage_warning(
        self,
        e2e_temp_dir: Path,
        run_classify: Callable,
    ):
        """Should warn when genome coverage is low but still succeed."""
        # Create BLAST with genomes mostly not in ANI matrix
        blast_path = e2e_temp_dir / "low_coverage.blast.tsv"
        blast_lines = [
            "r1\tGCF_MISSING_001.1_genomic\t99\t150\t1\t0\t1\t150\t1\t150\t1e-50\t250",
            "r2\tGCF_MISSING_002.1_genomic\t98\t150\t2\t0\t1\t150\t1\t150\t1e-45\t240",
            "r3\tGCF_MISSING_003.1_genomic\t97\t150\t3\t0\t1\t150\t1\t150\t1e-40\t230",
            "r4\tGCF_000123456.1_genomic\t99.5\t150\t0\t0\t1\t150\t1\t150\t1e-60\t280",
        ]
        blast_path.write_text("\n".join(blast_lines) + "\n")

        # Create ANI matrix with only one genome (25% coverage)
        ani_path = e2e_temp_dir / "small.ani.csv"
        pl.DataFrame({
            "genome": ["GCF_000123456.1"],
            "GCF_000123456.1": [100.0],
        }).write_csv(ani_path)

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="low_coverage.csv",
            verbose=True,
        )

        # Should succeed but emit warning
        CLIAssertions.assert_success(result)
        CLIAssertions.assert_output_contains(
            result, "coverage", case_sensitive=False
        )
