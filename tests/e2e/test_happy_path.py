"""
E2E happy path tests for metadarkmatter CLI.

Tests basic functionality with valid inputs and expected outputs.
These tests verify the primary use cases work correctly.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable

import polars as pl
import pytest

from tests.utils.assertions import (
    CLIAssertions,
    ClassificationAssertions,
    SummaryAssertions,
)


pytestmark = pytest.mark.e2e


class TestClassifyHappyPath:
    """Happy path tests for the classify command."""

    def test_basic_classification_csv(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """Basic classification should produce valid CSV output."""
        blast_path, ani_path = standard_dataset

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="basic.csv",
        )

        CLIAssertions.assert_success(result)
        CLIAssertions.assert_output_contains(result, "Classified")

        # Validate output file
        output_path = e2e_temp_dir / "basic.csv"
        df = ClassificationAssertions.assert_all_validations(output_path, min_rows=1)
        assert len(df) > 0

    def test_basic_classification_parquet(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """Classification with Parquet format should produce valid output."""
        blast_path, ani_path = standard_dataset

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="basic.parquet",
            output_format="parquet",
        )

        CLIAssertions.assert_success(result)

        output_path = e2e_temp_dir / "basic.parquet"
        df = ClassificationAssertions.assert_all_validations(output_path, min_rows=1)
        assert len(df) > 0

    def test_classification_with_summary(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """Classification with --summary should produce both output and summary."""
        blast_path, ani_path = standard_dataset

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="with_summary.csv",
            summary_name="summary.json",
        )

        CLIAssertions.assert_success(result)

        # Validate both files
        output_path = e2e_temp_dir / "with_summary.csv"
        summary_path = e2e_temp_dir / "summary.json"

        df = ClassificationAssertions.assert_all_validations(output_path)
        summary = SummaryAssertions.assert_all_validations(summary_path, df)

        assert summary["total_reads"] == len(df)

    def test_classification_fast_mode(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """Fast mode should produce valid output."""
        blast_path, ani_path = standard_dataset

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="fast.csv",
            fast=True,
        )

        CLIAssertions.assert_success(result)

        output_path = e2e_temp_dir / "fast.csv"
        ClassificationAssertions.assert_all_validations(output_path, min_rows=1)

    def test_classification_parallel_mode(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """Parallel mode should produce valid output."""
        blast_path, ani_path = standard_dataset

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="parallel.csv",
            parallel=True,
        )

        CLIAssertions.assert_success(result)

        output_path = e2e_temp_dir / "parallel.csv"
        ClassificationAssertions.assert_all_validations(output_path, min_rows=1)

    def test_classification_streaming_mode(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """Streaming mode should produce valid output."""
        blast_path, ani_path = standard_dataset

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="streaming.csv",
            streaming=True,
        )

        CLIAssertions.assert_success(result)

        output_path = e2e_temp_dir / "streaming.csv"
        ClassificationAssertions.assert_all_validations(output_path, min_rows=1)

    def test_classification_verbose(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """Verbose mode should show additional information."""
        blast_path, ani_path = standard_dataset

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="verbose.csv",
            verbose=True,
        )

        CLIAssertions.assert_success(result)
        # Verbose should show genome coverage info
        CLIAssertions.assert_output_contains(
            result, "coverage", case_sensitive=False
        )

    def test_creates_output_directory(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """Should create output directory if it doesn't exist."""
        blast_path, ani_path = standard_dataset
        nested_output = e2e_temp_dir / "nested" / "deep" / "output.csv"

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="nested/deep/output.csv",
        )

        CLIAssertions.assert_success(result)
        assert (e2e_temp_dir / "nested" / "deep" / "output.csv").exists()


class TestBatchHappyPath:
    """Happy path tests for the batch command."""

    def test_batch_basic(
        self,
        batch_dataset: tuple[Path, Path],
        run_batch: Callable,
        e2e_temp_dir: Path,
    ):
        """Basic batch processing should handle all files."""
        blast_dir, ani_path = batch_dataset

        result = run_batch(
            blast_dir=blast_dir,
            ani=ani_path,
            output_dir_name="batch_basic",
        )

        CLIAssertions.assert_success(result)
        CLIAssertions.assert_output_contains(result, "Batch Processing Complete")

        # Verify output files created
        output_dir = e2e_temp_dir / "batch_basic"
        output_files = list(output_dir.glob("*_classifications.csv"))
        assert len(output_files) == 3

    def test_batch_with_summaries(
        self,
        batch_dataset: tuple[Path, Path],
        run_batch: Callable,
        e2e_temp_dir: Path,
    ):
        """Batch should create summary files for each sample."""
        blast_dir, ani_path = batch_dataset

        result = run_batch(
            blast_dir=blast_dir,
            ani=ani_path,
            output_dir_name="batch_summaries",
        )

        CLIAssertions.assert_success(result)

        output_dir = e2e_temp_dir / "batch_summaries"
        summary_files = list(output_dir.glob("*_summary.json"))
        assert len(summary_files) == 3

        # Validate each summary
        for summary_path in summary_files:
            SummaryAssertions.assert_valid_summary_file(summary_path)

    def test_batch_parquet_format(
        self,
        batch_dataset: tuple[Path, Path],
        run_batch: Callable,
        e2e_temp_dir: Path,
    ):
        """Batch with Parquet format should work correctly."""
        blast_dir, ani_path = batch_dataset

        result = run_batch(
            blast_dir=blast_dir,
            ani=ani_path,
            output_dir_name="batch_parquet",
            output_format="parquet",
        )

        CLIAssertions.assert_success(result)

        output_dir = e2e_temp_dir / "batch_parquet"
        parquet_files = list(output_dir.glob("*_classifications.parquet"))
        assert len(parquet_files) == 3

    def test_batch_fast_mode(
        self,
        batch_dataset: tuple[Path, Path],
        run_batch: Callable,
        e2e_temp_dir: Path,
    ):
        """Batch with fast mode should work correctly."""
        blast_dir, ani_path = batch_dataset

        result = run_batch(
            blast_dir=blast_dir,
            ani=ani_path,
            output_dir_name="batch_fast",
            fast=True,
        )

        CLIAssertions.assert_success(result)

        output_dir = e2e_temp_dir / "batch_fast"
        output_files = list(output_dir.glob("*_classifications.csv"))
        assert len(output_files) == 3

    def test_batch_parallel_mode(
        self,
        batch_dataset: tuple[Path, Path],
        run_batch: Callable,
        e2e_temp_dir: Path,
    ):
        """Batch with parallel mode should work correctly."""
        blast_dir, ani_path = batch_dataset

        result = run_batch(
            blast_dir=blast_dir,
            ani=ani_path,
            output_dir_name="batch_parallel",
            parallel=True,
        )

        CLIAssertions.assert_success(result)
        CLIAssertions.assert_output_contains(result, "vectorized")


class TestOutputValidation:
    """Tests validating the content of classification outputs."""

    def test_all_taxonomic_calls_present(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """Standard dataset should produce all taxonomic call types."""
        blast_path, ani_path = standard_dataset

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="all_calls.csv",
            summary_name="all_calls_summary.json",
        )

        CLIAssertions.assert_success(result)

        output_path = e2e_temp_dir / "all_calls.csv"
        df = ClassificationAssertions.assert_all_validations(output_path)

        # Check that we have diversity in classifications
        unique_calls = set(df["taxonomic_call"].unique().to_list())
        assert len(unique_calls) >= 2, (
            f"Expected multiple classification types, got: {unique_calls}"
        )

    def test_csv_parquet_equivalence(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """CSV and Parquet outputs should be equivalent."""
        blast_path, ani_path = standard_dataset

        # Run with CSV
        csv_result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="equiv.csv",
            output_format="csv",
        )
        CLIAssertions.assert_success(csv_result)

        # Run with Parquet
        parquet_result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="equiv.parquet",
            output_format="parquet",
        )
        CLIAssertions.assert_success(parquet_result)

        # Compare outputs
        csv_df = pl.read_csv(e2e_temp_dir / "equiv.csv")
        parquet_df = pl.read_parquet(e2e_temp_dir / "equiv.parquet")

        assert len(csv_df) == len(parquet_df)
        assert set(csv_df.columns) == set(parquet_df.columns)

    def test_summary_matches_output(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """Summary statistics should match the classification output."""
        blast_path, ani_path = standard_dataset

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="summary_check.csv",
            summary_name="summary_check.json",
        )

        CLIAssertions.assert_success(result)

        output_path = e2e_temp_dir / "summary_check.csv"
        summary_path = e2e_temp_dir / "summary_check.json"

        df = ClassificationAssertions.assert_all_validations(output_path)
        summary = SummaryAssertions.assert_all_validations(summary_path, df)

        # Cross-validate
        assert summary["total_reads"] == len(df)
