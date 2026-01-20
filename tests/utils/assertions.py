"""
Custom assertion helpers for E2E testing.

Provides domain-specific assertions for validating classification outputs,
summary statistics, and CLI execution results. These helpers encapsulate
common validation patterns and provide clear, actionable error messages.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import TYPE_CHECKING, Any

import polars as pl

if TYPE_CHECKING:
    from click.testing import Result


# =============================================================================
# Valid Values for Classification Outputs
# =============================================================================

VALID_TAXONOMIC_CALLS = frozenset({
    "Known Species",
    "Novel Species",
    "Novel Genus",
    "Conserved Region",
    "Ambiguous",
})

REQUIRED_CLASSIFICATION_COLUMNS = frozenset({
    "read_id",
    "best_match_genome",
    "novelty_index",
    "placement_uncertainty",
    "taxonomic_call",
    "is_novel",
})

REQUIRED_SUMMARY_FIELDS = frozenset({
    "total_reads",
    "known_species",
    "novel_species",
    "novel_genus",
    "conserved_regions",
    "mean_novelty_index",
    "mean_placement_uncertainty",
})


# =============================================================================
# Classification Output Assertions
# =============================================================================


class ClassificationAssertions:
    """
    Assertions for validating classification output files.

    Provides methods to validate:
    - File structure and required columns
    - Value ranges for numeric fields
    - Categorical field validity
    - Logical consistency between fields
    """

    @staticmethod
    def assert_valid_output_file(
        path: Path,
        min_rows: int = 0,
    ) -> pl.DataFrame:
        """
        Assert that a classification output file is valid and readable.

        Args:
            path: Path to classification output file (CSV or Parquet)
            min_rows: Minimum expected number of rows

        Returns:
            The loaded DataFrame for further assertions

        Raises:
            AssertionError: If file is invalid or unreadable
        """
        assert path.exists(), f"Classification file does not exist: {path}"

        # Load based on extension
        if path.suffix == ".parquet":
            df = pl.read_parquet(path)
        elif path.suffix == ".csv":
            df = pl.read_csv(path)
        else:
            raise AssertionError(
                f"Unsupported file format: {path.suffix}. "
                f"Expected .csv or .parquet"
            )

        assert len(df) >= min_rows, (
            f"Expected at least {min_rows} rows, got {len(df)}"
        )

        return df

    @staticmethod
    def assert_required_columns(df: pl.DataFrame) -> None:
        """
        Assert that all required classification columns are present.

        Args:
            df: Classification DataFrame

        Raises:
            AssertionError: If any required column is missing
        """
        actual_columns = set(df.columns)
        missing = REQUIRED_CLASSIFICATION_COLUMNS - actual_columns

        assert not missing, (
            f"Missing required columns: {sorted(missing)}. "
            f"Got columns: {sorted(actual_columns)}"
        )

    @staticmethod
    def assert_valid_taxonomic_calls(df: pl.DataFrame) -> None:
        """
        Assert that all taxonomic_call values are valid categories.

        Args:
            df: Classification DataFrame with 'taxonomic_call' column

        Raises:
            AssertionError: If any invalid taxonomic call is found
        """
        if "taxonomic_call" not in df.columns:
            raise AssertionError("DataFrame missing 'taxonomic_call' column")

        actual_calls = set(df["taxonomic_call"].unique().to_list())
        invalid = actual_calls - VALID_TAXONOMIC_CALLS

        assert not invalid, (
            f"Invalid taxonomic calls found: {sorted(invalid)}. "
            f"Valid values are: {sorted(VALID_TAXONOMIC_CALLS)}"
        )

    @staticmethod
    def assert_novelty_index_range(
        df: pl.DataFrame,
        min_val: float = 0.0,
        max_val: float = 100.0,
    ) -> None:
        """
        Assert that novelty_index values are within valid range.

        Novelty index is defined as 100 - percent_identity, so valid
        range is [0, 100] where 0 means identical match.

        Args:
            df: Classification DataFrame with 'novelty_index' column
            min_val: Minimum allowed value (default 0.0)
            max_val: Maximum allowed value (default 100.0)

        Raises:
            AssertionError: If any value is out of range
        """
        if "novelty_index" not in df.columns:
            raise AssertionError("DataFrame missing 'novelty_index' column")

        col = df["novelty_index"]
        actual_min = col.min()
        actual_max = col.max()

        assert actual_min >= min_val, (
            f"novelty_index has values below {min_val}: min={actual_min}"
        )
        assert actual_max <= max_val, (
            f"novelty_index has values above {max_val}: max={actual_max}"
        )

    @staticmethod
    def assert_placement_uncertainty_range(
        df: pl.DataFrame,
        min_val: float = 0.0,
        max_val: float = 100.0,
    ) -> None:
        """
        Assert that placement_uncertainty values are within valid range.

        Placement uncertainty is derived from ANI values, so valid
        range is [0, 100] where 0 means unambiguous placement.

        Args:
            df: Classification DataFrame with 'placement_uncertainty' column
            min_val: Minimum allowed value (default 0.0)
            max_val: Maximum allowed value (default 100.0)

        Raises:
            AssertionError: If any value is out of range
        """
        if "placement_uncertainty" not in df.columns:
            raise AssertionError("DataFrame missing 'placement_uncertainty' column")

        col = df["placement_uncertainty"]
        actual_min = col.min()
        actual_max = col.max()

        assert actual_min >= min_val, (
            f"placement_uncertainty has values below {min_val}: min={actual_min}"
        )
        assert actual_max <= max_val, (
            f"placement_uncertainty has values above {max_val}: max={actual_max}"
        )

    @staticmethod
    def assert_is_novel_consistency(df: pl.DataFrame) -> None:
        """
        Assert that is_novel boolean is consistent with taxonomic_call.

        is_novel should be True when taxonomic_call is "Novel Species" or
        "Novel Genus", and False otherwise.

        Args:
            df: Classification DataFrame with 'is_novel' and 'taxonomic_call'

        Raises:
            AssertionError: If is_novel is inconsistent with taxonomic_call
        """
        if "is_novel" not in df.columns:
            raise AssertionError("DataFrame missing 'is_novel' column")
        if "taxonomic_call" not in df.columns:
            raise AssertionError("DataFrame missing 'taxonomic_call' column")

        novel_calls = {"Novel Species", "Novel Genus"}

        # Check each row for consistency
        inconsistent = df.filter(
            (pl.col("taxonomic_call").is_in(novel_calls) & ~pl.col("is_novel")) |
            (~pl.col("taxonomic_call").is_in(novel_calls) & pl.col("is_novel"))
        )

        assert len(inconsistent) == 0, (
            f"Found {len(inconsistent)} rows with inconsistent is_novel values. "
            f"is_novel should be True for Novel Species/Genus, False otherwise."
        )

    @staticmethod
    def assert_best_match_genome_not_empty(df: pl.DataFrame) -> None:
        """
        Assert that best_match_genome is never empty or null.

        Args:
            df: Classification DataFrame with 'best_match_genome' column

        Raises:
            AssertionError: If any genome is empty/null
        """
        if "best_match_genome" not in df.columns:
            raise AssertionError("DataFrame missing 'best_match_genome' column")

        empty_count = df.filter(
            pl.col("best_match_genome").is_null() |
            (pl.col("best_match_genome") == "") |
            (pl.col("best_match_genome") == "unknown")
        ).height

        assert empty_count == 0, (
            f"Found {empty_count} rows with empty/null/unknown best_match_genome"
        )

    @classmethod
    def assert_all_validations(
        cls,
        path: Path,
        min_rows: int = 1,
    ) -> pl.DataFrame:
        """
        Run all classification validations on an output file.

        This is the primary entry point for comprehensive validation.

        Args:
            path: Path to classification output file
            min_rows: Minimum expected number of rows

        Returns:
            The validated DataFrame

        Raises:
            AssertionError: If any validation fails
        """
        df = cls.assert_valid_output_file(path, min_rows=min_rows)
        cls.assert_required_columns(df)
        cls.assert_valid_taxonomic_calls(df)
        cls.assert_novelty_index_range(df)
        cls.assert_placement_uncertainty_range(df)
        cls.assert_is_novel_consistency(df)
        cls.assert_best_match_genome_not_empty(df)

        return df


# =============================================================================
# Summary Statistics Assertions
# =============================================================================


class SummaryAssertions:
    """
    Assertions for validating summary JSON files.

    Provides methods to validate:
    - JSON structure and required fields
    - Count consistency (sums to total)
    - Percentage validity
    - Correlation with source DataFrame
    """

    @staticmethod
    def assert_valid_summary_file(path: Path) -> dict[str, Any]:
        """
        Assert that a summary JSON file is valid and readable.

        Args:
            path: Path to summary JSON file

        Returns:
            The parsed JSON as a dictionary

        Raises:
            AssertionError: If file is invalid or unreadable
        """
        assert path.exists(), f"Summary file does not exist: {path}"
        assert path.suffix == ".json", (
            f"Expected .json file, got: {path.suffix}"
        )

        try:
            with open(path) as f:
                data = json.load(f)
        except json.JSONDecodeError as e:
            raise AssertionError(f"Invalid JSON in summary file: {e}")

        assert isinstance(data, dict), (
            f"Summary should be a JSON object, got: {type(data).__name__}"
        )

        return data

    @staticmethod
    def assert_required_fields(summary: dict[str, Any]) -> None:
        """
        Assert that all required summary fields are present.

        Args:
            summary: Summary dictionary

        Raises:
            AssertionError: If any required field is missing
        """
        missing = REQUIRED_SUMMARY_FIELDS - set(summary.keys())

        assert not missing, (
            f"Missing required summary fields: {sorted(missing)}"
        )

    @staticmethod
    def assert_counts_sum_to_total(summary: dict[str, Any]) -> None:
        """
        Assert that category counts sum to total_reads.

        Args:
            summary: Summary dictionary with count fields

        Raises:
            AssertionError: If counts don't sum to total
        """
        total = summary.get("total_reads", 0)
        known = summary.get("known_species", 0)
        novel_sp = summary.get("novel_species", 0)
        novel_gen = summary.get("novel_genus", 0)
        conserved = summary.get("conserved_regions", 0)
        ambiguous = summary.get("ambiguous", 0)
        unclassified = summary.get("unclassified", 0)

        sum_of_counts = known + novel_sp + novel_gen + conserved + ambiguous + unclassified

        assert sum_of_counts == total, (
            f"Category counts ({sum_of_counts}) don't sum to total_reads ({total}). "
            f"Breakdown: known={known}, novel_sp={novel_sp}, "
            f"novel_gen={novel_gen}, conserved={conserved}, "
            f"ambiguous={ambiguous}, unclassified={unclassified}"
        )

    @staticmethod
    def assert_percentages_valid(summary: dict[str, Any]) -> None:
        """
        Assert that percentage fields are in valid range [0, 100].

        Args:
            summary: Summary dictionary

        Raises:
            AssertionError: If any percentage is invalid
        """
        pct_fields = [
            "known_species_pct",
            "novel_species_pct",
            "novel_genus_pct",
            "novel_diversity_pct",
        ]

        for field in pct_fields:
            if field in summary:
                val = summary[field]
                assert 0.0 <= val <= 100.0, (
                    f"{field} out of range [0, 100]: {val}"
                )

    @staticmethod
    def assert_means_valid(summary: dict[str, Any]) -> None:
        """
        Assert that mean values are in valid ranges.

        Args:
            summary: Summary dictionary

        Raises:
            AssertionError: If any mean is invalid
        """
        if "mean_novelty_index" in summary:
            val = summary["mean_novelty_index"]
            if val is not None:
                assert 0.0 <= val <= 100.0, (
                    f"mean_novelty_index out of range [0, 100]: {val}"
                )

        if "mean_placement_uncertainty" in summary:
            val = summary["mean_placement_uncertainty"]
            if val is not None:
                assert 0.0 <= val <= 100.0, (
                    f"mean_placement_uncertainty out of range [0, 100]: {val}"
                )

    @staticmethod
    def assert_matches_classification(
        summary: dict[str, Any],
        df: pl.DataFrame,
    ) -> None:
        """
        Assert that summary statistics match the classification DataFrame.

        Args:
            summary: Summary dictionary
            df: Classification DataFrame that produced the summary

        Raises:
            AssertionError: If summary doesn't match DataFrame
        """
        # Check total
        assert summary["total_reads"] == len(df), (
            f"total_reads mismatch: summary={summary['total_reads']}, df={len(df)}"
        )

        # Check category counts
        call_counts = df.group_by("taxonomic_call").len().to_dict(as_series=False)
        count_dict = dict(zip(
            call_counts["taxonomic_call"],
            call_counts["len"],
            strict=False
        ))

        expected_known = count_dict.get("Known Species", 0)
        expected_novel_sp = count_dict.get("Novel Species", 0)
        expected_novel_gen = count_dict.get("Novel Genus", 0)
        expected_conserved = count_dict.get("Conserved Region", 0)

        assert summary["known_species"] == expected_known, (
            f"known_species mismatch: summary={summary['known_species']}, "
            f"expected={expected_known}"
        )
        assert summary["novel_species"] == expected_novel_sp, (
            f"novel_species mismatch: summary={summary['novel_species']}, "
            f"expected={expected_novel_sp}"
        )
        assert summary["novel_genus"] == expected_novel_gen, (
            f"novel_genus mismatch: summary={summary['novel_genus']}, "
            f"expected={expected_novel_gen}"
        )
        assert summary["conserved_regions"] == expected_conserved, (
            f"conserved_regions mismatch: summary={summary['conserved_regions']}, "
            f"expected={expected_conserved}"
        )

    @classmethod
    def assert_all_validations(
        cls,
        path: Path,
        df: pl.DataFrame | None = None,
    ) -> dict[str, Any]:
        """
        Run all summary validations on a summary file.

        Args:
            path: Path to summary JSON file
            df: Optional classification DataFrame to check consistency

        Returns:
            The validated summary dictionary

        Raises:
            AssertionError: If any validation fails
        """
        summary = cls.assert_valid_summary_file(path)
        cls.assert_required_fields(summary)
        cls.assert_counts_sum_to_total(summary)
        cls.assert_percentages_valid(summary)
        cls.assert_means_valid(summary)

        if df is not None:
            cls.assert_matches_classification(summary, df)

        return summary


# =============================================================================
# CLI Execution Assertions
# =============================================================================


class CLIAssertions:
    """
    Assertions for validating CLI execution results.

    Works with typer.testing.CliRunner results to validate:
    - Exit codes
    - Output content
    - Error messages
    """

    @staticmethod
    def assert_success(result: Result, message: str = "") -> None:
        """
        Assert that CLI command succeeded (exit code 0).

        Args:
            result: CliRunner invocation result
            message: Optional context message for failure

        Raises:
            AssertionError: If exit code is not 0
        """
        if result.exit_code != 0:
            context = f" - {message}" if message else ""
            raise AssertionError(
                f"Expected exit code 0, got {result.exit_code}{context}\n"
                f"Output:\n{result.stdout}"
            )

    @staticmethod
    def assert_failure(
        result: Result,
        expected_code: int = 1,
        message: str = "",
    ) -> None:
        """
        Assert that CLI command failed with expected exit code.

        Args:
            result: CliRunner invocation result
            expected_code: Expected non-zero exit code
            message: Optional context message

        Raises:
            AssertionError: If exit code doesn't match expected
        """
        if result.exit_code != expected_code:
            context = f" - {message}" if message else ""
            raise AssertionError(
                f"Expected exit code {expected_code}, got {result.exit_code}{context}\n"
                f"Output:\n{result.stdout}"
            )

    @staticmethod
    def assert_output_contains(
        result: Result,
        text: str,
        case_sensitive: bool = True,
    ) -> None:
        """
        Assert that CLI output contains expected text.

        Args:
            result: CliRunner invocation result
            text: Text that should appear in output
            case_sensitive: Whether comparison is case-sensitive

        Raises:
            AssertionError: If text not found in output
        """
        output = result.stdout
        search_text = text
        search_output = output

        if not case_sensitive:
            search_text = text.lower()
            search_output = output.lower()

        assert search_text in search_output, (
            f"Expected output to contain: '{text}'\n"
            f"Actual output:\n{output}"
        )

    @staticmethod
    def assert_output_not_contains(
        result: Result,
        text: str,
        case_sensitive: bool = True,
    ) -> None:
        """
        Assert that CLI output does NOT contain specified text.

        Args:
            result: CliRunner invocation result
            text: Text that should NOT appear in output
            case_sensitive: Whether comparison is case-sensitive

        Raises:
            AssertionError: If text found in output
        """
        output = result.stdout
        search_text = text
        search_output = output

        if not case_sensitive:
            search_text = text.lower()
            search_output = output.lower()

        assert search_text not in search_output, (
            f"Expected output to NOT contain: '{text}'\n"
            f"But it was found in output:\n{output}"
        )

    @staticmethod
    def assert_output_matches(
        result: Result,
        pattern: str,
        flags: int = 0,
    ) -> re.Match[str]:
        """
        Assert that CLI output matches a regex pattern.

        Args:
            result: CliRunner invocation result
            pattern: Regex pattern to match
            flags: Regex flags (e.g., re.IGNORECASE)

        Returns:
            The match object for extracting groups

        Raises:
            AssertionError: If pattern doesn't match
        """
        match = re.search(pattern, result.stdout, flags)

        assert match is not None, (
            f"Expected output to match pattern: '{pattern}'\n"
            f"Actual output:\n{result.stdout}"
        )

        return match

    @staticmethod
    def assert_file_created(
        result: Result,
        path: Path,
        min_size: int = 0,
    ) -> None:
        """
        Assert that CLI command created a file.

        Args:
            result: CliRunner invocation result (for context on failure)
            path: Path that should have been created
            min_size: Minimum file size in bytes

        Raises:
            AssertionError: If file doesn't exist or is too small
        """
        assert path.exists(), (
            f"Expected file to be created: {path}\n"
            f"CLI output:\n{result.stdout}"
        )

        if min_size > 0:
            actual_size = path.stat().st_size
            assert actual_size >= min_size, (
                f"File {path} is too small: {actual_size} bytes "
                f"(expected >= {min_size})"
            )

    @staticmethod
    def assert_classified_reads_count(
        result: Result,
        expected_count: int,
        tolerance: int = 0,
    ) -> None:
        """
        Assert that CLI output reports expected number of classified reads.

        Parses the "Classified N reads" message from output.

        Args:
            result: CliRunner invocation result
            expected_count: Expected number of reads
            tolerance: Allowed deviation from expected

        Raises:
            AssertionError: If count doesn't match
        """
        # Match patterns like "Classified 123 reads" or "Classified 1,234 reads"
        match = re.search(r"Classified\s+([\d,]+)\s+reads", result.stdout)

        assert match is not None, (
            f"Could not find 'Classified N reads' in output:\n{result.stdout}"
        )

        actual_count = int(match.group(1).replace(",", ""))

        if tolerance == 0:
            assert actual_count == expected_count, (
                f"Expected {expected_count} classified reads, got {actual_count}"
            )
        else:
            assert abs(actual_count - expected_count) <= tolerance, (
                f"Expected ~{expected_count} classified reads (tolerance={tolerance}), "
                f"got {actual_count}"
            )


# =============================================================================
# Convenience Functions
# =============================================================================


def validate_classification_output(
    path: Path,
    min_rows: int = 1,
) -> pl.DataFrame:
    """
    Convenience function to validate a classification output file.

    Args:
        path: Path to classification file
        min_rows: Minimum expected rows

    Returns:
        Validated DataFrame
    """
    return ClassificationAssertions.assert_all_validations(path, min_rows)


def validate_summary_output(
    path: Path,
    classification_df: pl.DataFrame | None = None,
) -> dict[str, Any]:
    """
    Convenience function to validate a summary output file.

    Args:
        path: Path to summary JSON file
        classification_df: Optional DataFrame for consistency check

    Returns:
        Validated summary dictionary
    """
    return SummaryAssertions.assert_all_validations(path, classification_df)
