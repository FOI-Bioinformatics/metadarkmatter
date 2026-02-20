"""Unit tests for custom exceptions module."""

import pytest

from metadarkmatter.core.exceptions import (
    ANIMatrixError,
    ANIMatrixNotSquareError,
    ANIMatrixRowColumnMismatchError,
    ANIMatrixValueError,
    BlastFileError,
    ConfigurationError,
    EmptyBlastFileError,
    GenomeCoverageWarning,
    InvalidThresholdError,
    MalformedBlastFileError,
    MetadarkmatterError,
    ProcessingModeError,
)


class TestMetadarkmatterError:
    """Tests for base exception class."""

    def test_basic_message(self):
        """Should create exception with just message."""
        error = MetadarkmatterError("Test error")
        assert error.message == "Test error"
        assert error.suggestion is None
        assert str(error) == "Test error"

    def test_message_with_suggestion(self):
        """Should include suggestion in full message."""
        error = MetadarkmatterError("Test error", suggestion="Try this fix")
        assert error.message == "Test error"
        assert error.suggestion == "Try this fix"
        assert "Suggestion: Try this fix" in str(error)


class TestANIMatrixErrors:
    """Tests for ANI matrix error classes."""

    def test_not_square_error(self):
        """Should format not-square error message."""
        error = ANIMatrixNotSquareError(rows=10, cols=8)
        assert "10 rows x 8 columns" in str(error)
        assert "Regenerate" in error.suggestion
        assert error.rows == 10
        assert error.cols == 8

    def test_row_column_mismatch_error(self):
        """Should format row/column mismatch error."""
        error = ANIMatrixRowColumnMismatchError(
            missing_in_rows={"G1", "G2"},
            missing_in_cols={"G3"},
        )
        assert "row and column genomes don't match" in str(error)
        assert "symmetric" in error.suggestion

    def test_value_error(self):
        """Should format value error with examples."""
        error = ANIMatrixValueError([
            ("G1", "G2", 105.0),
            ("G3", "G4", -5.0),
        ])
        assert "invalid values" in str(error).lower()
        assert "0-100" in str(error)
        assert "105" in str(error) or "-5" in str(error)


class TestBlastFileErrors:
    """Tests for BLAST file error classes."""

    def test_empty_blast_file_error(self):
        """Should format empty file error with suggestions."""
        error = EmptyBlastFileError("/path/to/file.blast.tsv")
        assert "empty" in str(error).lower()
        assert "file.blast.tsv" in str(error)
        assert "Query sequences" in error.suggestion

    def test_malformed_blast_file_error(self):
        """Should format malformed file error with details."""
        error = MalformedBlastFileError(
            path="/path/to/file.blast.tsv",
            expected_cols=12,
            actual_cols=5,
            line_num=42,
        )
        assert "line 42" in str(error)
        assert "12 columns" in str(error)
        assert "5" in str(error)
        assert "-outfmt 6" in error.suggestion


class TestGenomeCoverageWarning:
    """Tests for genome coverage warning."""

    def test_warning_message(self):
        """Should format coverage warning with stats."""
        warning = GenomeCoverageWarning(
            matched=10,
            total=50,
            coverage_pct=20.0,
            missing_examples=["G1", "G2", "G3"],
        )
        assert "10/50" in str(warning)
        assert "20.0%" in str(warning)
        assert "G1" in warning.suggestion

    def test_many_missing_genomes(self):
        """Should truncate long list of missing genomes."""
        missing = [f"G{i}" for i in range(20)]
        warning = GenomeCoverageWarning(
            matched=5,
            total=25,
            coverage_pct=20.0,
            missing_examples=missing,
        )
        assert "and 15 more" in warning.suggestion


class TestConfigurationErrors:
    """Tests for configuration error classes."""

    def test_invalid_threshold_error(self):
        """Should format threshold error with valid range."""
        error = InvalidThresholdError(
            param_name="bitscore_threshold",
            value=150.0,
            min_val=0.0,
            max_val=100.0,
        )
        assert "bitscore_threshold" in str(error)
        assert "150" in str(error)
        assert "[0.0, 100.0]" in str(error)

    def test_processing_mode_error(self):
        """Should format mode error with available options."""
        error = ProcessingModeError(["--streaming", "--other"])
        assert "streaming" in str(error)
        assert "other" in str(error)
        assert "Choose only one" in error.suggestion


class TestExceptionHierarchy:
    """Tests for exception class hierarchy."""

    def test_ani_error_is_metadarkmatter_error(self):
        """ANI errors should inherit from base class."""
        error = ANIMatrixNotSquareError(10, 8)
        assert isinstance(error, MetadarkmatterError)
        assert isinstance(error, ANIMatrixError)

    def test_blast_error_is_metadarkmatter_error(self):
        """BLAST errors should inherit from base class."""
        error = EmptyBlastFileError("/path")
        assert isinstance(error, MetadarkmatterError)
        assert isinstance(error, BlastFileError)

    def test_config_error_is_metadarkmatter_error(self):
        """Config errors should inherit from base class."""
        error = InvalidThresholdError("test", 0, 0, 100)
        assert isinstance(error, MetadarkmatterError)
        assert isinstance(error, ConfigurationError)
