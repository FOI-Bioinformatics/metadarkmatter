"""
Custom exceptions with actionable guidance.

Provides specific error types for common failure scenarios,
each with helpful suggestions for resolution.
"""

from __future__ import annotations


class MetadarkmatterError(Exception):
    """Base exception for metadarkmatter errors."""

    def __init__(self, message: str, suggestion: str | None = None):
        self.message = message
        self.suggestion = suggestion
        super().__init__(self.full_message)

    @property
    def full_message(self) -> str:
        if self.suggestion:
            return f"{self.message}\n\nSuggestion: {self.suggestion}"
        return self.message


class ANIMatrixError(MetadarkmatterError):
    """Base class for ANI matrix errors."""



class ANIMatrixNotSquareError(ANIMatrixError):
    """Raised when ANI matrix is not square."""

    def __init__(self, rows: int, cols: int):
        super().__init__(
            message=f"ANI matrix is not square: {rows} rows x {cols} columns",
            suggestion=(
                "Regenerate the ANI matrix using FastANI or pyani. "
                "Ensure all genome pairs are included and the output is symmetric."
            ),
        )
        self.rows = rows
        self.cols = cols


class ANIMatrixRowColumnMismatchError(ANIMatrixError):
    """Raised when row and column genome names don't match."""

    def __init__(self, missing_in_rows: set[str], missing_in_cols: set[str]):
        message = "ANI matrix row and column genomes don't match"
        details = []
        if missing_in_rows:
            details.append(f"Missing in rows: {', '.join(sorted(missing_in_rows)[:5])}")
        if missing_in_cols:
            details.append(f"Missing in columns: {', '.join(sorted(missing_in_cols)[:5])}")

        super().__init__(
            message=f"{message}. {'; '.join(details)}",
            suggestion=(
                "Check that the ANI matrix file has matching row and column headers. "
                "The matrix should be symmetric with all genomes present in both dimensions."
            ),
        )


class ANIMatrixValueError(ANIMatrixError):
    """Raised when ANI values are out of valid range."""

    def __init__(self, invalid_values: list[tuple[str, str, float]]):
        examples = invalid_values[:3]
        example_str = ", ".join(
            f"{g1}-{g2}: {val:.1f}" for g1, g2, val in examples
        )
        super().__init__(
            message=f"ANI matrix contains invalid values (must be 0-100): {example_str}",
            suggestion=(
                "ANI values must be between 0 and 100. Check your ANI calculation "
                "tool output. Negative values or values >100 indicate a malformed matrix."
            ),
        )


class BlastFileError(MetadarkmatterError):
    """Base class for BLAST file errors."""



class EmptyBlastFileError(BlastFileError):
    """Raised when BLAST file is empty or has no valid alignments."""

    def __init__(self, path: str):
        super().__init__(
            message=f"BLAST file is empty or contains no valid alignments: {path}",
            suggestion=(
                "Check that your BLAST search completed successfully. "
                "Verify that:\n"
                "  - Query sequences are in FASTA format\n"
                "  - Database was built correctly with makeblastdb\n"
                "  - BLAST parameters (-evalue, -perc_identity) aren't too strict\n"
                "  - Output format is -outfmt 6 (tabular)"
            ),
        )


class MalformedBlastFileError(BlastFileError):
    """Raised when BLAST file has incorrect format."""

    def __init__(self, path: str, expected_cols: int, actual_cols: int, line_num: int):
        super().__init__(
            message=(
                f"Malformed BLAST file '{path}' at line {line_num}: "
                f"expected {expected_cols} columns, got {actual_cols}"
            ),
            suggestion=(
                "BLAST file must be in -outfmt 6 format (12 tab-separated columns):\n"
                "  qseqid sseqid pident length mismatch gapopen qstart qend "
                "sstart send evalue bitscore\n\n"
                "Check that your BLAST command used -outfmt 6 and the file "
                "wasn't corrupted during transfer."
            ),
        )


class GenomeCoverageWarning(MetadarkmatterError):
    """Warning when genome coverage is low but not critical."""

    def __init__(
        self,
        matched: int,
        total: int,
        coverage_pct: float,
        missing_examples: list[str],
    ):
        missing_str = ", ".join(missing_examples[:5])
        if len(missing_examples) > 5:
            missing_str += f"... and {len(missing_examples) - 5} more"

        super().__init__(
            message=(
                f"Low genome coverage: only {matched}/{total} "
                f"({coverage_pct:.1f}%) of BLAST genomes are in the ANI matrix"
            ),
            suggestion=(
                f"Missing genomes: {missing_str}\n\n"
                f"To improve classification accuracy:\n"
                f"  1. Add missing genomes to your reference database\n"
                f"  2. Regenerate the ANI matrix to include all genomes\n"
                f"  3. Or use a BLAST database that matches your ANI matrix"
            ),
        )


class ConfigurationError(MetadarkmatterError):
    """Raised when configuration is invalid."""



class InvalidThresholdError(ConfigurationError):
    """Raised when a threshold parameter is out of valid range."""

    def __init__(self, param_name: str, value: float, min_val: float, max_val: float):
        super().__init__(
            message=f"{param_name} = {value} is out of valid range [{min_val}, {max_val}]",
            suggestion=f"Set {param_name} to a value between {min_val} and {max_val}.",
        )


class ProcessingModeError(ConfigurationError):
    """Raised when processing mode configuration is invalid."""

    def __init__(self, modes: list[str]):
        super().__init__(
            message=f"Multiple processing modes specified: {', '.join(modes)}",
            suggestion=(
                "Choose only one processing mode:\n"
                "  --fast      : Single-threaded optimized (~3x faster)\n"
                "  --parallel  : Polars vectorized (~16x faster, recommended)\n"
                "  --streaming : Memory-efficient for 100M+ alignments"
            ),
        )
