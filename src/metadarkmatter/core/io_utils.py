"""
I/O utilities for DataFrame serialization.

Provides consistent handling of output formats (CSV/Parquet) across the codebase.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import polars as pl

OutputFormat = Literal["csv", "parquet"]


def write_dataframe(
    df: pl.DataFrame,
    path: Path,
    output_format: OutputFormat = "csv",
) -> None:
    """
    Write DataFrame to file in specified format.

    For Parquet output, uses zstd compression for optimal size/speed tradeoff
    (approximately 10x smaller than CSV with fast decompression).

    Args:
        df: Polars DataFrame to write.
        path: Output file path.
        output_format: Output format - 'csv' or 'parquet'.

    Example:
        >>> df = pl.DataFrame({"a": [1, 2, 3]})
        >>> write_dataframe(df, Path("output.parquet"), "parquet")
    """
    if output_format == "parquet":
        df.write_parquet(path, compression="zstd")
    else:
        df.write_csv(path)


def write_dataframe_append(
    df: pl.DataFrame,
    path: Path,
    output_format: OutputFormat = "csv",
    is_first: bool = True,
) -> None:
    """
    Write or append DataFrame to file in specified format.

    For streaming/chunked processing where data is written incrementally.
    CSV appends use file append mode; Parquet requires read-concat-write.

    Args:
        df: Polars DataFrame to write or append.
        path: Output file path.
        output_format: Output format - 'csv' or 'parquet'.
        is_first: If True, create new file; if False, append to existing.

    Note:
        Parquet append is not natively supported, so this reads the existing
        file, concatenates, and rewrites. For very large files, consider
        partitioned output instead.
    """
    if output_format == "parquet":
        if is_first:
            df.write_parquet(path, compression="zstd")
        else:
            existing = pl.read_parquet(path)
            combined = pl.concat([existing, df])
            combined.write_parquet(path, compression="zstd")
    else:
        if is_first:
            df.write_csv(path)
        else:
            # Append to CSV without header
            with path.open("a") as f:
                f.write(df.write_csv(include_header=False))


def read_dataframe(path: Path) -> pl.DataFrame:
    """
    Read DataFrame from file, auto-detecting format from extension.

    Supports: .csv, .tsv, .parquet, .csv.gz, .tsv.gz

    Args:
        path: Input file path.

    Returns:
        Polars DataFrame.

    Raises:
        ValueError: If file extension is not recognized.
    """
    suffix = path.suffix.lower()
    name = path.name.lower()

    if suffix == ".parquet":
        return pl.read_parquet(path)
    if suffix == ".csv" or name.endswith(".csv.gz"):
        return pl.read_csv(path)
    if suffix == ".tsv" or name.endswith(".tsv.gz"):
        return pl.read_csv(path, separator="\t")
    msg = f"Unrecognized file format: {path}"
    raise ValueError(msg)
