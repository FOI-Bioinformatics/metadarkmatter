"""
I/O utilities for DataFrame serialization.

Provides consistent handling of output formats (CSV/Parquet) across the codebase.
"""

from __future__ import annotations

from pathlib import Path
from types import TracebackType
from typing import TYPE_CHECKING, Literal

import polars as pl

if TYPE_CHECKING:
    import pyarrow.parquet as pq

OutputFormat = Literal["csv", "parquet"]


class StreamingWriter:
    """Append DataFrame partitions to a single file in O(total rows).

    Replaces the naive parquet read-concat-rewrite append (which is O(n^2)
    across K partitions) with an incremental writer:

    - parquet: a ``pyarrow.parquet.ParquetWriter`` opened on the first
      partition and fed one row group per partition. Each partition is written
      exactly once; memory stays bounded to a single partition.
    - csv: native append mode (header on the first partition only).

    Schemas must match across partitions (the classifier emits a stable schema
    for a given config). Use as a context manager so the parquet file footer is
    always finalised::

        with StreamingWriter(path, "parquet") as w:
            for part in partitions:
                w.write(part)
    """

    def __init__(self, path: Path, output_format: OutputFormat = "parquet") -> None:
        self.path = path
        self.output_format = output_format
        self._pq_writer: pq.ParquetWriter | None = None
        self._csv_started = False

    def write(self, df: pl.DataFrame) -> None:
        """Append one partition."""
        if self.output_format == "parquet":
            import pyarrow.parquet as pq

            table = df.to_arrow()
            if self._pq_writer is None:
                self._pq_writer = pq.ParquetWriter(
                    self.path, table.schema, compression="zstd"
                )
            self._pq_writer.write_table(table)
        elif not self._csv_started:
            df.write_csv(self.path)
            self._csv_started = True
        else:
            with self.path.open("a") as f:
                f.write(df.write_csv(include_header=False))

    def close(self) -> None:
        if self._pq_writer is not None:
            self._pq_writer.close()
            self._pq_writer = None

    def __enter__(self) -> StreamingWriter:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        self.close()


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
