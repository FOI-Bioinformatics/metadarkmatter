"""
Unit tests for I/O utility functions.

Tests for write_dataframe, write_dataframe_append, and read_dataframe functions.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from metadarkmatter.core.io_utils import (
    read_dataframe,
    write_dataframe,
    write_dataframe_append,
)


@pytest.fixture
def sample_df() -> pl.DataFrame:
    """Create a sample DataFrame for testing."""
    return pl.DataFrame({
        "read_id": ["read1", "read2", "read3"],
        "genome": ["GCF_001", "GCF_002", "GCF_001"],
        "identity": [95.5, 88.2, 92.1],
    })


class TestWriteDataframe:
    """Tests for write_dataframe function."""

    def test_write_csv(self, sample_df: pl.DataFrame, tmp_path: Path) -> None:
        """Should write DataFrame to CSV format."""
        output = tmp_path / "output.csv"
        write_dataframe(sample_df, output, "csv")

        assert output.exists()
        loaded = pl.read_csv(output)
        assert loaded.shape == sample_df.shape
        assert loaded.columns == sample_df.columns

    def test_write_parquet(self, sample_df: pl.DataFrame, tmp_path: Path) -> None:
        """Should write DataFrame to Parquet format with zstd compression."""
        output = tmp_path / "output.parquet"
        write_dataframe(sample_df, output, "parquet")

        assert output.exists()
        loaded = pl.read_parquet(output)
        assert loaded.shape == sample_df.shape
        assert loaded.columns == sample_df.columns

    def test_parquet_smaller_than_csv(
        self, sample_df: pl.DataFrame, tmp_path: Path
    ) -> None:
        """Parquet with zstd compression should typically be smaller than CSV."""
        # Create larger DataFrame for meaningful size comparison
        large_df = pl.concat([sample_df] * 1000)

        csv_path = tmp_path / "output.csv"
        parquet_path = tmp_path / "output.parquet"

        write_dataframe(large_df, csv_path, "csv")
        write_dataframe(large_df, parquet_path, "parquet")

        # Parquet with compression should be smaller
        assert parquet_path.stat().st_size < csv_path.stat().st_size


class TestWriteDataframeAppend:
    """Tests for write_dataframe_append function."""

    def test_first_write_creates_file(
        self, sample_df: pl.DataFrame, tmp_path: Path
    ) -> None:
        """First write should create new file."""
        output = tmp_path / "output.csv"
        write_dataframe_append(sample_df, output, "csv", is_first=True)

        assert output.exists()
        loaded = pl.read_csv(output)
        assert loaded.shape == sample_df.shape

    def test_append_csv(self, sample_df: pl.DataFrame, tmp_path: Path) -> None:
        """Should append to existing CSV file."""
        output = tmp_path / "output.csv"

        # First write
        write_dataframe_append(sample_df, output, "csv", is_first=True)

        # Append
        write_dataframe_append(sample_df, output, "csv", is_first=False)

        loaded = pl.read_csv(output)
        assert loaded.height == sample_df.height * 2

    def test_append_parquet(self, sample_df: pl.DataFrame, tmp_path: Path) -> None:
        """Should append to existing Parquet file."""
        output = tmp_path / "output.parquet"

        # First write
        write_dataframe_append(sample_df, output, "parquet", is_first=True)

        # Append
        write_dataframe_append(sample_df, output, "parquet", is_first=False)

        loaded = pl.read_parquet(output)
        assert loaded.height == sample_df.height * 2

    def test_multiple_appends(self, sample_df: pl.DataFrame, tmp_path: Path) -> None:
        """Should handle multiple appends correctly."""
        output = tmp_path / "output.csv"

        write_dataframe_append(sample_df, output, "csv", is_first=True)
        write_dataframe_append(sample_df, output, "csv", is_first=False)
        write_dataframe_append(sample_df, output, "csv", is_first=False)

        loaded = pl.read_csv(output)
        assert loaded.height == sample_df.height * 3


class TestReadDataframe:
    """Tests for read_dataframe function."""

    def test_read_csv(self, sample_df: pl.DataFrame, tmp_path: Path) -> None:
        """Should read CSV file."""
        path = tmp_path / "input.csv"
        sample_df.write_csv(path)

        loaded = read_dataframe(path)
        assert loaded.shape == sample_df.shape

    def test_read_tsv(self, sample_df: pl.DataFrame, tmp_path: Path) -> None:
        """Should read TSV file."""
        path = tmp_path / "input.tsv"
        sample_df.write_csv(path, separator="\t")

        loaded = read_dataframe(path)
        assert loaded.shape == sample_df.shape

    def test_read_parquet(self, sample_df: pl.DataFrame, tmp_path: Path) -> None:
        """Should read Parquet file."""
        path = tmp_path / "input.parquet"
        sample_df.write_parquet(path)

        loaded = read_dataframe(path)
        assert loaded.shape == sample_df.shape

    def test_unrecognized_format_raises(self, tmp_path: Path) -> None:
        """Should raise ValueError for unrecognized format."""
        path = tmp_path / "input.xyz"
        path.touch()

        with pytest.raises(ValueError, match="Unrecognized file format"):
            read_dataframe(path)

    def test_auto_detect_format(self, sample_df: pl.DataFrame, tmp_path: Path) -> None:
        """Should auto-detect format from extension."""
        # Test CSV
        csv_path = tmp_path / "data.csv"
        sample_df.write_csv(csv_path)
        assert read_dataframe(csv_path).shape == sample_df.shape

        # Test Parquet
        parquet_path = tmp_path / "data.parquet"
        sample_df.write_parquet(parquet_path)
        assert read_dataframe(parquet_path).shape == sample_df.shape


class TestRoundtrip:
    """Test write and read roundtrip."""

    def test_csv_roundtrip_preserves_data(
        self, sample_df: pl.DataFrame, tmp_path: Path
    ) -> None:
        """CSV roundtrip should preserve data."""
        path = tmp_path / "roundtrip.csv"
        write_dataframe(sample_df, path, "csv")
        loaded = read_dataframe(path)

        # Compare values (types may differ slightly for CSV)
        assert loaded["read_id"].to_list() == sample_df["read_id"].to_list()
        assert loaded["genome"].to_list() == sample_df["genome"].to_list()

    def test_parquet_roundtrip_preserves_data(
        self, sample_df: pl.DataFrame, tmp_path: Path
    ) -> None:
        """Parquet roundtrip should preserve data exactly."""
        path = tmp_path / "roundtrip.parquet"
        write_dataframe(sample_df, path, "parquet")
        loaded = read_dataframe(path)

        # Parquet preserves types exactly
        assert loaded.equals(sample_df)
