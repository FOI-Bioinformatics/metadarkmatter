"""
Unit tests for CLI utility functions.

Tests for spinner_progress, extract_sample_name, and QuietConsole.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metadarkmatter.cli.utils import extract_sample_name


class TestExtractSampleName:
    """Tests for extract_sample_name function."""

    def test_simple_name(self) -> None:
        """Test extraction from simple filename."""
        path = Path("sample_001.csv")
        assert extract_sample_name(path) == "sample_001"

    def test_classifications_suffix(self) -> None:
        """Test removal of _classifications suffix."""
        path = Path("sample_001_classifications.csv")
        assert extract_sample_name(path) == "sample_001"

    def test_classifications_prefix(self) -> None:
        """Test removal of classifications_ prefix."""
        path = Path("classifications_sample_001.csv")
        assert extract_sample_name(path) == "sample_001"

    def test_kraken_suffix(self) -> None:
        """Test removal of .kraken suffix."""
        path = Path("sample_001.kraken")
        assert extract_sample_name(path) == "sample_001"

    def test_blast_tsv_suffix(self) -> None:
        """Test removal of .blast.tsv suffix."""
        path = Path("sample_001.blast.tsv")
        assert extract_sample_name(path) == "sample_001"

    def test_gzipped_file(self) -> None:
        """Test handling of .gz extension."""
        path = Path("sample_001.blast.tsv.gz")
        assert extract_sample_name(path) == "sample_001"

    def test_custom_suffixes(self) -> None:
        """Test removal of custom suffixes."""
        path = Path("sample_001_custom.csv")
        assert extract_sample_name(path, ("_custom",)) == "sample_001"

    def test_multiple_suffixes(self) -> None:
        """Test removal of multiple suffixes."""
        path = Path("sample_001_classifications.blast.tsv")
        # After removing _classifications, .blast, .tsv
        assert extract_sample_name(path) == "sample_001"

    def test_nested_path(self) -> None:
        """Test extraction from nested path."""
        path = Path("/data/results/sample_001_classifications.csv")
        assert extract_sample_name(path) == "sample_001"

    def test_no_suffix_to_remove(self) -> None:
        """Test when no suffixes need removal."""
        path = Path("my_sample.csv")
        assert extract_sample_name(path) == "my_sample"
