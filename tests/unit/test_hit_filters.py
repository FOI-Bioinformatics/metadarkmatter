"""
Unit tests for hit-level alignment filters.

Tests that evalue, percent identity, bitscore, read length, and query
coverage filters are applied correctly in VectorizedClassifier.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import polars as pl
import pytest

from metadarkmatter.core.classification.ani_matrix import ANIMatrix
from metadarkmatter.core.classification.classifiers.vectorized import VectorizedClassifier
from metadarkmatter.models.config import ScoringConfig


def _make_ani_matrix() -> ANIMatrix:
    """Build a minimal 3-genome ANI matrix for testing."""
    return ANIMatrix({
        "GCF_A": {"GCF_A": 100.0, "GCF_B": 96.0, "GCF_C": 80.0},
        "GCF_B": {"GCF_A": 96.0, "GCF_B": 100.0, "GCF_C": 81.0},
        "GCF_C": {"GCF_A": 80.0, "GCF_B": 81.0, "GCF_C": 100.0},
    })


def _write_blast_file(
    path: Path,
    rows: list[dict],
    include_qlen: bool = False,
) -> Path:
    """Write BLAST tabular rows to a TSV file.

    Each row dict must have the standard 12 BLAST columns.
    If include_qlen is True, a 13th column is written.
    """
    cols = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    ]
    if include_qlen:
        cols.append("qlen")
    df = pl.DataFrame(rows)
    df.select(cols).write_csv(path, separator="\t", include_header=False)
    return path


def _single_hit(
    *,
    read_id: str = "read_1",
    genome: str = "GCF_A",
    evalue: float = 1e-50,
    pident: float = 98.0,
    bitscore: float = 250.0,
    length: int = 300,
    qstart: int = 1,
    qend: int = 300,
    qlen: int | None = None,
) -> dict:
    """Build a single BLAST hit dict."""
    row = {
        "qseqid": read_id,
        "sseqid": f"{genome}|contig1",
        "pident": pident,
        "length": length,
        "mismatch": 1,
        "gapopen": 0,
        "qstart": qstart,
        "qend": qend,
        "sstart": 1,
        "send": length,
        "evalue": evalue,
        "bitscore": bitscore,
    }
    if qlen is not None:
        row["qlen"] = qlen
    return row


def _classify(config: ScoringConfig, rows: list[dict], include_qlen: bool = False) -> pl.DataFrame:
    """Helper: write hits to a temp file and classify."""
    ani = _make_ani_matrix()
    clf = VectorizedClassifier(ani_matrix=ani, config=config)
    with tempfile.TemporaryDirectory() as tmpdir:
        path = Path(tmpdir) / "hits.tsv"
        _write_blast_file(path, rows, include_qlen=include_qlen)
        return clf.classify_file(path)


class TestDefaultFiltersPassThrough:
    """Default config (all filters at 0) should not remove any hits."""

    def test_default_filters_pass_all_hits(self):
        config = ScoringConfig(min_alignment_length=0, min_alignment_fraction=0.0)
        rows = [
            _single_hit(read_id="r1", evalue=1e-80, pident=99.0, bitscore=500.0),
            _single_hit(read_id="r2", evalue=1e-5, pident=85.0, bitscore=80.0),
            _single_hit(read_id="r3", evalue=1e-40, pident=92.0, bitscore=300.0, genome="GCF_B"),
            _single_hit(read_id="r4", evalue=0.1, pident=70.0, bitscore=30.0, genome="GCF_C"),
        ]
        result = _classify(config, rows)
        assert len(result) == 4


class TestEvalueFilter:
    """Tests for max_evalue filter."""

    def test_evalue_filter_removes_high_evalue(self):
        config = ScoringConfig(
            max_evalue=1e-10,
            min_alignment_length=0,
            min_alignment_fraction=0.0,
        )
        rows = [_single_hit(evalue=1e-5)]
        result = _classify(config, rows)
        assert len(result) == 0

    def test_evalue_filter_keeps_low_evalue(self):
        config = ScoringConfig(
            max_evalue=1e-10,
            min_alignment_length=0,
            min_alignment_fraction=0.0,
        )
        rows = [_single_hit(evalue=1e-50)]
        result = _classify(config, rows)
        assert len(result) == 1

    def test_evalue_zero_disables_filter(self):
        config = ScoringConfig(
            max_evalue=0.0,
            min_alignment_length=0,
            min_alignment_fraction=0.0,
        )
        rows = [_single_hit(evalue=10.0)]
        result = _classify(config, rows)
        assert len(result) == 1


class TestPercentIdentityFilter:
    """Tests for min_percent_identity filter."""

    def test_identity_filter_removes_low_identity(self):
        config = ScoringConfig(
            min_percent_identity=90.0,
            min_alignment_length=0,
            min_alignment_fraction=0.0,
        )
        rows = [_single_hit(pident=85.0)]
        result = _classify(config, rows)
        assert len(result) == 0

    def test_identity_filter_keeps_high_identity(self):
        config = ScoringConfig(
            min_percent_identity=90.0,
            min_alignment_length=0,
            min_alignment_fraction=0.0,
        )
        rows = [_single_hit(pident=95.0)]
        result = _classify(config, rows)
        assert len(result) == 1


class TestBitscoreFilter:
    """Tests for min_bitscore filter."""

    def test_bitscore_filter_removes_low_bitscore(self):
        config = ScoringConfig(
            min_bitscore=100.0,
            min_alignment_length=0,
            min_alignment_fraction=0.0,
        )
        rows = [_single_hit(bitscore=50.0)]
        result = _classify(config, rows)
        assert len(result) == 0

    def test_bitscore_filter_keeps_high_bitscore(self):
        config = ScoringConfig(
            min_bitscore=100.0,
            min_alignment_length=0,
            min_alignment_fraction=0.0,
        )
        rows = [_single_hit(bitscore=250.0)]
        result = _classify(config, rows)
        assert len(result) == 1


class TestReadLengthFilter:
    """Tests for min_read_length filter."""

    def test_read_length_filter_without_qlen(self):
        """When qlen is absent, qend is used as proxy for read length."""
        config = ScoringConfig(
            min_read_length=200,
            min_alignment_length=0,
            min_alignment_fraction=0.0,
        )
        # qend=150 < 200 -> filtered out
        rows = [_single_hit(qend=150, length=150)]
        result = _classify(config, rows)
        assert len(result) == 0

    def test_read_length_filter_with_qlen(self):
        """When qlen is present, it takes priority over qend."""
        config = ScoringConfig(
            min_read_length=200,
            min_alignment_length=0,
            min_alignment_fraction=0.0,
        )
        # qlen=300 >= 200 -> kept (qend is small but qlen is authoritative)
        rows = [_single_hit(qend=150, length=150, qlen=300)]
        result = _classify(config, rows, include_qlen=True)
        assert len(result) == 1

    def test_read_length_filter_passes_large_reads(self):
        config = ScoringConfig(
            min_read_length=200,
            min_alignment_length=0,
            min_alignment_fraction=0.0,
        )
        rows = [_single_hit(qend=300)]
        result = _classify(config, rows)
        assert len(result) == 1


class TestQueryCoverageFilter:
    """Tests for min_query_coverage filter."""

    def test_query_coverage_filter_removes_low_coverage(self):
        config = ScoringConfig(
            min_query_coverage=80.0,
            min_alignment_length=0,
            min_alignment_fraction=0.0,
        )
        # aligned = qend - qstart + 1 = 100 - 50 + 1 = 51
        # read_length = qend = 300 (no qlen)
        # coverage = 51/300 * 100 = 17% < 80%
        rows = [_single_hit(qstart=50, qend=300, length=251)]
        # Need read_length = qend = 300, aligned = 300-50+1 = 251
        # coverage = 251/300*100 = 83.7% >= 80% -- that would pass
        # Let me make it fail: qstart=200, qend=300, length=101
        # aligned = 300-200+1 = 101, read_length = qend = 300
        # coverage = 101/300*100 = 33.7% < 80%
        rows = [_single_hit(qstart=200, qend=300, length=101)]
        result = _classify(config, rows)
        assert len(result) == 0

    def test_query_coverage_filter_keeps_high_coverage(self):
        config = ScoringConfig(
            min_query_coverage=80.0,
            min_alignment_length=0,
            min_alignment_fraction=0.0,
        )
        # aligned = 300 - 1 + 1 = 300, read_length = qend = 300
        # coverage = 300/300*100 = 100% >= 80%
        rows = [_single_hit(qstart=1, qend=300)]
        result = _classify(config, rows)
        assert len(result) == 1

    def test_query_coverage_with_qlen(self):
        """Query coverage should use qlen when available."""
        config = ScoringConfig(
            min_query_coverage=50.0,
            min_alignment_length=0,
            min_alignment_fraction=0.0,
        )
        # aligned = 150 - 1 + 1 = 150, qlen = 300
        # coverage = 150/300*100 = 50% >= 50%
        rows = [_single_hit(qstart=1, qend=150, length=150, qlen=300)]
        result = _classify(config, rows, include_qlen=True)
        assert len(result) == 1


class TestCombinedFilters:
    """Tests for multiple filters applied together."""

    def test_combined_filters_cumulative(self):
        """Multiple active filters should be applied cumulatively."""
        config = ScoringConfig(
            max_evalue=1e-10,
            min_percent_identity=90.0,
            min_bitscore=100.0,
            min_alignment_length=0,
            min_alignment_fraction=0.0,
        )
        rows = [
            # r1: evalue=1e-80 OK, pident=99 OK, bitscore=500 OK -> kept
            _single_hit(read_id="r1", evalue=1e-80, pident=99.0, bitscore=500.0),
            # r2: evalue=1e-5 FAIL -> removed
            _single_hit(read_id="r2", evalue=1e-5, pident=85.0, bitscore=80.0),
            # r3: evalue=1e-40 OK, pident=92 OK, bitscore=300 OK -> kept
            _single_hit(read_id="r3", evalue=1e-40, pident=92.0, bitscore=300.0, genome="GCF_B"),
            # r4: evalue=0.1 FAIL -> removed
            _single_hit(read_id="r4", evalue=0.1, pident=70.0, bitscore=30.0, genome="GCF_C"),
        ]
        result = _classify(config, rows)
        assert len(result) == 2
        read_ids = set(result["read_id"].to_list())
        assert "r1" in read_ids
        assert "r3" in read_ids

    def test_all_filters_disabled_passes_everything(self):
        """Config with all filters at 0 should pass every hit."""
        config = ScoringConfig(
            max_evalue=0.0,
            min_percent_identity=0.0,
            min_bitscore=0.0,
            min_read_length=0,
            min_query_coverage=0.0,
            min_alignment_length=0,
            min_alignment_fraction=0.0,
        )
        rows = [
            _single_hit(read_id="r1", evalue=1e-80, pident=99.0, bitscore=500.0),
            _single_hit(read_id="r2", evalue=1e-5, pident=85.0, bitscore=80.0),
            _single_hit(read_id="r3", evalue=1e-40, pident=92.0, bitscore=300.0, genome="GCF_B"),
            _single_hit(read_id="r4", evalue=0.1, pident=70.0, bitscore=30.0, genome="GCF_C"),
        ]
        result = _classify(config, rows)
        assert len(result) == 4
