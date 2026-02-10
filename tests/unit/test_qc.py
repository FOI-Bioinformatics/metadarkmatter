"""
Unit tests for QC metrics module.

Tests compute_pre_qc and compute_post_qc functions
including warning generation and metric calculations.
"""

from __future__ import annotations

import polars as pl
import pytest

from metadarkmatter.core.classification.qc import (
    QCMetrics,
    compute_post_qc,
    compute_pre_qc,
)
from metadarkmatter.models.config import ScoringConfig


class TestQCMetrics:
    """Tests for QCMetrics dataclass."""

    def test_default_values(self):
        """Default QCMetrics should have zero/empty values."""
        qc = QCMetrics()
        assert qc.total_alignments == 0
        assert qc.filtered_alignments == 0
        assert qc.warnings == []
        assert qc.classification_counts == {}

    def test_to_dict(self):
        """to_dict should produce JSON-serializable dict."""
        qc = QCMetrics(
            total_alignments=1000,
            filtered_alignments=200,
            filter_rate=0.2,
            warnings=["test warning"],
        )
        d = qc.to_dict()
        assert d["total_alignments"] == 1000
        assert d["filter_rate"] == 0.2
        assert d["warnings"] == ["test warning"]

    def test_to_dict_rounds_floats(self):
        """to_dict should round float values."""
        qc = QCMetrics(filter_rate=0.123456789)
        d = qc.to_dict()
        assert d["filter_rate"] == 0.1235


class TestPreQC:
    """Tests for compute_pre_qc function."""

    @pytest.fixture
    def config(self):
        return ScoringConfig()

    def test_basic_pre_qc(self, ani_matrix, config):
        """Should compute basic alignment statistics."""
        raw_df = pl.DataFrame({
            "qseqid": ["r1", "r1", "r2", "r3"],
            "sseqid": ["g1", "g2", "g1", "g1"],
            "pident": [99.0, 95.0, 90.0, 85.0],
            "length": [150, 120, 140, 130],
            "genome_name": ["GCF_000123456.1", "GCF_000789012.1", "GCF_000123456.1", "GCF_000123456.1"],
        })
        filtered_df = raw_df  # No filtering

        qc = compute_pre_qc(raw_df, filtered_df, ani_matrix, config)

        assert qc.total_alignments == 4
        assert qc.filtered_alignments == 0
        assert qc.filter_rate == 0.0
        assert qc.genomes_in_alignment == 2
        assert qc.genomes_in_ani == 3
        assert qc.total_reads == 3
        assert qc.mean_identity > 0

    def test_filter_rate(self, ani_matrix, config):
        """Should compute correct filter rate."""
        raw_df = pl.DataFrame({
            "qseqid": ["r1", "r2", "r3", "r4"],
            "sseqid": ["g1", "g2", "g3", "g4"],
            "pident": [99.0, 95.0, 90.0, 85.0],
            "length": [150, 120, 140, 130],
            "genome_name": ["GCF_000123456.1", "GCF_000789012.1", "GCF_000123456.1", "GCF_000123456.1"],
        })
        filtered_df = raw_df.head(2)  # Only 2 of 4 remain

        qc = compute_pre_qc(raw_df, filtered_df, ani_matrix, config)

        assert qc.filtered_alignments == 2
        assert qc.filter_rate == 0.5

    def test_genome_coverage(self, ani_matrix, config):
        """Should compute genome coverage against ANI matrix."""
        raw_df = pl.DataFrame({
            "qseqid": ["r1"],
            "sseqid": ["g1"],
            "pident": [99.0],
            "length": [150],
            "genome_name": ["GCF_000123456.1"],
        })

        qc = compute_pre_qc(raw_df, raw_df, ani_matrix, config)

        # Only 1 of 3 ANI genomes seen
        assert qc.genomes_in_alignment == 1
        assert qc.genome_coverage == pytest.approx(1 / 3, rel=0.01)
        assert len(qc.missing_genomes) == 2

    def test_single_hit_fraction(self, ani_matrix, config):
        """Should compute fraction of reads with single hits."""
        raw_df = pl.DataFrame({
            "qseqid": ["r1", "r2", "r3", "r3"],
            "sseqid": ["g1", "g1", "g1", "g2"],
            "pident": [99.0, 95.0, 90.0, 88.0],
            "length": [150, 150, 150, 150],
            "genome_name": ["GCF_000123456.1", "GCF_000123456.1", "GCF_000123456.1", "GCF_000789012.1"],
        })

        qc = compute_pre_qc(raw_df, raw_df, ani_matrix, config)

        # r1 and r2 are single-hit, r3 is multi-hit
        assert qc.total_reads == 3
        assert qc.single_hit_fraction == pytest.approx(2 / 3, rel=0.01)

    def test_high_filter_rate_warning(self, ani_matrix, config):
        """Should warn when filter rate exceeds 50%."""
        raw_df = pl.DataFrame({
            "qseqid": ["r1", "r2", "r3", "r4"],
            "sseqid": ["g1"] * 4,
            "pident": [99.0] * 4,
            "length": [150] * 4,
            "genome_name": ["GCF_000123456.1"] * 4,
        })
        filtered_df = raw_df.head(1)  # 75% filtered

        qc = compute_pre_qc(raw_df, filtered_df, ani_matrix, config)

        assert qc.filter_rate == 0.75
        assert any("filter rate" in w.lower() for w in qc.warnings)

    def test_low_genome_coverage_warning(self, config):
        """Should warn when genome coverage is low with >3 genomes."""
        from metadarkmatter.core.ani_placement import ANIMatrix

        # Need >3 genomes for the warning to trigger
        ani_dict = {
            f"G{i}": {f"G{j}": 100.0 if i == j else 85.0 for j in range(5)}
            for i in range(5)
        }
        large_ani = ANIMatrix(ani_dict)

        raw_df = pl.DataFrame({
            "qseqid": ["r1"],
            "sseqid": ["g1"],
            "pident": [99.0],
            "length": [150],
            "genome_name": ["G0"],
        })

        qc = compute_pre_qc(raw_df, raw_df, large_ani, config)

        assert qc.genome_coverage < 0.5
        assert any("genome coverage" in w.lower() for w in qc.warnings)

    def test_high_single_hit_warning(self, ani_matrix, config):
        """Should warn when single-hit fraction exceeds 80%."""
        raw_df = pl.DataFrame({
            "qseqid": [f"r{i}" for i in range(10)],
            "sseqid": ["g1"] * 10,
            "pident": [99.0] * 10,
            "length": [150] * 10,
            "genome_name": ["GCF_000123456.1"] * 10,
        })

        qc = compute_pre_qc(raw_df, raw_df, ani_matrix, config)

        assert qc.single_hit_fraction == 1.0
        assert any("single-hit" in w.lower() for w in qc.warnings)

    def test_empty_filtered_df(self, ani_matrix, config):
        """Should handle empty filtered DataFrame."""
        raw_df = pl.DataFrame({
            "qseqid": ["r1"],
            "sseqid": ["g1"],
            "pident": [99.0],
            "length": [150],
            "genome_name": ["GCF_000123456.1"],
        })
        filtered_df = raw_df.clear()

        qc = compute_pre_qc(raw_df, filtered_df, ani_matrix, config)

        assert qc.genomes_in_alignment == 0
        assert qc.genome_coverage == 0.0


class TestPostQC:
    """Tests for compute_post_qc function."""

    def test_classification_counts(self):
        """Should count classification categories."""
        qc = QCMetrics()
        result_df = pl.DataFrame({
            "taxonomic_call": ["Known Species", "Known Species", "Novel Species", "Ambiguous"],
            "confidence_score": [0.9, 0.8, 0.7, 0.3],
        })

        qc = compute_post_qc(qc, result_df)

        assert qc.classification_counts["Known Species"] == 2
        assert qc.classification_counts["Novel Species"] == 1
        assert qc.classification_counts["Ambiguous"] == 1

    def test_ambiguous_fraction(self):
        """Should compute ambiguous fraction."""
        qc = QCMetrics()
        result_df = pl.DataFrame({
            "taxonomic_call": ["Known Species", "Ambiguous", "Ambiguous", "Ambiguous"],
            "confidence_score": [0.9, 0.3, 0.2, 0.4],
        })

        qc = compute_post_qc(qc, result_df)

        assert qc.ambiguous_fraction == 0.75

    def test_novel_fraction(self):
        """Should compute fraction of novel reads."""
        qc = QCMetrics()
        result_df = pl.DataFrame({
            "taxonomic_call": ["Known Species", "Novel Species", "Novel Genus"],
            "confidence_score": [0.9, 0.7, 0.6],
        })

        qc = compute_post_qc(qc, result_df)

        assert qc.novel_fraction == pytest.approx(2 / 3, rel=0.01)

    def test_low_confidence_fraction(self):
        """Should compute fraction of low-confidence reads."""
        qc = QCMetrics()
        result_df = pl.DataFrame({
            "taxonomic_call": ["Known Species", "Known Species", "Known Species"],
            "confidence_score": [0.9, 0.3, 0.2],
        })

        qc = compute_post_qc(qc, result_df)

        assert qc.low_confidence_fraction == pytest.approx(2 / 3, rel=0.01)

    def test_high_ambiguous_warning(self):
        """Should warn when ambiguous fraction exceeds 50%."""
        qc = QCMetrics()
        result_df = pl.DataFrame({
            "taxonomic_call": ["Ambiguous"] * 6 + ["Known Species"] * 4,
            "confidence_score": [0.3] * 6 + [0.9] * 4,
        })

        qc = compute_post_qc(qc, result_df)

        assert any("ambiguous" in w.lower() for w in qc.warnings)

    def test_empty_result_df(self):
        """Should handle empty result DataFrame."""
        qc = QCMetrics(total_alignments=100)
        result_df = pl.DataFrame({
            "taxonomic_call": [],
            "confidence_score": [],
        }).cast({"taxonomic_call": pl.Utf8, "confidence_score": pl.Float64})

        qc = compute_post_qc(qc, result_df)

        # Pre-QC fields preserved
        assert qc.total_alignments == 100
        # Post-QC fields empty
        assert qc.classification_counts == {}
