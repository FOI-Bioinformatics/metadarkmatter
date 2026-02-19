"""Tests for family validation in VectorizedClassifier."""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from metadarkmatter.core.classification.ani_matrix import ANIMatrix
from metadarkmatter.core.classification.classifiers.vectorized import VectorizedClassifier
from metadarkmatter.models.config import ScoringConfig


def _make_ani_matrix(genomes: list[str], default_ani: float = 70.0) -> ANIMatrix:
    """Create a minimal ANIMatrix for testing."""
    n = len(genomes)
    arr = np.full((n, n), default_ani, dtype=np.float32)
    np.fill_diagonal(arr, 100.0)
    # Make within-family genomes related
    for i in range(n):
        for j in range(n):
            if i != j:
                arr[i, j] = 85.0
    matrix = ANIMatrix.__new__(ANIMatrix)
    matrix._genomes = tuple(genomes)
    matrix._genome_to_idx = {g: i for i, g in enumerate(genomes)}
    matrix._ani_array = arr
    matrix._default_ani = default_ani
    matrix._num_genomes = n
    return matrix


class TestFamilyValidationDisabled:
    """When target_family is None, output should have no family columns."""

    def test_no_family_columns_without_target_family(self, tmp_path):
        """Output should not contain family columns when target_family is None."""
        genomes = ["GCF_001", "GCF_002"]
        ani = _make_ani_matrix(genomes)
        config = ScoringConfig()  # No target_family

        blast_file = tmp_path / "test.tsv"
        blast_file.write_text(
            "read1\tGCF_001|contig1\t98.0\t150\t3\t0\t1\t150\t1\t150\t1e-50\t500\n"
            "read1\tGCF_002|contig1\t95.0\t150\t8\t0\t1\t150\t1\t150\t1e-40\t450\n"
        )

        classifier = VectorizedClassifier(ani, config=config)
        result = classifier.classify_file(blast_file)

        assert "family_bitscore_ratio" not in result.columns
        assert "external_best_genome" not in result.columns


class TestFamilyValidationEnabled:
    """When target_family is set, family validation should work."""

    def test_all_in_family_ratio_is_one(self, tmp_path):
        """When all hits are in-family, family_bitscore_ratio should be 1.0."""
        genomes = ["GCF_001", "GCF_002"]
        ani = _make_ani_matrix(genomes)
        config = ScoringConfig(target_family="f__TestFamily")

        blast_file = tmp_path / "test.tsv"
        blast_file.write_text(
            "read1\tGCF_001|c1\t98.0\t150\t3\t0\t1\t150\t1\t150\t1e-50\t500\n"
            "read1\tGCF_002|c1\t95.0\t150\t8\t0\t1\t150\t1\t150\t1e-40\t450\n"
        )

        classifier = VectorizedClassifier(ani, config=config)
        result = classifier.classify_file(blast_file)

        assert "family_bitscore_ratio" in result.columns
        assert result["family_bitscore_ratio"][0] == pytest.approx(1.0)
        assert result["taxonomic_call"][0] != "Off-target"

    def test_external_better_hit_marks_off_target(self, tmp_path):
        """Reads with much better external hits should be Off-target."""
        genomes = ["GCF_001"]
        ani = _make_ani_matrix(genomes)
        config = ScoringConfig(target_family="f__TestFamily", family_ratio_threshold=0.8)

        # read1: in-family=200, external=600 -> ratio=0.33 -> Off-target
        blast_file = tmp_path / "test.tsv"
        blast_file.write_text(
            "read1\tGCF_001|c1\t75.0\t150\t38\t0\t1\t150\t1\t150\t1e-10\t200\n"
            "read1\tEXTERNAL_001|c1\t95.0\t150\t8\t0\t1\t150\t1\t150\t1e-50\t600\n"
        )

        classifier = VectorizedClassifier(ani, config=config)
        result = classifier.classify_file(blast_file)

        assert result["taxonomic_call"][0] == "Off-target"
        assert result["family_bitscore_ratio"][0] == pytest.approx(200.0 / 600.0, abs=0.01)

    def test_no_in_family_hits_marks_off_target(self, tmp_path):
        """Reads with only external hits should be Off-target."""
        genomes = ["GCF_001"]
        ani = _make_ani_matrix(genomes)
        config = ScoringConfig(target_family="f__TestFamily")

        blast_file = tmp_path / "test.tsv"
        blast_file.write_text(
            "read1\tEXTERNAL_001|c1\t92.0\t150\t12\t0\t1\t150\t1\t150\t1e-45\t550\n"
        )

        classifier = VectorizedClassifier(ani, config=config)
        result = classifier.classify_file(blast_file)

        assert result["taxonomic_call"][0] == "Off-target"
        assert result["family_bitscore_ratio"][0] == 0.0

    def test_mixed_reads_correct_classification(self, tmp_path):
        """Mix of in-family and off-target reads should be classified correctly."""
        genomes = ["GCF_001", "GCF_002"]
        ani = _make_ani_matrix(genomes)
        config = ScoringConfig(target_family="f__TestFamily", family_ratio_threshold=0.8)

        blast_file = tmp_path / "test.tsv"
        lines = [
            # read1: all in-family -> normal classification
            "read1\tGCF_001|c1\t98.0\t150\t3\t0\t1\t150\t1\t150\t1e-50\t500",
            "read1\tGCF_002|c1\t95.0\t150\t8\t0\t1\t150\t1\t150\t1e-40\t450",
            # read2: external much better -> Off-target
            "read2\tGCF_001|c1\t75.0\t150\t38\t0\t1\t150\t1\t150\t1e-10\t200",
            "read2\tEXTERNAL_001|c1\t95.0\t150\t8\t0\t1\t150\t1\t150\t1e-50\t600",
        ]
        blast_file.write_text("\n".join(lines) + "\n")

        classifier = VectorizedClassifier(ani, config=config)
        result = classifier.classify_file(blast_file)

        # Should have 2 reads
        assert result.height == 2

        read1 = result.filter(pl.col("read_id") == "read1")
        assert read1["taxonomic_call"][0] != "Off-target"

        read2 = result.filter(pl.col("read_id") == "read2")
        assert read2["taxonomic_call"][0] == "Off-target"

    def test_family_columns_present(self, tmp_path):
        """Family validation columns should be present in output."""
        genomes = ["GCF_001"]
        ani = _make_ani_matrix(genomes)
        config = ScoringConfig(target_family="f__TestFamily")

        blast_file = tmp_path / "test.tsv"
        blast_file.write_text(
            "read1\tGCF_001|c1\t98.0\t150\t3\t0\t1\t150\t1\t150\t1e-50\t500\n"
        )

        classifier = VectorizedClassifier(ani, config=config)
        result = classifier.classify_file(blast_file)

        expected_cols = [
            "family_bitscore_ratio",
            "family_identity_gap",
            "in_family_hit_fraction",
            "external_best_genome",
            "external_best_identity",
        ]
        for col in expected_cols:
            assert col in result.columns, f"Missing column: {col}"
