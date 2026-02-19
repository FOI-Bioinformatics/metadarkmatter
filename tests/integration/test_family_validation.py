"""Integration test for family validation with mixed in-family and external hits."""
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
    arr = np.full((n, n), 85.0, dtype=np.float32)
    np.fill_diagonal(arr, 100.0)
    matrix = ANIMatrix.__new__(ANIMatrix)
    matrix._genomes = tuple(genomes)
    matrix._num_genomes = n
    matrix._genome_to_idx = {g: i for i, g in enumerate(genomes)}
    matrix._ani_array = arr
    matrix._default_ani = default_ani
    return matrix


class TestFamilyValidationIntegration:
    """Full pipeline integration test with mixed hits."""

    def test_mixed_hits_classifies_correctly(self, tmp_path):
        """Reads with external hits should be classified as Off-target."""
        family_genomes = ["GCF_001", "GCF_002"]
        ani = _make_ani_matrix(family_genomes)

        blast_lines = [
            # read1: best hit is in-family -> normal classification
            "read1\tGCF_001|c1\t98.0\t150\t3\t0\t1\t150\t1\t150\t1e-50\t500",
            "read1\tGCF_002|c1\t95.0\t150\t8\t0\t1\t150\t1\t150\t1e-40\t450",
            # read2: external hit much better -> Off-target (ratio = 200/600 = 0.33)
            "read2\tGCF_001|c1\t75.0\t150\t38\t0\t1\t150\t1\t150\t1e-10\t200",
            "read2\tEXTERNAL_001|c1\t95.0\t150\t8\t0\t1\t150\t1\t150\t1e-50\t600",
            # read3: only external hits -> Off-target
            "read3\tEXTERNAL_002|c1\t92.0\t150\t12\t0\t1\t150\t1\t150\t1e-45\t550",
        ]
        blast_file = tmp_path / "mixed.tsv"
        blast_file.write_text("\n".join(blast_lines) + "\n")

        config = ScoringConfig(
            target_family="f__TestFamily",
            family_ratio_threshold=0.8,
        )
        classifier = VectorizedClassifier(ani, config=config)
        result = classifier.classify_file(blast_file)

        # read1 should be classified normally (not Off-target)
        read1 = result.filter(pl.col("read_id") == "read1")
        assert read1["taxonomic_call"][0] != "Off-target"
        assert read1["family_bitscore_ratio"][0] == pytest.approx(1.0)

        # read2 should be Off-target
        read2 = result.filter(pl.col("read_id") == "read2")
        assert read2["taxonomic_call"][0] == "Off-target"
        assert read2["family_bitscore_ratio"][0] == pytest.approx(200.0 / 600.0, abs=0.01)

        # read3 should be Off-target (no in-family hits)
        read3 = result.filter(pl.col("read_id") == "read3")
        assert read3["taxonomic_call"][0] == "Off-target"
        assert read3["family_bitscore_ratio"][0] == 0.0

    def test_all_external_hits(self, tmp_path):
        """When all reads are external, all should be Off-target."""
        family_genomes = ["GCF_001"]
        ani = _make_ani_matrix(family_genomes)

        blast_lines = [
            "read1\tEXTERNAL_001|c1\t92.0\t150\t12\t0\t1\t150\t1\t150\t1e-45\t550",
            "read2\tEXTERNAL_002|c1\t88.0\t150\t18\t0\t1\t150\t1\t150\t1e-30\t400",
        ]
        blast_file = tmp_path / "all_external.tsv"
        blast_file.write_text("\n".join(blast_lines) + "\n")

        config = ScoringConfig(target_family="f__TestFamily")
        classifier = VectorizedClassifier(ani, config=config)
        result = classifier.classify_file(blast_file)

        assert result.height == 2
        assert all(c == "Off-target" for c in result["taxonomic_call"].to_list())

    def test_borderline_ratio(self, tmp_path):
        """Reads with ratio exactly at threshold should not be Off-target."""
        family_genomes = ["GCF_001"]
        ani = _make_ani_matrix(family_genomes)

        # ratio = 400/500 = 0.8 exactly at threshold
        blast_lines = [
            "read1\tGCF_001|c1\t90.0\t150\t15\t0\t1\t150\t1\t150\t1e-30\t400",
            "read1\tEXTERNAL_001|c1\t92.0\t150\t12\t0\t1\t150\t1\t150\t1e-40\t500",
        ]
        blast_file = tmp_path / "borderline.tsv"
        blast_file.write_text("\n".join(blast_lines) + "\n")

        config = ScoringConfig(target_family="f__TestFamily", family_ratio_threshold=0.8)
        classifier = VectorizedClassifier(ani, config=config)
        result = classifier.classify_file(blast_file)

        # 0.8 is NOT < 0.8, so should NOT be off-target
        assert result["taxonomic_call"][0] != "Off-target"

    def test_family_columns_output(self, tmp_path):
        """All family validation columns should be present and correct types."""
        family_genomes = ["GCF_001", "GCF_002"]
        ani = _make_ani_matrix(family_genomes)

        blast_lines = [
            "read1\tGCF_001|c1\t98.0\t150\t3\t0\t1\t150\t1\t150\t1e-50\t500",
            "read1\tEXTERNAL_001|c1\t85.0\t150\t23\t0\t1\t150\t1\t150\t1e-20\t300",
        ]
        blast_file = tmp_path / "output_test.tsv"
        blast_file.write_text("\n".join(blast_lines) + "\n")

        config = ScoringConfig(target_family="f__TestFamily")
        classifier = VectorizedClassifier(ani, config=config)
        result = classifier.classify_file(blast_file)

        assert "family_bitscore_ratio" in result.columns
        assert "family_identity_gap" in result.columns
        assert "in_family_hit_fraction" in result.columns
        assert "external_best_genome" in result.columns
        assert "external_best_identity" in result.columns
