"""Tests for family validation in VectorizedClassifier."""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from metadarkmatter.core.classification.ani_matrix import ANIMatrix
from metadarkmatter.core.classification.classifiers.vectorized import VectorizedClassifier
from metadarkmatter.core.classification.thresholds import apply_classification_thresholds
from metadarkmatter.models.classification import TaxonomicSummary
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


class TestThresholdsWithOffTarget:
    """Off-target reads should be preserved through threshold reclassification."""

    def test_off_target_preserved(self):
        """Off-target reads should not be reclassified."""
        df = pl.DataFrame({
            "read_id": ["read_ot", "read_normal"],
            "novelty_index": [15.0, 2.0],
            "placement_uncertainty": [0.5, 0.5],
            "num_ambiguous_hits": [1, 2],
            "identity_gap": [None, 3.0],
            "taxonomic_call": ["Off-target", "Known Species"],
        })
        config = ScoringConfig(target_family="f__Test")
        result = apply_classification_thresholds(df, config)

        ot_row = result.filter(pl.col("read_id") == "read_ot")
        assert ot_row["taxonomic_call"][0] == "Off-target"

        normal_row = result.filter(pl.col("read_id") == "read_normal")
        assert normal_row["taxonomic_call"][0] in [
            "Known Species", "Ambiguous", "Novel Species", "Unclassified",
        ]

    def test_no_off_target_unchanged(self):
        """Without Off-target reads, function should work identically."""
        df = pl.DataFrame({
            "novelty_index": [2.0, 10.0],
            "placement_uncertainty": [0.5, 0.5],
            "num_ambiguous_hits": [2, 3],
            "identity_gap": [3.0, 5.0],
            "taxonomic_call": ["Known Species", "Novel Species"],
        })
        config = ScoringConfig()
        result = apply_classification_thresholds(df, config)
        assert "Off-target" not in result["taxonomic_call"].to_list()


class TestBackwardCompatibility:
    """Verify identical output when family validation is disabled."""

    def test_output_identical_without_target_family(self, tmp_path):
        """Classification should produce identical core results."""
        genomes = ["GCF_001", "GCF_002", "GCF_003"]
        ani = _make_ani_matrix(genomes)

        blast_file = tmp_path / "test.tsv"
        lines = [
            "read1\tGCF_001|c1\t98.0\t150\t3\t0\t1\t150\t1\t150\t1e-50\t500",
            "read1\tGCF_002|c1\t95.0\t150\t8\t0\t1\t150\t1\t150\t1e-40\t450",
            "read2\tGCF_003|c1\t88.0\t150\t18\t0\t1\t150\t1\t150\t1e-30\t300",
        ]
        blast_file.write_text("\n".join(lines) + "\n")

        # Run without family validation
        config_without = ScoringConfig()
        classifier_without = VectorizedClassifier(ani, config=config_without)
        result_without = classifier_without.classify_file(blast_file)

        # Run with family validation but all hits in-family
        config_with = ScoringConfig(target_family="f__TestFamily")
        classifier_with = VectorizedClassifier(ani, config=config_with)
        result_with = classifier_with.classify_file(blast_file)

        # Core classification columns should match
        core_cols = [
            "read_id", "best_match_genome", "top_hit_identity",
            "novelty_index", "placement_uncertainty", "taxonomic_call",
        ]
        for col in core_cols:
            assert result_without[col].to_list() == result_with[col].to_list(), (
                f"Column {col} differs"
            )

        # Family columns should only exist in result_with
        assert "family_bitscore_ratio" not in result_without.columns
        assert "family_bitscore_ratio" in result_with.columns


class TestOffTargetDiversityStatus:
    """Off-target reads should get diversity_status='Off-target', not 'Uncertain'."""

    def test_off_target_diversity_status_in_classifier(self, tmp_path):
        """Off-target reads from the classifier should have diversity_status='Off-target'."""
        genomes = ["GCF_001"]
        ani = _make_ani_matrix(genomes)
        config = ScoringConfig(target_family="f__TestFamily", family_ratio_threshold=0.8)

        blast_file = tmp_path / "test.tsv"
        blast_file.write_text(
            "read1\tGCF_001|c1\t75.0\t150\t38\t0\t1\t150\t1\t150\t1e-10\t200\n"
            "read1\tEXTERNAL_001|c1\t95.0\t150\t8\t0\t1\t150\t1\t150\t1e-50\t600\n"
        )

        classifier = VectorizedClassifier(ani, config=config)
        result = classifier.classify_file(blast_file)

        assert result["taxonomic_call"][0] == "Off-target"
        assert result["diversity_status"][0] == "Off-target"

    def test_off_target_not_grouped_as_uncertain(self, tmp_path):
        """Off-target reads should not inflate the Uncertain count."""
        from metadarkmatter.models.classification import TAXONOMIC_TO_DIVERSITY
        assert TAXONOMIC_TO_DIVERSITY["Off-target"] == "Off-target"

    def test_in_family_reads_still_get_correct_diversity_status(self, tmp_path):
        """In-family reads should still get Known/Novel/Uncertain status normally."""
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

        # In-family reads should not have diversity_status = "Off-target"
        assert result["diversity_status"][0] != "Off-target"


class TestFamilyValidationSummary:
    """Tests for family validation in summary JSON."""

    def test_summary_includes_off_target_count(self):
        """TaxonomicSummary should include off_target count."""
        summary = TaxonomicSummary(
            total_reads=100,
            known_species=50,
            novel_species=20,
            novel_genus=5,
            conserved_regions=3,
            off_target=10,
            mean_novelty_index=8.0,
            mean_placement_uncertainty=1.2,
        )
        assert summary.off_target == 10
        data = summary.model_dump()
        assert data["off_target"] == 10
