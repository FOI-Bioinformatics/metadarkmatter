"""
Unit tests for ANI-weighted placement classifier.

Tests ANIMatrix, ANIWeightedClassifier, VectorizedClassifier including:
- ANI matrix operations and lookups
- Classification threshold logic
- Novelty index and placement uncertainty calculations
- Integration with BLAST parsing
"""

from __future__ import annotations

import polars as pl
import pytest

from metadarkmatter.core.ani_placement import (
    ANIMatrix,
    ANIWeightedClassifier,
    VectorizedClassifier,
)
from metadarkmatter.models.blast import BlastHit, BlastResult
from metadarkmatter.models.classification import TaxonomicCall
from metadarkmatter.models.config import ScoringConfig


class TestANIMatrix:
    """Tests for ANIMatrix class."""

    def test_create_from_dict(self, small_ani_dict):
        """Should create ANIMatrix from nested dictionary."""
        matrix = ANIMatrix(small_ani_dict)

        assert len(matrix) == 3

    def test_genomes_property(self, small_ani_dict):
        """genomes property should return set of genome names."""
        matrix = ANIMatrix(small_ani_dict)
        genomes = matrix.genomes

        assert isinstance(genomes, set)
        assert "GCF_000123456.1" in genomes
        assert "GCF_000789012.1" in genomes
        assert "GCA_000111222.1" in genomes

    def test_get_ani_same_genome(self, small_ani_dict):
        """get_ani with same genome should return 100.0."""
        matrix = ANIMatrix(small_ani_dict)

        ani = matrix.get_ani("GCF_000123456.1", "GCF_000123456.1")

        assert ani == 100.0

    def test_get_ani_different_genomes(self, small_ani_dict):
        """get_ani should return correct ANI value."""
        matrix = ANIMatrix(small_ani_dict)

        ani = matrix.get_ani("GCF_000123456.1", "GCF_000789012.1")

        assert ani == 95.5

    def test_get_ani_symmetric(self, small_ani_dict):
        """get_ani should be symmetric."""
        matrix = ANIMatrix(small_ani_dict)

        ani_ab = matrix.get_ani("GCF_000123456.1", "GCF_000789012.1")
        ani_ba = matrix.get_ani("GCF_000789012.1", "GCF_000123456.1")

        assert ani_ab == ani_ba

    def test_get_ani_unknown_genome(self, small_ani_dict):
        """get_ani with unknown genome should return 0.0."""
        matrix = ANIMatrix(small_ani_dict)

        ani = matrix.get_ani("GCF_000123456.1", "unknown_genome")

        assert ani == 0.0

    def test_get_ani_both_unknown(self, small_ani_dict):
        """get_ani with both genomes unknown should return 0.0."""
        matrix = ANIMatrix(small_ani_dict)

        ani = matrix.get_ani("unknown_1", "unknown_2")

        assert ani == 0.0

    def test_has_genome_true(self, small_ani_dict):
        """has_genome should return True for known genome."""
        matrix = ANIMatrix(small_ani_dict)

        assert matrix.has_genome("GCF_000123456.1") is True

    def test_has_genome_false(self, small_ani_dict):
        """has_genome should return False for unknown genome."""
        matrix = ANIMatrix(small_ani_dict)

        assert matrix.has_genome("unknown_genome") is False

    def test_get_genome_idx(self, small_ani_dict):
        """get_genome_idx should return integer index."""
        matrix = ANIMatrix(small_ani_dict)

        idx = matrix.get_genome_idx("GCF_000123456.1")

        assert isinstance(idx, int)
        assert 0 <= idx < len(matrix)

    def test_get_genome_idx_unknown(self, small_ani_dict):
        """get_genome_idx should return None for unknown genome."""
        matrix = ANIMatrix(small_ani_dict)

        idx = matrix.get_genome_idx("unknown_genome")

        assert idx is None

    def test_get_ani_by_idx(self, small_ani_dict):
        """get_ani_by_idx should work with integer indices."""
        matrix = ANIMatrix(small_ani_dict)

        idx1 = matrix.get_genome_idx("GCF_000123456.1")
        idx2 = matrix.get_genome_idx("GCF_000789012.1")

        ani = matrix.get_ani_by_idx(idx1, idx2)

        assert ani == 95.5

    def test_from_file(self, temp_ani_file):
        """from_file should load ANI matrix from CSV."""
        matrix = ANIMatrix.from_file(temp_ani_file)

        assert len(matrix) == 3
        assert matrix.has_genome("GCF_000123456.1")

    def test_memory_usage_bytes(self, small_ani_dict):
        """memory_usage_bytes should return positive integer."""
        matrix = ANIMatrix(small_ani_dict)

        mem = matrix.memory_usage_bytes()

        assert isinstance(mem, int)
        assert mem > 0


class TestANIWeightedClassifier:
    """Tests for ANIWeightedClassifier class."""

    def test_create_with_default_config(self, ani_matrix):
        """Should create classifier with default config."""
        classifier = ANIWeightedClassifier(ani_matrix)

        assert classifier.config.bitscore_threshold_pct == 95.0

    def test_create_with_custom_config(self, ani_matrix, custom_scoring_config):
        """Should create classifier with custom config."""
        classifier = ANIWeightedClassifier(ani_matrix, config=custom_scoring_config)

        assert classifier.config.bitscore_threshold_pct == 90.0

    def test_classify_empty_result(self, classifier):
        """classify_read with no hits should return None."""
        result = BlastResult(read_id="read_001", hits=())

        classification = classifier.classify_read(result)

        assert classification is None

    def test_classify_single_hit_known_species(self, classifier):
        """Single hit with high identity should be known species."""
        hit = BlastHit(
            qseqid="read_001",
            sseqid="GCF_000123456.1",
            pident=99.0,  # Novelty = 1.0 (< 2.0)
            length=150,
            mismatch=1,
            gapopen=0,
            qstart=1,
            qend=150,
            sstart=1000,
            send=1150,
            evalue=1e-50,
            bitscore=250.0,
        )
        result = BlastResult(read_id="read_001", hits=(hit,))

        classification = classifier.classify_read(result)

        assert classification is not None
        assert classification.taxonomic_call == TaxonomicCall.KNOWN_SPECIES
        assert classification.novelty_index == 1.0
        assert classification.placement_uncertainty == 0.0  # Single hit

    def test_classify_novel_species(self, ani_matrix):
        """Hit with moderate identity should be novel species.

        Uses high single_hit_uncertainty_threshold to test novelty classification
        independently of Rule 0b (single-hit inferred uncertainty).
        """
        config = ScoringConfig(single_hit_uncertainty_threshold=100.0)
        novelty_classifier = ANIWeightedClassifier(ani_matrix, config=config)

        hit = BlastHit(
            qseqid="read_001",
            sseqid="GCF_000123456.1",
            pident=90.0,  # Novelty = 10.0 (4-20)
            length=150,
            mismatch=15,
            gapopen=0,
            qstart=1,
            qend=150,
            sstart=1000,
            send=1150,
            evalue=1e-40,
            bitscore=200.0,
        )
        result = BlastResult(read_id="read_001", hits=(hit,))

        classification = novelty_classifier.classify_read(result)

        assert classification is not None
        assert classification.taxonomic_call == TaxonomicCall.NOVEL_SPECIES
        assert classification.novelty_index == 10.0

    def test_classify_novel_genus(self, ani_matrix):
        """Hit with low identity should be novel genus.

        Uses high single_hit_uncertainty_threshold to test novelty classification
        independently of Rule 0b (single-hit inferred uncertainty).
        """
        config = ScoringConfig(single_hit_uncertainty_threshold=100.0)
        novelty_classifier = ANIWeightedClassifier(ani_matrix, config=config)

        hit = BlastHit(
            qseqid="read_001",
            sseqid="GCF_000123456.1",
            pident=80.0,  # Novelty = 20.0 (20-25)
            length=150,
            mismatch=30,
            gapopen=0,
            qstart=1,
            qend=150,
            sstart=1000,
            send=1150,
            evalue=1e-30,
            bitscore=150.0,
        )
        result = BlastResult(read_id="read_001", hits=(hit,))

        classification = novelty_classifier.classify_read(result)

        assert classification is not None
        assert classification.taxonomic_call == TaxonomicCall.NOVEL_GENUS
        assert classification.novelty_index == 20.0

    def test_classify_with_ambiguous_hits(self, classifier):
        """Multiple ambiguous hits should affect placement uncertainty."""
        hits = (
            BlastHit(
                qseqid="read_001",
                sseqid="GCF_000123456.1",
                pident=99.0,
                length=150,
                mismatch=1,
                gapopen=0,
                qstart=1,
                qend=150,
                sstart=1000,
                send=1150,
                evalue=1e-50,
                bitscore=250.0,
            ),
            BlastHit(
                qseqid="read_001",
                sseqid="GCF_000789012.1",
                pident=98.0,
                length=150,
                mismatch=3,
                gapopen=0,
                qstart=1,
                qend=150,
                sstart=2000,
                send=2150,
                evalue=1e-48,
                bitscore=245.0,  # Within 95% of 250
            ),
        )
        result = BlastResult(read_id="read_001", hits=hits)

        classification = classifier.classify_read(result)

        assert classification is not None
        assert classification.num_ambiguous_hits == 2
        # Placement uncertainty = 100 - ANI(GCF_000123456.1, GCF_000789012.1) = 100 - 95.5 = 4.5
        assert classification.placement_uncertainty == pytest.approx(4.5, abs=0.1)

    def test_classify_ambiguous_high_uncertainty(self, classifier):
        """High placement uncertainty indicates ambiguous classification.

        Note: CONSERVED_REGION is only assigned when hits span multiple genera,
        which requires genus metadata (only available in Polars classifier).
        The standard classifier returns AMBIGUOUS for high uncertainty.
        """
        # Create hits with high uncertainty (low ANI between genomes)
        # AND large identity gap (> 2%) so identity-based ambiguity doesn't trigger
        hits = (
            BlastHit(
                qseqid="read_001",
                sseqid="GCF_000123456.1",
                pident=99.0,
                length=150,
                mismatch=1,
                gapopen=0,
                qstart=1,
                qend=150,
                sstart=1000,
                send=1150,
                evalue=1e-50,
                bitscore=250.0,
            ),
            BlastHit(
                qseqid="read_001",
                sseqid="GCA_000111222.1",  # ANI = 80.0 to GCF_000123456.1
                pident=95.0,  # Identity gap = 99.0 - 95.0 = 4.0 (> 2.0 threshold)
                length=150,
                mismatch=7,
                gapopen=0,
                qstart=1,
                qend=150,
                sstart=2000,
                send=2150,
                evalue=1e-45,
                bitscore=240.0,  # Within 95% of 250 (237.5)
            ),
        )
        result = BlastResult(read_id="read_001", hits=hits)

        classification = classifier.classify_read(result)

        # Uncertainty = 100 - 80.0 = 20.0 (>= 5.0 conserved threshold)
        # High uncertainty results in AMBIGUOUS (CONSERVED_REGION requires cross-genera check)
        assert classification.placement_uncertainty == 20.0
        assert classification.taxonomic_call == TaxonomicCall.AMBIGUOUS

class TestClassificationThresholds:
    """Tests for classification threshold boundaries."""

    @pytest.fixture
    def boundary_classifier(self, ani_matrix):
        """Classifier with default thresholds for boundary testing.

        Uses high single_hit_uncertainty_threshold to test novelty boundaries
        independently of Rule 0b (single-hit inferred uncertainty).
        """
        return ANIWeightedClassifier(
            ani_matrix,
            config=ScoringConfig(single_hit_uncertainty_threshold=100.0),
        )

    def test_known_species_boundary_novelty_2(self, boundary_classifier):
        """Novelty = 2.0 should still be known species."""
        hit = BlastHit(
            qseqid="read_001",
            sseqid="GCF_000123456.1",
            pident=98.0,  # Novelty = 2.0 (== novelty_known_max)
            length=150,
            mismatch=3,
            gapopen=0,
            qstart=1,
            qend=150,
            sstart=1000,
            send=1150,
            evalue=1e-50,
            bitscore=250.0,
        )
        result = BlastResult(read_id="read_001", hits=(hit,))

        classification = boundary_classifier.classify_read(result)

        # Novelty < 2.0 is known, so 2.0 should NOT be known species
        # It falls into the gap between known and novel species
        assert classification.novelty_index == 2.0

    def test_novel_species_boundary_novelty_5(self, boundary_classifier):
        """Novelty = 5.0 should be novel species."""
        hit = BlastHit(
            qseqid="read_001",
            sseqid="GCF_000123456.1",
            pident=95.0,  # Novelty = 5.0 (== novelty_novel_species_min)
            length=150,
            mismatch=7,
            gapopen=0,
            qstart=1,
            qend=150,
            sstart=1000,
            send=1150,
            evalue=1e-45,
            bitscore=230.0,
        )
        result = BlastResult(read_id="read_001", hits=(hit,))

        classification = boundary_classifier.classify_read(result)

        assert classification.novelty_index == 5.0
        assert classification.taxonomic_call == TaxonomicCall.NOVEL_SPECIES

    def test_novel_species_midrange_novelty_15(self, boundary_classifier):
        """Novelty = 15.0 should be novel species (mid-range, well below 20% boundary)."""
        hit = BlastHit(
            qseqid="read_001",
            sseqid="GCF_000123456.1",
            pident=85.0,  # Novelty = 15.0 (within 5-20% range for novel species)
            length=150,
            mismatch=22,
            gapopen=0,
            qstart=1,
            qend=150,
            sstart=1000,
            send=1150,
            evalue=1e-35,
            bitscore=180.0,
        )
        result = BlastResult(read_id="read_001", hits=(hit,))

        classification = boundary_classifier.classify_read(result)

        assert classification.novelty_index == 15.0
        # Novelty 15% is clearly Novel Species (threshold is 20% for genus)
        assert classification.taxonomic_call == TaxonomicCall.NOVEL_SPECIES

    def test_novel_species_boundary_novelty_19_9(self, boundary_classifier):
        """Novelty = 19.9 should be novel species (just below 20% boundary)."""
        hit = BlastHit(
            qseqid="read_001",
            sseqid="GCF_000123456.1",
            pident=80.1,  # Novelty = 19.9 (just below 20% boundary)
            length=150,
            mismatch=29,
            gapopen=0,
            qstart=1,
            qend=150,
            sstart=1000,
            send=1150,
            evalue=1e-30,
            bitscore=160.0,
        )
        result = BlastResult(read_id="read_001", hits=(hit,))

        classification = boundary_classifier.classify_read(result)

        assert classification.novelty_index == pytest.approx(19.9, rel=0.01)
        # Just below 20% is Novel Species
        assert classification.taxonomic_call == TaxonomicCall.NOVEL_SPECIES

    def test_novel_genus_boundary_novelty_20(self, boundary_classifier):
        """Novelty = 20.0 should be novel genus (exactly at boundary)."""
        hit = BlastHit(
            qseqid="read_001",
            sseqid="GCF_000123456.1",
            pident=80.0,  # Novelty = 20.0 (exactly at genus boundary)
            length=150,
            mismatch=30,
            gapopen=0,
            qstart=1,
            qend=150,
            sstart=1000,
            send=1150,
            evalue=1e-28,
            bitscore=150.0,
        )
        result = BlastResult(read_id="read_001", hits=(hit,))

        classification = boundary_classifier.classify_read(result)

        assert classification.novelty_index == 20.0
        # At 20% boundary, should be Novel Genus
        assert classification.taxonomic_call == TaxonomicCall.NOVEL_GENUS

    def test_novel_genus_boundary_novelty_25(self, boundary_classifier):
        """Novelty = 25.0 should be novel genus (max boundary)."""
        hit = BlastHit(
            qseqid="read_001",
            sseqid="GCF_000123456.1",
            pident=75.0,  # Novelty = 25.0 (== novelty_novel_genus_max)
            length=150,
            mismatch=37,
            gapopen=0,
            qstart=1,
            qend=150,
            sstart=1000,
            send=1150,
            evalue=1e-25,
            bitscore=140.0,
        )
        result = BlastResult(read_id="read_001", hits=(hit,))

        classification = boundary_classifier.classify_read(result)

        assert classification.novelty_index == 25.0
        assert classification.taxonomic_call == TaxonomicCall.NOVEL_GENUS


class TestVectorizedClassifier:
    """Tests for VectorizedClassifier class."""

    def test_create_vectorized_classifier(self, ani_matrix, default_scoring_config):
        """Should create VectorizedClassifier."""
        vectorized = VectorizedClassifier(ani_matrix, config=default_scoring_config)

        assert vectorized.ani_matrix == ani_matrix
        assert vectorized.config == default_scoring_config

    def test_classify_file(self, ani_matrix, default_scoring_config, temp_blast_file):
        """classify_file should return DataFrame."""
        vectorized = VectorizedClassifier(ani_matrix, config=default_scoring_config)

        df = vectorized.classify_file(temp_blast_file)

        assert isinstance(df, pl.DataFrame)
        assert "read_id" in df.columns
        assert "taxonomic_call" in df.columns
        assert len(df) > 0

    def test_classify_file_has_expected_columns(
        self, ani_matrix, default_scoring_config, temp_blast_file
    ):
        """VectorizedClassifier should produce DataFrame with expected columns."""
        vectorized = VectorizedClassifier(ani_matrix, config=default_scoring_config)

        df = vectorized.classify_file(temp_blast_file)

        assert isinstance(df, pl.DataFrame)
        assert len(df) > 0
        for col in ["read_id", "taxonomic_call", "novelty_index", "placement_uncertainty"]:
            assert col in df.columns

    def test_stream_to_file_csv(
        self, ani_matrix, default_scoring_config, temp_blast_file, temp_dir
    ):
        """stream_to_file should write CSV."""
        vectorized = VectorizedClassifier(ani_matrix, config=default_scoring_config)
        output_path = temp_dir / "streamed.csv"

        num_classified = vectorized.stream_to_file(
            temp_blast_file, output_path, output_format="csv"
        )

        assert output_path.exists()
        assert num_classified > 0

    def test_stream_to_file_parquet(
        self, ani_matrix, default_scoring_config, temp_blast_file, temp_dir
    ):
        """stream_to_file should write Parquet."""
        vectorized = VectorizedClassifier(ani_matrix, config=default_scoring_config)
        output_path = temp_dir / "streamed.parquet"

        num_classified = vectorized.stream_to_file(
            temp_blast_file, output_path, output_format="parquet"
        )

        assert output_path.exists()
        assert num_classified > 0

        # Verify content
        df = pl.read_parquet(output_path)
        assert len(df) == num_classified

    def test_empty_blast_file(self, ani_matrix, default_scoring_config, empty_blast_file):
        """Should raise ValueError for empty BLAST file during parsing."""
        vectorized = VectorizedClassifier(ani_matrix, config=default_scoring_config)

        # Parser raises ValueError during column detection for empty files
        with pytest.raises(ValueError, match="no data lines"):
            vectorized.classify_file(empty_blast_file)


class TestUncertaintyModes:
    """Tests for uncertainty mode calculation (max vs second)."""

    @pytest.fixture
    def uncertainty_test_hits(self):
        """Create hits for testing uncertainty modes.

        Scenario:
        - Best hit: GCF_000123456.1 (bitscore 250)
        - Second best: GCF_000789012.1 (bitscore 245, ANI 95.5 to best)
        - Third best: GCA_000111222.1 (bitscore 240, ANI 80.0 to best)

        Max mode: U = 100 - max(95.5, 80.0) = 100 - 95.5 = 4.5
        Second mode: U = 100 - 95.5 = 4.5 (same because second-best has max ANI)
        """
        return (
            BlastHit(
                qseqid="read_001",
                sseqid="GCF_000123456.1",
                pident=99.0,
                length=150,
                mismatch=1,
                gapopen=0,
                qstart=1,
                qend=150,
                sstart=1000,
                send=1150,
                evalue=1e-50,
                bitscore=250.0,
            ),
            BlastHit(
                qseqid="read_001",
                sseqid="GCF_000789012.1",  # ANI 95.5 to best
                pident=98.0,
                length=150,
                mismatch=3,
                gapopen=0,
                qstart=1,
                qend=150,
                sstart=2000,
                send=2150,
                evalue=1e-48,
                bitscore=245.0,
            ),
            BlastHit(
                qseqid="read_001",
                sseqid="GCA_000111222.1",  # ANI 80.0 to best
                pident=96.0,
                length=150,
                mismatch=6,
                gapopen=0,
                qstart=1,
                qend=150,
                sstart=3000,
                send=3150,
                evalue=1e-46,
                bitscore=240.0,
            ),
        )

    @pytest.fixture
    def uncertainty_test_hits_different(self):
        """Create hits where max and second modes give different results.

        Scenario:
        - Best hit: GCF_000123456.1 (bitscore 250)
        - Second best: GCA_000111222.1 (bitscore 245, ANI 80.0 to best) - LOWER ANI
        - Third best: GCF_000789012.1 (bitscore 240, ANI 95.5 to best) - HIGHER ANI

        95% threshold = 250 * 0.95 = 237.5, so all three hits are within threshold.

        Max mode: U = 100 - max(80.0, 95.5) = 100 - 95.5 = 4.5
        Second mode: U = 100 - 80.0 = 20.0 (uses second-best genome's ANI)
        """
        return (
            BlastHit(
                qseqid="read_001",
                sseqid="GCF_000123456.1",
                pident=99.0,
                length=150,
                mismatch=1,
                gapopen=0,
                qstart=1,
                qend=150,
                sstart=1000,
                send=1150,
                evalue=1e-50,
                bitscore=250.0,
            ),
            BlastHit(
                qseqid="read_001",
                sseqid="GCA_000111222.1",  # ANI 80.0 to best - second by bitscore
                pident=97.0,
                length=150,
                mismatch=4,
                gapopen=0,
                qstart=1,
                qend=150,
                sstart=3000,
                send=3150,
                evalue=1e-47,
                bitscore=245.0,
            ),
            BlastHit(
                qseqid="read_001",
                sseqid="GCF_000789012.1",  # ANI 95.5 to best - third by bitscore
                pident=95.0,
                length=150,
                mismatch=7,
                gapopen=0,
                qstart=1,
                qend=150,
                sstart=2000,
                send=2150,
                evalue=1e-44,
                bitscore=240.0,  # Within 95% threshold (237.5)
            ),
        )

    def test_max_mode_uses_maximum_ani(self, ani_matrix, uncertainty_test_hits):
        """Max mode should use maximum ANI among all competing genomes."""
        config = ScoringConfig(uncertainty_mode="max")
        classifier = ANIWeightedClassifier(ani_matrix, config=config)

        result = BlastResult(read_id="read_001", hits=uncertainty_test_hits)
        classification = classifier.classify_read(result)

        # Max ANI to any competitor is 95.5, so uncertainty = 4.5
        assert classification.placement_uncertainty == pytest.approx(4.5, abs=0.1)

    def test_second_mode_uses_second_best_ani(self, ani_matrix, uncertainty_test_hits):
        """Second mode should use ANI to second-best genome by bitscore."""
        config = ScoringConfig(uncertainty_mode="second")
        classifier = ANIWeightedClassifier(ani_matrix, config=config)

        result = BlastResult(read_id="read_001", hits=uncertainty_test_hits)
        classification = classifier.classify_read(result)

        # Second-best genome (by bitscore) is GCF_000789012.1 with ANI 95.5
        # So uncertainty = 100 - 95.5 = 4.5
        assert classification.placement_uncertainty == pytest.approx(4.5, abs=0.1)

    def test_modes_differ_when_second_has_lower_ani(
        self, ani_matrix, uncertainty_test_hits_different
    ):
        """Max and second modes should give different results when appropriate.

        When the second-best genome (by bitscore) has lower ANI than other
        competitors, max mode gives lower uncertainty than second mode.
        """
        config_max = ScoringConfig(uncertainty_mode="max")
        config_second = ScoringConfig(uncertainty_mode="second")

        classifier_max = ANIWeightedClassifier(ani_matrix, config=config_max)
        classifier_second = ANIWeightedClassifier(ani_matrix, config=config_second)

        result = BlastResult(read_id="read_001", hits=uncertainty_test_hits_different)

        classification_max = classifier_max.classify_read(result)
        classification_second = classifier_second.classify_read(result)

        # Max mode: uses max ANI (95.5), uncertainty = 4.5
        assert classification_max.placement_uncertainty == pytest.approx(4.5, abs=0.1)

        # Second mode: uses ANI to second-best genome (80.0), uncertainty = 20.0
        assert classification_second.placement_uncertainty == pytest.approx(20.0, abs=0.1)

        # Results should differ
        assert classification_max.placement_uncertainty < classification_second.placement_uncertainty

    def test_single_hit_uncertainty_zero_both_modes(self, ani_matrix):
        """Single hit should have zero uncertainty in both modes."""
        hit = BlastHit(
            qseqid="read_001",
            sseqid="GCF_000123456.1",
            pident=99.0,
            length=150,
            mismatch=1,
            gapopen=0,
            qstart=1,
            qend=150,
            sstart=1000,
            send=1150,
            evalue=1e-50,
            bitscore=250.0,
        )
        result = BlastResult(read_id="read_001", hits=(hit,))

        config_max = ScoringConfig(uncertainty_mode="max")
        config_second = ScoringConfig(uncertainty_mode="second")

        classifier_max = ANIWeightedClassifier(ani_matrix, config=config_max)
        classifier_second = ANIWeightedClassifier(ani_matrix, config=config_second)

        classification_max = classifier_max.classify_read(result)
        classification_second = classifier_second.classify_read(result)

        assert classification_max.placement_uncertainty == 0.0
        assert classification_second.placement_uncertainty == 0.0


class TestStreamingChunks:
    """Tests for streaming chunk writing."""

    @pytest.fixture
    def large_blast_file(self, temp_dir):
        """Create a BLAST file large enough to trigger multiple chunks."""
        blast_path = temp_dir / "large.blast.tsv"

        # Create 100 reads with multiple hits each
        rows = []
        for i in range(100):
            for hit in range(3):
                rows.append({
                    "qseqid": f"read_{i:04d}",
                    "sseqid": ["GCF_000123456.1", "GCF_000789012.1", "GCA_000111222.1"][hit],
                    "pident": [99.0, 95.5, 80.0][hit],
                    "length": 150,
                    "mismatch": 3,
                    "gapopen": 0,
                    "qstart": 1,
                    "qend": 150,
                    "sstart": 1000,
                    "send": 1150,
                    "evalue": 1e-50,
                    "bitscore": [250.0, 240.0, 200.0][hit],
                })

        df = pl.DataFrame(rows)
        df.write_csv(blast_path, separator="\t", include_header=False)
        return blast_path

    def test_stream_csv(
        self, ani_matrix, default_scoring_config, large_blast_file, temp_dir
    ):
        """Streaming should handle CSV format."""
        vectorized = VectorizedClassifier(ani_matrix, config=default_scoring_config)
        output_path = temp_dir / "streamed.csv"

        num_classified = vectorized.stream_to_file(
            large_blast_file,
            output_path,
            output_format="csv",
        )

        assert output_path.exists()
        assert num_classified == 100

        # Verify all data written correctly
        df = pl.read_csv(output_path)
        assert len(df) == 100

    def test_stream_parquet(
        self, ani_matrix, default_scoring_config, large_blast_file, temp_dir
    ):
        """Streaming should handle Parquet format."""
        vectorized = VectorizedClassifier(ani_matrix, config=default_scoring_config)
        output_path = temp_dir / "streamed.parquet"

        num_classified = vectorized.stream_to_file(
            large_blast_file,
            output_path,
            output_format="parquet",
        )

        assert output_path.exists()
        assert num_classified == 100

        df = pl.read_parquet(output_path)
        assert len(df) == 100

    def test_stream_with_progress_callback(
        self, ani_matrix, default_scoring_config, large_blast_file, temp_dir
    ):
        """Streaming should call progress callback."""
        vectorized = VectorizedClassifier(ani_matrix, config=default_scoring_config)
        output_path = temp_dir / "progress.csv"

        progress_calls = []

        def track_progress(rows: int, reads: int, elapsed: float):
            progress_calls.append((rows, reads, elapsed))

        num_classified = vectorized.stream_to_file(
            large_blast_file,
            output_path,
            output_format="csv",
            progress_callback=track_progress,
        )

        # Progress callback should have been invoked at least once
        assert len(progress_calls) > 0
