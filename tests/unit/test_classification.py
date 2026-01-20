"""
Unit tests for classification models.

Tests TaxonomicCall enum, ReadClassification, and TaxonomicSummary models
including computed fields and serialization.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metadarkmatter.models.classification import (
    ReadClassification,
    TaxonomicCall,
    TaxonomicSummary,
)


class TestTaxonomicCall:
    """Tests for TaxonomicCall enumeration."""

    def test_known_species_value(self):
        """Known Species should have correct string value."""
        assert TaxonomicCall.KNOWN_SPECIES.value == "Known Species"

    def test_novel_species_value(self):
        """Novel Species should have correct string value."""
        assert TaxonomicCall.NOVEL_SPECIES.value == "Novel Species"

    def test_novel_genus_value(self):
        """Novel Genus should have correct string value."""
        assert TaxonomicCall.NOVEL_GENUS.value == "Novel Genus"

    def test_conserved_region_value(self):
        """Conserved Region should have correct string value."""
        assert TaxonomicCall.CONSERVED_REGION.value == "Conserved Region"

    def test_enum_is_string(self):
        """TaxonomicCall should be a string enum."""
        assert isinstance(TaxonomicCall.KNOWN_SPECIES, str)
        assert TaxonomicCall.KNOWN_SPECIES == "Known Species"

    def test_all_values_unique(self):
        """All enum values should be unique."""
        values = [tc.value for tc in TaxonomicCall]
        assert len(values) == len(set(values))


class TestReadClassification:
    """Tests for ReadClassification Pydantic model."""

    @pytest.fixture
    def known_species_classification(self):
        """Classification result for known species."""
        return ReadClassification(
            read_id="read_001",
            best_match_genome="GCF_000123456.1",
            top_hit_identity=99.0,
            novelty_index=1.0,
            placement_uncertainty=0.3,
            num_ambiguous_hits=2,
            taxonomic_call=TaxonomicCall.KNOWN_SPECIES,
        )

    @pytest.fixture
    def novel_species_classification(self):
        """Classification result for novel species."""
        return ReadClassification(
            read_id="read_002",
            best_match_genome="GCF_000789012.1",
            top_hit_identity=90.0,
            novelty_index=10.0,
            placement_uncertainty=0.4,
            num_ambiguous_hits=3,
            taxonomic_call=TaxonomicCall.NOVEL_SPECIES,
        )

    @pytest.fixture
    def novel_genus_classification(self):
        """Classification result for novel genus."""
        return ReadClassification(
            read_id="read_003",
            best_match_genome="GCA_000111222.1",
            top_hit_identity=80.0,
            novelty_index=20.0,
            placement_uncertainty=1.5,
            num_ambiguous_hits=5,
            taxonomic_call=TaxonomicCall.NOVEL_GENUS,
        )

    def test_create_valid_classification(self, known_species_classification):
        """Should create ReadClassification from valid data."""
        assert known_species_classification.read_id == "read_001"
        assert known_species_classification.best_match_genome == "GCF_000123456.1"
        assert known_species_classification.top_hit_identity == 99.0
        assert known_species_classification.novelty_index == 1.0
        assert known_species_classification.placement_uncertainty == 0.3
        assert known_species_classification.num_ambiguous_hits == 2
        assert known_species_classification.taxonomic_call == TaxonomicCall.KNOWN_SPECIES

    def test_classification_is_frozen(self, known_species_classification):
        """ReadClassification should be immutable."""
        with pytest.raises(Exception):
            known_species_classification.read_id = "different_read"

    def test_is_novel_known_species(self, known_species_classification):
        """Known species should not be novel."""
        assert known_species_classification.is_novel is False

    def test_is_novel_novel_species(self, novel_species_classification):
        """Novel species should be novel."""
        assert novel_species_classification.is_novel is True

    def test_is_novel_novel_genus(self, novel_genus_classification):
        """Novel genus should be novel."""
        assert novel_genus_classification.is_novel is True

    def test_is_novel_conserved_region(self):
        """Conserved region should not be novel."""
        classification = ReadClassification(
            read_id="read_004",
            best_match_genome="GCF_000123456.1",
            top_hit_identity=95.0,
            novelty_index=5.0,
            placement_uncertainty=10.0,
            num_ambiguous_hits=10,
            taxonomic_call=TaxonomicCall.CONSERVED_REGION,
        )
        assert classification.is_novel is False

    def test_top_hit_identity_max_100(self):
        """top_hit_identity should not exceed 100."""
        with pytest.raises(Exception):
            ReadClassification(
                read_id="read_001",
                best_match_genome="GCF_000123456.1",
                top_hit_identity=105.0,  # Invalid
                novelty_index=0.0,
                placement_uncertainty=0.0,
                num_ambiguous_hits=1,
                taxonomic_call=TaxonomicCall.KNOWN_SPECIES,
            )

    def test_novelty_index_min_0(self):
        """novelty_index should not be negative."""
        with pytest.raises(Exception):
            ReadClassification(
                read_id="read_001",
                best_match_genome="GCF_000123456.1",
                top_hit_identity=100.0,
                novelty_index=-5.0,  # Invalid
                placement_uncertainty=0.0,
                num_ambiguous_hits=1,
                taxonomic_call=TaxonomicCall.KNOWN_SPECIES,
            )

    def test_placement_uncertainty_min_0(self):
        """placement_uncertainty should not be negative."""
        with pytest.raises(Exception):
            ReadClassification(
                read_id="read_001",
                best_match_genome="GCF_000123456.1",
                top_hit_identity=100.0,
                novelty_index=0.0,
                placement_uncertainty=-1.0,  # Invalid
                num_ambiguous_hits=1,
                taxonomic_call=TaxonomicCall.KNOWN_SPECIES,
            )

    def test_num_ambiguous_hits_min_0(self):
        """num_ambiguous_hits should not be negative."""
        with pytest.raises(Exception):
            ReadClassification(
                read_id="read_001",
                best_match_genome="GCF_000123456.1",
                top_hit_identity=100.0,
                novelty_index=0.0,
                placement_uncertainty=0.0,
                num_ambiguous_hits=-1,  # Invalid
                taxonomic_call=TaxonomicCall.KNOWN_SPECIES,
            )

    def test_to_dict_method(self, known_species_classification):
        """to_dict should return dictionary with all fields."""
        result = known_species_classification.to_dict()

        assert result["read_id"] == "read_001"
        assert result["best_match_genome"] == "GCF_000123456.1"
        assert result["top_hit_identity"] == 99.0
        assert result["novelty_index"] == 1.0
        assert result["placement_uncertainty"] == 0.3
        assert result["num_ambiguous_hits"] == 2
        assert result["taxonomic_call"] == "Known Species"
        assert result["is_novel"] is False

    def test_to_dict_novel_species(self, novel_species_classification):
        """to_dict should correctly set is_novel for novel species."""
        result = novel_species_classification.to_dict()

        assert result["taxonomic_call"] == "Novel Species"
        assert result["is_novel"] is True


class TestTaxonomicSummary:
    """Tests for TaxonomicSummary aggregation model."""

    @pytest.fixture
    def sample_summary(self):
        """Sample summary with known distribution."""
        return TaxonomicSummary(
            total_reads=1000,
            known_species=600,
            novel_species=150,
            novel_genus=50,
            conserved_regions=200,
            mean_novelty_index=5.5,
            mean_placement_uncertainty=1.2,
            genome_hit_counts={
                "GCF_000123456.1": 300,
                "GCF_000789012.1": 200,
                "GCA_000111222.1": 100,
            },
        )

    def test_create_valid_summary(self, sample_summary):
        """Should create TaxonomicSummary from valid data."""
        assert sample_summary.total_reads == 1000
        assert sample_summary.known_species == 600
        assert sample_summary.novel_species == 150
        assert sample_summary.novel_genus == 50
        assert sample_summary.conserved_regions == 200

    def test_known_species_pct(self, sample_summary):
        """known_species_pct should calculate correctly."""
        # 600 / 1000 * 100 = 60.0%
        assert sample_summary.known_species_pct == 60.0

    def test_novel_species_pct(self, sample_summary):
        """novel_species_pct should calculate correctly."""
        # 150 / 1000 * 100 = 15.0%
        assert sample_summary.novel_species_pct == 15.0

    def test_novel_genus_pct(self, sample_summary):
        """novel_genus_pct should calculate correctly."""
        # 50 / 1000 * 100 = 5.0%
        assert sample_summary.novel_genus_pct == 5.0

    def test_novel_diversity_pct(self, sample_summary):
        """novel_diversity_pct should sum novel species and genus."""
        # 15.0% + 5.0% = 20.0%
        assert sample_summary.novel_diversity_pct == 20.0

    def test_percentages_with_zero_reads(self):
        """Percentages should be 0 when total_reads is 0."""
        summary = TaxonomicSummary(
            total_reads=0,
            known_species=0,
            novel_species=0,
            novel_genus=0,
            conserved_regions=0,
            mean_novelty_index=0.0,
            mean_placement_uncertainty=0.0,
        )

        assert summary.known_species_pct == 0.0
        assert summary.novel_species_pct == 0.0
        assert summary.novel_genus_pct == 0.0
        assert summary.novel_diversity_pct == 0.0

    def test_empty_genome_hit_counts(self):
        """Should handle empty genome_hit_counts."""
        summary = TaxonomicSummary(
            total_reads=100,
            known_species=50,
            novel_species=25,
            novel_genus=10,
            conserved_regions=15,
            mean_novelty_index=3.0,
            mean_placement_uncertainty=0.8,
        )

        assert summary.genome_hit_counts == {}

    def test_to_json(self, sample_summary, temp_dir):
        """to_json should write valid JSON file."""
        output_path = temp_dir / "summary.json"
        sample_summary.to_json(output_path)

        assert output_path.exists()
        content = output_path.read_text()
        assert "total_reads" in content
        assert "1000" in content

    def test_from_json(self, sample_summary, temp_dir):
        """from_json should load summary from JSON file."""
        output_path = temp_dir / "summary.json"
        sample_summary.to_json(output_path)

        loaded = TaxonomicSummary.from_json(output_path)

        assert loaded.total_reads == sample_summary.total_reads
        assert loaded.known_species == sample_summary.known_species
        assert loaded.novel_species == sample_summary.novel_species
        assert loaded.novel_genus == sample_summary.novel_genus
        assert loaded.mean_novelty_index == sample_summary.mean_novelty_index

    def test_json_roundtrip_preserves_data(self, sample_summary, temp_dir):
        """JSON serialization should preserve all data."""
        output_path = temp_dir / "summary.json"
        sample_summary.to_json(output_path)
        loaded = TaxonomicSummary.from_json(output_path)

        # Computed fields should match
        assert loaded.known_species_pct == sample_summary.known_species_pct
        assert loaded.novel_diversity_pct == sample_summary.novel_diversity_pct
        assert loaded.genome_hit_counts == sample_summary.genome_hit_counts

    def test_genome_hit_counts_preserved(self, sample_summary, temp_dir):
        """genome_hit_counts should survive JSON roundtrip."""
        output_path = temp_dir / "summary.json"
        sample_summary.to_json(output_path)
        loaded = TaxonomicSummary.from_json(output_path)

        assert "GCF_000123456.1" in loaded.genome_hit_counts
        assert loaded.genome_hit_counts["GCF_000123456.1"] == 300


class TestClassificationScenarios:
    """Tests for realistic classification scenarios."""

    def test_high_confidence_known_species(self):
        """High identity, low uncertainty should be known species."""
        classification = ReadClassification(
            read_id="read_high_conf",
            best_match_genome="GCF_000123456.1",
            top_hit_identity=99.5,
            novelty_index=0.5,
            placement_uncertainty=0.1,
            num_ambiguous_hits=1,
            taxonomic_call=TaxonomicCall.KNOWN_SPECIES,
        )

        assert classification.is_novel is False
        assert classification.taxonomic_call == TaxonomicCall.KNOWN_SPECIES

    def test_novel_species_boundary(self):
        """Novelty at species boundary should classify correctly."""
        # Novelty = 5.0 is minimum for novel species
        classification = ReadClassification(
            read_id="read_boundary",
            best_match_genome="GCF_000789012.1",
            top_hit_identity=95.0,  # Novelty = 5.0
            novelty_index=5.0,
            placement_uncertainty=0.3,
            num_ambiguous_hits=2,
            taxonomic_call=TaxonomicCall.NOVEL_SPECIES,
        )

        assert classification.is_novel is True

    def test_novel_genus_boundary(self):
        """Novelty at genus boundary should classify correctly."""
        # Novelty = 15.0 is minimum for novel genus
        classification = ReadClassification(
            read_id="read_genus",
            best_match_genome="GCA_000111222.1",
            top_hit_identity=85.0,  # Novelty = 15.0
            novelty_index=15.0,
            placement_uncertainty=1.8,
            num_ambiguous_hits=4,
            taxonomic_call=TaxonomicCall.NOVEL_GENUS,
        )

        assert classification.is_novel is True

    def test_conserved_region_high_uncertainty(self):
        """High uncertainty should indicate conserved region."""
        classification = ReadClassification(
            read_id="read_conserved",
            best_match_genome="GCF_000123456.1",
            top_hit_identity=95.0,
            novelty_index=5.0,
            placement_uncertainty=15.0,  # High uncertainty
            num_ambiguous_hits=20,
            taxonomic_call=TaxonomicCall.CONSERVED_REGION,
        )

        assert classification.is_novel is False
        assert classification.taxonomic_call == TaxonomicCall.CONSERVED_REGION
