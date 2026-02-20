"""
Tests for protein-level read classification.

Tests the protein constants, alignment mode in ScoringConfig, and the
integration with the classifier for protein-mode classification.
"""

from __future__ import annotations

import pytest

from metadarkmatter.core.constants import (
    CONFIDENCE_IDENTITY_GAP_THRESHOLDS,
    CONFIDENCE_IDENTITY_SCORE_BASE,
    CONFIDENCE_IDENTITY_SCORE_RANGE,
    CONFIDENCE_MARGIN_DIVISOR_KNOWN,
    CONFIDENCE_MARGIN_DIVISOR_NOVEL_GENUS,
    CONFIDENCE_MARGIN_DIVISOR_NOVEL_SPECIES,
    calculate_confidence_score,
)
from metadarkmatter.core.protein_constants import (
    PROTEIN_CONFIDENCE_IDENTITY_GAP_THRESHOLDS,
    PROTEIN_CONFIDENCE_IDENTITY_SCORE_BASE,
    PROTEIN_CONFIDENCE_IDENTITY_SCORE_RANGE,
    PROTEIN_CONFIDENCE_MARGIN_DIVISOR_KNOWN,
    PROTEIN_CONFIDENCE_MARGIN_DIVISOR_NOVEL_GENUS,
    PROTEIN_CONFIDENCE_MARGIN_DIVISOR_NOVEL_SPECIES,
    PROTEIN_NOVELTY_KNOWN_MAX,
    PROTEIN_NOVELTY_NOVEL_GENUS_MAX,
    PROTEIN_NOVELTY_NOVEL_GENUS_MIN,
    PROTEIN_NOVELTY_NOVEL_SPECIES_MAX,
    PROTEIN_NOVELTY_NOVEL_SPECIES_MIN,
    PROTEIN_UNCERTAINTY_CONFIDENT_MAX,
    PROTEIN_UNCERTAINTY_CONSERVED_MIN,
)
from metadarkmatter.models.config import ScoringConfig


def _calculate_protein_confidence_score(
    novelty_index: float,
    placement_uncertainty: float,
    num_ambiguous_hits: int,
    identity_gap: float | None,
    top_hit_identity: float,
    taxonomic_call: str,
    novelty_known_max: float = PROTEIN_NOVELTY_KNOWN_MAX,
    novelty_novel_species_max: float = PROTEIN_NOVELTY_NOVEL_SPECIES_MAX,
    novelty_novel_genus_max: float = PROTEIN_NOVELTY_NOVEL_GENUS_MAX,
    uncertainty_confident_max: float = PROTEIN_UNCERTAINTY_CONFIDENT_MAX,
) -> float:
    """Test helper: calculate confidence score with protein-specific parameters."""
    return calculate_confidence_score(
        novelty_index=novelty_index,
        placement_uncertainty=placement_uncertainty,
        num_ambiguous_hits=num_ambiguous_hits,
        identity_gap=identity_gap,
        top_hit_identity=top_hit_identity,
        taxonomic_call=taxonomic_call,
        novelty_known_max=novelty_known_max,
        novelty_novel_species_max=novelty_novel_species_max,
        novelty_novel_genus_max=novelty_novel_genus_max,
        uncertainty_confident_max=uncertainty_confident_max,
        margin_divisor_known=PROTEIN_CONFIDENCE_MARGIN_DIVISOR_KNOWN,
        margin_divisor_novel_species=PROTEIN_CONFIDENCE_MARGIN_DIVISOR_NOVEL_SPECIES,
        margin_divisor_novel_genus=PROTEIN_CONFIDENCE_MARGIN_DIVISOR_NOVEL_GENUS,
        identity_gap_thresholds=PROTEIN_CONFIDENCE_IDENTITY_GAP_THRESHOLDS,
        identity_score_base=PROTEIN_CONFIDENCE_IDENTITY_SCORE_BASE,
        identity_score_range=PROTEIN_CONFIDENCE_IDENTITY_SCORE_RANGE,
    )


class TestProteinConstants:
    """Test protein threshold constants."""

    def test_protein_novelty_thresholds_are_wider_than_nucleotide(self) -> None:
        """Protein thresholds should be wider than nucleotide thresholds."""
        # Protein Known Species: N < 10% (vs N < 5% for nucleotide)
        assert PROTEIN_NOVELTY_KNOWN_MAX == 10.0

        # Protein Novel Species: 10-25% (vs 5-20% for nucleotide)
        assert PROTEIN_NOVELTY_NOVEL_SPECIES_MIN == 10.0
        assert PROTEIN_NOVELTY_NOVEL_SPECIES_MAX == 25.0

        # Protein Novel Genus: 25-40% (vs 20-25% for nucleotide)
        assert PROTEIN_NOVELTY_NOVEL_GENUS_MIN == 25.0
        assert PROTEIN_NOVELTY_NOVEL_GENUS_MAX == 40.0

    def test_protein_uncertainty_thresholds(self) -> None:
        """Test protein uncertainty thresholds."""
        # Protein confident: U < 5% (vs U < 2% for nucleotide)
        assert PROTEIN_UNCERTAINTY_CONFIDENT_MAX == 5.0

        # Protein conserved region: U >= 10%
        assert PROTEIN_UNCERTAINTY_CONSERVED_MIN == 10.0

    def test_threshold_boundaries_are_continuous(self) -> None:
        """Threshold boundaries should be continuous (no gaps)."""
        # Known -> Novel Species boundary
        assert PROTEIN_NOVELTY_KNOWN_MAX == PROTEIN_NOVELTY_NOVEL_SPECIES_MIN

        # Novel Species -> Novel Genus boundary
        assert PROTEIN_NOVELTY_NOVEL_SPECIES_MAX == PROTEIN_NOVELTY_NOVEL_GENUS_MIN


class TestScoringConfigAlignmentMode:
    """Test alignment mode in ScoringConfig."""

    def test_default_alignment_mode_is_nucleotide(self) -> None:
        """Default alignment mode should be nucleotide."""
        config = ScoringConfig()
        assert config.alignment_mode == "nucleotide"

    def test_can_create_protein_mode_config(self) -> None:
        """Should be able to create protein mode config."""
        config = ScoringConfig(alignment_mode="protein")
        assert config.alignment_mode == "protein"

    def test_for_protein_mode_factory(self) -> None:
        """Factory method should create protein mode config."""
        config = ScoringConfig.for_protein_mode()
        assert config.alignment_mode == "protein"

    def test_get_effective_thresholds_nucleotide(self) -> None:
        """Nucleotide mode should return standard thresholds."""
        config = ScoringConfig(alignment_mode="nucleotide")
        eff = config.get_effective_thresholds()

        # Should match config's own thresholds
        assert eff["novelty_known_max"] == config.novelty_known_max
        assert eff["novelty_novel_species_max"] == config.novelty_novel_species_max
        assert eff["novelty_novel_genus_max"] == config.novelty_novel_genus_max
        assert eff["uncertainty_known_max"] == config.uncertainty_known_max

    def test_get_effective_thresholds_protein(self) -> None:
        """Protein mode should return protein-specific thresholds."""
        config = ScoringConfig(alignment_mode="protein")
        eff = config.get_effective_thresholds()

        # Should match protein constants
        assert eff["novelty_known_max"] == PROTEIN_NOVELTY_KNOWN_MAX
        assert eff["novelty_novel_species_min"] == PROTEIN_NOVELTY_NOVEL_SPECIES_MIN
        assert eff["novelty_novel_species_max"] == PROTEIN_NOVELTY_NOVEL_SPECIES_MAX
        assert eff["novelty_novel_genus_min"] == PROTEIN_NOVELTY_NOVEL_GENUS_MIN
        assert eff["novelty_novel_genus_max"] == PROTEIN_NOVELTY_NOVEL_GENUS_MAX
        assert eff["uncertainty_known_max"] == PROTEIN_UNCERTAINTY_CONFIDENT_MAX
        assert eff["uncertainty_conserved_min"] == PROTEIN_UNCERTAINTY_CONSERVED_MIN

    def test_protein_mode_ignores_explicit_novelty_thresholds(self) -> None:
        """Protein mode should use protein constants even if config has custom values."""
        config = ScoringConfig(
            alignment_mode="protein",
            novelty_known_max=3.0,  # Custom value, should be ignored
        )
        eff = config.get_effective_thresholds()

        # Should use protein constant, not config value
        assert eff["novelty_known_max"] == PROTEIN_NOVELTY_KNOWN_MAX
        assert eff["novelty_known_max"] != config.novelty_known_max


class TestProteinConfidenceScore:
    """Test protein confidence score calculation."""

    def test_known_species_high_confidence(self) -> None:
        """Known species with clear margin should have high confidence."""
        score = _calculate_protein_confidence_score(
            novelty_index=3.0,  # Well below 10% threshold
            placement_uncertainty=1.0,  # Low uncertainty
            num_ambiguous_hits=1,
            identity_gap=15.0,  # Large gap
            top_hit_identity=97.0,
            taxonomic_call="Known Species",
        )
        assert score >= 70.0

    def test_novel_species_moderate_confidence(self) -> None:
        """Novel species in middle of range should have moderate confidence."""
        score = _calculate_protein_confidence_score(
            novelty_index=17.5,  # Middle of 10-25% range
            placement_uncertainty=2.0,
            num_ambiguous_hits=2,
            identity_gap=8.0,
            top_hit_identity=82.5,
            taxonomic_call="Novel Species",
        )
        assert 40.0 <= score <= 80.0

    def test_novel_genus_moderate_confidence(self) -> None:
        """Novel genus classification should have moderate confidence."""
        score = _calculate_protein_confidence_score(
            novelty_index=32.5,  # Middle of 25-40% range
            placement_uncertainty=3.0,
            num_ambiguous_hits=3,
            identity_gap=5.0,
            top_hit_identity=67.5,
            taxonomic_call="Novel Genus",
        )
        assert 30.0 <= score <= 70.0

    def test_ambiguous_low_confidence(self) -> None:
        """Ambiguous classifications should have low confidence."""
        score = _calculate_protein_confidence_score(
            novelty_index=15.0,
            placement_uncertainty=8.0,  # High uncertainty
            num_ambiguous_hits=10,
            identity_gap=0.5,  # Small gap
            top_hit_identity=85.0,
            taxonomic_call="Ambiguous",
        )
        assert score <= 40.0


class TestProteinVsNucleotideClassification:
    """Test that protein mode classifies differently than nucleotide mode."""

    def test_same_identity_different_classification(self) -> None:
        """Same identity value should classify differently based on mode."""
        # 92% identity: N = 8%
        # Nucleotide: N > 5% -> Novel Species
        # Protein: N < 10% -> Known Species
        nuc_config = ScoringConfig(alignment_mode="nucleotide")
        prot_config = ScoringConfig(alignment_mode="protein")

        nuc_eff = nuc_config.get_effective_thresholds()
        prot_eff = prot_config.get_effective_thresholds()

        novelty = 8.0  # 92% identity

        # In nucleotide mode, this is above known species threshold
        assert novelty >= nuc_eff["novelty_known_max"]

        # In protein mode, this is below known species threshold
        assert novelty < prot_eff["novelty_known_max"]

    def test_protein_novel_species_thresholds_shifted_higher(self) -> None:
        """Protein mode should have higher novel species thresholds.

        Protein identity is more conserved than nucleotide identity,
        so the thresholds are shifted to account for slower divergence.
        """
        nuc_config = ScoringConfig(alignment_mode="nucleotide")
        prot_config = ScoringConfig(alignment_mode="protein")

        nuc_eff = nuc_config.get_effective_thresholds()
        prot_eff = prot_config.get_effective_thresholds()

        # Protein thresholds should be shifted higher
        assert prot_eff["novelty_novel_species_min"] > nuc_eff["novelty_novel_species_min"]
        assert prot_eff["novelty_novel_species_max"] > nuc_eff["novelty_novel_species_max"]

    def test_wider_novel_genus_range_in_protein(self) -> None:
        """Protein mode should have wider novel genus range."""
        nuc_config = ScoringConfig(alignment_mode="nucleotide")
        prot_config = ScoringConfig(alignment_mode="protein")

        nuc_eff = nuc_config.get_effective_thresholds()
        prot_eff = prot_config.get_effective_thresholds()

        nuc_range = nuc_eff["novelty_novel_genus_max"] - nuc_eff["novelty_novel_genus_min"]
        prot_range = prot_eff["novelty_novel_genus_max"] - prot_eff["novelty_novel_genus_min"]

        # Protein range should be wider
        assert prot_range > nuc_range


class TestConfidenceScoreScalingParameters:
    """Test confidence score scaling parameters in effective thresholds."""

    def test_nucleotide_mode_includes_scaling_parameters(self) -> None:
        """Nucleotide mode should include confidence score scaling parameters."""
        config = ScoringConfig(alignment_mode="nucleotide")
        eff = config.get_effective_thresholds()

        # Should have all scaling parameters
        assert "margin_divisor_known" in eff
        assert "margin_divisor_novel_species" in eff
        assert "margin_divisor_novel_genus" in eff
        assert "identity_gap_thresholds" in eff
        assert "identity_score_base" in eff
        assert "identity_score_range" in eff

        # Should match nucleotide constants
        assert eff["margin_divisor_known"] == CONFIDENCE_MARGIN_DIVISOR_KNOWN
        assert eff["margin_divisor_novel_species"] == CONFIDENCE_MARGIN_DIVISOR_NOVEL_SPECIES
        assert eff["margin_divisor_novel_genus"] == CONFIDENCE_MARGIN_DIVISOR_NOVEL_GENUS
        assert eff["identity_gap_thresholds"] == CONFIDENCE_IDENTITY_GAP_THRESHOLDS
        assert eff["identity_score_base"] == CONFIDENCE_IDENTITY_SCORE_BASE
        assert eff["identity_score_range"] == CONFIDENCE_IDENTITY_SCORE_RANGE

    def test_protein_mode_includes_scaling_parameters(self) -> None:
        """Protein mode should include protein-specific scaling parameters."""
        config = ScoringConfig(alignment_mode="protein")
        eff = config.get_effective_thresholds()

        # Should have all scaling parameters
        assert "margin_divisor_known" in eff
        assert "margin_divisor_novel_species" in eff
        assert "margin_divisor_novel_genus" in eff
        assert "identity_gap_thresholds" in eff
        assert "identity_score_base" in eff
        assert "identity_score_range" in eff

        # Should match protein constants
        assert eff["margin_divisor_known"] == PROTEIN_CONFIDENCE_MARGIN_DIVISOR_KNOWN
        assert eff["margin_divisor_novel_species"] == PROTEIN_CONFIDENCE_MARGIN_DIVISOR_NOVEL_SPECIES
        assert eff["margin_divisor_novel_genus"] == PROTEIN_CONFIDENCE_MARGIN_DIVISOR_NOVEL_GENUS
        assert eff["identity_gap_thresholds"] == PROTEIN_CONFIDENCE_IDENTITY_GAP_THRESHOLDS
        assert eff["identity_score_base"] == PROTEIN_CONFIDENCE_IDENTITY_SCORE_BASE
        assert eff["identity_score_range"] == PROTEIN_CONFIDENCE_IDENTITY_SCORE_RANGE

    def test_protein_scaling_parameters_are_wider(self) -> None:
        """Protein scaling parameters should be wider than nucleotide."""
        nuc_config = ScoringConfig(alignment_mode="nucleotide")
        prot_config = ScoringConfig(alignment_mode="protein")

        nuc_eff = nuc_config.get_effective_thresholds()
        prot_eff = prot_config.get_effective_thresholds()

        # Protein margin divisors should be larger (wider scaling)
        assert prot_eff["margin_divisor_known"] > nuc_eff["margin_divisor_known"]
        assert prot_eff["margin_divisor_novel_species"] > nuc_eff["margin_divisor_novel_species"]
        assert prot_eff["margin_divisor_novel_genus"] > nuc_eff["margin_divisor_novel_genus"]

        # Protein identity gap thresholds should be higher
        assert prot_eff["identity_gap_thresholds"][0] > nuc_eff["identity_gap_thresholds"][0]

        # Protein identity score base should be lower (accommodate lower identities)
        assert prot_eff["identity_score_base"] < nuc_eff["identity_score_base"]

        # Protein identity score range should be wider
        assert prot_eff["identity_score_range"] > nuc_eff["identity_score_range"]


class TestConfidenceScoreConsistency:
    """Test that confidence scores are consistent between direct call and delegated call."""

    def test_protein_confidence_score_delegates_correctly(self) -> None:
        """Protein confidence score should delegate to parameterized nucleotide function."""
        # Calculate using protein-specific function
        protein_score = _calculate_protein_confidence_score(
            novelty_index=5.0,
            placement_uncertainty=2.0,
            num_ambiguous_hits=2,
            identity_gap=8.0,
            top_hit_identity=95.0,
            taxonomic_call="Known Species",
        )

        # Calculate using parameterized base function with protein parameters
        direct_score = calculate_confidence_score(
            novelty_index=5.0,
            placement_uncertainty=2.0,
            num_ambiguous_hits=2,
            identity_gap=8.0,
            top_hit_identity=95.0,
            taxonomic_call="Known Species",
            novelty_known_max=PROTEIN_NOVELTY_KNOWN_MAX,
            novelty_novel_species_max=PROTEIN_NOVELTY_NOVEL_SPECIES_MAX,
            novelty_novel_genus_max=PROTEIN_NOVELTY_NOVEL_GENUS_MAX,
            uncertainty_confident_max=PROTEIN_UNCERTAINTY_CONFIDENT_MAX,
            margin_divisor_known=PROTEIN_CONFIDENCE_MARGIN_DIVISOR_KNOWN,
            margin_divisor_novel_species=PROTEIN_CONFIDENCE_MARGIN_DIVISOR_NOVEL_SPECIES,
            margin_divisor_novel_genus=PROTEIN_CONFIDENCE_MARGIN_DIVISOR_NOVEL_GENUS,
            identity_gap_thresholds=PROTEIN_CONFIDENCE_IDENTITY_GAP_THRESHOLDS,
            identity_score_base=PROTEIN_CONFIDENCE_IDENTITY_SCORE_BASE,
            identity_score_range=PROTEIN_CONFIDENCE_IDENTITY_SCORE_RANGE,
        )

        # Should produce identical results
        assert protein_score == direct_score

    def test_nucleotide_and_protein_produce_different_scores(self) -> None:
        """Same inputs with different modes should produce different confidence scores."""
        # Use inputs that would produce different scores due to scaling
        inputs = {
            "novelty_index": 8.0,  # Known species in protein, novel species in nucleotide
            "placement_uncertainty": 3.0,
            "num_ambiguous_hits": 3,
            "identity_gap": 3.0,  # Above protein threshold (2.0), below nucleotide threshold (5.0)
            "top_hit_identity": 92.0,
            "taxonomic_call": "Known Species",
        }

        # Nucleotide mode with nucleotide defaults
        nuc_score = calculate_confidence_score(**inputs)

        # Protein mode with protein scaling
        prot_score = _calculate_protein_confidence_score(**inputs)

        # Scores should be different due to different scaling parameters
        # Note: The exact difference depends on the specific inputs and thresholds
        # but they should not be identical
        assert abs(nuc_score - prot_score) > 0.1
