"""
Unit tests for configuration models.

Tests ScoringConfig, BlastConfig, Bowtie2Config, KrakenConfig, and GlobalConfig
including threshold validation and YAML serialization.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from pydantic import ValidationError

from metadarkmatter.models.config import (
    BlastConfig,
    Bowtie2Config,
    GlobalConfig,
    ScoringConfig,
)


class TestScoringConfig:
    """Tests for ScoringConfig model."""

    def test_default_values(self):
        """Default ScoringConfig should have expected values.

        Thresholds are based on literature:
        - 96% ANI = prokaryotic species boundary (Jain et al. 2018)
        - Novelty N = 100 - pident, so N < 4% = pident > 96% (same species)
        """
        config = ScoringConfig()

        assert config.bitscore_threshold_pct == 95.0
        # Novelty thresholds: continuous boundaries at species cutoff
        # Note: 20% novelty threshold accounts for read-level identity being
        # 10-20% lower than genome-level ANI
        assert config.novelty_known_max == 4.0  # pident > 96%
        assert config.novelty_novel_species_min == 4.0  # pident <= 96%
        assert config.novelty_novel_species_max == 20.0  # pident >= 80%
        assert config.novelty_novel_genus_min == 20.0  # pident < 80%
        assert config.novelty_novel_genus_max == 25.0
        # Uncertainty thresholds: based on competing genome ANI
        assert config.uncertainty_known_max == 1.5  # ANI > 98.5%
        assert config.uncertainty_novel_species_max == 1.5
        assert config.uncertainty_novel_genus_max == 1.5
        assert config.uncertainty_conserved_min == 5.0
        # Coverage weighting is on by default
        assert config.coverage_weight_mode == "linear"
        assert config.min_alignment_fraction == 0.3

    def test_custom_bitscore_threshold(self):
        """Should accept custom bitscore threshold."""
        config = ScoringConfig(bitscore_threshold_pct=90.0)

        assert config.bitscore_threshold_pct == 90.0

    def test_bitscore_threshold_max_100(self):
        """bitscore_threshold_pct should not exceed 100."""
        with pytest.raises(ValidationError) as exc_info:
            ScoringConfig(bitscore_threshold_pct=110.0)

        assert "bitscore_threshold_pct" in str(exc_info.value)

    def test_bitscore_threshold_min_0(self):
        """bitscore_threshold_pct should not be negative."""
        with pytest.raises(ValidationError) as exc_info:
            ScoringConfig(bitscore_threshold_pct=-5.0)

        assert "bitscore_threshold_pct" in str(exc_info.value)

    def test_novelty_species_min_less_than_known_max_rejected(self):
        """novelty_novel_species_min must be >= novelty_known_max.

        This ensures continuous boundaries: known < threshold <= novel.
        Equal values are allowed for gapless classification at the boundary.
        """
        with pytest.raises(ValidationError) as exc_info:
            ScoringConfig(
                novelty_known_max=6.0,
                novelty_novel_species_min=5.0,  # Must be >= 6.0
            )

        assert "novelty_novel_species_min" in str(exc_info.value)

    def test_novelty_genus_min_equals_species_max(self):
        """novelty_novel_genus_min must equal novelty_novel_species_max."""
        with pytest.raises(ValidationError) as exc_info:
            ScoringConfig(
                novelty_novel_species_max=15.0,
                novelty_novel_genus_min=16.0,  # Must be 15.0
            )

        assert "novelty_novel_genus_min" in str(exc_info.value)

    def test_uncertainty_species_max_le_genus_max(self):
        """uncertainty_novel_species_max must be <= uncertainty_novel_genus_max."""
        with pytest.raises(ValidationError) as exc_info:
            ScoringConfig(
                uncertainty_novel_species_max=3.0,
                uncertainty_novel_genus_max=2.0,  # Must be >= 3.0
            )

        assert "uncertainty_novel_species_max" in str(exc_info.value)

    def test_valid_custom_thresholds(self):
        """Valid custom thresholds should be accepted."""
        config = ScoringConfig(
            bitscore_threshold_pct=90.0,
            novelty_known_max=3.0,
            novelty_novel_species_min=6.0,
            novelty_novel_species_max=18.0,
            novelty_novel_genus_min=18.0,  # Must equal species_max
            novelty_novel_genus_max=30.0,
            uncertainty_known_max=1.0,
            uncertainty_novel_species_max=1.0,
            uncertainty_novel_genus_max=3.0,
            uncertainty_conserved_min=6.0,
        )

        assert config.novelty_known_max == 3.0
        assert config.novelty_novel_genus_min == 18.0

    def test_config_is_frozen(self):
        """ScoringConfig should be immutable."""
        config = ScoringConfig()

        with pytest.raises(ValidationError):
            config.bitscore_threshold_pct = 50.0

    def test_negative_novelty_rejected(self):
        """Negative novelty thresholds should be rejected."""
        with pytest.raises(ValidationError):
            ScoringConfig(novelty_known_max=-1.0)

    def test_negative_uncertainty_rejected(self):
        """Negative uncertainty thresholds should be rejected."""
        with pytest.raises(ValidationError):
            ScoringConfig(uncertainty_known_max=-0.5)

    def test_uncertainty_mode_default(self):
        """Default uncertainty_mode should be 'second'."""
        config = ScoringConfig()
        assert config.uncertainty_mode == "second"

    def test_uncertainty_mode_max(self):
        """Should accept 'max' uncertainty mode."""
        config = ScoringConfig(uncertainty_mode="max")
        assert config.uncertainty_mode == "max"

    def test_uncertainty_mode_second(self):
        """Should accept 'second' uncertainty mode."""
        config = ScoringConfig(uncertainty_mode="second")
        assert config.uncertainty_mode == "second"

    def test_uncertainty_mode_invalid_rejected(self):
        """Invalid uncertainty_mode should be rejected."""
        with pytest.raises(ValidationError):
            ScoringConfig(uncertainty_mode="invalid")


class TestBlastConfig:
    """Tests for BlastConfig model."""

    def test_default_values(self):
        """Default BlastConfig should have expected values."""
        config = BlastConfig()

        assert config.num_threads == 4
        assert config.max_target_seqs == 500
        assert config.evalue == 1e-5
        assert config.perc_identity == 0.0
        assert "6 qseqid sseqid" in config.outfmt

    def test_custom_threads(self):
        """Should accept custom thread count."""
        config = BlastConfig(num_threads=16)

        assert config.num_threads == 16

    def test_threads_min_1(self):
        """num_threads should be at least 1."""
        with pytest.raises(ValidationError):
            BlastConfig(num_threads=0)

    def test_max_target_seqs_min_1(self):
        """max_target_seqs should be at least 1."""
        with pytest.raises(ValidationError):
            BlastConfig(max_target_seqs=0)

    def test_evalue_positive(self):
        """evalue should be positive."""
        with pytest.raises(ValidationError):
            BlastConfig(evalue=0.0)

    def test_perc_identity_range(self):
        """perc_identity should be between 0 and 100."""
        config = BlastConfig(perc_identity=50.0)
        assert config.perc_identity == 50.0

        with pytest.raises(ValidationError):
            BlastConfig(perc_identity=150.0)

    def test_config_is_frozen(self):
        """BlastConfig should be immutable."""
        config = BlastConfig()

        with pytest.raises(ValidationError):
            config.num_threads = 8


class TestBowtie2Config:
    """Tests for Bowtie2Config model."""

    def test_default_values(self):
        """Default Bowtie2Config should have expected values."""
        config = Bowtie2Config()

        assert config.num_threads == 4
        assert config.mode == "local"
        assert config.very_sensitive is True
        assert config.max_alignments == 500

    def test_mode_local(self):
        """Should accept 'local' mode."""
        config = Bowtie2Config(mode="local")
        assert config.mode == "local"

    def test_mode_end_to_end(self):
        """Should accept 'end-to-end' mode."""
        config = Bowtie2Config(mode="end-to-end")
        assert config.mode == "end-to-end"

    def test_invalid_mode_rejected(self):
        """Invalid mode should be rejected."""
        with pytest.raises(ValidationError):
            Bowtie2Config(mode="invalid")

    def test_max_alignments_min_1(self):
        """max_alignments should be at least 1."""
        with pytest.raises(ValidationError):
            Bowtie2Config(max_alignments=0)

    def test_config_is_frozen(self):
        """Bowtie2Config should be immutable."""
        config = Bowtie2Config()

        with pytest.raises(ValidationError):
            config.mode = "end-to-end"


class TestGlobalConfig:
    """Tests for GlobalConfig settings model."""

    def test_default_values(self):
        """Default GlobalConfig should have expected values."""
        config = GlobalConfig()

        assert config.project_name == "metadarkmatter_analysis"
        assert config.num_threads == 4
        assert config.verbose is False

    def test_nested_scoring_config(self):
        """Should include nested ScoringConfig."""
        config = GlobalConfig()

        assert isinstance(config.scoring, ScoringConfig)
        assert config.scoring.bitscore_threshold_pct == 95.0

    def test_nested_blast_config(self):
        """Should include nested BlastConfig."""
        config = GlobalConfig()

        assert isinstance(config.blast, BlastConfig)
        assert config.blast.num_threads == 4

    def test_nested_bowtie2_config(self):
        """Should include nested Bowtie2Config."""
        config = GlobalConfig()

        assert isinstance(config.bowtie2, Bowtie2Config)
        assert config.bowtie2.mode == "local"

    def test_kraken_config_optional(self):
        """KrakenConfig should be optional (None by default)."""
        config = GlobalConfig()

        assert config.kraken is None

    def test_custom_project_name(self):
        """Should accept custom project name."""
        config = GlobalConfig(project_name="my_analysis")

        assert config.project_name == "my_analysis"

    def test_output_dir_created(self, temp_dir):
        """output_dir should be created if it doesn't exist."""
        new_dir = temp_dir / "new_output"
        config = GlobalConfig(output_dir=new_dir)

        assert new_dir.exists()
        assert config.output_dir == new_dir

    def test_to_yaml(self, temp_dir):
        """to_yaml should write valid YAML file."""
        config = GlobalConfig(
            project_name="test_project",
            output_dir=temp_dir / "output",
            num_threads=8,
        )
        yaml_path = temp_dir / "config.yaml"
        config.to_yaml(yaml_path)

        assert yaml_path.exists()
        content = yaml_path.read_text()
        assert "test_project" in content
        assert "num_threads: 8" in content

    def test_from_yaml(self, temp_dir):
        """from_yaml should load config from YAML file."""
        # Create a config and save it
        original = GlobalConfig(
            project_name="yaml_test",
            output_dir=temp_dir / "yaml_output",
            num_threads=12,
            verbose=True,
        )
        yaml_path = temp_dir / "config.yaml"
        original.to_yaml(yaml_path)

        # Load it back
        loaded = GlobalConfig.from_yaml(yaml_path)

        assert loaded.project_name == "yaml_test"
        assert loaded.num_threads == 12
        assert loaded.verbose is True

    def test_yaml_roundtrip_preserves_nested_config(self, temp_dir):
        """YAML roundtrip should preserve nested config values."""
        original = GlobalConfig(
            project_name="nested_test",
            output_dir=temp_dir / "nested_output",
        )
        yaml_path = temp_dir / "config.yaml"
        original.to_yaml(yaml_path)
        loaded = GlobalConfig.from_yaml(yaml_path)

        assert loaded.scoring.bitscore_threshold_pct == original.scoring.bitscore_threshold_pct
        assert loaded.blast.max_target_seqs == original.blast.max_target_seqs
        assert loaded.bowtie2.mode == original.bowtie2.mode


class TestScoringConfigBoundaries:
    """Edge case tests for ScoringConfig threshold boundaries."""

    def test_novelty_boundaries_are_continuous(self):
        """Novel species max should equal novel genus min."""
        config = ScoringConfig()

        # Ensure no gap in classification boundaries
        assert config.novelty_novel_species_max == config.novelty_novel_genus_min

    def test_known_species_and_novel_species_continuous(self):
        """Boundaries should be continuous for gapless classification.

        The defaults use novelty_known_max == novelty_novel_species_min (both 4.0)
        to ensure every read is classified without gaps.
        Classification logic: known if N < 4.0, novel species if N >= 4.0.
        """
        config = ScoringConfig()

        # Continuous boundary at 4.0 (96% identity = species boundary)
        gap = config.novelty_novel_species_min - config.novelty_known_max
        assert gap == 0.0  # Continuous, no gap

    def test_uncertainty_progression(self):
        """Uncertainty thresholds should increase: known < novel_species <= novel_genus < conserved."""
        config = ScoringConfig()

        assert config.uncertainty_known_max <= config.uncertainty_novel_species_max
        assert config.uncertainty_novel_species_max <= config.uncertainty_novel_genus_max
        assert config.uncertainty_novel_genus_max < config.uncertainty_conserved_min

    def test_extreme_but_valid_thresholds(self):
        """Extreme but valid threshold values should be accepted."""
        config = ScoringConfig(
            novelty_known_max=0.1,
            novelty_novel_species_min=0.2,
            novelty_novel_species_max=50.0,
            novelty_novel_genus_min=50.0,
            novelty_novel_genus_max=99.0,
        )

        assert config.novelty_known_max == 0.1
        assert config.novelty_novel_genus_max == 99.0

    def test_zero_thresholds_where_valid(self):
        """Zero thresholds should be valid where semantically meaningful."""
        # novelty_known_max = 0 means only 100% identity is "known"
        config = ScoringConfig(
            novelty_known_max=0.0,
            novelty_novel_species_min=1.0,
            novelty_novel_species_max=15.0,
            novelty_novel_genus_min=15.0,
        )

        assert config.novelty_known_max == 0.0
