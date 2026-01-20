"""
Unit tests for AAI matrix and genus-level classification.

Tests AAIMatrix, AAI-based classification thresholds, and integration
with the ANI-weighted classifier for genus-level decisions.

AAI (Average Amino Acid Identity) is used for genus-level classification
because ANI becomes unreliable below approximately 80% identity.

AAI genus boundaries (Riesco & Trujillo 2024):
- Same genus: AAI > 65%
- Genus boundary: AAI 58-65%
- Different genus: AAI < 58%
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metadarkmatter.core.aai_matrix_builder import (
    AAIMatrix,
    extract_genome_from_protein_header,
    extract_genome_name_from_filename,
)
from metadarkmatter.core.constants import (
    AAI_DEFAULT_UNRELATED,
    AAI_GENUS_BOUNDARY_HIGH,
    AAI_GENUS_BOUNDARY_LOW,
)


# =============================================================================
# AAI Matrix Test Data Fixtures
# =============================================================================


@pytest.fixture
def small_aai_dict() -> dict[str, dict[str, float]]:
    """Small AAI matrix as nested dictionary for testing."""
    return {
        "GCF_000123456.1": {
            "GCF_000123456.1": 100.0,
            "GCF_000789012.1": 72.5,  # Same genus (AAI > 65%)
            "GCA_000111222.1": 55.0,  # Different genus (AAI < 58%)
        },
        "GCF_000789012.1": {
            "GCF_000123456.1": 72.5,
            "GCF_000789012.1": 100.0,
            "GCA_000111222.1": 52.0,  # Different genus
        },
        "GCA_000111222.1": {
            "GCF_000123456.1": 55.0,
            "GCF_000789012.1": 52.0,
            "GCA_000111222.1": 100.0,
        },
    }


@pytest.fixture
def genus_boundary_aai_dict() -> dict[str, dict[str, float]]:
    """AAI matrix with values in the genus boundary zone (58-65%)."""
    return {
        "GCF_000001.1": {
            "GCF_000001.1": 100.0,
            "GCF_000002.1": 62.0,  # In boundary zone
            "GCF_000003.1": 59.0,  # In boundary zone, lower end
        },
        "GCF_000002.1": {
            "GCF_000001.1": 62.0,
            "GCF_000002.1": 100.0,
            "GCF_000003.1": 70.0,  # Same genus
        },
        "GCF_000003.1": {
            "GCF_000001.1": 59.0,
            "GCF_000002.1": 70.0,
            "GCF_000003.1": 100.0,
        },
    }


@pytest.fixture
def temp_aai_file(tmp_path: Path, small_aai_dict: dict) -> Path:
    """Create temporary AAI matrix CSV file."""
    aai_file = tmp_path / "aai_matrix.csv"
    genomes = sorted(small_aai_dict.keys())

    # Write CSV with genome column and genome columns
    lines = ["genome," + ",".join(genomes)]
    for g1 in genomes:
        row = [g1] + [str(small_aai_dict[g1].get(g2, 0.0)) for g2 in genomes]
        lines.append(",".join(row))

    aai_file.write_text("\n".join(lines))
    return aai_file


# =============================================================================
# AAIMatrix Tests
# =============================================================================


class TestAAIMatrix:
    """Tests for AAIMatrix class."""

    def test_create_from_dict(self, small_aai_dict):
        """Should create AAIMatrix from nested dictionary."""
        matrix = AAIMatrix(small_aai_dict)

        assert len(matrix) == 3

    def test_genomes_property(self, small_aai_dict):
        """genomes property should return set of genome names."""
        matrix = AAIMatrix(small_aai_dict)
        genomes = matrix.genomes

        assert isinstance(genomes, set)
        assert "GCF_000123456.1" in genomes
        assert "GCF_000789012.1" in genomes
        assert "GCA_000111222.1" in genomes

    def test_get_aai_same_genome(self, small_aai_dict):
        """get_aai with same genome should return 100.0."""
        matrix = AAIMatrix(small_aai_dict)

        aai = matrix.get_aai("GCF_000123456.1", "GCF_000123456.1")

        assert aai == 100.0

    def test_get_aai_same_genus(self, small_aai_dict):
        """get_aai should return correct value for same genus pair."""
        matrix = AAIMatrix(small_aai_dict)

        aai = matrix.get_aai("GCF_000123456.1", "GCF_000789012.1")

        assert aai == 72.5
        assert aai > AAI_GENUS_BOUNDARY_HIGH  # Same genus

    def test_get_aai_different_genus(self, small_aai_dict):
        """get_aai should return correct value for different genus pair."""
        matrix = AAIMatrix(small_aai_dict)

        aai = matrix.get_aai("GCF_000123456.1", "GCA_000111222.1")

        assert aai == 55.0
        assert aai < AAI_GENUS_BOUNDARY_LOW  # Different genus

    def test_get_aai_symmetric(self, small_aai_dict):
        """get_aai should be symmetric."""
        matrix = AAIMatrix(small_aai_dict)

        aai_ab = matrix.get_aai("GCF_000123456.1", "GCF_000789012.1")
        aai_ba = matrix.get_aai("GCF_000789012.1", "GCF_000123456.1")

        assert aai_ab == aai_ba

    def test_get_aai_unknown_genome(self, small_aai_dict):
        """get_aai with unknown genome should return 0.0."""
        matrix = AAIMatrix(small_aai_dict)

        aai = matrix.get_aai("GCF_000123456.1", "unknown_genome")

        assert aai == 0.0

    def test_get_aai_both_unknown(self, small_aai_dict):
        """get_aai with both genomes unknown should return 0.0."""
        matrix = AAIMatrix(small_aai_dict)

        aai = matrix.get_aai("unknown_1", "unknown_2")

        assert aai == 0.0

    def test_has_genome_true(self, small_aai_dict):
        """has_genome should return True for known genome."""
        matrix = AAIMatrix(small_aai_dict)

        assert matrix.has_genome("GCF_000123456.1") is True

    def test_has_genome_false(self, small_aai_dict):
        """has_genome should return False for unknown genome."""
        matrix = AAIMatrix(small_aai_dict)

        assert matrix.has_genome("unknown_genome") is False

    def test_get_genome_idx(self, small_aai_dict):
        """get_genome_idx should return integer index."""
        matrix = AAIMatrix(small_aai_dict)

        idx = matrix.get_genome_idx("GCF_000123456.1")

        assert isinstance(idx, int)
        assert 0 <= idx < len(matrix)

    def test_get_genome_idx_unknown(self, small_aai_dict):
        """get_genome_idx should return None for unknown genome."""
        matrix = AAIMatrix(small_aai_dict)

        idx = matrix.get_genome_idx("unknown_genome")

        assert idx is None

    def test_get_aai_by_idx(self, small_aai_dict):
        """get_aai_by_idx should work with integer indices."""
        matrix = AAIMatrix(small_aai_dict)

        idx1 = matrix.get_genome_idx("GCF_000123456.1")
        idx2 = matrix.get_genome_idx("GCF_000789012.1")

        aai = matrix.get_aai_by_idx(idx1, idx2)

        assert aai == 72.5

    def test_from_file(self, temp_aai_file):
        """from_file should load AAI matrix from CSV."""
        matrix = AAIMatrix.from_file(temp_aai_file)

        assert len(matrix) == 3
        assert matrix.has_genome("GCF_000123456.1")

    def test_memory_usage_bytes(self, small_aai_dict):
        """memory_usage_bytes should return positive integer."""
        matrix = AAIMatrix(small_aai_dict)

        mem = matrix.memory_usage_bytes()

        assert isinstance(mem, int)
        assert mem > 0

    def test_default_aai_for_missing_values(self, small_aai_dict):
        """Missing AAI values (0.0) should return default_aai."""
        # Create matrix with a missing entry
        sparse_dict = {
            "GCF_000001.1": {
                "GCF_000001.1": 100.0,
                "GCF_000002.1": 0.0,  # Missing value
            },
            "GCF_000002.1": {
                "GCF_000001.1": 0.0,  # Missing value
                "GCF_000002.1": 100.0,
            },
        }
        matrix = AAIMatrix(sparse_dict, default_aai=50.0)

        # Should return default_aai instead of 0.0
        aai = matrix.get_aai("GCF_000001.1", "GCF_000002.1")
        assert aai == 50.0


# =============================================================================
# AAI Threshold Tests
# =============================================================================


class TestAAIThresholds:
    """Tests for AAI-based genus classification thresholds."""

    def test_same_genus_threshold(self, small_aai_dict):
        """AAI > 65% should indicate same genus."""
        matrix = AAIMatrix(small_aai_dict)

        aai = matrix.get_aai("GCF_000123456.1", "GCF_000789012.1")

        # 72.5% > 65% = same genus
        assert aai > AAI_GENUS_BOUNDARY_HIGH

    def test_different_genus_threshold(self, small_aai_dict):
        """AAI < 58% should indicate different genus."""
        matrix = AAIMatrix(small_aai_dict)

        aai = matrix.get_aai("GCF_000123456.1", "GCA_000111222.1")

        # 55.0% < 58% = different genus
        assert aai < AAI_GENUS_BOUNDARY_LOW

    def test_boundary_zone_identification(self, genus_boundary_aai_dict):
        """AAI 58-65% should be identified as genus boundary zone."""
        matrix = AAIMatrix(genus_boundary_aai_dict)

        aai = matrix.get_aai("GCF_000001.1", "GCF_000002.1")

        # 62.0% is in boundary zone
        assert AAI_GENUS_BOUNDARY_LOW <= aai <= AAI_GENUS_BOUNDARY_HIGH

    def test_uncertainty_calculation_same_genus(self, small_aai_dict):
        """AAI uncertainty for same genus should be < 35%."""
        matrix = AAIMatrix(small_aai_dict)

        aai = matrix.get_aai("GCF_000123456.1", "GCF_000789012.1")
        uncertainty = 100.0 - aai

        # 72.5% AAI -> 27.5% uncertainty
        assert uncertainty < 35.0  # Same genus threshold

    def test_uncertainty_calculation_different_genus(self, small_aai_dict):
        """AAI uncertainty for different genus should be > 42%."""
        matrix = AAIMatrix(small_aai_dict)

        aai = matrix.get_aai("GCF_000123456.1", "GCA_000111222.1")
        uncertainty = 100.0 - aai

        # 55.0% AAI -> 45.0% uncertainty
        assert uncertainty > 42.0  # Different genus threshold


# =============================================================================
# Helper Function Tests
# =============================================================================


class TestHelperFunctions:
    """Tests for AAI helper functions."""

    def test_extract_genome_from_protein_header_pipe_format(self):
        """Should extract genome from pipe-separated format."""
        header = ">GCF_000123456.1|WP_012345678.1"
        genome = extract_genome_from_protein_header(header)
        assert genome == "GCF_000123456.1"

    def test_extract_genome_from_protein_header_prodigal_format(self):
        """Should extract genome from Prodigal format."""
        header = ">GCF_000123456.1_1 # 1 # 500 # 1 # ID=1_1"
        genome = extract_genome_from_protein_header(header)
        assert genome == "GCF_000123456.1"

    def test_extract_genome_from_protein_header_no_prefix(self):
        """Should handle header without '>' prefix."""
        header = "GCF_000123456.1|protein_id"
        genome = extract_genome_from_protein_header(header)
        assert genome == "GCF_000123456.1"

    def test_extract_genome_name_from_filename_full_ncbi(self):
        """Should extract genome from full NCBI filename."""
        path = Path("/data/GCF_000123456.1_ASM123v1_protein.faa")
        genome = extract_genome_name_from_filename(path)
        assert genome == "GCF_000123456.1"

    def test_extract_genome_name_from_filename_simple(self):
        """Should extract genome from simple filename."""
        path = Path("/data/GCF_000123456.1.faa")
        genome = extract_genome_name_from_filename(path)
        assert genome == "GCF_000123456.1"

    def test_extract_genome_name_from_filename_gzipped(self):
        """Should handle gzipped filename."""
        path = Path("/data/GCF_000123456.1_protein.faa.gz")
        genome = extract_genome_name_from_filename(path)
        assert genome == "GCF_000123456.1"


# =============================================================================
# AAI Constants Tests
# =============================================================================


class TestAAIConstants:
    """Tests for AAI-related constants."""

    def test_genus_boundary_low(self):
        """Lower genus boundary should be 58%."""
        assert AAI_GENUS_BOUNDARY_LOW == 58.0

    def test_genus_boundary_high(self):
        """Upper genus boundary should be 65%."""
        assert AAI_GENUS_BOUNDARY_HIGH == 65.0

    def test_default_unrelated(self):
        """Default AAI for unrelated genomes should be 50%."""
        assert AAI_DEFAULT_UNRELATED == 50.0

    def test_boundary_zone_range(self):
        """Boundary zone should span 7% (58-65%)."""
        zone_size = AAI_GENUS_BOUNDARY_HIGH - AAI_GENUS_BOUNDARY_LOW
        assert zone_size == 7.0
