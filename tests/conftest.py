"""
Shared pytest fixtures for metadarkmatter tests.

Provides reusable test data, temporary files, and mock objects
for unit and integration testing.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
import pytest


# =============================================================================
# BLAST Test Data Fixtures
# =============================================================================


@pytest.fixture
def valid_blast_line() -> str:
    """Single valid BLAST tabular output line."""
    return "read_001\tGCF_000123456.1_ASM123v1_genomic\t98.5\t150\t2\t0\t1\t150\t1000\t1150\t1e-50\t250.0"


@pytest.fixture
def valid_blast_hit_dict() -> dict[str, Any]:
    """Valid BLAST hit as dictionary."""
    return {
        "qseqid": "read_001",
        "sseqid": "GCF_000123456.1_ASM123v1_genomic",
        "pident": 98.5,
        "length": 150,
        "mismatch": 2,
        "gapopen": 0,
        "qstart": 1,
        "qend": 150,
        "sstart": 1000,
        "send": 1150,
        "evalue": 1e-50,
        "bitscore": 250.0,
    }


@pytest.fixture
def multi_hit_blast_data() -> list[dict[str, Any]]:
    """Multiple BLAST hits for a single read with varying bitscores."""
    return [
        {
            "qseqid": "read_001",
            "sseqid": "GCF_000123456.1_ASM123v1_genomic",
            "pident": 98.5,
            "length": 150,
            "mismatch": 2,
            "gapopen": 0,
            "qstart": 1,
            "qend": 150,
            "sstart": 1000,
            "send": 1150,
            "evalue": 1e-50,
            "bitscore": 250.0,
        },
        {
            "qseqid": "read_001",
            "sseqid": "GCF_000789012.1_ASM789v1_genomic",
            "pident": 95.0,
            "length": 150,
            "mismatch": 7,
            "gapopen": 0,
            "qstart": 1,
            "qend": 150,
            "sstart": 2000,
            "send": 2150,
            "evalue": 1e-45,
            "bitscore": 240.0,
        },
        {
            "qseqid": "read_001",
            "sseqid": "GCA_000111222.1_ASM111v1_genomic",
            "pident": 85.0,
            "length": 150,
            "mismatch": 22,
            "gapopen": 0,
            "qstart": 1,
            "qend": 150,
            "sstart": 3000,
            "send": 3150,
            "evalue": 1e-30,
            "bitscore": 180.0,
        },
    ]


# =============================================================================
# ANI Matrix Test Data Fixtures
# =============================================================================


@pytest.fixture
def small_ani_dict() -> dict[str, dict[str, float]]:
    """Small ANI matrix as nested dictionary (3 genomes)."""
    return {
        "GCF_000123456.1": {
            "GCF_000123456.1": 100.0,
            "GCF_000789012.1": 95.5,
            "GCA_000111222.1": 80.0,
        },
        "GCF_000789012.1": {
            "GCF_000123456.1": 95.5,
            "GCF_000789012.1": 100.0,
            "GCA_000111222.1": 82.0,
        },
        "GCA_000111222.1": {
            "GCF_000123456.1": 80.0,
            "GCF_000789012.1": 82.0,
            "GCA_000111222.1": 100.0,
        },
    }


@pytest.fixture
def ani_matrix_df(small_ani_dict: dict) -> pl.DataFrame:
    """ANI matrix as Polars DataFrame."""
    genomes = list(small_ani_dict.keys())
    data = {"genome": genomes}
    for genome in genomes:
        data[genome] = [small_ani_dict[genome][g] for g in genomes]
    return pl.DataFrame(data)


# =============================================================================
# Temporary File Fixtures
# =============================================================================


@pytest.fixture
def temp_dir() -> Path:
    """Temporary directory that is cleaned up after tests."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def temp_blast_file(temp_dir: Path, multi_hit_blast_data: list[dict]) -> Path:
    """Temporary BLAST tabular file."""
    blast_path = temp_dir / "test.blast.tsv"

    # Add more reads for realistic testing
    all_hits = []
    for i in range(5):
        for hit in multi_hit_blast_data:
            new_hit = hit.copy()
            new_hit["qseqid"] = f"read_{i:03d}"
            all_hits.append(new_hit)

    df = pl.DataFrame(all_hits)
    df.write_csv(blast_path, separator="\t", include_header=False)

    return blast_path


@pytest.fixture
def temp_ani_file(temp_dir: Path, ani_matrix_df: pl.DataFrame) -> Path:
    """Temporary ANI matrix CSV file."""
    ani_path = temp_dir / "test.ani.csv"
    ani_matrix_df.write_csv(ani_path)
    return ani_path


@pytest.fixture
def empty_blast_file(temp_dir: Path) -> Path:
    """Empty BLAST file."""
    blast_path = temp_dir / "empty.blast.tsv"
    blast_path.touch()
    return blast_path


@pytest.fixture
def malformed_blast_file(temp_dir: Path) -> Path:
    """BLAST file with malformed lines."""
    blast_path = temp_dir / "malformed.blast.tsv"
    blast_path.write_text(
        "read_001\tGCF_000123456.1\t98.5\n"  # Too few fields
        "read_002\tGCF_000123456.1\t98.5\t150\t2\t0\t1\t150\t1000\t1150\t1e-50\t250.0\n"  # Valid
        "\n"  # Empty line
        "# Comment line\n"
    )
    return blast_path


# =============================================================================
# Configuration Fixtures
# =============================================================================


@pytest.fixture
def default_scoring_config():
    """Default scoring configuration."""
    from metadarkmatter.models.config import ScoringConfig

    return ScoringConfig()


@pytest.fixture
def custom_scoring_config():
    """Custom scoring configuration with modified thresholds."""
    from metadarkmatter.models.config import ScoringConfig

    return ScoringConfig(
        bitscore_threshold_pct=90.0,
        novelty_known_max=3.0,
        novelty_novel_species_min=6.0,
        novelty_novel_species_max=18.0,
        novelty_novel_genus_min=18.0,
        novelty_novel_genus_max=30.0,
    )


# =============================================================================
# Classifier Fixtures
# =============================================================================


@pytest.fixture
def ani_matrix(small_ani_dict: dict):
    """ANIMatrix instance for testing."""
    from metadarkmatter.core.ani_placement import ANIMatrix

    return ANIMatrix(small_ani_dict)


@pytest.fixture
def classifier(ani_matrix, default_scoring_config):
    """ANIWeightedClassifier instance for testing."""
    from metadarkmatter.core.ani_placement import ANIWeightedClassifier

    return ANIWeightedClassifier(ani_matrix=ani_matrix, config=default_scoring_config)


# =============================================================================
# Genome Name Extraction Test Cases
# =============================================================================


@pytest.fixture
def genome_name_test_cases() -> list[tuple[str, str]]:
    """Test cases for genome name extraction: (sseqid, expected_genome_name)."""
    return [
        # RefSeq format
        ("GCF_000123456.1_ASM123v1_genomic", "GCF_000123456.1"),
        ("GCF_000123456.2", "GCF_000123456.2"),
        # GenBank format
        ("GCA_000789012.1_ASM789v1_genomic", "GCA_000789012.1"),
        ("GCA_000789012.3", "GCA_000789012.3"),
        # RefSeq with NZ prefix
        ("NZ_CP012345.1", "NZ_CP012345.1"),
        ("NZ_ABCD01000001.1", "NZ_ABCD01000001.1"),
        # Fallback cases
        ("custom_genome_name", "custom_genome_name"),
        ("genome with spaces suffix", "genome"),
        # Edge cases
        ("", "unknown"),
        ("   ", "unknown"),
    ]


# =============================================================================
# Classification Scenario Fixtures
# =============================================================================


@pytest.fixture
def known_species_scenario() -> dict:
    """Scenario parameters for known species classification."""
    return {
        "top_hit_identity": 99.0,  # Novelty = 1.0 (< 2.0)
        "placement_uncertainty": 0.3,  # < 0.5
        "expected_call": "Known Species",
    }


@pytest.fixture
def novel_species_scenario() -> dict:
    """Scenario parameters for novel species classification."""
    return {
        "top_hit_identity": 90.0,  # Novelty = 10.0 (5-15)
        "placement_uncertainty": 0.3,  # < 0.5
        "expected_call": "Novel Species",
    }


@pytest.fixture
def novel_genus_scenario() -> dict:
    """Scenario parameters for novel genus classification."""
    return {
        "top_hit_identity": 80.0,  # Novelty = 20.0 (15-25)
        "placement_uncertainty": 1.5,  # < 2.0
        "expected_call": "Novel Genus",
    }


@pytest.fixture
def conserved_region_scenario() -> dict:
    """Scenario parameters for conserved region classification."""
    return {
        "top_hit_identity": 95.0,  # Any novelty
        "placement_uncertainty": 10.0,  # >= 5.0
        "expected_call": "Conserved Region",
    }


# =============================================================================
# CLI Testing Fixtures
# =============================================================================


@pytest.fixture
def cli_runner():
    """Typer CLI runner for testing commands."""
    from typer.testing import CliRunner

    return CliRunner()


@pytest.fixture
def temp_classification_file(temp_dir: Path) -> Path:
    """Temporary classification CSV file with test data."""
    data = {
        "read_id": [f"read_{i:03d}" for i in range(10)],
        "best_match_genome": ["GCF_000123456.1"] * 5 + ["GCF_000789012.1"] * 5,
        "taxonomic_call": (
            ["Known Species"] * 3 + ["Novel Species"] * 2 +
            ["Novel Genus"] * 2 + ["Conserved Region"] * 3
        ),
        "novelty_index": [1.0, 2.0, 3.0, 8.0, 12.0, 18.0, 22.0, 5.0, 6.0, 7.0],
        "placement_uncertainty": [0.2, 0.3, 0.4, 0.5, 0.6, 1.0, 1.5, 8.0, 9.0, 10.0],
        "top_hit_identity": [99.0, 98.0, 97.0, 92.0, 88.0, 82.0, 78.0, 95.0, 94.0, 93.0],
    }
    df = pl.DataFrame(data)
    output_path = temp_dir / "classifications.csv"
    df.write_csv(output_path)
    return output_path


@pytest.fixture
def temp_metadata_file(temp_dir: Path) -> Path:
    """Temporary genome metadata TSV file."""
    data = {
        "accession": ["GCF_000123456.1", "GCF_000789012.1", "GCA_000111222.1"],
        "species": ["Francisella tularensis", "Francisella novicida", "Francisella philomiragia"],
        "genus": ["Francisella", "Francisella", "Francisella"],
        "family": ["Francisellaceae", "Francisellaceae", "Francisellaceae"],
    }
    df = pl.DataFrame(data)
    output_path = temp_dir / "genome_metadata.tsv"
    df.write_csv(output_path, separator="\t")
    return output_path
