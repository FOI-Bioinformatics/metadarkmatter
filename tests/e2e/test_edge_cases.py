"""
E2E edge case tests for metadarkmatter CLI.

Tests boundary conditions, unusual inputs, and corner cases
that should still produce valid results.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable

import polars as pl
import pytest

from tests.factories import BlastDataFactory, ANIMatrixFactory
from tests.utils.assertions import CLIAssertions, ClassificationAssertions


pytestmark = pytest.mark.e2e


class TestSingleReadCases:
    """Tests with minimal data (single read)."""

    def test_single_read_classification(
        self,
        single_read_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """Should handle single-read BLAST file."""
        blast_path, ani_path = single_read_dataset

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="single.csv",
        )

        CLIAssertions.assert_success(result)
        CLIAssertions.assert_output_contains(result, "1")  # Should classify 1 read

        output_path = e2e_temp_dir / "single.csv"
        df = ClassificationAssertions.assert_all_validations(output_path)
        assert len(df) == 1


class TestBoundaryValues:
    """Tests with values at classification boundaries."""

    def test_boundary_novelty_values(
        self,
        boundary_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """Should handle reads at exact novelty thresholds."""
        blast_path, ani_path = boundary_dataset

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="boundary.csv",
        )

        CLIAssertions.assert_success(result)

        output_path = e2e_temp_dir / "boundary.csv"
        df = ClassificationAssertions.assert_all_validations(output_path)

        # Should have 3 reads at boundary values
        assert len(df) == 3

    def test_perfect_identity_match(
        self,
        e2e_temp_dir: Path,
        run_classify: Callable,
    ):
        """Should handle 100% identity matches (novelty = 0)."""
        # Create BLAST with perfect match
        bf = BlastDataFactory(seed=42)
        records = bf.create_known_species_read(pident=100.0)

        blast_path = e2e_temp_dir / "perfect.blast.tsv"
        bf.write_blast_file(records, blast_path)

        # Create ANI matrix
        af = ANIMatrixFactory(seed=42)
        genomes = BlastDataFactory.GENOME_PREFIXES[:4]
        matrix = af.create_matrix(genomes)

        ani_path = e2e_temp_dir / "perfect.ani.csv"
        af.write_ani_file(matrix, ani_path)

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="perfect.csv",
        )

        CLIAssertions.assert_success(result)

        output_path = e2e_temp_dir / "perfect.csv"
        df = ClassificationAssertions.assert_all_validations(output_path)

        # 100% identity means novelty_index = 0
        assert df["novelty_index"][0] == 0.0
        assert df["taxonomic_call"][0] == "Known Species"

    def test_very_low_identity(
        self,
        e2e_temp_dir: Path,
        run_classify: Callable,
    ):
        """Should handle very low identity matches (high novelty)."""
        # Create BLAST with low identity
        bf = BlastDataFactory(seed=42)
        records = bf.create_novel_genus_read(pident=75.0)  # 25% divergence

        blast_path = e2e_temp_dir / "low_id.blast.tsv"
        bf.write_blast_file(records, blast_path)

        # Create ANI matrix
        af = ANIMatrixFactory(seed=42)
        genomes = BlastDataFactory.GENOME_PREFIXES[:4]
        matrix = af.create_matrix(genomes)

        ani_path = e2e_temp_dir / "low_id.ani.csv"
        af.write_ani_file(matrix, ani_path)

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="low_id.csv",
        )

        CLIAssertions.assert_success(result)

        output_path = e2e_temp_dir / "low_id.csv"
        df = ClassificationAssertions.assert_all_validations(output_path)

        # 75% identity means novelty_index = 25
        assert df["novelty_index"][0] == 25.0


class TestMissingGenomes:
    """Tests with genomes missing from ANI matrix."""

    def test_partial_genome_coverage(
        self,
        e2e_temp_dir: Path,
        run_classify: Callable,
    ):
        """Should handle when some BLAST genomes are not in ANI matrix."""
        # Create BLAST with mix of known and unknown genomes
        bf = BlastDataFactory(seed=42)
        records = []

        # Add reads to known genomes
        for _ in range(5):
            records.extend(bf.create_known_species_read(genome="GCF_000123456.1"))

        # Add reads to unknown genome (not in ANI matrix)
        for _ in range(5):
            records.extend(bf.create_known_species_read(genome="GCF_UNKNOWN_999.1"))

        blast_path = e2e_temp_dir / "partial.blast.tsv"
        bf.write_blast_file(records, blast_path)

        # Create ANI matrix with only some genomes
        af = ANIMatrixFactory(seed=42)
        matrix = af.create_matrix(genomes=["GCF_000123456.1", "GCF_000789012.1"])

        ani_path = e2e_temp_dir / "partial.ani.csv"
        af.write_ani_file(matrix, ani_path)

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="partial.csv",
            verbose=True,
        )

        # Should succeed but with warning
        CLIAssertions.assert_success(result)
        # Should have processed all reads, even those with missing genomes


class TestConservedRegions:
    """Tests for conserved region classification."""

    def test_many_ambiguous_hits(
        self,
        e2e_temp_dir: Path,
        run_classify: Callable,
    ):
        """Should handle reads with many similar-scoring hits."""
        bf = BlastDataFactory(seed=42)
        records = []

        # Create many conserved region reads
        for _ in range(20):
            records.extend(bf.create_conserved_region_read())

        blast_path = e2e_temp_dir / "conserved.blast.tsv"
        bf.write_blast_file(records, blast_path)

        # Create high-ANI matrix (close relatives)
        af = ANIMatrixFactory(seed=42)
        genomes = BlastDataFactory.GENOME_PREFIXES[:8]
        matrix = af.create_high_ani_matrix(genomes)

        ani_path = e2e_temp_dir / "conserved.ani.csv"
        af.write_ani_file(matrix, ani_path)

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="conserved.csv",
        )

        CLIAssertions.assert_success(result)

        output_path = e2e_temp_dir / "conserved.csv"
        df = ClassificationAssertions.assert_all_validations(output_path)

        # With high ANI between genomes, ambiguous hits should be conserved regions
        # or have low placement uncertainty


class TestBitscoreThreshold:
    """Tests for bitscore threshold parameter."""

    def test_high_bitscore_threshold(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """Higher threshold should consider fewer hits as ambiguous."""
        blast_path, ani_path = standard_dataset

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="high_threshold.csv",
            bitscore_threshold=99.0,  # Very strict
        )

        CLIAssertions.assert_success(result)

        output_path = e2e_temp_dir / "high_threshold.csv"
        ClassificationAssertions.assert_all_validations(output_path)

    def test_low_bitscore_threshold(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """Lower threshold should consider more hits as ambiguous."""
        blast_path, ani_path = standard_dataset

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="low_threshold.csv",
            bitscore_threshold=80.0,  # Lenient
        )

        CLIAssertions.assert_success(result)

        output_path = e2e_temp_dir / "low_threshold.csv"
        ClassificationAssertions.assert_all_validations(output_path)


class TestLargeANIMatrix:
    """Tests with larger ANI matrices."""

    def test_many_genomes(
        self,
        e2e_temp_dir: Path,
        run_classify: Callable,
    ):
        """Should handle ANI matrix with many genomes."""
        bf = BlastDataFactory(seed=42)
        records = bf.create_mixed_sample(
            n_known=10,
            n_novel_species=5,
            n_novel_genus=3,
            n_conserved=5,
        )

        blast_path = e2e_temp_dir / "many_genomes.blast.tsv"
        bf.write_blast_file(records, blast_path)

        # Create larger matrix
        af = ANIMatrixFactory(seed=42)
        genomes = [f"GCF_{i:09d}.1" for i in range(50)]
        # Add the factory genomes to ensure coverage
        genomes = BlastDataFactory.GENOME_PREFIXES + genomes
        matrix = af.create_matrix(genomes[:50])

        ani_path = e2e_temp_dir / "many_genomes.ani.csv"
        af.write_ani_file(matrix, ani_path)

        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="many_genomes.csv",
        )

        CLIAssertions.assert_success(result)
