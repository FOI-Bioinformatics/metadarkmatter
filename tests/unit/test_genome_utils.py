"""Tests for genome utility functions.

Tests accession extraction from filenames, including standard NCBI
filenames and custom prefixed filenames with embedded taxonomy.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

from metadarkmatter.core.genome_utils import extract_accession_from_filename
from metadarkmatter.cli.mapping import _parse_structured_filename, app


# =============================================================================
# extract_accession_from_filename tests
# =============================================================================


class TestExtractAccessionFromFilename:
    """Tests for extract_accession_from_filename()."""

    def test_standard_ncbi_filename(self) -> None:
        """Should extract accession from standard NCBI download filenames."""
        assert (
            extract_accession_from_filename("GCF_000005845.2_ASM584v2_genomic.fna")
            == "GCF_000005845.2"
        )

    def test_genbank_accession(self) -> None:
        """Should extract GCA accessions the same as GCF."""
        assert (
            extract_accession_from_filename("GCA_000001405.15_GRCh38_genomic.fna")
            == "GCA_000001405.15"
        )

    def test_bare_accession(self) -> None:
        """Should extract accession from minimal filenames."""
        assert (
            extract_accession_from_filename("GCF_000195955.2.fna")
            == "GCF_000195955.2"
        )

    def test_prefixed_filename_family_genus_species(self) -> None:
        """Should extract accession from Family_Genus_species_Accession format."""
        assert (
            extract_accession_from_filename(
                "Francisellaceae_Allofrancisella_frigidaquae_GCA_000710735.1.fasta"
            )
            == "GCA_000710735.1"
        )

    def test_prefixed_filename_with_numeric_species(self) -> None:
        """Should extract accession when species contains numbers."""
        assert (
            extract_accession_from_filename(
                "Francisellaceae_CAJXRW01_sp913060525_GCA_913060525.1.fasta"
            )
            == "GCA_913060525.1"
        )

    def test_prefixed_filename_refseq(self) -> None:
        """Should extract GCF accession from prefixed filename."""
        assert (
            extract_accession_from_filename(
                "Francisellaceae_Francisella_tularensis_GCF_000195955.2.fasta"
            )
            == "GCF_000195955.2"
        )

    def test_gzipped_fasta_extension(self) -> None:
        """Should handle .fasta.gz extension."""
        assert (
            extract_accession_from_filename(
                "Francisellaceae_Francisella_tularensis_GCF_000195955.2.fasta.gz"
            )
            == "GCF_000195955.2"
        )

    def test_fna_gz_extension(self) -> None:
        """Should handle .fna.gz extension."""
        assert (
            extract_accession_from_filename("GCF_000195955.2_genomic.fna.gz")
            == "GCF_000195955.2"
        )

    def test_no_accession_returns_stem(self) -> None:
        """Should return whole stem when no GCF/GCA pattern is found."""
        assert (
            extract_accession_from_filename("custom_genome.fasta")
            == "custom_genome"
        )

    def test_multiple_accessions_returns_first(self) -> None:
        """Should return the first GCF/GCA match when multiple are present."""
        result = extract_accession_from_filename(
            "GCF_000111111.1_extra_GCA_000222222.1.fna"
        )
        assert result == "GCF_000111111.1"

    def test_fa_extension(self) -> None:
        """Should handle .fa extension."""
        assert (
            extract_accession_from_filename("GCA_000123456.1.fa")
            == "GCA_000123456.1"
        )


# =============================================================================
# _parse_structured_filename tests
# =============================================================================


class TestParseStructuredFilename:
    """Tests for _parse_structured_filename()."""

    def test_full_taxonomy_filename(self) -> None:
        """Should parse Family_Genus_species_Accession format."""
        result = _parse_structured_filename(
            "Francisellaceae_Allofrancisella_frigidaquae_GCA_000710735.1.fasta"
        )
        assert result["accession"] == "GCA_000710735.1"
        assert result["family"] == "Francisellaceae"
        assert result["genus"] == "Allofrancisella"
        assert result["species"] == "Allofrancisella frigidaquae"

    def test_alphanumeric_genus_and_species(self) -> None:
        """Should handle non-standard genus/species names from GTDB."""
        result = _parse_structured_filename(
            "Francisellaceae_CAJXRW01_sp913060525_GCA_913060525.1.fasta"
        )
        assert result["accession"] == "GCA_913060525.1"
        assert result["family"] == "Francisellaceae"
        assert result["genus"] == "CAJXRW01"
        assert result["species"] == "CAJXRW01 sp913060525"

    def test_family_override(self) -> None:
        """Should use family_override instead of parsing from filename."""
        result = _parse_structured_filename(
            "Francisella_tularensis_GCF_000195955.2.fasta",
            family_override="Francisellaceae",
        )
        assert result["family"] == "Francisellaceae"
        assert result["genus"] == "Francisella"
        assert result["species"] == "Francisella tularensis"

    def test_standard_ncbi_filename(self) -> None:
        """Should handle standard NCBI filename with empty taxonomy."""
        result = _parse_structured_filename(
            "GCF_000195955.2_ASM584v2_genomic.fna"
        )
        assert result["accession"] == "GCF_000195955.2"
        # Standard NCBI filenames have the accession first, so prefix
        # is empty or contains assembly name, not taxonomy
        assert result["family"] == ""

    def test_bare_accession_filename(self) -> None:
        """Should handle bare accession with no prefix."""
        result = _parse_structured_filename("GCF_000195955.2.fna")
        assert result["accession"] == "GCF_000195955.2"
        assert result["family"] == ""
        assert result["genus"] == ""
        assert result["species"] == ""

    def test_no_accession_raises_error(self) -> None:
        """Should raise ValueError when no GCF/GCA accession is found."""
        with pytest.raises(ValueError, match="Cannot extract"):
            _parse_structured_filename("custom_genome.fasta")


# =============================================================================
# generate-metadata CLI command tests
# =============================================================================


class TestGenerateMetadata:
    """Tests for the generate-metadata CLI command."""

    @pytest.fixture
    def runner(self) -> CliRunner:
        return CliRunner()

    @pytest.fixture
    def genome_dir(self, tmp_path: Path) -> Path:
        """Create a directory with structured genome filenames."""
        gdir = tmp_path / "genomes"
        gdir.mkdir()

        # Create genome files with structured names
        filenames = [
            "Francisellaceae_Francisella_tularensis_GCF_000195955.2.fasta",
            "Francisellaceae_Francisella_tularensis_GCF_000123456.1.fasta",
            "Francisellaceae_Francisella_novicida_GCF_000789012.1.fasta",
            "Francisellaceae_Allofrancisella_frigidaquae_GCA_000710735.1.fasta",
        ]
        for fn in filenames:
            (gdir / fn).write_text(">contig1\nATCG\n")

        return gdir

    def test_generates_metadata_from_structured_filenames(
        self, runner: CliRunner, genome_dir: Path, tmp_path: Path
    ) -> None:
        """Should produce metadata TSV with correct taxonomy columns."""
        output = tmp_path / "metadata.tsv"

        result = runner.invoke(
            app,
            [
                "generate-metadata",
                "--genomes", str(genome_dir),
                "--pattern", "*.fasta",
                "--output", str(output),
            ],
        )

        assert result.exit_code == 0, result.output
        assert output.exists()

        df = pl.read_csv(output, separator="\t")
        assert len(df) == 4
        assert set(df.columns) == {
            "accession", "species", "genus", "family",
            "representative", "gtdb_taxonomy",
        }

        # Check taxonomy parsing
        row = df.filter(pl.col("accession") == "GCA_000710735.1")
        assert row["genus"][0] == "Allofrancisella"
        assert row["species"][0] == "Allofrancisella frigidaquae"
        assert row["family"][0] == "Francisellaceae"

    def test_representative_assignment_without_gtdb(
        self, runner: CliRunner, genome_dir: Path, tmp_path: Path
    ) -> None:
        """Should assign first accession per species as representative."""
        output = tmp_path / "metadata.tsv"

        result = runner.invoke(
            app,
            [
                "generate-metadata",
                "--genomes", str(genome_dir),
                "--pattern", "*.fasta",
                "--output", str(output),
            ],
        )

        assert result.exit_code == 0, result.output
        df = pl.read_csv(output, separator="\t")

        # F. tularensis has two genomes; representative should be first alphabetically
        tularensis = df.filter(
            pl.col("species") == "Francisella tularensis"
        ).sort("accession")
        assert len(tularensis) == 2
        # GCF_000123456.1 < GCF_000195955.2 alphabetically
        expected_rep = "GCF_000123456.1"
        for rep in tularensis["representative"].to_list():
            assert rep == expected_rep

    def test_representative_assignment_with_gtdb_metadata(
        self, runner: CliRunner, genome_dir: Path, tmp_path: Path
    ) -> None:
        """Should use GTDB metadata for representative assignment."""
        # Create a GTDB metadata file
        gtdb_file = tmp_path / "gtdb_metadata.tsv"
        gtdb_df = pl.DataFrame({
            "accession": [
                "GCF_000195955.2", "GCF_000123456.1",
                "GCF_000789012.1", "GCA_000710735.1",
            ],
            "species": [
                "Francisella tularensis", "Francisella tularensis",
                "Francisella novicida", "Allofrancisella frigidaquae",
            ],
            "representative": [
                "GCF_000195955.2", "GCF_000195955.2",
                "GCF_000789012.1", "GCA_000710735.1",
            ],
            "gtdb_taxonomy": [
                "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria",
                "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria",
                "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria",
                "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria",
            ],
        })
        gtdb_df.write_csv(gtdb_file, separator="\t")

        output = tmp_path / "metadata.tsv"

        result = runner.invoke(
            app,
            [
                "generate-metadata",
                "--genomes", str(genome_dir),
                "--pattern", "*.fasta",
                "--gtdb-metadata", str(gtdb_file),
                "--output", str(output),
            ],
        )

        assert result.exit_code == 0, result.output
        df = pl.read_csv(output, separator="\t")

        # F. tularensis representative should be GCF_000195955.2 (from GTDB)
        tularensis = df.filter(
            pl.col("species") == "Francisella tularensis"
        )
        for rep in tularensis["representative"].to_list():
            assert rep == "GCF_000195955.2"

        # Should also include gtdb_taxonomy
        assert all(
            "Bacteria" in t
            for t in df["gtdb_taxonomy"].to_list()
            if t
        )

    def test_family_override(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should use --family override instead of parsing from filename."""
        gdir = tmp_path / "genomes"
        gdir.mkdir()
        # Filename without family prefix
        (gdir / "Francisella_tularensis_GCF_000195955.2.fasta").write_text(
            ">c\nA\n"
        )

        output = tmp_path / "metadata.tsv"

        result = runner.invoke(
            app,
            [
                "generate-metadata",
                "--genomes", str(gdir),
                "--pattern", "*.fasta",
                "--family", "Francisellaceae",
                "--output", str(output),
            ],
        )

        assert result.exit_code == 0, result.output
        df = pl.read_csv(output, separator="\t")
        assert df["family"][0] == "Francisellaceae"
        assert df["genus"][0] == "Francisella"

    def test_no_matching_files(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should exit with error when no files match pattern."""
        gdir = tmp_path / "empty"
        gdir.mkdir()

        output = tmp_path / "metadata.tsv"

        result = runner.invoke(
            app,
            [
                "generate-metadata",
                "--genomes", str(gdir),
                "--pattern", "*.fasta",
                "--output", str(output),
            ],
        )

        assert result.exit_code == 1
        assert "No files matching" in result.output

    def test_quiet_mode(
        self, runner: CliRunner, genome_dir: Path, tmp_path: Path
    ) -> None:
        """Should suppress output in quiet mode."""
        output = tmp_path / "metadata.tsv"

        result = runner.invoke(
            app,
            [
                "generate-metadata",
                "--genomes", str(genome_dir),
                "--pattern", "*.fasta",
                "--output", str(output),
                "--quiet",
            ],
        )

        assert result.exit_code == 0
