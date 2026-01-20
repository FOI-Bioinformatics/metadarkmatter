"""Unit tests for genome accession models.

Tests for GenomeAccession and AccessionList models.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metadarkmatter.models.genomes import AccessionList, GenomeAccession


class TestGenomeAccession:
    """Tests for GenomeAccession model."""

    def test_create_accession(self):
        """Should create accession with all fields."""
        accession = GenomeAccession(
            accession="GCF_000005845.2",
            gtdb_taxonomy="d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli",
            species="Escherichia coli",
            genome_size=4641652,
        )
        assert accession.accession == "GCF_000005845.2"
        assert accession.species == "Escherichia coli"
        assert accession.genome_size == 4641652

    def test_frozen(self):
        """Should be immutable."""
        accession = GenomeAccession(
            accession="GCF_000005845.2",
            gtdb_taxonomy="d__Bacteria",
            species="E. coli",
        )
        with pytest.raises(Exception):  # ValidationError or AttributeError
            accession.accession = "GCF_999999999.1"

    def test_extract_genus(self):
        """Should extract genus from taxonomy."""
        accession = GenomeAccession(
            accession="GCF_000005845.2",
            gtdb_taxonomy="d__Bacteria;p__Proteobacteria;g__Escherichia;s__Escherichia coli",
            species="Escherichia coli",
        )
        assert accession.genus == "Escherichia"

    def test_extract_genus_missing(self):
        """Should return empty string if no genus in taxonomy."""
        accession = GenomeAccession(
            accession="GCF_000005845.2",
            gtdb_taxonomy="d__Bacteria;p__Proteobacteria",
            species="Unknown",
        )
        assert accession.genus == ""

    def test_extract_family(self):
        """Should extract family from taxonomy."""
        accession = GenomeAccession(
            accession="GCF_000005845.2",
            gtdb_taxonomy="d__Bacteria;f__Enterobacteriaceae;g__Escherichia",
            species="Escherichia coli",
        )
        assert accession.family == "Enterobacteriaceae"

    def test_extract_family_missing(self):
        """Should return empty string if no family in taxonomy."""
        accession = GenomeAccession(
            accession="GCF_000005845.2",
            gtdb_taxonomy="d__Bacteria;g__Escherichia",
            species="Unknown",
        )
        assert accession.family == ""

    def test_none_genome_size(self):
        """Should allow None genome_size."""
        accession = GenomeAccession(
            accession="GCF_000005845.2",
            gtdb_taxonomy="d__Bacteria",
            species="E. coli",
            genome_size=None,
        )
        assert accession.genome_size is None


class TestAccessionList:
    """Tests for AccessionList model."""

    def test_create_empty_list(self):
        """Should create empty accession list."""
        acc_list = AccessionList(taxon="f__Enterobacteriaceae")
        assert acc_list.taxon == "f__Enterobacteriaceae"
        assert acc_list.total_count == 0
        assert acc_list.accessions == []

    def test_create_with_accessions(self):
        """Should create list with accessions."""
        accessions = [
            GenomeAccession(
                accession="GCF_000005845.2",
                gtdb_taxonomy="g__Escherichia;s__Escherichia coli",
                species="Escherichia coli",
            ),
            GenomeAccession(
                accession="GCF_000006765.1",
                gtdb_taxonomy="g__Salmonella;s__Salmonella enterica",
                species="Salmonella enterica",
            ),
        ]

        acc_list = AccessionList(
            taxon="f__Enterobacteriaceae",
            accessions=accessions,
            genus_counts={"Escherichia": 1, "Salmonella": 1},
        )

        assert acc_list.total_count == 2
        assert acc_list.genus_count == 2

    def test_total_count(self):
        """Should return correct total count."""
        accessions = [
            GenomeAccession(accession=f"GCF_{i}", gtdb_taxonomy="", species="")
            for i in range(5)
        ]
        acc_list = AccessionList(taxon="test", accessions=accessions)
        assert acc_list.total_count == 5

    def test_species_count(self):
        """Should return unique species count."""
        accessions = [
            GenomeAccession(accession="GCF_1", gtdb_taxonomy="", species="Species A"),
            GenomeAccession(accession="GCF_2", gtdb_taxonomy="", species="Species A"),
            GenomeAccession(accession="GCF_3", gtdb_taxonomy="", species="Species B"),
        ]
        acc_list = AccessionList(taxon="test", accessions=accessions)
        assert acc_list.species_count == 2

    def test_genus_count(self):
        """Should return genus count from genus_counts dict."""
        acc_list = AccessionList(
            taxon="test",
            genus_counts={"Genus1": 5, "Genus2": 3, "Genus3": 2},
        )
        assert acc_list.genus_count == 3

    def test_get_accession_strings(self):
        """Should return list of accession strings."""
        accessions = [
            GenomeAccession(accession="GCF_000005845.2", gtdb_taxonomy="", species=""),
            GenomeAccession(accession="GCF_000006765.1", gtdb_taxonomy="", species=""),
        ]
        acc_list = AccessionList(taxon="test", accessions=accessions)

        acc_strings = acc_list.get_accession_strings()
        assert acc_strings == ["GCF_000005845.2", "GCF_000006765.1"]


class TestAccessionListTSV:
    """Tests for TSV read/write operations."""

    def test_to_tsv_basic(self, tmp_path):
        """Should write basic TSV file."""
        accessions = [
            GenomeAccession(
                accession="GCF_000005845.2",
                gtdb_taxonomy="d__Bacteria;g__Escherichia;s__Escherichia coli",
                species="Escherichia coli",
            ),
        ]
        acc_list = AccessionList(taxon="f__Test", accessions=accessions)

        output_path = tmp_path / "accessions.tsv"
        acc_list.to_tsv(output_path)

        assert output_path.exists()
        content = output_path.read_text()
        assert "GCF_000005845.2" in content
        assert "Escherichia coli" in content

    def test_to_tsv_with_metadata(self, tmp_path):
        """Should include metadata columns when requested."""
        accessions = [
            GenomeAccession(
                accession="GCF_000005845.2",
                gtdb_taxonomy="d__Bacteria;f__Enterobacteriaceae;g__Escherichia",
                species="Escherichia coli",
                genome_size=4641652,
            ),
        ]
        acc_list = AccessionList(taxon="f__Test", accessions=accessions)

        output_path = tmp_path / "accessions.tsv"
        acc_list.to_tsv(output_path, include_metadata=True)

        content = output_path.read_text()
        assert "genome_size" in content
        assert "genus" in content
        assert "family" in content
        assert "4641652" in content
        assert "Escherichia" in content

    def test_from_tsv(self, tmp_path):
        """Should read TSV file back into AccessionList."""
        # Write a TSV file
        tsv_content = "accession\tgtdb_taxonomy\tspecies\n"
        tsv_content += "GCF_000005845.2\td__Bacteria;g__Escherichia\tEscherichia coli\n"
        tsv_content += "GCF_000006765.1\td__Bacteria;g__Salmonella\tSalmonella enterica\n"

        tsv_path = tmp_path / "input.tsv"
        tsv_path.write_text(tsv_content)

        acc_list = AccessionList.from_tsv(tsv_path)

        assert acc_list.total_count == 2
        assert acc_list.accessions[0].accession == "GCF_000005845.2"
        assert acc_list.accessions[1].species == "Salmonella enterica"

    def test_roundtrip(self, tmp_path):
        """Should preserve data through write/read cycle."""
        original = AccessionList(
            taxon="f__Enterobacteriaceae",
            accessions=[
                GenomeAccession(
                    accession="GCF_000005845.2",
                    gtdb_taxonomy="d__Bacteria;g__Escherichia;s__Escherichia coli",
                    species="Escherichia coli",
                    genome_size=4641652,
                ),
                GenomeAccession(
                    accession="GCF_000006765.1",
                    gtdb_taxonomy="d__Bacteria;g__Salmonella;s__Salmonella enterica",
                    species="Salmonella enterica",
                    genome_size=4857432,
                ),
            ],
            genus_counts={"Escherichia": 1, "Salmonella": 1},
        )

        tsv_path = tmp_path / "roundtrip.tsv"
        original.to_tsv(tsv_path, include_metadata=True)
        loaded = AccessionList.from_tsv(tsv_path)

        assert loaded.total_count == original.total_count
        assert loaded.accessions[0].accession == original.accessions[0].accession
        assert loaded.accessions[1].species == original.accessions[1].species
