"""Unit tests for GTDB API client.

Tests the GTDBClient, taxon validation, and response parsing
without making actual API requests.
"""

from __future__ import annotations

import pytest

from metadarkmatter.clients.gtdb import (
    GTDBClient,
    GTDBGenome,
    GTDBQueryResult,
    InvalidTaxonFormatError,
)


class TestTaxonValidation:
    """Tests for GTDB taxon format validation."""

    def test_valid_family_taxon(self):
        """Should accept valid family taxon."""
        assert GTDBClient.validate_taxon_format("f__Enterobacteriaceae") is True

    def test_valid_genus_taxon(self):
        """Should accept valid genus taxon."""
        assert GTDBClient.validate_taxon_format("g__Escherichia") is True

    def test_valid_species_taxon(self):
        """Should accept valid species taxon with spaces."""
        assert GTDBClient.validate_taxon_format("s__Escherichia coli") is True

    def test_valid_domain_taxon(self):
        """Should accept valid domain taxon."""
        assert GTDBClient.validate_taxon_format("d__Bacteria") is True

    def test_valid_phylum_taxon(self):
        """Should accept valid phylum taxon."""
        assert GTDBClient.validate_taxon_format("p__Proteobacteria") is True

    def test_valid_class_taxon(self):
        """Should accept valid class taxon."""
        assert GTDBClient.validate_taxon_format("c__Gammaproteobacteria") is True

    def test_valid_order_taxon(self):
        """Should accept valid order taxon."""
        assert GTDBClient.validate_taxon_format("o__Enterobacterales") is True

    def test_invalid_missing_prefix(self):
        """Should reject taxon without prefix."""
        assert GTDBClient.validate_taxon_format("Enterobacteriaceae") is False

    def test_invalid_wrong_prefix(self):
        """Should reject taxon with wrong prefix."""
        assert GTDBClient.validate_taxon_format("x__Enterobacteriaceae") is False

    def test_invalid_empty_name(self):
        """Should reject taxon with empty name."""
        assert GTDBClient.validate_taxon_format("f__") is False

    def test_invalid_empty_string(self):
        """Should reject empty string."""
        assert GTDBClient.validate_taxon_format("") is False


class TestInvalidTaxonFormatError:
    """Tests for InvalidTaxonFormatError exception."""

    def test_error_message_includes_taxon(self):
        """Should include the invalid taxon in error message."""
        error = InvalidTaxonFormatError("BadTaxon")
        assert "BadTaxon" in str(error)
        assert error.taxon == "BadTaxon"

    def test_error_suggestion_includes_format_examples(self):
        """Should include format examples in suggestion."""
        error = InvalidTaxonFormatError("BadTaxon")
        assert "f__" in str(error.suggestion)
        assert "g__" in str(error.suggestion)


class TestGTDBGenome:
    """Tests for GTDBGenome dataclass."""

    def test_create_genome(self):
        """Should create genome with all fields."""
        genome = GTDBGenome(
            accession="GCF_000005845.2",
            gtdb_taxonomy="d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli",
            species="Escherichia coli",
            genome_size=4641652,
        )
        assert genome.accession == "GCF_000005845.2"
        assert genome.species == "Escherichia coli"
        assert genome.genome_size == 4641652

    def test_genome_is_frozen(self):
        """Should be immutable."""
        genome = GTDBGenome(
            accession="GCF_000005845.2",
            gtdb_taxonomy="d__Bacteria",
            species="Escherichia coli",
        )
        with pytest.raises(AttributeError):
            genome.accession = "GCF_999999999.1"

    def test_genome_without_size(self):
        """Should allow None genome_size."""
        genome = GTDBGenome(
            accession="GCF_000005845.2",
            gtdb_taxonomy="d__Bacteria",
            species="Escherichia coli",
            genome_size=None,
        )
        assert genome.genome_size is None


class TestGTDBQueryResult:
    """Tests for GTDBQueryResult dataclass."""

    def test_create_result(self):
        """Should create result with all fields."""
        genomes = (
            GTDBGenome(
                accession="GCF_000005845.2",
                gtdb_taxonomy="g__Escherichia;s__Escherichia coli",
                species="Escherichia coli",
            ),
            GTDBGenome(
                accession="GCF_000006765.1",
                gtdb_taxonomy="g__Salmonella;s__Salmonella enterica",
                species="Salmonella enterica",
            ),
        )

        result = GTDBQueryResult(
            taxon="f__Enterobacteriaceae",
            genomes=genomes,
            total_count=2,
            genus_counts={"Escherichia": 1, "Salmonella": 1},
            species_counts={"Escherichia coli": 1, "Salmonella enterica": 1},
        )

        assert result.taxon == "f__Enterobacteriaceae"
        assert result.total_count == 2
        assert len(result.genomes) == 2
        assert result.genus_counts["Escherichia"] == 1

    def test_result_is_frozen(self):
        """Should be immutable."""
        result = GTDBQueryResult(
            taxon="f__Test",
            genomes=(),
            total_count=0,
            genus_counts={},
            species_counts={},
        )
        with pytest.raises(AttributeError):
            result.taxon = "f__Other"


class TestGTDBClientParsing:
    """Tests for GTDBClient parsing methods."""

    def test_extract_species_from_taxonomy(self):
        """Should extract species from GTDB taxonomy string."""
        taxonomy = "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli"
        species = GTDBClient._extract_species_from_taxonomy(taxonomy)
        assert species == "Escherichia coli"

    def test_extract_species_missing(self):
        """Should return empty string if no species."""
        taxonomy = "d__Bacteria;p__Proteobacteria"
        species = GTDBClient._extract_species_from_taxonomy(taxonomy)
        assert species == ""

    def test_extract_genus_from_taxonomy(self):
        """Should extract genus from GTDB taxonomy string."""
        taxonomy = "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli"
        genus = GTDBClient._extract_genus_from_taxonomy(taxonomy)
        assert genus == "Escherichia"

    def test_extract_genus_missing(self):
        """Should return empty string if no genus."""
        taxonomy = "d__Bacteria;p__Proteobacteria"
        genus = GTDBClient._extract_genus_from_taxonomy(taxonomy)
        assert genus == ""


class TestGTDBClientInit:
    """Tests for GTDBClient initialization."""

    def test_default_timeout(self):
        """Should use default timeout."""
        client = GTDBClient()
        assert client.timeout == 60.0

    def test_custom_timeout(self):
        """Should accept custom timeout."""
        client = GTDBClient(timeout=30.0)
        assert client.timeout == 30.0

    def test_context_manager(self):
        """Should work as context manager."""
        with GTDBClient() as client:
            assert client is not None
        # After exit, client should be closed
        assert client._client is None


class TestGTDBClientQueryValidation:
    """Tests for GTDBClient query validation."""

    def test_query_invalid_taxon_raises(self):
        """Should raise InvalidTaxonFormatError for invalid taxon."""
        client = GTDBClient()
        with pytest.raises(InvalidTaxonFormatError) as exc_info:
            client.query_genomes("InvalidTaxon")
        assert exc_info.value.taxon == "InvalidTaxon"
