"""
Integration tests for download commands with live API access.

These tests hit real GTDB and NCBI APIs and should be skipped in CI.
Run manually for verification against production services.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from metadarkmatter.clients.gtdb import GTDBClient, GTDBQueryResult
from metadarkmatter.models.genomes import AccessionList

# Skip all tests in CI
pytestmark = [
    pytest.mark.integration,
    pytest.mark.slow,
    pytest.mark.skipif(
        os.environ.get("CI") == "true",
        reason="Skip live API tests in CI",
    ),
]


# =============================================================================
# Live GTDB API Tests
# =============================================================================


class TestLiveGTDBAPI:
    """Integration tests hitting real GTDB API."""

    def test_live_query_francisella(self):
        """Query real GTDB API for g__Francisella."""
        with GTDBClient(timeout=30.0) as client:
            result = client.query_genomes("g__Francisella")

        assert isinstance(result, GTDBQueryResult)
        assert result.total_count >= 15, (
            f"Expected at least 15 Francisella genomes, got {result.total_count}"
        )
        assert result.taxon == "g__Francisella"

        # Verify genome data structure
        for genome in result.genomes:
            assert genome.accession.startswith("GC")
            assert "Francisella" in genome.gtdb_taxonomy
            assert genome.species  # Species should not be empty

    def test_live_query_francisellaceae(self):
        """Query real GTDB API for f__Francisellaceae."""
        with GTDBClient(timeout=30.0) as client:
            result = client.query_genomes("f__Francisellaceae")

        assert isinstance(result, GTDBQueryResult)
        assert result.total_count >= 40, (
            f"Expected at least 40 Francisellaceae genomes, got {result.total_count}"
        )

        # Family should have multiple genera
        assert len(result.genus_counts) >= 5, (
            f"Expected at least 5 genera, got {len(result.genus_counts)}"
        )

        # Francisella should be present
        assert "Francisella" in result.genus_counts

    def test_live_query_small_genus(self):
        """Query a smaller genus to verify response parsing."""
        with GTDBClient(timeout=30.0) as client:
            result = client.query_genomes("g__Allofrancisella")

        assert isinstance(result, GTDBQueryResult)
        assert result.total_count >= 1

        for genome in result.genomes:
            assert "Allofrancisella" in genome.gtdb_taxonomy

    def test_live_query_unknown_taxon(self):
        """Should handle unknown taxon gracefully."""
        from metadarkmatter.clients.gtdb import GTDBAPIError

        with GTDBClient(timeout=30.0) as client:
            with pytest.raises(GTDBAPIError) as exc_info:
                client.query_genomes("g__FakeTaxonDoesNotExist")

        assert exc_info.value.status_code == 404

    def test_live_query_all_genomes(self):
        """Query with all genomes (not just representatives)."""
        with GTDBClient(timeout=30.0) as client:
            reps_only = client.query_genomes("g__Francisella", representatives_only=True)
            all_genomes = client.query_genomes("g__Francisella", representatives_only=False)

        # All genomes should be >= representative genomes
        assert all_genomes.total_count >= reps_only.total_count


# =============================================================================
# Live CLI Tests
# =============================================================================


class TestLiveCLI:
    """Integration tests for CLI commands with live APIs."""

    def test_live_cli_list_command(self, tmp_path: Path):
        """Full CLI list command against real API."""
        from typer.testing import CliRunner

        from metadarkmatter.cli.main import app

        runner = CliRunner()
        output_path = tmp_path / "live_accessions.tsv"

        result = runner.invoke(
            app,
            [
                "download", "genomes", "list",
                "g__Francisella",
                "--output", str(output_path),
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert output_path.exists()

        # Verify TSV content
        acc_list = AccessionList.from_tsv(output_path)
        assert acc_list.total_count >= 15

    def test_live_cli_list_family(self, tmp_path: Path):
        """Full CLI list command for family against real API."""
        from typer.testing import CliRunner

        from metadarkmatter.cli.main import app

        runner = CliRunner()
        output_path = tmp_path / "live_family.tsv"

        result = runner.invoke(
            app,
            [
                "download", "genomes", "list",
                "f__Francisellaceae",
                "--output", str(output_path),
                "--verbose",
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert output_path.exists()

        acc_list = AccessionList.from_tsv(output_path)
        assert acc_list.total_count >= 40

        # Extract genera from results
        genera = set()
        for acc in acc_list.accessions:
            for part in acc.gtdb_taxonomy.split(";"):
                if part.startswith("g__"):
                    genera.add(part[3:])

        assert len(genera) >= 5, f"Expected multiple genera: {genera}"

    def test_live_cli_list_with_metadata(self, tmp_path: Path):
        """CLI list with metadata columns."""
        from typer.testing import CliRunner

        from metadarkmatter.cli.main import app

        runner = CliRunner()
        output_path = tmp_path / "metadata.tsv"

        result = runner.invoke(
            app,
            [
                "download", "genomes", "list",
                "g__Francisella",
                "--output", str(output_path),
                "--include-metadata",
            ],
        )

        assert result.exit_code == 0
        content = output_path.read_text()

        # Check for metadata columns
        assert "genus" in content.lower() or "family" in content.lower()

    def test_live_cli_invalid_taxon(self):
        """CLI should error on invalid taxon format."""
        from typer.testing import CliRunner

        from metadarkmatter.cli.main import app

        runner = CliRunner()

        result = runner.invoke(
            app,
            [
                "download", "genomes", "list",
                "Francisella",  # Missing prefix
                "--output", "/tmp/test.tsv",
            ],
        )

        assert result.exit_code != 0
        assert "invalid" in result.output.lower() or "format" in result.output.lower()


# =============================================================================
# Validation Tests
# =============================================================================


class TestLiveDataValidation:
    """Tests to validate data quality from live API."""

    def test_all_accessions_have_valid_format(self):
        """All accessions should match NCBI format."""
        import re

        pattern = re.compile(r"^GC[AF]_\d{9}\.\d+$")

        with GTDBClient(timeout=30.0) as client:
            result = client.query_genomes("g__Francisella")

        for genome in result.genomes:
            assert pattern.match(genome.accession), (
                f"Invalid accession format: {genome.accession}"
            )

    def test_all_species_have_proper_naming(self):
        """Species names should be properly formatted."""
        with GTDBClient(timeout=30.0) as client:
            result = client.query_genomes("g__Francisella")

        for genome in result.genomes:
            # Species should start with genus name
            assert genome.species.startswith("Francisella"), (
                f"Species should start with genus: {genome.species}"
            )

    def test_taxonomy_strings_complete(self):
        """Taxonomy strings should have complete hierarchy."""
        with GTDBClient(timeout=30.0) as client:
            result = client.query_genomes("g__Francisella")

        for genome in result.genomes:
            taxonomy = genome.gtdb_taxonomy

            # Should have domain at minimum
            assert taxonomy.startswith("d__"), f"Missing domain: {taxonomy}"

            # Should have genus
            assert "g__" in taxonomy, f"Missing genus: {taxonomy}"

            # Should have species for most
            assert "s__" in taxonomy, f"Missing species: {taxonomy}"
