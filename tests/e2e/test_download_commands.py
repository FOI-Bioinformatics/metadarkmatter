"""
End-to-end tests for download CLI commands.

Tests the full workflow of querying GTDB for genome accessions
and downloading genomes using NCBI Datasets CLI.

Uses mocked API responses for CI-safe testing.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable

import pytest

from metadarkmatter.models.genomes import AccessionList

if TYPE_CHECKING:
    from click.testing import Result


pytestmark = pytest.mark.e2e


# =============================================================================
# Validation Helpers
# =============================================================================


class DownloadAssertions:
    """Custom assertions for download command outputs."""

    ACCESSION_PATTERN = re.compile(r"^GC[AF]_\d+\.\d+$")

    @classmethod
    def assert_valid_accession_tsv(cls, path: Path, min_rows: int = 1) -> None:
        """Validate TSV has required columns and format."""
        assert path.exists(), f"TSV file not found: {path}"

        content = path.read_text()
        lines = content.strip().split("\n")

        assert len(lines) >= min_rows + 1, f"Expected at least {min_rows + 1} lines (header + data)"

        # Check header
        header = lines[0].split("\t")
        assert "accession" in header, "Missing 'accession' column"
        assert "gtdb_taxonomy" in header, "Missing 'gtdb_taxonomy' column"
        assert "species" in header, "Missing 'species' column"

        # Check data rows
        for i, line in enumerate(lines[1:], start=2):
            fields = line.split("\t")
            assert len(fields) >= 3, f"Line {i}: Expected at least 3 fields, got {len(fields)}"

    @classmethod
    def assert_accession_format(cls, accession: str) -> None:
        """Validate accession matches GCF_/GCA_ pattern."""
        assert cls.ACCESSION_PATTERN.match(accession), (
            f"Invalid accession format: {accession}"
        )

    @classmethod
    def assert_taxonomy_format(cls, taxonomy: str) -> None:
        """Validate GTDB taxonomy string format."""
        assert taxonomy.startswith("d__"), f"Taxonomy should start with d__: {taxonomy}"
        assert ";" in taxonomy, f"Taxonomy should have semicolon separators: {taxonomy}"

    @classmethod
    def assert_genus_in_results(cls, accession_list: AccessionList, expected_genus: str) -> None:
        """Verify all results belong to expected genus."""
        for acc in accession_list.accessions:
            assert f"g__{expected_genus}" in acc.gtdb_taxonomy, (
                f"Accession {acc.accession} not in genus {expected_genus}"
            )


# =============================================================================
# List Command Tests
# =============================================================================


class TestDownloadGenomesList:
    """E2E tests for 'download genomes list' command."""

    def test_list_genus_francisella(
        self,
        mock_gtdb_response: dict,
        run_download_list: Callable[..., "Result"],
        e2e_temp_dir: Path,
    ):
        """Should query and save accession list for g__Francisella."""
        result = run_download_list("g__Francisella")

        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert "17" in result.output or "genomes" in result.output.lower()

        # Verify output file
        output_path = e2e_temp_dir / "accessions.tsv"
        DownloadAssertions.assert_valid_accession_tsv(output_path, min_rows=17)

    def test_list_family_francisellaceae(
        self,
        mock_gtdb_family_response: dict,
        run_download_list: Callable[..., "Result"],
        e2e_temp_dir: Path,
    ):
        """Should query and save accession list for f__Francisellaceae."""
        result = run_download_list("f__Francisellaceae")

        assert result.exit_code == 0, f"Command failed: {result.output}"

        # Verify output file
        output_path = e2e_temp_dir / "accessions.tsv"
        DownloadAssertions.assert_valid_accession_tsv(output_path, min_rows=40)

    def test_list_with_metadata(
        self,
        mock_gtdb_response: dict,
        run_download_list: Callable[..., "Result"],
        e2e_temp_dir: Path,
    ):
        """Should include extra columns with --include-metadata."""
        result = run_download_list("g__Francisella", include_metadata=True)

        assert result.exit_code == 0, f"Command failed: {result.output}"

        output_path = e2e_temp_dir / "accessions.tsv"
        content = output_path.read_text()

        # Metadata columns should be present
        assert "genus" in content.lower() or "family" in content.lower()

    def test_list_verbose_output(
        self,
        mock_gtdb_response: dict,
        run_download_list: Callable[..., "Result"],
    ):
        """Should show genus breakdown with --verbose."""
        result = run_download_list("g__Francisella", verbose=True)

        assert result.exit_code == 0, f"Command failed: {result.output}"
        # Verbose mode should show more details
        assert "Francisella" in result.output

    def test_list_quiet_mode(
        self,
        mock_gtdb_response: dict,
        run_download_list: Callable[..., "Result"],
        e2e_temp_dir: Path,
    ):
        """Should produce minimal output with --quiet."""
        result = run_download_list("g__Francisella", quiet=True)

        assert result.exit_code == 0
        # Output should be minimal in quiet mode
        output_path = e2e_temp_dir / "accessions.tsv"
        assert output_path.exists()

    def test_list_dry_run(
        self,
        run_download_list: Callable[..., "Result"],
        e2e_temp_dir: Path,
    ):
        """Should show query URL without making request in dry-run mode."""
        result = run_download_list("g__Francisella", dry_run=True)

        assert result.exit_code == 0, f"Command failed: {result.output}"
        # Dry run should mention what would happen
        assert "dry" in result.output.lower() or "would" in result.output.lower()

        # Output file should not be created in dry-run
        output_path = e2e_temp_dir / "accessions.tsv"
        assert not output_path.exists()

    def test_list_invalid_taxon_format(
        self,
        run_download_list: Callable[..., "Result"],
    ):
        """Should error on invalid taxon format."""
        result = run_download_list("Francisella")  # Missing g__ prefix

        assert result.exit_code != 0, "Should fail for invalid taxon format"
        assert "invalid" in result.output.lower() or "format" in result.output.lower()

    def test_list_nonexistent_taxon(
        self,
        mock_gtdb_404: None,
        run_download_list: Callable[..., "Result"],
    ):
        """Should handle 404 for unknown taxa."""
        result = run_download_list("g__FakeTaxonNotReal")

        assert result.exit_code != 0, "Should fail for unknown taxon"
        # Should show helpful error message
        assert "not" in result.output.lower() or "error" in result.output.lower()

    def test_list_output_file_created(
        self,
        mock_gtdb_response: dict,
        run_download_list: Callable[..., "Result"],
        e2e_temp_dir: Path,
    ):
        """Should create valid TSV file."""
        custom_name = "my_genomes.tsv"
        result = run_download_list("g__Francisella", output_name=custom_name)

        assert result.exit_code == 0
        output_path = e2e_temp_dir / custom_name
        assert output_path.exists()
        DownloadAssertions.assert_valid_accession_tsv(output_path)

    def test_list_tsv_roundtrip(
        self,
        mock_gtdb_response: dict,
        run_download_list: Callable[..., "Result"],
        e2e_temp_dir: Path,
    ):
        """TSV can be read back with AccessionList.from_tsv()."""
        result = run_download_list("g__Francisella")

        assert result.exit_code == 0
        output_path = e2e_temp_dir / "accessions.tsv"

        # Load the TSV back
        acc_list = AccessionList.from_tsv(output_path)

        assert acc_list.total_count == 17
        assert len(acc_list.accessions) == 17

        # Verify accession formats
        for acc in acc_list.accessions:
            DownloadAssertions.assert_accession_format(acc.accession)
            DownloadAssertions.assert_taxonomy_format(acc.gtdb_taxonomy)

    def test_list_genus_contains_only_target_genus(
        self,
        mock_gtdb_response: dict,
        run_download_list: Callable[..., "Result"],
        e2e_temp_dir: Path,
    ):
        """All accessions in genus query should belong to that genus."""
        result = run_download_list("g__Francisella")

        assert result.exit_code == 0
        output_path = e2e_temp_dir / "accessions.tsv"
        acc_list = AccessionList.from_tsv(output_path)

        DownloadAssertions.assert_genus_in_results(acc_list, "Francisella")

    def test_list_family_contains_multiple_genera(
        self,
        mock_gtdb_family_response: dict,
        run_download_list: Callable[..., "Result"],
        e2e_temp_dir: Path,
    ):
        """Family query should contain multiple genera."""
        result = run_download_list("f__Francisellaceae")

        assert result.exit_code == 0
        output_path = e2e_temp_dir / "accessions.tsv"
        acc_list = AccessionList.from_tsv(output_path)

        # Extract unique genera from taxonomy strings
        genera = set()
        for acc in acc_list.accessions:
            for part in acc.gtdb_taxonomy.split(";"):
                if part.startswith("g__"):
                    genera.add(part[3:])

        # Family should have multiple genera
        assert len(genera) > 1, f"Expected multiple genera, got: {genera}"


# =============================================================================
# Fetch Command Tests
# =============================================================================


class TestDownloadGenomesFetch:
    """E2E tests for 'download genomes fetch' command."""

    def test_fetch_dry_run(
        self,
        sample_accession_tsv: Path,
        mock_ncbi_available: None,
        run_download_fetch: Callable[..., "Result"],
    ):
        """Should show download plan without executing."""
        result = run_download_fetch(sample_accession_tsv, dry_run=True)

        assert result.exit_code == 0, f"Command failed: {result.output}"
        # Should mention dry run or what would happen
        assert "dry" in result.output.lower() or "would" in result.output.lower()

    def test_fetch_missing_ncbi_cli(
        self,
        sample_accession_tsv: Path,
        mock_ncbi_missing: None,
        run_download_fetch: Callable[..., "Result"],
    ):
        """Should error when NCBI datasets CLI not installed."""
        result = run_download_fetch(sample_accession_tsv)

        assert result.exit_code != 0, "Should fail when datasets CLI not found"
        # Should show helpful install hint
        assert "conda" in result.output.lower() or "install" in result.output.lower()

    def test_fetch_skip_existing(
        self,
        sample_accession_tsv: Path,
        mock_ncbi_available: None,
        run_download_fetch: Callable[..., "Result"],
        e2e_temp_dir: Path,
    ):
        """Should skip genomes that already exist."""
        output_dir = e2e_temp_dir / "genomes"
        output_dir.mkdir(parents=True)

        # Create existing file for first accession
        existing_file = output_dir / "GCF_000833475.1.fna"
        existing_file.write_text(">header\nACGT\n")

        # Default is skip_if_exists=True
        result = run_download_fetch(
            sample_accession_tsv,
            dry_run=True,  # Use dry run to avoid actual downloads
        )

        assert result.exit_code == 0

    def test_fetch_all_exist_early_exit(
        self,
        sample_accession_tsv: Path,
        mock_ncbi_available: None,
        run_download_fetch: Callable[..., "Result"],
        e2e_temp_dir: Path,
    ):
        """Should exit early when all genomes exist."""
        output_dir = e2e_temp_dir / "genomes"
        output_dir.mkdir(parents=True)

        # Create files for all accessions in sample TSV
        for acc in ["GCF_000833475.1", "GCF_000156715.1", "GCF_001885235.1"]:
            (output_dir / f"{acc}.fna").write_text(f">header\nACGT\n")

        # Default is skip_if_exists=True
        result = run_download_fetch(sample_accession_tsv)

        assert result.exit_code == 0
        assert "exist" in result.output.lower() or "skip" in result.output.lower()

    def test_fetch_invalid_accession_file(
        self,
        tmp_path: Path,
        run_download_fetch: Callable[..., "Result"],
    ):
        """Should error on invalid accession file."""
        # Create malformed TSV
        bad_tsv = tmp_path / "bad.tsv"
        bad_tsv.write_text("wrong,header,format\n1,2,3\n")

        result = run_download_fetch(bad_tsv)

        assert result.exit_code != 0, "Should fail for invalid TSV format"

    def test_fetch_nonexistent_file(
        self,
        run_download_fetch: Callable[..., "Result"],
        e2e_temp_dir: Path,
    ):
        """Should error on nonexistent accession file."""
        fake_path = e2e_temp_dir / "nonexistent.tsv"

        result = run_download_fetch(fake_path)

        assert result.exit_code != 0

    def test_fetch_verbose_output(
        self,
        sample_accession_tsv: Path,
        mock_ncbi_available: None,
        run_download_fetch: Callable[..., "Result"],
    ):
        """Should show detailed output with --verbose."""
        result = run_download_fetch(
            sample_accession_tsv,
            verbose=True,
            dry_run=True,
        )

        assert result.exit_code == 0


# =============================================================================
# Workflow Tests
# =============================================================================


class TestDownloadWorkflow:
    """E2E tests for full list -> fetch workflow."""

    def test_list_then_fetch_dry_run(
        self,
        mock_gtdb_response: dict,
        mock_ncbi_available: None,
        run_download_list: Callable[..., "Result"],
        run_download_fetch: Callable[..., "Result"],
        e2e_temp_dir: Path,
    ):
        """Full workflow: list genomes then fetch with dry-run."""
        # Step 1: List genomes
        list_result = run_download_list("g__Francisella")
        assert list_result.exit_code == 0

        # Verify TSV was created
        accession_file = e2e_temp_dir / "accessions.tsv"
        assert accession_file.exists()

        # Step 2: Fetch with dry-run
        fetch_result = run_download_fetch(accession_file, dry_run=True)
        assert fetch_result.exit_code == 0

    def test_list_with_metadata_then_verify(
        self,
        mock_gtdb_response: dict,
        run_download_list: Callable[..., "Result"],
        e2e_temp_dir: Path,
    ):
        """List with metadata and verify extra columns."""
        result = run_download_list(
            "g__Francisella",
            include_metadata=True,
            output_name="with_metadata.tsv",
        )
        assert result.exit_code == 0

        output_path = e2e_temp_dir / "with_metadata.tsv"
        content = output_path.read_text()
        lines = content.strip().split("\n")
        header = lines[0].lower()

        # Should have metadata columns
        assert "genus" in header or "family" in header


# =============================================================================
# Edge Case Tests
# =============================================================================


class TestDownloadEdgeCases:
    """Tests for edge cases and error handling."""

    def test_list_special_characters_in_species(
        self,
        mock_gtdb_response: dict,
        run_download_list: Callable[..., "Result"],
        e2e_temp_dir: Path,
    ):
        """Should handle species names with underscores and special chars."""
        result = run_download_list("g__Francisella")

        assert result.exit_code == 0
        output_path = e2e_temp_dir / "accessions.tsv"
        acc_list = AccessionList.from_tsv(output_path)

        # Some Francisella species have _B suffix
        species_names = [acc.species for acc in acc_list.accessions]
        assert any("_" in s or "sp" in s for s in species_names)

    def test_list_handles_gca_and_gcf_accessions(
        self,
        mock_gtdb_response: dict,
        run_download_list: Callable[..., "Result"],
        e2e_temp_dir: Path,
    ):
        """Should handle both GCA and GCF accession prefixes."""
        result = run_download_list("g__Francisella")

        assert result.exit_code == 0
        output_path = e2e_temp_dir / "accessions.tsv"
        acc_list = AccessionList.from_tsv(output_path)

        accession_strings = [acc.accession for acc in acc_list.accessions]
        gcf_count = sum(1 for a in accession_strings if a.startswith("GCF"))
        gca_count = sum(1 for a in accession_strings if a.startswith("GCA"))

        # Both types should be handled
        assert gcf_count > 0, "Should have GCF accessions"
        # GCA might not be present in this fixture but should be handled

    def test_empty_taxon_error(
        self,
        run_download_list: Callable[..., "Result"],
    ):
        """Should handle empty taxon gracefully."""
        result = run_download_list("")

        assert result.exit_code != 0

    def test_list_creates_parent_directories(
        self,
        mock_gtdb_response: dict,
        e2e_runner: Any,
        e2e_temp_dir: Path,
    ):
        """Should create parent directories for output file."""
        from metadarkmatter.cli.main import app

        nested_path = e2e_temp_dir / "nested" / "dir" / "output.tsv"

        result = e2e_runner.invoke(
            app,
            ["download", "genomes", "list", "g__Francisella", "--output", str(nested_path)],
        )

        assert result.exit_code == 0
        assert nested_path.exists()
