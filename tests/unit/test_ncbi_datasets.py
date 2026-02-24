"""Unit tests for NCBI Datasets CLI wrapper.

Tests command building logic without requiring the actual tool
to be installed.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metadarkmatter.external.ncbi_datasets import (
    DownloadOutcome,
    DownloadReport,
    NCBIDatasets,
    _alternative_versions,
    _parse_accession_version,
)


class TestNCBIDatasetsCommandBuilding:
    """Tests for NCBIDatasets command building."""

    def test_basic_download_command(self):
        """Should build basic download command."""
        ncbi = NCBIDatasets()
        # Mock executable to avoid needing datasets installed
        NCBIDatasets._executable_cache["datasets"] = Path("/usr/bin/datasets")

        cmd = ncbi.build_download_command(
            accessions=["GCF_000005845.2"],
            output_file=Path("/tmp/output.zip"),
        )

        assert "/usr/bin/datasets" in cmd
        assert "download" in cmd
        assert "genome" in cmd
        assert "accession" in cmd
        assert "GCF_000005845.2" in cmd
        assert "--filename" in cmd
        assert "/tmp/output.zip" in " ".join(cmd)
        assert "--include" in cmd
        assert "genome" in cmd

    def test_multiple_accessions(self):
        """Should include all accessions in command."""
        ncbi = NCBIDatasets()
        NCBIDatasets._executable_cache["datasets"] = Path("/usr/bin/datasets")

        accessions = ["GCF_000005845.2", "GCF_000006765.1", "GCF_000009045.1"]
        cmd = ncbi.build_download_command(
            accessions=accessions,
            output_file=Path("/tmp/output.zip"),
        )

        for acc in accessions:
            assert acc in cmd

    def test_without_sequence(self):
        """Should not include --include genome when include_sequence=False."""
        ncbi = NCBIDatasets()
        NCBIDatasets._executable_cache["datasets"] = Path("/usr/bin/datasets")

        cmd = ncbi.build_download_command(
            accessions=["GCF_000005845.2"],
            output_file=Path("/tmp/output.zip"),
            include_sequence=False,
        )

        # Should not have "--include" and "genome" together
        cmd_str = " ".join(cmd)
        assert "--include genome" not in cmd_str

    def test_build_command_requires_accessions(self):
        """Should raise error when accessions is empty."""
        ncbi = NCBIDatasets()
        NCBIDatasets._executable_cache["datasets"] = Path("/usr/bin/datasets")

        with pytest.raises(ValueError, match="accessions"):
            ncbi.build_command(accessions=[], output_file=Path("/tmp/out.zip"))

    def test_build_command_requires_output_file(self):
        """Should raise error when output_file is missing."""
        ncbi = NCBIDatasets()
        NCBIDatasets._executable_cache["datasets"] = Path("/usr/bin/datasets")

        with pytest.raises(ValueError, match="output_file"):
            ncbi.build_command(accessions=["GCF_000005845.2"])


class TestNCBIDatasetsAttributes:
    """Tests for NCBIDatasets class attributes."""

    def test_tool_name(self):
        """Should have correct tool name."""
        assert NCBIDatasets.TOOL_NAME == "datasets"

    def test_install_hint(self):
        """Should have conda install hint."""
        assert "conda" in NCBIDatasets.INSTALL_HINT
        assert "conda-forge" in NCBIDatasets.INSTALL_HINT
        assert "ncbi-datasets-cli" in NCBIDatasets.INSTALL_HINT


class TestNCBIDatasetsDryRun:
    """Tests for dry run mode."""

    def test_download_batch_dry_run(self):
        """Should return dry run result without executing."""
        ncbi = NCBIDatasets()
        NCBIDatasets._executable_cache["datasets"] = Path("/usr/bin/datasets")

        result = ncbi.download_batch(
            accessions=["GCF_000005845.2"],
            output_dir=Path("/tmp/genomes"),
            dry_run=True,
        )

        assert result.success is True
        assert "[dry-run]" in result.stdout
        assert result.elapsed_seconds == 0.0

    def test_download_genomes_dry_run(self, tmp_path):
        """Should return report with all outcomes marked as success."""
        ncbi = NCBIDatasets()
        NCBIDatasets._executable_cache["datasets"] = Path("/usr/bin/datasets")

        accessions = ["GCF_000005845.2", "GCF_000006765.1"]
        output_dir = tmp_path / "genomes"

        report = ncbi.download_genomes(
            accessions=accessions,
            output_dir=output_dir,
            batch_size=100,
            dry_run=True,
        )

        assert isinstance(report, DownloadReport)
        assert len(report.succeeded) == 2
        assert len(report.failed) == 0


class TestNCBIDatasetsSkipExisting:
    """Tests for skip-if-exists functionality."""

    def test_skip_existing_fna_files(self, tmp_path):
        """Should skip accessions with existing .fna files."""
        ncbi = NCBIDatasets()
        NCBIDatasets._executable_cache["datasets"] = Path("/usr/bin/datasets")

        # Create existing file
        existing_fna = tmp_path / "GCF_000005845.2.fna"
        existing_fna.write_text(">header\nACGT\n")

        accessions = ["GCF_000005845.2", "GCF_000006765.1"]

        report = ncbi.download_genomes(
            accessions=accessions,
            output_dir=tmp_path,
            skip_if_exists=True,
            dry_run=True,  # Use dry run to avoid actual downloads
        )

        # Should only process the non-existing accession
        assert isinstance(report, DownloadReport)
        assert len(report.outcomes) == 1

    def test_skip_existing_fna_gz_files(self, tmp_path):
        """Should skip accessions with existing .fna.gz files."""
        ncbi = NCBIDatasets()
        NCBIDatasets._executable_cache["datasets"] = Path("/usr/bin/datasets")

        # Create existing compressed file
        existing_gz = tmp_path / "GCF_000005845.2.fna.gz"
        existing_gz.write_bytes(b"\x1f\x8b")  # gzip magic bytes

        accessions = ["GCF_000005845.2", "GCF_000006765.1"]

        report = ncbi.download_genomes(
            accessions=accessions,
            output_dir=tmp_path,
            skip_if_exists=True,
            dry_run=True,
        )

        assert isinstance(report, DownloadReport)
        assert len(report.outcomes) == 1

    def test_all_exist_returns_early(self, tmp_path):
        """Should return empty report when all accessions exist."""
        ncbi = NCBIDatasets()
        NCBIDatasets._executable_cache["datasets"] = Path("/usr/bin/datasets")

        # Create all existing files
        (tmp_path / "GCF_000005845.2.fna").write_text(">h\nA\n")
        (tmp_path / "GCF_000006765.1.fna").write_text(">h\nC\n")

        accessions = ["GCF_000005845.2", "GCF_000006765.1"]

        report = ncbi.download_genomes(
            accessions=accessions,
            output_dir=tmp_path,
            skip_if_exists=True,
        )

        assert isinstance(report, DownloadReport)
        assert len(report.outcomes) == 0
        assert report.elapsed_seconds == 0.0


class TestNCBIDatasetsBatching:
    """Tests for batch processing."""

    def test_respects_batch_size(self, tmp_path):
        """Should process in batches of specified size."""
        ncbi = NCBIDatasets()
        NCBIDatasets._executable_cache["datasets"] = Path("/usr/bin/datasets")

        # Create more accessions than batch size
        accessions = [f"GCF_00000000{i}.1" for i in range(10)]

        report = ncbi.download_genomes(
            accessions=accessions,
            output_dir=tmp_path,
            batch_size=3,
            dry_run=True,
        )

        assert isinstance(report, DownloadReport)
        assert len(report.succeeded) == 10


# =============================================================================
# Accession version parsing tests
# =============================================================================


class TestParseAccessionVersion:
    """Tests for _parse_accession_version()."""

    def test_gcf_version_2(self):
        """Should parse GCF accession with version 2."""
        result = _parse_accession_version("GCF_000195955.2")
        assert result == ("GCF_000195955", 2)

    def test_gca_version_1(self):
        """Should parse GCA accession with version 1."""
        result = _parse_accession_version("GCA_000710735.1")
        assert result == ("GCA_000710735", 1)

    def test_high_version_number(self):
        """Should parse multi-digit version numbers."""
        result = _parse_accession_version("GCF_000123456.15")
        assert result == ("GCF_000123456", 15)

    def test_invalid_format_returns_none(self):
        """Should return None for non-accession strings."""
        assert _parse_accession_version("custom_genome") is None

    def test_no_version_returns_none(self):
        """Should return None when no version dot is present."""
        assert _parse_accession_version("GCF_000195955") is None


class TestAlternativeVersions:
    """Tests for _alternative_versions()."""

    def test_version_1_tries_higher(self):
        """For version .1, should try .2, .3, .4 (no lower versions available)."""
        result = _alternative_versions("GCA_000710735.1", max_attempts=3)
        assert result == [
            "GCA_000710735.2",
            "GCA_000710735.3",
            "GCA_000710735.4",
        ]

    def test_version_2_interleaves(self):
        """For version .2, should interleave: .3, .1, .4."""
        result = _alternative_versions("GCF_000195955.2", max_attempts=3)
        assert result == [
            "GCF_000195955.3",
            "GCF_000195955.1",
            "GCF_000195955.4",
        ]

    def test_version_3_includes_lower(self):
        """For version .3 with enough attempts, should include lower versions."""
        result = _alternative_versions("GCF_000195955.3", max_attempts=5)
        # Interleaved: .4, .2, .5, .1, .6
        assert result == [
            "GCF_000195955.4",
            "GCF_000195955.2",
            "GCF_000195955.5",
            "GCF_000195955.1",
            "GCF_000195955.6",
        ]

    def test_max_attempts_limits_output(self):
        """Should return at most max_attempts candidates."""
        result = _alternative_versions("GCA_000710735.1", max_attempts=2)
        assert len(result) == 2

    def test_invalid_accession_returns_empty(self):
        """Should return empty list for non-accession strings."""
        assert _alternative_versions("custom_genome") == []

    def test_no_negative_versions(self):
        """Should not produce version 0 or negative versions."""
        result = _alternative_versions("GCA_000710735.1", max_attempts=5)
        for v in result:
            parsed = _parse_accession_version(v)
            assert parsed is not None
            assert parsed[1] >= 1


# =============================================================================
# DownloadOutcome / DownloadReport tests
# =============================================================================


class TestDownloadOutcome:
    """Tests for DownloadOutcome dataclass."""

    def test_simple_success(self):
        o = DownloadOutcome(accession="GCF_000005845.2", success=True)
        assert o.success
        assert not o.was_version_recovered

    def test_version_recovered(self):
        o = DownloadOutcome(
            accession="GCF_000005845.1",
            success=True,
            resolved_version="GCF_000005845.2",
        )
        assert o.success
        assert o.was_version_recovered

    def test_failure(self):
        o = DownloadOutcome(
            accession="GCF_000005845.1",
            success=False,
            reason="not found in NCBI",
        )
        assert not o.success
        assert not o.was_version_recovered
        assert "not found" in o.reason


class TestDownloadReport:
    """Tests for DownloadReport dataclass."""

    def test_empty_report(self):
        r = DownloadReport()
        assert len(r.succeeded) == 0
        assert len(r.failed) == 0
        assert len(r.recovered) == 0

    def test_mixed_outcomes(self):
        r = DownloadReport(outcomes=[
            DownloadOutcome(accession="GCF_000001.1", success=True),
            DownloadOutcome(accession="GCF_000002.1", success=False, reason="not found"),
            DownloadOutcome(
                accession="GCF_000003.1", success=True,
                resolved_version="GCF_000003.2",
            ),
        ])
        assert len(r.succeeded) == 2
        assert len(r.failed) == 1
        assert len(r.recovered) == 1
        assert r.failed[0].accession == "GCF_000002.1"
        assert r.recovered[0].accession == "GCF_000003.1"

    def test_write_failures_tsv(self, tmp_path):
        """Should write TSV with failed accession details."""
        r = DownloadReport(outcomes=[
            DownloadOutcome(accession="GCF_000001.1", success=True),
            DownloadOutcome(
                accession="GCF_000002.1", success=False,
                reason="not found in NCBI",
            ),
            DownloadOutcome(
                accession="GCA_000003.1", success=False,
                reason="not found in NCBI (tried: GCA_000003.1, GCA_000003.2)",
            ),
        ])

        out = tmp_path / "failures.tsv"
        r.write_failures_tsv(out)

        assert out.exists()
        lines = out.read_text().strip().split("\n")
        assert lines[0] == "accession\treason"
        assert len(lines) == 3  # header + 2 failures
        assert "GCF_000002.1" in lines[1]
        assert "GCA_000003.1" in lines[2]

    def test_write_failures_tsv_no_failures(self, tmp_path):
        """Should not write file when there are no failures."""
        r = DownloadReport(outcomes=[
            DownloadOutcome(accession="GCF_000001.1", success=True),
        ])

        out = tmp_path / "failures.tsv"
        r.write_failures_tsv(out)
        assert not out.exists()
