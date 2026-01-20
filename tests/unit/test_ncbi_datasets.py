"""Unit tests for NCBI Datasets CLI wrapper.

Tests command building logic without requiring the actual tool
to be installed.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metadarkmatter.external.ncbi_datasets import NCBIDatasets


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
        """Should show commands in dry run mode."""
        ncbi = NCBIDatasets()
        NCBIDatasets._executable_cache["datasets"] = Path("/usr/bin/datasets")

        accessions = ["GCF_000005845.2", "GCF_000006765.1"]
        output_dir = tmp_path / "genomes"

        result = ncbi.download_genomes(
            accessions=accessions,
            output_dir=output_dir,
            batch_size=100,
            dry_run=True,
        )

        assert result.success is True
        assert "Would execute" in result.stdout or "[dry-run]" in result.stdout


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

        result = ncbi.download_genomes(
            accessions=accessions,
            output_dir=tmp_path,
            skip_if_exists=True,
            dry_run=True,  # Use dry run to avoid actual downloads
        )

        # Should only process the non-existing accession
        assert result.success is True

    def test_skip_existing_fna_gz_files(self, tmp_path):
        """Should skip accessions with existing .fna.gz files."""
        ncbi = NCBIDatasets()
        NCBIDatasets._executable_cache["datasets"] = Path("/usr/bin/datasets")

        # Create existing compressed file
        existing_gz = tmp_path / "GCF_000005845.2.fna.gz"
        existing_gz.write_bytes(b"\x1f\x8b")  # gzip magic bytes

        accessions = ["GCF_000005845.2", "GCF_000006765.1"]

        result = ncbi.download_genomes(
            accessions=accessions,
            output_dir=tmp_path,
            skip_if_exists=True,
            dry_run=True,
        )

        assert result.success is True

    def test_all_exist_returns_early(self, tmp_path):
        """Should return early when all accessions exist."""
        ncbi = NCBIDatasets()
        NCBIDatasets._executable_cache["datasets"] = Path("/usr/bin/datasets")

        # Create all existing files
        (tmp_path / "GCF_000005845.2.fna").write_text(">h\nA\n")
        (tmp_path / "GCF_000006765.1.fna").write_text(">h\nC\n")

        accessions = ["GCF_000005845.2", "GCF_000006765.1"]

        result = ncbi.download_genomes(
            accessions=accessions,
            output_dir=tmp_path,
            skip_if_exists=True,
        )

        assert result.success is True
        assert "already exist" in result.stdout


class TestNCBIDatasetsBatching:
    """Tests for batch processing."""

    def test_respects_batch_size(self, tmp_path):
        """Should process in batches of specified size."""
        ncbi = NCBIDatasets()
        NCBIDatasets._executable_cache["datasets"] = Path("/usr/bin/datasets")

        # Create more accessions than batch size
        accessions = [f"GCF_00000000{i}.1" for i in range(10)]

        result = ncbi.download_genomes(
            accessions=accessions,
            output_dir=tmp_path,
            batch_size=3,
            dry_run=True,
        )

        assert result.success is True
