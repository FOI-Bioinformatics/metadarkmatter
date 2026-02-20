"""Tests for CLI mapping (util) commands.

Tests the generate-mapping and validate-mapping CLI commands
using Typer's CliRunner with mocked ContigIdMapping operations.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import polars as pl
import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.mapping import app


@pytest.fixture
def runner() -> CliRunner:
    """Create a Typer CLI runner."""
    return CliRunner()


# =============================================================================
# generate-mapping command tests
# =============================================================================


class TestGenerateMapping:
    """Tests for the generate-mapping CLI command."""

    def test_successful_generation(self, runner: CliRunner, tmp_path: Path) -> None:
        """Should generate mapping from genome directory and save to TSV."""
        genomes_dir = tmp_path / "genomes"
        genomes_dir.mkdir()
        # Create a dummy genome file so the directory is not empty
        (genomes_dir / "GCF_000123456.1.fna").write_text(">contig1\nATCG\n")
        output_file = tmp_path / "mapping.tsv"

        mock_mapping = MagicMock()
        mock_mapping.__len__ = MagicMock(return_value=5)

        with patch(
            "metadarkmatter.cli.mapping.ContigIdMapping.from_genome_dir",
            return_value=mock_mapping,
        ) as mock_from_dir:
            result = runner.invoke(
                app,
                [
                    "generate-mapping",
                    "--genomes", str(genomes_dir),
                    "--output", str(output_file),
                ],
            )

        assert result.exit_code == 0
        assert "5" in result.output  # contig count
        mock_from_dir.assert_called_once()
        mock_mapping.to_tsv.assert_called_once_with(output_file)

    def test_custom_pattern(self, runner: CliRunner, tmp_path: Path) -> None:
        """Should pass custom glob pattern to ContigIdMapping."""
        genomes_dir = tmp_path / "genomes"
        genomes_dir.mkdir()
        (genomes_dir / "genome.fa").write_text(">c1\nATCG\n")
        output_file = tmp_path / "mapping.tsv"

        mock_mapping = MagicMock()
        mock_mapping.__len__ = MagicMock(return_value=3)

        with patch(
            "metadarkmatter.cli.mapping.ContigIdMapping.from_genome_dir",
            return_value=mock_mapping,
        ) as mock_from_dir:
            result = runner.invoke(
                app,
                [
                    "generate-mapping",
                    "--genomes", str(genomes_dir),
                    "--output", str(output_file),
                    "--pattern", "*.fa",
                ],
            )

        assert result.exit_code == 0
        mock_from_dir.assert_called_once_with(genomes_dir, pattern="*.fa")

    def test_quiet_mode(self, runner: CliRunner, tmp_path: Path) -> None:
        """Should suppress output when quiet flag is set."""
        genomes_dir = tmp_path / "genomes"
        genomes_dir.mkdir()
        (genomes_dir / "GCF_000123456.1.fna").write_text(">c\nA\n")
        output_file = tmp_path / "mapping.tsv"

        mock_mapping = MagicMock()
        mock_mapping.__len__ = MagicMock(return_value=1)

        with patch(
            "metadarkmatter.cli.mapping.ContigIdMapping.from_genome_dir",
            return_value=mock_mapping,
        ):
            result = runner.invoke(
                app,
                [
                    "generate-mapping",
                    "--genomes", str(genomes_dir),
                    "--output", str(output_file),
                    "--quiet",
                ],
            )

        assert result.exit_code == 0

    def test_file_not_found_error(self, runner: CliRunner, tmp_path: Path) -> None:
        """Should handle FileNotFoundError with exit code 1."""
        genomes_dir = tmp_path / "genomes"
        genomes_dir.mkdir()
        (genomes_dir / "dummy.fna").write_text(">c\nA\n")
        output_file = tmp_path / "mapping.tsv"

        with patch(
            "metadarkmatter.cli.mapping.ContigIdMapping.from_genome_dir",
            side_effect=FileNotFoundError("No genome files found"),
        ):
            result = runner.invoke(
                app,
                [
                    "generate-mapping",
                    "--genomes", str(genomes_dir),
                    "--output", str(output_file),
                ],
            )

        assert result.exit_code == 1
        assert "Error" in result.output

    def test_value_error(self, runner: CliRunner, tmp_path: Path) -> None:
        """Should handle ValueError with exit code 1."""
        genomes_dir = tmp_path / "genomes"
        genomes_dir.mkdir()
        (genomes_dir / "dummy.fna").write_text(">c\nA\n")
        output_file = tmp_path / "mapping.tsv"

        with patch(
            "metadarkmatter.cli.mapping.ContigIdMapping.from_genome_dir",
            side_effect=ValueError("No genome files matched pattern"),
        ):
            result = runner.invoke(
                app,
                [
                    "generate-mapping",
                    "--genomes", str(genomes_dir),
                    "--output", str(output_file),
                ],
            )

        assert result.exit_code == 1
        assert "Error" in result.output


# =============================================================================
# validate-mapping command tests
# =============================================================================


class TestValidateMapping:
    """Tests for the validate-mapping CLI command."""

    def _create_mapping_tsv(self, path: Path) -> None:
        """Helper to create a valid mapping TSV file."""
        df = pl.DataFrame({
            "original_contig_id": ["contig_1", "contig_2", "contig_3"],
            "genome_accession": ["GCF_000123456.1", "GCF_000123456.1", "GCF_000789012.1"],
        })
        df.write_csv(path, separator="\t")

    def test_validate_valid_mapping(self, runner: CliRunner, tmp_path: Path) -> None:
        """Should validate a properly formatted mapping file."""
        mapping_file = tmp_path / "mapping.tsv"
        self._create_mapping_tsv(mapping_file)

        result = runner.invoke(app, ["validate-mapping", str(mapping_file)])

        assert result.exit_code == 0
        assert "Valid mapping file" in result.output
        assert "3" in result.output  # 3 entries

    def test_validate_shows_sample_entries(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should display sample entries from the mapping."""
        mapping_file = tmp_path / "mapping.tsv"
        self._create_mapping_tsv(mapping_file)

        result = runner.invoke(app, ["validate-mapping", str(mapping_file)])

        assert result.exit_code == 0
        assert "Sample entries" in result.output

    def test_validate_with_blast_file_full_coverage(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should report full coverage when all BLAST IDs are mapped."""
        mapping_file = tmp_path / "mapping.tsv"
        self._create_mapping_tsv(mapping_file)

        # Create a BLAST-format file with mapped subject IDs
        blast_file = tmp_path / "results.blast.tsv"
        blast_lines = [
            "read_1\tcontig_1\t99.0\t150\t1\t0\t1\t150\t100\t250\t1e-50\t200",
            "read_2\tcontig_2\t98.0\t150\t2\t0\t1\t150\t200\t350\t1e-45\t190",
        ]
        blast_file.write_text("\n".join(blast_lines) + "\n")

        result = runner.invoke(
            app,
            ["validate-mapping", str(mapping_file), "--blast", str(blast_file)],
        )

        assert result.exit_code == 0
        assert "100.0%" in result.output
        assert "Mapped: 2" in result.output

    def test_validate_with_blast_file_partial_coverage(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should warn about low coverage when some BLAST IDs are unmapped."""
        mapping_file = tmp_path / "mapping.tsv"
        self._create_mapping_tsv(mapping_file)

        # Create a BLAST file with some unmapped subject IDs
        blast_file = tmp_path / "results.blast.tsv"
        blast_lines = [
            "read_1\tcontig_1\t99.0\t150\t1\t0\t1\t150\t100\t250\t1e-50\t200",
            "read_2\tunknown_id\t98.0\t150\t2\t0\t1\t150\t200\t350\t1e-45\t190",
        ]
        blast_file.write_text("\n".join(blast_lines) + "\n")

        result = runner.invoke(
            app,
            ["validate-mapping", str(mapping_file), "--blast", str(blast_file)],
        )

        assert result.exit_code == 0
        assert "50.0%" in result.output
        assert "Warning" in result.output  # Low coverage warning

    def test_validate_shows_unmapped_ids(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should list unmapped IDs when count is small."""
        mapping_file = tmp_path / "mapping.tsv"
        self._create_mapping_tsv(mapping_file)

        blast_file = tmp_path / "results.blast.tsv"
        blast_lines = [
            "read_1\tcontig_1\t99.0\t150\t1\t0\t1\t150\t100\t250\t1e-50\t200",
            "read_2\tmissing_1\t98.0\t150\t2\t0\t1\t150\t200\t350\t1e-45\t190",
            "read_3\tmissing_2\t97.0\t150\t3\t0\t1\t150\t300\t450\t1e-40\t180",
        ]
        blast_file.write_text("\n".join(blast_lines) + "\n")

        result = runner.invoke(
            app,
            ["validate-mapping", str(mapping_file), "--blast", str(blast_file)],
        )

        assert result.exit_code == 0
        assert "Unmapped IDs" in result.output

    def test_validate_invalid_mapping_file(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should handle errors from malformed mapping files."""
        mapping_file = tmp_path / "bad_mapping.tsv"
        # Write a file with wrong column names
        mapping_file.write_text("wrong_col\tanother_col\nval1\tval2\n")

        result = runner.invoke(app, ["validate-mapping", str(mapping_file)])

        assert result.exit_code == 1
        assert "Error" in result.output

    def test_validate_quiet_mode(self, runner: CliRunner, tmp_path: Path) -> None:
        """Should suppress output in quiet mode."""
        mapping_file = tmp_path / "mapping.tsv"
        self._create_mapping_tsv(mapping_file)

        result = runner.invoke(
            app,
            ["validate-mapping", str(mapping_file), "--quiet"],
        )

        assert result.exit_code == 0
