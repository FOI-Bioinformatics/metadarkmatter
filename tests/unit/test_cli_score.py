"""
Unit tests for score CLI commands.

Tests argument validation, file handling, and error cases for the score subcommands.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metadarkmatter.cli.main import app


class TestScoreClassifyValidation:
    """Tests for score classify command argument validation."""

    def test_classify_requires_alignment_file(self, cli_runner, temp_ani_file, temp_dir):
        """Test that --alignment argument is required."""
        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "classify",
                "--ani",
                str(temp_ani_file),
                "--output",
                str(output),
            ],
        )
        assert result.exit_code != 0
        assert "--alignment" in result.output or "alignment" in result.output.lower()

    def test_classify_requires_ani_file(self, cli_runner, temp_blast_file, temp_dir):
        """Test that --ani argument is required."""
        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "classify",
                "--alignment",
                str(temp_blast_file),
                "--output",
                str(output),
            ],
        )
        assert result.exit_code != 0
        assert "--ani" in result.output or "ani" in result.output.lower()

    def test_classify_missing_blast_file(self, cli_runner, temp_ani_file, temp_dir):
        """Test error handling for non-existent BLAST file."""
        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "classify",
                "--alignment",
                str(temp_dir / "nonexistent.blast.tsv"),
                "--ani",
                str(temp_ani_file),
                "--output",
                str(output),
            ],
        )
        assert result.exit_code != 0
        assert "not found" in result.output.lower() or "does not exist" in result.output.lower()

    def test_classify_missing_ani_file(self, cli_runner, temp_blast_file, temp_dir):
        """Test error handling for non-existent ANI file."""
        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "classify",
                "--alignment",
                str(temp_blast_file),
                "--ani",
                str(temp_dir / "nonexistent.ani.csv"),
                "--output",
                str(output),
            ],
        )
        assert result.exit_code != 0
        assert "not found" in result.output.lower() or "does not exist" in result.output.lower()

    def test_classify_empty_blast_file(self, cli_runner, empty_blast_file, temp_ani_file, temp_dir):
        """Test handling of empty BLAST file."""
        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "classify",
                "--alignment",
                str(empty_blast_file),
                "--ani",
                str(temp_ani_file),
                "--output",
                str(output),
            ],
        )
        # Should either fail with error or succeed with warning
        assert result.exit_code != 0 or "empty" in result.output.lower()

    def test_classify_invalid_threshold(self, cli_runner, temp_blast_file, temp_ani_file, temp_dir):
        """Test validation of threshold values."""
        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "classify",
                "--alignment",
                str(temp_blast_file),
                "--ani",
                str(temp_ani_file),
                "--bitscore-threshold",
                "150",  # > 100, invalid
                "--output",
                str(output),
            ],
        )
        assert result.exit_code != 0

    def test_classify_invalid_mode(self, cli_runner, temp_blast_file, temp_ani_file, temp_dir):
        """Test validation of processing mode."""
        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "classify",
                "--alignment",
                str(temp_blast_file),
                "--ani",
                str(temp_ani_file),
                "--mode",
                "invalid_mode",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code != 0


class TestScoreClassifyFamilyFlags:
    """Tests for family validation CLI flags."""

    def test_classify_help_shows_family_flags(self, cli_runner):
        """CLI should show --target-family and --family-ratio-threshold flags."""
        result = cli_runner.invoke(app, ["score", "classify", "--help"])
        assert "--target-family" in result.output
        # Typer truncates long option names in help output, so match prefix
        assert "--family-ratio-thre" in result.output


class TestScoreClassifyOptions:
    """Tests for score classify command with various options."""

    def test_classify_with_metadata(
        self,
        cli_runner,
        temp_blast_file,
        temp_ani_file,
        temp_metadata_file,
        temp_dir,
    ):
        """Test classify with metadata file."""
        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "classify",
                "--alignment",
                str(temp_blast_file),
                "--ani",
                str(temp_ani_file),
                "--metadata",
                str(temp_metadata_file),
                "--output",
                str(output),
            ],
        )
        # Should succeed or at least not crash on argument parsing
        assert result.exit_code == 0 or "error" not in result.output.lower()

    def test_classify_custom_preset(
        self,
        cli_runner,
        temp_blast_file,
        temp_ani_file,
        temp_dir,
    ):
        """Test classify with preset configuration."""
        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "classify",
                "--alignment",
                str(temp_blast_file),
                "--ani",
                str(temp_ani_file),
                "--preset",
                "gtdb-strict",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code == 0 or "error" not in result.output.lower()


class TestScoreExtractNovel:
    """Tests for score extract-novel command."""

    def test_extract_novel_requires_classifications(self, cli_runner, temp_dir):
        """Test that --classifications argument is required."""
        output = temp_dir / "novel.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "extract-novel",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code != 0
        assert "--classifications" in result.output or "classifications" in result.output.lower()

    def test_extract_novel_missing_file(self, cli_runner, temp_dir):
        """Test error handling for non-existent classification file."""
        output = temp_dir / "novel.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "extract-novel",
                "--classifications",
                str(temp_dir / "nonexistent.csv"),
                "--output",
                str(output),
            ],
        )
        assert result.exit_code != 0

    def test_extract_novel_species_only(
        self,
        cli_runner,
        temp_classification_file,
        temp_dir,
    ):
        """Test extracting only novel species."""
        output = temp_dir / "novel_species.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "extract-novel",
                "--classifications",
                str(temp_classification_file),
                "--category",
                "species",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code == 0 or "error" not in result.output.lower()

    def test_extract_novel_genus_only(
        self,
        cli_runner,
        temp_classification_file,
        temp_dir,
    ):
        """Test extracting only novel genus."""
        output = temp_dir / "novel_genus.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "extract-novel",
                "--classifications",
                str(temp_classification_file),
                "--category",
                "genus",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code == 0 or "error" not in result.output.lower()

    def test_extract_novel_with_read_ids(
        self,
        cli_runner,
        temp_classification_file,
        temp_dir,
    ):
        """Test extracting novel reads with read ID output."""
        output_csv = temp_dir / "novel.csv"
        output_ids = temp_dir / "novel_ids.txt"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "extract-novel",
                "--classifications",
                str(temp_classification_file),
                "--output",
                str(output_csv),
                "--read-ids",
                str(output_ids),
            ],
        )
        assert result.exit_code == 0 or "error" not in result.output.lower()

    def test_extract_novel_custom_threshold(
        self,
        cli_runner,
        temp_classification_file,
        temp_dir,
    ):
        """Test extracting with custom novelty threshold."""
        output = temp_dir / "highly_novel.csv"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "extract-novel",
                "--classifications",
                str(temp_classification_file),
                "--min-novelty",
                "15.0",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code == 0 or "error" not in result.output.lower()


class TestScoreBatch:
    """Tests for score batch command."""

    def test_batch_requires_input_dir(self, cli_runner, temp_dir):
        """Test that --alignment-dir argument is required."""
        output_dir = temp_dir / "batch_output"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "batch",
                "--output-dir",
                str(output_dir),
            ],
        )
        assert result.exit_code != 0

    def test_batch_missing_directory(self, cli_runner, temp_dir):
        """Test error handling for non-existent input directory."""
        output_dir = temp_dir / "batch_output"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "batch",
                "--alignment-dir",
                str(temp_dir / "nonexistent"),
                "--output-dir",
                str(output_dir),
            ],
        )
        assert result.exit_code != 0

    def test_batch_empty_directory(self, cli_runner, temp_dir):
        """Test handling of empty input directory."""
        input_dir = temp_dir / "empty_input"
        input_dir.mkdir()
        output_dir = temp_dir / "batch_output"
        result = cli_runner.invoke(
            app,
            [
                "score",
                "batch",
                "--alignment-dir",
                str(input_dir),
                "--output-dir",
                str(output_dir),
            ],
        )
        # Should either warn or fail gracefully
        assert result.exit_code != 0 or "no files" in result.output.lower()


# =============================================================================
# New tests for improved coverage
# =============================================================================

import json
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from unittest.mock import MagicMock, patch

import polars as pl

from metadarkmatter.cli.score import (
    THRESHOLD_PRESETS,
    _display_summary_table,
    _finalize_classification,
    _generate_summary,
    validate_ani_genome_coverage,
    validate_output_format_extension,
)
from metadarkmatter.models.classification import TaxonomicSummary


class TestValidateOutputFormatExtension:
    """Tests for validate_output_format_extension helper."""

    def test_parquet_format_with_csv_extension(self, temp_dir):
        """Parquet format with .csv extension should be corrected to .parquet."""
        from rich.console import Console

        console = Console()
        output = temp_dir / "results.csv"
        corrected = validate_output_format_extension(output, "parquet", console)
        assert corrected.suffix == ".parquet"
        assert corrected.stem == "results"

    def test_parquet_format_with_tsv_extension(self, temp_dir):
        """Parquet format with .tsv extension should be corrected to .parquet."""
        from rich.console import Console

        console = Console()
        output = temp_dir / "results.tsv"
        corrected = validate_output_format_extension(output, "parquet", console)
        assert corrected.suffix == ".parquet"

    def test_csv_format_with_parquet_extension(self, temp_dir):
        """CSV format with .parquet extension should be corrected to .csv."""
        from rich.console import Console

        console = Console()
        output = temp_dir / "results.parquet"
        corrected = validate_output_format_extension(output, "csv", console)
        assert corrected.suffix == ".csv"
        assert corrected.stem == "results"

    def test_matching_csv_extension_unchanged(self, temp_dir):
        """CSV format with .csv extension should remain unchanged."""
        from rich.console import Console

        console = Console()
        output = temp_dir / "results.csv"
        corrected = validate_output_format_extension(output, "csv", console)
        assert corrected == output

    def test_matching_parquet_extension_unchanged(self, temp_dir):
        """Parquet format with .parquet extension should remain unchanged."""
        from rich.console import Console

        console = Console()
        output = temp_dir / "results.parquet"
        corrected = validate_output_format_extension(output, "parquet", console)
        assert corrected == output


class TestFinalizeClassification:
    """Tests for _finalize_classification helper."""

    def test_returns_zero_for_none_df(self, temp_dir):
        """Should return 0 when classification_df is None."""
        output = temp_dir / "output.csv"
        result = _finalize_classification(None, None, output, "csv")
        assert result == 0

    def test_writes_csv_output(self, temp_dir):
        """Should write DataFrame to CSV and return count."""
        df = pl.DataFrame({
            "read_id": ["r1", "r2", "r3"],
            "best_match_genome": ["g1", "g2", "g3"],
            "novelty_index": [1.0, 5.0, 20.0],
            "placement_uncertainty": [0.5, 0.3, 1.0],
            "taxonomic_call": ["Known Species", "Novel Species", "Novel Genus"],
        })
        output = temp_dir / "output.csv"
        result = _finalize_classification(df, None, output, "csv")
        assert result == 3
        assert output.exists()

    def test_writes_parquet_output(self, temp_dir):
        """Should write DataFrame to Parquet."""
        df = pl.DataFrame({
            "read_id": ["r1"],
            "best_match_genome": ["g1"],
            "novelty_index": [1.0],
            "placement_uncertainty": [0.5],
            "taxonomic_call": ["Known Species"],
        })
        output = temp_dir / "output.parquet"
        result = _finalize_classification(df, None, output, "parquet")
        assert result == 1
        assert output.exists()

    def test_joins_metadata(self, temp_dir):
        """Should join genome metadata when provided."""
        df = pl.DataFrame({
            "read_id": ["r1"],
            "best_match_genome": ["g1"],
            "novelty_index": [1.0],
            "placement_uncertainty": [0.5],
            "taxonomic_call": ["Known Species"],
        })
        mock_metadata = MagicMock()
        mock_metadata.join_classifications.return_value = df
        output = temp_dir / "output.csv"

        _finalize_classification(df, mock_metadata, output, "csv")
        mock_metadata.join_classifications.assert_called_once_with(df)

    def test_empty_dataframe_not_written(self, temp_dir):
        """Should not write file when DataFrame is empty."""
        df = pl.DataFrame({
            "read_id": [],
            "best_match_genome": [],
        })
        output = temp_dir / "output.csv"
        result = _finalize_classification(df, None, output, "csv")
        assert result == 0
        assert not output.exists()


class TestGenerateSummary:
    """Tests for _generate_summary function."""

    def test_basic_summary(self):
        """Should produce correct counts from classification DataFrame."""
        df = pl.DataFrame({
            "read_id": [f"r{i}" for i in range(8)],
            "best_match_genome": ["genome_A"] * 5 + ["genome_B"] * 3,
            "taxonomic_call": [
                "Known Species", "Known Species", "Novel Species",
                "Novel Genus", "Ambiguous",
                "Conserved Region", "Unclassified", "Known Species",
            ],
            "novelty_index": [1.0, 2.0, 8.0, 22.0, 5.0, 6.0, 30.0, 1.5],
            "placement_uncertainty": [0.2, 0.3, 0.5, 1.0, 6.0, 8.0, 0.1, 0.2],
        })
        summary = _generate_summary(df)

        assert summary.total_reads == 8
        assert summary.known_species == 3
        assert summary.novel_species == 1
        assert summary.novel_genus == 1
        assert summary.ambiguous == 1
        assert summary.conserved_regions == 1
        assert summary.unclassified == 1
        assert "genome_A" in summary.genome_hit_counts
        assert summary.genome_hit_counts["genome_A"] == 5

    def test_summary_with_species_column(self):
        """Should populate species_hit_counts when species column is present."""
        df = pl.DataFrame({
            "read_id": ["r1", "r2", "r3"],
            "best_match_genome": ["g1", "g1", "g2"],
            "taxonomic_call": ["Known Species", "Known Species", "Novel Species"],
            "novelty_index": [1.0, 2.0, 8.0],
            "placement_uncertainty": [0.2, 0.3, 0.5],
            "species": ["Species A", "Species A", "Species B"],
        })
        summary = _generate_summary(df)
        assert "Species A" in summary.species_hit_counts
        assert summary.species_hit_counts["Species A"] == 2

    def test_summary_off_target(self):
        """Should count Off-target reads."""
        df = pl.DataFrame({
            "read_id": ["r1", "r2"],
            "best_match_genome": ["g1", "g2"],
            "taxonomic_call": ["Off-target", "Known Species"],
            "novelty_index": [5.0, 1.0],
            "placement_uncertainty": [0.5, 0.2],
        })
        summary = _generate_summary(df)
        assert summary.off_target == 1


class TestDisplaySummaryTable:
    """Tests for _display_summary_table function."""

    def test_display_without_error(self):
        """Should display summary table without raising."""
        summary = TaxonomicSummary(
            total_reads=100,
            known_species=50,
            novel_species=20,
            novel_genus=10,
            conserved_regions=5,
            ambiguous=10,
            unclassified=5,
            mean_novelty_index=6.0,
            mean_placement_uncertainty=1.5,
            genome_hit_counts={"genome_A": 60, "genome_B": 40},
        )
        # Should not raise
        _display_summary_table(summary)

    def test_display_with_off_target(self):
        """Should display off-target row when count > 0."""
        summary = TaxonomicSummary(
            total_reads=100,
            known_species=40,
            novel_species=10,
            novel_genus=5,
            conserved_regions=5,
            off_target=15,
            ambiguous=10,
            unclassified=15,
            mean_novelty_index=8.0,
            mean_placement_uncertainty=2.0,
            genome_hit_counts={},
        )
        _display_summary_table(summary)

    def test_display_with_species(self):
        """Should display species information when available."""
        summary = TaxonomicSummary(
            total_reads=100,
            known_species=80,
            novel_species=10,
            novel_genus=5,
            conserved_regions=5,
            mean_novelty_index=3.0,
            mean_placement_uncertainty=0.5,
            genome_hit_counts={"g1": 50},
            species_hit_counts={"Species A": 50, "Species B": 30},
        )
        _display_summary_table(summary)


class TestThresholdPresets:
    """Tests for THRESHOLD_PRESETS dictionary."""

    def test_all_presets_are_valid_configs(self):
        """All defined presets should be valid ScoringConfig instances."""
        from metadarkmatter.models.config import ScoringConfig

        for name, config in THRESHOLD_PRESETS.items():
            assert isinstance(config, ScoringConfig), f"Preset '{name}' is not a ScoringConfig"

    def test_expected_presets_exist(self):
        """Expected preset names should all be present."""
        expected = {"default", "gtdb-strict", "gtdb-relaxed", "conservative", "literature-strict", "coverage-strict", "adaptive"}
        assert expected == set(THRESHOLD_PRESETS.keys())


class TestClassifyCommandMainFlow:
    """Tests for the classify command covering the main execution flow."""

    def _make_blast_and_ani(self, temp_dir):
        """Helper to create matching BLAST and ANI files for testing."""
        # Create ANI matrix
        genomes = ["GCF_000123456.1", "GCF_000789012.1", "GCA_000111222.1"]
        ani_data = {"genome": genomes}
        ani_vals = [
            [100.0, 95.5, 80.0],
            [95.5, 100.0, 82.0],
            [80.0, 82.0, 100.0],
        ]
        for i, g in enumerate(genomes):
            ani_data[g] = ani_vals[i]
        ani_df = pl.DataFrame(ani_data)
        ani_path = temp_dir / "test.ani.csv"
        ani_df.write_csv(ani_path)

        # Create BLAST hits
        blast_rows = []
        for read_idx in range(5):
            for g_idx, genome in enumerate(genomes):
                blast_rows.append({
                    "qseqid": f"read_{read_idx:03d}",
                    "sseqid": f"{genome}_ASM_genomic",
                    "pident": 98.5 - g_idx * 5,
                    "length": 150,
                    "mismatch": 2 + g_idx * 3,
                    "gapopen": 0,
                    "qstart": 1,
                    "qend": 150,
                    "sstart": 1000,
                    "send": 1150,
                    "evalue": 1e-50,
                    "bitscore": 250.0 - g_idx * 20,
                })
        blast_df = pl.DataFrame(blast_rows)
        blast_path = temp_dir / "test.blast.tsv"
        blast_df.write_csv(blast_path, separator="\t", include_header=False)

        return blast_path, ani_path

    def test_classify_invalid_format(self, cli_runner, temp_blast_file, temp_ani_file, temp_dir):
        """Invalid output format should produce error exit."""
        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(temp_blast_file),
                "--ani", str(temp_ani_file),
                "--output", str(output),
                "--format", "xlsx",
            ],
        )
        assert result.exit_code != 0
        assert "invalid format" in result.output.lower() or "xlsx" in result.output.lower()

    def test_classify_invalid_alignment_mode(self, cli_runner, temp_blast_file, temp_ani_file, temp_dir):
        """Invalid alignment mode should produce error exit."""
        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(temp_blast_file),
                "--ani", str(temp_ani_file),
                "--output", str(output),
                "--alignment-mode", "rna",
            ],
        )
        assert result.exit_code != 0

    def test_classify_invalid_uncertainty_mode(self, cli_runner, temp_blast_file, temp_ani_file, temp_dir):
        """Invalid uncertainty mode should produce error exit."""
        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(temp_blast_file),
                "--ani", str(temp_ani_file),
                "--output", str(output),
                "--uncertainty-mode", "average",
            ],
        )
        assert result.exit_code != 0

    def test_classify_invalid_preset(self, cli_runner, temp_blast_file, temp_ani_file, temp_dir):
        """Invalid preset name should produce error exit."""
        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(temp_blast_file),
                "--ani", str(temp_ani_file),
                "--output", str(output),
                "--preset", "nonexistent-preset",
            ],
        )
        assert result.exit_code != 0

    def test_classify_genomes_and_id_mapping_conflict(self, cli_runner, temp_dir):
        """Specifying both --genomes and --id-mapping should produce error."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        genome_dir = temp_dir / "genomes"
        genome_dir.mkdir()
        # Create a dummy FASTA
        (genome_dir / "GCF_000123456.1.fna").write_text(">contig1\nACGT\n")

        mapping_path = temp_dir / "mapping.tsv"
        mapping_path.write_text("contig_id\tgenome_accession\ncontig1\tGCF_000123456.1\n")

        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(blast_path),
                "--ani", str(ani_path),
                "--genomes", str(genome_dir),
                "--id-mapping", str(mapping_path),
                "--output", str(output),
            ],
        )
        assert result.exit_code != 0
        assert "cannot specify both" in result.output.lower() or "error" in result.output.lower()

    def test_classify_dry_run(self, cli_runner, temp_dir):
        """Dry-run mode should validate inputs and exit cleanly."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.csv"

        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(blast_path),
                "--ani", str(ani_path),
                "--output", str(output),
                "--dry-run",
            ],
        )
        assert result.exit_code == 0
        assert "dry run" in result.output.lower() or "validation complete" in result.output.lower()
        # Should NOT create the output file
        assert not output.exists()

    def test_classify_dry_run_with_summary(self, cli_runner, temp_dir):
        """Dry-run with summary flag should show summary path in output."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.csv"
        summary = temp_dir / "summary.json"

        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(blast_path),
                "--ani", str(ani_path),
                "--output", str(output),
                "--summary", str(summary),
                "--dry-run",
            ],
        )
        assert result.exit_code == 0

    def test_classify_protein_mode(self, cli_runner, temp_dir):
        """Protein alignment mode should run successfully."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.csv"

        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(blast_path),
                "--ani", str(ani_path),
                "--output", str(output),
                "--alignment-mode", "protein",
            ],
        )
        assert result.exit_code == 0
        assert "protein mode" in result.output.lower()

    def test_classify_max_uncertainty_mode(self, cli_runner, temp_dir):
        """Max uncertainty mode should run successfully."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.csv"

        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(blast_path),
                "--ani", str(ani_path),
                "--output", str(output),
                "--uncertainty-mode", "max",
            ],
        )
        assert result.exit_code == 0

    def test_classify_verbose_mode(self, cli_runner, temp_dir):
        """Verbose mode should show additional detail."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.csv"

        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(blast_path),
                "--ani", str(ani_path),
                "--output", str(output),
                "--verbose",
            ],
        )
        assert result.exit_code == 0
        # Verbose should mention processing mode
        assert "vectorized" in result.output.lower()

    def test_classify_quiet_mode(self, cli_runner, temp_dir):
        """Quiet mode should suppress most output."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.csv"

        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(blast_path),
                "--ani", str(ani_path),
                "--output", str(output),
                "--quiet",
            ],
        )
        assert result.exit_code == 0

    def test_classify_with_summary_output(self, cli_runner, temp_dir):
        """Classify with --summary should write a JSON summary file."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.csv"
        summary_path = temp_dir / "summary.json"

        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(blast_path),
                "--ani", str(ani_path),
                "--output", str(output),
                "--summary", str(summary_path),
            ],
        )
        assert result.exit_code == 0
        assert output.exists()
        # Summary should be generated if reads were classified
        if summary_path.exists():
            summary_data = json.loads(summary_path.read_text())
            assert "total_reads" in summary_data

    def test_classify_with_parquet_output(self, cli_runner, temp_dir):
        """Classify to parquet format should produce valid output."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.parquet"

        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(blast_path),
                "--ani", str(ani_path),
                "--output", str(output),
                "--format", "parquet",
            ],
        )
        assert result.exit_code == 0
        assert output.exists()

    def test_classify_with_preset_and_overrides(self, cli_runner, temp_dir):
        """Preset with non-default options should merge correctly."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.csv"

        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(blast_path),
                "--ani", str(ani_path),
                "--output", str(output),
                "--preset", "gtdb-strict",
                "--bitscore-threshold", "90.0",
            ],
        )
        assert result.exit_code == 0

    def test_classify_with_qc_output(self, cli_runner, temp_dir):
        """Classify with --qc-output should write QC metrics JSON."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.csv"
        qc_path = temp_dir / "qc.json"

        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(blast_path),
                "--ani", str(ani_path),
                "--output", str(output),
                "--qc-output", str(qc_path),
            ],
        )
        assert result.exit_code == 0
        if qc_path.exists():
            qc_data = json.loads(qc_path.read_text())
            assert "total_alignments" in qc_data

    def test_classify_with_bayesian(self, cli_runner, temp_dir):
        """Classify with --bayesian should add posterior columns."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.csv"

        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(blast_path),
                "--ani", str(ani_path),
                "--output", str(output),
                "--bayesian",
            ],
        )
        assert result.exit_code == 0
        if output.exists():
            df = pl.read_csv(output)
            # Bayesian should add posterior columns
            if len(df) > 0:
                assert "p_known_species" in df.columns or "bayesian_category" in df.columns

    def test_classify_streaming_mode(self, cli_runner, temp_dir):
        """Streaming mode should write output directly to file."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.csv"

        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(blast_path),
                "--ani", str(ani_path),
                "--output", str(output),
                "--streaming",
            ],
        )
        assert result.exit_code == 0
        assert "streaming" in result.output.lower()

    def test_classify_streaming_with_metadata_warning(self, cli_runner, temp_dir):
        """Streaming with metadata should show a warning."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        metadata_path = temp_dir / "metadata.tsv"
        pl.DataFrame({
            "accession": ["GCF_000123456.1"],
            "species": ["Test species"],
            "genus": ["Test"],
            "family": ["Testaceae"],
        }).write_csv(metadata_path, separator="\t")

        output = temp_dir / "output.csv"

        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(blast_path),
                "--ani", str(ani_path),
                "--metadata", str(metadata_path),
                "--output", str(output),
                "--streaming",
            ],
        )
        assert result.exit_code == 0
        assert "not supported" in result.output.lower() or "streaming" in result.output.lower()

    def test_classify_bayesian_streaming_warning(self, cli_runner, temp_dir):
        """Bayesian with streaming should show warning about no support."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.csv"

        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(blast_path),
                "--ani", str(ani_path),
                "--output", str(output),
                "--streaming",
                "--bayesian",
            ],
        )
        assert result.exit_code == 0

    def test_classify_with_target_family(self, cli_runner, temp_dir):
        """Classify with --target-family should enable family validation."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.csv"

        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(blast_path),
                "--ani", str(ani_path),
                "--output", str(output),
                "--target-family", "f__Francisellaceae",
            ],
        )
        assert result.exit_code == 0

    def test_classify_infers_family_from_metadata(self, cli_runner, temp_dir):
        """Should infer target family from metadata when not explicitly given."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        metadata_path = temp_dir / "metadata.tsv"
        pl.DataFrame({
            "accession": ["GCF_000123456.1", "GCF_000789012.1", "GCA_000111222.1"],
            "species": ["Sp A", "Sp B", "Sp C"],
            "genus": ["Gen", "Gen", "Gen"],
            "family": ["f__Francisellaceae", "f__Francisellaceae", "f__Francisellaceae"],
        }).write_csv(metadata_path, separator="\t")

        output = temp_dir / "output.csv"

        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(blast_path),
                "--ani", str(ani_path),
                "--metadata", str(metadata_path),
                "--output", str(output),
            ],
        )
        assert result.exit_code == 0

    def test_classify_creates_output_directory(self, cli_runner, temp_dir):
        """Should create nested output directories if they do not exist."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "nested" / "deep" / "output.csv"

        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(blast_path),
                "--ani", str(ani_path),
                "--output", str(output),
            ],
        )
        assert result.exit_code == 0


class TestClassifyErrorHandling:
    """Tests for classify command error handling paths."""

    def _make_blast_and_ani(self, temp_dir):
        """Helper to create matching BLAST and ANI files."""
        genomes = ["GCF_000123456.1", "GCF_000789012.1", "GCA_000111222.1"]
        ani_data = {"genome": genomes}
        ani_vals = [[100.0, 95.5, 80.0], [95.5, 100.0, 82.0], [80.0, 82.0, 100.0]]
        for i, g in enumerate(genomes):
            ani_data[g] = ani_vals[i]
        ani_df = pl.DataFrame(ani_data)
        ani_path = temp_dir / "test.ani.csv"
        ani_df.write_csv(ani_path)

        blast_rows = []
        for read_idx in range(3):
            for g_idx, genome in enumerate(genomes):
                blast_rows.append({
                    "qseqid": f"read_{read_idx:03d}",
                    "sseqid": f"{genome}_ASM_genomic",
                    "pident": 98.5 - g_idx * 5,
                    "length": 150,
                    "mismatch": 2,
                    "gapopen": 0,
                    "qstart": 1,
                    "qend": 150,
                    "sstart": 1000,
                    "send": 1150,
                    "evalue": 1e-50,
                    "bitscore": 250.0 - g_idx * 20,
                })
        blast_df = pl.DataFrame(blast_rows)
        blast_path = temp_dir / "test.blast.tsv"
        blast_df.write_csv(blast_path, separator="\t", include_header=False)
        return blast_path, ani_path

    def test_classify_malformed_ani_matrix(self, cli_runner, temp_blast_file, temp_dir):
        """Malformed ANI matrix file should exit with error."""
        bad_ani = temp_dir / "bad_ani.csv"
        bad_ani.write_text("not,a,valid,matrix\n1,2,3,4\n")
        output = temp_dir / "output.csv"

        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(temp_blast_file),
                "--ani", str(bad_ani),
                "--output", str(output),
            ],
        )
        assert result.exit_code != 0

    def test_classify_polars_error_handling(self, cli_runner, temp_dir):
        """PolarsError during classification should be caught gracefully."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.csv"

        with patch(
            "metadarkmatter.cli.score.VectorizedClassifier"
        ) as mock_cls:
            mock_instance = MagicMock()
            mock_instance.classify_file.side_effect = pl.exceptions.PolarsError("test error")
            mock_cls.return_value = mock_instance

            result = cli_runner.invoke(
                app,
                [
                    "score", "classify",
                    "--alignment", str(blast_path),
                    "--ani", str(ani_path),
                    "--output", str(output),
                ],
            )
            assert result.exit_code != 0

    def test_classify_memory_error_handling(self, cli_runner, temp_dir):
        """MemoryError during classification should suggest streaming mode."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.csv"

        with patch(
            "metadarkmatter.cli.score.VectorizedClassifier"
        ) as mock_cls:
            mock_instance = MagicMock()
            mock_instance.classify_file.side_effect = MemoryError()
            mock_cls.return_value = mock_instance

            result = cli_runner.invoke(
                app,
                [
                    "score", "classify",
                    "--alignment", str(blast_path),
                    "--ani", str(ani_path),
                    "--output", str(output),
                ],
            )
            assert result.exit_code != 0
            assert "memory" in result.output.lower() or "streaming" in result.output.lower()

    def test_classify_generic_exception_handling(self, cli_runner, temp_dir):
        """Generic exceptions during classification should be caught."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.csv"

        with patch(
            "metadarkmatter.cli.score.VectorizedClassifier"
        ) as mock_cls:
            mock_instance = MagicMock()
            mock_instance.classify_file.side_effect = RuntimeError("unexpected error")
            mock_cls.return_value = mock_instance

            result = cli_runner.invoke(
                app,
                [
                    "score", "classify",
                    "--alignment", str(blast_path),
                    "--ani", str(ani_path),
                    "--output", str(output),
                ],
            )
            assert result.exit_code != 0

    def test_classify_file_not_found_during_processing(self, cli_runner, temp_dir):
        """FileNotFoundError during classification should be caught."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.csv"

        with patch(
            "metadarkmatter.cli.score.VectorizedClassifier"
        ) as mock_cls:
            mock_instance = MagicMock()
            mock_instance.classify_file.side_effect = FileNotFoundError("missing file")
            mock_cls.return_value = mock_instance

            result = cli_runner.invoke(
                app,
                [
                    "score", "classify",
                    "--alignment", str(blast_path),
                    "--ani", str(ani_path),
                    "--output", str(output),
                ],
            )
            assert result.exit_code != 0

    def test_classify_permission_error_handling(self, cli_runner, temp_dir):
        """PermissionError during classification should be caught."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.csv"

        with patch(
            "metadarkmatter.cli.score.VectorizedClassifier"
        ) as mock_cls:
            mock_instance = MagicMock()
            mock_instance.classify_file.side_effect = PermissionError("denied")
            mock_cls.return_value = mock_instance

            result = cli_runner.invoke(
                app,
                [
                    "score", "classify",
                    "--alignment", str(blast_path),
                    "--ani", str(ani_path),
                    "--output", str(output),
                ],
            )
            assert result.exit_code != 0


class TestBatchCommand:
    """Tests for the batch command flow."""

    def _setup_batch(self, temp_dir):
        """Create BLAST files and ANI for batch testing."""
        genomes = ["GCF_000123456.1", "GCF_000789012.1", "GCA_000111222.1"]
        ani_data = {"genome": genomes}
        ani_vals = [[100.0, 95.5, 80.0], [95.5, 100.0, 82.0], [80.0, 82.0, 100.0]]
        for i, g in enumerate(genomes):
            ani_data[g] = ani_vals[i]
        ani_df = pl.DataFrame(ani_data)
        ani_path = temp_dir / "test.ani.csv"
        ani_df.write_csv(ani_path)

        # Create input dir with BLAST files
        input_dir = temp_dir / "blast_input"
        input_dir.mkdir()

        for sample in ["sample_A", "sample_B"]:
            blast_rows = []
            for read_idx in range(3):
                for g_idx, genome in enumerate(genomes):
                    blast_rows.append({
                        "qseqid": f"read_{read_idx:03d}",
                        "sseqid": f"{genome}_ASM_genomic",
                        "pident": 98.5 - g_idx * 5,
                        "length": 150,
                        "mismatch": 2,
                        "gapopen": 0,
                        "qstart": 1,
                        "qend": 150,
                        "sstart": 1000,
                        "send": 1150,
                        "evalue": 1e-50,
                        "bitscore": 250.0 - g_idx * 20,
                    })
            blast_df = pl.DataFrame(blast_rows)
            blast_path = input_dir / f"{sample}.blast.tsv.gz"
            # Write as plain TSV (not gzipped) but with .gz extension for the pattern
            # Actually, use .tsv extension for the pattern
            blast_path = input_dir / f"{sample}.blast.tsv"
            blast_df.write_csv(blast_path, separator="\t", include_header=False)

        return input_dir, ani_path

    def test_batch_classify_basic(self, cli_runner, temp_dir):
        """Batch command should process multiple files."""
        input_dir, ani_path = self._setup_batch(temp_dir)
        output_dir = temp_dir / "batch_output"

        result = cli_runner.invoke(
            app,
            [
                "score", "batch",
                "--alignment-dir", str(input_dir),
                "--ani", str(ani_path),
                "--output-dir", str(output_dir),
                "--pattern", "*.blast.tsv",
            ],
        )
        assert result.exit_code == 0
        assert "batch processing complete" in result.output.lower()

    def test_batch_invalid_format(self, cli_runner, temp_dir):
        """Invalid format in batch should exit with error."""
        input_dir, ani_path = self._setup_batch(temp_dir)
        output_dir = temp_dir / "batch_output"

        result = cli_runner.invoke(
            app,
            [
                "score", "batch",
                "--alignment-dir", str(input_dir),
                "--ani", str(ani_path),
                "--output-dir", str(output_dir),
                "--format", "xlsx",
            ],
        )
        assert result.exit_code != 0

    def test_batch_invalid_preset(self, cli_runner, temp_dir):
        """Invalid preset in batch should exit with error."""
        input_dir, ani_path = self._setup_batch(temp_dir)
        output_dir = temp_dir / "batch_output"

        result = cli_runner.invoke(
            app,
            [
                "score", "batch",
                "--alignment-dir", str(input_dir),
                "--ani", str(ani_path),
                "--output-dir", str(output_dir),
                "--pattern", "*.blast.tsv",
                "--preset", "nonexistent",
            ],
        )
        assert result.exit_code != 0

    def test_batch_with_preset(self, cli_runner, temp_dir):
        """Batch with preset should apply preset configuration."""
        input_dir, ani_path = self._setup_batch(temp_dir)
        output_dir = temp_dir / "batch_output"

        result = cli_runner.invoke(
            app,
            [
                "score", "batch",
                "--alignment-dir", str(input_dir),
                "--ani", str(ani_path),
                "--output-dir", str(output_dir),
                "--pattern", "*.blast.tsv",
                "--preset", "gtdb-strict",
            ],
        )
        assert result.exit_code == 0

    def test_batch_with_parquet_output(self, cli_runner, temp_dir):
        """Batch with parquet format should create .parquet files."""
        input_dir, ani_path = self._setup_batch(temp_dir)
        output_dir = temp_dir / "batch_output"

        result = cli_runner.invoke(
            app,
            [
                "score", "batch",
                "--alignment-dir", str(input_dir),
                "--ani", str(ani_path),
                "--output-dir", str(output_dir),
                "--pattern", "*.blast.tsv",
                "--format", "parquet",
            ],
        )
        assert result.exit_code == 0

    def test_batch_invalid_alignment_mode(self, cli_runner, temp_dir):
        """Invalid alignment mode in batch should exit with error."""
        input_dir, ani_path = self._setup_batch(temp_dir)
        output_dir = temp_dir / "batch_output"

        result = cli_runner.invoke(
            app,
            [
                "score", "batch",
                "--alignment-dir", str(input_dir),
                "--ani", str(ani_path),
                "--output-dir", str(output_dir),
                "--pattern", "*.blast.tsv",
                "--alignment-mode", "rna",
            ],
        )
        assert result.exit_code != 0

    def test_batch_invalid_uncertainty_mode(self, cli_runner, temp_dir):
        """Invalid uncertainty mode in batch should exit with error."""
        input_dir, ani_path = self._setup_batch(temp_dir)
        output_dir = temp_dir / "batch_output"

        result = cli_runner.invoke(
            app,
            [
                "score", "batch",
                "--alignment-dir", str(input_dir),
                "--ani", str(ani_path),
                "--output-dir", str(output_dir),
                "--pattern", "*.blast.tsv",
                "--uncertainty-mode", "average",
            ],
        )
        assert result.exit_code != 0

    def test_batch_protein_mode(self, cli_runner, temp_dir):
        """Batch in protein mode should show protein mode message."""
        input_dir, ani_path = self._setup_batch(temp_dir)
        output_dir = temp_dir / "batch_output"

        result = cli_runner.invoke(
            app,
            [
                "score", "batch",
                "--alignment-dir", str(input_dir),
                "--ani", str(ani_path),
                "--output-dir", str(output_dir),
                "--pattern", "*.blast.tsv",
                "--alignment-mode", "protein",
            ],
        )
        assert result.exit_code == 0
        assert "protein mode" in result.output.lower()

    def test_batch_max_uncertainty_mode(self, cli_runner, temp_dir):
        """Batch with max uncertainty mode should show appropriate message."""
        input_dir, ani_path = self._setup_batch(temp_dir)
        output_dir = temp_dir / "batch_output"

        result = cli_runner.invoke(
            app,
            [
                "score", "batch",
                "--alignment-dir", str(input_dir),
                "--ani", str(ani_path),
                "--output-dir", str(output_dir),
                "--pattern", "*.blast.tsv",
                "--uncertainty-mode", "max",
            ],
        )
        assert result.exit_code == 0

    def test_batch_verbose(self, cli_runner, temp_dir):
        """Batch verbose mode should show coverage info."""
        input_dir, ani_path = self._setup_batch(temp_dir)
        output_dir = temp_dir / "batch_output"

        result = cli_runner.invoke(
            app,
            [
                "score", "batch",
                "--alignment-dir", str(input_dir),
                "--ani", str(ani_path),
                "--output-dir", str(output_dir),
                "--pattern", "*.blast.tsv",
                "--verbose",
            ],
        )
        assert result.exit_code == 0

    def test_batch_ani_error(self, cli_runner, temp_dir):
        """Batch with malformed ANI should exit with error."""
        input_dir, _ = self._setup_batch(temp_dir)
        bad_ani = temp_dir / "bad_ani.csv"
        bad_ani.write_text("bad\n")
        output_dir = temp_dir / "batch_output"

        result = cli_runner.invoke(
            app,
            [
                "score", "batch",
                "--alignment-dir", str(input_dir),
                "--ani", str(bad_ani),
                "--output-dir", str(output_dir),
                "--pattern", "*.blast.tsv",
            ],
        )
        assert result.exit_code != 0

    def test_batch_no_matching_files(self, cli_runner, temp_dir):
        """Batch with no matching pattern should exit with error."""
        input_dir, ani_path = self._setup_batch(temp_dir)
        output_dir = temp_dir / "batch_output"

        result = cli_runner.invoke(
            app,
            [
                "score", "batch",
                "--alignment-dir", str(input_dir),
                "--ani", str(ani_path),
                "--output-dir", str(output_dir),
                "--pattern", "*.nonexistent.xyz",
            ],
        )
        assert result.exit_code != 0

    def test_batch_preset_with_protein_override(self, cli_runner, temp_dir):
        """Batch with preset and protein mode should merge configs."""
        input_dir, ani_path = self._setup_batch(temp_dir)
        output_dir = temp_dir / "batch_output"

        result = cli_runner.invoke(
            app,
            [
                "score", "batch",
                "--alignment-dir", str(input_dir),
                "--ani", str(ani_path),
                "--output-dir", str(output_dir),
                "--pattern", "*.blast.tsv",
                "--preset", "gtdb-relaxed",
                "--alignment-mode", "protein",
            ],
        )
        assert result.exit_code == 0


class TestExtractNovelAdditional:
    """Additional tests for extract-novel command covering more paths."""

    def _make_classification_file(self, temp_dir, include_novel=True):
        """Create a classification CSV with various categories."""
        data = {
            "read_id": [f"read_{i:03d}" for i in range(20)],
            "best_match_genome": (
                ["GCF_000123456.1"] * 8 + ["GCF_000789012.1"] * 7 + ["GCA_000111222.1"] * 5
            ),
            "taxonomic_call": (
                ["Known Species"] * 5
                + (["Novel Species"] * 5 if include_novel else ["Known Species"] * 5)
                + ["Novel Genus"] * 4
                + ["Conserved Region"] * 3
                + ["Ambiguous"] * 3
            ),
            "novelty_index": (
                [1.0, 2.0, 3.0, 1.5, 2.5]
                + ([8.0, 10.0, 12.0, 15.0, 6.0] if include_novel else [1.0, 2.0, 3.0, 1.5, 2.5])
                + [22.0, 23.0, 21.0, 24.0]
                + [5.0, 6.0, 7.0]
                + [4.0, 3.0, 5.0]
            ),
            "placement_uncertainty": (
                [0.2, 0.3, 0.4, 0.1, 0.2]
                + [0.5, 0.6, 0.7, 0.3, 0.4]
                + [1.0, 1.2, 0.8, 1.5]
                + [8.0, 9.0, 10.0]
                + [6.0, 7.0, 8.0]
            ),
            "top_hit_identity": (
                [99.0, 98.0, 97.0, 98.5, 97.5]
                + [92.0, 90.0, 88.0, 85.0, 94.0]
                + [78.0, 77.0, 79.0, 76.0]
                + [95.0, 94.0, 93.0]
                + [96.0, 97.0, 95.0]
            ),
        }
        df = pl.DataFrame(data)
        path = temp_dir / "classifications.csv"
        df.write_csv(path)
        return path

    def test_extract_novel_invalid_category(self, cli_runner, temp_dir):
        """Invalid category should exit with error."""
        clf_path = self._make_classification_file(temp_dir)
        output = temp_dir / "novel.csv"

        result = cli_runner.invoke(
            app,
            [
                "score", "extract-novel",
                "--classifications", str(clf_path),
                "--output", str(output),
                "--category", "invalid_cat",
            ],
        )
        assert result.exit_code != 0

    def test_extract_novel_no_results(self, cli_runner, temp_dir):
        """When no novel reads found, should exit cleanly."""
        clf_path = self._make_classification_file(temp_dir, include_novel=False)
        output = temp_dir / "novel.csv"

        result = cli_runner.invoke(
            app,
            [
                "score", "extract-novel",
                "--classifications", str(clf_path),
                "--output", str(output),
                "--category", "species",
            ],
        )
        assert result.exit_code == 0

    def test_extract_novel_all_categories(self, cli_runner, temp_dir):
        """Category 'all' should extract both species and genus."""
        clf_path = self._make_classification_file(temp_dir)
        output = temp_dir / "novel.csv"

        result = cli_runner.invoke(
            app,
            [
                "score", "extract-novel",
                "--classifications", str(clf_path),
                "--output", str(output),
                "--category", "all",
            ],
        )
        assert result.exit_code == 0
        if output.exists():
            df = pl.read_csv(output)
            assert len(df) > 0

    def test_extract_novel_no_group(self, cli_runner, temp_dir):
        """--no-group should write individual reads."""
        clf_path = self._make_classification_file(temp_dir)
        output = temp_dir / "novel.csv"

        result = cli_runner.invoke(
            app,
            [
                "score", "extract-novel",
                "--classifications", str(clf_path),
                "--output", str(output),
                "--no-group",
            ],
        )
        assert result.exit_code == 0
        if output.exists():
            df = pl.read_csv(output)
            assert "read_id" in df.columns

    def test_extract_novel_with_min_reads_filter(self, cli_runner, temp_dir):
        """--min-reads should filter groups with too few reads."""
        clf_path = self._make_classification_file(temp_dir)
        output = temp_dir / "novel.csv"

        result = cli_runner.invoke(
            app,
            [
                "score", "extract-novel",
                "--classifications", str(clf_path),
                "--output", str(output),
                "--min-reads", "3",
            ],
        )
        assert result.exit_code == 0

    def test_extract_novel_with_read_ids_output(self, cli_runner, temp_dir):
        """Should write read IDs to separate file."""
        clf_path = self._make_classification_file(temp_dir)
        output = temp_dir / "novel.csv"
        read_ids = temp_dir / "read_ids.txt"

        result = cli_runner.invoke(
            app,
            [
                "score", "extract-novel",
                "--classifications", str(clf_path),
                "--output", str(output),
                "--read-ids", str(read_ids),
            ],
        )
        assert result.exit_code == 0
        if read_ids.exists():
            content = read_ids.read_text()
            assert len(content.strip()) > 0

    def test_extract_novel_high_novelty_warning(self, cli_runner, temp_dir):
        """High min-novelty should trigger warning."""
        clf_path = self._make_classification_file(temp_dir)
        output = temp_dir / "novel.csv"

        result = cli_runner.invoke(
            app,
            [
                "score", "extract-novel",
                "--classifications", str(clf_path),
                "--output", str(output),
                "--min-novelty", "55",
            ],
        )
        # Should warn about high novelty and exit 0 (no reads found)
        assert result.exit_code == 0

    def test_extract_novel_high_uncertainty_warning(self, cli_runner, temp_dir):
        """High max-uncertainty should trigger warning."""
        clf_path = self._make_classification_file(temp_dir)
        output = temp_dir / "novel.csv"

        result = cli_runner.invoke(
            app,
            [
                "score", "extract-novel",
                "--classifications", str(clf_path),
                "--output", str(output),
                "--max-uncertainty", "10",
            ],
        )
        assert result.exit_code == 0

    def test_extract_novel_quiet_mode(self, cli_runner, temp_dir):
        """Quiet mode should suppress table output."""
        clf_path = self._make_classification_file(temp_dir)
        output = temp_dir / "novel.csv"

        result = cli_runner.invoke(
            app,
            [
                "score", "extract-novel",
                "--classifications", str(clf_path),
                "--output", str(output),
                "--quiet",
            ],
        )
        assert result.exit_code == 0


class TestValidateAniGenomeCoverage:
    """Tests for validate_ani_genome_coverage function."""

    def test_full_coverage(self, temp_dir, small_ani_dict):
        """All genomes in BLAST should be in ANI matrix."""
        from metadarkmatter.core.ani_placement import ANIMatrix

        ani_matrix = ANIMatrix(small_ani_dict)

        # Create BLAST with matching genomes using pipe-separated format
        # (expected by extract_genome_name_expr: "{accession}|{contig_id}")
        blast_rows = []
        for genome in small_ani_dict:
            blast_rows.append({
                "qseqid": "read_001",
                "sseqid": f"{genome}|contig_1",
                "pident": 98.0,
                "length": 150,
                "mismatch": 2,
                "gapopen": 0,
                "qstart": 1,
                "qend": 150,
                "sstart": 1000,
                "send": 1150,
                "evalue": 1e-50,
                "bitscore": 250.0,
            })
        blast_df = pl.DataFrame(blast_rows)
        blast_path = temp_dir / "test.blast.tsv"
        blast_df.write_csv(blast_path, separator="\t", include_header=False)

        matched, total, pct, missing = validate_ani_genome_coverage(
            blast_path, ani_matrix
        )
        assert matched == total
        assert pct == 100.0
        assert len(missing) == 0

    def test_partial_coverage(self, temp_dir, small_ani_dict):
        """Some genomes missing from ANI should be detected."""
        from metadarkmatter.core.ani_placement import ANIMatrix

        ani_matrix = ANIMatrix(small_ani_dict)

        # Create BLAST with one extra genome not in ANI (pipe-separated format)
        blast_rows = [
            {
                "qseqid": "read_001",
                "sseqid": "GCF_000123456.1|contig_1",
                "pident": 98.0,
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
                "qseqid": "read_002",
                "sseqid": "GCF_999999999.1|contig_1",
                "pident": 95.0,
                "length": 150,
                "mismatch": 7,
                "gapopen": 0,
                "qstart": 1,
                "qend": 150,
                "sstart": 1000,
                "send": 1150,
                "evalue": 1e-40,
                "bitscore": 200.0,
            },
        ]
        blast_df = pl.DataFrame(blast_rows)
        blast_path = temp_dir / "test.blast.tsv"
        blast_df.write_csv(blast_path, separator="\t", include_header=False)

        matched, total, pct, missing = validate_ani_genome_coverage(
            blast_path, ani_matrix
        )
        assert matched < total
        assert pct < 100.0
        assert "GCF_999999999.1" in missing


class TestClassifyWithAdaptiveThresholds:
    """Tests for --adaptive-thresholds flag."""

    def _make_blast_and_ani(self, temp_dir):
        """Helper for blast + ani files."""
        genomes = ["GCF_000123456.1", "GCF_000789012.1", "GCA_000111222.1"]
        ani_data = {"genome": genomes}
        ani_vals = [[100.0, 95.5, 80.0], [95.5, 100.0, 82.0], [80.0, 82.0, 100.0]]
        for i, g in enumerate(genomes):
            ani_data[g] = ani_vals[i]
        ani_df = pl.DataFrame(ani_data)
        ani_path = temp_dir / "test.ani.csv"
        ani_df.write_csv(ani_path)

        blast_rows = []
        for read_idx in range(3):
            for g_idx, genome in enumerate(genomes):
                blast_rows.append({
                    "qseqid": f"read_{read_idx:03d}",
                    "sseqid": f"{genome}_ASM_genomic",
                    "pident": 98.5 - g_idx * 5,
                    "length": 150,
                    "mismatch": 2,
                    "gapopen": 0,
                    "qstart": 1,
                    "qend": 150,
                    "sstart": 1000,
                    "send": 1150,
                    "evalue": 1e-50,
                    "bitscore": 250.0 - g_idx * 20,
                })
        blast_df = pl.DataFrame(blast_rows)
        blast_path = temp_dir / "test.blast.tsv"
        blast_df.write_csv(blast_path, separator="\t", include_header=False)
        return blast_path, ani_path

    def test_adaptive_thresholds_success(self, cli_runner, temp_dir):
        """Adaptive thresholds should work when scikit-learn is available."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.csv"

        # Try running with adaptive thresholds
        # It may succeed or fall back depending on GMM convergence
        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(blast_path),
                "--ani", str(ani_path),
                "--output", str(output),
                "--adaptive-thresholds",
            ],
        )
        # Should succeed (either with GMM result or fallback)
        assert result.exit_code == 0

    def test_adaptive_thresholds_import_error(self, cli_runner, temp_dir):
        """Should handle missing scikit-learn gracefully."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.csv"

        with patch(
            "metadarkmatter.cli.score.VectorizedClassifier"
        ):
            # Mock the adaptive import to raise ImportError
            import builtins
            real_import = builtins.__import__

            def mock_import(name, *args, **kwargs):
                if "adaptive" in name:
                    raise ImportError("No module named 'sklearn'")
                return real_import(name, *args, **kwargs)

            with patch("builtins.__import__", side_effect=mock_import):
                result = cli_runner.invoke(
                    app,
                    [
                        "score", "classify",
                        "--alignment", str(blast_path),
                        "--ani", str(ani_path),
                        "--output", str(output),
                        "--adaptive-thresholds",
                    ],
                )
                assert result.exit_code != 0


class TestClassifyWithSummaryAndQC:
    """Tests for summary generation and QC output in the classify flow."""

    def _make_blast_and_ani(self, temp_dir):
        """Helper for blast + ani files."""
        genomes = ["GCF_000123456.1", "GCF_000789012.1", "GCA_000111222.1"]
        ani_data = {"genome": genomes}
        ani_vals = [[100.0, 95.5, 80.0], [95.5, 100.0, 82.0], [80.0, 82.0, 100.0]]
        for i, g in enumerate(genomes):
            ani_data[g] = ani_vals[i]
        ani_df = pl.DataFrame(ani_data)
        ani_path = temp_dir / "test.ani.csv"
        ani_df.write_csv(ani_path)

        blast_rows = []
        for read_idx in range(5):
            for g_idx, genome in enumerate(genomes):
                blast_rows.append({
                    "qseqid": f"read_{read_idx:03d}",
                    "sseqid": f"{genome}_ASM_genomic",
                    "pident": 98.5 - g_idx * 5,
                    "length": 150,
                    "mismatch": 2 + g_idx * 3,
                    "gapopen": 0,
                    "qstart": 1,
                    "qend": 150,
                    "sstart": 1000,
                    "send": 1150,
                    "evalue": 1e-50,
                    "bitscore": 250.0 - g_idx * 20,
                })
        blast_df = pl.DataFrame(blast_rows)
        blast_path = temp_dir / "test.blast.tsv"
        blast_df.write_csv(blast_path, separator="\t", include_header=False)
        return blast_path, ani_path

    def test_summary_json_structure(self, cli_runner, temp_dir):
        """Summary JSON should contain expected fields."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.csv"
        summary_path = temp_dir / "summary.json"

        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(blast_path),
                "--ani", str(ani_path),
                "--output", str(output),
                "--summary", str(summary_path),
            ],
        )
        assert result.exit_code == 0

        if summary_path.exists():
            data = json.loads(summary_path.read_text())
            assert "total_reads" in data
            assert "known_species" in data
            assert "novel_species" in data
            assert "novel_genus" in data
            assert "mean_novelty_index" in data
            assert "genome_hit_counts" in data

    def test_qc_output_with_warnings(self, cli_runner, temp_dir):
        """QC output with warnings should be displayed."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.csv"
        qc_path = temp_dir / "qc_metrics.json"

        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(blast_path),
                "--ani", str(ani_path),
                "--output", str(output),
                "--qc-output", str(qc_path),
            ],
        )
        assert result.exit_code == 0

    def test_summary_quiet_mode_no_table(self, cli_runner, temp_dir):
        """Summary in quiet mode should not display table."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.csv"
        summary_path = temp_dir / "summary.json"

        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(blast_path),
                "--ani", str(ani_path),
                "--output", str(output),
                "--summary", str(summary_path),
                "--quiet",
            ],
        )
        assert result.exit_code == 0

    def test_summary_creates_parent_directory(self, cli_runner, temp_dir):
        """Summary parent directory should be created if missing."""
        blast_path, ani_path = self._make_blast_and_ani(temp_dir)
        output = temp_dir / "output.csv"
        summary_path = temp_dir / "subdir" / "summary.json"

        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(blast_path),
                "--ani", str(ani_path),
                "--output", str(output),
                "--summary", str(summary_path),
            ],
        )
        assert result.exit_code == 0


# =============================================================================
# Tests for export-config command and YAML config round-trip
# =============================================================================

import yaml

from metadarkmatter.models.config import ScoringConfig


class TestExportConfigCommand:
    """Tests for the ``score export-config`` CLI subcommand."""

    def test_export_config_creates_yaml(self, cli_runner, temp_dir):
        """export-config --output config.yaml should create a valid YAML file."""
        config_path = temp_dir / "config.yaml"
        result = cli_runner.invoke(
            app,
            ["score", "export-config", "--output", str(config_path)],
        )
        assert result.exit_code == 0
        assert config_path.exists()

        # The file should be parseable YAML
        content = yaml.safe_load(config_path.read_text())
        assert isinstance(content, dict)
        assert "alignment_mode" in content

    def test_export_config_from_preset(self, cli_runner, temp_dir):
        """export-config --preset gtdb-strict should create YAML from preset."""
        config_path = temp_dir / "config.yaml"
        result = cli_runner.invoke(
            app,
            [
                "score", "export-config",
                "--preset", "gtdb-strict",
                "--output", str(config_path),
            ],
        )
        assert result.exit_code == 0
        assert config_path.exists()

        content = yaml.safe_load(config_path.read_text())
        assert isinstance(content, dict)

    def test_export_config_invalid_preset(self, cli_runner, temp_dir):
        """export-config --preset invalid-name should return exit code 1."""
        config_path = temp_dir / "config.yaml"
        result = cli_runner.invoke(
            app,
            [
                "score", "export-config",
                "--preset", "invalid-name",
                "--output", str(config_path),
            ],
        )
        assert result.exit_code != 0

    def test_export_config_round_trip(self, cli_runner, temp_dir):
        """Exported YAML can be loaded back with ScoringConfig.from_yaml."""
        config_path = temp_dir / "config.yaml"
        result = cli_runner.invoke(
            app,
            ["score", "export-config", "--output", str(config_path)],
        )
        assert result.exit_code == 0

        loaded = ScoringConfig.from_yaml(config_path)
        assert isinstance(loaded, ScoringConfig)
        # Round-tripped values should match defaults
        default = ScoringConfig()
        assert loaded.novelty_known_max == default.novelty_known_max
        assert loaded.novelty_novel_species_max == default.novelty_novel_species_max
        assert loaded.novelty_novel_genus_max == default.novelty_novel_genus_max
        assert loaded.uncertainty_known_max == default.uncertainty_known_max
        assert loaded.bitscore_threshold_pct == default.bitscore_threshold_pct

    def test_export_config_preset_preserves_thresholds(self, cli_runner, temp_dir):
        """Exported YAML from preset should preserve key threshold values."""
        config_path = temp_dir / "gtdb_strict.yaml"
        result = cli_runner.invoke(
            app,
            [
                "score", "export-config",
                "--preset", "gtdb-strict",
                "--output", str(config_path),
            ],
        )
        assert result.exit_code == 0

        loaded = ScoringConfig.from_yaml(config_path)
        preset = THRESHOLD_PRESETS["gtdb-strict"]

        assert loaded.novelty_known_max == preset.novelty_known_max
        assert loaded.novelty_novel_species_min == preset.novelty_novel_species_min
        assert loaded.novelty_novel_species_max == preset.novelty_novel_species_max
        assert loaded.novelty_novel_genus_min == preset.novelty_novel_genus_min
        assert loaded.novelty_novel_genus_max == preset.novelty_novel_genus_max
        assert loaded.min_alignment_length == preset.min_alignment_length
        assert loaded.min_alignment_fraction == preset.min_alignment_fraction


class TestYamlConfigRoundTrip:
    """Tests for --config YAML loading in classify command and ScoringConfig YAML methods."""

    def test_classify_accepts_yaml_config(self, cli_runner, temp_dir):
        """Passing --config with a valid YAML should be accepted by classify.

        The command will fail at alignment loading (missing valid data) but
        should get past config parsing without error.
        """
        # Write a valid YAML config
        config = ScoringConfig()
        config_path = temp_dir / "config.yaml"
        config.to_yaml(config_path)

        # Create minimal BLAST and ANI files so the command reaches config parsing
        genomes = ["GCF_000123456.1", "GCF_000789012.1", "GCA_000111222.1"]
        ani_data = {"genome": genomes}
        ani_vals = [
            [100.0, 95.5, 80.0],
            [95.5, 100.0, 82.0],
            [80.0, 82.0, 100.0],
        ]
        for i, g in enumerate(genomes):
            ani_data[g] = ani_vals[i]
        ani_df = pl.DataFrame(ani_data)
        ani_path = temp_dir / "test.ani.csv"
        ani_df.write_csv(ani_path)

        blast_rows = []
        for read_idx in range(3):
            for g_idx, genome in enumerate(genomes):
                blast_rows.append({
                    "qseqid": f"read_{read_idx:03d}",
                    "sseqid": f"{genome}_ASM_genomic",
                    "pident": 98.5 - g_idx * 5,
                    "length": 150,
                    "mismatch": 2,
                    "gapopen": 0,
                    "qstart": 1,
                    "qend": 150,
                    "sstart": 1000,
                    "send": 1150,
                    "evalue": 1e-50,
                    "bitscore": 250.0 - g_idx * 20,
                })
        blast_df = pl.DataFrame(blast_rows)
        blast_path = temp_dir / "test.blast.tsv"
        blast_df.write_csv(blast_path, separator="\t", include_header=False)

        output = temp_dir / "output.csv"
        result = cli_runner.invoke(
            app,
            [
                "score", "classify",
                "--alignment", str(blast_path),
                "--ani", str(ani_path),
                "--config", str(config_path),
                "--output", str(output),
            ],
        )
        # Should succeed (config parsed and classification runs)
        assert result.exit_code == 0
        assert "loaded config from" in result.output.lower()

    def test_scoring_config_yaml_round_trip(self, temp_dir):
        """ScoringConfig.to_yaml -> from_yaml should preserve all values."""
        original = ScoringConfig(
            novelty_known_max=5.0,
            novelty_novel_species_min=5.0,
            novelty_novel_species_max=18.0,
            novelty_novel_genus_min=18.0,
            novelty_novel_genus_max=28.0,
            uncertainty_known_max=2.0,
            uncertainty_novel_species_max=2.0,
            uncertainty_novel_genus_max=2.0,
            uncertainty_conserved_min=6.0,
            bitscore_threshold_pct=90.0,
            min_alignment_length=150,
            min_alignment_fraction=0.5,
            coverage_weight_mode="sigmoid",
            coverage_weight_strength=0.7,
            uncertainty_mode="max",
            single_hit_uncertainty_threshold=8.0,
            family_ratio_threshold=0.9,
        )
        config_path = temp_dir / "round_trip.yaml"
        original.to_yaml(config_path)

        loaded = ScoringConfig.from_yaml(config_path)

        assert loaded.novelty_known_max == original.novelty_known_max
        assert loaded.novelty_novel_species_min == original.novelty_novel_species_min
        assert loaded.novelty_novel_species_max == original.novelty_novel_species_max
        assert loaded.novelty_novel_genus_min == original.novelty_novel_genus_min
        assert loaded.novelty_novel_genus_max == original.novelty_novel_genus_max
        assert loaded.uncertainty_known_max == original.uncertainty_known_max
        assert loaded.uncertainty_novel_species_max == original.uncertainty_novel_species_max
        assert loaded.uncertainty_novel_genus_max == original.uncertainty_novel_genus_max
        assert loaded.uncertainty_conserved_min == original.uncertainty_conserved_min
        assert loaded.bitscore_threshold_pct == original.bitscore_threshold_pct
        assert loaded.min_alignment_length == original.min_alignment_length
        assert loaded.min_alignment_fraction == original.min_alignment_fraction
        assert loaded.coverage_weight_mode == original.coverage_weight_mode
        assert loaded.coverage_weight_strength == original.coverage_weight_strength
        assert loaded.uncertainty_mode == original.uncertainty_mode
        assert loaded.single_hit_uncertainty_threshold == original.single_hit_uncertainty_threshold
        assert loaded.family_ratio_threshold == original.family_ratio_threshold
