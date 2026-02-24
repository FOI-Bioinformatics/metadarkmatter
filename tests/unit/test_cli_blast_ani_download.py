"""
Unit tests for BLAST, ANI, and Download CLI commands.

Tests argument validation, external tool mocking, file handling,
and error cases for the blast, ani, and download subcommands.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import polars as pl
import pytest

from metadarkmatter.cli.main import app
from metadarkmatter.external.base import ToolExecutionError, ToolResult
from metadarkmatter.external.ncbi_datasets import DownloadOutcome, DownloadReport


# =============================================================================
# Shared helpers
# =============================================================================


def _make_tool_result(
    command: tuple[str, ...] = ("mock_tool",),
    return_code: int = 0,
    stdout: str = "",
    stderr: str = "",
    elapsed_seconds: float = 1.5,
) -> ToolResult:
    """Create a ToolResult for mocking external tool calls."""
    return ToolResult(
        command=command,
        return_code=return_code,
        stdout=stdout,
        stderr=stderr,
        elapsed_seconds=elapsed_seconds,
    )


def _make_download_report(
    accessions: list[str] | None = None,
    elapsed_seconds: float = 1.5,
) -> DownloadReport:
    """Create a DownloadReport for mocking download_genomes calls."""
    outcomes = []
    for acc in accessions or []:
        outcomes.append(DownloadOutcome(accession=acc, success=True))
    return DownloadReport(outcomes=outcomes, elapsed_seconds=elapsed_seconds)


# =============================================================================
# BLAST CLI Tests
# =============================================================================


class TestBlastMakedbCommand:
    """Tests for 'blast makedb' CLI command."""

    def test_makedb_help_shows_options(self, cli_runner):
        """Help text should display key options."""
        result = cli_runner.invoke(app, ["blast", "makedb", "--help"])
        assert result.exit_code == 0
        assert "--genomes" in result.output
        assert "--output" in result.output
        assert "--title" in result.output

    def test_makedb_requires_genomes(self, cli_runner, temp_dir):
        """Command fails when --genomes is missing."""
        result = cli_runner.invoke(
            app,
            ["blast", "makedb", "--output", str(temp_dir / "db")],
        )
        assert result.exit_code != 0

    def test_makedb_requires_output(self, cli_runner, temp_dir):
        """Command fails when --output is missing."""
        result = cli_runner.invoke(
            app,
            ["blast", "makedb", "--genomes", str(temp_dir)],
        )
        assert result.exit_code != 0

    def test_makedb_nonexistent_genomes_path(self, cli_runner, temp_dir):
        """Command fails for a path that does not exist."""
        result = cli_runner.invoke(
            app,
            [
                "blast",
                "makedb",
                "--genomes",
                str(temp_dir / "no_such_dir"),
                "--output",
                str(temp_dir / "db"),
            ],
        )
        assert result.exit_code != 0

    @patch("metadarkmatter.cli.blast.MakeBlastDb.check_available", return_value=False)
    def test_makedb_tool_not_found(self, mock_check, cli_runner, temp_dir):
        """Command exits with error when makeblastdb is not installed."""
        genome_dir = temp_dir / "genomes"
        genome_dir.mkdir()
        (genome_dir / "genome1.fna").write_text(">contig1\nACGT\n")

        result = cli_runner.invoke(
            app,
            [
                "blast",
                "makedb",
                "--genomes",
                str(genome_dir),
                "--output",
                str(temp_dir / "db"),
            ],
        )
        assert result.exit_code != 0
        assert "not found" in result.output.lower() or "error" in result.output.lower()

    @patch("metadarkmatter.cli.blast.MakeBlastDb.check_available", return_value=True)
    @patch("metadarkmatter.cli.blast._concatenate_genomes", return_value=(3, 15))
    @patch("metadarkmatter.cli.blast.MakeBlastDb.run_or_raise")
    def test_makedb_happy_path_directory(
        self,
        mock_run,
        mock_concat,
        mock_check,
        cli_runner,
        temp_dir,
    ):
        """Happy path: makedb from a genome directory succeeds."""
        genome_dir = temp_dir / "genomes"
        genome_dir.mkdir()
        (genome_dir / "genome1.fna").write_text(">c1\nACGT\n")

        output_prefix = temp_dir / "blastdb" / "pangenome"
        mock_run.return_value = _make_tool_result(
            command=("makeblastdb", "-in", "pangenome.fasta", "-out", str(output_prefix)),
        )

        # Create fake db files so the summary section can stat them
        output_prefix.parent.mkdir(parents=True, exist_ok=True)
        for suffix in [".nhr", ".nin", ".nsq"]:
            (output_prefix.with_suffix(suffix)).write_bytes(b"\x00" * 100)

        result = cli_runner.invoke(
            app,
            [
                "blast",
                "makedb",
                "--genomes",
                str(genome_dir),
                "--output",
                str(output_prefix),
            ],
        )
        assert result.exit_code == 0
        assert "database ready" in result.output.lower() or "database created" in result.output.lower()
        mock_concat.assert_called_once()
        mock_run.assert_called_once()

    @patch("metadarkmatter.cli.blast.MakeBlastDb.check_available", return_value=True)
    @patch("metadarkmatter.cli.blast.MakeBlastDb.run_or_raise")
    def test_makedb_happy_path_single_fasta(
        self,
        mock_run,
        mock_check,
        cli_runner,
        temp_dir,
    ):
        """Happy path: makedb from a single FASTA file succeeds."""
        fasta = temp_dir / "pangenome.fasta"
        fasta.write_text(">contig1\nACGT\n>contig2\nTGCA\n")

        output_prefix = temp_dir / "blastdb" / "pangenome"
        mock_run.return_value = _make_tool_result()

        output_prefix.parent.mkdir(parents=True, exist_ok=True)
        for suffix in [".nhr", ".nin", ".nsq"]:
            (output_prefix.with_suffix(suffix)).write_bytes(b"\x00" * 200)

        result = cli_runner.invoke(
            app,
            [
                "blast",
                "makedb",
                "--genomes",
                str(fasta),
                "--output",
                str(output_prefix),
            ],
        )
        assert result.exit_code == 0
        mock_run.assert_called_once()

    @patch("metadarkmatter.cli.blast.MakeBlastDb.check_available", return_value=True)
    @patch("metadarkmatter.cli.blast._concatenate_genomes", return_value=(0, 0))
    def test_makedb_no_genomes_found(
        self,
        mock_concat,
        mock_check,
        cli_runner,
        temp_dir,
    ):
        """Command fails when the genome directory has no matching files."""
        genome_dir = temp_dir / "genomes"
        genome_dir.mkdir()

        result = cli_runner.invoke(
            app,
            [
                "blast",
                "makedb",
                "--genomes",
                str(genome_dir),
                "--output",
                str(temp_dir / "db"),
            ],
        )
        assert result.exit_code != 0
        assert "no genome" in result.output.lower() or "error" in result.output.lower()

    @patch("metadarkmatter.cli.blast.MakeBlastDb.check_available", return_value=True)
    def test_makedb_skip_if_exists(self, mock_check, cli_runner, temp_dir):
        """Command exits gracefully when database already exists and --skip-if-exists."""
        output_prefix = temp_dir / "db"
        for suffix in [".nhr", ".nin", ".nsq"]:
            output_prefix.with_suffix(suffix).write_bytes(b"\x00")

        fasta = temp_dir / "genome.fasta"
        fasta.write_text(">c1\nACGT\n")

        result = cli_runner.invoke(
            app,
            [
                "blast",
                "makedb",
                "--genomes",
                str(fasta),
                "--output",
                str(output_prefix),
                "--skip-if-exists",
            ],
        )
        assert result.exit_code == 0
        assert "skipping" in result.output.lower()

    @patch("metadarkmatter.cli.blast.MakeBlastDb.check_available", return_value=True)
    @patch("metadarkmatter.cli.blast.MakeBlastDb.run")
    def test_makedb_dry_run(self, mock_run, mock_check, cli_runner, temp_dir):
        """Dry run shows command without executing."""
        fasta = temp_dir / "genomes.fasta"
        fasta.write_text(">c1\nACGT\n")

        mock_run.return_value = _make_tool_result(
            command=("makeblastdb", "-in", str(fasta)),
        )

        result = cli_runner.invoke(
            app,
            [
                "blast",
                "makedb",
                "--genomes",
                str(fasta),
                "--output",
                str(temp_dir / "db"),
                "--dry-run",
            ],
        )
        assert result.exit_code == 0
        assert "dry run" in result.output.lower()

    @patch("metadarkmatter.cli.blast.MakeBlastDb.check_available", return_value=True)
    @patch(
        "metadarkmatter.cli.blast.MakeBlastDb.run_or_raise",
        side_effect=ToolExecutionError("makeblastdb", ["makeblastdb"], 1, "segfault"),
    )
    def test_makedb_tool_execution_error(
        self,
        mock_run,
        mock_check,
        cli_runner,
        temp_dir,
    ):
        """Command reports failure when makeblastdb returns non-zero."""
        fasta = temp_dir / "genomes.fasta"
        fasta.write_text(">c1\nACGT\n")

        result = cli_runner.invoke(
            app,
            [
                "blast",
                "makedb",
                "--genomes",
                str(fasta),
                "--output",
                str(temp_dir / "db"),
            ],
        )
        assert result.exit_code != 0
        assert "failed" in result.output.lower()


class TestBlastAlignCommand:
    """Tests for 'blast align' CLI command."""

    def test_align_help_shows_options(self, cli_runner):
        """Help text should display key options."""
        result = cli_runner.invoke(app, ["blast", "align", "--help"])
        assert result.exit_code == 0
        assert "--query" in result.output
        assert "--database" in result.output
        assert "--output" in result.output
        assert "--threads" in result.output

    def test_align_requires_query(self, cli_runner, temp_dir):
        """Command fails when --query is missing."""
        result = cli_runner.invoke(
            app,
            [
                "blast",
                "align",
                "--database",
                str(temp_dir / "db"),
                "--output",
                str(temp_dir / "out.tsv"),
            ],
        )
        assert result.exit_code != 0

    def test_align_nonexistent_query(self, cli_runner, temp_dir):
        """Command fails when query file does not exist."""
        result = cli_runner.invoke(
            app,
            [
                "blast",
                "align",
                "--query",
                str(temp_dir / "missing.fasta"),
                "--database",
                str(temp_dir / "db"),
                "--output",
                str(temp_dir / "out.tsv"),
            ],
        )
        assert result.exit_code != 0

    @patch("metadarkmatter.cli.blast.BlastN.check_available", return_value=False)
    def test_align_tool_not_found(self, mock_check, cli_runner, temp_dir):
        """Command exits with error when blastn is not installed."""
        query = temp_dir / "reads.fasta"
        query.write_text(">read1\nACGT\n")

        result = cli_runner.invoke(
            app,
            [
                "blast",
                "align",
                "--query",
                str(query),
                "--database",
                str(temp_dir / "db"),
                "--output",
                str(temp_dir / "out.tsv"),
            ],
        )
        assert result.exit_code != 0
        assert "not found" in result.output.lower()

    @patch("metadarkmatter.cli.blast.BlastN.check_available", return_value=True)
    def test_align_database_not_found(self, mock_check, cli_runner, temp_dir):
        """Command exits with error when database files are missing."""
        query = temp_dir / "reads.fasta"
        query.write_text(">read1\nACGT\n")

        result = cli_runner.invoke(
            app,
            [
                "blast",
                "align",
                "--query",
                str(query),
                "--database",
                str(temp_dir / "nonexistent_db"),
                "--output",
                str(temp_dir / "out.tsv"),
            ],
        )
        assert result.exit_code != 0
        assert "database not found" in result.output.lower() or "not found" in result.output.lower()

    @patch("metadarkmatter.cli.blast.BlastN.check_available", return_value=True)
    @patch("metadarkmatter.cli.blast.BlastN.run_or_raise")
    def test_align_happy_path(self, mock_run, mock_check, cli_runner, temp_dir):
        """Happy path: alignment runs and produces output."""
        query = temp_dir / "reads.fasta"
        query.write_text(">read1\nACGT\n")

        # Create fake database files
        db_prefix = temp_dir / "blastdb"
        for suffix in [".nhr", ".nin", ".nsq"]:
            db_prefix.with_suffix(suffix).write_bytes(b"\x00")

        output_tsv = temp_dir / "results.tsv"
        blast_line = "read1\tGCF_000001.1|c1\t98.5\t150\t2\t0\t1\t150\t100\t250\t1e-50\t250.0\t150\n"

        def _write_output(**kwargs):
            """Side effect that writes a mock BLAST output file."""
            out_path = kwargs.get("output", output_tsv)
            out_path.write_text(blast_line)
            return _make_tool_result(command=("blastn",))

        mock_run.side_effect = _write_output

        result = cli_runner.invoke(
            app,
            [
                "blast",
                "align",
                "--query",
                str(query),
                "--database",
                str(db_prefix),
                "--output",
                str(output_tsv),
                "--no-compress",
            ],
        )
        assert result.exit_code == 0
        assert "alignment complete" in result.output.lower() or "complete" in result.output.lower()
        mock_run.assert_called_once()

    @patch("metadarkmatter.cli.blast.BlastN.check_available", return_value=True)
    def test_align_skip_if_exists(self, mock_check, cli_runner, temp_dir):
        """Command exits gracefully when output file already exists and --skip-if-exists."""
        query = temp_dir / "reads.fasta"
        query.write_text(">read1\nACGT\n")

        db_prefix = temp_dir / "db"
        for suffix in [".nhr", ".nin", ".nsq"]:
            db_prefix.with_suffix(suffix).write_bytes(b"\x00")

        output = temp_dir / "out.tsv"
        output.write_text("existing data")

        result = cli_runner.invoke(
            app,
            [
                "blast",
                "align",
                "--query",
                str(query),
                "--database",
                str(db_prefix),
                "--output",
                str(output),
                "--skip-if-exists",
                "--no-compress",
            ],
        )
        assert result.exit_code == 0
        assert "skipping" in result.output.lower()

    @patch("metadarkmatter.cli.blast.BlastN.check_available", return_value=True)
    @patch("metadarkmatter.cli.blast.BlastN.run")
    def test_align_dry_run(self, mock_run, mock_check, cli_runner, temp_dir):
        """Dry run shows the command that would be executed."""
        query = temp_dir / "reads.fasta"
        query.write_text(">read1\nACGT\n")

        db_prefix = temp_dir / "db"
        for suffix in [".nhr", ".nin", ".nsq"]:
            db_prefix.with_suffix(suffix).write_bytes(b"\x00")

        mock_run.return_value = _make_tool_result(
            command=("blastn", "-query", str(query)),
        )

        result = cli_runner.invoke(
            app,
            [
                "blast",
                "align",
                "--query",
                str(query),
                "--database",
                str(db_prefix),
                "--output",
                str(temp_dir / "out.tsv"),
                "--dry-run",
            ],
        )
        assert result.exit_code == 0
        assert "dry run" in result.output.lower()

    @patch("metadarkmatter.cli.blast.BlastN.check_available", return_value=True)
    @patch(
        "metadarkmatter.cli.blast.BlastN.run_or_raise",
        side_effect=ToolExecutionError("blastn", ["blastn"], 1, "error details"),
    )
    def test_align_tool_execution_error(
        self,
        mock_run,
        mock_check,
        cli_runner,
        temp_dir,
    ):
        """Command reports failure when blastn returns non-zero."""
        query = temp_dir / "reads.fasta"
        query.write_text(">read1\nACGT\n")

        db_prefix = temp_dir / "db"
        for suffix in [".nhr", ".nin", ".nsq"]:
            db_prefix.with_suffix(suffix).write_bytes(b"\x00")

        result = cli_runner.invoke(
            app,
            [
                "blast",
                "align",
                "--query",
                str(query),
                "--database",
                str(db_prefix),
                "--output",
                str(temp_dir / "out.tsv"),
                "--no-compress",
            ],
        )
        assert result.exit_code != 0
        assert "failed" in result.output.lower()


class TestBlastHelperFunctions:
    """Tests for blast module helper functions."""

    @patch("metadarkmatter.cli.blast.MakeBlastDb.check_available", return_value=True)
    def test_check_makeblastdb_available_true(self, mock_check):
        """Returns (True, '') when makeblastdb is available."""
        from metadarkmatter.cli.blast import _check_makeblastdb_available

        ok, missing = _check_makeblastdb_available()
        assert ok is True
        assert missing == ""

    @patch("metadarkmatter.cli.blast.MakeBlastDb.check_available", return_value=False)
    def test_check_makeblastdb_available_false(self, mock_check):
        """Returns (False, ...) when makeblastdb is not available."""
        from metadarkmatter.cli.blast import _check_makeblastdb_available

        ok, missing = _check_makeblastdb_available()
        assert ok is False
        assert "makeblastdb" in missing.lower()

    @patch("metadarkmatter.cli.blast.BlastN.check_available", return_value=True)
    def test_check_blastn_available_true(self, mock_check):
        """Returns (True, '') when blastn is available."""
        from metadarkmatter.cli.blast import _check_blastn_available

        ok, missing = _check_blastn_available()
        assert ok is True
        assert missing == ""

    @patch("metadarkmatter.cli.blast.BlastN.check_available", return_value=False)
    def test_check_blastn_available_false(self, mock_check):
        """Returns (False, ...) when blastn is not available."""
        from metadarkmatter.cli.blast import _check_blastn_available

        ok, missing = _check_blastn_available()
        assert ok is False
        assert "blastn" in missing.lower()


# =============================================================================
# ANI CLI Tests
# =============================================================================


class TestANIComputeCommand:
    """Tests for 'ani compute' CLI command."""

    def test_compute_help_shows_options(self, cli_runner):
        """Help text should display key options."""
        result = cli_runner.invoke(app, ["ani", "compute", "--help"])
        assert result.exit_code == 0
        assert "--genomes" in result.output
        assert "--output" in result.output
        assert "--backend" in result.output
        assert "--threads" in result.output

    def test_compute_requires_genomes(self, cli_runner, temp_dir):
        """Command fails when --genomes is missing."""
        result = cli_runner.invoke(
            app,
            ["ani", "compute", "--output", str(temp_dir / "ani.csv")],
        )
        assert result.exit_code != 0

    def test_compute_requires_output(self, cli_runner, temp_dir):
        """Command fails when --output is missing."""
        result = cli_runner.invoke(
            app,
            ["ani", "compute", "--genomes", str(temp_dir)],
        )
        assert result.exit_code != 0

    def test_compute_nonexistent_genome_dir(self, cli_runner, temp_dir):
        """Command fails for a nonexistent genome directory."""
        result = cli_runner.invoke(
            app,
            [
                "ani",
                "compute",
                "--genomes",
                str(temp_dir / "no_such_dir"),
                "--output",
                str(temp_dir / "ani.csv"),
            ],
        )
        assert result.exit_code != 0

    @patch("metadarkmatter.external.skani.Skani.check_available", return_value=False)
    @patch("metadarkmatter.external.fastani.FastANI.check_available", return_value=False)
    def test_compute_no_backend_available(
        self, mock_fastani, mock_skani, cli_runner, temp_dir
    ):
        """Command exits with error when no ANI backend is found."""
        genome_dir = temp_dir / "genomes"
        genome_dir.mkdir()
        (genome_dir / "g1.fna").write_text(">c\nACGT\n")

        result = cli_runner.invoke(
            app,
            [
                "ani",
                "compute",
                "--genomes",
                str(genome_dir),
                "--output",
                str(temp_dir / "ani.csv"),
            ],
        )
        assert result.exit_code != 0
        assert "no ani backend" in result.output.lower() or "error" in result.output.lower()

    @patch("metadarkmatter.external.skani.Skani.check_available", return_value=True)
    def test_compute_no_genome_files(self, mock_check, cli_runner, temp_dir):
        """Command fails when no genome files match the pattern."""
        genome_dir = temp_dir / "genomes"
        genome_dir.mkdir()
        # Create a file that does not match *.fna
        (genome_dir / "readme.txt").write_text("hello")

        result = cli_runner.invoke(
            app,
            [
                "ani",
                "compute",
                "--genomes",
                str(genome_dir),
                "--output",
                str(temp_dir / "ani.csv"),
                "--backend",
                "skani",
                "--genome-pattern",
                "*.fna",
            ],
        )
        assert result.exit_code != 0
        assert "no genome" in result.output.lower() or "error" in result.output.lower()

    @patch("metadarkmatter.external.skani.Skani.check_available", return_value=True)
    @patch("metadarkmatter.external.skani.Skani.run_or_raise")
    @patch("metadarkmatter.cli.ani.parse_skani_output")
    @patch("metadarkmatter.cli.ani.ani_dict_to_csv", return_value=3)
    def test_compute_happy_path_skani(
        self,
        mock_csv,
        mock_parse,
        mock_run,
        mock_check,
        cli_runner,
        temp_dir,
    ):
        """Happy path: skani compute succeeds."""
        genome_dir = temp_dir / "genomes"
        genome_dir.mkdir()
        for i in range(3):
            (genome_dir / f"genome_{i}.fna").write_text(f">contig{i}\nACGTACGT\n")

        output_csv = temp_dir / "ani.csv"

        mock_run.return_value = _make_tool_result(command=("skani", "triangle"))
        mock_parse.return_value = {
            "genome_0": {"genome_0": 100.0, "genome_1": 95.0, "genome_2": 80.0},
            "genome_1": {"genome_0": 95.0, "genome_1": 100.0, "genome_2": 82.0},
            "genome_2": {"genome_0": 80.0, "genome_1": 82.0, "genome_2": 100.0},
        }

        # Create the output file so stat works in the summary
        def _write_csv(ani_dict, out_path, compress=False):
            pl.DataFrame(
                {
                    "genome": ["g0", "g1", "g2"],
                    "g0": [100.0, 95.0, 80.0],
                    "g1": [95.0, 100.0, 82.0],
                    "g2": [80.0, 82.0, 100.0],
                }
            ).write_csv(out_path)
            return 3

        mock_csv.side_effect = _write_csv

        result = cli_runner.invoke(
            app,
            [
                "ani",
                "compute",
                "--genomes",
                str(genome_dir),
                "--output",
                str(output_csv),
                "--backend",
                "skani",
            ],
        )
        assert result.exit_code == 0
        assert "computed successfully" in result.output.lower() or "complete" in result.output.lower()
        mock_run.assert_called_once()
        mock_parse.assert_called_once()

    @patch("metadarkmatter.external.fastani.FastANI.check_available", return_value=True)
    @patch("metadarkmatter.external.skani.Skani.check_available", return_value=False)
    @patch("metadarkmatter.external.fastani.FastANI.run_or_raise")
    @patch("metadarkmatter.cli.ani.parse_fastani_output")
    @patch("metadarkmatter.cli.ani.ani_dict_to_csv", return_value=2)
    @patch("metadarkmatter.external.fastani.create_genome_list_file")
    def test_compute_happy_path_fastani(
        self,
        mock_list,
        mock_csv,
        mock_parse,
        mock_run,
        mock_skani_check,
        mock_fastani_check,
        cli_runner,
        temp_dir,
    ):
        """Happy path: fastANI compute succeeds."""
        genome_dir = temp_dir / "genomes"
        genome_dir.mkdir()
        for i in range(2):
            (genome_dir / f"genome_{i}.fna").write_text(f">contig{i}\nACGT\n")

        output_csv = temp_dir / "ani.csv"

        mock_run.return_value = _make_tool_result(command=("fastANI",))
        mock_parse.return_value = {
            "genome_0": {"genome_0": 100.0, "genome_1": 95.0},
            "genome_1": {"genome_0": 95.0, "genome_1": 100.0},
        }

        def _write_csv(ani_dict, out_path, compress=False):
            pl.DataFrame(
                {
                    "genome": ["g0", "g1"],
                    "g0": [100.0, 95.0],
                    "g1": [95.0, 100.0],
                }
            ).write_csv(out_path)
            return 2

        mock_csv.side_effect = _write_csv

        result = cli_runner.invoke(
            app,
            [
                "ani",
                "compute",
                "--genomes",
                str(genome_dir),
                "--output",
                str(output_csv),
                "--backend",
                "fastani",
            ],
        )
        assert result.exit_code == 0
        mock_run.assert_called_once()

    @patch("metadarkmatter.external.skani.Skani.check_available", return_value=True)
    def test_compute_dry_run(self, mock_check, cli_runner, temp_dir):
        """Dry run shows command without executing."""
        genome_dir = temp_dir / "genomes"
        genome_dir.mkdir()
        (genome_dir / "g1.fna").write_text(">c\nACGT\n")

        result = cli_runner.invoke(
            app,
            [
                "ani",
                "compute",
                "--genomes",
                str(genome_dir),
                "--output",
                str(temp_dir / "ani.csv"),
                "--backend",
                "skani",
                "--dry-run",
            ],
        )
        assert result.exit_code == 0
        assert "dry run" in result.output.lower() or "would be" in result.output.lower()

    @patch("metadarkmatter.external.skani.Skani.check_available", return_value=True)
    @patch(
        "metadarkmatter.external.skani.Skani.run_or_raise",
        side_effect=ToolExecutionError("skani", ["skani"], 1, "error msg"),
    )
    def test_compute_tool_execution_error(
        self,
        mock_run,
        mock_check,
        cli_runner,
        temp_dir,
    ):
        """Command reports failure when the backend tool returns non-zero."""
        genome_dir = temp_dir / "genomes"
        genome_dir.mkdir()
        (genome_dir / "g1.fna").write_text(">c\nACGT\n")

        result = cli_runner.invoke(
            app,
            [
                "ani",
                "compute",
                "--genomes",
                str(genome_dir),
                "--output",
                str(temp_dir / "ani.csv"),
                "--backend",
                "skani",
            ],
        )
        assert result.exit_code != 0
        assert "failed" in result.output.lower()


class TestANIValidateCommand:
    """Tests for 'ani validate' CLI command."""

    def test_validate_help_shows_options(self, cli_runner):
        """Help text should display key options."""
        result = cli_runner.invoke(app, ["ani", "validate", "--help"])
        assert result.exit_code == 0
        assert "--ani" in result.output
        assert "--blast" in result.output

    def test_validate_requires_both_files(self, cli_runner, temp_dir):
        """Command fails when either --ani or --blast is missing."""
        result_no_ani = cli_runner.invoke(
            app,
            [
                "ani",
                "validate",
                "--blast",
                str(temp_dir / "blast.tsv"),
            ],
        )
        assert result_no_ani.exit_code != 0

        result_no_blast = cli_runner.invoke(
            app,
            [
                "ani",
                "validate",
                "--ani",
                str(temp_dir / "ani.csv"),
            ],
        )
        assert result_no_blast.exit_code != 0

    def test_validate_nonexistent_ani_file(self, cli_runner, temp_dir):
        """Command fails when ANI file does not exist."""
        blast_file = temp_dir / "blast.tsv"
        blast_file.write_text(
            "read1\tGCF_000001.1|c1\t98.0\t150\t3\t0\t1\t150\t100\t250\t1e-50\t250\n"
        )

        result = cli_runner.invoke(
            app,
            [
                "ani",
                "validate",
                "--ani",
                str(temp_dir / "missing.csv"),
                "--blast",
                str(blast_file),
            ],
        )
        assert result.exit_code != 0

    def test_validate_nonexistent_blast_file(self, cli_runner, temp_dir, temp_ani_file):
        """Command fails when BLAST file does not exist."""
        result = cli_runner.invoke(
            app,
            [
                "ani",
                "validate",
                "--ani",
                str(temp_ani_file),
                "--blast",
                str(temp_dir / "missing_blast.tsv"),
            ],
        )
        assert result.exit_code != 0


class TestANIHelperFunctions:
    """Tests for ani module helper functions."""

    def test_find_genome_files_primary_pattern(self, temp_dir):
        """Finds genome files matching the primary pattern."""
        from metadarkmatter.cli.ani import _find_genome_files

        for i in range(3):
            (temp_dir / f"g{i}.fna").write_text(f">c{i}\nACGT\n")

        files = _find_genome_files(temp_dir, "*.fna")
        assert len(files) == 3

    def test_find_genome_files_fallback_pattern(self, temp_dir):
        """Falls back to alternative patterns when primary matches nothing."""
        from metadarkmatter.cli.ani import _find_genome_files

        (temp_dir / "genome.fasta").write_text(">c\nACGT\n")

        files = _find_genome_files(temp_dir, "*.xyz")
        assert len(files) == 1
        assert files[0].name == "genome.fasta"

    def test_find_genome_files_empty_directory(self, temp_dir):
        """Returns empty list for a directory with no genome files."""
        from metadarkmatter.cli.ani import _find_genome_files

        files = _find_genome_files(temp_dir, "*.fna")
        assert files == []

    def test_detect_backend_skani_preferred(self):
        """Auto-detect prefers skani over fastANI."""
        from metadarkmatter.cli.ani import ANIBackend, _detect_backend

        with (
            patch("metadarkmatter.external.skani.Skani.check_available", return_value=True),
            patch("metadarkmatter.external.fastani.FastANI.check_available", return_value=True),
        ):
            backend = _detect_backend()
            assert backend == ANIBackend.SKANI

    def test_detect_backend_fastani_fallback(self):
        """Falls back to fastANI when skani is not available."""
        from metadarkmatter.cli.ani import ANIBackend, _detect_backend

        with (
            patch("metadarkmatter.external.skani.Skani.check_available", return_value=False),
            patch("metadarkmatter.external.fastani.FastANI.check_available", return_value=True),
        ):
            backend = _detect_backend()
            assert backend == ANIBackend.FASTANI

    def test_detect_backend_none_exits(self):
        """Raises typer.Exit when no backend is available."""
        from click.exceptions import Exit

        from metadarkmatter.cli.ani import _detect_backend

        with (
            patch("metadarkmatter.external.skani.Skani.check_available", return_value=False),
            patch("metadarkmatter.external.fastani.FastANI.check_available", return_value=False),
            pytest.raises(Exit),
        ):
            _detect_backend()


class TestANISensitivityPresets:
    """Tests that sensitivity presets affect parameters."""

    @patch("metadarkmatter.external.skani.Skani.check_available", return_value=True)
    def test_sensitive_preset_dry_run(self, mock_check, cli_runner, temp_dir):
        """Sensitive mode shows lower min_af and screen values in dry run."""
        genome_dir = temp_dir / "genomes"
        genome_dir.mkdir()
        (genome_dir / "g1.fna").write_text(">c\nACGT\n")

        result = cli_runner.invoke(
            app,
            [
                "ani",
                "compute",
                "--genomes",
                str(genome_dir),
                "--output",
                str(temp_dir / "ani.csv"),
                "--backend",
                "skani",
                "--sensitivity",
                "sensitive",
                "--dry-run",
            ],
        )
        assert result.exit_code == 0
        # Sensitive defaults: min_af=5, screen=50
        assert "5" in result.output or "50" in result.output

    @patch("metadarkmatter.external.skani.Skani.check_available", return_value=True)
    def test_default_preset_dry_run(self, mock_check, cli_runner, temp_dir):
        """Default mode shows standard min_af and screen values in dry run."""
        genome_dir = temp_dir / "genomes"
        genome_dir.mkdir()
        (genome_dir / "g1.fna").write_text(">c\nACGT\n")

        result = cli_runner.invoke(
            app,
            [
                "ani",
                "compute",
                "--genomes",
                str(genome_dir),
                "--output",
                str(temp_dir / "ani.csv"),
                "--backend",
                "skani",
                "--sensitivity",
                "default",
                "--dry-run",
            ],
        )
        assert result.exit_code == 0
        # Default: min_af=15, screen=80
        assert "15" in result.output or "80" in result.output


# =============================================================================
# Download CLI Tests
# =============================================================================


class TestDownloadGenomesListCommand:
    """Tests for 'download genomes list' CLI command."""

    def test_list_help_shows_options(self, cli_runner):
        """Help text should display key options."""
        result = cli_runner.invoke(app, ["download", "genomes", "list", "--help"])
        assert result.exit_code == 0
        assert "--output" in result.output

    def test_list_requires_taxon(self, cli_runner, temp_dir):
        """Command fails when taxon argument is missing."""
        result = cli_runner.invoke(
            app,
            [
                "download",
                "genomes",
                "list",
                "--output",
                str(temp_dir / "genomes.tsv"),
            ],
        )
        assert result.exit_code != 0

    def test_list_requires_output(self, cli_runner):
        """Command fails when --output is missing."""
        result = cli_runner.invoke(
            app,
            ["download", "genomes", "list", "f__Francisellaceae"],
        )
        assert result.exit_code != 0

    def test_list_dry_run(self, cli_runner, temp_dir):
        """Dry run shows query details without making requests."""
        result = cli_runner.invoke(
            app,
            [
                "download",
                "genomes",
                "list",
                "f__Francisellaceae",
                "--output",
                str(temp_dir / "genomes.tsv"),
                "--dry-run",
            ],
        )
        assert result.exit_code == 0
        assert "dry run" in result.output.lower()

    @patch("metadarkmatter.cli.download.GTDBClient")
    def test_list_invalid_taxon_format(self, mock_client_cls, cli_runner, temp_dir):
        """Command reports error for an invalid GTDB taxon format."""
        from metadarkmatter.clients.gtdb import InvalidTaxonFormatError

        mock_client = MagicMock()
        mock_client.query_genomes.side_effect = InvalidTaxonFormatError("BadTaxon")
        mock_client_cls.return_value = mock_client

        result = cli_runner.invoke(
            app,
            [
                "download",
                "genomes",
                "list",
                "BadTaxon",
                "--output",
                str(temp_dir / "genomes.tsv"),
            ],
        )
        assert result.exit_code != 0

    @patch("metadarkmatter.cli.download.GTDBClient")
    def test_list_api_error(self, mock_client_cls, cli_runner, temp_dir):
        """Command reports error when GTDB API returns an error."""
        from metadarkmatter.clients.gtdb import GTDBAPIError

        mock_client = MagicMock()
        mock_client.query_genomes.side_effect = GTDBAPIError("Server error", status_code=500)
        mock_client_cls.return_value = mock_client

        result = cli_runner.invoke(
            app,
            [
                "download",
                "genomes",
                "list",
                "f__Francisellaceae",
                "--output",
                str(temp_dir / "genomes.tsv"),
            ],
        )
        assert result.exit_code != 0

    @patch("metadarkmatter.cli.download.GTDBClient")
    def test_list_no_genomes_found(self, mock_client_cls, cli_runner, temp_dir):
        """Command reports error when GTDB returns zero genomes."""
        from metadarkmatter.clients.gtdb import GTDBQueryResult

        mock_client = MagicMock()
        mock_client.query_genomes.return_value = GTDBQueryResult(
            taxon="f__EmptyFamily",
            genomes=(),
            total_count=0,
            genus_counts={},
            species_counts={},
        )
        mock_client_cls.return_value = mock_client

        result = cli_runner.invoke(
            app,
            [
                "download",
                "genomes",
                "list",
                "f__EmptyFamily",
                "--output",
                str(temp_dir / "genomes.tsv"),
            ],
        )
        assert result.exit_code != 0
        assert "no genomes" in result.output.lower()

    @patch("metadarkmatter.cli.download.GTDBClient")
    def test_list_happy_path(self, mock_client_cls, cli_runner, temp_dir):
        """Happy path: list genomes queries GTDB and writes TSV."""
        from metadarkmatter.clients.gtdb import GTDBGenome, GTDBQueryResult

        genomes = (
            GTDBGenome(
                accession="GCF_000001.1",
                gtdb_taxonomy="d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Thiotrichales;f__Francisellaceae;g__Francisella;s__Francisella tularensis",
                species="Francisella tularensis",
                genome_size=1900000,
            ),
            GTDBGenome(
                accession="GCF_000002.1",
                gtdb_taxonomy="d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Thiotrichales;f__Francisellaceae;g__Francisella;s__Francisella novicida",
                species="Francisella novicida",
                genome_size=1800000,
            ),
        )
        query_result = GTDBQueryResult(
            taxon="f__Francisellaceae",
            genomes=genomes,
            total_count=2,
            genus_counts={"Francisella": 2},
            species_counts={"Francisella tularensis": 1, "Francisella novicida": 1},
        )

        mock_client = MagicMock()
        mock_client.query_genomes.return_value = query_result
        mock_client_cls.return_value = mock_client

        output = temp_dir / "genomes.tsv"

        result = cli_runner.invoke(
            app,
            [
                "download",
                "genomes",
                "list",
                "f__Francisellaceae",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code == 0
        assert "query complete" in result.output.lower() or "complete" in result.output.lower()
        assert output.exists()
        # Verify accession list content
        df = pl.read_csv(output, separator="\t")
        assert len(df) == 2
        assert "GCF_000001.1" in df["accession"].to_list()
        # Verify metadata file was also created
        metadata_path = temp_dir / "genome_metadata.tsv"
        assert metadata_path.exists()

    @patch("metadarkmatter.cli.download.GTDBClient")
    def test_list_verbose_shows_genus_table(self, mock_client_cls, cli_runner, temp_dir):
        """Verbose mode shows genus breakdown table."""
        from metadarkmatter.clients.gtdb import GTDBGenome, GTDBQueryResult

        genomes = (
            GTDBGenome(
                accession="GCF_000001.1",
                gtdb_taxonomy="d__Bacteria;f__Francisellaceae;g__Francisella;s__Francisella tularensis",
                species="Francisella tularensis",
            ),
        )
        mock_client = MagicMock()
        mock_client.query_genomes.return_value = GTDBQueryResult(
            taxon="f__Francisellaceae",
            genomes=genomes,
            total_count=1,
            genus_counts={"Francisella": 1},
            species_counts={"Francisella tularensis": 1},
        )
        mock_client_cls.return_value = mock_client

        result = cli_runner.invoke(
            app,
            [
                "download",
                "genomes",
                "list",
                "f__Francisellaceae",
                "--output",
                str(temp_dir / "genomes.tsv"),
                "--verbose",
            ],
        )
        assert result.exit_code == 0
        assert "francisella" in result.output.lower()


class TestDownloadGenomesFetchCommand:
    """Tests for 'download genomes fetch' CLI command."""

    def test_fetch_help_shows_options(self, cli_runner):
        """Help text should display key options."""
        result = cli_runner.invoke(app, ["download", "genomes", "fetch", "--help"])
        assert result.exit_code == 0
        assert "--accessions" in result.output
        assert "--output-dir" in result.output

    def test_fetch_requires_accessions(self, cli_runner, temp_dir):
        """Command fails when --accessions is missing."""
        result = cli_runner.invoke(
            app,
            [
                "download",
                "genomes",
                "fetch",
                "--output-dir",
                str(temp_dir / "genomes"),
            ],
        )
        assert result.exit_code != 0

    def test_fetch_requires_output_dir(self, cli_runner, temp_dir):
        """Command fails when --output-dir is missing."""
        acc_file = temp_dir / "acc.tsv"
        acc_file.write_text("accession\tgtdb_taxonomy\tspecies\nGCF_000001.1\ttax\tsp\n")

        result = cli_runner.invoke(
            app,
            [
                "download",
                "genomes",
                "fetch",
                "--accessions",
                str(acc_file),
            ],
        )
        assert result.exit_code != 0

    def test_fetch_nonexistent_accessions_file(self, cli_runner, temp_dir):
        """Command fails when accessions file does not exist."""
        result = cli_runner.invoke(
            app,
            [
                "download",
                "genomes",
                "fetch",
                "--accessions",
                str(temp_dir / "missing.tsv"),
                "--output-dir",
                str(temp_dir / "out"),
            ],
        )
        assert result.exit_code != 0

    @patch("metadarkmatter.cli.download.NCBIDatasets")
    def test_fetch_ncbi_tool_not_found(self, mock_ncbi_cls, cli_runner, temp_dir):
        """Command exits with error when NCBI datasets CLI is not installed."""
        mock_ncbi = MagicMock()
        mock_ncbi.check_available.return_value = False
        mock_ncbi_cls.return_value = mock_ncbi

        acc_file = temp_dir / "acc.tsv"
        acc_file.write_text("accession\tgtdb_taxonomy\tspecies\nGCF_000001.1\ttax\tsp\n")

        result = cli_runner.invoke(
            app,
            [
                "download",
                "genomes",
                "fetch",
                "--accessions",
                str(acc_file),
                "--output-dir",
                str(temp_dir / "out"),
            ],
        )
        assert result.exit_code != 0
        assert "not found" in result.output.lower() or "error" in result.output.lower()

    @patch("metadarkmatter.cli.download.NCBIDatasets")
    def test_fetch_dry_run(self, mock_ncbi_cls, cli_runner, temp_dir):
        """Dry run shows what would be downloaded without making requests."""
        mock_ncbi = MagicMock()
        mock_ncbi.check_available.return_value = True
        mock_ncbi_cls.return_value = mock_ncbi

        acc_file = temp_dir / "acc.tsv"
        pl.DataFrame(
            {
                "accession": ["GCF_000001.1", "GCF_000002.1"],
                "gtdb_taxonomy": ["tax1", "tax2"],
                "species": ["sp1", "sp2"],
            }
        ).write_csv(acc_file, separator="\t")

        result = cli_runner.invoke(
            app,
            [
                "download",
                "genomes",
                "fetch",
                "--accessions",
                str(acc_file),
                "--output-dir",
                str(temp_dir / "out"),
                "--dry-run",
            ],
        )
        assert result.exit_code == 0
        assert "dry run" in result.output.lower()

    @patch("metadarkmatter.cli.download.NCBIDatasets")
    def test_fetch_all_already_exist(self, mock_ncbi_cls, cli_runner, temp_dir):
        """Command exits gracefully when all genomes are already downloaded."""
        mock_ncbi = MagicMock()
        mock_ncbi.check_available.return_value = True
        mock_ncbi_cls.return_value = mock_ncbi

        acc_file = temp_dir / "acc.tsv"
        pl.DataFrame(
            {
                "accession": ["GCF_000001.1"],
                "gtdb_taxonomy": ["tax"],
                "species": ["sp"],
            }
        ).write_csv(acc_file, separator="\t")

        output_dir = temp_dir / "genomes"
        output_dir.mkdir()
        (output_dir / "GCF_000001.1.fna").write_text(">c\nACGT\n")

        result = cli_runner.invoke(
            app,
            [
                "download",
                "genomes",
                "fetch",
                "--accessions",
                str(acc_file),
                "--output-dir",
                str(output_dir),
            ],
        )
        assert result.exit_code == 0
        assert "already exist" in result.output.lower()

    @patch("metadarkmatter.cli.download.NCBIDatasets")
    def test_fetch_empty_accession_list(self, mock_ncbi_cls, cli_runner, temp_dir):
        """Command exits when accession list has zero entries."""
        mock_ncbi = MagicMock()
        mock_ncbi.check_available.return_value = True
        mock_ncbi_cls.return_value = mock_ncbi

        acc_file = temp_dir / "acc.tsv"
        pl.DataFrame(
            {
                "accession": [],
                "gtdb_taxonomy": [],
                "species": [],
            }
        ).write_csv(acc_file, separator="\t")

        result = cli_runner.invoke(
            app,
            [
                "download",
                "genomes",
                "fetch",
                "--accessions",
                str(acc_file),
                "--output-dir",
                str(temp_dir / "out"),
            ],
        )
        assert result.exit_code != 0
        assert "no accessions" in result.output.lower() or "0" in result.output

    @patch("metadarkmatter.cli.download.NCBIDatasets")
    def test_fetch_happy_path(self, mock_ncbi_cls, cli_runner, temp_dir):
        """Happy path: fetching genomes invokes download and reports summary."""
        mock_ncbi = MagicMock()
        mock_ncbi.check_available.return_value = True
        mock_ncbi_cls.return_value = mock_ncbi

        acc_file = temp_dir / "acc.tsv"
        pl.DataFrame(
            {
                "accession": ["GCF_000001.1", "GCF_000002.1"],
                "gtdb_taxonomy": ["tax1", "tax2"],
                "species": ["sp1", "sp2"],
            }
        ).write_csv(acc_file, separator="\t")

        output_dir = temp_dir / "genomes"

        # Simulate the download creating genome files and returning DownloadReport
        def _mock_download(**kwargs):
            od = kwargs.get("output_dir", output_dir)
            od.mkdir(parents=True, exist_ok=True)
            (od / "GCF_000001.1.fna").write_text(">c1\nACGTACGT\n")
            (od / "GCF_000002.1.fna").write_text(">c2\nTGCATGCA\n")
            return _make_download_report(
                accessions=["GCF_000001.1", "GCF_000002.1"],
                elapsed_seconds=3.0,
            )

        mock_ncbi.download_genomes.side_effect = _mock_download

        result = cli_runner.invoke(
            app,
            [
                "download",
                "genomes",
                "fetch",
                "--accessions",
                str(acc_file),
                "--output-dir",
                str(output_dir),
                "--redownload",
            ],
        )
        assert result.exit_code == 0
        assert "download complete" in result.output.lower() or "complete" in result.output.lower()
        mock_ncbi.download_genomes.assert_called_once()

    @patch("metadarkmatter.cli.download.NCBIDatasets")
    def test_fetch_download_failure(self, mock_ncbi_cls, cli_runner, temp_dir):
        """Command reports failure when download raises an exception."""
        mock_ncbi = MagicMock()
        mock_ncbi.check_available.return_value = True
        mock_ncbi.download_genomes.side_effect = RuntimeError("network timeout")
        mock_ncbi_cls.return_value = mock_ncbi

        acc_file = temp_dir / "acc.tsv"
        pl.DataFrame(
            {
                "accession": ["GCF_000001.1"],
                "gtdb_taxonomy": ["tax"],
                "species": ["sp"],
            }
        ).write_csv(acc_file, separator="\t")

        result = cli_runner.invoke(
            app,
            [
                "download",
                "genomes",
                "fetch",
                "--accessions",
                str(acc_file),
                "--output-dir",
                str(temp_dir / "out"),
                "--redownload",
            ],
        )
        assert result.exit_code != 0
        assert "failed" in result.output.lower() or "error" in result.output.lower()

    @patch("metadarkmatter.cli.download.NCBIDatasets")
    def test_fetch_skip_already_downloaded(self, mock_ncbi_cls, cli_runner, temp_dir):
        """With skip_if_exists (default), only missing genomes are downloaded."""
        mock_ncbi = MagicMock()
        mock_ncbi.check_available.return_value = True
        mock_ncbi_cls.return_value = mock_ncbi

        acc_file = temp_dir / "acc.tsv"
        pl.DataFrame(
            {
                "accession": ["GCF_000001.1", "GCF_000002.1"],
                "gtdb_taxonomy": ["tax1", "tax2"],
                "species": ["sp1", "sp2"],
            }
        ).write_csv(acc_file, separator="\t")

        output_dir = temp_dir / "genomes"
        output_dir.mkdir()
        (output_dir / "GCF_000001.1.fna").write_text(">c\nACGT\n")

        def _mock_download(**kwargs):
            od = kwargs.get("output_dir", output_dir)
            (od / "GCF_000002.1.fna").write_text(">c2\nTGCA\n")
            return _make_download_report(
                accessions=["GCF_000002.1"],
            )

        mock_ncbi.download_genomes.side_effect = _mock_download

        result = cli_runner.invoke(
            app,
            [
                "download",
                "genomes",
                "fetch",
                "--accessions",
                str(acc_file),
                "--output-dir",
                str(output_dir),
            ],
        )
        assert result.exit_code == 0
        # Verify skipping message appears
        assert "skipping" in result.output.lower() or "1" in result.output
