"""
Unit tests for CLI modules with low coverage.

Tests for kraken2, map, proteins, blastx, and aai CLI commands.
Uses typer.testing.CliRunner with mocked external tool calls.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.main import app
from metadarkmatter.external.base import ToolExecutionError, ToolResult


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _make_tool_result(
    command: tuple[str, ...] = ("mock",),
    return_code: int = 0,
    stdout: str = "",
    stderr: str = "",
    elapsed: float = 1.5,
) -> ToolResult:
    """Create a ToolResult for mocking external tool calls."""
    return ToolResult(
        command=command,
        return_code=return_code,
        stdout=stdout,
        stderr=stderr,
        elapsed_seconds=elapsed,
    )


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


# ===========================================================================
# Kraken2 CLI Tests
# ===========================================================================


class TestKraken2Classify:
    """Tests for kraken2 classify command."""

    def test_classify_tool_not_available(self, runner: CliRunner) -> None:
        """Should exit with error when kraken2 is not installed."""
        with patch(
            "metadarkmatter.cli.kraken2.Kraken2.check_available",
            return_value=False,
        ):
            result = runner.invoke(
                app,
                [
                    "kraken2",
                    "classify",
                    "--reads-1", "/fake/reads_R1.fastq.gz",
                    "--kraken-db", "/fake/krakendb",
                    "--output", "/fake/output",
                ],
            )
        assert result.exit_code != 0

    def test_classify_missing_reads_file(self, runner: CliRunner) -> None:
        """Should fail for non-existent reads file."""
        result = runner.invoke(
            app,
            [
                "kraken2",
                "classify",
                "--reads-1", "/nonexistent/reads_R1.fastq.gz",
                "--kraken-db", "/fake/krakendb",
                "--output", "/fake/output",
            ],
        )
        assert result.exit_code != 0

    def test_classify_happy_path(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should run classify successfully with mocked Kraken2."""
        reads = tmp_path / "sample_R1.fastq.gz"
        reads.write_text("@read1\nACGT\n+\nIIII\n")
        kraken_db = tmp_path / "krakendb"
        kraken_db.mkdir()
        output_dir = tmp_path / "kraken_out"

        mock_result = _make_tool_result(
            command=("kraken2", "--db", str(kraken_db)),
        )

        with (
            patch(
                "metadarkmatter.cli.kraken2.Kraken2.check_available",
                return_value=True,
            ),
            patch(
                "metadarkmatter.cli.kraken2.Kraken2.run_or_raise",
                return_value=mock_result,
            ),
            patch(
                "metadarkmatter.cli.kraken2.KrakenReport",
            ),
        ):
            result = runner.invoke(
                app,
                [
                    "kraken2",
                    "classify",
                    "--reads-1", str(reads),
                    "--kraken-db", str(kraken_db),
                    "--output", str(output_dir),
                    "--quiet",
                ],
            )
        assert result.exit_code == 0

    def test_classify_dry_run(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Dry run should show command without executing."""
        reads = tmp_path / "sample_R1.fastq.gz"
        reads.write_text("@read1\nACGT\n+\nIIII\n")
        kraken_db = tmp_path / "krakendb"
        kraken_db.mkdir()
        output_dir = tmp_path / "kraken_out"

        mock_result = _make_tool_result()

        with (
            patch(
                "metadarkmatter.cli.kraken2.Kraken2.check_available",
                return_value=True,
            ),
            patch(
                "metadarkmatter.cli.kraken2.Kraken2.run",
                return_value=mock_result,
            ),
        ):
            result = runner.invoke(
                app,
                [
                    "kraken2",
                    "classify",
                    "--reads-1", str(reads),
                    "--kraken-db", str(kraken_db),
                    "--output", str(output_dir),
                    "--dry-run",
                ],
            )
        assert result.exit_code == 0
        assert "DRY RUN" in result.output

    def test_classify_skip_if_exists(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should skip when output files already exist."""
        reads = tmp_path / "sample_R1.fastq.gz"
        reads.write_text("@read1\nACGT\n+\nIIII\n")
        kraken_db = tmp_path / "krakendb"
        kraken_db.mkdir()
        output_dir = tmp_path / "kraken_out"
        output_dir.mkdir()
        # Create existing outputs
        (output_dir / "sample.kraken").write_text("data")
        (output_dir / "sample.kreport").write_text("data")

        with patch(
            "metadarkmatter.cli.kraken2.Kraken2.check_available",
            return_value=True,
        ):
            result = runner.invoke(
                app,
                [
                    "kraken2",
                    "classify",
                    "--reads-1", str(reads),
                    "--kraken-db", str(kraken_db),
                    "--output", str(output_dir),
                    "--sample-name", "sample",
                    "--skip-if-exists",
                ],
            )
        assert result.exit_code == 0

    def test_classify_tool_execution_error(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should handle ToolExecutionError from kraken2."""
        reads = tmp_path / "sample_R1.fastq.gz"
        reads.write_text("@read1\nACGT\n+\nIIII\n")
        kraken_db = tmp_path / "krakendb"
        kraken_db.mkdir()
        output_dir = tmp_path / "kraken_out"

        with (
            patch(
                "metadarkmatter.cli.kraken2.Kraken2.check_available",
                return_value=True,
            ),
            patch(
                "metadarkmatter.cli.kraken2.Kraken2.run_or_raise",
                side_effect=ToolExecutionError(
                    "kraken2", ["kraken2"], 1, "kraken2 failed"
                ),
            ),
        ):
            result = runner.invoke(
                app,
                [
                    "kraken2",
                    "classify",
                    "--reads-1", str(reads),
                    "--kraken-db", str(kraken_db),
                    "--output", str(output_dir),
                    "--quiet",
                ],
            )
        assert result.exit_code != 0


class TestKraken2Extract:
    """Tests for kraken2 extract command."""

    def test_extract_tool_not_available(self, runner: CliRunner) -> None:
        """Should exit with error when krakentools is not installed."""
        with patch(
            "metadarkmatter.cli.kraken2.ExtractKrakenReads.check_available",
            return_value=False,
        ):
            result = runner.invoke(
                app,
                [
                    "kraken2",
                    "extract",
                    "--kraken-output", "/fake/sample.kraken",
                    "--reads-1", "/fake/reads.fastq.gz",
                    "--taxid", "262",
                    "--output", "/fake/out",
                ],
            )
        assert result.exit_code != 0

    def test_extract_missing_kreport(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should fail when kreport file is missing."""
        kraken_out = tmp_path / "sample.kraken"
        kraken_out.write_text("C\tread1\t262\t150\t262:150")
        reads = tmp_path / "reads_R1.fastq.gz"
        reads.write_text("@read1\nACGT\n+\nIIII\n")
        output_dir = tmp_path / "extracted"

        with patch(
            "metadarkmatter.cli.kraken2.ExtractKrakenReads.check_available",
            return_value=True,
        ):
            result = runner.invoke(
                app,
                [
                    "kraken2",
                    "extract",
                    "--kraken-output", str(kraken_out),
                    "--reads-1", str(reads),
                    "--taxid", "262",
                    "--output", str(output_dir),
                ],
            )
        assert result.exit_code != 0

    def test_extract_happy_path(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should successfully extract reads with mocked tools."""
        kraken_out = tmp_path / "sample.kraken"
        kraken_out.write_text("C\tread1\t262\t150\t262:150")
        kreport = tmp_path / "sample.kreport"
        kreport.write_text("100.00\t1000\t0\tR\t1\troot\n")
        reads = tmp_path / "reads_R1.fastq.gz"
        reads.write_text("@read1\nACGT\n+\nIIII\n")
        output_dir = tmp_path / "extracted"

        mock_result = _make_tool_result()

        def fake_run_or_raise(**kwargs):
            # Create the temp output file that extraction would create
            output_1 = kwargs.get("output_1")
            if output_1 and isinstance(output_1, Path):
                output_1.parent.mkdir(parents=True, exist_ok=True)
                output_1.write_text("@read1\nACGT\n+\nIIII\n")
            return mock_result

        with (
            patch(
                "metadarkmatter.cli.kraken2.ExtractKrakenReads.check_available",
                return_value=True,
            ),
            patch(
                "metadarkmatter.cli.kraken2.ExtractKrakenReads.run_or_raise",
                side_effect=fake_run_or_raise,
            ),
            patch(
                "metadarkmatter.cli.kraken2.KrakenReport",
            ),
        ):
            result = runner.invoke(
                app,
                [
                    "kraken2",
                    "extract",
                    "--kraken-output", str(kraken_out),
                    "--reads-1", str(reads),
                    "--taxid", "262",
                    "--output", str(output_dir),
                    "--quiet",
                ],
            )
        assert result.exit_code == 0

    def test_extract_dry_run(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Dry run should show command without executing."""
        kraken_out = tmp_path / "sample.kraken"
        kraken_out.write_text("C\tread1\t262\t150\t262:150")
        kreport = tmp_path / "sample.kreport"
        kreport.write_text("100.00\t1000\t0\tR\t1\troot\n")
        reads = tmp_path / "reads_R1.fastq.gz"
        reads.write_text("@read1\nACGT\n+\nIIII\n")
        output_dir = tmp_path / "extracted"

        mock_result = _make_tool_result()

        with (
            patch(
                "metadarkmatter.cli.kraken2.ExtractKrakenReads.check_available",
                return_value=True,
            ),
            patch(
                "metadarkmatter.cli.kraken2.ExtractKrakenReads.run",
                return_value=mock_result,
            ),
        ):
            result = runner.invoke(
                app,
                [
                    "kraken2",
                    "extract",
                    "--kraken-output", str(kraken_out),
                    "--reads-1", str(reads),
                    "--taxid", "262",
                    "--output", str(output_dir),
                    "--dry-run",
                ],
            )
        assert result.exit_code == 0
        assert "DRY RUN" in result.output

    def test_extract_skip_if_exists(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should skip when extracted reads already exist."""
        kraken_out = tmp_path / "sample.kraken"
        kraken_out.write_text("C\tread1\t262\t150\t262:150")
        kreport = tmp_path / "sample.kreport"
        kreport.write_text("100.00\t1000\t0\tR\t1\troot\n")
        reads = tmp_path / "reads_R1.fastq.gz"
        reads.write_text("@read1\nACGT\n+\nIIII\n")
        output_dir = tmp_path / "extracted"
        output_dir.mkdir()
        # Create existing extracted file
        (output_dir / "sample_taxid262_R1.fastq.gz").write_text("data")

        with patch(
            "metadarkmatter.cli.kraken2.ExtractKrakenReads.check_available",
            return_value=True,
        ):
            result = runner.invoke(
                app,
                [
                    "kraken2",
                    "extract",
                    "--kraken-output", str(kraken_out),
                    "--reads-1", str(reads),
                    "--taxid", "262",
                    "--output", str(output_dir),
                    "--skip-if-exists",
                ],
            )
        assert result.exit_code == 0

    def test_extract_tool_execution_error(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should handle ToolExecutionError from extract tool."""
        kraken_out = tmp_path / "sample.kraken"
        kraken_out.write_text("C\tread1\t262\t150\t262:150")
        kreport = tmp_path / "sample.kreport"
        kreport.write_text("100.00\t1000\t0\tR\t1\troot\n")
        reads = tmp_path / "reads_R1.fastq.gz"
        reads.write_text("@read1\nACGT\n+\nIIII\n")
        output_dir = tmp_path / "extracted"

        with (
            patch(
                "metadarkmatter.cli.kraken2.ExtractKrakenReads.check_available",
                return_value=True,
            ),
            patch(
                "metadarkmatter.cli.kraken2.ExtractKrakenReads.run_or_raise",
                side_effect=ToolExecutionError(
                    "extract_kraken_reads.py",
                    ["extract_kraken_reads.py"],
                    1,
                    "extraction failed",
                ),
            ),
            patch(
                "metadarkmatter.cli.kraken2.KrakenReport",
            ),
        ):
            result = runner.invoke(
                app,
                [
                    "kraken2",
                    "extract",
                    "--kraken-output", str(kraken_out),
                    "--reads-1", str(reads),
                    "--taxid", "262",
                    "--output", str(output_dir),
                    "--quiet",
                ],
            )
        assert result.exit_code != 0


# ===========================================================================
# Map CLI Tests
# ===========================================================================


class TestMapReads:
    """Tests for map reads command."""

    def test_map_invalid_mode(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should fail with an invalid alignment mode."""
        reads = tmp_path / "reads_R1.fastq"
        reads.write_text("@read1\nACGT\n+\nIIII\n")
        genomes = tmp_path / "genomes"
        genomes.mkdir()
        output = tmp_path / "output.bam"

        with patch(
            "metadarkmatter.cli.map._check_tools_available",
            return_value=(True, []),
        ):
            result = runner.invoke(
                app,
                [
                    "map",
                    "reads",
                    "--reads-1", str(reads),
                    "--genomes", str(genomes),
                    "--output", str(output),
                    "--mode", "invalid_mode",
                ],
            )
        assert result.exit_code != 0

    def test_map_tools_not_available(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should fail when required mapping tools are missing."""
        reads = tmp_path / "reads_R1.fastq"
        reads.write_text("@read1\nACGT\n+\nIIII\n")
        genomes = tmp_path / "genomes"
        genomes.mkdir()
        output = tmp_path / "output.bam"

        with patch(
            "metadarkmatter.cli.map._check_tools_available",
            return_value=(False, ["bowtie2", "samtools"]),
        ):
            result = runner.invoke(
                app,
                [
                    "map",
                    "reads",
                    "--reads-1", str(reads),
                    "--genomes", str(genomes),
                    "--output", str(output),
                ],
            )
        assert result.exit_code != 0

    def test_map_skip_if_exists(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should skip when output BAM already exists."""
        reads = tmp_path / "reads_R1.fastq"
        reads.write_text("@read1\nACGT\n+\nIIII\n")
        genomes = tmp_path / "genomes"
        genomes.mkdir()
        output = tmp_path / "output.bam"
        output.write_text("fake bam")

        with patch(
            "metadarkmatter.cli.map._check_tools_available",
            return_value=(True, []),
        ):
            result = runner.invoke(
                app,
                [
                    "map",
                    "reads",
                    "--reads-1", str(reads),
                    "--genomes", str(genomes),
                    "--output", str(output),
                    "--skip-if-exists",
                ],
            )
        assert result.exit_code == 0

    def test_map_dry_run_with_directory(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Dry run should show commands for genome directory input."""
        reads = tmp_path / "reads_R1.fastq"
        reads.write_text("@read1\nACGT\n+\nIIII\n")
        genomes = tmp_path / "genomes"
        genomes.mkdir()
        (genomes / "genome1.fna").write_text(">contig1\nACGT\n")
        output = tmp_path / "output.bam"

        mock_result = _make_tool_result()

        with (
            patch(
                "metadarkmatter.cli.map._check_tools_available",
                return_value=(True, []),
            ),
            patch(
                "metadarkmatter.cli.map.concatenate_genomes",
                return_value=["genome1"],
            ),
            patch(
                "metadarkmatter.cli.map.Bowtie2Build.run",
                return_value=mock_result,
            ),
            patch(
                "metadarkmatter.cli.map.Bowtie2.run",
                return_value=mock_result,
            ),
        ):
            result = runner.invoke(
                app,
                [
                    "map",
                    "reads",
                    "--reads-1", str(reads),
                    "--genomes", str(genomes),
                    "--output", str(output),
                    "--dry-run",
                ],
            )
        assert result.exit_code == 0
        assert "DRY RUN" in result.output

    def test_map_dry_run_with_existing_index(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Dry run should work with pre-built index."""
        reads = tmp_path / "reads_R1.fastq"
        reads.write_text("@read1\nACGT\n+\nIIII\n")
        pangenome = tmp_path / "pangenome.fasta"
        pangenome.write_text(">contig1\nACGT\n")
        index = tmp_path / "pangenome_db"
        output = tmp_path / "output.bam"

        mock_result = _make_tool_result()

        with (
            patch(
                "metadarkmatter.cli.map._check_tools_available",
                return_value=(True, []),
            ),
            patch(
                "metadarkmatter.cli.map.Bowtie2.run",
                return_value=mock_result,
            ),
        ):
            result = runner.invoke(
                app,
                [
                    "map",
                    "reads",
                    "--reads-1", str(reads),
                    "--genomes", str(pangenome),
                    "--index", str(index),
                    "--output", str(output),
                    "--dry-run",
                ],
            )
        assert result.exit_code == 0

    def test_map_happy_path_with_index(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should run mapping successfully with pre-built index."""
        reads = tmp_path / "reads_R1.fastq"
        reads.write_text("@read1\nACGT\n+\nIIII\n")
        pangenome = tmp_path / "pangenome.fasta"
        pangenome.write_text(">contig1\nACGT\n")
        index = tmp_path / "pangenome_db"
        output = tmp_path / "output.bam"

        mock_result = _make_tool_result()
        mock_sam_results = (mock_result, mock_result, mock_result)

        with (
            patch(
                "metadarkmatter.cli.map._check_tools_available",
                return_value=(True, []),
            ),
            patch(
                "metadarkmatter.cli.map.Bowtie2.run_or_raise",
                return_value=mock_result,
            ),
            patch(
                "metadarkmatter.cli.map.Samtools.sam_to_sorted_bam",
                return_value=mock_sam_results,
            ),
        ):
            result = runner.invoke(
                app,
                [
                    "map",
                    "reads",
                    "--reads-1", str(reads),
                    "--genomes", str(pangenome),
                    "--index", str(index),
                    "--output", str(output),
                    "--quiet",
                ],
            )
        assert result.exit_code == 0

    def test_map_concatenate_error(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should handle FileNotFoundError from concatenate_genomes."""
        reads = tmp_path / "reads_R1.fastq"
        reads.write_text("@read1\nACGT\n+\nIIII\n")
        genomes = tmp_path / "genomes"
        genomes.mkdir()
        output = tmp_path / "output.bam"

        with (
            patch(
                "metadarkmatter.cli.map._check_tools_available",
                return_value=(True, []),
            ),
            patch(
                "metadarkmatter.cli.map.concatenate_genomes",
                side_effect=FileNotFoundError("No genome files found"),
            ),
        ):
            result = runner.invoke(
                app,
                [
                    "map",
                    "reads",
                    "--reads-1", str(reads),
                    "--genomes", str(genomes),
                    "--output", str(output),
                    "--quiet",
                ],
            )
        assert result.exit_code != 0

    def test_map_bowtie2_alignment_error(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should handle ToolExecutionError from Bowtie2 alignment."""
        reads = tmp_path / "reads_R1.fastq"
        reads.write_text("@read1\nACGT\n+\nIIII\n")
        pangenome = tmp_path / "pangenome.fasta"
        pangenome.write_text(">contig1\nACGT\n")
        index = tmp_path / "pangenome_db"
        output = tmp_path / "output.bam"

        with (
            patch(
                "metadarkmatter.cli.map._check_tools_available",
                return_value=(True, []),
            ),
            patch(
                "metadarkmatter.cli.map.Bowtie2.run_or_raise",
                side_effect=ToolExecutionError(
                    "bowtie2", ["bowtie2"], 1, "alignment failed"
                ),
            ),
        ):
            result = runner.invoke(
                app,
                [
                    "map",
                    "reads",
                    "--reads-1", str(reads),
                    "--genomes", str(pangenome),
                    "--index", str(index),
                    "--output", str(output),
                    "--quiet",
                ],
            )
        assert result.exit_code != 0


class TestMapDeriveOutputName:
    """Tests for _derive_output_name helper."""

    def test_common_suffixes(self) -> None:
        from metadarkmatter.cli.map import _derive_output_name

        assert _derive_output_name(Path("sample_R1.fastq.gz")) == "sample"
        assert _derive_output_name(Path("sample_1.fastq.gz")) == "sample"
        assert _derive_output_name(Path("sample.fastq.gz")) == "sample"
        assert _derive_output_name(Path("sample_family_R1.fastq")) == "sample"

    def test_fallback(self) -> None:
        from metadarkmatter.cli.map import _derive_output_name

        assert _derive_output_name(Path("myfile.txt")) == "myfile"


class TestMapCheckToolsAvailable:
    """Tests for _check_tools_available helper."""

    def test_all_available(self) -> None:
        from metadarkmatter.cli.map import _check_tools_available

        with (
            patch("metadarkmatter.cli.map.Bowtie2.check_available", return_value=True),
            patch("metadarkmatter.cli.map.Bowtie2Build.check_available", return_value=True),
            patch("metadarkmatter.cli.map.Samtools.check_available", return_value=True),
        ):
            ok, missing = _check_tools_available()
        assert ok is True
        assert missing == []

    def test_none_available(self) -> None:
        from metadarkmatter.cli.map import _check_tools_available

        with (
            patch("metadarkmatter.cli.map.Bowtie2.check_available", return_value=False),
            patch("metadarkmatter.cli.map.Bowtie2Build.check_available", return_value=False),
            patch("metadarkmatter.cli.map.Samtools.check_available", return_value=False),
        ):
            ok, missing = _check_tools_available()
        assert ok is False
        assert len(missing) == 3


# ===========================================================================
# Proteins CLI Tests
# ===========================================================================


class TestProteinsPredict:
    """Tests for proteins predict command."""

    def test_predict_tool_not_available(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should exit when Prodigal is not installed."""
        genomes = tmp_path / "genomes"
        genomes.mkdir()
        (genomes / "genome1.fna").write_text(">c1\nACGT\n")

        with patch(
            "metadarkmatter.external.prodigal.Prodigal.check_available",
            return_value=False,
        ):
            result = runner.invoke(
                app,
                [
                    "proteins",
                    "predict",
                    "--genomes", str(genomes),
                ],
            )
        assert result.exit_code != 0

    def test_predict_no_genome_files(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should fail when no genome files match the pattern."""
        genomes = tmp_path / "genomes"
        genomes.mkdir()

        with patch(
            "metadarkmatter.external.prodigal.Prodigal.check_available",
            return_value=True,
        ):
            result = runner.invoke(
                app,
                [
                    "proteins",
                    "predict",
                    "--genomes", str(genomes),
                ],
            )
        assert result.exit_code != 0

    def test_predict_dry_run(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Dry run should list genomes that would be processed."""
        genomes = tmp_path / "genomes"
        genomes.mkdir()
        (genomes / "GCF_000001.fna").write_text(">c1\nACGT\n")
        (genomes / "GCF_000002.fna").write_text(">c2\nTGCA\n")

        with patch(
            "metadarkmatter.external.prodigal.Prodigal.check_available",
            return_value=True,
        ):
            result = runner.invoke(
                app,
                [
                    "proteins",
                    "predict",
                    "--genomes", str(genomes),
                    "--dry-run",
                ],
            )
        assert result.exit_code == 0
        assert "DRY RUN" in result.output

    def test_predict_missing_only_all_exist(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should exit cleanly when all proteins already exist."""
        genomes = tmp_path / "genomes"
        genomes.mkdir()
        (genomes / "GCF_000001.fna").write_text(">c1\nACGT\n")
        (genomes / "GCF_000001.faa").write_text(">p1\nM\n")

        with patch(
            "metadarkmatter.external.prodigal.Prodigal.check_available",
            return_value=True,
        ):
            result = runner.invoke(
                app,
                [
                    "proteins",
                    "predict",
                    "--genomes", str(genomes),
                    "--missing-only",
                ],
            )
        assert result.exit_code == 0
        assert "already have protein files" in result.output

    def test_predict_happy_path(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should predict proteins for all genomes."""
        genomes = tmp_path / "genomes"
        genomes.mkdir()
        (genomes / "GCF_000001.fna").write_text(">c1\nACGT\n")

        def fake_predict(genome_file, output_dir, procedure):
            # Create a fake protein file
            accession = genome_file.stem
            protein_file = output_dir / f"{accession}.faa"
            protein_file.write_text(">prot1\nMAAA\n>prot2\nMBBB\n")
            return (genome_file, 2, None)

        with (
            patch(
                "metadarkmatter.external.prodigal.Prodigal.check_available",
                return_value=True,
            ),
            patch(
                "metadarkmatter.cli.proteins._predict_single_genome",
                side_effect=fake_predict,
            ),
        ):
            result = runner.invoke(
                app,
                [
                    "proteins",
                    "predict",
                    "--genomes", str(genomes),
                    "--quiet",
                ],
            )
        assert result.exit_code == 0

    def test_predict_with_errors(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should report failed genomes in the summary."""
        genomes = tmp_path / "genomes"
        genomes.mkdir()
        (genomes / "GCF_000001.fna").write_text(">c1\nACGT\n")

        def fake_predict_error(genome_file, output_dir, procedure):
            return (genome_file, 0, "Prodigal failed: too short")

        with (
            patch(
                "metadarkmatter.external.prodigal.Prodigal.check_available",
                return_value=True,
            ),
            patch(
                "metadarkmatter.cli.proteins._predict_single_genome",
                side_effect=fake_predict_error,
            ),
        ):
            result = runner.invoke(
                app,
                [
                    "proteins",
                    "predict",
                    "--genomes", str(genomes),
                    "--quiet",
                ],
            )
        assert result.exit_code == 0


class TestProteinHelpers:
    """Tests for proteins module helper functions."""

    def test_find_genome_files_primary_pattern(self, tmp_path: Path) -> None:
        from metadarkmatter.cli.proteins import _find_genome_files

        (tmp_path / "g1.fna").write_text(">c1\nACGT\n")
        (tmp_path / "g2.fna").write_text(">c2\nTGCA\n")
        files = _find_genome_files(tmp_path, "*.fna")
        assert len(files) == 2

    def test_find_genome_files_fallback_pattern(self, tmp_path: Path) -> None:
        from metadarkmatter.cli.proteins import _find_genome_files

        (tmp_path / "g1.fasta").write_text(">c1\nACGT\n")
        files = _find_genome_files(tmp_path, "*.xyz")
        assert len(files) == 1

    def test_find_missing_proteins(self, tmp_path: Path) -> None:
        from metadarkmatter.cli.proteins import _find_missing_proteins

        genome_dir = tmp_path / "genomes"
        genome_dir.mkdir()
        (genome_dir / "GCF_001.fna").write_text(">c1\nACGT\n")
        (genome_dir / "GCF_002.fna").write_text(">c2\nACGT\n")
        # Only GCF_001 has proteins
        (genome_dir / "GCF_001.faa").write_text(">p1\nM\n")

        missing, existing = _find_missing_proteins(genome_dir, genome_dir, "*.fna")
        assert len(missing) == 1
        assert len(existing) == 1
        assert missing[0].stem == "GCF_002"


# ===========================================================================
# BLASTX CLI Tests
# ===========================================================================


class TestBlastxMakedb:
    """Tests for blastx makedb command."""

    def test_makedb_tool_not_available(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should exit when Diamond is not installed."""
        proteins = tmp_path / "proteins.faa"
        proteins.write_text(">p1\nMAAA\n")
        output = tmp_path / "db" / "panproteome"

        with patch(
            "metadarkmatter.cli.blastx.Diamond.check_available",
            return_value=False,
        ):
            result = runner.invoke(
                app,
                [
                    "blastx",
                    "makedb",
                    "--proteins", str(proteins),
                    "--output", str(output),
                ],
            )
        assert result.exit_code != 0

    def test_makedb_skip_if_exists(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should skip when database already exists."""
        proteins = tmp_path / "proteins.faa"
        proteins.write_text(">p1\nMAAA\n")
        output = tmp_path / "panproteome"
        (tmp_path / "panproteome.dmnd").write_text("fake db")

        with patch(
            "metadarkmatter.cli.blastx.Diamond.check_available",
            return_value=True,
        ):
            result = runner.invoke(
                app,
                [
                    "blastx",
                    "makedb",
                    "--proteins", str(proteins),
                    "--output", str(output),
                    "--skip-if-exists",
                ],
            )
        assert result.exit_code == 0

    def test_makedb_dry_run_single_file(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Dry run with a single protein FASTA file."""
        proteins = tmp_path / "proteins.faa"
        proteins.write_text(">p1\nMAAA\n")
        output = tmp_path / "panproteome"

        mock_result = _make_tool_result()

        with (
            patch(
                "metadarkmatter.cli.blastx.Diamond.check_available",
                return_value=True,
            ),
            patch(
                "metadarkmatter.cli.blastx.Diamond.run",
                return_value=mock_result,
            ),
        ):
            result = runner.invoke(
                app,
                [
                    "blastx",
                    "makedb",
                    "--proteins", str(proteins),
                    "--output", str(output),
                    "--dry-run",
                ],
            )
        assert result.exit_code == 0

    def test_makedb_dry_run_directory(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Dry run with a directory of protein FASTA files."""
        protein_dir = tmp_path / "proteins"
        protein_dir.mkdir()
        (protein_dir / "GCF_001.faa").write_text(">p1\nMAAA\n")
        (protein_dir / "GCF_002.faa").write_text(">p2\nMBBB\n")
        output = tmp_path / "panproteome"

        mock_result = _make_tool_result()

        with (
            patch(
                "metadarkmatter.cli.blastx.Diamond.check_available",
                return_value=True,
            ),
            patch(
                "metadarkmatter.cli.blastx.Diamond.run",
                return_value=mock_result,
            ),
        ):
            result = runner.invoke(
                app,
                [
                    "blastx",
                    "makedb",
                    "--proteins", str(protein_dir),
                    "--output", str(output),
                    "--dry-run",
                ],
            )
        assert result.exit_code == 0

    def test_makedb_happy_path_single_file(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should create database from a single protein FASTA."""
        proteins = tmp_path / "proteins.faa"
        proteins.write_text(">p1\nMAAA\n")
        output = tmp_path / "panproteome"

        mock_result = _make_tool_result()

        def fake_run_or_raise(**kwargs):
            # Create the expected .dmnd file
            db_path = output.with_suffix(".dmnd")
            db_path.write_text("fake diamond db")
            return mock_result

        with (
            patch(
                "metadarkmatter.cli.blastx.Diamond.check_available",
                return_value=True,
            ),
            patch(
                "metadarkmatter.cli.blastx.Diamond.run_or_raise",
                side_effect=fake_run_or_raise,
            ),
        ):
            result = runner.invoke(
                app,
                [
                    "blastx",
                    "makedb",
                    "--proteins", str(proteins),
                    "--output", str(output),
                    "--quiet",
                ],
            )
        assert result.exit_code == 0

    def test_makedb_happy_path_directory(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should concatenate and create database from protein directory."""
        protein_dir = tmp_path / "proteins"
        protein_dir.mkdir()
        (protein_dir / "GCF_001.faa").write_text(">p1\nMAAA\n")
        output = tmp_path / "panproteome"

        mock_result = _make_tool_result()

        def fake_run_or_raise(**kwargs):
            db_path = output.with_suffix(".dmnd")
            db_path.write_text("fake diamond db")
            return mock_result

        with (
            patch(
                "metadarkmatter.cli.blastx.Diamond.check_available",
                return_value=True,
            ),
            patch(
                "metadarkmatter.cli.blastx.Diamond.run_or_raise",
                side_effect=fake_run_or_raise,
            ),
        ):
            result = runner.invoke(
                app,
                [
                    "blastx",
                    "makedb",
                    "--proteins", str(protein_dir),
                    "--output", str(output),
                    "--quiet",
                ],
            )
        assert result.exit_code == 0

    def test_makedb_empty_directory(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should fail when protein directory is empty."""
        protein_dir = tmp_path / "proteins"
        protein_dir.mkdir()
        output = tmp_path / "panproteome"

        with patch(
            "metadarkmatter.cli.blastx.Diamond.check_available",
            return_value=True,
        ):
            result = runner.invoke(
                app,
                [
                    "blastx",
                    "makedb",
                    "--proteins", str(protein_dir),
                    "--output", str(output),
                ],
            )
        assert result.exit_code != 0

    def test_makedb_tool_execution_error(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should handle ToolExecutionError from Diamond makedb."""
        proteins = tmp_path / "proteins.faa"
        proteins.write_text(">p1\nMAAA\n")
        output = tmp_path / "panproteome"

        with (
            patch(
                "metadarkmatter.cli.blastx.Diamond.check_available",
                return_value=True,
            ),
            patch(
                "metadarkmatter.cli.blastx.Diamond.run_or_raise",
                side_effect=ToolExecutionError(
                    "diamond", ["diamond", "makedb"], 1, "makedb failed"
                ),
            ),
        ):
            result = runner.invoke(
                app,
                [
                    "blastx",
                    "makedb",
                    "--proteins", str(proteins),
                    "--output", str(output),
                    "--quiet",
                ],
            )
        assert result.exit_code != 0


class TestBlastxAlign:
    """Tests for blastx align command."""

    def test_align_tool_not_available(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should exit when Diamond is not installed."""
        query = tmp_path / "reads.fastq.gz"
        query.write_text("@read1\nACGT\n+\nIIII\n")
        db = tmp_path / "panproteome"

        with patch(
            "metadarkmatter.cli.blastx.Diamond.check_available",
            return_value=False,
        ):
            result = runner.invoke(
                app,
                [
                    "blastx",
                    "align",
                    "--query", str(query),
                    "--database", str(db),
                    "--output", str(tmp_path / "results.tsv.gz"),
                ],
            )
        assert result.exit_code != 0

    def test_align_missing_database(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should fail when database file is missing."""
        query = tmp_path / "reads.fastq.gz"
        query.write_text("@read1\nACGT\n+\nIIII\n")
        db = tmp_path / "nonexistent_db"
        output = tmp_path / "results.tsv.gz"

        with patch(
            "metadarkmatter.cli.blastx.Diamond.check_available",
            return_value=True,
        ):
            result = runner.invoke(
                app,
                [
                    "blastx",
                    "align",
                    "--query", str(query),
                    "--database", str(db),
                    "--output", str(output),
                ],
            )
        assert result.exit_code != 0

    def test_align_skip_if_exists(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should skip when output already exists."""
        query = tmp_path / "reads.fastq.gz"
        query.write_text("@read1\nACGT\n+\nIIII\n")
        db = tmp_path / "panproteome"
        (tmp_path / "panproteome.dmnd").write_text("fake db")
        output = tmp_path / "results.tsv.gz"
        output.write_text("existing results")

        with patch(
            "metadarkmatter.cli.blastx.Diamond.check_available",
            return_value=True,
        ):
            result = runner.invoke(
                app,
                [
                    "blastx",
                    "align",
                    "--query", str(query),
                    "--database", str(db),
                    "--output", str(output),
                    "--skip-if-exists",
                ],
            )
        assert result.exit_code == 0

    def test_align_dry_run(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Dry run should show command without executing."""
        query = tmp_path / "reads.fastq.gz"
        query.write_text("@read1\nACGT\n+\nIIII\n")
        db = tmp_path / "panproteome"
        (tmp_path / "panproteome.dmnd").write_text("fake db")
        output = tmp_path / "results.tsv.gz"

        mock_result = _make_tool_result()

        with (
            patch(
                "metadarkmatter.cli.blastx.Diamond.check_available",
                return_value=True,
            ),
            patch(
                "metadarkmatter.cli.blastx.Diamond.run",
                return_value=mock_result,
            ),
        ):
            result = runner.invoke(
                app,
                [
                    "blastx",
                    "align",
                    "--query", str(query),
                    "--database", str(db),
                    "--output", str(output),
                    "--dry-run",
                ],
            )
        assert result.exit_code == 0

    def test_align_happy_path(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should run alignment and produce compressed output."""
        query = tmp_path / "reads.fastq.gz"
        query.write_text("@read1\nACGT\n+\nIIII\n")
        db = tmp_path / "panproteome"
        (tmp_path / "panproteome.dmnd").write_text("fake db")
        output = tmp_path / "results.tsv.gz"

        mock_result = _make_tool_result()

        def fake_run_or_raise(**kwargs):
            # Create the temp output file that Diamond would create
            temp_out = kwargs.get("output")
            if temp_out and isinstance(temp_out, Path):
                temp_out.parent.mkdir(parents=True, exist_ok=True)
                temp_out.write_text(
                    "read1\tGCF_001|p1\t85.0\t50\t7\t0\t1\t150\t1\t50\t1e-10\t100.0\n"
                )
            return mock_result

        with (
            patch(
                "metadarkmatter.cli.blastx.Diamond.check_available",
                return_value=True,
            ),
            patch(
                "metadarkmatter.cli.blastx.Diamond.run_or_raise",
                side_effect=fake_run_or_raise,
            ),
        ):
            result = runner.invoke(
                app,
                [
                    "blastx",
                    "align",
                    "--query", str(query),
                    "--database", str(db),
                    "--output", str(output),
                    "--quiet",
                ],
            )
        assert result.exit_code == 0

    def test_align_tool_execution_error(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should handle ToolExecutionError from Diamond blastx."""
        query = tmp_path / "reads.fastq.gz"
        query.write_text("@read1\nACGT\n+\nIIII\n")
        db = tmp_path / "panproteome"
        (tmp_path / "panproteome.dmnd").write_text("fake db")
        output = tmp_path / "results.tsv.gz"

        with (
            patch(
                "metadarkmatter.cli.blastx.Diamond.check_available",
                return_value=True,
            ),
            patch(
                "metadarkmatter.cli.blastx.Diamond.run_or_raise",
                side_effect=ToolExecutionError(
                    "diamond", ["diamond", "blastx"], 1, "blastx failed"
                ),
            ),
        ):
            result = runner.invoke(
                app,
                [
                    "blastx",
                    "align",
                    "--query", str(query),
                    "--database", str(db),
                    "--output", str(output),
                    "--quiet",
                ],
            )
        assert result.exit_code != 0


class TestBlastxConcatenateProteins:
    """Tests for _concatenate_proteins helper."""

    def test_concatenate_basic(self, tmp_path: Path) -> None:
        from metadarkmatter.cli.blastx import _concatenate_proteins

        protein_dir = tmp_path / "proteins"
        protein_dir.mkdir()
        (protein_dir / "GCF_001.faa").write_text(">prot_A desc\nMAAA\n>prot_B\nMBBB\n")
        (protein_dir / "GCF_002.faa").write_text(">prot_C\nMCCC\n")

        output = tmp_path / "panproteome.faa"
        genome_count, protein_count = _concatenate_proteins(protein_dir, output)

        assert genome_count == 2
        assert protein_count == 3
        content = output.read_text()
        assert ">GCF_001|prot_A" in content
        assert ">GCF_002|prot_C" in content

    def test_concatenate_empty_dir(self, tmp_path: Path) -> None:
        from metadarkmatter.cli.blastx import _concatenate_proteins

        protein_dir = tmp_path / "empty"
        protein_dir.mkdir()
        output = tmp_path / "out.faa"
        genome_count, protein_count = _concatenate_proteins(protein_dir, output)
        assert genome_count == 0
        assert protein_count == 0


class TestBlastxCheckDiamond:
    """Tests for _check_diamond_available helper."""

    def test_available(self) -> None:
        from metadarkmatter.cli.blastx import _check_diamond_available

        with patch(
            "metadarkmatter.cli.blastx.Diamond.check_available",
            return_value=True,
        ):
            ok, missing = _check_diamond_available()
        assert ok is True
        assert missing == ""

    def test_not_available(self) -> None:
        from metadarkmatter.cli.blastx import _check_diamond_available

        with patch(
            "metadarkmatter.cli.blastx.Diamond.check_available",
            return_value=False,
        ):
            ok, missing = _check_diamond_available()
        assert ok is False
        assert missing == "diamond"


# ===========================================================================
# AAI CLI Tests
# ===========================================================================


class TestAAICompute:
    """Tests for aai compute command."""

    def test_compute_tool_not_available(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should exit when Diamond is not installed."""
        genomes = tmp_path / "genomes"
        genomes.mkdir()
        (genomes / "GCF_001.faa").write_text(">p1\nMAAA\n")

        with patch(
            "metadarkmatter.external.diamond.Diamond.check_available",
            return_value=False,
        ):
            result = runner.invoke(
                app,
                [
                    "aai",
                    "compute",
                    "--genomes", str(genomes),
                    "--output", str(tmp_path / "aai_matrix.csv"),
                ],
            )
        assert result.exit_code != 0

    def test_compute_no_protein_files(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should fail when no protein files are found."""
        genomes = tmp_path / "genomes"
        genomes.mkdir()

        with patch(
            "metadarkmatter.external.diamond.Diamond.check_available",
            return_value=True,
        ):
            result = runner.invoke(
                app,
                [
                    "aai",
                    "compute",
                    "--genomes", str(genomes),
                    "--output", str(tmp_path / "aai_matrix.csv"),
                ],
            )
        assert result.exit_code != 0

    def test_compute_dry_run(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Dry run should show steps without execution."""
        genomes = tmp_path / "genomes"
        genomes.mkdir()
        (genomes / "GCF_001.faa").write_text(">p1\nMAAA\n")
        (genomes / "GCF_002.faa").write_text(">p2\nMBBB\n")

        with patch(
            "metadarkmatter.external.diamond.Diamond.check_available",
            return_value=True,
        ):
            result = runner.invoke(
                app,
                [
                    "aai",
                    "compute",
                    "--genomes", str(genomes),
                    "--output", str(tmp_path / "aai_matrix.csv"),
                    "--dry-run",
                ],
            )
        assert result.exit_code == 0
        assert "DRY RUN" in result.output

    def test_compute_happy_path(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should compute AAI matrix successfully."""
        genomes = tmp_path / "genomes"
        genomes.mkdir()
        (genomes / "GCF_001.faa").write_text(">p1\nMAAA\n")
        (genomes / "GCF_002.faa").write_text(">p2\nMBBB\n")
        output = tmp_path / "aai_matrix.csv"

        # Create a fake AAI matrix output that compute_aai_matrix would produce
        fake_aai_dict = {
            "GCF_001": {"GCF_001": 100.0, "GCF_002": 75.0},
            "GCF_002": {"GCF_001": 75.0, "GCF_002": 100.0},
        }

        def fake_compute(**kwargs):
            # Write a real CSV so the summary code can read it
            out_path = kwargs.get("output_path", output)
            out_path.write_text(
                "genome,GCF_001,GCF_002\n"
                "GCF_001,100.0,75.0\n"
                "GCF_002,75.0,100.0\n"
            )
            return fake_aai_dict

        with (
            patch(
                "metadarkmatter.external.diamond.Diamond.check_available",
                return_value=True,
            ),
            patch(
                "metadarkmatter.core.aai_matrix_builder.compute_aai_matrix",
                side_effect=fake_compute,
            ),
        ):
            result = runner.invoke(
                app,
                [
                    "aai",
                    "compute",
                    "--genomes", str(genomes),
                    "--output", str(output),
                    "--quiet",
                ],
            )
        assert result.exit_code == 0

    def test_compute_tool_execution_error(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should handle ToolExecutionError from Diamond."""
        genomes = tmp_path / "genomes"
        genomes.mkdir()
        (genomes / "GCF_001.faa").write_text(">p1\nMAAA\n")

        with (
            patch(
                "metadarkmatter.external.diamond.Diamond.check_available",
                return_value=True,
            ),
            patch(
                "metadarkmatter.core.aai_matrix_builder.compute_aai_matrix",
                side_effect=ToolExecutionError(
                    "diamond", ["diamond", "blastp"], 1, "diamond failed"
                ),
            ),
        ):
            result = runner.invoke(
                app,
                [
                    "aai",
                    "compute",
                    "--genomes", str(genomes),
                    "--output", str(tmp_path / "aai_matrix.csv"),
                    "--quiet",
                ],
            )
        assert result.exit_code != 0


class TestAAIValidate:
    """Tests for aai validate command."""

    def _create_matrix_csv(self, path: Path, genomes: list[str]) -> None:
        """Helper to write a CSV matrix file."""
        header = "genome," + ",".join(genomes)
        rows = []
        for g1 in genomes:
            values = []
            for g2 in genomes:
                values.append("100.0" if g1 == g2 else "75.0")
            rows.append(f"{g1}," + ",".join(values))
        path.write_text(header + "\n" + "\n".join(rows) + "\n")

    def test_validate_good_coverage(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should report good coverage when both matrices match."""
        genomes = ["GCF_001", "GCF_002", "GCF_003"]
        ani_path = tmp_path / "ani.csv"
        aai_path = tmp_path / "aai.csv"
        self._create_matrix_csv(ani_path, genomes)
        self._create_matrix_csv(aai_path, genomes)

        result = runner.invoke(
            app,
            [
                "aai",
                "validate",
                "--aai", str(aai_path),
                "--ani", str(ani_path),
            ],
        )
        assert result.exit_code == 0
        assert "100.0%" in result.output

    def test_validate_low_coverage(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should warn and return non-zero exit for low coverage."""
        ani_genomes = ["GCF_001", "GCF_002", "GCF_003", "GCF_004"]
        aai_genomes = ["GCF_001"]
        ani_path = tmp_path / "ani.csv"
        aai_path = tmp_path / "aai.csv"
        self._create_matrix_csv(ani_path, ani_genomes)
        self._create_matrix_csv(aai_path, aai_genomes)

        result = runner.invoke(
            app,
            [
                "aai",
                "validate",
                "--aai", str(aai_path),
                "--ani", str(ani_path),
            ],
        )
        assert result.exit_code != 0

    def test_validate_moderate_coverage(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Should report moderate coverage (50-90%)."""
        ani_genomes = ["GCF_001", "GCF_002", "GCF_003"]
        aai_genomes = ["GCF_001", "GCF_002"]
        ani_path = tmp_path / "ani.csv"
        aai_path = tmp_path / "aai.csv"
        self._create_matrix_csv(ani_path, ani_genomes)
        self._create_matrix_csv(aai_path, aai_genomes)

        result = runner.invoke(
            app,
            [
                "aai",
                "validate",
                "--aai", str(aai_path),
                "--ani", str(ani_path),
            ],
        )
        assert result.exit_code == 0
        assert "66.7%" in result.output

    def test_validate_verbose(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """Verbose should show missing genomes and AAI statistics."""
        ani_genomes = ["GCF_001", "GCF_002", "GCF_003"]
        aai_genomes = ["GCF_001", "GCF_002", "GCF_003"]
        ani_path = tmp_path / "ani.csv"
        aai_path = tmp_path / "aai.csv"
        self._create_matrix_csv(ani_path, ani_genomes)
        self._create_matrix_csv(aai_path, aai_genomes)

        result = runner.invoke(
            app,
            [
                "aai",
                "validate",
                "--aai", str(aai_path),
                "--ani", str(ani_path),
                "--verbose",
            ],
        )
        assert result.exit_code == 0


class TestAAIFindProteinFiles:
    """Tests for _find_protein_files helper."""

    def test_primary_pattern(self, tmp_path: Path) -> None:
        from metadarkmatter.cli.aai import _find_protein_files

        (tmp_path / "g1.faa").write_text(">p1\nM\n")
        (tmp_path / "g2.faa").write_text(">p2\nM\n")
        files = _find_protein_files(tmp_path, "*.faa")
        assert len(files) == 2

    def test_fallback_pattern(self, tmp_path: Path) -> None:
        from metadarkmatter.cli.aai import _find_protein_files

        (tmp_path / "g1.fasta").write_text(">p1\nM\n")
        files = _find_protein_files(tmp_path, "*.xyz")
        assert len(files) == 1

    def test_no_files(self, tmp_path: Path) -> None:
        from metadarkmatter.cli.aai import _find_protein_files

        files = _find_protein_files(tmp_path, "*.faa")
        assert len(files) == 0
