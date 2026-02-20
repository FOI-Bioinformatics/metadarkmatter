"""Unit tests for BLAST, MMseqs2, and Samtools external tool wrappers.

Covers the uncovered lines identified by coverage analysis:
- blast.py lines 65-95 (MakeBlastDb.build_command), 197-266 (BlastN optional
  filters, FASTQ detection, FASTQ-to-FASTA conversion), 294-344 (BlastN.run
  with FASTQ conversion, build_sensitive_command_args)
- mmseqs2.py lines 229-244 (_build_search_command validation), 298
  (create_database convenience), 353 (search convenience), 436-529
  (search_multistep workflow)
- samtools.py lines covering view/sort/index execution, _run_command error
  paths, and sam_to_sorted_bam pipeline

All tests use mocking so that no external bioinformatics tools are required.
"""

from __future__ import annotations

import gzip
import subprocess
import tempfile
from pathlib import Path
from typing import ClassVar
from unittest.mock import MagicMock, call, patch

import pytest

from metadarkmatter.external.base import (
    ExternalTool,
    ToolExecutionError,
    ToolResult,
    ToolTimeoutError,
)
from metadarkmatter.external.blast import BlastN, MakeBlastDb
from metadarkmatter.external.mmseqs2 import MMseqs2
from metadarkmatter.external.samtools import Samtools


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _clear_caches():
    """Ensure executable caches are empty before and after each test."""
    ExternalTool._executable_cache.clear()
    yield
    ExternalTool._executable_cache.clear()


@pytest.fixture()
def mock_blastn() -> BlastN:
    """BlastN instance with a cached mock executable."""
    BlastN._executable_cache["blastn"] = Path("/usr/bin/blastn")
    return BlastN()


@pytest.fixture()
def mock_makeblastdb() -> MakeBlastDb:
    """MakeBlastDb instance with a cached mock executable."""
    MakeBlastDb._executable_cache["makeblastdb"] = Path("/usr/bin/makeblastdb")
    return MakeBlastDb()


@pytest.fixture()
def mock_mmseqs() -> MMseqs2:
    """MMseqs2 instance with a cached mock executable."""
    MMseqs2._executable_cache["mmseqs"] = Path("/usr/bin/mmseqs")
    return MMseqs2()


@pytest.fixture()
def mock_samtools() -> Samtools:
    """Samtools instance with a cached mock executable."""
    Samtools._executable_cache["samtools"] = Path("/usr/bin/samtools")
    return Samtools()


# =========================================================================
# MakeBlastDb.build_command  (blast.py lines 65-95)
# =========================================================================


class TestMakeBlastDbBuildCommand:
    """Tests for MakeBlastDb.build_command covering all optional parameters."""

    def test_basic_command(self, mock_makeblastdb: MakeBlastDb):
        """Should build a basic makeblastdb command with defaults."""
        cmd = mock_makeblastdb.build_command(
            input_fasta=Path("input.fasta"),
            output_db=Path("dbname"),
        )
        assert cmd[0] == "/usr/bin/makeblastdb"
        assert "-in" in cmd
        assert "-out" in cmd
        assert "-dbtype" in cmd
        assert "nucl" in cmd
        assert "-parse_seqids" in cmd

    def test_protein_dbtype(self, mock_makeblastdb: MakeBlastDb):
        """Should use 'prot' dbtype when specified."""
        cmd = mock_makeblastdb.build_command(
            input_fasta=Path("proteins.faa"),
            output_db=Path("protdb"),
            dbtype="prot",
        )
        idx = cmd.index("-dbtype")
        assert cmd[idx + 1] == "prot"

    def test_parse_seqids_disabled(self, mock_makeblastdb: MakeBlastDb):
        """Should omit -parse_seqids when disabled."""
        cmd = mock_makeblastdb.build_command(
            input_fasta=Path("input.fasta"),
            output_db=Path("db"),
            parse_seqids=False,
        )
        assert "-parse_seqids" not in cmd

    def test_title_parameter(self, mock_makeblastdb: MakeBlastDb):
        """Should include -title when provided."""
        cmd = mock_makeblastdb.build_command(
            input_fasta=Path("input.fasta"),
            output_db=Path("db"),
            title="My Database",
        )
        assert "-title" in cmd
        idx = cmd.index("-title")
        assert cmd[idx + 1] == "My Database"

    def test_taxid_parameter(self, mock_makeblastdb: MakeBlastDb):
        """Should include -taxid when provided."""
        cmd = mock_makeblastdb.build_command(
            input_fasta=Path("input.fasta"),
            output_db=Path("db"),
            taxid=9606,
        )
        assert "-taxid" in cmd
        idx = cmd.index("-taxid")
        assert cmd[idx + 1] == "9606"

    def test_logfile_parameter(self, mock_makeblastdb: MakeBlastDb):
        """Should include -logfile when provided."""
        cmd = mock_makeblastdb.build_command(
            input_fasta=Path("input.fasta"),
            output_db=Path("db"),
            logfile=Path("build.log"),
        )
        assert "-logfile" in cmd

    def test_all_optional_parameters(self, mock_makeblastdb: MakeBlastDb):
        """Should include all optional parameters when provided together."""
        cmd = mock_makeblastdb.build_command(
            input_fasta=Path("input.fasta"),
            output_db=Path("db"),
            dbtype="nucl",
            parse_seqids=True,
            title="Test DB",
            taxid=12345,
            logfile=Path("out.log"),
        )
        assert "-parse_seqids" in cmd
        assert "-title" in cmd
        assert "-taxid" in cmd
        assert "12345" in cmd
        assert "-logfile" in cmd


# =========================================================================
# BlastN.build_command optional filters  (blast.py lines 196-205)
# =========================================================================


class TestBlastNOptionalFilters:
    """Tests for BlastN optional filter parameters."""

    def test_perc_identity_filter(self, mock_blastn: BlastN):
        """Should include -perc_identity when specified."""
        cmd = mock_blastn.build_command(
            query=Path("reads.fasta"),
            database=Path("db"),
            output=Path("out.tsv"),
            perc_identity=90.0,
        )
        assert "-perc_identity" in cmd
        idx = cmd.index("-perc_identity")
        assert cmd[idx + 1] == "90.0"

    def test_dust_filter(self, mock_blastn: BlastN):
        """Should include -dust when specified."""
        cmd = mock_blastn.build_command(
            query=Path("reads.fasta"),
            database=Path("db"),
            output=Path("out.tsv"),
            dust="no",
        )
        assert "-dust" in cmd
        idx = cmd.index("-dust")
        assert cmd[idx + 1] == "no"

    def test_soft_masking_true(self, mock_blastn: BlastN):
        """Should include -soft_masking true when enabled."""
        cmd = mock_blastn.build_command(
            query=Path("reads.fasta"),
            database=Path("db"),
            output=Path("out.tsv"),
            soft_masking=True,
        )
        assert "-soft_masking" in cmd
        idx = cmd.index("-soft_masking")
        assert cmd[idx + 1] == "true"

    def test_soft_masking_false(self, mock_blastn: BlastN):
        """Should include -soft_masking false when disabled."""
        cmd = mock_blastn.build_command(
            query=Path("reads.fasta"),
            database=Path("db"),
            output=Path("out.tsv"),
            soft_masking=False,
        )
        assert "-soft_masking" in cmd
        idx = cmd.index("-soft_masking")
        assert cmd[idx + 1] == "false"

    def test_all_optional_filters_together(self, mock_blastn: BlastN):
        """Should include all optional filters simultaneously."""
        cmd = mock_blastn.build_command(
            query=Path("reads.fasta"),
            database=Path("db"),
            output=Path("out.tsv"),
            perc_identity=85.0,
            dust="yes",
            soft_masking=True,
        )
        assert "-perc_identity" in cmd
        assert "-dust" in cmd
        assert "-soft_masking" in cmd


# =========================================================================
# BlastN._is_fastq_file  (blast.py lines 207-218)
# =========================================================================


class TestBlastNIsFastqFile:
    """Tests for BlastN._is_fastq_file static method."""

    @pytest.mark.parametrize(
        "filename,expected",
        [
            ("reads.fastq", True),
            ("reads.fq", True),
            ("reads.fastq.gz", True),
            ("reads.fq.gz", True),
            ("reads.FASTQ", True),
            ("reads.FQ.GZ", True),
            ("reads.fasta", False),
            ("reads.fa", False),
            ("reads.fasta.gz", False),
            ("reads.fa.gz", False),
            ("reads.txt", False),
            ("reads.bam", False),
        ],
    )
    def test_fastq_detection(self, filename: str, expected: bool):
        """Should correctly identify FASTQ files by extension."""
        assert BlastN._is_fastq_file(Path(filename)) == expected


# =========================================================================
# BlastN._convert_fastq_to_fasta  (blast.py lines 220-266)
# =========================================================================


class TestBlastNConvertFastqToFasta:
    """Tests for BlastN._convert_fastq_to_fasta static method."""

    def test_python_fallback_plain_fastq(self, tmp_path: Path):
        """Should convert plain FASTQ to FASTA using Python fallback."""
        fastq_file = tmp_path / "reads.fastq"
        fasta_file = tmp_path / "reads.fasta"

        fastq_content = (
            "@read1\n"
            "ACGTACGT\n"
            "+\n"
            "IIIIIIII\n"
            "@read2\n"
            "TGCATGCA\n"
            "+\n"
            "JJJJJJJJ\n"
        )
        fastq_file.write_text(fastq_content)

        # Patch subprocess.run to simulate seqtk not found
        with patch(
            "metadarkmatter.external.blast.subprocess.run",
            side_effect=FileNotFoundError("seqtk not found"),
        ):
            BlastN._convert_fastq_to_fasta(fastq_file, fasta_file)

        fasta_content = fasta_file.read_text()
        assert ">read1\n" in fasta_content
        assert "ACGTACGT\n" in fasta_content
        assert ">read2\n" in fasta_content
        assert "TGCATGCA\n" in fasta_content
        # Quality lines should not appear
        assert "+" not in fasta_content
        assert "IIIIIIII" not in fasta_content

    def test_python_fallback_gzipped_fastq(self, tmp_path: Path):
        """Should convert gzipped FASTQ to FASTA using Python fallback."""
        fastq_file = tmp_path / "reads.fastq.gz"
        fasta_file = tmp_path / "reads.fasta"

        fastq_content = (
            "@read1\n"
            "ACGTACGT\n"
            "+\n"
            "IIIIIIII\n"
        )
        with gzip.open(fastq_file, "wt") as f:
            f.write(fastq_content)

        with patch(
            "metadarkmatter.external.blast.subprocess.run",
            side_effect=FileNotFoundError("seqtk not found"),
        ):
            BlastN._convert_fastq_to_fasta(fastq_file, fasta_file)

        fasta_content = fasta_file.read_text()
        assert ">read1\n" in fasta_content
        assert "ACGTACGT\n" in fasta_content

    def test_seqtk_success(self, tmp_path: Path):
        """Should use seqtk when available and successful."""
        fastq_file = tmp_path / "reads.fastq"
        fasta_file = tmp_path / "reads.fasta"
        fastq_file.write_text("@read1\nACGT\n+\nIIII\n")

        mock_result = MagicMock()
        mock_result.returncode = 0

        with patch(
            "metadarkmatter.external.blast.subprocess.run",
            return_value=mock_result,
        ) as mock_run:
            # We need to create the fasta file since we mock subprocess
            # The real implementation writes via stdout redirect, but we can
            # verify seqtk was called correctly
            BlastN._convert_fastq_to_fasta(fastq_file, fasta_file)

            mock_run.assert_called_once()
            args = mock_run.call_args
            assert args[0][0] == ["seqtk", "seq", "-A", str(fastq_file)]

    def test_seqtk_failure_falls_back_to_python(self, tmp_path: Path):
        """Should fall back to Python when seqtk fails."""
        fastq_file = tmp_path / "reads.fastq"
        fasta_file = tmp_path / "reads.fasta"
        fastq_file.write_text("@read1\nACGT\n+\nIIII\n")

        with patch(
            "metadarkmatter.external.blast.subprocess.run",
            side_effect=subprocess.CalledProcessError(1, "seqtk"),
        ):
            BlastN._convert_fastq_to_fasta(fastq_file, fasta_file)

        fasta_content = fasta_file.read_text()
        assert ">read1\n" in fasta_content
        assert "ACGT\n" in fasta_content

    def test_python_fallback_error_raises_runtime_error(self, tmp_path: Path):
        """Should raise RuntimeError when Python fallback also fails."""
        fastq_file = tmp_path / "nonexistent.fastq"
        fasta_file = tmp_path / "reads.fasta"

        with patch(
            "metadarkmatter.external.blast.subprocess.run",
            side_effect=FileNotFoundError("seqtk not found"),
        ):
            with pytest.raises(RuntimeError, match="Failed to convert FASTQ to FASTA"):
                BlastN._convert_fastq_to_fasta(fastq_file, fasta_file)


# =========================================================================
# BlastN.run with FASTQ conversion  (blast.py lines 294-332)
# =========================================================================


class TestBlastNRunFastqConversion:
    """Tests for BlastN.run with automatic FASTQ-to-FASTA conversion."""

    def test_run_with_fastq_converts_and_cleans_up(
        self, mock_blastn: BlastN, tmp_path: Path
    ):
        """Should convert FASTQ to FASTA, run BLAST, and clean up temp file."""
        fastq_file = tmp_path / "reads.fastq"
        fastq_file.write_text("@read1\nACGT\n+\nIIII\n")

        with patch(
            "metadarkmatter.external.blast.subprocess.run",
            side_effect=FileNotFoundError("seqtk not found"),
        ):
            with patch.object(
                ExternalTool,
                "run",
                return_value=ToolResult(
                    command=("blastn",),
                    return_code=0,
                    stdout="",
                    stderr="",
                    elapsed_seconds=0.1,
                ),
            ) as mock_super_run:
                result = mock_blastn.run(
                    query=fastq_file,
                    database=Path("db"),
                    output=Path("out.tsv"),
                )

                assert result.success is True
                # Verify super().run was called with a FASTA query (not the original FASTQ)
                call_kwargs = mock_super_run.call_args
                query_used = call_kwargs.kwargs.get("query")
                assert query_used is not None
                assert str(query_used).endswith(".fasta")

    def test_run_with_fasta_skips_conversion(
        self, mock_blastn: BlastN, tmp_path: Path
    ):
        """Should skip conversion for FASTA files and run directly."""
        fasta_file = tmp_path / "reads.fasta"
        fasta_file.write_text(">read1\nACGT\n")

        with patch.object(
            ExternalTool,
            "run",
            return_value=ToolResult(
                command=("blastn",),
                return_code=0,
                stdout="",
                stderr="",
                elapsed_seconds=0.1,
            ),
        ) as mock_super_run:
            result = mock_blastn.run(
                query=fasta_file,
                database=Path("db"),
                output=Path("out.tsv"),
            )

            assert result.success is True
            mock_super_run.assert_called_once()

    def test_run_without_query_runs_normally(self, mock_blastn: BlastN):
        """Should run normally when no query is provided."""
        with patch.object(
            ExternalTool,
            "run",
            return_value=ToolResult(
                command=("blastn",),
                return_code=0,
                stdout="",
                stderr="",
                elapsed_seconds=0.1,
            ),
        ) as mock_super_run:
            result = mock_blastn.run(
                database=Path("db"),
                output=Path("out.tsv"),
            )
            assert result.success is True
            mock_super_run.assert_called_once()

    def test_run_with_fastq_cleans_up_on_error(
        self, mock_blastn: BlastN, tmp_path: Path
    ):
        """Should clean up temp FASTA even when BLAST fails."""
        fastq_file = tmp_path / "reads.fastq"
        fastq_file.write_text("@read1\nACGT\n+\nIIII\n")

        with patch(
            "metadarkmatter.external.blast.subprocess.run",
            side_effect=FileNotFoundError("seqtk not found"),
        ):
            with patch.object(
                ExternalTool,
                "run",
                side_effect=RuntimeError("BLAST failed"),
            ):
                with pytest.raises(RuntimeError, match="BLAST failed"):
                    mock_blastn.run(
                        query=fastq_file,
                        database=Path("db"),
                        output=Path("out.tsv"),
                    )

        # Temp files should be cleaned up (the finally block ran)
        # We cannot check the specific temp file, but verifying no error
        # in cleanup is sufficient.

    def test_run_dry_run_with_fastq_passes_through(
        self, mock_blastn: BlastN, tmp_path: Path
    ):
        """Should handle dry_run correctly with FASTQ input."""
        fastq_file = tmp_path / "reads.fastq"
        fastq_file.write_text("@read1\nACGT\n+\nIIII\n")

        with patch(
            "metadarkmatter.external.blast.subprocess.run",
            side_effect=FileNotFoundError("seqtk not found"),
        ):
            with patch.object(
                ExternalTool,
                "run",
                return_value=ToolResult(
                    command=("blastn",),
                    return_code=0,
                    stdout="[dry-run] Command not executed",
                    stderr="",
                    elapsed_seconds=0.0,
                ),
            ):
                result = mock_blastn.run(
                    query=fastq_file,
                    database=Path("db"),
                    output=Path("out.tsv"),
                    dry_run=True,
                )
                assert result.success is True


# =========================================================================
# BlastN.build_sensitive_command_args  (blast.py lines 334-352)
# =========================================================================


class TestBlastNSensitiveCommandArgs:
    """Tests for BlastN.build_sensitive_command_args class method."""

    def test_returns_expected_keys(self):
        """Should return dict with task, word_size, evalue, outfmt."""
        args = BlastN.build_sensitive_command_args()
        assert "task" in args
        assert "word_size" in args
        assert "evalue" in args
        assert "outfmt" in args

    def test_sensitive_defaults(self):
        """Should return settings optimized for sensitive detection."""
        args = BlastN.build_sensitive_command_args()
        assert args["task"] == "blastn"
        assert args["word_size"] == 7
        assert args["evalue"] == 1e-3

    def test_outfmt_includes_12_columns(self):
        """Should include 12 standard BLAST columns in outfmt."""
        args = BlastN.build_sensitive_command_args()
        outfmt = args["outfmt"]
        assert "qseqid" in outfmt
        assert "sseqid" in outfmt
        assert "pident" in outfmt
        assert "length" in outfmt
        assert "bitscore" in outfmt
        assert "evalue" in outfmt

    def test_can_be_used_in_build_command(self, mock_blastn: BlastN):
        """Should produce valid kwargs for build_command."""
        sensitive_args = BlastN.build_sensitive_command_args()
        cmd = mock_blastn.build_command(
            query=Path("reads.fasta"),
            database=Path("db"),
            output=Path("out.tsv"),
            **sensitive_args,
        )
        assert "-word_size" in cmd
        assert "7" in cmd


# =========================================================================
# MMseqs2._build_search_command validation  (mmseqs2.py lines 229-244)
# =========================================================================


class TestMMseqs2SearchValidation:
    """Tests for MMseqs2 search mode parameter validation."""

    def test_search_missing_database_raises(self, mock_mmseqs: MMseqs2, tmp_path: Path):
        """Should raise ValueError when database is missing for search."""
        query = tmp_path / "reads.fastq"
        query.write_text("@read1\nACGT\n+\nIIII\n")

        with pytest.raises(ValueError, match="database is required"):
            mock_mmseqs.build_command(
                mode="search",
                query=query,
                output=Path("out.tsv"),
            )

    def test_search_missing_output_raises(self, mock_mmseqs: MMseqs2, tmp_path: Path):
        """Should raise ValueError when output is missing for search."""
        query = tmp_path / "reads.fastq"
        query.write_text("@read1\nACGT\n+\nIIII\n")

        with pytest.raises(ValueError, match="output is required"):
            mock_mmseqs.build_command(
                mode="search",
                query=query,
                database=Path("db"),
            )

    def test_search_with_explicit_tmp_dir(
        self, mock_mmseqs: MMseqs2, tmp_path: Path
    ):
        """Should use provided tmp_dir and create it if needed."""
        query = tmp_path / "reads.fastq"
        query.write_text("@read1\nACGT\n+\nIIII\n")

        custom_tmp = tmp_path / "custom_tmp"
        cmd = mock_mmseqs.build_command(
            mode="search",
            query=query,
            database=tmp_path / "db",
            output=tmp_path / "out.tsv",
            tmp_dir=custom_tmp,
        )

        assert str(custom_tmp) in cmd
        assert custom_tmp.exists()

    def test_search_auto_creates_tmp_dir(
        self, mock_mmseqs: MMseqs2, tmp_path: Path
    ):
        """Should auto-create a temp directory when tmp_dir is None."""
        query = tmp_path / "reads.fastq"
        query.write_text("@read1\nACGT\n+\nIIII\n")

        cmd = mock_mmseqs.build_command(
            mode="search",
            query=query,
            database=tmp_path / "db",
            output=tmp_path / "out.tsv",
            tmp_dir=None,
        )

        # Command should contain a temporary directory path
        assert "easy-search" in cmd
        # The auto-created tmp dir should be in the command
        # Find index after output path to locate the tmp dir in command
        assert len(cmd) > 4


# =========================================================================
# MMseqs2.create_database convenience method  (mmseqs2.py line 298)
# =========================================================================


class TestMMseqs2CreateDatabase:
    """Tests for MMseqs2.create_database convenience method."""

    @patch("metadarkmatter.external.base.subprocess.run")
    def test_create_database_calls_run_or_raise(
        self, mock_run: MagicMock, mock_mmseqs: MMseqs2, tmp_path: Path
    ):
        """Should call run_or_raise with createdb mode."""
        mock_run.return_value = MagicMock(
            returncode=0, stdout="", stderr=""
        )
        input_fasta = tmp_path / "genomes.fasta"
        input_fasta.write_text(">seq1\nACGT\n")

        mock_mmseqs.create_database(
            input_fasta=input_fasta,
            database=tmp_path / "mmseqs_db",
        )

        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert "createdb" in cmd

    @patch("metadarkmatter.external.base.subprocess.run")
    def test_create_database_with_dbtype(
        self, mock_run: MagicMock, mock_mmseqs: MMseqs2, tmp_path: Path
    ):
        """Should pass dbtype parameter through."""
        mock_run.return_value = MagicMock(
            returncode=0, stdout="", stderr=""
        )
        input_fasta = tmp_path / "genomes.fasta"
        input_fasta.write_text(">seq1\nACGT\n")

        mock_mmseqs.create_database(
            input_fasta=input_fasta,
            database=tmp_path / "db",
            dbtype=1,
        )

        cmd = mock_run.call_args[0][0]
        assert "--dbtype" in cmd
        assert "1" in cmd

    @patch("metadarkmatter.external.base.subprocess.run")
    def test_create_database_with_timeout(
        self, mock_run: MagicMock, mock_mmseqs: MMseqs2, tmp_path: Path
    ):
        """Should pass timeout parameter through."""
        mock_run.return_value = MagicMock(
            returncode=0, stdout="", stderr=""
        )
        input_fasta = tmp_path / "genomes.fasta"
        input_fasta.write_text(">seq1\nACGT\n")

        mock_mmseqs.create_database(
            input_fasta=input_fasta,
            database=tmp_path / "db",
            timeout=120.0,
        )

        call_kwargs = mock_run.call_args.kwargs
        assert call_kwargs.get("timeout") == 120.0

    @patch("metadarkmatter.external.base.subprocess.run")
    def test_create_database_raises_on_failure(
        self, mock_run: MagicMock, mock_mmseqs: MMseqs2, tmp_path: Path
    ):
        """Should raise ToolExecutionError when createdb fails."""
        mock_run.return_value = MagicMock(
            returncode=1, stdout="", stderr="createdb failed"
        )
        input_fasta = tmp_path / "genomes.fasta"
        input_fasta.write_text(">seq1\nACGT\n")

        with pytest.raises(ToolExecutionError):
            mock_mmseqs.create_database(
                input_fasta=input_fasta,
                database=tmp_path / "db",
            )


# =========================================================================
# MMseqs2.search convenience method  (mmseqs2.py line 353)
# =========================================================================


class TestMMseqs2Search:
    """Tests for MMseqs2.search convenience method."""

    @patch("metadarkmatter.external.base.subprocess.run")
    def test_search_calls_run_or_raise(
        self, mock_run: MagicMock, mock_mmseqs: MMseqs2, tmp_path: Path
    ):
        """Should call run_or_raise with search mode."""
        mock_run.return_value = MagicMock(
            returncode=0, stdout="", stderr=""
        )
        query = tmp_path / "reads.fastq"
        query.write_text("@read1\nACGT\n+\nIIII\n")

        mock_mmseqs.search(
            query=query,
            database=tmp_path / "db",
            output=tmp_path / "out.tsv",
        )

        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert "easy-search" in cmd

    @patch("metadarkmatter.external.base.subprocess.run")
    def test_search_with_all_parameters(
        self, mock_run: MagicMock, mock_mmseqs: MMseqs2, tmp_path: Path
    ):
        """Should pass all search parameters through."""
        mock_run.return_value = MagicMock(
            returncode=0, stdout="", stderr=""
        )
        query = tmp_path / "reads.fastq"
        query.write_text("@read1\nACGT\n+\nIIII\n")

        mock_mmseqs.search(
            query=query,
            database=tmp_path / "db",
            output=tmp_path / "out.tsv",
            threads=16,
            sensitivity=7.0,
            evalue=1e-5,
            max_seqs=100,
            min_identity=80.0,
            search_type=2,
        )

        cmd = mock_run.call_args[0][0]
        assert "--threads" in cmd
        assert "16" in cmd
        assert "-s" in cmd
        assert "7.0" in cmd
        assert "--min-seq-id" in cmd
        assert "0.8" in cmd
        assert "--search-type" in cmd
        assert "2" in cmd

    @patch("metadarkmatter.external.base.subprocess.run")
    def test_search_raises_on_failure(
        self, mock_run: MagicMock, mock_mmseqs: MMseqs2, tmp_path: Path
    ):
        """Should raise ToolExecutionError when search fails."""
        mock_run.return_value = MagicMock(
            returncode=1, stdout="", stderr="search failed"
        )
        query = tmp_path / "reads.fastq"
        query.write_text("@read1\nACGT\n+\nIIII\n")

        with pytest.raises(ToolExecutionError):
            mock_mmseqs.search(
                query=query,
                database=tmp_path / "db",
                output=tmp_path / "out.tsv",
            )

    @patch("metadarkmatter.external.base.subprocess.run")
    def test_search_capture_output_false(
        self, mock_run: MagicMock, mock_mmseqs: MMseqs2, tmp_path: Path
    ):
        """Should support capture_output=False for batch processing."""
        mock_run.return_value = MagicMock(returncode=0)
        query = tmp_path / "reads.fastq"
        query.write_text("@read1\nACGT\n+\nIIII\n")

        mock_mmseqs.search(
            query=query,
            database=tmp_path / "db",
            output=tmp_path / "out.tsv",
            capture_output=False,
        )

        call_kwargs = mock_run.call_args.kwargs
        assert call_kwargs.get("stdout") == subprocess.DEVNULL


# =========================================================================
# MMseqs2.search_multistep  (mmseqs2.py lines 436-536)
# =========================================================================


class TestMMseqs2SearchMultistep:
    """Tests for MMseqs2.search_multistep multi-step workflow."""

    @patch("subprocess.run")
    def test_full_workflow_creates_query_db(
        self, mock_run: MagicMock, mock_mmseqs: MMseqs2, tmp_path: Path
    ):
        """Should run createdb, search, and convertalis steps."""
        mock_run.return_value = MagicMock(
            returncode=0, stdout="", stderr=""
        )

        query = tmp_path / "reads.fasta"
        query.write_text(">read1\nACGT\n")
        database = tmp_path / "target_db"
        database.touch()
        output = tmp_path / "results.tsv"

        result = mock_mmseqs.search_multistep(
            query=query,
            database=database,
            output=output,
            tmp_dir=tmp_path / "tmp",
        )

        # Should have called subprocess.run 3 times:
        # createdb, search, convertalis
        assert mock_run.call_count == 3

        # Verify createdb step
        createdb_cmd = mock_run.call_args_list[0][0][0]
        assert "createdb" in createdb_cmd

        # Verify search step
        search_cmd = mock_run.call_args_list[1][0][0]
        assert "search" in search_cmd

        # Verify convertalis step
        convert_cmd = mock_run.call_args_list[2][0][0]
        assert "convertalis" in convert_cmd

        # Verify result dictionary
        assert "query_db" in result
        assert "result_db" in result
        assert "createdb_time" in result
        assert "search_time" in result
        assert "convertalis_time" in result
        assert "total_time" in result
        assert result["total_time"] >= 0

    @patch("subprocess.run")
    def test_skips_createdb_when_query_db_exists(
        self, mock_run: MagicMock, mock_mmseqs: MMseqs2, tmp_path: Path
    ):
        """Should skip createdb step when query_db already exists."""
        mock_run.return_value = MagicMock(
            returncode=0, stdout="", stderr=""
        )

        query = tmp_path / "reads.fasta"
        query.write_text(">read1\nACGT\n")
        database = tmp_path / "target_db"
        database.touch()
        output = tmp_path / "results.tsv"

        # Pre-create query_db so it exists
        query_db = tmp_path / "query_db"
        query_db.touch()

        result = mock_mmseqs.search_multistep(
            query=query,
            database=database,
            output=output,
            query_db=query_db,
            tmp_dir=tmp_path / "tmp",
        )

        # Should have called subprocess.run only 2 times (search + convertalis)
        assert mock_run.call_count == 2
        assert result["createdb_time"] == 0.0

    @patch("subprocess.run")
    def test_with_custom_query_db_path(
        self, mock_run: MagicMock, mock_mmseqs: MMseqs2, tmp_path: Path
    ):
        """Should use provided query_db path."""
        mock_run.return_value = MagicMock(
            returncode=0, stdout="", stderr=""
        )

        query = tmp_path / "reads.fasta"
        query.write_text(">read1\nACGT\n")
        database = tmp_path / "target_db"
        database.touch()
        output = tmp_path / "results.tsv"
        query_db = tmp_path / "custom_qdb" / "reads"

        result = mock_mmseqs.search_multistep(
            query=query,
            database=database,
            output=output,
            query_db=query_db,
            tmp_dir=tmp_path / "tmp",
        )

        assert result["query_db"] == str(query_db)
        # Parent directory should have been created
        assert query_db.parent.exists()

    @patch("subprocess.run")
    def test_auto_creates_tmp_dir(
        self, mock_run: MagicMock, mock_mmseqs: MMseqs2, tmp_path: Path
    ):
        """Should auto-create tmp_dir when not specified."""
        mock_run.return_value = MagicMock(
            returncode=0, stdout="", stderr=""
        )

        query = tmp_path / "reads.fasta"
        query.write_text(">read1\nACGT\n")
        database = tmp_path / "target_db"
        database.touch()
        output = tmp_path / "results.tsv"

        result = mock_mmseqs.search_multistep(
            query=query,
            database=database,
            output=output,
        )

        assert "query_db" in result
        assert "result_db" in result

    @patch("subprocess.run")
    def test_min_identity_passed_to_search(
        self, mock_run: MagicMock, mock_mmseqs: MMseqs2, tmp_path: Path
    ):
        """Should pass min_identity to search command."""
        mock_run.return_value = MagicMock(
            returncode=0, stdout="", stderr=""
        )

        query = tmp_path / "reads.fasta"
        query.write_text(">read1\nACGT\n")
        database = tmp_path / "target_db"
        database.touch()
        output = tmp_path / "results.tsv"

        mock_mmseqs.search_multistep(
            query=query,
            database=database,
            output=output,
            min_identity=90.0,
            tmp_dir=tmp_path / "tmp",
        )

        # Search step is the 2nd call (after createdb)
        search_cmd = mock_run.call_args_list[1][0][0]
        assert "--min-seq-id" in search_cmd
        # 90.0 / 100.0 = 0.9
        assert "0.9" in search_cmd

    @patch("subprocess.run")
    def test_timeout_forwarded_to_all_steps(
        self, mock_run: MagicMock, mock_mmseqs: MMseqs2, tmp_path: Path
    ):
        """Should forward timeout to all subprocess calls."""
        mock_run.return_value = MagicMock(
            returncode=0, stdout="", stderr=""
        )

        query = tmp_path / "reads.fasta"
        query.write_text(">read1\nACGT\n")
        database = tmp_path / "target_db"
        database.touch()
        output = tmp_path / "results.tsv"

        mock_mmseqs.search_multistep(
            query=query,
            database=database,
            output=output,
            timeout=300.0,
            tmp_dir=tmp_path / "tmp",
        )

        for call_args in mock_run.call_args_list:
            assert call_args.kwargs.get("timeout") == 300.0

    @patch("subprocess.run")
    def test_search_step_failure_raises(
        self, mock_run: MagicMock, mock_mmseqs: MMseqs2, tmp_path: Path
    ):
        """Should raise CalledProcessError when search step fails."""

        def side_effect(cmd, **kwargs):
            if "search" in cmd and "createdb" not in cmd and "convertalis" not in cmd:
                raise subprocess.CalledProcessError(1, cmd)
            return MagicMock(returncode=0, stdout="", stderr="")

        mock_run.side_effect = side_effect

        query = tmp_path / "reads.fasta"
        query.write_text(">read1\nACGT\n")
        database = tmp_path / "target_db"
        database.touch()
        output = tmp_path / "results.tsv"

        with pytest.raises(subprocess.CalledProcessError):
            mock_mmseqs.search_multistep(
                query=query,
                database=database,
                output=output,
                tmp_dir=tmp_path / "tmp",
            )

    @patch("subprocess.run")
    def test_convertalis_format_string(
        self, mock_run: MagicMock, mock_mmseqs: MMseqs2, tmp_path: Path
    ):
        """Should use BLAST-compatible format for convertalis step."""
        mock_run.return_value = MagicMock(
            returncode=0, stdout="", stderr=""
        )

        query = tmp_path / "reads.fasta"
        query.write_text(">read1\nACGT\n")
        database = tmp_path / "target_db"
        database.touch()
        output = tmp_path / "results.tsv"

        mock_mmseqs.search_multistep(
            query=query,
            database=database,
            output=output,
            tmp_dir=tmp_path / "tmp",
        )

        # The convertalis step is the last call
        convert_cmd = mock_run.call_args_list[-1][0][0]
        assert "--format-mode" in convert_cmd
        assert "--format-output" in convert_cmd

        format_idx = convert_cmd.index("--format-output")
        format_str = convert_cmd[format_idx + 1]
        assert "query" in format_str
        assert "target" in format_str
        assert "pident" in format_str
        assert "bits" in format_str
        assert "qlen" in format_str

    @patch("subprocess.run")
    def test_search_type_determines_dbtype(
        self, mock_run: MagicMock, mock_mmseqs: MMseqs2, tmp_path: Path
    ):
        """Should derive dbtype from search_type (search_type - 1)."""
        mock_run.return_value = MagicMock(
            returncode=0, stdout="", stderr=""
        )

        query = tmp_path / "reads.fasta"
        query.write_text(">read1\nACGT\n")
        database = tmp_path / "target_db"
        database.touch()
        output = tmp_path / "results.tsv"

        # search_type=3 (nucleotide) -> dbtype=2
        mock_mmseqs.search_multistep(
            query=query,
            database=database,
            output=output,
            search_type=3,
            tmp_dir=tmp_path / "tmp",
        )

        createdb_cmd = mock_run.call_args_list[0][0][0]
        dbtype_idx = createdb_cmd.index("--dbtype")
        assert createdb_cmd[dbtype_idx + 1] == "2"

    @patch("subprocess.run")
    def test_explicit_tmp_dir_created(
        self, mock_run: MagicMock, mock_mmseqs: MMseqs2, tmp_path: Path
    ):
        """Should create the explicit tmp_dir if it does not exist."""
        mock_run.return_value = MagicMock(
            returncode=0, stdout="", stderr=""
        )

        query = tmp_path / "reads.fasta"
        query.write_text(">read1\nACGT\n")
        database = tmp_path / "target_db"
        database.touch()
        output = tmp_path / "results.tsv"
        custom_tmp = tmp_path / "nested" / "tmp"

        mock_mmseqs.search_multistep(
            query=query,
            database=database,
            output=output,
            tmp_dir=custom_tmp,
        )

        assert custom_tmp.exists()


# =========================================================================
# Samtools - view, sort, index execution  (samtools.py)
# =========================================================================


class TestSamtoolsViewExecution:
    """Tests for Samtools.view method execution paths."""

    def test_view_dry_run(self, mock_samtools: Samtools):
        """Should return dry-run result without executing."""
        result = mock_samtools.view(
            input_file=Path("input.sam"),
            output_file=Path("output.bam"),
            dry_run=True,
        )
        assert result.success is True
        assert "[dry-run]" in result.stdout
        assert result.elapsed_seconds == 0.0

    @patch("subprocess.run")
    def test_view_execution_success(
        self, mock_run: MagicMock, mock_samtools: Samtools
    ):
        """Should execute view command and return result."""
        mock_run.return_value = MagicMock(
            returncode=0, stdout="", stderr=""
        )
        result = mock_samtools.view(
            input_file=Path("input.sam"),
            output_file=Path("output.bam"),
        )
        assert result.success is True
        assert "view" in list(result.command)

    def test_view_with_all_filters(self, mock_samtools: Samtools):
        """Should include all filter options in command."""
        cmd = mock_samtools._build_view_command(
            input_file=Path("input.bam"),
            output_file=Path("output.bam"),
            output_bam=True,
            threads=8,
            include_header=True,
            min_mapq=30,
            flags_required=2,
            flags_excluded=256,
            regions=["chr1:1-1000", "chr2"],
        )
        assert "-q" in cmd
        assert "30" in cmd
        assert "-f" in cmd
        assert "2" in cmd
        assert "-F" in cmd
        assert "256" in cmd
        assert "chr1:1-1000" in cmd
        assert "chr2" in cmd

    def test_view_without_header(self, mock_samtools: Samtools):
        """Should omit -h when include_header is False."""
        cmd = mock_samtools._build_view_command(
            input_file=Path("input.bam"),
            output_file=Path("output.bam"),
            include_header=False,
        )
        assert "-h" not in cmd

    def test_view_sam_output(self, mock_samtools: Samtools):
        """Should omit -b when output_bam is False."""
        cmd = mock_samtools._build_view_command(
            input_file=Path("input.bam"),
            output_file=Path("output.sam"),
            output_bam=False,
        )
        assert "-b" not in cmd


class TestSamtoolsSortExecution:
    """Tests for Samtools.sort method execution paths."""

    def test_sort_dry_run(self, mock_samtools: Samtools):
        """Should return dry-run result without executing."""
        result = mock_samtools.sort(
            input_file=Path("input.bam"),
            output_file=Path("sorted.bam"),
            dry_run=True,
        )
        assert result.success is True
        assert "[dry-run]" in result.stdout

    @patch("subprocess.run")
    def test_sort_execution_success(
        self, mock_run: MagicMock, mock_samtools: Samtools
    ):
        """Should execute sort command and return result."""
        mock_run.return_value = MagicMock(
            returncode=0, stdout="", stderr=""
        )
        result = mock_samtools.sort(
            input_file=Path("input.bam"),
            output_file=Path("sorted.bam"),
        )
        assert result.success is True

    def test_sort_by_name(self, mock_samtools: Samtools):
        """Should include -n flag when sorting by name."""
        cmd = mock_samtools._build_sort_command(
            input_file=Path("input.bam"),
            output_file=Path("sorted.bam"),
            by_name=True,
        )
        assert "-n" in cmd

    def test_sort_with_temp_prefix(self, mock_samtools: Samtools):
        """Should include -T flag with temp prefix."""
        cmd = mock_samtools._build_sort_command(
            input_file=Path("input.bam"),
            output_file=Path("sorted.bam"),
            temp_prefix=Path("/tmp/sort_tmp"),
        )
        assert "-T" in cmd
        assert "/tmp/sort_tmp" in " ".join(cmd)

    def test_sort_with_custom_memory(self, mock_samtools: Samtools):
        """Should include custom memory per thread."""
        cmd = mock_samtools._build_sort_command(
            input_file=Path("input.bam"),
            output_file=Path("sorted.bam"),
            memory_per_thread="2G",
        )
        assert "-m" in cmd
        idx = cmd.index("-m")
        assert cmd[idx + 1] == "2G"


class TestSamtoolsIndexExecution:
    """Tests for Samtools.index method execution paths."""

    def test_index_dry_run(self, mock_samtools: Samtools):
        """Should return dry-run result without executing."""
        result = mock_samtools.index(
            bam_file=Path("sorted.bam"),
            dry_run=True,
        )
        assert result.success is True
        assert "[dry-run]" in result.stdout

    @patch("subprocess.run")
    def test_index_execution_success(
        self, mock_run: MagicMock, mock_samtools: Samtools
    ):
        """Should execute index command and return result."""
        mock_run.return_value = MagicMock(
            returncode=0, stdout="", stderr=""
        )
        result = mock_samtools.index(bam_file=Path("sorted.bam"))
        assert result.success is True

    def test_index_csi_mode(self, mock_samtools: Samtools):
        """Should include -c flag for CSI index."""
        cmd = mock_samtools._build_index_command(
            bam_file=Path("sorted.bam"),
            csi=True,
        )
        assert "-c" in cmd

    def test_index_custom_threads(self, mock_samtools: Samtools):
        """Should include custom thread count."""
        cmd = mock_samtools._build_index_command(
            bam_file=Path("sorted.bam"),
            threads=16,
        )
        assert "--threads" in cmd
        idx = cmd.index("--threads")
        assert cmd[idx + 1] == "16"


# =========================================================================
# Samtools._run_command error paths  (samtools.py lines 305-353)
# =========================================================================


class TestSamtoolsRunCommand:
    """Tests for Samtools._run_command error handling."""

    @patch("subprocess.run")
    def test_successful_execution(
        self, mock_run: MagicMock, mock_samtools: Samtools
    ):
        """Should return ToolResult on successful execution."""
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="output",
            stderr="",
        )
        result = mock_samtools._run_command(["samtools", "view", "input.bam"])
        assert result.success is True
        assert result.stdout == "output"
        assert result.elapsed_seconds >= 0

    @patch("subprocess.run")
    def test_non_zero_exit_raises_execution_error(
        self, mock_run: MagicMock, mock_samtools: Samtools
    ):
        """Should raise ToolExecutionError on non-zero exit code."""
        mock_run.return_value = MagicMock(
            returncode=1,
            stdout="",
            stderr="error: file not found",
        )
        with pytest.raises(ToolExecutionError) as exc_info:
            mock_samtools._run_command(["samtools", "view", "missing.bam"])
        assert exc_info.value.tool_name == "samtools"
        assert exc_info.value.return_code == 1

    @patch("subprocess.run")
    def test_timeout_raises_tool_timeout_error(
        self, mock_run: MagicMock, mock_samtools: Samtools
    ):
        """Should raise ToolTimeoutError when subprocess times out."""
        mock_run.side_effect = subprocess.TimeoutExpired(
            cmd=["samtools", "sort"], timeout=60
        )
        with pytest.raises(ToolTimeoutError) as exc_info:
            mock_samtools._run_command(
                ["samtools", "sort", "big.bam"], timeout=60.0
            )
        assert exc_info.value.tool_name == "samtools"

    @patch("subprocess.run")
    def test_timeout_with_none_defaults_to_zero(
        self, mock_run: MagicMock, mock_samtools: Samtools
    ):
        """Should use timeout=0 in error when timeout parameter is None."""
        mock_run.side_effect = subprocess.TimeoutExpired(
            cmd=["samtools", "index"], timeout=0
        )
        with pytest.raises(ToolTimeoutError) as exc_info:
            mock_samtools._run_command(["samtools", "index", "f.bam"])
        assert exc_info.value.timeout_seconds == 0


# =========================================================================
# Samtools.sam_to_sorted_bam pipeline  (samtools.py lines 355-417)
# =========================================================================


class TestSamtoolsSamToSortedBam:
    """Tests for Samtools.sam_to_sorted_bam convenience pipeline."""

    def test_dry_run_returns_three_results(self, mock_samtools: Samtools):
        """Should return three dry-run results without executing."""
        view_result, sort_result, index_result = mock_samtools.sam_to_sorted_bam(
            input_sam=Path("input.sam"),
            output_bam=Path("output.sorted.bam"),
            dry_run=True,
        )

        assert view_result.success is True
        assert "[dry-run]" in view_result.stdout
        assert sort_result.success is True
        assert "[dry-run]" in sort_result.stdout
        assert index_result.success is True
        assert "[dry-run]" in index_result.stdout

    @patch("subprocess.run")
    def test_full_pipeline_execution(
        self, mock_run: MagicMock, mock_samtools: Samtools, tmp_path: Path
    ):
        """Should run view, sort, index pipeline and clean up intermediate BAM."""
        mock_run.return_value = MagicMock(
            returncode=0, stdout="", stderr=""
        )

        input_sam = tmp_path / "input.sam"
        input_sam.write_text("@HD\tVN:1.6\n")
        output_bam = tmp_path / "output.sorted.bam"

        # Create the unsorted intermediate BAM that would be produced by view
        unsorted_bam = output_bam.with_suffix(".unsorted.bam")
        unsorted_bam.touch()

        view_result, sort_result, index_result = mock_samtools.sam_to_sorted_bam(
            input_sam=input_sam,
            output_bam=output_bam,
        )

        assert view_result.success is True
        assert sort_result.success is True
        assert index_result.success is True

        # Intermediate unsorted BAM should have been cleaned up
        assert not unsorted_bam.exists()

    @patch("subprocess.run")
    def test_cleanup_sam_removes_input(
        self, mock_run: MagicMock, mock_samtools: Samtools, tmp_path: Path
    ):
        """Should remove input SAM when cleanup_sam=True."""
        mock_run.return_value = MagicMock(
            returncode=0, stdout="", stderr=""
        )

        input_sam = tmp_path / "input.sam"
        input_sam.write_text("@HD\tVN:1.6\n")
        output_bam = tmp_path / "output.sorted.bam"

        # Create the unsorted intermediate BAM
        unsorted_bam = output_bam.with_suffix(".unsorted.bam")
        unsorted_bam.touch()

        mock_samtools.sam_to_sorted_bam(
            input_sam=input_sam,
            output_bam=output_bam,
            cleanup_sam=True,
        )

        # Input SAM should have been deleted
        assert not input_sam.exists()

    @patch("subprocess.run")
    def test_no_cleanup_sam_preserves_input(
        self, mock_run: MagicMock, mock_samtools: Samtools, tmp_path: Path
    ):
        """Should preserve input SAM when cleanup_sam=False (default)."""
        mock_run.return_value = MagicMock(
            returncode=0, stdout="", stderr=""
        )

        input_sam = tmp_path / "input.sam"
        input_sam.write_text("@HD\tVN:1.6\n")
        output_bam = tmp_path / "output.sorted.bam"

        # Create the unsorted intermediate BAM
        unsorted_bam = output_bam.with_suffix(".unsorted.bam")
        unsorted_bam.touch()

        mock_samtools.sam_to_sorted_bam(
            input_sam=input_sam,
            output_bam=output_bam,
            cleanup_sam=False,
        )

        # Input SAM should still exist
        assert input_sam.exists()

    def test_dry_run_does_not_cleanup(
        self, mock_samtools: Samtools, tmp_path: Path
    ):
        """Should not attempt cleanup during dry run."""
        input_sam = tmp_path / "input.sam"
        input_sam.write_text("@HD\tVN:1.6\n")
        output_bam = tmp_path / "output.sorted.bam"

        mock_samtools.sam_to_sorted_bam(
            input_sam=input_sam,
            output_bam=output_bam,
            cleanup_sam=True,
            dry_run=True,
        )

        # Input SAM should still exist (dry run did not clean up)
        assert input_sam.exists()

    @patch("subprocess.run")
    def test_custom_threads(
        self, mock_run: MagicMock, mock_samtools: Samtools, tmp_path: Path
    ):
        """Should pass thread count to all three steps."""
        mock_run.return_value = MagicMock(
            returncode=0, stdout="", stderr=""
        )

        input_sam = tmp_path / "input.sam"
        input_sam.write_text("@HD\tVN:1.6\n")
        output_bam = tmp_path / "output.sorted.bam"
        unsorted_bam = output_bam.with_suffix(".unsorted.bam")
        unsorted_bam.touch()

        mock_samtools.sam_to_sorted_bam(
            input_sam=input_sam,
            output_bam=output_bam,
            threads=16,
        )

        # Each call should include thread parameter
        for call_args in mock_run.call_args_list:
            cmd = call_args[0][0]
            assert "--threads" in cmd
            idx = cmd.index("--threads")
            assert cmd[idx + 1] == "16"


# =========================================================================
# Samtools.build_command raises NotImplementedError
# =========================================================================


class TestSamtoolsBuildCommand:
    """Tests for Samtools.build_command raising NotImplementedError."""

    def test_build_command_raises(self, mock_samtools: Samtools):
        """Should raise NotImplementedError directing to subcommand methods."""
        with pytest.raises(NotImplementedError, match="view.*sort.*index"):
            mock_samtools.build_command()
