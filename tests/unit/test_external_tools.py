"""Unit tests for external tool wrappers.

Tests the command building logic without requiring the actual tools
to be installed.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metadarkmatter.external.base import (
    ExternalTool,
    ToolExecutionError,
    ToolNotFoundError,
    ToolResult,
    ToolTimeoutError,
)


class TestToolResult:
    """Tests for ToolResult dataclass."""

    def test_create_successful_result(self):
        """Should create result with success=True for exit code 0."""
        result = ToolResult(
            command=("echo", "hello"),
            return_code=0,
            stdout="hello\n",
            stderr="",
            elapsed_seconds=0.1,
        )
        assert result.success is True
        assert result.return_code == 0
        assert result.command_string == "echo hello"

    def test_create_failed_result(self):
        """Should create result with success=False for non-zero exit code."""
        result = ToolResult(
            command=("false",),
            return_code=1,
            stdout="",
            stderr="error message",
            elapsed_seconds=0.05,
        )
        assert result.success is False
        assert result.return_code == 1

    def test_result_is_frozen(self):
        """Should be immutable."""
        result = ToolResult(
            command=("test",),
            return_code=0,
            stdout="",
            stderr="",
            elapsed_seconds=0.0,
        )
        with pytest.raises(AttributeError):
            result.return_code = 1


class TestToolNotFoundError:
    """Tests for ToolNotFoundError exception."""

    def test_error_message_basic(self):
        """Should include tool name in message."""
        error = ToolNotFoundError("kraken2")
        assert "kraken2" in str(error)
        assert "not found" in str(error).lower()

    def test_error_with_install_hint(self):
        """Should include install hint in suggestion."""
        error = ToolNotFoundError("kraken2", "conda install -c bioconda kraken2")
        assert "conda install" in str(error)


class TestToolExecutionError:
    """Tests for ToolExecutionError exception."""

    def test_error_includes_command(self):
        """Should include command in message."""
        error = ToolExecutionError(
            tool_name="blastn",
            command=["blastn", "-query", "input.fa", "-db", "database"],
            return_code=1,
            stderr="BLAST error: database not found",
        )
        assert "blastn" in str(error)
        assert "database not found" in str(error)
        assert error.return_code == 1

    def test_error_truncates_long_stderr(self):
        """Should truncate very long error output."""
        long_stderr = "x" * 1000
        error = ToolExecutionError(
            tool_name="test",
            command=["test"],
            return_code=1,
            stderr=long_stderr,
        )
        assert len(str(error)) < len(long_stderr)
        assert "truncated" in str(error).lower()


class TestToolTimeoutError:
    """Tests for ToolTimeoutError exception."""

    def test_error_includes_timeout(self):
        """Should include timeout value in message."""
        error = ToolTimeoutError(
            tool_name="slowtool",
            timeout_seconds=60.0,
            command=["slowtool", "arg"],
        )
        assert "60" in str(error)
        assert "timeout" in str(error).lower()


class TestKraken2CommandBuilding:
    """Tests for Kraken2 command building."""

    def test_basic_command(self):
        """Should build basic Kraken2 command."""
        from metadarkmatter.external.kraken import Kraken2

        kraken = Kraken2()

        # Mock get_executable to avoid needing kraken2 installed
        Kraken2._executable_cache["kraken2"] = Path("/usr/bin/kraken2")

        cmd = kraken.build_command(
            reads_1=Path("reads_R1.fastq"),
            database=Path("/db/kraken"),
            output=Path("output.kraken"),
            report=Path("output.kreport"),
        )

        assert "/usr/bin/kraken2" in cmd
        assert "--db" in cmd
        assert "/db/kraken" in " ".join(cmd)
        assert "--output" in cmd
        assert "--report" in cmd

    def test_paired_end_command(self):
        """Should include --paired flag for paired-end reads."""
        from metadarkmatter.external.kraken import Kraken2

        kraken = Kraken2()
        Kraken2._executable_cache["kraken2"] = Path("/usr/bin/kraken2")

        cmd = kraken.build_command(
            reads_1=Path("reads_R1.fastq"),
            reads_2=Path("reads_R2.fastq"),
            database=Path("/db/kraken"),
            output=Path("output.kraken"),
            report=Path("output.kreport"),
        )

        assert "--paired" in cmd
        assert "reads_R1.fastq" in " ".join(cmd)
        assert "reads_R2.fastq" in " ".join(cmd)

    def test_confidence_threshold(self):
        """Should include confidence when specified."""
        from metadarkmatter.external.kraken import Kraken2

        kraken = Kraken2()
        Kraken2._executable_cache["kraken2"] = Path("/usr/bin/kraken2")

        cmd = kraken.build_command(
            reads_1=Path("reads.fastq"),
            database=Path("/db"),
            output=Path("out.kraken"),
            report=Path("out.kreport"),
            confidence=0.5,
        )

        assert "--confidence" in cmd
        assert "0.5" in cmd


class TestExtractKrakenReadsCommandBuilding:
    """Tests for ExtractKrakenReads command building."""

    def test_basic_command(self):
        """Should build basic extraction command."""
        from metadarkmatter.external.kraken import ExtractKrakenReads

        extractor = ExtractKrakenReads()
        ExtractKrakenReads._executable_cache["extract_kraken_reads.py"] = Path(
            "/usr/bin/extract_kraken_reads.py"
        )

        cmd = extractor.build_command(
            kraken_output=Path("sample.kraken"),
            kraken_report=Path("sample.kreport"),
            reads_1=Path("sample_R1.fastq"),
            output_1=Path("extracted_R1.fastq"),
            taxid=1224,
        )

        assert "-k" in cmd
        assert "-r" in cmd
        assert "-t" in cmd
        assert "1224" in cmd

    def test_include_children_flag(self):
        """Should include --include-children by default."""
        from metadarkmatter.external.kraken import ExtractKrakenReads

        extractor = ExtractKrakenReads()
        ExtractKrakenReads._executable_cache["extract_kraken_reads.py"] = Path(
            "/usr/bin/extract_kraken_reads.py"
        )

        cmd = extractor.build_command(
            kraken_output=Path("sample.kraken"),
            kraken_report=Path("sample.kreport"),
            reads_1=Path("sample_R1.fastq"),
            output_1=Path("extracted_R1.fastq"),
            taxid=1224,
            include_children=True,
        )

        assert "--include-children" in cmd


class TestBowtie2CommandBuilding:
    """Tests for Bowtie2 command building."""

    def test_basic_command(self):
        """Should build basic Bowtie2 alignment command."""
        from metadarkmatter.external.bowtie2 import Bowtie2

        bowtie = Bowtie2()
        Bowtie2._executable_cache["bowtie2"] = Path("/usr/bin/bowtie2")

        cmd = bowtie.build_command(
            index_prefix=Path("index"),
            reads_1=Path("reads_R1.fastq"),
            output_sam=Path("output.sam"),
        )

        assert "-x" in cmd
        assert "-S" in cmd
        assert "index" in " ".join(cmd)

    def test_local_mode(self):
        """Should include --local flag for local alignment."""
        from metadarkmatter.external.bowtie2 import Bowtie2

        bowtie = Bowtie2()
        Bowtie2._executable_cache["bowtie2"] = Path("/usr/bin/bowtie2")

        cmd = bowtie.build_command(
            index_prefix=Path("index"),
            reads_1=Path("reads.fastq"),
            output_sam=Path("output.sam"),
            mode="local",
        )

        assert "--local" in cmd

    def test_max_alignments(self):
        """Should include -k flag for multiple alignments."""
        from metadarkmatter.external.bowtie2 import Bowtie2

        bowtie = Bowtie2()
        Bowtie2._executable_cache["bowtie2"] = Path("/usr/bin/bowtie2")

        cmd = bowtie.build_command(
            index_prefix=Path("index"),
            reads_1=Path("reads.fastq"),
            output_sam=Path("output.sam"),
            max_alignments=100,
        )

        assert "-k" in cmd
        assert "100" in cmd


class TestBowtie2BuildCommandBuilding:
    """Tests for Bowtie2Build command building."""

    def test_basic_command(self):
        """Should build basic index building command."""
        from metadarkmatter.external.bowtie2 import Bowtie2Build

        builder = Bowtie2Build()
        Bowtie2Build._executable_cache["bowtie2-build"] = Path("/usr/bin/bowtie2-build")

        cmd = builder.build_command(
            reference=Path("reference.fasta"),
            index_prefix=Path("index"),
        )

        assert "reference.fasta" in " ".join(cmd)
        assert "index" in " ".join(cmd)
        assert "--threads" in cmd


class TestSamtoolsCommands:
    """Tests for Samtools command building."""

    def test_view_command(self):
        """Should build samtools view command."""
        from metadarkmatter.external.samtools import Samtools

        samtools = Samtools()
        Samtools._executable_cache["samtools"] = Path("/usr/bin/samtools")

        cmd = samtools._build_view_command(
            input_file=Path("input.sam"),
            output_file=Path("output.bam"),
            output_bam=True,
        )

        assert "view" in cmd
        assert "-b" in cmd
        assert "-o" in cmd

    def test_sort_command(self):
        """Should build samtools sort command."""
        from metadarkmatter.external.samtools import Samtools

        samtools = Samtools()
        Samtools._executable_cache["samtools"] = Path("/usr/bin/samtools")

        cmd = samtools._build_sort_command(
            input_file=Path("input.bam"),
            output_file=Path("sorted.bam"),
        )

        assert "sort" in cmd
        assert "-o" in cmd

    def test_index_command(self):
        """Should build samtools index command."""
        from metadarkmatter.external.samtools import Samtools

        samtools = Samtools()
        Samtools._executable_cache["samtools"] = Path("/usr/bin/samtools")

        cmd = samtools._build_index_command(
            bam_file=Path("sorted.bam"),
        )

        assert "index" in cmd
        assert "sorted.bam" in " ".join(cmd)


class TestBlastNCommandBuilding:
    """Tests for BlastN command building."""

    def test_basic_command(self):
        """Should build basic blastn command."""
        from metadarkmatter.external.blast import BlastN

        blastn = BlastN()
        BlastN._executable_cache["blastn"] = Path("/usr/bin/blastn")

        cmd = blastn.build_command(
            query=Path("reads.fasta"),
            database=Path("database"),
            output=Path("results.tsv"),
        )

        assert "-query" in cmd
        assert "-db" in cmd
        assert "-out" in cmd
        assert "-outfmt" in cmd

    def test_sensitivity_parameters(self):
        """Should include sensitivity parameters."""
        from metadarkmatter.external.blast import BlastN

        blastn = BlastN()
        BlastN._executable_cache["blastn"] = Path("/usr/bin/blastn")

        cmd = blastn.build_command(
            query=Path("reads.fasta"),
            database=Path("database"),
            output=Path("results.tsv"),
            word_size=7,
            evalue=1e-3,
        )

        assert "-word_size" in cmd
        assert "7" in cmd
        assert "-evalue" in cmd


class TestConcatenateGenomes:
    """Tests for genome concatenation function."""

    def test_concatenate_creates_file(self, tmp_path):
        """Should create concatenated pangenome file."""
        from metadarkmatter.external.bowtie2 import concatenate_genomes

        # Create test genome files
        genome_dir = tmp_path / "genomes"
        genome_dir.mkdir()

        (genome_dir / "genome1.fna").write_text(">contig1\nACGT\n")
        (genome_dir / "genome2.fna").write_text(">contig2\nTGCA\n")

        output_path = tmp_path / "pangenome.fasta"
        genome_names = concatenate_genomes(
            genome_dir=genome_dir,
            output_path=output_path,
            pattern="*.fna",
        )

        assert output_path.exists()
        assert len(genome_names) == 2
        assert "genome1" in genome_names
        assert "genome2" in genome_names

    def test_concatenate_prefixes_headers(self, tmp_path):
        """Should prefix headers with genome names using standardized pipe format.

        The concatenation function rewrites headers to the standardized format:
        {accession}|{original_contig_id}

        This enables reliable genome identification from multi-contig draft genomes.
        """
        from metadarkmatter.external.bowtie2 import concatenate_genomes

        genome_dir = tmp_path / "genomes"
        genome_dir.mkdir()

        (genome_dir / "mygenome.fna").write_text(">contig1\nACGT\n")

        output_path = tmp_path / "pangenome.fasta"
        concatenate_genomes(
            genome_dir=genome_dir,
            output_path=output_path,
        )

        content = output_path.read_text()
        # Standardized format uses pipe delimiter: {accession}|{contig_id}
        assert ">mygenome|contig1" in content

    def test_concatenate_no_files_raises(self, tmp_path):
        """Should raise error when no files match pattern."""
        from metadarkmatter.external.bowtie2 import concatenate_genomes

        genome_dir = tmp_path / "empty"
        genome_dir.mkdir()

        with pytest.raises(FileNotFoundError):
            concatenate_genomes(
                genome_dir=genome_dir,
                output_path=tmp_path / "pangenome.fasta",
            )


class TestDiamondCommandBuilding:
    """Tests for Diamond command building."""

    def test_makedb_command(self, tmp_path):
        """Should build Diamond makedb command."""
        from metadarkmatter.external.diamond import Diamond

        diamond = Diamond()
        Diamond._executable_cache["diamond"] = Path("/usr/bin/diamond")

        # Create test file for validation
        input_fasta = tmp_path / "proteins.faa"
        input_fasta.write_text(">prot1\nMKTL\n")

        cmd = diamond.build_command(
            mode="makedb",
            input_fasta=input_fasta,
            database=tmp_path / "proteins.dmnd",
        )

        assert "makedb" in cmd
        assert "--in" in cmd
        assert "-d" in cmd
        assert "proteins.faa" in " ".join(cmd)

    def test_blastp_command(self, tmp_path):
        """Should build Diamond blastp command."""
        from metadarkmatter.external.diamond import Diamond

        diamond = Diamond()
        Diamond._executable_cache["diamond"] = Path("/usr/bin/diamond")

        # Create test file for validation
        query = tmp_path / "query.faa"
        query.write_text(">prot1\nMKTL\n")

        cmd = diamond.build_command(
            mode="blastp",
            query=query,
            database=tmp_path / "proteins.dmnd",
            output=tmp_path / "results.tsv",
            threads=8,
            evalue=1e-10,
        )

        assert "blastp" in cmd
        assert "-q" in cmd
        assert "-d" in cmd
        assert "-o" in cmd
        assert "-p" in cmd
        assert "8" in cmd
        assert "1e-10" in cmd

    def test_blastp_sensitive_mode(self, tmp_path):
        """Should include --sensitive flag by default."""
        from metadarkmatter.external.diamond import Diamond

        diamond = Diamond()
        Diamond._executable_cache["diamond"] = Path("/usr/bin/diamond")

        query = tmp_path / "query.faa"
        query.write_text(">prot1\nMKTL\n")

        cmd = diamond.build_command(
            mode="blastp",
            query=query,
            database=tmp_path / "db.dmnd",
            output=tmp_path / "out.tsv",
            sensitive=True,
        )

        assert "--sensitive" in cmd

    def test_blastx_command(self, tmp_path):
        """Should build Diamond blastx command."""
        from metadarkmatter.external.diamond import Diamond

        diamond = Diamond()
        Diamond._executable_cache["diamond"] = Path("/usr/bin/diamond")

        # Create test file for validation
        query = tmp_path / "reads.fastq"
        query.write_text("@read1\nACGT\n+\nIIII\n")

        cmd = diamond.build_command(
            mode="blastx",
            query=query,
            database=tmp_path / "proteins.dmnd",
            output=tmp_path / "results.tsv",
            threads=16,
            evalue=1e-5,
            min_identity=30.0,
        )

        assert "blastx" in cmd
        assert "-q" in cmd
        assert "-d" in cmd
        assert "-o" in cmd
        assert "-p" in cmd
        assert "16" in cmd
        assert "--id" in cmd
        assert "30.0" in cmd

    def test_blastx_sensitive_mode(self, tmp_path):
        """Should include --sensitive flag by default for blastx."""
        from metadarkmatter.external.diamond import Diamond

        diamond = Diamond()
        Diamond._executable_cache["diamond"] = Path("/usr/bin/diamond")

        query = tmp_path / "reads.fastq"
        query.write_text("@read1\nACGT\n+\nIIII\n")

        cmd = diamond.build_command(
            mode="blastx",
            query=query,
            database=tmp_path / "db.dmnd",
            output=tmp_path / "out.tsv",
            sensitive=True,
        )

        assert "--sensitive" in cmd

    def test_blastx_default_outfmt(self, tmp_path):
        """Should use 12-column tabular format by default for blastx."""
        from metadarkmatter.external.diamond import Diamond

        diamond = Diamond()
        Diamond._executable_cache["diamond"] = Path("/usr/bin/diamond")

        query = tmp_path / "reads.fastq"
        query.write_text("@read1\nACGT\n+\nIIII\n")

        cmd = diamond.build_command(
            mode="blastx",
            query=query,
            database=tmp_path / "db.dmnd",
            output=tmp_path / "out.tsv",
        )

        assert "--outfmt" in cmd
        outfmt_idx = cmd.index("--outfmt")
        # Diamond outfmt is split into separate arguments after --outfmt
        # Format is: --outfmt 6 qseqid sseqid pident length ...
        assert cmd[outfmt_idx + 1] == "6"
        # Check that key fields are included as separate arguments
        assert "qseqid" in cmd
        assert "sseqid" in cmd
        assert "pident" in cmd
        assert "bitscore" in cmd

    def test_blastx_custom_outfmt(self, tmp_path):
        """Should allow custom output format for blastx."""
        from metadarkmatter.external.diamond import Diamond

        diamond = Diamond()
        Diamond._executable_cache["diamond"] = Path("/usr/bin/diamond")

        query = tmp_path / "reads.fastq"
        query.write_text("@read1\nACGT\n+\nIIII\n")

        custom_fmt = "6 qseqid sseqid pident"
        cmd = diamond.build_command(
            mode="blastx",
            query=query,
            database=tmp_path / "db.dmnd",
            output=tmp_path / "out.tsv",
            outfmt=custom_fmt,
        )

        assert "--outfmt" in cmd
        outfmt_idx = cmd.index("--outfmt")
        # Custom format is split into separate arguments
        assert cmd[outfmt_idx + 1] == "6"
        assert "qseqid" in cmd
        assert "sseqid" in cmd
        assert "pident" in cmd

    def test_blastx_max_target_seqs(self, tmp_path):
        """Should include max_target_seqs parameter for blastx."""
        from metadarkmatter.external.diamond import Diamond

        diamond = Diamond()
        Diamond._executable_cache["diamond"] = Path("/usr/bin/diamond")

        query = tmp_path / "reads.fastq"
        query.write_text("@read1\nACGT\n+\nIIII\n")

        cmd = diamond.build_command(
            mode="blastx",
            query=query,
            database=tmp_path / "db.dmnd",
            output=tmp_path / "out.tsv",
            max_target_seqs=100,
        )

        assert "-k" in cmd
        assert "100" in cmd

    def test_unknown_mode_raises_error(self, tmp_path):
        """Should raise ValueError for unknown mode."""
        from metadarkmatter.external.diamond import Diamond

        diamond = Diamond()
        Diamond._executable_cache["diamond"] = Path("/usr/bin/diamond")

        with pytest.raises(ValueError, match="Unknown Diamond mode"):
            diamond.build_command(
                mode="invalid",
                query=tmp_path / "query.faa",
                database=tmp_path / "db.dmnd",
                output=tmp_path / "out.tsv",
            )


class TestMMseqs2CommandBuilding:
    """Tests for MMseqs2 command building."""

    def test_createdb_command(self, tmp_path):
        """Should build MMseqs2 createdb command."""
        from metadarkmatter.external.mmseqs2 import MMseqs2

        mmseqs = MMseqs2()
        MMseqs2._executable_cache["mmseqs"] = Path("/usr/bin/mmseqs")

        # Create test file for validation
        input_fasta = tmp_path / "genomes.fasta"
        input_fasta.write_text(">seq1\nACGT\n")

        cmd = mmseqs.build_command(
            mode="createdb",
            input_fasta=input_fasta,
            database=tmp_path / "mmseqs_db",
        )

        assert "createdb" in cmd
        assert "genomes.fasta" in " ".join(cmd)
        assert "mmseqs_db" in " ".join(cmd)

    def test_createdb_with_dbtype(self, tmp_path):
        """Should include dbtype parameter when specified."""
        from metadarkmatter.external.mmseqs2 import MMseqs2

        mmseqs = MMseqs2()
        MMseqs2._executable_cache["mmseqs"] = Path("/usr/bin/mmseqs")

        input_fasta = tmp_path / "genomes.fasta"
        input_fasta.write_text(">seq1\nACGT\n")

        cmd = mmseqs.build_command(
            mode="createdb",
            input_fasta=input_fasta,
            database=tmp_path / "mmseqs_db",
            dbtype=1,  # nucleotide
        )

        assert "--dbtype" in cmd
        assert "1" in cmd

    def test_search_command(self, tmp_path):
        """Should build MMseqs2 easy-search command."""
        from metadarkmatter.external.mmseqs2 import MMseqs2

        mmseqs = MMseqs2()
        MMseqs2._executable_cache["mmseqs"] = Path("/usr/bin/mmseqs")

        query = tmp_path / "reads.fastq"
        query.write_text("@read1\nACGT\n+\nIIII\n")

        cmd = mmseqs.build_command(
            mode="search",
            query=query,
            database=tmp_path / "mmseqs_db",
            output=tmp_path / "results.tsv",
            threads=8,
        )

        assert "easy-search" in cmd
        assert "reads.fastq" in " ".join(cmd)
        assert "--threads" in cmd
        assert "8" in cmd

    def test_search_sensitivity_parameter(self, tmp_path):
        """Should include sensitivity parameter."""
        from metadarkmatter.external.mmseqs2 import MMseqs2

        mmseqs = MMseqs2()
        MMseqs2._executable_cache["mmseqs"] = Path("/usr/bin/mmseqs")

        query = tmp_path / "reads.fastq"
        query.write_text("@read1\nACGT\n+\nIIII\n")

        cmd = mmseqs.build_command(
            mode="search",
            query=query,
            database=tmp_path / "db",
            output=tmp_path / "out.tsv",
            sensitivity=7.0,
        )

        assert "-s" in cmd
        assert "7.0" in cmd

    def test_search_evalue_parameter(self, tmp_path):
        """Should include evalue parameter."""
        from metadarkmatter.external.mmseqs2 import MMseqs2

        mmseqs = MMseqs2()
        MMseqs2._executable_cache["mmseqs"] = Path("/usr/bin/mmseqs")

        query = tmp_path / "reads.fastq"
        query.write_text("@read1\nACGT\n+\nIIII\n")

        cmd = mmseqs.build_command(
            mode="search",
            query=query,
            database=tmp_path / "db",
            output=tmp_path / "out.tsv",
            evalue=1e-5,
        )

        assert "-e" in cmd
        assert "1e-05" in cmd or "0.00001" in cmd

    def test_search_max_seqs_parameter(self, tmp_path):
        """Should include max-seqs parameter."""
        from metadarkmatter.external.mmseqs2 import MMseqs2

        mmseqs = MMseqs2()
        MMseqs2._executable_cache["mmseqs"] = Path("/usr/bin/mmseqs")

        query = tmp_path / "reads.fastq"
        query.write_text("@read1\nACGT\n+\nIIII\n")

        cmd = mmseqs.build_command(
            mode="search",
            query=query,
            database=tmp_path / "db",
            output=tmp_path / "out.tsv",
            max_seqs=100,
        )

        assert "--max-seqs" in cmd
        assert "100" in cmd

    def test_search_min_identity_parameter(self, tmp_path):
        """Should include min-seq-id parameter when specified."""
        from metadarkmatter.external.mmseqs2 import MMseqs2

        mmseqs = MMseqs2()
        MMseqs2._executable_cache["mmseqs"] = Path("/usr/bin/mmseqs")

        query = tmp_path / "reads.fastq"
        query.write_text("@read1\nACGT\n+\nIIII\n")

        cmd = mmseqs.build_command(
            mode="search",
            query=query,
            database=tmp_path / "db",
            output=tmp_path / "out.tsv",
            min_identity=75.0,
        )

        assert "--min-seq-id" in cmd
        # Should convert from 0-100 scale to 0-1 scale
        assert "0.75" in cmd

    def test_search_blast_compatible_format(self, tmp_path):
        """Should use BLAST-compatible output format."""
        from metadarkmatter.external.mmseqs2 import MMseqs2

        mmseqs = MMseqs2()
        MMseqs2._executable_cache["mmseqs"] = Path("/usr/bin/mmseqs")

        query = tmp_path / "reads.fastq"
        query.write_text("@read1\nACGT\n+\nIIII\n")

        cmd = mmseqs.build_command(
            mode="search",
            query=query,
            database=tmp_path / "db",
            output=tmp_path / "out.tsv",
        )

        assert "--format-mode" in cmd
        assert "0" in cmd
        assert "--format-output" in cmd

        # Check for BLAST-compatible column names
        format_idx = cmd.index("--format-output")
        format_str = cmd[format_idx + 1]
        assert "query" in format_str
        assert "target" in format_str
        assert "pident" in format_str
        assert "evalue" in format_str
        assert "bits" in format_str

    def test_search_type_parameter(self, tmp_path):
        """Should include search-type parameter."""
        from metadarkmatter.external.mmseqs2 import MMseqs2

        mmseqs = MMseqs2()
        MMseqs2._executable_cache["mmseqs"] = Path("/usr/bin/mmseqs")

        query = tmp_path / "reads.fastq"
        query.write_text("@read1\nACGT\n+\nIIII\n")

        cmd = mmseqs.build_command(
            mode="search",
            query=query,
            database=tmp_path / "db",
            output=tmp_path / "out.tsv",
            search_type=3,  # nucleotide
        )

        assert "--search-type" in cmd
        assert "3" in cmd

    def test_search_default_parameters(self, tmp_path):
        """Should use sensible defaults when parameters not specified."""
        from metadarkmatter.external.mmseqs2 import MMseqs2

        mmseqs = MMseqs2()
        MMseqs2._executable_cache["mmseqs"] = Path("/usr/bin/mmseqs")

        query = tmp_path / "reads.fastq"
        query.write_text("@read1\nACGT\n+\nIIII\n")

        cmd = mmseqs.build_command(
            mode="search",
            query=query,
            database=tmp_path / "db",
            output=tmp_path / "out.tsv",
        )

        # Default sensitivity: 5.7
        assert "-s" in cmd
        assert "5.7" in cmd

        # Default evalue: 1e-3
        assert "-e" in cmd
        assert "0.001" in cmd or "1e-03" in cmd

        # Default max_seqs: 500
        assert "--max-seqs" in cmd
        assert "500" in cmd

        # Default threads: 4
        assert "--threads" in cmd
        assert "4" in cmd

    def test_unknown_mode_raises_error(self, tmp_path):
        """Should raise ValueError for unknown mode."""
        from metadarkmatter.external.mmseqs2 import MMseqs2

        mmseqs = MMseqs2()
        MMseqs2._executable_cache["mmseqs"] = Path("/usr/bin/mmseqs")

        with pytest.raises(ValueError, match="Unknown MMseqs2 mode"):
            mmseqs.build_command(
                mode="invalid",
                query=tmp_path / "query.fasta",
                database=tmp_path / "db",
                output=tmp_path / "out.tsv",
            )

    def test_createdb_missing_input_raises_error(self):
        """Should raise ValueError when input_fasta is missing for createdb."""
        from metadarkmatter.external.mmseqs2 import MMseqs2

        mmseqs = MMseqs2()
        MMseqs2._executable_cache["mmseqs"] = Path("/usr/bin/mmseqs")

        with pytest.raises(ValueError, match="input_fasta is required"):
            mmseqs.build_command(
                mode="createdb",
                database=Path("db"),
            )

    def test_createdb_missing_database_raises_error(self, tmp_path):
        """Should raise ValueError when database is missing for createdb."""
        from metadarkmatter.external.mmseqs2 import MMseqs2

        mmseqs = MMseqs2()
        MMseqs2._executable_cache["mmseqs"] = Path("/usr/bin/mmseqs")

        input_fasta = tmp_path / "genomes.fasta"
        input_fasta.write_text(">seq1\nACGT\n")

        with pytest.raises(ValueError, match="database is required"):
            mmseqs.build_command(
                mode="createdb",
                input_fasta=input_fasta,
            )

    def test_search_missing_query_raises_error(self):
        """Should raise ValueError when query is missing for search."""
        from metadarkmatter.external.mmseqs2 import MMseqs2

        mmseqs = MMseqs2()
        MMseqs2._executable_cache["mmseqs"] = Path("/usr/bin/mmseqs")

        with pytest.raises(ValueError, match="query is required"):
            mmseqs.build_command(
                mode="search",
                database=Path("db"),
                output=Path("out.tsv"),
            )


class TestDryRunMode:
    """Tests for dry-run execution mode."""

    def test_dry_run_returns_command_without_execution(self):
        """Should return result without executing command."""
        from metadarkmatter.external.kraken import Kraken2

        kraken = Kraken2()
        Kraken2._executable_cache["kraken2"] = Path("/usr/bin/kraken2")

        result = kraken.run(
            reads_1=Path("reads.fastq"),
            database=Path("/db"),
            output=Path("out.kraken"),
            report=Path("out.kreport"),
            dry_run=True,
        )

        assert result.success is True
        assert "[dry-run]" in result.stdout
        assert result.elapsed_seconds == 0.0
        assert "kraken2" in result.command_string
