"""
BLAST+ wrapper classes.

Provides Python interfaces for:
- MakeBlastDb: Creating BLAST databases
- BlastN: Nucleotide-nucleotide alignment
"""

from __future__ import annotations

import gzip
import os
import subprocess
import tempfile
from pathlib import Path

from metadarkmatter.external.base import ExternalTool, validate_path_safe


class MakeBlastDb(ExternalTool):
    """Wrapper for makeblastdb database construction.

    Creates a BLAST database from FASTA sequences for use
    with blastn and other BLAST+ tools.

    Example:
        >>> builder = MakeBlastDb()
        >>> result = builder.run(
        ...     input_fasta=Path("pangenome.fasta"),
        ...     output_db=Path("pangenome_blastdb"),
        ...     dbtype="nucl",
        ... )
    """

    TOOL_NAME = "makeblastdb"
    INSTALL_HINT = "conda install -c bioconda blast"

    def build_command(
        self,
        *,
        input_fasta: Path,
        output_db: Path,
        dbtype: str = "nucl",
        parse_seqids: bool = True,
        title: str | None = None,
        taxid: int | None = None,
        logfile: Path | None = None,
    ) -> list[str]:
        """Build makeblastdb command.

        Args:
            input_fasta: Input FASTA file.
            output_db: Output database prefix.
            dbtype: Database type - "nucl" or "prot".
            parse_seqids: Parse sequence IDs for retrieval.
            title: Title for the database.
            taxid: Taxonomy ID to assign to sequences.
            logfile: Path for log output.

        Returns:
            Command as list of strings.
        """
        # Validate input paths for safety (resolve and check for unsafe chars)
        # Note: must_exist=False allows command building for dry-run/testing
        input_fasta = validate_path_safe(input_fasta, must_exist=False)
        output_db = validate_path_safe(output_db, must_exist=False)
        if logfile:
            logfile = validate_path_safe(logfile, must_exist=False)

        exe = str(self.get_executable())
        cmd = [exe]

        # Input
        cmd.extend(["-in", str(input_fasta)])

        # Output
        cmd.extend(["-out", str(output_db)])

        # Database type
        cmd.extend(["-dbtype", dbtype])

        # Options
        if parse_seqids:
            cmd.append("-parse_seqids")

        if title:
            cmd.extend(["-title", title])

        if taxid is not None:
            cmd.extend(["-taxid", str(taxid)])

        if logfile:
            cmd.extend(["-logfile", str(logfile)])

        return cmd


class BlastN(ExternalTool):
    """Wrapper for blastn nucleotide alignment.

    Performs nucleotide-nucleotide sequence alignment using
    the BLAST algorithm. Optimized for metagenomic read mapping.

    Example:
        >>> blastn = BlastN()
        >>> result = blastn.run(
        ...     query=Path("reads.fasta"),
        ...     database=Path("pangenome_blastdb"),
        ...     output=Path("blast_results.tsv"),
        ...     word_size=7,
        ...     max_target_seqs=10,
        ... )
    """

    TOOL_NAME = "blastn"
    INSTALL_HINT = "conda install -c bioconda blast"

    # Standard BLAST tabular format columns
    OUTFMT_6_COLUMNS = (
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
    )

    def build_command(
        self,
        *,
        query: Path,
        database: Path,
        output: Path,
        outfmt: str = "6 qseqid sseqid pident length bitscore",
        task: str = "blastn",
        word_size: int = 7,
        evalue: float = 1e-3,
        max_target_seqs: int = 10,
        threads: int = 4,
        perc_identity: float | None = None,
        dust: str | None = None,
        soft_masking: bool | None = None,
    ) -> list[str]:
        """Build blastn command.

        Args:
            query: Query sequences (FASTA/FASTQ).
            database: BLAST database prefix.
            output: Output file path.
            outfmt: Output format string.
            task: BLAST task - "blastn", "megablast", etc.
            word_size: Word size for seed matches (smaller = more sensitive).
            evalue: E-value threshold.
            max_target_seqs: Maximum target sequences per query.
            threads: Number of threads.
            perc_identity: Minimum percent identity filter.
            dust: DUST filtering ("yes", "no", or parameters).
            soft_masking: Use soft masking for filtering.

        Returns:
            Command as list of strings.
        """
        # Validate input paths for safety (resolve and check for unsafe chars)
        # Note: must_exist=False allows command building for dry-run/testing
        query = validate_path_safe(query, must_exist=False)
        database = validate_path_safe(database, must_exist=False)
        output = validate_path_safe(output, must_exist=False)

        exe = str(self.get_executable())
        cmd = [exe]

        # Required inputs
        cmd.extend(["-query", str(query)])
        cmd.extend(["-db", str(database)])
        cmd.extend(["-out", str(output)])

        # Output format
        cmd.extend(["-outfmt", outfmt])

        # Algorithm parameters
        cmd.extend(["-task", task])
        cmd.extend(["-word_size", str(word_size)])
        cmd.extend(["-evalue", str(evalue)])
        cmd.extend(["-max_target_seqs", str(max_target_seqs)])

        # Threading
        cmd.extend(["-num_threads", str(threads)])

        # Optional filters
        if perc_identity is not None:
            cmd.extend(["-perc_identity", str(perc_identity)])

        if dust is not None:
            cmd.extend(["-dust", dust])

        if soft_masking is not None:
            cmd.extend(["-soft_masking", "true" if soft_masking else "false"])

        return cmd

    @staticmethod
    def _is_fastq_file(file_path: Path) -> bool:
        """Check if file is FASTQ based on extension.

        Args:
            file_path: Path to check.

        Returns:
            True if file appears to be FASTQ format.
        """
        name_lower = file_path.name.lower()
        return name_lower.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz'))

    @staticmethod
    def _convert_fastq_to_fasta(fastq_path: Path, fasta_path: Path) -> None:
        """Convert FASTQ to FASTA using seqtk or fallback to Python.

        Args:
            fastq_path: Input FASTQ file (can be gzipped).
            fasta_path: Output FASTA file.

        Raises:
            RuntimeError: If conversion fails.
        """
        # Try seqtk first (fastest)
        try:
            cmd = ["seqtk", "seq", "-A", str(fastq_path)]
            with open(fasta_path, 'w') as f:
                result = subprocess.run(
                    cmd,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                )
            return
        except (subprocess.CalledProcessError, FileNotFoundError):
            # seqtk not available or failed, use Python fallback
            pass

        # Python fallback - handle gzipped and plain FASTQ
        try:
            open_func = gzip.open if str(fastq_path).endswith('.gz') else open

            with open_func(fastq_path, 'rt') as fin, open(fasta_path, 'w') as fout:
                line_num = 0
                for line in fin:
                    line_num += 1
                    mod = line_num % 4

                    if mod == 1:  # Header line
                        # Convert @header to >header
                        fout.write('>' + line[1:])
                    elif mod == 2:  # Sequence line
                        fout.write(line)
                    # Skip quality header (mod==3) and quality scores (mod==0)

        except Exception as e:
            msg = f"Failed to convert FASTQ to FASTA: {e}"
            raise RuntimeError(msg) from e

    def run(
        self,
        *,
        timeout: float | None = None,
        dry_run: bool = False,
        capture_output: bool = True,
        **kwargs: object,
    ):
        """Execute BLASTN with automatic FASTQ to FASTA conversion.

        If the query file is FASTQ format, it will be automatically converted
        to FASTA before running BLAST. The temporary FASTA file is cleaned up
        after execution.

        Args:
            timeout: Maximum execution time in seconds.
            dry_run: If True, return command without execution.
            capture_output: If True, capture stdout/stderr.
            **kwargs: Arguments passed to build_command(), including 'query'.

        Returns:
            ToolResult with command, exit code, and output.

        Raises:
            RuntimeError: If FASTQ conversion fails.
        """
        query = kwargs.get('query')
        temp_fasta = None

        # Check if we need to convert FASTQ to FASTA
        if query and isinstance(query, Path) and self._is_fastq_file(query):
            # Create temporary FASTA file
            fd, temp_fasta_str = tempfile.mkstemp(suffix='.fasta', prefix='blast_query_')
            os.close(fd)
            temp_fasta = Path(temp_fasta_str)

            try:
                # Convert FASTQ to FASTA
                self._convert_fastq_to_fasta(query, temp_fasta)

                # Replace query with temporary FASTA
                kwargs['query'] = temp_fasta

                # Run BLAST with FASTA
                result = super().run(
                    timeout=timeout,
                    dry_run=dry_run,
                    capture_output=capture_output,
                    **kwargs,
                )

                return result

            finally:
                # Clean up temporary file
                if temp_fasta and temp_fasta.exists():
                    temp_fasta.unlink()
        else:
            # Already FASTA or no query provided, run normally
            return super().run(
                timeout=timeout,
                dry_run=dry_run,
                capture_output=capture_output,
                **kwargs,
            )

    @classmethod
    def build_sensitive_command_args(cls) -> dict[str, object]:
        """Return parameters optimized for sensitive remote homology detection.

        These settings are recommended for detecting divergent sequences
        in metagenomic dark matter analysis.

        Returns:
            Dictionary of kwargs for build_command().
        """
        return {
            "task": "blastn",
            "word_size": 7,
            "evalue": 1e-3,
            "outfmt": (
                "6 qseqid sseqid pident length mismatch gapopen "
                "qstart qend sstart send evalue bitscore"
            ),
        }
