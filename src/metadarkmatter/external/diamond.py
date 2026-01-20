"""
Diamond wrapper for protein alignment and AAI computation.

Diamond is a fast protein aligner, used here for all-vs-all protein BLAST
(blastp) to compute Average Amino Acid Identity (AAI) between genomes.
AAI provides a more reliable metric than ANI for genus-level classification,
as ANI becomes unreliable below approximately 80% identity.

AAI genus boundaries based on Riesco & Trujillo 2024:
- Same genus: AAI > 65%
- Genus boundary zone: AAI 58-65%
- Different genus: AAI < 58%
"""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

from metadarkmatter.external.base import ExternalTool, validate_path_safe


class Diamond(ExternalTool):
    """Wrapper for Diamond protein aligner.

    Diamond provides fast protein-level BLAST alignment, suitable for
    computing AAI matrices from annotated genome protein files (.faa).

    The makedb command creates a Diamond database from protein sequences.
    The blastp command performs all-vs-all protein alignment for AAI calculation.

    Example:
        >>> diamond = Diamond()
        >>> # Create database
        >>> diamond.run(
        ...     mode="makedb",
        ...     input_fasta=Path("proteins.faa"),
        ...     database=Path("proteins.dmnd"),
        ... )
        >>> # Run all-vs-all blastp
        >>> diamond.run(
        ...     mode="blastp",
        ...     query=Path("proteins.faa"),
        ...     database=Path("proteins.dmnd"),
        ...     output=Path("hits.tsv"),
        ...     threads=16,
        ... )
    """

    TOOL_NAME: ClassVar[str] = "diamond"
    INSTALL_HINT: ClassVar[str] = "conda install -c bioconda diamond"

    def build_command(
        self,
        *,
        mode: str,
        # makedb parameters
        input_fasta: Path | None = None,
        database: Path | None = None,
        # blastp/blastx parameters
        query: Path | None = None,
        output: Path | None = None,
        threads: int = 4,
        evalue: float = 1e-5,
        max_target_seqs: int = 500,
        min_identity: float = 30.0,
        sensitive: bool = True,
        outfmt: str | None = None,
    ) -> list[str]:
        """Build Diamond command.

        Supports three modes:
        - 'makedb': Create Diamond database from protein FASTA
        - 'blastp': Run protein-vs-protein BLAST alignment
        - 'blastx': Run translated nucleotide-vs-protein BLAST alignment

        Args:
            mode: Operation mode ('makedb', 'blastp', or 'blastx')
            input_fasta: Input FASTA file (for makedb)
            database: Diamond database path (for both modes)
            query: Query FASTA file (for blastp/blastx)
            output: Output file path (for blastp/blastx)
            threads: Number of CPU threads (default: 4)
            evalue: E-value threshold (default: 1e-5)
            max_target_seqs: Maximum targets per query (default: 500)
            min_identity: Minimum percent identity filter (default: 30%)
            sensitive: Use sensitive mode for divergent sequences (default: True)
            outfmt: Custom output format (default: 6-column for AAI, 12-column for blastx)

        Returns:
            Command as list of strings.

        Raises:
            ValueError: If required parameters are missing for the mode.
        """
        exe = str(self.get_executable())

        if mode == "makedb":
            return self._build_makedb_command(
                exe=exe,
                input_fasta=input_fasta,
                database=database,
            )
        elif mode == "blastp":
            return self._build_blastp_command(
                exe=exe,
                query=query,
                database=database,
                output=output,
                threads=threads,
                evalue=evalue,
                max_target_seqs=max_target_seqs,
                min_identity=min_identity,
                sensitive=sensitive,
                outfmt=outfmt,
            )
        elif mode == "blastx":
            return self._build_blastx_command(
                exe=exe,
                query=query,
                database=database,
                output=output,
                threads=threads,
                evalue=evalue,
                max_target_seqs=max_target_seqs,
                min_identity=min_identity,
                sensitive=sensitive,
                outfmt=outfmt,
            )
        else:
            msg = f"Unknown Diamond mode: {mode}. Use 'makedb', 'blastp', or 'blastx'."
            raise ValueError(msg)

    def _build_makedb_command(
        self,
        exe: str,
        input_fasta: Path | None,
        database: Path | None,
    ) -> list[str]:
        """Build Diamond makedb command.

        Args:
            exe: Path to Diamond executable
            input_fasta: Input protein FASTA file
            database: Output database path

        Returns:
            Command as list of strings.

        Raises:
            ValueError: If required parameters are missing.
        """
        if input_fasta is None:
            msg = "input_fasta is required for makedb mode"
            raise ValueError(msg)
        if database is None:
            msg = "database is required for makedb mode"
            raise ValueError(msg)

        input_fasta = validate_path_safe(input_fasta, must_exist=True)
        database = validate_path_safe(database, must_exist=False)

        return [
            exe,
            "makedb",
            "--in",
            str(input_fasta),
            "-d",
            str(database),
        ]

    def _build_blastp_command(
        self,
        exe: str,
        query: Path | None,
        database: Path | None,
        output: Path | None,
        threads: int,
        evalue: float,
        max_target_seqs: int,
        min_identity: float,
        sensitive: bool,
        outfmt: str | None,
    ) -> list[str]:
        """Build Diamond blastp command.

        Args:
            exe: Path to Diamond executable
            query: Query protein FASTA file
            database: Diamond database path
            output: Output file path
            threads: Number of CPU threads
            evalue: E-value threshold
            max_target_seqs: Maximum targets per query
            min_identity: Minimum percent identity filter
            sensitive: Use sensitive mode
            outfmt: Custom output format

        Returns:
            Command as list of strings.

        Raises:
            ValueError: If required parameters are missing.
        """
        if query is None:
            msg = "query is required for blastp mode"
            raise ValueError(msg)
        if database is None:
            msg = "database is required for blastp mode"
            raise ValueError(msg)
        if output is None:
            msg = "output is required for blastp mode"
            raise ValueError(msg)

        query = validate_path_safe(query, must_exist=True)
        database = validate_path_safe(database, must_exist=False)
        output = validate_path_safe(output, must_exist=False)

        cmd = [
            exe,
            "blastp",
            "-q",
            str(query),
            "-d",
            str(database),
            "-o",
            str(output),
            "-p",
            str(threads),
            "-e",
            str(evalue),
            "-k",
            str(max_target_seqs),
            "--id",
            str(min_identity),
        ]

        # Add sensitive mode for divergent sequences (AAI often < 70%)
        if sensitive:
            cmd.append("--sensitive")

        # Output format for AAI computation:
        # qseqid sseqid pident length qlen slen
        # This provides the data needed for reciprocal best hit filtering
        # and AAI calculation (identity weighted by alignment)
        # Diamond expects each format field as a separate argument
        if outfmt is None:
            outfmt = "6 qseqid sseqid pident length qlen slen"
        cmd.append("--outfmt")
        cmd.extend(outfmt.split())

        return cmd

    def _build_blastx_command(
        self,
        exe: str,
        query: Path | None,
        database: Path | None,
        output: Path | None,
        threads: int,
        evalue: float,
        max_target_seqs: int,
        min_identity: float,
        sensitive: bool,
        outfmt: str | None,
    ) -> list[str]:
        """Build Diamond blastx command.

        BLASTX translates nucleotide query sequences in all six reading frames
        and aligns them against a protein database.

        Args:
            exe: Path to Diamond executable
            query: Query nucleotide FASTA/FASTQ file
            database: Diamond protein database path
            output: Output file path
            threads: Number of CPU threads
            evalue: E-value threshold
            max_target_seqs: Maximum targets per query
            min_identity: Minimum percent identity filter
            sensitive: Use sensitive mode
            outfmt: Custom output format

        Returns:
            Command as list of strings.

        Raises:
            ValueError: If required parameters are missing.
        """
        if query is None:
            msg = "query is required for blastx mode"
            raise ValueError(msg)
        if database is None:
            msg = "database is required for blastx mode"
            raise ValueError(msg)
        if output is None:
            msg = "output is required for blastx mode"
            raise ValueError(msg)

        query = validate_path_safe(query, must_exist=True)
        database = validate_path_safe(database, must_exist=False)
        output = validate_path_safe(output, must_exist=False)

        cmd = [
            exe,
            "blastx",
            "-q",
            str(query),
            "-d",
            str(database),
            "-o",
            str(output),
            "-p",
            str(threads),
            "-e",
            str(evalue),
            "-k",
            str(max_target_seqs),
            "--id",
            str(min_identity),
        ]

        # Add sensitive mode for divergent sequences
        if sensitive:
            cmd.append("--sensitive")

        # Output format for read classification:
        # Standard 12-column tabular format matching BLASTN output
        # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        # Diamond expects each format field as a separate argument
        if outfmt is None:
            outfmt = (
                "6 qseqid sseqid pident length mismatch gapopen "
                "qstart qend sstart send evalue bitscore"
            )
        cmd.append("--outfmt")
        cmd.extend(outfmt.split())

        return cmd

    def make_database(
        self,
        input_fasta: Path,
        database: Path,
        timeout: float | None = None,
    ) -> None:
        """Create Diamond database from protein FASTA.

        Convenience method for database creation.

        Args:
            input_fasta: Input protein FASTA file (.faa)
            database: Output database path (will add .dmnd extension)
            timeout: Optional timeout in seconds

        Raises:
            ToolExecutionError: If Diamond fails
        """
        self.run_or_raise(
            mode="makedb",
            input_fasta=input_fasta,
            database=database,
            timeout=timeout,
        )

    def blastp(
        self,
        query: Path,
        database: Path,
        output: Path,
        threads: int = 4,
        evalue: float = 1e-5,
        min_identity: float = 30.0,
        sensitive: bool = True,
        timeout: float | None = None,
        capture_output: bool = True,
    ) -> None:
        """Run Diamond blastp alignment.

        Convenience method for protein BLAST.

        Args:
            query: Query protein FASTA file
            database: Diamond database path
            output: Output file path
            threads: Number of CPU threads
            evalue: E-value threshold
            min_identity: Minimum percent identity filter
            sensitive: Use sensitive mode for divergent sequences
            timeout: Optional timeout in seconds
            capture_output: If True, capture stdout/stderr. Set to False for
                batch processing to prevent pipe buffer blocking.

        Raises:
            ToolExecutionError: If Diamond fails
        """
        self.run_or_raise(
            mode="blastp",
            query=query,
            database=database,
            output=output,
            threads=threads,
            evalue=evalue,
            min_identity=min_identity,
            sensitive=sensitive,
            timeout=timeout,
            capture_output=capture_output,
        )

    def blastx(
        self,
        query: Path,
        database: Path,
        output: Path,
        threads: int = 4,
        evalue: float = 1e-5,
        min_identity: float = 30.0,
        max_target_seqs: int = 500,
        sensitive: bool = True,
        outfmt: str | None = None,
        timeout: float | None = None,
    ) -> None:
        """Run Diamond BLASTX alignment (DNA query vs protein database).

        BLASTX translates nucleotide query sequences in all six reading frames
        and aligns them against a protein database. This enables detection of
        divergent homology that may not be detectable at the nucleotide level.

        For highly divergent taxa (genus-level novelty), nucleotide-based
        alignment may miss homology that is detectable at the protein level
        because amino acid sequences are more conserved than nucleotide sequences.

        Args:
            query: Query nucleotide FASTA/FASTQ file (DNA reads)
            database: Diamond protein database path (from makedb)
            output: Output file path for tabular results
            threads: Number of CPU threads (default: 4)
            evalue: E-value threshold (default: 1e-5)
            min_identity: Minimum percent identity filter (default: 30%)
            max_target_seqs: Maximum targets per query (default: 500)
            sensitive: Use sensitive mode for divergent sequences (default: True)
            outfmt: Custom output format (default: 12-column tabular)
            timeout: Optional timeout in seconds

        Raises:
            ToolExecutionError: If Diamond fails

        Example:
            >>> diamond = Diamond()
            >>> diamond.blastx(
            ...     query=Path("reads.fastq"),
            ...     database=Path("proteins.dmnd"),
            ...     output=Path("blastx_results.tsv"),
            ...     threads=16,
            ... )
        """
        self.run_or_raise(
            mode="blastx",
            query=query,
            database=database,
            output=output,
            threads=threads,
            evalue=evalue,
            min_identity=min_identity,
            max_target_seqs=max_target_seqs,
            sensitive=sensitive,
            outfmt=outfmt,
            timeout=timeout,
        )
