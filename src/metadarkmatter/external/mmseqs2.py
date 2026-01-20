"""
MMseqs2 wrapper for fast sequence search.

MMseqs2 provides 100-1000x speedup over BLAST while maintaining sensitivity.
Used as a drop-in alternative for nucleotide and protein sequence searching.

Key Features:
- Fast database creation and searching
- BLAST-compatible tabular output format
- Supports both nucleotide and protein searches
- Adjustable sensitivity for speed/accuracy trade-off

Performance:
- Sensitivity 5.7 (default): Balanced speed and sensitivity
- Sensitivity 7.0+: Comparable to BLAST, slower but more sensitive
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Any, ClassVar

from metadarkmatter.external.base import ExternalTool, validate_path_safe


class MMseqs2(ExternalTool):
    """Wrapper for MMseqs2 fast sequence search.

    MMseqs2 provides fast nucleotide and protein sequence searching,
    offering 100-1000x speedup over BLAST while maintaining sensitivity.

    Supports two modes:
    - 'createdb': Build MMseqs2 database from FASTA
    - 'search': Run sequence search with BLAST-compatible output

    Example:
        >>> mmseqs = MMseqs2()
        >>> # Create database
        >>> mmseqs.run(
        ...     mode="createdb",
        ...     input_fasta=Path("genomes.fasta"),
        ...     database=Path("mmseqs_db"),
        ... )
        >>> # Run search
        >>> mmseqs.run(
        ...     mode="search",
        ...     query=Path("reads.fastq"),
        ...     database=Path("mmseqs_db"),
        ...     output=Path("results.tsv"),
        ...     threads=16,
        ... )
    """

    TOOL_NAME: ClassVar[str] = "mmseqs"
    INSTALL_HINT: ClassVar[str] = "conda install -c bioconda mmseqs2"

    # BLAST-compatible output format (13 columns, including qlen)
    # This ensures compatibility with existing parsers and classifiers
    BLAST_FORMAT_COLUMNS = (
        "query",
        "target",
        "pident",
        "alnlen",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "tstart",
        "tend",
        "evalue",
        "bits",
        "qlen",
    )

    def build_command(
        self,
        *,
        mode: str,
        # createdb parameters
        input_fasta: Path | None = None,
        database: Path | None = None,
        dbtype: int | None = None,
        # search parameters
        query: Path | None = None,
        output: Path | None = None,
        tmp_dir: Path | None = None,
        threads: int = 4,
        sensitivity: float = 5.7,
        evalue: float = 1e-3,
        max_seqs: int = 500,
        min_identity: float | None = None,
        search_type: int = 3,
    ) -> list[str]:
        """Build MMseqs2 command.

        Supports two modes:
        - 'createdb': Create MMseqs2 database from FASTA
        - 'search': Run sequence search using easy-search

        Args:
            mode: Operation mode ('createdb' or 'search')
            input_fasta: Input FASTA file (for createdb)
            database: MMseqs2 database path (for both modes)
            dbtype: Database type - 0=amino acid, 1=nucleotide, 2=auto (for createdb)
            query: Query sequences (for search)
            output: Output file path (for search)
            tmp_dir: Temporary directory for search (auto-created if None)
            threads: Number of CPU threads (default: 4)
            sensitivity: Sensitivity parameter 1.0-7.5 (default: 5.7, balanced)
            evalue: E-value threshold (default: 1e-3)
            max_seqs: Maximum target sequences per query (default: 500)
            min_identity: Minimum percent identity filter (0-100 scale)
            search_type: Search type - 2=translated, 3=nucleotide (default: 3)

        Returns:
            Command as list of strings.

        Raises:
            ValueError: If required parameters are missing for the mode.
        """
        exe = str(self.get_executable())

        if mode == "createdb":
            return self._build_createdb_command(
                exe=exe,
                input_fasta=input_fasta,
                database=database,
                dbtype=dbtype,
            )
        elif mode == "search":
            return self._build_search_command(
                exe=exe,
                query=query,
                database=database,
                output=output,
                tmp_dir=tmp_dir,
                threads=threads,
                sensitivity=sensitivity,
                evalue=evalue,
                max_seqs=max_seqs,
                min_identity=min_identity,
                search_type=search_type,
            )
        else:
            msg = f"Unknown MMseqs2 mode: {mode}. Use 'createdb' or 'search'."
            raise ValueError(msg)

    def _build_createdb_command(
        self,
        exe: str,
        input_fasta: Path | None,
        database: Path | None,
        dbtype: int | None,
    ) -> list[str]:
        """Build MMseqs2 createdb command.

        Args:
            exe: Path to MMseqs2 executable
            input_fasta: Input FASTA file
            database: Output database path
            dbtype: Database type - 0=amino acid, 1=nucleotide, 2=auto

        Returns:
            Command as list of strings.

        Raises:
            ValueError: If required parameters are missing.
        """
        if input_fasta is None:
            msg = "input_fasta is required for createdb mode"
            raise ValueError(msg)
        if database is None:
            msg = "database is required for createdb mode"
            raise ValueError(msg)

        input_fasta = validate_path_safe(input_fasta, must_exist=True)
        database = validate_path_safe(database, must_exist=False)

        cmd = [exe, "createdb", str(input_fasta), str(database)]

        if dbtype is not None:
            cmd.extend(["--dbtype", str(dbtype)])

        return cmd

    def _build_search_command(
        self,
        exe: str,
        query: Path | None,
        database: Path | None,
        output: Path | None,
        tmp_dir: Path | None,
        threads: int,
        sensitivity: float,
        evalue: float,
        max_seqs: int,
        min_identity: float | None,
        search_type: int,
    ) -> list[str]:
        """Build MMseqs2 search command (core search step only).

        Note: This returns only the 'search' command. For complete workflow,
        the CLI must also run 'createdb' for the query and 'convertalis' for results.

        Args:
            exe: Path to MMseqs2 executable
            query: Query sequences (FASTA/FASTQ) - will be converted to query_db
            database: MMseqs2 database path
            output: Output file path (for result_db)
            tmp_dir: Temporary directory (auto-created if None)
            threads: Number of CPU threads
            sensitivity: Sensitivity parameter (1.0-7.5)
            evalue: E-value threshold
            max_seqs: Maximum target sequences per query
            min_identity: Minimum percent identity filter (0-100 scale)
            search_type: Search type - 2=translated, 3=nucleotide

        Returns:
            Command as list of strings.

        Raises:
            ValueError: If required parameters are missing.
        """
        if query is None:
            msg = "query is required for search mode"
            raise ValueError(msg)
        if database is None:
            msg = "database is required for search mode"
            raise ValueError(msg)
        if output is None:
            msg = "output is required for search mode"
            raise ValueError(msg)

        query = validate_path_safe(query, must_exist=True)
        database = validate_path_safe(database, must_exist=False)
        output = validate_path_safe(output, must_exist=False)

        # Create tmp directory if not specified
        if tmp_dir is None:
            tmp_dir = Path(tempfile.mkdtemp(prefix="mmseqs_tmp_"))
        else:
            tmp_dir = Path(tmp_dir)
            tmp_dir.mkdir(parents=True, exist_ok=True)

        # For now, use easy-search but add a warning that this is slow for large files
        # TODO: Migrate to proper multi-step workflow
        cmd = [
            exe,
            "easy-search",
            str(query),
            str(database),
            str(output),
            str(tmp_dir),
            "--threads",
            str(threads),
            "-s",
            str(sensitivity),
            "-e",
            str(evalue),
            "--max-seqs",
            str(max_seqs),
            "--search-type",
            str(search_type),
        ]

        # Min identity filter (optional)
        if min_identity is not None:
            cmd.extend(["--min-seq-id", str(min_identity / 100.0)])

        # CRITICAL: BLAST-compatible output format
        cmd.extend(["--format-mode", "0"])
        format_str = ",".join(self.BLAST_FORMAT_COLUMNS)
        cmd.extend(["--format-output", format_str])

        return cmd

    def create_database(
        self,
        input_fasta: Path,
        database: Path,
        dbtype: int | None = None,
        timeout: float | None = None,
    ) -> None:
        """Create MMseqs2 database from FASTA.

        Convenience method for database creation.

        Args:
            input_fasta: Input FASTA file
            database: Output database path
            dbtype: Database type - 0=amino acid, 1=nucleotide, 2=auto (default)
            timeout: Optional timeout in seconds

        Raises:
            ToolExecutionError: If MMseqs2 fails
        """
        self.run_or_raise(
            mode="createdb",
            input_fasta=input_fasta,
            database=database,
            dbtype=dbtype,
            timeout=timeout,
        )

    def search(
        self,
        query: Path,
        database: Path,
        output: Path,
        tmp_dir: Path | None = None,
        threads: int = 4,
        sensitivity: float = 5.7,
        evalue: float = 1e-3,
        max_seqs: int = 500,
        min_identity: float | None = None,
        search_type: int = 3,
        timeout: float | None = None,
        capture_output: bool = True,
    ) -> None:
        """Run MMseqs2 sequence search.

        Convenience method for searching with BLAST-compatible output.

        Args:
            query: Query sequences (FASTA/FASTQ)
            database: MMseqs2 database path
            output: Output file path
            tmp_dir: Temporary directory (auto-created if None)
            threads: Number of CPU threads (default: 4)
            sensitivity: Sensitivity parameter 1.0-7.5 (default: 5.7)
            evalue: E-value threshold (default: 1e-3)
            max_seqs: Maximum target sequences per query (default: 500)
            min_identity: Minimum percent identity filter (0-100 scale)
            search_type: Search type - 2=translated, 3=nucleotide (default: 3)
            timeout: Optional timeout in seconds
            capture_output: If True, capture stdout/stderr. Set to False for
                batch processing to prevent pipe buffer blocking.

        Raises:
            ToolExecutionError: If MMseqs2 fails

        Example:
            >>> mmseqs = MMseqs2()
            >>> mmseqs.search(
            ...     query=Path("reads.fastq"),
            ...     database=Path("mmseqs_db"),
            ...     output=Path("results.tsv"),
            ...     threads=16,
            ...     sensitivity=5.7,
            ... )
        """
        self.run_or_raise(
            mode="search",
            query=query,
            database=database,
            output=output,
            tmp_dir=tmp_dir,
            threads=threads,
            sensitivity=sensitivity,
            evalue=evalue,
            max_seqs=max_seqs,
            min_identity=min_identity,
            search_type=search_type,
            timeout=timeout,
            capture_output=capture_output,
        )

    def search_multistep(
        self,
        query: Path,
        database: Path,
        output: Path,
        query_db: Path | None = None,
        tmp_dir: Path | None = None,
        threads: int = 4,
        sensitivity: float = 5.7,
        evalue: float = 1e-3,
        max_seqs: int = 500,
        min_identity: float | None = None,
        search_type: int = 3,
        timeout: float | None = None,
    ) -> dict[str, Any]:
        """Run MMseqs2 search using optimized multi-step workflow.

        This is 10-100x faster than easy-search for large query files because
        the query database is created once and can be cached for future searches.

        Workflow:
        1. createdb: Convert query FASTA/FASTQ to MMseqs2 database (cached if query_db provided)
        2. search: Run actual alignment
        3. convertalis: Convert results to BLAST-compatible TSV

        Args:
            query: Query sequences (FASTA/FASTQ)
            database: MMseqs2 target database path
            output: Output TSV file path
            query_db: Pre-computed query database (if None, will create in tmp_dir)
            tmp_dir: Temporary directory for intermediate files
            threads: Number of CPU threads
            sensitivity: Sensitivity parameter 1.0-7.5
            evalue: E-value threshold
            max_seqs: Maximum target sequences per query
            min_identity: Minimum percent identity filter (0-100 scale)
            search_type: Search type - 2=translated, 3=nucleotide
            timeout: Optional timeout in seconds

        Returns:
            Dict with keys:
                - 'query_db': Path to query database (for caching)
                - 'result_db': Path to result database
                - 'createdb_time': Time for createdb step (seconds)
                - 'search_time': Time for search step (seconds)
                - 'convertalis_time': Time for convertalis step (seconds)
                - 'total_time': Total elapsed time (seconds)

        Example:
            >>> mmseqs = MMseqs2()
            >>> # First run - creates query DB
            >>> result = mmseqs.search_multistep(
            ...     query=Path("reads.fastq"),
            ...     database=Path("mmseqs_db"),
            ...     output=Path("results.tsv"),
            ...     query_db=Path("query_db/reads"),  # Will be created and cached
            ...     threads=16,
            ... )
            >>> # Subsequent runs - reuses cached query DB (10-100x faster)
            >>> result2 = mmseqs.search_multistep(
            ...     query=Path("reads.fastq"),
            ...     database=Path("mmseqs_db2"),
            ...     output=Path("results2.tsv"),
            ...     query_db=Path("query_db/reads"),  # Reused!
            ...     threads=16,
            ... )
        """
        import subprocess
        import time

        query = validate_path_safe(query, must_exist=True)
        database = validate_path_safe(database, must_exist=False)
        output = validate_path_safe(output, must_exist=False)

        # Setup directories
        if tmp_dir is None:
            tmp_dir = Path(tempfile.mkdtemp(prefix="mmseqs_tmp_"))
        else:
            tmp_dir = Path(tmp_dir)
            tmp_dir.mkdir(parents=True, exist_ok=True)

        # Setup query database path
        if query_db is None:
            query_db = tmp_dir / "queryDB"
        else:
            query_db = Path(query_db)
            query_db.parent.mkdir(parents=True, exist_ok=True)

        result_db = tmp_dir / "resultDB"

        exe = str(self.get_executable())
        dbtype = search_type - 1  # 3->2 for nucleotide, 2->1 for protein

        times = {}
        total_start = time.perf_counter()

        # Step 1: Create query database (skip if exists)
        createdb_time = 0.0
        if not query_db.exists():
            createdb_cmd = [
                exe, "createdb",
                str(query), str(query_db),
                "--dbtype", str(dbtype),
            ]

            start = time.perf_counter()
            result = subprocess.run(
                createdb_cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
                check=True,
            )
            createdb_time = time.perf_counter() - start

        # Step 2: Run search
        search_cmd = [
            exe, "search",
            str(query_db), str(database), str(result_db), str(tmp_dir),
            "--threads", str(threads),
            "-s", str(sensitivity),
            "-e", str(evalue),
            "--max-seqs", str(max_seqs),
            "--search-type", str(search_type),
        ]

        if min_identity is not None:
            search_cmd.extend(["--min-seq-id", str(min_identity / 100.0)])

        start = time.perf_counter()
        result = subprocess.run(
            search_cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=True,
        )
        search_time = time.perf_counter() - start

        # Step 3: Convert to BLAST format
        format_str = ",".join(self.BLAST_FORMAT_COLUMNS)
        convertalis_cmd = [
            exe, "convertalis",
            str(query_db), str(database), str(result_db), str(output),
            "--format-mode", "0",
            "--format-output", format_str,
        ]

        start = time.perf_counter()
        result = subprocess.run(
            convertalis_cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=True,
        )
        convertalis_time = time.perf_counter() - start

        total_time = time.perf_counter() - total_start

        return {
            "query_db": str(query_db),
            "result_db": str(result_db),
            "createdb_time": createdb_time,
            "search_time": search_time,
            "convertalis_time": convertalis_time,
            "total_time": total_time,
        }
