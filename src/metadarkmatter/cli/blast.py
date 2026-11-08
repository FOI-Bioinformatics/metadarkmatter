"""
BLAST command for competitive alignment of reads against reference genomes.

Wraps BLAST+ tools (makeblastdb, blastn) for building databases and running
high-sensitivity nucleotide alignments optimized for detecting divergent sequences.
"""

from __future__ import annotations

import contextlib
import gzip
from pathlib import Path

import typer
from rich.console import Console
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
)

from metadarkmatter.cli.utils import QuietConsole, spinner_progress
from metadarkmatter.core.runtime import is_dry_run
from metadarkmatter.external import ToolExecutionError
from metadarkmatter.external.blast import BlastN, MakeBlastDb

app = typer.Typer(
    name="blast",
    help="BLAST competitive alignment for ANI-weighted placement",
    no_args_is_help=True,
)

console = Console()


def _check_makeblastdb_available() -> tuple[bool, str]:
    """Check if makeblastdb is available in PATH.

    Returns:
        Tuple of (is_available, missing_tool_name). If available,
        returns (True, ""). If not, returns (False, "makeblastdb (BLAST+)").
    """
    if MakeBlastDb.check_available():
        return True, ""
    return False, "makeblastdb (BLAST+)"


def _check_blastn_available() -> tuple[bool, str]:
    """Check if blastn is available in PATH.

    Returns:
        Tuple of (is_available, missing_tool_name). If available,
        returns (True, ""). If not, returns (False, "blastn (BLAST+)").
    """
    if BlastN.check_available():
        return True, ""
    return False, "blastn (BLAST+)"


# makeblastdb -parse_seqids rejects local IDs longer than this.
_MAX_PARSE_SEQID_LEN = 50

# Default seed word size per BLAST task when --word-size is not given.
# blastn keeps the sensitive word size 7; dc-megablast only accepts 11 or 12.
_TASK_DEFAULT_WORD_SIZE = {
    "blastn": 7,
    "blastn-short": 7,
    "dc-megablast": 11,
    "megablast": 28,
}


def blast_db_exists(prefix: Path) -> bool:
    """True when a finished BLAST nucleotide database exists at ``prefix``.

    Accepts a single-volume DB (``<prefix>.nhr``), a multi-volume alias
    (``<prefix>.nal``) or a first volume (``<prefix>.00.nhr``). Suffixes are
    appended to the full prefix string so dotted prefixes such as
    ``db/pangenome.v2`` are not truncated by ``Path.with_suffix``.
    """
    return any(
        Path(f"{prefix}{suffix}").is_file() for suffix in (".nhr", ".nal", ".00.nhr")
    )


def _fasta_ids_fit_parse_seqids(fasta: Path) -> bool:
    """True when every FASTA ID is short enough for ``makeblastdb -parse_seqids``.

    Standardized ``{accession}|{contig_id}`` headers built from assembler
    contig names (e.g. SPAdes ``NODE_..._length_..._cov_...``) easily exceed
    the 50-character local-ID limit, which makes makeblastdb fail outright.
    """
    opener = gzip.open if str(fasta).endswith(".gz") else open
    try:
        with opener(fasta, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    parts = line[1:].split(maxsplit=1)
                    if parts and len(parts[0]) > _MAX_PARSE_SEQID_LEN:
                        return False
    except OSError:
        return True
    return True


def _sniff_sequence_suffix(path: Path) -> str:
    """Return ``.fastq`` or ``.fasta`` from the first record marker of a file."""
    opener = gzip.open if str(path).endswith((".gz", ".gzip")) else open
    try:
        with opener(path, "rt") as fh:
            for line in fh:
                stripped = line.strip()
                if stripped:
                    return ".fastq" if stripped.startswith("@") else ".fasta"
    except (OSError, EOFError):
        pass
    inner = Path(Path(path).stem).suffix.lower()
    return ".fastq" if inner in (".fastq", ".fq") else ".fasta"


def _concatenate_genomes(
    genome_dir: Path,
    output_fasta: Path,
    contig_mapping_path: Path,
    pattern: str = "*.fna",
) -> tuple[int, int]:
    """Concatenate genome files with standardized headers and create contig mapping.

    Rewrites FASTA headers to format: {accession}|{original_contig_id}
    Creates a contig mapping TSV for tracking which contigs belong to which genome.

    Returns:
        Tuple of (genome_count, contig_count)
    """
    from metadarkmatter.core.genome_utils import concatenate_genomes_with_mapping

    return concatenate_genomes_with_mapping(
        genome_dir=genome_dir,
        output_fasta=output_fasta,
        contig_mapping_path=contig_mapping_path,
        pattern=pattern,
    )


@app.command(name="makedb")
def make_database(
    genomes: Path = typer.Option(
        ...,
        "--genomes", "-g",
        help="Path to genome FASTA file or directory of genome files",
        exists=True,
    ),
    output: Path = typer.Option(
        ...,
        "--output", "-o",
        help="Output database prefix (e.g., 'blastdb/pangenome')",
    ),
    title: str | None = typer.Option(
        None,
        "--title", "-t",
        help="Database title (default: derived from output name)",
    ),
    genome_pattern: str = typer.Option(
        "*.fna",
        "--genome-pattern", "-p",
        help="Glob pattern for genome files when --genomes is a directory",
    ),
    skip_if_exists: bool = typer.Option(
        False,
        "--skip-if-exists",
        help="Skip database creation if output files already exist",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose", "-v",
        help="Enable verbose output",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet", "-q",
        help="Suppress progress output",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Show commands without executing",
    ),
) -> None:
    """
    Create a BLAST nucleotide database from reference genomes.

    The database is used for competitive alignment with 'metadarkmatter blast align'.
    Accepts either a single FASTA file or a directory of genome files.

    Example:

        # From directory of genomes
        metadarkmatter blast makedb \\
            --genomes reference_genomes/ \\
            --output blastdb/pangenome

        # From single FASTA
        metadarkmatter blast makedb \\
            --genomes pangenome.fasta \\
            --output blastdb/pangenome
    """
    out = QuietConsole(console, quiet=quiet)

    # makeblastdb cannot take a database name containing whitespace
    # ("Please provide a database name using -out"); fail before writing.
    if any(ch.isspace() for ch in str(output)):
        console.print(
            f"[red]Error: the BLAST database path must not contain spaces: {output}[/red]"
        )
        raise typer.Exit(code=1)

    out.print("\n[bold blue]Metadarkmatter BLAST Database Builder[/bold blue]\n")

    dry_run = dry_run or is_dry_run()

    # Check tool availability (a dry run only prints the plan)
    tool_ok, missing = _check_makeblastdb_available()
    if not tool_ok and not dry_run:
        console.print(f"[red]Error: Required tool not found: {missing}[/red]")
        console.print("\n[dim]Install with: conda install -c bioconda blast[/dim]")
        raise typer.Exit(code=1) from None

    # Setup output directory
    if not dry_run:
        output.parent.mkdir(parents=True, exist_ok=True)

    # Check skip-if-exists
    if skip_if_exists and blast_db_exists(output):
        out.print(f"[yellow]Skipping: BLAST database already exists at {output}[/yellow]")
        raise typer.Exit(code=0)

    # Determine input FASTA
    if genomes.is_dir():
        out.print(f"[bold]Step 1:[/bold] Concatenating genomes from {genomes}...")
        pangenome_fasta = output.parent / f"{output.name}_pangenome.fasta"
        contig_mapping_path = output.parent / f"{output.name}_contig_mapping.tsv"

        if dry_run:
            console.print(
                f"  [dim]Would concatenate genomes matching "
                f"'{genome_pattern}' to {pangenome_fasta}[/dim]"
            )
            console.print(
                f"  [dim]Would create contig mapping at {contig_mapping_path}[/dim]"
            )
            input_fasta = pangenome_fasta
            genome_count = len(list(genomes.glob(genome_pattern)))
        else:
            genome_count, contig_count = _concatenate_genomes(
                genomes, pangenome_fasta, contig_mapping_path, genome_pattern
            )
            if genome_count == 0:
                console.print(
                    f"[red]Error: No genome files found matching "
                    f"'{genome_pattern}' in {genomes}[/red]"
                )
                raise typer.Exit(code=1) from None
            out.print(
                f"  [green]Concatenated {genome_count} genomes "
                f"({contig_count} contigs)[/green]"
            )
            out.print(f"  [green]Created contig mapping: {contig_mapping_path}[/green]")
            input_fasta = pangenome_fasta
    else:
        out.print(f"[bold]Input:[/bold] Using {genomes}")
        input_fasta = genomes
        genome_count = 1

    # Build database
    out.print("\n[bold]Step 2:[/bold] Building BLAST database...")

    if title is None:
        title = output.name

    builder = MakeBlastDb()

    if dry_run:
        if not tool_ok:
            console.print(f"[yellow]Note: {missing} not found; required for a real run.[/yellow]")
            console.print("\n[green]Dry run complete. No files were created.[/green]")
            raise typer.Exit(code=0)
        result = builder.run(
            input_fasta=input_fasta,
            output_db=output,
            dbtype="nucl",
            parse_seqids=True,
            title=title,
            dry_run=True,
        )
        console.print(f"\n[dim]Command: {result.command_string}[/dim]")
        console.print("\n[green]Dry run complete. No files were created.[/green]")
        raise typer.Exit(code=0)

    parse_seqids = _fasta_ids_fit_parse_seqids(input_fasta)
    if not parse_seqids:
        out.print(
            f"  [yellow]Note: some sequence IDs exceed {_MAX_PARSE_SEQID_LEN} characters; "
            "building without -parse_seqids (sseqids in alignments are unchanged).[/yellow]"
        )

    with spinner_progress("Building BLAST database...", console, quiet):
        try:
            result = builder.run_or_raise(
                input_fasta=input_fasta,
                output_db=output,
                dbtype="nucl",
                parse_seqids=parse_seqids,
                title=title,
            )
        except ToolExecutionError as e:
            console.print(f"\n[red]makeblastdb failed:[/red]\n{e.message}")
            raise typer.Exit(code=1) from None

    out.print(f"  [green]Database created ({result.elapsed_seconds:.1f}s)[/green]")

    if verbose:
        out.print(f"\n[dim]Command: {result.command_string}[/dim]")
        if result.stderr.strip():
            out.print(f"[dim]Tool output: {result.stderr.strip()[:200]}[/dim]")

    # Summary
    out.print("\n[bold green]BLAST database ready![/bold green]")
    out.print("\n[bold]Output files:[/bold]")
    for suffix in [".nhr", ".nin", ".nsq", ".nal"]:
        db_file = Path(f"{output}{suffix}")
        if db_file.exists():
            size_mb = db_file.stat().st_size / (1024 * 1024)
            out.print(f"  {db_file} ({size_mb:.1f} MB)")

    if verbose:
        out.print(f"\n[dim]Database title: {title}[/dim]")
        out.print(f"[dim]Genomes included: {genome_count}[/dim]")

    out.print()


@app.command(name="align")
def align_reads(
    query: Path | None = typer.Option(
        None,
        "--query",
        help="Single-end query reads (FASTA/FASTQ). Use this OR --query-1/--query-2.",
        exists=True,
        dir_okay=False,
    ),
    query_1: Path | None = typer.Option(
        None,
        "--query-1", "-1",
        help="Forward (R1) reads for paired-end data (FASTA/FASTQ, gzipped OK)",
        exists=True,
        dir_okay=False,
    ),
    query_2: Path | None = typer.Option(
        None,
        "--query-2", "-2",
        help="Reverse (R2) reads for paired-end data (FASTA/FASTQ, gzipped OK)",
        exists=True,
        dir_okay=False,
    ),
    database: Path = typer.Option(
        ...,
        "--database", "-d",
        help="BLAST database prefix (from 'blast makedb')",
    ),
    output: Path = typer.Option(
        ...,
        "--output", "-o",
        help="Output file for BLAST results (TSV, optionally .gz)",
    ),
    word_size: int | None = typer.Option(
        None,
        "--word-size", "-w",
        help=(
            "Word size for seeds (smaller = more sensitive). Default depends on "
            "--task: 7 for blastn, 11 for dc-megablast, 28 for megablast."
        ),
        min=4,
        max=28,
    ),
    evalue: float = typer.Option(
        1e-5,
        "--evalue", "-e",
        help="E-value threshold (default 1e-5 for 75%+ identity)",
    ),
    max_targets: int = typer.Option(
        500,
        "--max-targets", "-m",
        help="Max target sequences per query (high for competitive alignment)",
        min=1,
        max=10000,
    ),
    min_identity: float | None = typer.Option(
        None,
        "--min-identity",
        help="Minimum percent identity filter (optional)",
        min=0.0,
        max=100.0,
    ),
    task: str = typer.Option(
        "blastn",
        "--task",
        help="BLAST task (blastn, megablast, dc-megablast)",
    ),
    threads: int = typer.Option(
        4,
        "--threads", "-p",
        help="Number of threads",
        min=1,
    ),
    compress: bool | None = typer.Option(
        None,
        "--compress/--no-compress",
        help=(
            "Compress output with gzip (default: auto, compress only when --output "
            "ends in .gz; --compress appends .gz)"
        ),
    ),
    skip_if_exists: bool = typer.Option(
        False,
        "--skip-if-exists",
        help="Skip alignment if output file already exists",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose", "-v",
        help="Enable verbose output",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet", "-q",
        help="Suppress progress output",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Show commands without executing",
    ),
) -> None:
    """
    Run BLASTN competitive alignment for ANI-weighted placement.

    Aligns query reads against a reference database using parameters optimized
    for detecting divergent sequences. Output is in tabular format (outfmt 6)
    with columns: qseqid, sseqid, pident, length, mismatch, gapopen, qstart,
    qend, sstart, send, evalue, bitscore, qlen.

    Defaults: --task blastn with --word-size following the task (blastn 7,
    dc-megablast 11, megablast 28), --evalue 1e-5, --max-targets 500.

    Supports single-end (--query) and paired-end (--query-1, --query-2) input.
    For paired-end data both mates are aligned in one BLAST run and written to
    one output file; mate IDs get /1 and /2 suffixes unless already present
    (the same convention as 'mmseqs2 search'), so mates are classified as
    separate reads and each read's hits stay contiguous.

    Example:

        # Basic alignment
        metadarkmatter blast align \\
            --query extracted_reads.fasta \\
            --database blastdb/pangenome \\
            --output sample.blast.tsv.gz

        # Paired-end alignment (both mates, one combined output)
        metadarkmatter blast align \\
            --query-1 extracted_R1.fastq.gz \\
            --query-2 extracted_R2.fastq.gz \\
            --database blastdb/pangenome \\
            --output sample.blast.tsv.gz

        # Faster alignment against closely related references
        metadarkmatter blast align \\
            --query extracted_reads.fasta \\
            --database blastdb/pangenome \\
            --output sample.blast.tsv.gz \\
            --task dc-megablast \\
            --threads 16
    """
    out = QuietConsole(console, quiet=quiet)

    out.print("\n[bold blue]Metadarkmatter BLAST Alignment[/bold blue]\n")

    dry_run = dry_run or is_dry_run()

    # Validate input options (mutual exclusivity)
    has_single_end = query is not None
    has_paired_end = query_1 is not None or query_2 is not None
    if has_single_end and has_paired_end:
        console.print(
            "[red]Error: Cannot use --query with --query-1/--query-2. "
            "Choose one input mode.[/red]"
        )
        raise typer.Exit(code=1) from None
    if not has_single_end and not has_paired_end:
        console.print(
            "[red]Error: Must provide either --query (single-end) or "
            "--query-1 and --query-2 (paired-end).[/red]"
        )
        raise typer.Exit(code=1) from None
    if has_paired_end and (query_1 is None or query_2 is None):
        console.print(
            "[red]Error: Paired-end mode requires both --query-1 and --query-2.[/red]"
        )
        raise typer.Exit(code=1) from None
    if query_1 is not None and query_2 is not None and (
        _sniff_sequence_suffix(query_1) != _sniff_sequence_suffix(query_2)
    ):
        console.print(
            "[red]Error: --query-1 and --query-2 must use the same format "
            "(both FASTA or both FASTQ).[/red]"
        )
        raise typer.Exit(code=1) from None

    # Check tool availability
    tool_ok, missing = _check_blastn_available()
    if not tool_ok and not dry_run:
        console.print(f"[red]Error: Required tool not found: {missing}[/red]")
        console.print("\n[dim]Install with: conda install -c bioconda blast[/dim]")
        raise typer.Exit(code=1) from None

    # Validate database exists (a dry run may plan against a DB built earlier
    # in the same pipeline, so only warn there)
    if not blast_db_exists(database):
        if not dry_run:
            console.print(f"[red]Error: BLAST database not found at {database}[/red]")
            console.print(
                "[dim]Run 'metadarkmatter blast makedb' first to create the database[/dim]"
            )
            raise typer.Exit(code=1) from None
        console.print(f"[yellow]Note: BLAST database not found at {database}[/yellow]")

    if word_size is None:
        word_size = _TASK_DEFAULT_WORD_SIZE.get(task, 7)

    # Setup output
    if not dry_run:
        output.parent.mkdir(parents=True, exist_ok=True)

    # Compression follows the output extension unless --compress is explicit.
    should_compress = bool(compress) or output.suffix == ".gz"
    if should_compress and output.suffix != ".gz":
        output = Path(str(output) + ".gz")

    # Check skip-if-exists
    if skip_if_exists and output.exists():
        out.print(f"[yellow]Skipping: Output file already exists: {output}[/yellow]")
        raise typer.Exit(code=0)

    # Show parameters
    out.print("[bold]Input:[/bold]")
    if has_paired_end:
        out.print("  Paired-end reads")
        out.print(f"  R1:       {query_1}")
        out.print(f"  R2:       {query_2}")
    else:
        out.print(f"  Query:    {query}")
    out.print(f"  Database: {database}")

    out.print("\n[bold]Parameters:[/bold]")
    out.print(f"  Task:         {task}")
    out.print(f"  Word size:    {word_size}")
    out.print(f"  E-value:      {evalue}")
    out.print(f"  Max targets:  {max_targets}")
    out.print(f"  Threads:      {threads}")
    if min_identity:
        out.print(f"  Min identity: {min_identity}%")

    # Build BLAST command
    blastn = BlastN()

    # Output format: standard tabular with all 13 columns (including qlen)
    outfmt = (
        "6 qseqid sseqid pident length mismatch gapopen "
        "qstart qend sstart send evalue bitscore qlen"
    )

    # Prepare the query file. Gzipped input is decompressed (BLASTN does not
    # read gzip); paired-end input is combined into one file with distinct
    # /1 and /2 mate IDs so each mate is a separate read in the output.
    temp_query: Path | None = None

    def _cleanup_temp_query() -> None:
        """Ensure temp query file is cleaned up."""
        if temp_query and temp_query.exists():
            with contextlib.suppress(OSError):
                temp_query.unlink()  # Best effort cleanup

    query_is_gzipped = query is not None and str(query).endswith((".gz", ".gzip"))
    query_for_blast: Path
    if has_paired_end:
        assert query_1 is not None and query_2 is not None
        if dry_run:
            query_for_blast = Path("<paired R1+R2>")
        else:
            import os
            import tempfile

            from metadarkmatter.core.paired_reads import combine_mates

            out.print("  [dim]Combining paired-end reads (mates tagged /1 and /2)...[/dim]")
            fd, temp_query_path = tempfile.mkstemp(suffix=_sniff_sequence_suffix(query_1))
            os.close(fd)
            temp_query = Path(temp_query_path)
            try:
                n_reads = combine_mates(query_1, query_2, temp_query)
            except Exception:
                _cleanup_temp_query()
                raise
            out.print(f"  [dim]Combined {n_reads:,} reads from paired-end files[/dim]")
            query_for_blast = temp_query
    elif query_is_gzipped and not dry_run:
        assert query is not None
        import os
        import shutil
        import tempfile

        out.print("  [dim]Decompressing query file...[/dim]")
        # Keep the real format in the temp suffix: BlastN converts FASTQ to
        # FASTA based on it, and treating FASTA as FASTQ drops records.
        fd, temp_query_path = tempfile.mkstemp(suffix=_sniff_sequence_suffix(query))
        os.close(fd)
        temp_query = Path(temp_query_path)
        try:
            with gzip.open(query, "rb") as f_in, temp_query.open("wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        except Exception:
            _cleanup_temp_query()
            raise
        query_for_blast = temp_query
    else:
        assert query is not None
        query_for_blast = query

    # For compressed output, we write to temp then compress
    # Uncompressed temp file next to the output; a name derived from the
    # output itself so it can never be an unrelated existing file.
    temp_output = output.with_name(output.name + ".partial") if should_compress else output

    if dry_run:
        # Clean up temp query if created
        _cleanup_temp_query()
        if not tool_ok:
            console.print(f"[yellow]Note: {missing} not found; required for a real run.[/yellow]")
            console.print("\n[green]Dry run complete. No files were created.[/green]")
            raise typer.Exit(code=0)
        result = blastn.run(
            query=query_for_blast,  # Show original query (or paired placeholder) in dry run
            database=database,
            output=temp_output,
            outfmt=outfmt,
            task=task,
            word_size=word_size,
            evalue=evalue,
            max_target_seqs=max_targets,
            threads=threads,
            perc_identity=min_identity,
            dry_run=True,
        )
        console.print(f"\n[dim]Command: {result.command_string}[/dim]")
        if has_paired_end:
            console.print(
                f"[dim](R1 {query_1} and R2 {query_2} will be combined into one temp "
                "query with /1 and /2 mate suffixes before alignment)[/dim]"
            )
        elif query_is_gzipped:
            console.print("[dim](Query will be decompressed to temp file before alignment)[/dim]")
        if should_compress:
            console.print(f"[dim]Then: gzip {temp_output} -> {output}[/dim]")
        console.print("\n[green]Dry run complete. No files were created.[/green]")
        raise typer.Exit(code=0)

    # Run BLAST
    out.print("\n[bold]Running BLASTN alignment...[/bold]")

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TimeElapsedColumn(),
        console=console if not quiet else None,
    ) as progress:
        progress.add_task(
            description=f"Aligning reads (threads={threads})...",
            total=None,
        )

        try:
            result = blastn.run_or_raise(
                query=query_for_blast,
                database=database,
                output=temp_output,
                outfmt=outfmt,
                task=task,
                word_size=word_size,
                evalue=evalue,
                max_target_seqs=max_targets,
                threads=threads,
                perc_identity=min_identity,
            )
        except ToolExecutionError as e:
            # Clean up temp query file on error
            _cleanup_temp_query()
            console.print(f"\n[red]BLASTN failed:[/red]\n{e.message}")
            raise typer.Exit(code=1) from None

    # Clean up temp query file
    _cleanup_temp_query()

    out.print(f"  [green]Alignment complete ({result.elapsed_seconds:.1f}s)[/green]")

    # Compress output if needed
    if should_compress and temp_output.exists():
        out.print("  Compressing output...")
        with temp_output.open("rb") as src_in, gzip.open(output, "wb") as gz_out:
            gz_out.writelines(src_in)
        temp_output.unlink()

    # Count alignments
    final_output = output if should_compress else temp_output
    if final_output.exists():
        if final_output.suffix == ".gz":
            with gzip.open(final_output, "rt") as f:
                alignment_count = sum(1 for _ in f)
        else:
            with final_output.open() as f:
                alignment_count = sum(1 for _ in f)

        out.print(f"  [dim]Total alignments: {alignment_count:,}[/dim]")

    # Summary
    out.print("\n[bold green]BLAST alignment complete![/bold green]")
    out.print("\n[bold]Output:[/bold]")

    file_size = final_output.stat().st_size / (1024 * 1024)
    out.print(f"  {final_output} ({file_size:.1f} MB)")

    if verbose:
        out.print(f"\n[dim]Command: {result.command_string}[/dim]")
        out.print("[dim]Output format: tabular (outfmt 6)[/dim]")
        out.print(
            "[dim]Columns: qseqid, sseqid, pident, length, mismatch, gapopen, "
            "qstart, qend, sstart, send, evalue, bitscore[/dim]"
        )
        if result.stderr.strip():
            out.print(f"[dim]Tool messages: {result.stderr.strip()[:200]}[/dim]")

    out.print(
        f"\n[dim]Next step: metadarkmatter score classify "
        f"--alignment {final_output} --ani <ani_matrix.csv> "
        f"--output classifications.csv[/dim]"
    )
    out.print()
