"""
BLAST command for competitive alignment of reads against reference genomes.

Wraps BLAST+ tools (makeblastdb, blastn) for building databases and running
high-sensitivity nucleotide alignments optimized for detecting divergent sequences.
"""

from __future__ import annotations

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

    out.print("\n[bold blue]Metadarkmatter BLAST Database Builder[/bold blue]\n")

    # Check tool availability
    tool_ok, missing = _check_makeblastdb_available()
    if not tool_ok:
        console.print(f"[red]Error: Required tool not found: {missing}[/red]")
        console.print("\n[dim]Install with: conda install -c bioconda blast[/dim]")
        raise typer.Exit(code=1) from None

    # Setup output directory
    output.parent.mkdir(parents=True, exist_ok=True)

    # Check skip-if-exists
    db_files = [output.with_suffix(s) for s in [".nhr", ".nin", ".nsq"]]
    if skip_if_exists and all(f.exists() for f in db_files):
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

    with spinner_progress("Building BLAST database...", console, quiet):
        try:
            result = builder.run_or_raise(
                input_fasta=input_fasta,
                output_db=output,
                dbtype="nucl",
                parse_seqids=True,
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
    for suffix in [".nhr", ".nin", ".nsq"]:
        db_file = output.with_suffix(suffix)
        if db_file.exists():
            size_mb = db_file.stat().st_size / (1024 * 1024)
            out.print(f"  {db_file} ({size_mb:.1f} MB)")

    if verbose:
        out.print(f"\n[dim]Database title: {title}[/dim]")
        out.print(f"[dim]Genomes included: {genome_count}[/dim]")

    out.print()


@app.command(name="align")
def align_reads(
    query: Path = typer.Option(
        ...,
        "--query", "-q",
        help="Query reads file (FASTA/FASTQ)",
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
    word_size: int = typer.Option(
        7,
        "--word-size", "-w",
        help="Word size for seeds (smaller = more sensitive, default 7)",
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
    compress: bool = typer.Option(
        True,
        "--compress/--no-compress",
        help="Compress output with gzip (default: auto based on .gz extension)",
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
    qend, sstart, send, evalue, bitscore.

    The default parameters (word_size=7, evalue=1e-3) are optimized for
    detecting remote homology in metagenomic dark matter analysis.

    Example:

        # Basic alignment
        metadarkmatter blast align \\
            --query extracted_reads.fasta \\
            --database blastdb/pangenome \\
            --output sample.blast.tsv.gz

        # High-sensitivity alignment
        metadarkmatter blast align \\
            --query extracted_reads.fasta \\
            --database blastdb/pangenome \\
            --output sample.blast.tsv.gz \\
            --word-size 7 \\
            --evalue 1e-3 \\
            --max-targets 100 \\
            --threads 16
    """
    out = QuietConsole(console, quiet=quiet)

    out.print("\n[bold blue]Metadarkmatter BLAST Alignment[/bold blue]\n")

    # Check tool availability
    tool_ok, missing = _check_blastn_available()
    if not tool_ok:
        console.print(f"[red]Error: Required tool not found: {missing}[/red]")
        console.print("\n[dim]Install with: conda install -c bioconda blast[/dim]")
        raise typer.Exit(code=1) from None

    # Validate database exists
    db_files = [database.with_suffix(s) for s in [".nhr", ".nin", ".nsq"]]
    if not all(f.exists() for f in db_files):
        # Try without suffix (user may have provided full path)
        db_check = Path(str(database) + ".nhr")
        if not db_check.exists():
            console.print(f"[red]Error: BLAST database not found at {database}[/red]")
            console.print(
                "[dim]Run 'metadarkmatter blast makedb' first to create the database[/dim]"
            )
            raise typer.Exit(code=1) from None

    # Setup output
    output.parent.mkdir(parents=True, exist_ok=True)

    # Determine if output should be compressed
    should_compress = compress or output.suffix == ".gz"
    if should_compress and output.suffix != ".gz":
        output = Path(str(output) + ".gz")

    # Check skip-if-exists
    if skip_if_exists and output.exists():
        out.print(f"[yellow]Skipping: Output file already exists: {output}[/yellow]")
        raise typer.Exit(code=0)

    # Show parameters
    out.print("[bold]Input:[/bold]")
    out.print(f"  Query:    {query}")
    out.print(f"  Database: {database}")

    out.print("\n[bold]Parameters:[/bold]")
    out.print(f"  Word size:    {word_size}")
    out.print(f"  E-value:      {evalue}")
    out.print(f"  Max targets:  {max_targets}")
    out.print(f"  Threads:      {threads}")
    if min_identity:
        out.print(f"  Min identity: {min_identity}%")

    # Build BLAST command
    blastn = BlastN()

    # Output format: standard tabular with all 12 columns
    outfmt = (
        "6 qseqid sseqid pident length mismatch gapopen "
        "qstart qend sstart send evalue bitscore"
    )

    # Handle gzipped query input (BLASTN doesn't support gzip directly)
    query_is_gzipped = str(query).endswith(".gz") or str(query).endswith(".gzip")
    temp_query: Path | None = None

    def _cleanup_temp_query() -> None:
        """Ensure temp query file is cleaned up."""
        if temp_query and temp_query.exists():
            try:
                temp_query.unlink()
            except OSError:
                pass  # Best effort cleanup

    if query_is_gzipped:
        import shutil
        import tempfile

        out.print("  [dim]Decompressing query file...[/dim]")
        # Create temp file for decompressed query
        _, temp_query_path = tempfile.mkstemp(suffix=".fastq")
        temp_query = Path(temp_query_path)
        try:
            with gzip.open(query, "rb") as f_in, temp_query.open("wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        except Exception:
            _cleanup_temp_query()
            raise
        query_for_blast = temp_query
    else:
        query_for_blast = query

    # For compressed output, we write to temp then compress
    temp_output = output.with_suffix("").with_suffix(".tsv") if should_compress else output

    if dry_run:
        # Clean up temp query if created
        _cleanup_temp_query()
        result = blastn.run(
            query=query,  # Show original query in dry run
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
        if query_is_gzipped:
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
        with temp_output.open("rb") as f_in, gzip.open(output, "wb") as f_out:
            f_out.writelines(f_in)
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
        f"--blast {final_output} --ani <ani_matrix.csv> "
        f"--output classifications.csv[/dim]"
    )
    out.print()
