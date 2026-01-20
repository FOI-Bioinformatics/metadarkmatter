"""
MMseqs2 command for fast sequence search (BLAST alternative).

Wraps MMseqs2 tools for building databases and running high-speed nucleotide
alignments optimized for large-scale metagenomic datasets. Provides 100-1000x
speedup over BLAST while maintaining sensitivity.
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
from metadarkmatter.external.mmseqs2 import MMseqs2

app = typer.Typer(
    name="mmseqs2",
    help="MMseqs2 fast sequence search (BLAST alternative)",
    no_args_is_help=True,
)

console = Console()


def _check_mmseqs2_available() -> tuple[bool, str]:
    """Check if mmseqs is available in PATH.

    Returns:
        Tuple of (is_available, missing_tool_name). If available,
        returns (True, ""). If not, returns (False, "mmseqs2").
    """
    if MMseqs2.check_available():
        return True, ""
    return False, "mmseqs2"


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
        help="Output database path (e.g., 'mmseqs_db/pangenome')",
    ),
    genome_pattern: str = typer.Option(
        "*.fna",
        "--genome-pattern", "-p",
        help="Glob pattern for genome files when --genomes is a directory",
    ),
    dbtype: int | None = typer.Option(
        None,
        "--dbtype",
        help="Database type: 0=protein, 1=nucleotide, 2=auto (default: auto)",
        min=0,
        max=2,
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
        "--quiet",
        help="Suppress progress output",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Show commands without executing",
    ),
) -> None:
    """
    Create an MMseqs2 nucleotide database from reference genomes.

    The database is used for fast sequence search with 'metadarkmatter mmseqs2 search'.
    Accepts either a single FASTA file or a directory of genome files.

    MMseqs2 provides 100-1000x speedup over BLAST while maintaining sensitivity.

    Example:

        # From directory of genomes
        metadarkmatter mmseqs2 makedb \\
            --genomes reference_genomes/ \\
            --output mmseqs_db/pangenome

        # From single FASTA
        metadarkmatter mmseqs2 makedb \\
            --genomes pangenome.fasta \\
            --output mmseqs_db/pangenome
    """
    out = QuietConsole(console, quiet=quiet)

    out.print("\n[bold blue]Metadarkmatter MMseqs2 Database Builder[/bold blue]\n")

    # Check tool availability
    tool_ok, missing = _check_mmseqs2_available()
    if not tool_ok:
        console.print(f"[red]Error: Required tool not found: {missing}[/red]")
        console.print("\n[dim]Install with: conda install -c bioconda mmseqs2[/dim]")
        raise typer.Exit(code=1) from None

    # Setup output directory
    output.parent.mkdir(parents=True, exist_ok=True)

    # Check skip-if-exists
    # MMseqs2 database files have different extensions than BLAST
    db_files = [output, output.with_suffix(".dbtype"), output.with_suffix(".index")]
    if skip_if_exists and all(f.exists() for f in db_files):
        out.print(f"[yellow]Skipping: MMseqs2 database already exists at {output}[/yellow]")
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
    out.print("\n[bold]Step 2:[/bold] Building MMseqs2 database...")

    mmseqs = MMseqs2()

    if dry_run:
        result = mmseqs.run(
            mode="createdb",
            input_fasta=input_fasta,
            database=output,
            dbtype=dbtype,
            dry_run=True,
        )
        console.print(f"\n[dim]Command: {result.command_string}[/dim]")
        console.print("\n[green]Dry run complete. No files were created.[/green]")
        raise typer.Exit(code=0)

    with spinner_progress("Building MMseqs2 database...", console, quiet):
        try:
            result = mmseqs.run_or_raise(
                mode="createdb",
                input_fasta=input_fasta,
                database=output,
                dbtype=dbtype,
            )
        except ToolExecutionError as e:
            console.print(f"\n[red]MMseqs2 createdb failed:[/red]\n{e.message}")
            raise typer.Exit(code=1) from None

    out.print(f"  [green]Database created ({result.elapsed_seconds:.1f}s)[/green]")

    if verbose:
        out.print(f"\n[dim]Command: {result.command_string}[/dim]")
        if result.stderr.strip():
            out.print(f"[dim]Tool output: {result.stderr.strip()[:200]}[/dim]")

    # Summary
    out.print("\n[bold green]MMseqs2 database ready![/bold green]")
    out.print("\n[bold]Output files:[/bold]")

    # MMseqs2 creates multiple files with different extensions
    db_extensions = ["", ".dbtype", ".index", ".lookup", "_h", "_h.dbtype", "_h.index"]
    for ext in db_extensions:
        db_file = Path(str(output) + ext)
        if db_file.exists():
            size_mb = db_file.stat().st_size / (1024 * 1024)
            out.print(f"  {db_file} ({size_mb:.1f} MB)")

    if verbose:
        out.print(f"[dim]Genomes included: {genome_count}[/dim]")

    out.print()


@app.command(name="search")
def search_reads(
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
        help="MMseqs2 database path (from 'mmseqs2 makedb')",
    ),
    output: Path = typer.Option(
        ...,
        "--output", "-o",
        help="Output file for search results (TSV, optionally .gz)",
    ),
    query_db: Path | None = typer.Option(
        None,
        "--query-db",
        help="Pre-computed query database for caching (optional, speeds up repeated searches)",
    ),
    sensitivity: float = typer.Option(
        5.7,
        "--sensitivity", "-s",
        help="Sensitivity: 1.0 (fast) to 7.5 (sensitive), default 5.7 (balanced)",
        min=1.0,
        max=7.5,
    ),
    evalue: float = typer.Option(
        1e-3,
        "--evalue", "-e",
        help="E-value threshold (default 1e-3)",
    ),
    max_seqs: int = typer.Option(
        500,
        "--max-seqs", "-m",
        help="Max target sequences per query (high for competitive alignment)",
        min=1,
        max=10000,
    ),
    min_identity: float | None = typer.Option(
        None,
        "--min-identity",
        help="Minimum percent identity filter (0-100 scale, optional)",
        min=0.0,
        max=100.0,
    ),
    search_type: int = typer.Option(
        3,
        "--search-type",
        help="Search type: 2=translated, 3=nucleotide (default: 3)",
        min=2,
        max=3,
    ),
    threads: int = typer.Option(
        4,
        "--threads", "-p",
        help="Number of threads",
        min=1,
    ),
    tmp_dir: Path | None = typer.Option(
        None,
        "--tmp-dir",
        help="Temporary directory for MMseqs2 (auto-created if not specified)",
    ),
    compress: bool = typer.Option(
        True,
        "--compress/--no-compress",
        help="Compress output with gzip (default: auto based on .gz extension)",
    ),
    skip_if_exists: bool = typer.Option(
        False,
        "--skip-if-exists",
        help="Skip search if output file already exists",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose", "-v",
        help="Enable verbose output",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet",
        help="Suppress progress output",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Show commands without executing",
    ),
) -> None:
    """
    Run MMseqs2 fast sequence search with BLAST-compatible output.

    Searches query reads against an MMseqs2 database using parameters optimized
    for detecting divergent sequences. Output is in BLAST-compatible tabular format
    with columns: query, target, pident, alnlen, mismatch, gapopen, qstart, qend,
    tstart, tend, evalue, bits.

    MMseqs2 provides 100-1000x speedup over BLAST. Default sensitivity (5.7) provides
    a good balance between speed and accuracy for metagenomic dark matter analysis.

    Example:

        # Basic search
        metadarkmatter mmseqs2 search \\
            --query extracted_reads.fasta \\
            --database mmseqs_db/pangenome \\
            --output sample.mmseqs2.tsv.gz

        # High-sensitivity search
        metadarkmatter mmseqs2 search \\
            --query extracted_reads.fasta \\
            --database mmseqs_db/pangenome \\
            --output sample.mmseqs2.tsv.gz \\
            --sensitivity 7.0 \\
            --evalue 1e-5 \\
            --max-seqs 100 \\
            --threads 16
    """
    out = QuietConsole(console, quiet=quiet)

    out.print("\n[bold blue]Metadarkmatter MMseqs2 Search[/bold blue]\n")

    # Check tool availability
    tool_ok, missing = _check_mmseqs2_available()
    if not tool_ok:
        console.print(f"[red]Error: Required tool not found: {missing}[/red]")
        console.print("\n[dim]Install with: conda install -c bioconda mmseqs2[/dim]")
        raise typer.Exit(code=1) from None

    # Validate database exists
    if not database.exists():
        console.print(f"[red]Error: MMseqs2 database not found at {database}[/red]")
        console.print(
            "[dim]Run 'metadarkmatter mmseqs2 makedb' first to create the database[/dim]"
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
    out.print(f"  Sensitivity:  {sensitivity} (1.0=fast, 7.5=sensitive)")
    out.print(f"  E-value:      {evalue}")
    out.print(f"  Max targets:  {max_seqs}")
    out.print(f"  Threads:      {threads}")
    if min_identity:
        out.print(f"  Min identity: {min_identity}%")

    # Build MMseqs2 command
    mmseqs = MMseqs2()

    # For compressed output, we write to temp then compress
    temp_output = output.with_suffix("").with_suffix(".tsv") if should_compress else output

    if dry_run:
        console.print("\n[dim]Would run the following MMseqs2 workflow:[/dim]")
        console.print(f"[dim]  1. Create query DB: {query} -> queryDB[/dim]")
        console.print(f"[dim]  2. Run search: queryDB + {database} -> resultDB[/dim]")
        console.print(f"[dim]  3. Convert to TSV: resultDB -> {temp_output}[/dim]")
        if should_compress:
            console.print(f"[dim]  4. Compress: {temp_output} -> {output}[/dim]")
        console.print("\n[green]Dry run complete. No files were created.[/green]")
        raise typer.Exit(code=0)

    # Run MMseqs2 multi-step search (optimized workflow)
    out.print("\n[bold]Running optimized MMseqs2 workflow...[/bold]")

    try:
        # Use multistep workflow for better performance
        result_dict = mmseqs.search_multistep(
            query=query,
            database=database,
            output=temp_output,
            query_db=query_db,  # Will be cached if provided
            tmp_dir=tmp_dir,
            threads=threads,
            sensitivity=sensitivity,
            evalue=evalue,
            max_seqs=max_seqs,
            min_identity=min_identity,
            search_type=search_type,
            timeout=None,
        )

        # Report timing
        if result_dict["createdb_time"] > 0:
            out.print(f"  [green]✓ Query DB created ({result_dict['createdb_time']:.1f}s)[/green]")
        else:
            out.print(f"  [green]✓ Query DB reused (cached at {result_dict['query_db']})[/green]")

        out.print(f"  [green]✓ Search complete ({result_dict['search_time']:.1f}s)[/green]")
        out.print(f"  [green]✓ Results converted ({result_dict['convertalis_time']:.1f}s)[/green]")
        out.print(f"  [bold green]Total: {result_dict['total_time']:.1f}s[/bold green]")

        # Save query_db path for future caching
        if query_db is None and verbose:
            out.print(f"\n[dim]Tip: Use --query-db {result_dict['query_db']} to cache the query database[/dim]")
            out.print(f"[dim]This speeds up subsequent searches with the same query file.[/dim]")

    except Exception as e:
        console.print(f"\n[red]MMseqs2 search failed:[/red]\n{str(e)}")
        raise typer.Exit(code=1) from None

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
    out.print("\n[bold green]MMseqs2 search complete![/bold green]")
    out.print("\n[bold]Output:[/bold]")

    file_size = final_output.stat().st_size / (1024 * 1024)
    out.print(f"  {final_output} ({file_size:.1f} MB)")

    if verbose:
        out.print(f"\n[dim]Command: {result.command_string}[/dim]")
        out.print("[dim]Output format: BLAST-compatible tabular[/dim]")
        out.print(
            "[dim]Columns: query, target, pident, alnlen, mismatch, gapopen, "
            "qstart, qend, tstart, tend, evalue, bits[/dim]"
        )

    out.print(
        f"\n[dim]Next step: metadarkmatter score classify "
        f"--blast {final_output} --ani <ani_matrix.csv> "
        f"--output classifications.csv[/dim]"
    )
    out.print()
