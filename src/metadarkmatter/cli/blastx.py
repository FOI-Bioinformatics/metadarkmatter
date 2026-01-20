"""
BLASTX commands for protein-level read classification.

Wraps Diamond BLASTX for translated nucleotide-vs-protein alignment,
enabling detection of divergent homology that may not be detectable
at the nucleotide level.

For highly divergent taxa (genus-level novelty), nucleotide-based alignment
may miss homology that is detectable at the protein level because amino acid
sequences are more conserved than nucleotide sequences.
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
from metadarkmatter.external.diamond import Diamond

app = typer.Typer(
    name="blastx",
    help="Protein-level read classification using Diamond BLASTX",
    no_args_is_help=True,
)

console = Console()


def _check_diamond_available() -> tuple[bool, str]:
    """Check if Diamond is available in PATH.

    Returns:
        Tuple of (is_available, missing_tool_name). If available,
        returns (True, ""). If not, returns (False, "diamond").
    """
    if Diamond.check_available():
        return True, ""
    return False, "diamond"


def _concatenate_proteins(
    protein_dir: Path,
    output_fasta: Path,
    pattern: str = "*.faa",
) -> tuple[int, int]:
    """Concatenate protein FASTA files with genome-prefixed headers.

    Rewrites FASTA headers to format: {accession}|{original_protein_id}
    This matches the header format used for nucleotide genomes.

    Args:
        protein_dir: Directory containing protein FASTA files
        pattern: Glob pattern for protein files (default: *.faa)
        output_fasta: Output concatenated FASTA file

    Returns:
        Tuple of (genome_count, protein_count)
    """
    protein_files = sorted(protein_dir.glob(pattern))
    genome_count = 0
    protein_count = 0

    with output_fasta.open("w") as out_f:
        for protein_file in protein_files:
            # Extract accession from filename (e.g., GCA_000001.1.faa -> GCA_000001.1)
            accession = protein_file.stem
            if accession.endswith("_protein"):
                accession = accession[:-8]  # Remove _protein suffix

            genome_count += 1

            # Handle gzipped files
            if str(protein_file).endswith(".gz"):
                open_func = gzip.open
                mode = "rt"
            else:
                open_func = open
                mode = "r"

            with open_func(protein_file, mode) as in_f:
                for line in in_f:
                    if line.startswith(">"):
                        # Extract original protein ID (first word after >)
                        original_id = line[1:].split()[0]
                        # Rewrite header with genome prefix
                        out_f.write(f">{accession}|{original_id}\n")
                        protein_count += 1
                    else:
                        out_f.write(line)

    return genome_count, protein_count


@app.command(name="makedb")
def make_database(
    proteins: Path = typer.Option(
        ...,
        "--proteins", "-p",
        help="Path to protein FASTA file or directory of .faa files",
        exists=True,
    ),
    output: Path = typer.Option(
        ...,
        "--output", "-o",
        help="Output database prefix (e.g., 'blastdb/panproteome')",
    ),
    protein_pattern: str = typer.Option(
        "*.faa",
        "--protein-pattern",
        help="Glob pattern for protein files when --proteins is a directory",
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
    Create a Diamond protein database from reference genome proteins.

    The database is used for BLASTX alignment with 'metadarkmatter blastx align'.
    Accepts either a single protein FASTA file or a directory of .faa files.

    Protein files can be obtained from NCBI using:

        datasets download genome accession <accession> --include protein

    Example:

        # From directory of protein files
        metadarkmatter blastx makedb \\
            --proteins reference_proteins/ \\
            --output blastdb/panproteome

        # From single FASTA
        metadarkmatter blastx makedb \\
            --proteins all_proteins.faa \\
            --output blastdb/panproteome
    """
    out = QuietConsole(console, quiet=quiet)

    out.print("\n[bold blue]Metadarkmatter Diamond Protein Database Builder[/bold blue]\n")

    # Check tool availability
    tool_ok, missing = _check_diamond_available()
    if not tool_ok:
        console.print(f"[red]Error: Required tool not found: {missing}[/red]")
        console.print("\n[dim]Install with: conda install -c bioconda diamond[/dim]")
        raise typer.Exit(code=1) from None

    # Setup output directory
    output.parent.mkdir(parents=True, exist_ok=True)

    # Check skip-if-exists
    db_file = output.with_suffix(".dmnd")
    if skip_if_exists and db_file.exists():
        out.print(f"[yellow]Skipping: Diamond database already exists at {db_file}[/yellow]")
        raise typer.Exit(code=0)

    # Determine input FASTA
    if proteins.is_dir():
        out.print(f"[bold]Step 1:[/bold] Concatenating proteins from {proteins}...")
        panproteome_fasta = output.parent / f"{output.name}_panproteome.faa"

        if dry_run:
            console.print(
                f"  [dim]Would concatenate proteins matching "
                f"'{protein_pattern}' to {panproteome_fasta}[/dim]"
            )
            input_fasta = panproteome_fasta
            genome_count = len(list(proteins.glob(protein_pattern)))
        else:
            genome_count, protein_count = _concatenate_proteins(
                proteins, panproteome_fasta, protein_pattern
            )
            if genome_count == 0:
                console.print(
                    f"[red]Error: No protein files found matching "
                    f"'{protein_pattern}' in {proteins}[/red]"
                )
                raise typer.Exit(code=1) from None
            out.print(
                f"  [green]Concatenated {genome_count} genomes "
                f"({protein_count:,} proteins)[/green]"
            )
            input_fasta = panproteome_fasta
    else:
        out.print(f"[bold]Input:[/bold] Using {proteins}")
        input_fasta = proteins
        genome_count = 1

    # Build database
    out.print("\n[bold]Step 2:[/bold] Building Diamond protein database...")

    diamond = Diamond()

    if dry_run:
        result = diamond.run(
            mode="makedb",
            input_fasta=input_fasta,
            database=output,
            dry_run=True,
        )
        console.print(f"\n[dim]Command: {result.command_string}[/dim]")
        console.print("\n[green]Dry run complete. No files were created.[/green]")
        raise typer.Exit(code=0)

    with spinner_progress("Building Diamond database...", console, quiet):
        try:
            result = diamond.run_or_raise(
                mode="makedb",
                input_fasta=input_fasta,
                database=output,
            )
        except ToolExecutionError as e:
            console.print(f"\n[red]Diamond makedb failed:[/red]\n{e.message}")
            raise typer.Exit(code=1) from None

    out.print(f"  [green]Database created ({result.elapsed_seconds:.1f}s)[/green]")

    if verbose:
        out.print(f"\n[dim]Command: {result.command_string}[/dim]")
        if result.stderr.strip():
            out.print(f"[dim]Tool output: {result.stderr.strip()[:200]}[/dim]")

    # Summary
    out.print("\n[bold green]Diamond protein database ready![/bold green]")
    out.print("\n[bold]Output files:[/bold]")

    if db_file.exists():
        size_mb = db_file.stat().st_size / (1024 * 1024)
        out.print(f"  {db_file} ({size_mb:.1f} MB)")

    if verbose:
        out.print(f"\n[dim]Genomes included: {genome_count}[/dim]")

    out.print(
        f"\n[dim]Next step: metadarkmatter blastx align "
        f"--query <reads.fastq.gz> --database {output} "
        f"--output results.blastx.tsv.gz[/dim]"
    )
    out.print()


@app.command(name="align")
def align_reads(
    query: Path = typer.Option(
        ...,
        "--query", "-q",
        help="Query reads file (FASTA/FASTQ, optionally gzipped)",
        exists=True,
        dir_okay=False,
    ),
    database: Path = typer.Option(
        ...,
        "--database", "-d",
        help="Diamond database prefix (from 'blastx makedb')",
    ),
    output: Path = typer.Option(
        ...,
        "--output", "-o",
        help="Output file for BLASTX results (TSV, optionally .gz)",
    ),
    evalue: float = typer.Option(
        1e-5,
        "--evalue", "-e",
        help="E-value threshold (default 1e-5)",
    ),
    max_targets: int = typer.Option(
        500,
        "--max-targets", "-k",
        help="Max target sequences per query (high for competitive alignment)",
        min=1,
        max=10000,
    ),
    min_identity: float = typer.Option(
        30.0,
        "--min-identity",
        help="Minimum percent identity filter (default 30%)",
        min=0.0,
        max=100.0,
    ),
    sensitive: bool = typer.Option(
        True,
        "--sensitive/--fast",
        help="Use sensitive mode for divergent sequences (default: sensitive)",
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
    Run Diamond BLASTX for protein-level read classification.

    BLASTX translates nucleotide query sequences in all six reading frames
    and aligns them against a protein database. This enables detection of
    divergent homology that may not be detectable at the nucleotide level.

    Output is in tabular format (outfmt 6) with columns: qseqid, sseqid,
    pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue,
    bitscore. This format is compatible with the existing BLAST parsers.

    Example:

        # Basic alignment
        metadarkmatter blastx align \\
            --query extracted_reads.fastq.gz \\
            --database blastdb/panproteome \\
            --output sample.blastx.tsv.gz

        # High-sensitivity alignment
        metadarkmatter blastx align \\
            --query extracted_reads.fastq.gz \\
            --database blastdb/panproteome \\
            --output sample.blastx.tsv.gz \\
            --sensitive \\
            --max-targets 500 \\
            --threads 16
    """
    out = QuietConsole(console, quiet=quiet)

    out.print("\n[bold blue]Metadarkmatter Diamond BLASTX Alignment[/bold blue]\n")

    # Check tool availability
    tool_ok, missing = _check_diamond_available()
    if not tool_ok:
        console.print(f"[red]Error: Required tool not found: {missing}[/red]")
        console.print("\n[dim]Install with: conda install -c bioconda diamond[/dim]")
        raise typer.Exit(code=1) from None

    # Validate database exists
    db_file = database.with_suffix(".dmnd")
    if not db_file.exists():
        console.print(f"[red]Error: Diamond database not found at {db_file}[/red]")
        console.print(
            "[dim]Run 'metadarkmatter blastx makedb' first to create the database[/dim]"
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
    out.print(f"  E-value:      {evalue}")
    out.print(f"  Max targets:  {max_targets}")
    out.print(f"  Min identity: {min_identity}%")
    out.print(f"  Sensitive:    {sensitive}")
    out.print(f"  Threads:      {threads}")

    diamond = Diamond()

    # For compressed output, we write to temp then compress
    temp_output = output.with_suffix("").with_suffix(".tsv") if should_compress else output

    if dry_run:
        result = diamond.run(
            mode="blastx",
            query=query,
            database=database,
            output=temp_output,
            threads=threads,
            evalue=evalue,
            max_target_seqs=max_targets,
            min_identity=min_identity,
            sensitive=sensitive,
            dry_run=True,
        )
        console.print(f"\n[dim]Command: {result.command_string}[/dim]")
        if should_compress:
            console.print(f"[dim]Then: gzip {temp_output} -> {output}[/dim]")
        console.print("\n[green]Dry run complete. No files were created.[/green]")
        raise typer.Exit(code=0)

    # Run BLASTX
    out.print("\n[bold]Running Diamond BLASTX alignment...[/bold]")

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
            result = diamond.run_or_raise(
                mode="blastx",
                query=query,
                database=database,
                output=temp_output,
                threads=threads,
                evalue=evalue,
                max_target_seqs=max_targets,
                min_identity=min_identity,
                sensitive=sensitive,
            )
        except ToolExecutionError as e:
            console.print(f"\n[red]Diamond BLASTX failed:[/red]\n{e.message}")
            raise typer.Exit(code=1) from None

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
    out.print("\n[bold green]BLASTX alignment complete![/bold green]")
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
        f"--alignment-mode protein --output classifications.csv[/dim]"
    )
    out.print()
