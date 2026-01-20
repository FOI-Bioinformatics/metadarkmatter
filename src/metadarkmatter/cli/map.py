"""
Map command for competitive read mapping to reference genomes.

Wraps Bowtie2 and samtools to perform competitive mapping of extracted
reads against a pangenome of reference sequences.
"""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console

from metadarkmatter.cli.utils import QuietConsole, spinner_progress
from metadarkmatter.external import (
    Bowtie2,
    Bowtie2Build,
    Samtools,
    ToolExecutionError,
    concatenate_genomes,
)

app = typer.Typer(
    name="map",
    help="Competitive mapping of reads to reference genomes",
    no_args_is_help=True,
)

console = Console()


def _check_tools_available() -> tuple[bool, list[str]]:
    """Check if required tools are available."""
    missing = []

    if not Bowtie2.check_available():
        missing.append("bowtie2")

    if not Bowtie2Build.check_available():
        missing.append("bowtie2-build")

    if not Samtools.check_available():
        missing.append("samtools")

    return len(missing) == 0, missing


def _derive_output_name(reads_1: Path) -> str:
    """Derive output name from reads file."""
    name = reads_1.name
    for suffix in [
        "_family_R1.fastq", "_family_R1.fq",
        "_R1.fastq.gz", "_R1.fq.gz", "_1.fastq.gz", "_1.fq.gz",
        "_R1.fastq", "_R1.fq", "_1.fastq", "_1.fq",
        ".fastq.gz", ".fq.gz", ".fastq", ".fq",
    ]:
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return name.split(".")[0]


@app.command(name="reads")
def map_reads(
    reads_1: Path = typer.Option(
        ...,
        "--reads-1", "-1",
        help="Forward reads file (FASTQ, optionally gzipped)",
        exists=True,
        dir_okay=False,
    ),
    reads_2: Path | None = typer.Option(
        None,
        "--reads-2", "-2",
        help="Reverse reads file for paired-end data",
        exists=True,
        dir_okay=False,
    ),
    genomes: Path = typer.Option(
        ...,
        "--genomes", "-g",
        help="Reference genome directory (FASTA files) or pangenome FASTA",
        exists=True,
    ),
    output: Path = typer.Option(
        ...,
        "--output", "-o",
        help="Output BAM file path",
    ),
    index: Path | None = typer.Option(
        None,
        "--index", "-i",
        help="Existing Bowtie2 index prefix (skip index building)",
    ),
    genome_pattern: str = typer.Option(
        "*.fna",
        "--genome-pattern",
        help="Glob pattern for genome files in directory",
    ),
    mode: str = typer.Option(
        "local",
        "--mode", "-m",
        help="Alignment mode: 'local' or 'end-to-end'",
    ),
    very_sensitive: bool = typer.Option(
        True,
        "--very-sensitive/--fast",
        help="Use very-sensitive preset for divergent sequences (default: very-sensitive)",
    ),
    max_alignments: int = typer.Option(
        500,
        "--max-alignments", "-k",
        help="Report up to this many alignments per read (high value captures all competing hits)",
        min=1,
    ),
    no_unal: bool = typer.Option(
        True,
        "--no-unal/--include-unal",
        help="Exclude unaligned reads from output (default: exclude)",
    ),
    no_discordant: bool = typer.Option(
        False,
        "--no-discordant",
        help="Do not report discordant paired-end alignments",
    ),
    skip_if_exists: bool = typer.Option(
        False,
        "--skip-if-exists",
        help="Skip processing if output BAM already exists",
    ),
    threads: int = typer.Option(
        4,
        "--threads", "-p",
        help="Number of threads",
        min=1,
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
    Competitively map reads against reference genomes.

    This command maps extracted reads against a pangenome constructed from
    multiple reference genomes. Competitive mapping forces reads to align
    to their best matching genome, enabling detection of novel diversity
    through recruitment plot analysis.

    If --genomes is a directory, all genome FASTA files will be concatenated
    into a pangenome with genome-prefixed sequence headers.

    Example:

        # Map to reference genome directory
        metadarkmatter map reads \\
            --reads-1 sample_family_R1.fastq \\
            --reads-2 sample_family_R2.fastq \\
            --genomes ./reference_genomes/ \\
            --output sample_mapped.bam

        # Map to existing pangenome with pre-built index
        metadarkmatter map reads \\
            --reads-1 sample_family_R1.fastq \\
            --genomes pangenome.fasta \\
            --index pangenome_db \\
            --output sample_mapped.bam
    """
    out = QuietConsole(console, quiet=quiet)

    out.print("\n[bold blue]Metadarkmatter Competitive Mapping[/bold blue]\n")

    # Validate mode
    if mode not in ("local", "end-to-end"):
        console.print(f"[red]Error: Invalid mode '{mode}'. Use 'local' or 'end-to-end'.[/red]")
        raise typer.Exit(code=1) from None

    # Check tool availability
    tools_ok, missing_tools = _check_tools_available()
    if not tools_ok:
        console.print("[red]Error: Required tools not found:[/red]")
        for tool in missing_tools:
            console.print(f"  - {tool}")
        console.print("\n[dim]Install with: conda install -c bioconda bowtie2 samtools[/dim]")
        raise typer.Exit(code=1) from None

    # Ensure output has .bam extension
    if output.suffix.lower() != ".bam":
        output = output.with_suffix(".bam")

    # Check skip-if-exists
    if skip_if_exists and output.exists():
        out.print(f"[yellow]Skipping: Output file already exists: {output}[/yellow]")
        raise typer.Exit(code=0)

    # Setup paths
    output.parent.mkdir(parents=True, exist_ok=True)
    output_sam = output.with_suffix(".sam")

    # Determine pangenome and index paths
    pangenome_path: Path
    index_prefix: Path

    if genomes.is_dir():
        # Concatenate genomes from directory
        pangenome_path = output.parent / f"{_derive_output_name(reads_1)}_pangenome.fasta"
        index_prefix = pangenome_path.with_suffix("")
        need_concatenate = True
        need_build_index = index is None
    else:
        # Use provided FASTA as pangenome
        pangenome_path = genomes
        if index is not None:
            index_prefix = index
            need_build_index = False
        else:
            index_prefix = pangenome_path.with_suffix("")
            need_build_index = True
        need_concatenate = False

    # Override if user provided an index
    if index is not None:
        index_prefix = index
        need_build_index = False
        need_concatenate = False

    # Dry run mode
    if dry_run:
        console.print(
            "[bold cyan]DRY RUN MODE[/bold cyan] - Commands shown but not executed\n"
        )

        console.print("[bold]Input Files:[/bold]")
        console.print(f"  Forward reads: {reads_1}")
        if reads_2:
            console.print(f"  Reverse reads: {reads_2}")
        console.print(f"  Genomes: {genomes}")

        console.print("\n[bold]Parameters:[/bold]")
        console.print(f"  Alignment mode: {mode}")
        console.print(f"  Very sensitive: {very_sensitive}")
        console.print(f"  Max alignments: {max_alignments}")
        console.print(f"  No unaligned: {no_unal}")
        console.print(f"  Threads: {threads}")

        console.print("\n[bold]Output Files:[/bold]")
        if need_concatenate:
            console.print(f"  Pangenome: {pangenome_path}")
        console.print(f"  Index prefix: {index_prefix}")
        console.print(f"  Output BAM: {output}")

    # Step 1: Concatenate genomes (if needed)
    genome_names: list[str] = []

    if need_concatenate:
        out.print("[bold]Step 1:[/bold] Concatenating reference genomes...")

        if dry_run:
            console.print(
                f"  [dim]Would concatenate genomes matching "
                f"'{genome_pattern}' from {genomes}[/dim]"
            )
        else:
            with spinner_progress("Concatenating genomes...", console, quiet):
                try:
                    genome_names = concatenate_genomes(
                        genome_dir=genomes,
                        output_path=pangenome_path,
                        pattern=genome_pattern,
                        prefix_headers=True,
                    )
                except FileNotFoundError as e:
                    console.print(f"\n[red]Error: {e}[/red]")
                    raise typer.Exit(code=1) from None

            out.print(f"  [green]Concatenated {len(genome_names)} genomes[/green]")
            if verbose:
                out.print(f"  [dim]Pangenome: {pangenome_path}[/dim]")
    else:
        out.print(f"[bold]Step 1:[/bold] [dim]Using existing pangenome: {pangenome_path}[/dim]")

    # Step 2: Build Bowtie2 index (if needed)
    if need_build_index:
        out.print("\n[bold]Step 2:[/bold] Building Bowtie2 index...")

        builder = Bowtie2Build()

        if dry_run:
            result = builder.run(
                reference=pangenome_path,
                index_prefix=index_prefix,
                threads=threads,
                dry_run=True,
            )
            console.print(f"  [dim]Command: {result.command_string}[/dim]")
        else:
            with spinner_progress("Building index...", console, quiet):
                try:
                    result = builder.run_or_raise(
                        reference=pangenome_path,
                        index_prefix=index_prefix,
                        threads=threads,
                    )
                except ToolExecutionError as e:
                    console.print(f"\n[red]Index building failed:[/red]\n{e.message}")
                    raise typer.Exit(code=1) from None

            out.print(f"  [green]Index built ({result.elapsed_seconds:.1f}s)[/green]")

            if verbose:
                out.print(f"  [dim]Command: {result.command_string}[/dim]")
    else:
        out.print(f"\n[bold]Step 2:[/bold] [dim]Using existing index: {index_prefix}[/dim]")

    # Step 3: Run Bowtie2 alignment
    out.print("\n[bold]Step 3:[/bold] Aligning reads with Bowtie2...")

    aligner = Bowtie2()

    if dry_run:
        result = aligner.run(
            index_prefix=index_prefix,
            reads_1=reads_1,
            reads_2=reads_2,
            output_sam=output_sam,
            threads=threads,
            mode=mode,
            very_sensitive=very_sensitive,
            max_alignments=max_alignments,
            no_unal=no_unal,
            no_discordant=no_discordant,
            dry_run=True,
        )
        console.print(f"  [dim]Command: {result.command_string}[/dim]")
    else:
        with spinner_progress("Aligning reads...", console, quiet):
            try:
                result = aligner.run_or_raise(
                    index_prefix=index_prefix,
                    reads_1=reads_1,
                    reads_2=reads_2,
                    output_sam=output_sam,
                    threads=threads,
                    mode=mode,
                    very_sensitive=very_sensitive,
                    max_alignments=max_alignments,
                    no_unal=no_unal,
                    no_discordant=no_discordant,
                )
            except ToolExecutionError as e:
                console.print(f"\n[red]Alignment failed:[/red]\n{e.message}")
                raise typer.Exit(code=1) from None

        out.print(f"  [green]Alignment complete ({result.elapsed_seconds:.1f}s)[/green]")

        if verbose:
            out.print(f"  [dim]Command: {result.command_string}[/dim]")

        # Parse alignment stats from stderr
        if verbose and result.stderr:
            for line in result.stderr.split("\n"):
                if "aligned" in line.lower() or "overall" in line.lower():
                    out.print(f"  [dim]{line.strip()}[/dim]")

    # Step 4: Convert SAM to sorted BAM
    out.print("\n[bold]Step 4:[/bold] Converting to sorted BAM...")

    samtools = Samtools()

    if dry_run:
        console.print("  [dim]Would run: samtools view -> sort -> index[/dim]")
        console.print("\n[green]Dry run complete. No files were created.[/green]")
        raise typer.Exit(code=0)

    with spinner_progress("Converting and indexing BAM...", console, quiet):
        try:
            view_result, sort_result, index_result = samtools.sam_to_sorted_bam(
                input_sam=output_sam,
                output_bam=output,
                threads=threads,
                cleanup_sam=True,
            )
        except ToolExecutionError as e:
            console.print(f"\n[red]BAM conversion failed:[/red]\n{e.message}")
            raise typer.Exit(code=1) from None

    total_time = (
        view_result.elapsed_seconds
        + sort_result.elapsed_seconds
        + index_result.elapsed_seconds
    )
    out.print(f"  [green]BAM conversion complete ({total_time:.1f}s)[/green]")

    # Summary
    out.print("\n[bold green]Mapping complete![/bold green]")
    out.print("\n[bold]Output files:[/bold]")
    if need_concatenate and pangenome_path.exists():
        out.print(f"  Pangenome:  {pangenome_path}")
    out.print(f"  BAM file:   {output}")
    out.print(f"  BAM index:  {output}.bai")
    out.print()
