"""
Kraken2 command for taxonomic classification and read extraction.

Provides two subcommands:
- classify: Run Kraken2 taxonomic classification
- extract: Extract reads for a specific taxid from Kraken2 output
"""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console

from metadarkmatter.cli.utils import (
    QuietConsole,
    extract_sample_name,
    extract_sample_name_from_reads,
    spinner_progress,
)
from metadarkmatter.external import (
    ExtractKrakenReads,
    Kraken2,
    KrakenReport,
    ToolExecutionError,
)

app = typer.Typer(
    name="kraken2",
    help="Kraken2 taxonomic classification and read extraction",
    no_args_is_help=True,
)

console = Console()


@app.command(name="classify")
def classify_reads(
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
    kraken_db: Path = typer.Option(
        ...,
        "--kraken-db", "-k",
        help="Path to Kraken2 database directory",
        exists=True,
        file_okay=False,
    ),
    output: Path = typer.Option(
        ...,
        "--output", "-o",
        help="Output directory for Kraken2 files",
    ),
    sample_name: str | None = typer.Option(
        None,
        "--sample-name", "-n",
        help="Sample name (default: derived from reads filename)",
    ),
    confidence: float = typer.Option(
        0.0,
        "--confidence",
        help="Kraken2 confidence threshold (0-1)",
        min=0.0,
        max=1.0,
    ),
    minimum_hit_groups: int = typer.Option(
        2,
        "--minimum-hit-groups",
        help="Minimum hit groups for classification",
        min=1,
    ),
    skip_if_exists: bool = typer.Option(
        False,
        "--skip-if-exists",
        help="Skip classification if output files already exist",
    ),
    threads: int = typer.Option(
        4,
        "--threads", "-p",
        help="Number of threads for Kraken2",
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
    Run Kraken2 taxonomic classification on metagenomic reads.

    This command classifies reads against a Kraken2 database and generates
    output files that can be used with 'extract reads' to extract reads
    belonging to specific taxa.

    Example:

        # Classify reads
        metadarkmatter extract classify \\
            --reads-1 sample_R1.fastq.gz \\
            --reads-2 sample_R2.fastq.gz \\
            --kraken-db /path/to/kraken_db \\
            --output ./kraken_output/

    Output files:
        - {sample}.kraken: Per-read classification output
        - {sample}.kreport: Hierarchical taxonomic report
    """
    out = QuietConsole(console, quiet=quiet)

    out.print("\n[bold blue]Metadarkmatter Kraken2 Classification[/bold blue]\n")

    # Check tool availability
    if not Kraken2.check_available():
        console.print("[red]Error: kraken2 not found[/red]")
        console.print("\n[dim]Install with: conda install -c bioconda kraken2[/dim]")
        raise typer.Exit(code=1) from None

    # Derive sample name
    if sample_name is None:
        sample_name = extract_sample_name_from_reads(reads_1)

    # Setup output paths
    output.mkdir(parents=True, exist_ok=True)

    kraken_output = output / f"{sample_name}.kraken"
    kraken_report = output / f"{sample_name}.kreport"

    # Check skip-if-exists
    if skip_if_exists and kraken_output.exists() and kraken_report.exists():
        out.print(f"[yellow]Skipping: Kraken2 output already exists for {sample_name}[/yellow]")
        out.print(f"  {kraken_output}")
        out.print(f"  {kraken_report}")
        raise typer.Exit(code=0)

    # Show parameters
    out.print("[bold]Input:[/bold]")
    out.print(f"  Forward reads: {reads_1}")
    if reads_2:
        out.print(f"  Reverse reads: {reads_2}")
    out.print(f"  Kraken2 DB:    {kraken_db}")

    out.print("\n[bold]Parameters:[/bold]")
    out.print(f"  Confidence:    {confidence}")
    out.print(f"  Threads:       {threads}")

    # Run Kraken2
    kraken = Kraken2()

    if dry_run:
        console.print("\n[bold cyan]DRY RUN MODE[/bold cyan]\n")
        result = kraken.run(
            reads_1=reads_1,
            reads_2=reads_2,
            database=kraken_db,
            output=kraken_output,
            report=kraken_report,
            threads=threads,
            confidence=confidence,
            dry_run=True,
        )
        console.print(f"[dim]Command: {result.command_string}[/dim]")
        console.print("\n[green]Dry run complete. No files were created.[/green]")
        raise typer.Exit(code=0)

    out.print("\n[bold]Running Kraken2 classification...[/bold]")

    with spinner_progress("Classifying reads...", console, quiet):
        try:
            result = kraken.run_or_raise(
                reads_1=reads_1,
                reads_2=reads_2,
                database=kraken_db,
                output=kraken_output,
                report=kraken_report,
                threads=threads,
                confidence=confidence,
            )
        except ToolExecutionError as e:
            console.print(f"\n[red]Kraken2 failed:[/red]\n{e.message}")
            raise typer.Exit(code=1) from None

    out.print(f"  [green]Classification complete ({result.elapsed_seconds:.1f}s)[/green]")

    if verbose:
        out.print(f"\n[dim]Command: {result.command_string}[/dim]")
        if result.stderr.strip():
            out.print(f"[dim]Tool output: {result.stderr.strip()[:300]}[/dim]")

    # Parse report for summary
    try:
        KrakenReport(kraken_report)  # Validates report format
        if verbose:
            out.print("[dim]Report parsed successfully[/dim]")
    except Exception as e:
        if verbose:
            out.print(f"\n[dim]Could not parse report: {e}[/dim]")

    # Summary
    out.print("\n[bold green]Classification complete![/bold green]")
    out.print("\n[bold]Output files:[/bold]")
    out.print(f"  Kraken output: {kraken_output}")
    out.print(f"  Kraken report: {kraken_report}")

    out.print(
        f"\n[dim]Next step: metadarkmatter kraken2 extract "
        f"--kraken-output {kraken_output} --taxid <TAXID> --output <dir>[/dim]"
    )
    out.print()


@app.command(name="extract")
def extract_reads(
    kraken_output: Path = typer.Option(
        ...,
        "--kraken-output", "-k",
        help="Kraken2 output file (.kraken) from 'kraken2 classify'",
        exists=True,
        dir_okay=False,
    ),
    reads_1: Path = typer.Option(
        ...,
        "--reads-1", "-1",
        help="Original forward reads file (FASTQ)",
        exists=True,
        dir_okay=False,
    ),
    reads_2: Path | None = typer.Option(
        None,
        "--reads-2", "-2",
        help="Original reverse reads file for paired-end data",
        exists=True,
        dir_okay=False,
    ),
    taxid: int = typer.Option(
        ...,
        "--taxid", "-t",
        help="Target taxonomic ID to extract (e.g., family TaxID)",
        min=1,
    ),
    output: Path = typer.Option(
        ...,
        "--output", "-o",
        help="Output directory for extracted reads",
    ),
    sample_name: str | None = typer.Option(
        None,
        "--sample-name", "-n",
        help="Sample name for output files (default: derived from kraken output)",
    ),
    include_children: bool = typer.Option(
        True,
        "--include-children/--exact-match",
        help="Include reads from child taxa (default) or exact match only",
    ),
    skip_if_exists: bool = typer.Option(
        False,
        "--skip-if-exists",
        help="Skip extraction if output files already exist",
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
    compress: bool = typer.Option(
        True,
        "--compress/--no-compress",
        help="Compress output files with gzip (default: compress)",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Show commands without executing",
    ),
) -> None:
    """
    Extract reads belonging to a specific taxon from Kraken2 output.

    This command uses KrakenTools to extract reads classified to a specific
    taxonomic ID (and optionally its children) from existing Kraken2 output.

    Run 'kraken2 classify' first to generate the required Kraken2 files.

    Example:

        # Extract all Francisellaceae reads (TaxID 262)
        metadarkmatter kraken2 extract \\
            --kraken-output ./kraken_output/sample.kraken \\
            --reads-1 sample_R1.fastq.gz \\
            --reads-2 sample_R2.fastq.gz \\
            --taxid 262 \\
            --output ./extracted/

        # Extract without child taxa (exact match only)
        metadarkmatter kraken2 extract \\
            --kraken-output ./kraken_output/sample.kraken \\
            --reads-1 sample_R1.fastq.gz \\
            --taxid 262 \\
            --output ./extracted/ \\
            --exact-match

    Output files:
        - {sample}_taxid{taxid}_R1.fastq.gz: Extracted forward reads (compressed)
        - {sample}_taxid{taxid}_R2.fastq.gz: Extracted reverse reads (if paired)

    Use --no-compress for uncompressed output if needed.
    """
    out = QuietConsole(console, quiet=quiet)

    out.print("\n[bold blue]Metadarkmatter Read Extraction[/bold blue]\n")

    # Check tool availability
    if not ExtractKrakenReads.check_available():
        console.print("[red]Error: extract_kraken_reads.py not found[/red]")
        console.print("\n[dim]Install with: conda install -c bioconda krakentools[/dim]")
        raise typer.Exit(code=1) from None

    # Derive sample name from kraken output if not provided
    if sample_name is None:
        sample_name = extract_sample_name(kraken_output)

    # Find corresponding kreport file
    kraken_report = kraken_output.with_suffix(".kreport")
    if not kraken_report.exists():
        # Try alternative naming
        kraken_report = kraken_output.parent / f"{sample_name}.kreport"
        if not kraken_report.exists():
            console.print("[red]Error: Kraken report not found[/red]")
            console.print(f"  Expected: {kraken_output.with_suffix('.kreport')}")
            console.print(f"  Or: {kraken_output.parent / f'{sample_name}.kreport'}")
            raise typer.Exit(code=1) from None

    # Setup output paths
    output.mkdir(parents=True, exist_ok=True)

    # Determine output filenames (compressed or not)
    suffix = ".fastq.gz" if compress else ".fastq"
    extracted_r1 = output / f"{sample_name}_taxid{taxid}_R1{suffix}"
    extracted_r2 = output / f"{sample_name}_taxid{taxid}_R2{suffix}" if reads_2 else None

    # Intermediate uncompressed files (used during extraction)
    temp_r1 = output / f"{sample_name}_taxid{taxid}_R1.fastq"
    temp_r2 = output / f"{sample_name}_taxid{taxid}_R2.fastq" if reads_2 else None

    # Check skip-if-exists
    if skip_if_exists:
        output_exists = extracted_r1.exists()
        if reads_2 and extracted_r2:
            output_exists = output_exists and extracted_r2.exists()

        if output_exists:
            out.print("[yellow]Skipping: Output files already exist[/yellow]")
            out.print(f"  {extracted_r1}")
            if extracted_r2:
                out.print(f"  {extracted_r2}")
            raise typer.Exit(code=0)

    # Show parameters
    out.print("[bold]Input:[/bold]")
    out.print(f"  Kraken output: {kraken_output}")
    out.print(f"  Kraken report: {kraken_report}")
    out.print(f"  Forward reads: {reads_1}")
    if reads_2:
        out.print(f"  Reverse reads: {reads_2}")

    out.print("\n[bold]Parameters:[/bold]")
    out.print(f"  Target TaxID:     {taxid}")
    out.print(f"  Include children: {include_children}")

    # Parse Kraken report to show expected reads
    if not dry_run:
        try:
            report = KrakenReport(kraken_report)
            clade_reads = report.get_clade_reads(taxid)
            out.print(f"  Expected reads:   {clade_reads:,}")
        except Exception as e:
            if verbose:
                out.print(f"  [dim]Could not estimate reads: {e}[/dim]")

    # Run extraction
    extractor = ExtractKrakenReads()

    if dry_run:
        console.print("\n[bold cyan]DRY RUN MODE[/bold cyan]\n")
        result = extractor.run(
            kraken_output=kraken_output,
            kraken_report=kraken_report,
            reads_1=reads_1,
            reads_2=reads_2,
            taxid=taxid,
            output_1=temp_r1,
            output_2=temp_r2,
            include_children=include_children,
            dry_run=True,
        )
        console.print(f"[dim]Command: {result.command_string}[/dim]")
        if compress:
            console.print(f"[dim]Then compress: gzip {temp_r1}[/dim]")
        console.print("\n[green]Dry run complete. No files were created.[/green]")
        raise typer.Exit(code=0)

    out.print(f"\n[bold]Extracting reads for TaxID {taxid}...[/bold]")

    with spinner_progress("Extracting reads...", console, quiet):
        try:
            result = extractor.run_or_raise(
                kraken_output=kraken_output,
                kraken_report=kraken_report,
                reads_1=reads_1,
                reads_2=reads_2,
                taxid=taxid,
                output_1=temp_r1,
                output_2=temp_r2,
                include_children=include_children,
            )
        except ToolExecutionError as e:
            console.print(f"\n[red]Extraction failed:[/red]\n{e.message}")
            raise typer.Exit(code=1) from None

    out.print(f"  [green]Extraction complete ({result.elapsed_seconds:.1f}s)[/green]")

    if verbose:
        out.print(f"\n[dim]Command: {result.command_string}[/dim]")
        if result.stderr.strip():
            out.print(f"[dim]Tool output: {result.stderr.strip()[:300]}[/dim]")

    # Compress output files if requested
    if compress:
        import gzip
        import shutil

        out.print("  [dim]Compressing output files...[/dim]")

        for temp_file, final_file in [(temp_r1, extracted_r1), (temp_r2, extracted_r2)]:
            if temp_file and temp_file.exists():
                with temp_file.open("rb") as f_in, gzip.open(final_file, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
                temp_file.unlink()  # Remove uncompressed file

        out.print("  [green]Compression complete[/green]")

    # Count extracted reads
    read_count = 0
    if extracted_r1.exists():
        # Quick line count (4 lines per FASTQ record)
        if compress:
            import gzip
            with gzip.open(extracted_r1, "rt") as f:
                line_count = sum(1 for _ in f)
        else:
            with extracted_r1.open() as f:
                line_count = sum(1 for _ in f)
        read_count = line_count // 4
        out.print(f"  [dim]Extracted reads: {read_count:,}[/dim]")

    # Summary
    out.print("\n[bold green]Extraction complete![/bold green]")
    out.print("\n[bold]Output files:[/bold]")
    out.print(f"  Extracted R1: {extracted_r1}")
    if extracted_r2:
        out.print(f"  Extracted R2: {extracted_r2}")

    out.print(
        f"\n[dim]Next step: metadarkmatter blast align --query {extracted_r1} "
        f"--database <blastdb> --output <output.tsv.gz>[/dim]"
    )
    out.print()
