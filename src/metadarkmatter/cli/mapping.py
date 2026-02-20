"""
Utility commands for metadarkmatter.

Provides helper commands for data preparation and troubleshooting.
"""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console

from metadarkmatter.cli.utils import QuietConsole, spinner_progress
from metadarkmatter.core.id_mapping import ContigIdMapping

app = typer.Typer(
    name="util",
    help="Utility commands for data preparation",
    no_args_is_help=True,
)

console = Console()


@app.command(name="generate-mapping")
def generate_mapping(
    genomes: Annotated[
        Path,
        typer.Option(
            "--genomes",
            "-g",
            help="Directory containing genome FASTA files",
            exists=True,
            file_okay=False,
            resolve_path=True,
        ),
    ],
    output: Annotated[
        Path,
        typer.Option(
            "--output",
            "-o",
            help="Output TSV file for ID mapping",
            dir_okay=False,
            resolve_path=True,
        ),
    ],
    pattern: Annotated[
        str,
        typer.Option(
            "--pattern",
            "-p",
            help="Glob pattern for genome files",
        ),
    ] = "*.fna",
    quiet: Annotated[
        bool,
        typer.Option(
            "--quiet",
            "-q",
            help="Suppress progress output",
        ),
    ] = False,
) -> None:
    """
    Generate ID mapping from genome FASTA directory.

    Creates a TSV file mapping original contig IDs (e.g., NZ_CP007557.1)
    to their parent genome accessions (e.g., GCF_000195955.2).

    This mapping is used to transform external BLAST, Bowtie2, or Kraken2
    results to use metadarkmatter's standardized ID format.

    Example:
        metadarkmatter util generate-mapping --genomes genomes/ --output id_mapping.tsv
    """
    qc = QuietConsole(console, quiet)

    qc.print(f"[bold]Generating ID mapping from:[/bold] {genomes}")

    try:
        with spinner_progress(
            "Scanning genome files...", console, quiet
        ) as _progress:
            mapping = ContigIdMapping.from_genome_dir(genomes, pattern=pattern)

        # Save to file
        mapping.to_tsv(output)

        qc.print(f"\n[green]Success![/green] Mapped {len(mapping):,} contigs")
        qc.print(f"[dim]Output:[/dim] {output}")

    except FileNotFoundError as e:
        console.print(f"[red]Error:[/red] {e}", style="bold")
        raise typer.Exit(1) from e
    except ValueError as e:
        console.print(f"[red]Error:[/red] {e}", style="bold")
        raise typer.Exit(1) from e


@app.command(name="validate-mapping")
def validate_mapping(
    mapping_file: Annotated[
        Path,
        typer.Argument(
            help="Path to ID mapping TSV file",
            exists=True,
            resolve_path=True,
        ),
    ],
    blast_file: Annotated[
        Path | None,
        typer.Option(
            "--blast",
            "-b",
            help="Optional BLAST TSV to check coverage",
            exists=True,
            resolve_path=True,
        ),
    ] = None,
    quiet: Annotated[
        bool,
        typer.Option(
            "--quiet",
            "-q",
            help="Suppress progress output",
        ),
    ] = False,
) -> None:
    """
    Validate an ID mapping file.

    Checks that the mapping file is properly formatted and optionally
    validates coverage against a BLAST results file.

    Example:
        metadarkmatter util validate-mapping id_mapping.tsv --blast results.blast.tsv
    """
    import polars as pl

    qc = QuietConsole(console, quiet)

    try:
        mapping = ContigIdMapping.from_tsv(mapping_file)
        qc.print(f"[green]Valid mapping file[/green] with {len(mapping):,} entries")

        # Show sample entries
        sample_entries = list(mapping.contig_to_accession.items())[:5]
        qc.print("\n[bold]Sample entries:[/bold]")
        for contig, accession in sample_entries:
            qc.print(f"  {contig} -> {accession}")

        # Validate against BLAST file if provided
        if blast_file:
            qc.print(f"\n[bold]Checking coverage against:[/bold] {blast_file}")

            # Read unique sseqids from BLAST file
            df = pl.scan_csv(
                blast_file,
                separator="\t",
                has_header=False,
            ).select(pl.col("column_2").alias("sseqid")).unique().collect()

            unique_ids = set(df["sseqid"].to_list())
            mapped = sum(1 for sid in unique_ids if sid in mapping)
            unmapped = unique_ids - set(mapping.contig_to_accession.keys())

            coverage = (mapped / len(unique_ids)) * 100 if unique_ids else 0

            qc.print(f"  Unique subject IDs: {len(unique_ids):,}")
            qc.print(f"  Mapped: {mapped:,} ({coverage:.1f}%)")
            qc.print(f"  Unmapped: {len(unmapped):,}")

            if unmapped and len(unmapped) <= 10:
                qc.print("\n[bold]Unmapped IDs:[/bold]")
                for uid in sorted(unmapped)[:10]:
                    qc.print(f"  {uid}")

            if coverage < 95:
                console.print(
                    f"\n[yellow]Warning:[/yellow] Low coverage ({coverage:.1f}%). "
                    "Some IDs may not be transformed correctly."
                )

    except Exception as e:
        console.print(f"[red]Error:[/red] {e}", style="bold")
        raise typer.Exit(1) from e
