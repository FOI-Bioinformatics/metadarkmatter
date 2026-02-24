"""
Download commands for GTDB genomes.

Provides commands to query GTDB for representative genome accessions
and download genome FASTAs from NCBI.
"""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.table import Table

from metadarkmatter.cli.utils import QuietConsole, spinner_progress
from metadarkmatter.clients.gtdb import (
    GTDBAPIError,
    GTDBClient,
    InvalidTaxonFormatError,
)
from metadarkmatter.external import NCBIDatasets, ToolNotFoundError
from metadarkmatter.models.genomes import AccessionList, GenomeAccession

app = typer.Typer(
    name="download",
    help="Download reference genomes from GTDB/NCBI",
    no_args_is_help=True,
)

genomes_app = typer.Typer(
    name="genomes",
    help="Download representative genomes from GTDB taxonomy",
    no_args_is_help=True,
)
app.add_typer(genomes_app, name="genomes")

console = Console()


@genomes_app.command(name="list")
def list_genomes(
    taxon: str = typer.Argument(
        ...,
        help="GTDB taxon (e.g., f__Enterobacteriaceae, g__Escherichia)",
    ),
    output: Path = typer.Option(
        ...,
        "--output", "-o",
        help="Output TSV file for accession list",
    ),
    representatives_only: bool = typer.Option(
        True,
        "--reps-only/--all-genomes",
        help="Only include representative genomes (default) or all genomes",
    ),
    include_metadata: bool = typer.Option(
        False,
        "--include-metadata",
        help="Include additional metadata columns (genome_size, etc.)",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose", "-v",
        help="Show detailed breakdown by genus",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet", "-q",
        help="Suppress progress output",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Show what would be queried without making requests",
    ),
) -> None:
    """
    Query GTDB API for representative genome accessions.

    Retrieves a list of genome accessions from the Genome Taxonomy Database
    (GTDB) for a specified taxonomic group. Results are saved to a TSV file
    that can be used with the 'fetch' command to download genome FASTAs.

    Taxon format uses GTDB prefixes:
        d__ = domain (e.g., d__Bacteria)
        p__ = phylum (e.g., p__Proteobacteria)
        c__ = class (e.g., c__Gammaproteobacteria)
        o__ = order (e.g., o__Enterobacterales)
        f__ = family (e.g., f__Enterobacteriaceae)
        g__ = genus (e.g., g__Escherichia)
        s__ = species (e.g., s__Escherichia coli)

    Example:

        # Get representative genomes for a bacterial family
        metadarkmatter download genomes list f__Enterobacteriaceae \\
            --output enterobacteriaceae_genomes.tsv

        # Get all genomes (not just representatives) for a genus
        metadarkmatter download genomes list g__Escherichia \\
            --output escherichia_all.tsv \\
            --all-genomes

        # Include metadata and show detailed breakdown
        metadarkmatter download genomes list f__Vibrionaceae \\
            --output vibrionaceae_genomes.tsv \\
            --include-metadata \\
            --verbose
    """
    out = QuietConsole(console, quiet=quiet)

    out.print("\n[bold blue]GTDB Genome Query[/bold blue]\n")
    out.print(f"Querying GTDB API for [bold]{taxon}[/bold]...")

    if dry_run:
        console.print(f"\n[dim]Would query: https://gtdb-api.ecogenomic.org/taxon/{taxon}/genomes[/dim]")
        console.print(f"[dim]Parameters: sp_reps_only={str(representatives_only).lower()}[/dim]")
        console.print("\n[green]Dry run complete. No requests were made.[/green]")
        raise typer.Exit(code=0)

    # Query GTDB API
    client = GTDBClient()

    try:
        with spinner_progress("Querying GTDB API...", console, quiet):
            result = client.query_genomes(taxon, representatives_only=representatives_only)

        # When downloading all genomes, also query for representatives
        # to build species-to-representative mapping
        reps_result = None
        if not representatives_only:
            with spinner_progress("Querying GTDB for species representatives...", console, quiet):
                reps_result = client.query_genomes(taxon, representatives_only=True)
    except InvalidTaxonFormatError as e:
        console.print(f"\n[red]Error: {e.message}[/red]")
        console.print(f"\n[dim]{e.suggestion}[/dim]")
        raise typer.Exit(code=1) from None
    except GTDBAPIError as e:
        console.print(f"\n[red]Error: {e.message}[/red]")
        console.print(f"\n[dim]{e.suggestion}[/dim]")
        raise typer.Exit(code=1) from None
    finally:
        client.close()

    if result.total_count == 0:
        console.print(f"\n[yellow]No genomes found for taxon '{taxon}'[/yellow]")
        console.print("[dim]Check that the taxon exists in GTDB and is spelled correctly.[/dim]")
        raise typer.Exit(code=1) from None

    # Display summary
    out.print("\n[bold]Summary:[/bold]")
    mode = "representative" if representatives_only else "all"
    out.print(
        f"  Found [green]{result.total_count:,}[/green] {mode} genomes "
        f"in [bold]{taxon}[/bold]"
    )
    out.print(f"  - {len(result.genus_counts):,} genera")
    out.print(f"  - {len(result.species_counts):,} species")

    # Show top genera in verbose mode
    if verbose and result.genus_counts:
        out.print()
        table = Table(title="Genomes by Genus (top 15)")
        table.add_column("Genus", style="cyan")
        table.add_column("Count", justify="right")

        sorted_genera = sorted(
            result.genus_counts.items(),
            key=lambda x: -x[1]
        )[:15]

        for genus, count in sorted_genera:
            table.add_row(genus, str(count))

        if len(result.genus_counts) > 15:
            table.add_row(
                f"... and {len(result.genus_counts) - 15} more",
                "",
                style="dim"
            )

        out.console.print(table)

    # Build representative mapping from reps query
    representative_map: dict[str, str] = {}
    rep_accession_set: set[str] = set()
    if reps_result is not None:
        for g in reps_result.genomes:
            if g.species:
                representative_map[g.species] = g.accession
                rep_accession_set.add(g.accession)
        out.print(
            f"  Identified [green]{len(representative_map):,}[/green] "
            f"species representatives"
        )

    # Convert to AccessionList
    accessions = [
        GenomeAccession(
            accession=g.accession,
            gtdb_taxonomy=g.gtdb_taxonomy,
            species=g.species,
            genome_size=g.genome_size,
            is_representative=g.accession in rep_accession_set if rep_accession_set else True,
        )
        for g in result.genomes
    ]

    accession_list = AccessionList(
        taxon=taxon,
        accessions=accessions,
        genus_counts=result.genus_counts,
        representative_map=representative_map,
    )

    # Save to TSV
    output.parent.mkdir(parents=True, exist_ok=True)
    accession_list.to_tsv(output, include_metadata=include_metadata)

    # Save genome metadata for pipeline use
    metadata_path = output.parent / "genome_metadata.tsv"
    accession_list.to_metadata_tsv(metadata_path)

    out.print("\n[bold green]Query complete![/bold green]")
    out.print(f"\nSaved accession list to: [bold]{output}[/bold]")
    out.print(f"Saved genome metadata to: [bold]{metadata_path}[/bold]")
    out.print()


@genomes_app.command(name="fetch")
def fetch_genomes(
    accessions: Path = typer.Option(
        ...,
        "--accessions", "-a",
        help="TSV file with accession list (from 'list' command)",
        exists=True,
        dir_okay=False,
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output-dir", "-o",
        help="Output directory for genome FASTA files",
    ),
    include_protein: bool = typer.Option(
        False,
        "--include-protein",
        help="Also download protein FASTA files (.faa) for AAI computation",
    ),
    skip_if_exists: bool = typer.Option(
        True,
        "--skip-if-exists/--redownload",
        help="Skip already downloaded genomes (default) or redownload all",
    ),
    decompress: bool = typer.Option(
        True,
        "--decompress/--keep-compressed",
        help="Decompress downloaded .gz files (default) or keep compressed",
    ),
    batch_size: int = typer.Option(
        100,
        "--batch-size",
        help="Number of accessions per NCBI request (max 1000)",
        min=1,
        max=1000,
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose", "-v",
        help="Show detailed progress",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet", "-q",
        help="Suppress progress output",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Show what would be downloaded without making requests",
    ),
) -> None:
    """
    Download genome FASTAs from NCBI using accession list.

    Uses the NCBI Datasets CLI to download genome assemblies from NCBI
    based on a list of accessions generated by the 'list' command.

    The NCBI Datasets CLI must be installed separately:
        conda install -c conda-forge ncbi-datasets-cli

    Example:

        # Download genomes from an accession list
        metadarkmatter download genomes fetch \\
            --accessions enterobacteriaceae_genomes.tsv \\
            --output-dir ./reference_genomes/

        # Download genomes with protein files for AAI computation
        metadarkmatter download genomes fetch \\
            --accessions enterobacteriaceae_genomes.tsv \\
            --output-dir ./reference_genomes/ \\
            --include-protein

        # Force redownload of all genomes
        metadarkmatter download genomes fetch \\
            --accessions enterobacteriaceae_genomes.tsv \\
            --output-dir ./reference_genomes/ \\
            --redownload

        # Keep compressed files (faster for large datasets)
        metadarkmatter download genomes fetch \\
            --accessions enterobacteriaceae_genomes.tsv \\
            --output-dir ./reference_genomes/ \\
            --keep-compressed
    """
    out = QuietConsole(console, quiet=quiet)

    out.print("\n[bold blue]NCBI Genome Download[/bold blue]\n")

    # Check NCBI Datasets CLI availability
    ncbi = NCBIDatasets()
    if not ncbi.check_available():
        console.print("[red]Error: NCBI Datasets CLI not found[/red]")
        console.print("\n[dim]Install with: conda install -c conda-forge ncbi-datasets-cli[/dim]")
        raise typer.Exit(code=1) from None

    # Load accession list
    out.print(f"Reading accession list: [bold]{accessions}[/bold]")

    try:
        accession_list = AccessionList.from_tsv(accessions)
    except Exception as e:
        console.print(f"[red]Error reading accession file: {e}[/red]")
        raise typer.Exit(code=1) from None

    out.print(f"  {accession_list.total_count:,} accessions to process")

    if accession_list.taxon:
        out.print(f"  Taxon: {accession_list.taxon}")

    if accession_list.total_count == 0:
        console.print("[yellow]No accessions found in file[/yellow]")
        raise typer.Exit(code=1) from None

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check for existing files if skip_if_exists
    accessions_to_download = accession_list.get_accession_strings()

    if skip_if_exists:
        existing_count = 0
        filtered_accessions = []

        for acc in accessions_to_download:
            fasta_path = output_dir / f"{acc}.fna"
            fasta_gz_path = output_dir / f"{acc}.fna.gz"
            if fasta_path.exists() or fasta_gz_path.exists():
                existing_count += 1
            else:
                filtered_accessions.append(acc)

        if existing_count > 0:
            out.print(f"  [dim]Skipping {existing_count:,} already downloaded genomes[/dim]")

        accessions_to_download = filtered_accessions

    if not accessions_to_download:
        out.print(f"\n[green]All genomes already exist in {output_dir}[/green]")
        raise typer.Exit(code=0)

    out.print(f"  {len(accessions_to_download):,} genomes to download")

    # Dry run mode
    if dry_run:
        console.print("\n[bold cyan]DRY RUN MODE[/bold cyan]")
        console.print(
            f"\n[dim]Would download {len(accessions_to_download):,} genomes "
            f"to {output_dir}[/dim]"
        )
        console.print(f"[dim]Batch size: {batch_size}[/dim]")
        console.print(f"[dim]Decompress: {decompress}[/dim]")
        console.print(f"[dim]Include protein: {include_protein}[/dim]")

        # Show sample of accessions
        sample = accessions_to_download[:5]
        console.print("\n[dim]Sample accessions:[/dim]")
        for acc in sample:
            console.print(f"  [dim]{acc}[/dim]")
        if len(accessions_to_download) > 5:
            console.print(f"  [dim]... and {len(accessions_to_download) - 5} more[/dim]")

        console.print("\n[green]Dry run complete. No files were downloaded.[/green]")
        raise typer.Exit(code=0)

    # Download genomes with progress
    out.print("\n[bold]Downloading genomes from NCBI...[/bold]")

    total = len(accessions_to_download)

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        TimeElapsedColumn(),
        console=console if not quiet else None,
    ) as progress:
        task = progress.add_task("Downloading...", total=total)

        def update_progress(completed: int, _total: int) -> None:
            progress.update(task, completed=completed)

        try:
            report = ncbi.download_genomes(
                accessions=accessions_to_download,
                output_dir=output_dir,
                batch_size=batch_size,
                decompress=decompress,
                skip_if_exists=False,  # Already filtered above
                include_protein=include_protein,
                retry_versions=True,
                progress_callback=update_progress,
            )
        except ToolNotFoundError as e:
            console.print(f"\n[red]Error: {e.message}[/red]")
            console.print(f"\n[dim]{e.suggestion}[/dim]")
            raise typer.Exit(code=1) from None
        except Exception as e:
            console.print(f"\n[red]Download failed: {e}[/red]")
            raise typer.Exit(code=1) from None

    # Count downloaded files
    genome_files = list(output_dir.glob("*.fna")) + list(output_dir.glob("*.fna.gz"))
    protein_files = list(output_dir.glob("*.faa")) + list(output_dir.glob("*.faa.gz"))
    all_downloaded = genome_files + protein_files
    total_size = sum(f.stat().st_size for f in all_downloaded)

    # Format size
    if total_size >= 1024 ** 3:
        size_str = f"{total_size / (1024 ** 3):.1f} GB"
    elif total_size >= 1024 ** 2:
        size_str = f"{total_size / (1024 ** 2):.1f} MB"
    else:
        size_str = f"{total_size / 1024:.1f} KB"

    # Summary from structured report
    n_succeeded = len(report.succeeded)
    n_failed = len(report.failed)
    n_recovered = len(report.recovered)

    out.print("\n[bold green]Download complete![/bold green]")
    out.print("\n[bold]Summary:[/bold]")
    out.print(f"  Requested: {len(accessions_to_download):,} genomes")
    out.print(f"  Downloaded: {n_succeeded:,} genomes")
    if n_recovered > 0:
        out.print(
            f"  [green]Recovered via version retry: {n_recovered:,} genomes[/green]"
        )
    out.print(f"  Genome files (.fna): {len(genome_files):,}")
    if n_failed > 0:
        out.print(
            f"  [yellow]Unavailable from NCBI: {n_failed:,} "
            f"(suppressed or removed assemblies)[/yellow]"
        )

    # Show version-recovered accessions
    if n_recovered > 0 and verbose:
        out.print("\n[dim]Recovered via alternative assembly versions:[/dim]")
        for o in report.recovered[:15]:
            out.print(f"[dim]  {o.accession} -> {o.resolved_version}[/dim]")
        if n_recovered > 15:
            out.print(f"[dim]  ... and {n_recovered - 15} more[/dim]")

    if include_protein:
        out.print(f"  Protein files (.faa): {len(protein_files):,}")

        # Warn if protein files are missing
        if len(protein_files) < len(genome_files):
            missing_count = len(genome_files) - len(protein_files)

            # Identify which genomes are missing proteins
            genome_accessions = {f.stem for f in genome_files}
            protein_accessions = {f.stem for f in protein_files}
            missing_accessions = sorted(genome_accessions - protein_accessions)

            out.print(f"\n[yellow]Warning: {missing_count} genomes missing protein files[/yellow]")
            out.print("[dim]This is common for GenBank-only assemblies or draft genomes.[/dim]")

            # Try to load metadata to show species information
            metadata_path = output_dir.parent / "genome_metadata.tsv"
            if not metadata_path.exists():
                metadata_path = output_dir / "genome_metadata.tsv"

            if metadata_path.exists() and missing_accessions:
                try:
                    import polars as pl
                    metadata = pl.read_csv(metadata_path, separator="\t")

                    missing_info = metadata.filter(
                        pl.col("accession").is_in(missing_accessions)
                    ).select(["accession", "species"])

                    if len(missing_info) > 0:
                        out.print("\n[dim]Missing protein files for:[/dim]")
                        for row in missing_info.head(10).iter_rows(named=True):
                            acc = row["accession"]
                            species = row["species"]
                            out.print(f"[dim]  - {acc:<20} {species}[/dim]")
                        if len(missing_info) > 10:
                            out.print(f"[dim]  ... and {len(missing_info) - 10} more[/dim]")
                except Exception:
                    if verbose and missing_accessions:
                        out.print("\n[dim]Missing protein files for:[/dim]")
                        for acc in missing_accessions[:10]:
                            out.print(f"[dim]  - {acc}[/dim]")
                        if len(missing_accessions) > 10:
                            out.print(f"[dim]  ... and {len(missing_accessions) - 10} more[/dim]")
            elif verbose and missing_accessions:
                out.print("\n[dim]Missing protein files for:[/dim]")
                for acc in missing_accessions[:10]:
                    out.print(f"[dim]  - {acc}[/dim]")
                if len(missing_accessions) > 10:
                    out.print(f"[dim]  ... and {len(missing_accessions) - 10} more[/dim]")

            out.print("\n[dim]To generate missing protein files with Prodigal:[/dim]")
            out.print(f"[dim]  cd {output_dir}[/dim]")
            out.print("[dim]  for fna in *.fna; do[/dim]")
            out.print("[dim]    base=$(basename $fna .fna)[/dim]")
            out.print("[dim]    [ ! -f ${base}.faa ] && prodigal -i $fna -a ${base}.faa -q[/dim]")
            out.print("[dim]  done[/dim]")

            if not verbose and not metadata_path.exists():
                out.print("\n[dim]Use --verbose to see which genomes are affected.[/dim]")

    out.print(f"  Output directory: [bold]{output_dir}[/bold]")
    out.print(f"  Total size: {size_str}")

    # Write failure report
    if n_failed > 0:
        failures_path = output_dir / "download_failures.tsv"
        report.write_failures_tsv(failures_path)
        out.print(f"\n  Failed accessions written to: [bold]{failures_path}[/bold]")

    if verbose:
        out.print(f"\n[dim]Elapsed time: {report.elapsed_seconds:.1f}s[/dim]")

    out.print()
