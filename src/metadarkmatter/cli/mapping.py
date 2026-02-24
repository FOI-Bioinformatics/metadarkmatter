"""
Utility commands for metadarkmatter.

Provides helper commands for data preparation and troubleshooting.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console

import polars as pl

from metadarkmatter.cli.utils import QuietConsole, spinner_progress
from metadarkmatter.core.genome_utils import extract_accession_from_filename
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


@app.command(name="select-representatives")
def select_representatives(
    metadata: Annotated[
        Path,
        typer.Option(
            "--metadata",
            "-m",
            help="Path to genome_metadata.tsv file",
            exists=True,
            dir_okay=False,
            resolve_path=True,
        ),
    ],
    output: Annotated[
        Path,
        typer.Option(
            "--output",
            "-o",
            help="Output path for updated metadata TSV",
            dir_okay=False,
            resolve_path=True,
        ),
    ],
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
    Assign or update species representative genomes in metadata.

    Groups genomes by species name and assigns a representative for each
    species. If a genome already has a self-referencing representative
    entry, it is kept. Otherwise, the first accession alphabetically is
    selected.

    This is useful when:
    - Metadata was created with --reps-only and more genomes were added later
    - Custom (non-GTDB) genomes need representative assignments
    - Representative assignments need to be regenerated

    Example:
        metadarkmatter util select-representatives \\
            --metadata genome_metadata.tsv \\
            --output genome_metadata_updated.tsv
    """
    qc = QuietConsole(console, quiet)

    qc.print("[bold]Selecting species representatives[/bold]")

    try:
        df = pl.read_csv(metadata, separator="\t")

        required = {"accession", "species"}
        missing = required - set(df.columns)
        if missing:
            console.print(f"[red]Error:[/red] Missing required columns: {missing}")
            raise typer.Exit(1)

        # Group genomes by species
        has_existing_reps = "representative" in df.columns

        # Build species -> representative mapping
        species_to_rep: dict[str, str] = {}

        if has_existing_reps:
            # Preserve existing self-referencing representatives
            for row in df.iter_rows(named=True):
                species = row["species"]
                accession = row["accession"]
                rep = row.get("representative", "")
                if species and rep == accession:
                    species_to_rep[species] = accession

        # For species without a representative, pick first alphabetically
        species_groups = df.group_by("species").agg(
            pl.col("accession").sort().first().alias("first_accession")
        )
        for row in species_groups.iter_rows(named=True):
            species = row["species"]
            if species and species not in species_to_rep:
                species_to_rep[species] = row["first_accession"]

        # Apply mapping to DataFrame
        representative_col = [
            species_to_rep.get(species, accession)
            for species, accession in zip(
                df["species"].to_list(),
                df["accession"].to_list(),
                strict=True,
            )
        ]

        if "representative" in df.columns:
            df = df.drop("representative")

        # Insert representative column after family (or after genus if no family)
        insert_after = "family" if "family" in df.columns else "genus"
        cols = df.columns
        insert_idx = cols.index(insert_after) + 1
        before_cols = cols[:insert_idx]
        after_cols = cols[insert_idx:]

        df = df.select(
            [pl.col(c) for c in before_cols]
            + [pl.Series("representative", representative_col)]
            + [pl.col(c) for c in after_cols]
        )

        output.parent.mkdir(parents=True, exist_ok=True)
        df.write_csv(output, separator="\t")

        num_species = len(species_to_rep)
        num_genomes = len(df)
        qc.print(
            f"\n[green]Success![/green] Assigned {num_species:,} representatives "
            f"for {num_genomes:,} genomes"
        )
        qc.print(f"[dim]Output:[/dim] {output}")

    except typer.Exit:
        raise
    except Exception as e:
        console.print(f"[red]Error:[/red] {e}", style="bold")
        raise typer.Exit(1) from e


def _parse_structured_filename(
    filename: str,
    family_override: str | None = None,
) -> dict[str, str]:
    """Parse taxonomy from a structured genome filename.

    Expected pattern: {Family}_{Genus}_{species}_{Accession}.{ext}

    The accession is located by regex search for GCF_/GCA_ patterns,
    and everything before it is split into family, genus, and species
    tokens.

    Args:
        filename: Genome filename (basename only)
        family_override: If set, use this instead of parsing family
            from the filename

    Returns:
        Dict with keys: accession, species, genus, family

    Raises:
        ValueError: If no GCF/GCA accession can be found in the filename
    """
    accession = extract_accession_from_filename(filename)

    # Check that we actually found a real accession (not the whole stem)
    if not re.match(r"GC[FA]_\d+\.\d+$", accession):
        raise ValueError(
            f"Cannot extract GCF/GCA accession from filename: {filename}"
        )

    # Strip extension to get stem
    stem = filename
    for suffix in [".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta"]:
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break

    # Find where the accession starts in the stem
    acc_match = re.search(r"GC[FA]_\d+\.\d+", stem)
    if not acc_match:
        raise ValueError(
            f"Cannot locate accession position in filename: {filename}"
        )

    prefix = stem[: acc_match.start()].rstrip("_")

    if not prefix:
        # No prefix before accession (standard NCBI filename)
        return {
            "accession": accession,
            "species": "",
            "genus": "",
            "family": family_override or "",
        }

    parts = prefix.split("_")

    if family_override:
        # All tokens are genus + species
        genus = parts[0] if parts else ""
        species_epithet = "_".join(parts[1:]) if len(parts) > 1 else ""
        family = family_override
    elif len(parts) >= 3:
        # Family_Genus_species...
        family = parts[0]
        genus = parts[1]
        species_epithet = "_".join(parts[2:])
    elif len(parts) == 2:
        # Genus_species (no family)
        family = ""
        genus = parts[0]
        species_epithet = parts[1]
    else:
        # Single token before accession
        family = ""
        genus = parts[0]
        species_epithet = ""

    # Construct full species name: "Genus epithet"
    if species_epithet:
        full_species = f"{genus} {species_epithet}"
    else:
        full_species = genus

    return {
        "accession": accession,
        "species": full_species,
        "genus": genus,
        "family": family,
    }


@app.command(name="generate-metadata")
def generate_metadata(
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
            help="Output path for genome_metadata.tsv",
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
    gtdb_metadata: Annotated[
        Path | None,
        typer.Option(
            "--gtdb-metadata",
            help="GTDB metadata file for representative assignment",
            exists=True,
            dir_okay=False,
            resolve_path=True,
        ),
    ] = None,
    family: Annotated[
        str | None,
        typer.Option(
            "--family",
            help="Override family name (instead of parsing from filename)",
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
    Generate genome_metadata.tsv from a directory of genome FASTA files.

    Parses structured filenames to extract taxonomy. Expected filename
    pattern: {Family}_{Genus}_{species}_{Accession}.{ext}

    Example filenames:
        Francisellaceae_Allofrancisella_frigidaquae_GCA_000710735.1.fasta
        Francisellaceae_Francisella_tularensis_GCF_000195955.2.fasta

    Representative assignment:
        Without --gtdb-metadata: first accession per species (alphabetically)
        With --gtdb-metadata: cross-references GTDB representatives

    Example:
        metadarkmatter util generate-metadata \\
            --genomes genomes/ --pattern "*.fasta" --output genome_metadata.tsv
    """
    qc = QuietConsole(console, quiet)

    qc.print(f"[bold]Generating metadata from:[/bold] {genomes}")

    try:
        genome_files = sorted(genomes.glob(pattern))
        if not genome_files:
            console.print(
                f"[red]Error:[/red] No files matching '{pattern}' "
                f"in {genomes}"
            )
            raise typer.Exit(1)

        # Parse each filename
        records: list[dict[str, str]] = []
        parse_errors: list[str] = []

        with spinner_progress(
            "Parsing genome filenames...", console, quiet
        ) as _progress:
            for gf in genome_files:
                try:
                    info = _parse_structured_filename(gf.name, family)
                    records.append(info)
                except ValueError as e:
                    parse_errors.append(str(e))

        if parse_errors:
            for err in parse_errors:
                console.print(f"[yellow]Warning:[/yellow] {err}")

        if not records:
            console.print(
                "[red]Error:[/red] No genomes could be parsed. "
                "Check filename format and glob pattern."
            )
            raise typer.Exit(1)

        df = pl.DataFrame(records)

        # Load GTDB metadata for representative cross-referencing
        gtdb_rep_map: dict[str, str] = {}
        gtdb_tax_map: dict[str, str] = {}
        if gtdb_metadata is not None:
            qc.print(
                f"[dim]Cross-referencing with GTDB metadata:[/dim] "
                f"{gtdb_metadata}"
            )
            gtdb_df = pl.read_csv(gtdb_metadata, separator="\t")

            if "representative" in gtdb_df.columns:
                for row in gtdb_df.iter_rows(named=True):
                    acc = row.get("accession", "")
                    rep = row.get("representative", "")
                    tax = row.get("gtdb_taxonomy", "")
                    if acc:
                        if rep:
                            gtdb_rep_map[acc] = rep
                        if tax:
                            gtdb_tax_map[acc] = tax

        # Assign representatives
        if gtdb_rep_map:
            # Use GTDB representative mapping
            representative_col = [
                gtdb_rep_map.get(row["accession"], row["accession"])
                for row in records
            ]
            gtdb_taxonomy_col = [
                gtdb_tax_map.get(row["accession"], "")
                for row in records
            ]
        else:
            # First accession per species alphabetically
            species_to_rep: dict[str, str] = {}
            species_groups = df.group_by("species").agg(
                pl.col("accession").sort().first().alias("first_accession")
            )
            for row in species_groups.iter_rows(named=True):
                sp = row["species"]
                if sp:
                    species_to_rep[sp] = row["first_accession"]

            representative_col = [
                species_to_rep.get(row["species"], row["accession"])
                for row in records
            ]
            gtdb_taxonomy_col = [
                gtdb_tax_map.get(row["accession"], "")
                for row in records
            ]

        # Build output DataFrame
        out_df = pl.DataFrame({
            "accession": df["accession"],
            "species": df["species"],
            "genus": df["genus"],
            "family": df["family"],
            "representative": representative_col,
            "gtdb_taxonomy": gtdb_taxonomy_col,
        })

        output.parent.mkdir(parents=True, exist_ok=True)
        out_df.write_csv(output, separator="\t")

        num_genomes = len(out_df)
        num_species = out_df["species"].n_unique()
        num_reps = len(
            set(r for r in representative_col if r in set(df["accession"]))
        )
        qc.print(
            f"\n[green]Success![/green] {num_genomes:,} genomes, "
            f"{num_species:,} species, {num_reps:,} representatives"
        )
        qc.print(f"[dim]Output:[/dim] {output}")

        if parse_errors:
            qc.print(
                f"\n[yellow]Warning:[/yellow] {len(parse_errors):,} files "
                "could not be parsed (see warnings above)"
            )

    except typer.Exit:
        raise
    except Exception as e:
        console.print(f"[red]Error:[/red] {e}", style="bold")
        raise typer.Exit(1) from e
