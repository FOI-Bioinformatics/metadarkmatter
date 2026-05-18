"""
Validate input files for the metadarkmatter pipeline.

Provides subcommands for checking BLAST output, classification results,
and ANI matrices before running downstream analyses.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import polars as pl
import typer
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from metadarkmatter.core.sequence_validation import (
    ValidationResult,
    validate_fasta,
    validate_fastq,
)

app = typer.Typer(
    name="validate",
    help="Validate input files for the metadarkmatter pipeline",
    no_args_is_help=True,
)

console = Console()

BLAST_EXPECTED_COLUMNS = 12
BLAST_COLUMN_NAMES = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]

REQUIRED_CLASSIFICATION_COLUMNS = [
    "read_id",
    "best_match_genome",
    "novelty_index",
    "placement_uncertainty",
    "taxonomic_call",
]

VALID_TAXONOMIC_CALLS = {
    "Known Species",
    "Novel Species",
    "Novel Genus",
    "Ambiguous",
    "Conserved Region Hit",
}


def _read_lines_sample(
    path: Path, max_lines: int = 1000
) -> list[str]:
    """Read up to max_lines from a file, handling gzip transparently."""
    lines: list[str] = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i >= max_lines:
                break
            lines.append(line.rstrip("\n"))
    return lines


@app.command()
def blast(
    input: Path = typer.Option(
        ..., "--input", "-i", help="BLAST tabular output file (.tsv or .tsv.gz)"
    ),
) -> None:
    """Validate BLAST tabular output format (outfmt 6)."""
    if not input.exists():
        console.print(f"[red]File not found:[/red] {input}")
        raise typer.Exit(code=1)

    issues: list[str] = []

    # File size
    file_size_mb = input.stat().st_size / (1024 * 1024)
    console.print(f"File: {input}")
    console.print(f"Size: {file_size_mb:.2f} MB")

    # Sample lines
    lines = _read_lines_sample(input, max_lines=1000)
    if not lines:
        console.print("[red]File is empty[/red]")
        raise typer.Exit(code=1)

    console.print(f"Lines sampled: {len(lines)}")

    # Check column count on first 10 lines
    for i, line in enumerate(lines[:10]):
        cols = line.split("\t")
        if len(cols) != BLAST_EXPECTED_COLUMNS:
            issues.append(
                f"Line {i + 1}: expected {BLAST_EXPECTED_COLUMNS} columns, "
                f"found {len(cols)}"
            )

    # Check numeric columns (columns 3-12, 0-indexed 2-11)
    for i, line in enumerate(lines[:10]):
        cols = line.split("\t")
        if len(cols) != BLAST_EXPECTED_COLUMNS:
            continue
        for col_idx in range(2, BLAST_EXPECTED_COLUMNS):
            try:
                float(cols[col_idx])
            except ValueError:
                issues.append(
                    f"Line {i + 1}, column {col_idx + 1} ({BLAST_COLUMN_NAMES[col_idx]}): "
                    f"expected numeric, got '{cols[col_idx]}'"
                )

    # Try reading with Polars
    try:
        df = pl.read_csv(
            input,
            separator="\t",
            has_header=False,
            n_rows=100,
            new_columns=BLAST_COLUMN_NAMES,
        )
        schema_table = Table(title="Detected Schema")
        schema_table.add_column("Column")
        schema_table.add_column("Type")
        for col_name in df.columns:
            schema_table.add_row(col_name, str(df[col_name].dtype))
        console.print(schema_table)
    except Exception as e:
        issues.append(f"Polars schema error: {e}")

    # Report
    if issues:
        console.print(Panel("\n".join(issues), title="Issues Found", border_style="red"))
        raise typer.Exit(code=1)
    else:
        console.print("[green]BLAST file is valid[/green]")


@app.command()
def classifications(
    input: Path = typer.Option(
        ..., "--input", "-i", help="Classification output CSV"
    ),
) -> None:
    """Validate classification output file."""
    if not input.exists():
        console.print(f"[red]File not found:[/red] {input}")
        raise typer.Exit(code=1)

    issues: list[str] = []

    try:
        df = pl.read_csv(input)
    except Exception as e:
        console.print(f"[red]Failed to read CSV:[/red] {e}")
        raise typer.Exit(code=1)

    console.print(f"File: {input}")
    console.print(f"Total reads: {len(df)}")

    # Check required columns
    missing = [c for c in REQUIRED_CLASSIFICATION_COLUMNS if c not in df.columns]
    if missing:
        issues.append(f"Missing required columns: {', '.join(missing)}")

    # Check taxonomic_call values
    if "taxonomic_call" in df.columns:
        unique_calls = set(df["taxonomic_call"].unique().to_list())
        invalid = unique_calls - VALID_TAXONOMIC_CALLS
        if invalid:
            issues.append(f"Invalid taxonomic_call values: {', '.join(sorted(invalid))}")

        # Distribution table
        dist = df.group_by("taxonomic_call").len().sort("len", descending=True)
        dist_table = Table(title="Classification Distribution")
        dist_table.add_column("Category")
        dist_table.add_column("Count", justify="right")
        dist_table.add_column("Percent", justify="right")
        for row in dist.iter_rows():
            pct = 100.0 * row[1] / len(df)
            dist_table.add_row(str(row[0]), str(row[1]), f"{pct:.1f}%")
        console.print(dist_table)

    # Report
    if issues:
        console.print(Panel("\n".join(issues), title="Issues Found", border_style="red"))
        raise typer.Exit(code=1)
    else:
        console.print("[green]Classification file is valid[/green]")


@app.command()
def ani(
    input: Path = typer.Option(
        ..., "--input", "-i", help="ANI matrix CSV"
    ),
) -> None:
    """Validate ANI matrix format and values."""
    if not input.exists():
        console.print(f"[red]File not found:[/red] {input}")
        raise typer.Exit(code=1)

    issues: list[str] = []

    try:
        df = pl.read_csv(input)
    except Exception as e:
        console.print(f"[red]Failed to read CSV:[/red] {e}")
        raise typer.Exit(code=1)

    # First column is genome labels; remaining columns are the matrix
    genome_labels = df.columns[1:]
    n_rows = len(df)
    n_cols = len(genome_labels)

    console.print(f"File: {input}")
    console.print(f"Genomes: {n_rows}")

    # Check square
    if n_rows != n_cols:
        issues.append(f"Matrix is not square: {n_rows} rows x {n_cols} columns")

    # Extract numeric matrix
    try:
        matrix = df.select(genome_labels).to_numpy()
    except Exception as e:
        issues.append(f"Could not extract numeric matrix: {e}")
        console.print(Panel("\n".join(issues), title="Issues Found", border_style="red"))
        raise typer.Exit(code=1)

    # Value range
    min_val = float(matrix.min())
    max_val = float(matrix.max())
    console.print(f"Value range: {min_val:.2f} - {max_val:.2f}")

    if min_val < 0 or max_val > 100:
        issues.append(
            f"Values outside 0-100 range: min={min_val:.2f}, max={max_val:.2f}"
        )

    # Diagonal check (should be ~100)
    n_check = min(n_rows, n_cols)
    for i in range(n_check):
        diag_val = float(matrix[i, i])
        if abs(diag_val - 100.0) > 1.0:
            issues.append(
                f"Diagonal [{i}][{i}] ({genome_labels[i]}): "
                f"expected ~100, got {diag_val:.2f}"
            )

    # Symmetry check
    if n_rows == n_cols:
        max_asym = 0.0
        asym_pair = ("", "")
        for i in range(n_rows):
            for j in range(i + 1, n_cols):
                diff = abs(float(matrix[i, j]) - float(matrix[j, i]))
                if diff > max_asym:
                    max_asym = diff
                    asym_pair = (genome_labels[i], genome_labels[j])
        console.print(f"Max asymmetry: {max_asym:.4f}")
        if max_asym > 1.0:
            issues.append(
                f"Asymmetry exceeds 1.0: {max_asym:.4f} "
                f"between {asym_pair[0]} and {asym_pair[1]}"
            )

    # Summary table
    summary = Table(title="ANI Matrix Summary")
    summary.add_column("Metric")
    summary.add_column("Value", justify="right")
    summary.add_row("Genomes", str(n_rows))
    summary.add_row("Min ANI", f"{min_val:.2f}")
    summary.add_row("Max ANI", f"{max_val:.2f}")
    summary.add_row("Mean ANI", f"{float(matrix.mean()):.2f}")
    console.print(summary)

    # Report
    if issues:
        console.print(Panel("\n".join(issues), title="Issues Found", border_style="red"))
        raise typer.Exit(code=1)
    else:
        console.print("[green]ANI matrix is valid[/green]")


def _report_sequence_result(result: ValidationResult, kind: str) -> None:
    """Render a ValidationResult to the console and exit on failure."""
    console.print(f"File: {result.path}")
    console.print(f"Records sampled: {result.record_count}")
    if result.sample_record_id:
        console.print(f"First record id: {result.sample_record_id}")
    if result.issues:
        console.print(
            Panel(
                "\n".join(result.issues),
                title=f"{kind} issues found",
                border_style="red",
            )
        )
        raise typer.Exit(code=1)
    console.print(f"[green]{kind} file is valid[/green]")


@app.command()
def fasta(
    input: Path = typer.Option(
        ..., "--input", "-i", help="FASTA file (.fa/.fasta/.fna, optionally .gz)"
    ),
    max_records: int = typer.Option(
        0,
        "--max-records",
        help="Stop after this many records (0 = read whole file).",
        min=0,
    ),
    sequence_type: str = typer.Option(
        "nucleotide",
        "--sequence-type",
        "-t",
        help="Allowed-character set: 'nucleotide' or 'protein'.",
    ),
) -> None:
    """Validate a FASTA file (streaming, dependency-free)."""
    if sequence_type not in ("nucleotide", "protein"):
        console.print(
            f"[red]--sequence-type must be 'nucleotide' or 'protein', "
            f"got '{sequence_type}'[/red]"
        )
        raise typer.Exit(code=2)

    cap = max_records if max_records > 0 else None
    result = validate_fasta(input, max_records=cap, sequence_type=sequence_type)
    _report_sequence_result(result, "FASTA")


@app.command()
def fastq(
    input: Path = typer.Option(
        ..., "--input", "-i", help="FASTQ file (.fq/.fastq, optionally .gz)"
    ),
    max_records: int = typer.Option(
        0,
        "--max-records",
        help="Stop after this many records (0 = read whole file).",
        min=0,
    ),
) -> None:
    """Validate a FASTQ file (4-line record structure)."""
    cap = max_records if max_records > 0 else None
    result = validate_fastq(input, max_records=cap)
    _report_sequence_result(result, "FASTQ")
