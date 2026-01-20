"""
Read recruitment data extraction from BAM files.

Extracts alignment information for recruitment plot visualization
without requiring pysam as a dependency.
"""

from __future__ import annotations

import re
import subprocess
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from metadarkmatter.external.base import validate_path_safe


@dataclass(frozen=True)
class AlignmentRecord:
    """Single alignment record from BAM file."""

    read_id: str
    genome_name: str
    contig: str
    position: int
    mapq: int
    cigar: str
    percent_identity: float
    alignment_length: int


def extract_genome_name(contig: str) -> str:
    """Extract genome name from contig header.

    Assumes format: GenomeName_ContigID or GenomeName_scaffold_N
    """
    # Handle common genome prefixes
    if contig.startswith(("GCF_", "GCA_")):
        # NCBI accession format: GCF_000123456.1_ASM123v1_scaffold_1
        match = re.match(r"(GC[AF]_\d+\.\d+)", contig)
        if match:
            return match.group(1)

    # Generic format: take first part before underscore
    parts = contig.split("_")
    if len(parts) >= 2:
        # Check if second part looks like a version or number
        if parts[1].isdigit() or parts[1].startswith(("v", "scaffold", "contig")):
            return parts[0]
        # Otherwise might be GenomeName_Species format
        return "_".join(parts[:2])

    return contig


def calculate_percent_identity(cigar: str, md_tag: str | None = None) -> float:
    """Calculate percent identity from CIGAR string.

    This is an approximation based on CIGAR operations.
    More accurate calculation would use the MD tag if available.
    """
    # Parse CIGAR
    matches = 0
    total = 0

    for match in re.finditer(r"(\d+)([MIDNSHP=X])", cigar):
        length = int(match.group(1))
        op = match.group(2)

        if op in ("M", "="):
            matches += length
            total += length
        elif op == "X" or op in ("I", "D"):
            total += length
        # S, H, N, P don't affect identity calculation

    if total == 0:
        return 0.0

    return 100.0 * matches / total


def parse_sam_line(line: str) -> AlignmentRecord | None:
    """Parse a single SAM line into an AlignmentRecord."""
    if line.startswith("@"):
        return None

    fields = line.strip().split("\t")
    if len(fields) < 11:
        return None

    # Skip unmapped reads
    flag = int(fields[1])
    if flag & 4:  # Unmapped
        return None

    read_id = fields[0]
    contig = fields[2]
    position = int(fields[3])
    mapq = int(fields[4])
    cigar = fields[5]

    if cigar == "*":
        return None

    # Calculate alignment length from CIGAR
    alignment_length = 0
    for match in re.finditer(r"(\d+)([MIDNSHP=X])", cigar):
        length = int(match.group(1))
        op = match.group(2)
        if op in ("M", "=", "X", "I", "D"):
            alignment_length += length

    # Calculate percent identity
    percent_identity = calculate_percent_identity(cigar)

    # Extract genome name
    genome_name = extract_genome_name(contig)

    return AlignmentRecord(
        read_id=read_id,
        genome_name=genome_name,
        contig=contig,
        position=position,
        mapq=mapq,
        cigar=cigar,
        percent_identity=percent_identity,
        alignment_length=alignment_length,
    )


def stream_bam_alignments(
    bam_path: Path,
    min_mapq: int = 0,
    min_identity: float = 0.0,
    samtools_path: str = "samtools",
) -> Iterator[AlignmentRecord]:
    """Stream alignments from BAM file using samtools.

    Args:
        bam_path: Path to BAM file.
        min_mapq: Minimum mapping quality filter.
        min_identity: Minimum percent identity filter.
        samtools_path: Path to samtools executable.

    Yields:
        AlignmentRecord for each alignment passing filters.
    """
    # Validate path for safety
    bam_path = validate_path_safe(bam_path, must_exist=True)

    cmd = [samtools_path, "view", str(bam_path)]

    if min_mapq > 0:
        cmd.extend(["-q", str(min_mapq)])

    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if process.stdout is None:
        return

    for line in process.stdout:
        record = parse_sam_line(line)
        if record is None:
            continue

        if record.percent_identity >= min_identity:
            yield record

    process.wait()


def load_recruitment_data(
    bam_path: Path,
    min_mapq: int = 0,
    min_identity: float = 70.0,
    max_records: int | None = None,
    samtools_path: str = "samtools",
) -> pl.DataFrame:
    """Load recruitment data from BAM file into Polars DataFrame.

    Args:
        bam_path: Path to BAM file.
        min_mapq: Minimum mapping quality filter.
        min_identity: Minimum percent identity filter.
        max_records: Maximum records to load (for sampling).
        samtools_path: Path to samtools executable.

    Returns:
        Polars DataFrame with columns:
        - read_id: Read identifier
        - genome_name: Reference genome name
        - contig: Full contig name
        - position: Alignment position
        - mapq: Mapping quality
        - percent_identity: Alignment percent identity
        - alignment_length: Length of alignment
    """
    records = []

    for count, record in enumerate(
        stream_bam_alignments(
            bam_path,
            min_mapq=min_mapq,
            min_identity=min_identity,
            samtools_path=samtools_path,
        ),
        start=1,
    ):
        records.append({
            "read_id": record.read_id,
            "genome_name": record.genome_name,
            "contig": record.contig,
            "position": record.position,
            "mapq": record.mapq,
            "percent_identity": record.percent_identity,
            "alignment_length": record.alignment_length,
        })

        if max_records and count >= max_records:
            break

    return pl.DataFrame(records)


def get_contig_lengths(bam_path: Path, samtools_path: str = "samtools") -> dict[str, int]:
    """Extract contig lengths from BAM header.

    Args:
        bam_path: Path to BAM file.
        samtools_path: Path to samtools executable.

    Returns:
        Dictionary mapping contig name to length.
    """
    # Validate path for safety
    bam_path = validate_path_safe(bam_path, must_exist=True)

    cmd = [samtools_path, "view", "-H", str(bam_path)]

    result = subprocess.run(cmd, capture_output=True, text=True)

    lengths = {}
    for line in result.stdout.split("\n"):
        if line.startswith("@SQ"):
            parts = dict(p.split(":", 1) for p in line.split("\t")[1:] if ":" in p)
            if "SN" in parts and "LN" in parts:
                lengths[parts["SN"]] = int(parts["LN"])

    return lengths


def aggregate_by_genome(df: pl.DataFrame) -> pl.DataFrame:
    """Aggregate recruitment data by genome.

    Args:
        df: Recruitment DataFrame from load_recruitment_data.

    Returns:
        DataFrame with per-genome statistics:
        - genome_name: Genome identifier
        - num_reads: Number of aligned reads
        - mean_identity: Mean percent identity
        - median_identity: Median percent identity
        - min_identity: Minimum percent identity
        - max_identity: Maximum percent identity
    """
    return (
        df.group_by("genome_name")
        .agg(
            pl.len().alias("num_reads"),
            pl.col("percent_identity").mean().alias("mean_identity"),
            pl.col("percent_identity").median().alias("median_identity"),
            pl.col("percent_identity").min().alias("min_identity"),
            pl.col("percent_identity").max().alias("max_identity"),
        )
        .sort("num_reads", descending=True)
    )
