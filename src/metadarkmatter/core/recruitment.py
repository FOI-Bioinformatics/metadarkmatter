"""
Read recruitment data extraction from BAM files.

Extracts alignment information for recruitment plot visualization
without requiring pysam as a dependency.
"""

from __future__ import annotations

import re
import subprocess
import tempfile
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


def _cigar_counts(cigar: str) -> tuple[int, int, int, int, int]:
    """Return (M_or_eq_len, mismatch_X_len, ins_len, del_len, eq_len) from a CIGAR."""
    m = x = ins = dele = eq = 0
    for match in re.finditer(r"(\d+)([MIDNSHP=X])", cigar):
        length = int(match.group(1))
        op = match.group(2)
        if op == "M":
            m += length
        elif op == "=":
            eq += length
        elif op == "X":
            x += length
        elif op == "I":
            ins += length
        elif op == "D":
            dele += length
    return m, x, ins, dele, eq


def calculate_percent_identity(
    cigar: str,
    md_tag: str | None = None,
    nm_tag: int | None = None,
) -> float:
    """Calculate percent identity of an alignment.

    Preference order:

    1. ``NM`` tag: identity = 1 - NM / aligned_length, where aligned_length is
       the number of alignment columns (M/=/X + I + D) and NM is the edit
       distance (mismatches + inserted + deleted bases).
    2. ``MD`` tag: mismatches counted from MD, plus I and D from the CIGAR.
    3. CIGAR only: ``=``/``X`` operations are exact; plain ``M`` cannot
       distinguish matches from mismatches and is counted as matching, so the
       result is an upper bound.
    """
    m, x, ins, dele, eq = _cigar_counts(cigar)
    aligned = m + x + eq + ins + dele
    if aligned == 0:
        return 0.0

    if nm_tag is not None:
        return max(0.0, 100.0 * (1.0 - nm_tag / aligned))

    if md_tag:
        md_clean = re.sub(r"\^[A-Za-z]+", "", md_tag)
        mismatches = len(re.findall(r"[A-Za-z]", md_clean))
        return max(0.0, 100.0 * (1.0 - (mismatches + ins + dele) / aligned))

    matches = m + eq
    return 100.0 * matches / aligned


def parse_sam_line(line: str) -> AlignmentRecord | None:
    """Parse a single SAM line into an AlignmentRecord."""
    if line.startswith("@"):
        return None

    fields = line.strip().split("\t")
    if len(fields) < 11:
        return None

    # Skip unmapped (0x4), secondary (0x100) and supplementary (0x800)
    # alignments so each read is counted once.
    flag = int(fields[1])
    if flag & 0x904:
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

    # Optional NM / MD tags for an accurate identity
    nm_tag: int | None = None
    md_tag: str | None = None
    for tag in fields[11:]:
        if tag.startswith("NM:i:"):
            try:
                nm_tag = int(tag[5:])
            except ValueError:
                nm_tag = None
        elif tag.startswith("MD:Z:"):
            md_tag = tag[5:]

    percent_identity = calculate_percent_identity(cigar, md_tag=md_tag, nm_tag=nm_tag)

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

    with tempfile.TemporaryFile(mode="w+") as stderr_file:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=stderr_file,
            text=True,
        )
        completed = False
        returncode: int | None = None
        try:
            if process.stdout is None:
                return

            for line in process.stdout:
                record = parse_sam_line(line)
                if record is None:
                    continue

                if record.percent_identity >= min_identity:
                    yield record
            completed = True
        finally:
            # Also runs when the consumer stops iterating early (generator
            # close): stop samtools, release the pipe and reap the process.
            if not completed and process.poll() is None:
                process.terminate()
            close = getattr(process.stdout, "close", None)
            if close is not None:
                close()
            try:
                returncode = process.wait(timeout=10)
            except subprocess.TimeoutExpired:  # pragma: no cover
                process.kill()
                returncode = process.wait()

        if completed and isinstance(returncode, int) and returncode != 0:
            stderr_file.seek(0)
            msg = f"samtools view failed (exit {returncode}): {stderr_file.read().strip()}"
            raise RuntimeError(msg)


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

    df = pl.DataFrame(records)
    if df.height:
        try:
            lengths = get_contig_lengths(bam_path, samtools_path=samtools_path)
        except (OSError, ValueError):
            lengths = {}
        df = add_genome_positions(df, lengths)
    return df


def add_genome_positions(df: pl.DataFrame, contig_lengths: dict[str, int]) -> pl.DataFrame:
    """Add ``genome_position``: the position along the concatenated genome.

    ``position`` is local to a contig, so the contigs of a draft genome would
    all be drawn over the same range. Each contig is offset by the summed
    lengths of the genome's earlier contigs (BAM header order; contigs
    missing from the header use their largest observed position).
    """
    if df.is_empty() or not {"genome_name", "contig", "position"} <= set(df.columns):
        return df
    observed = dict(
        df.group_by("contig").agg(pl.col("position").max()).iter_rows()
    )
    genome_of = dict(df.select("contig", "genome_name").unique().iter_rows())
    order = list(contig_lengths) + sorted(c for c in observed if c not in contig_lengths)
    offsets: dict[str, int] = {}
    running: dict[str, int] = {}
    for contig in order:
        genome = genome_of.get(contig)
        if genome is None:
            continue  # contig without reads: no offset needed
        offsets[contig] = running.get(genome, 0)
        running[genome] = offsets[contig] + int(
            contig_lengths.get(contig) or observed.get(contig) or 0
        )
    return df.with_columns(
        (pl.col("position") + pl.col("contig").replace_strict(offsets, default=0))
        .alias("genome_position")
    )


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
