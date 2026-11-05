"""
Paired-end read helpers shared by the BLAST and MMseqs2 aligners.

Neither aligner has native paired-end support, so mates are combined into one
query file and searched independently. Mates usually share an ID (Illumina/SRA
style); after combining, the classifier would merge both mates into a single
read. Each mate therefore gets a distinct ID by appending ``/1`` or ``/2``
unless it already carries that suffix. Compression is detected per file, and
every record is written with a trailing newline so files without one cannot
merge records.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import IO


def open_reads(path: Path) -> IO[str]:
    """Open a FASTA/FASTQ file as text, detecting gzip from its magic bytes."""
    with path.open("rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return path.open("r")


def _mate_id(header_line: str, mate: int) -> str:
    """Return a header line whose ID ends in /``mate`` (unchanged if already suffixed)."""
    marker = header_line[0]
    body = header_line[1:].rstrip("\r\n")
    parts = body.split(maxsplit=1)
    if not parts:
        return f"{marker}{body}\n"
    read_id = parts[0]
    # Only "/1" / "/2" count as existing mate suffixes: SRA spot IDs such as
    # "SRR123.1" are identical across mates and must still be disambiguated.
    if not read_id.endswith(f"/{mate}"):
        read_id = f"{read_id}/{mate}"
    rest = f" {parts[1]}" if len(parts) > 1 else ""
    return f"{marker}{read_id}{rest}\n"


def _write_mate_records(reads_file: Path, mate: int, out_f: IO[str]) -> int:
    """Copy one mate file to ``out_f`` with distinct mate IDs; return record count."""
    count = 0
    with open_reads(reads_file) as in_f:
        first = ""
        for candidate in in_f:
            if candidate.strip():
                first = candidate
                break
        if not first:
            return 0
        if first.startswith("@"):
            # FASTQ: 4-line records
            line: str | None = first
            while line:
                if not line.strip():
                    line = in_f.readline()
                    continue
                seq = in_f.readline()
                plus = in_f.readline()
                qual = in_f.readline()
                out_f.write(_mate_id(line, mate))
                for part in (seq, plus, qual):
                    out_f.write(part if part.endswith("\n") else part + "\n")
                count += 1
                line = in_f.readline()
        else:
            # FASTA: header + one or more sequence lines
            for raw in [first, *in_f]:
                if raw.startswith(">"):
                    out_f.write(_mate_id(raw, mate))
                    count += 1
                elif raw.strip():
                    out_f.write(raw if raw.endswith("\n") else raw + "\n")
    return count


def combine_mates(reads_1: Path, reads_2: Path, output_path: Path) -> int:
    """Write R1 then R2 to ``output_path`` with distinct mate IDs.

    Records of one mate are contiguous and R1 precedes R2. Returns the total
    number of reads written.
    """
    read_count = 0
    with output_path.open("w") as out_f:
        for mate, reads_file in ((1, reads_1), (2, reads_2)):
            read_count += _write_mate_records(reads_file, mate, out_f)
    return read_count
