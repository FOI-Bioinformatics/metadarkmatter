"""
Lightweight FASTA and FASTQ format validators.

These validators are intentionally streaming and dependency-free so they
can run as a first-stop check on large input files without loading them
into memory or invoking BioPython. The aim is to catch the common forms
of malformed input before they cause cryptic failures in BLAST,
MMseqs2, or Kraken2: truncated files, missing header lines, mixed
formats, and bad quality strings.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from io import TextIOBase
from pathlib import Path
from typing import IO

# Permitted characters per record type. We allow ambiguity codes (IUPAC)
# and gap characters since assemblies and aligned references can contain
# them; sequencing reads will not in practice.
_NUCLEOTIDE_CHARS = set("ACGTUNRYSWKMBDHVacgtunryswkmbdhv-.*")
_PROTEIN_CHARS = set("ACDEFGHIKLMNPQRSTVWYBZXJUOacdefghiklmnpqrstvwybzxjuo-.*")


@dataclass
class ValidationResult:
    """Outcome of a single-file validation."""

    path: Path
    record_count: int = 0
    issues: list[str] = field(default_factory=list)
    sample_record_id: str | None = None

    @property
    def ok(self) -> bool:
        return not self.issues


def _open_text(path: Path) -> IO[str]:
    """Open a possibly-gzipped file in text mode."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def validate_fasta(
    path: Path,
    *,
    max_records: int | None = None,
    sequence_type: str = "nucleotide",
) -> ValidationResult:
    """Validate a FASTA file.

    Streams the file line by line and reports the first issues it finds
    (capped at 20 distinct issues). A record is defined as a header line
    starting with ``>`` followed by one or more sequence lines.

    Args:
        path: FASTA file (gzipped supported via the ``.gz`` suffix).
        max_records: Stop after this many records. ``None`` reads the
            whole file.
        sequence_type: ``"nucleotide"`` (default) or ``"protein"``,
            chooses the allowed-character set.

    Returns:
        A :class:`ValidationResult`. ``ok`` is true when no issues were
        found.
    """
    allowed = _NUCLEOTIDE_CHARS if sequence_type == "nucleotide" else _PROTEIN_CHARS
    result = ValidationResult(path=path)
    in_record = False
    current_header: str | None = None
    current_has_sequence = False

    if not path.exists():
        result.issues.append(f"File does not exist: {path}")
        return result

    try:
        handle = _open_text(path)
    except OSError as exc:
        result.issues.append(f"Cannot open file: {exc}")
        return result

    with handle:
        for lineno, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None and not current_has_sequence:
                    result.issues.append(
                        f"Line {lineno - 1}: header '{current_header[:40]}' "
                        "has no sequence body"
                    )
                current_header = line[1:].strip()
                if not current_header:
                    result.issues.append(f"Line {lineno}: empty FASTA header")
                if result.sample_record_id is None:
                    result.sample_record_id = current_header.split()[0] if current_header else None
                current_has_sequence = False
                in_record = True
                result.record_count += 1
                if max_records is not None and result.record_count > max_records:
                    break
            else:
                if not in_record:
                    result.issues.append(
                        f"Line {lineno}: sequence data before any '>' header"
                    )
                    in_record = True  # avoid repeating
                current_has_sequence = True
                bad_chars = set(line) - allowed
                if bad_chars and len(result.issues) < 20:
                    result.issues.append(
                        f"Line {lineno}: invalid {sequence_type} characters: "
                        f"{sorted(bad_chars)[:5]}"
                    )

    if current_header is not None and not current_has_sequence:
        result.issues.append(
            f"Header '{current_header[:40]}' at end of file has no sequence body"
        )
    if result.record_count == 0 and not result.issues:
        result.issues.append("No FASTA records found")
    return result


def validate_fastq(
    path: Path,
    *,
    max_records: int | None = None,
) -> ValidationResult:
    """Validate a FASTQ file (4-line record structure).

    Streams the file and checks that records follow the
    ``@header / sequence / + / quality`` layout, that sequence and
    quality strings are equal length, and that quality scores fall in
    the printable-ASCII range (33-126, covering Phred+33 and Phred+64).

    Args:
        path: FASTQ file (gzipped supported via the ``.gz`` suffix).
        max_records: Stop after this many records.

    Returns:
        :class:`ValidationResult`.
    """
    result = ValidationResult(path=path)

    if not path.exists():
        result.issues.append(f"File does not exist: {path}")
        return result

    try:
        handle = _open_text(path)
    except OSError as exc:
        result.issues.append(f"Cannot open file: {exc}")
        return result

    with handle:
        record_lines: list[str] = []
        lineno = 0
        for raw_line in handle:
            lineno += 1
            line = raw_line.rstrip("\r\n")
            record_lines.append(line)
            if len(record_lines) == 4:
                header, seq, plus, qual = record_lines
                record_lines = []
                result.record_count += 1

                if not header.startswith("@"):
                    result.issues.append(
                        f"Line {lineno - 3}: header must start with '@', got "
                        f"'{header[:20]}'"
                    )
                if not plus.startswith("+"):
                    result.issues.append(
                        f"Line {lineno - 1}: third line must start with '+', got "
                        f"'{plus[:20]}'"
                    )
                if len(seq) != len(qual):
                    result.issues.append(
                        f"Record {result.record_count}: sequence length "
                        f"({len(seq)}) != quality length ({len(qual)})"
                    )
                bad_q = [c for c in qual if not 33 <= ord(c) <= 126]
                if bad_q and len(result.issues) < 20:
                    result.issues.append(
                        f"Record {result.record_count}: quality string contains "
                        f"non-printable characters"
                    )
                if result.sample_record_id is None and header.startswith("@"):
                    result.sample_record_id = header[1:].split()[0]

                if max_records is not None and result.record_count >= max_records:
                    break

    if record_lines:
        result.issues.append(
            f"Trailing partial record: {len(record_lines)} lines after the "
            "last complete record (FASTQ requires multiples of 4 lines)"
        )
    if result.record_count == 0 and not result.issues:
        result.issues.append("No FASTQ records found")
    return result
