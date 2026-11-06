"""Streaming extraction of selected reads from FASTQ/FASTA files.

Used by ``score extract-novel --reads-1 ... --fastq-out ...`` to turn a set of
classified read IDs back into sequence files for assembly. Inputs are streamed
record by record (gzip handled transparently), so memory use is bounded by the
ID set, never by the size of the FASTQ.

Read-ID matching normalises both sides: only the first whitespace-delimited
token of a header is used and a trailing ``/1`` or ``/2`` mate suffix is
removed. This reconciles IDs classified from a paired MMseqs2 search (which
appends ``/1`` and ``/2`` to the mates) with FASTQ headers that may or may not
carry the suffix. A pair is selected when either mate's ID is selected, so the
two mate files stay in sync for paired assembly.
"""

from __future__ import annotations

import gzip
import re
from collections.abc import Iterable, Iterator, Mapping
from contextlib import ExitStack
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO, NamedTuple

from metadarkmatter.core.exceptions import ConfigurationError

_MATE_SUFFIX = ("/1", "/2")
_SEQ_SUFFIXES = (
    ".fastq", ".fq", ".fasta", ".fa", ".fna",
)
# Maximum number of per-lineage files kept open at once (the rest are closed
# and re-opened in append mode on demand; gzip members concatenate validly).
_MAX_OPEN_LINEAGE_FILES = 64


def normalize_read_id(raw: str) -> str:
    """Return the canonical read ID: first token, no leading marker, no /1 or /2."""
    token = raw.strip().split(None, 1)[0] if raw.strip() else ""
    if token[:1] in ("@", ">"):
        token = token[1:]
    if token.endswith(_MATE_SUFFIX):
        token = token[:-2]
    return token


class SeqRecord(NamedTuple):
    """One sequence record; ``qual`` is ``None`` for FASTA."""

    header: bytes  # without the leading '@' / '>' and without newline
    seq: bytes
    qual: bytes | None

    @property
    def read_id(self) -> str:
        return normalize_read_id(self.header.decode("utf-8", errors="replace"))


def _open_binary(path: Path) -> IO[bytes]:
    with path.open("rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rb")  # type: ignore[return-value]
    return path.open("rb")


def detect_format(path: Path) -> str:
    """Return ``"fastq"`` or ``"fasta"`` from the first byte of the file."""
    with _open_binary(path) as fh:
        first = fh.read(1)
    if first == b"@":
        return "fastq"
    if first == b">":
        return "fasta"
    if first == b"":
        return "fastq"  # empty file: treat as FASTQ, yields no records
    raise ConfigurationError(
        f"{path} is not FASTQ or FASTA (first character {first!r}).",
        suggestion="Provide an uncompressed or gzip-compressed FASTQ/FASTA file.",
    )


def iter_records(path: Path) -> Iterator[SeqRecord]:
    """Stream records from a FASTQ or FASTA file (optionally gzipped)."""
    fmt = detect_format(path)
    with _open_binary(path) as fh:
        if fmt == "fastq":
            while True:
                head = fh.readline()
                if not head:
                    return
                if not head.strip():
                    continue
                seq = fh.readline().rstrip(b"\r\n")
                fh.readline()  # '+' separator
                qual = fh.readline().rstrip(b"\r\n")
                if not head.startswith(b"@") or (not qual and seq):
                    raise ConfigurationError(f"Malformed FASTQ record near {head[:60]!r} in {path}.")
                yield SeqRecord(head[1:].rstrip(b"\r\n"), seq, qual)
        else:
            header: bytes | None = None
            parts: list[bytes] = []
            for line in fh:
                if line.startswith(b">"):
                    if header is not None:
                        yield SeqRecord(header, b"".join(parts), None)
                    header = line[1:].rstrip(b"\r\n")
                    parts = []
                elif header is not None:
                    parts.append(line.strip())
            if header is not None:
                yield SeqRecord(header, b"".join(parts), None)


def format_record(rec: SeqRecord, fmt: str) -> bytes:
    """Serialise a record as ``fmt`` ('fastq' or 'fasta')."""
    if fmt == "fastq":
        if rec.qual is None:
            raise ConfigurationError("Cannot write FASTQ from FASTA input (no qualities).")
        return b"@" + rec.header + b"\n" + rec.seq + b"\n+\n" + rec.qual + b"\n"
    return b">" + rec.header + b"\n" + rec.seq + b"\n"


def mate_output_path(path: Path, mate: int | None) -> Path:
    """``novel.fastq.gz`` -> ``novel_R1.fastq.gz`` for ``mate=1`` (unchanged for None)."""
    if mate is None:
        return path
    name = path.name
    gz = name.endswith(".gz")
    core = name[:-3] if gz else name
    suffix = ""
    for s in _SEQ_SUFFIXES:
        if core.endswith(s):
            suffix = s
            core = core[: -len(s)]
            break
    return path.with_name(f"{core}_R{mate}{suffix}{'.gz' if gz else ''}")


def _open_out(path: Path, append: bool = False) -> IO[bytes]:
    path.parent.mkdir(parents=True, exist_ok=True)
    mode = "ab" if append else "wb"
    if path.name.endswith(".gz"):
        return gzip.open(path, mode)  # type: ignore[return-value]
    return path.open(mode)


class _LineageWriters:
    """Per-lineage output files with a bounded number of open handles."""

    def __init__(self, directory: Path, lineages: Iterable[str], fmt: str, paired: bool):
        self._fmt = fmt
        self._paired = paired
        self._dir = directory
        self._handles: dict[Path, IO[bytes]] = {}
        self._seen: set[Path] = set()
        self.paths: dict[str, Path] = {}
        taken: set[str] = set()
        for lineage in sorted(set(lineages)):
            safe = re.sub(r"[^A-Za-z0-9._-]+", "_", lineage).strip("._") or "lineage"
            name, n = safe, 0
            # A dedupe suffix may itself equal another lineage's name, so
            # loop until the name is unused.
            while name in taken:
                n += 1
                name = f"{safe}_{n}"
            taken.add(name)
            self.paths[lineage] = directory / name

    def path_for(self, lineage: str, mate: int | None) -> Path:
        base = self.paths[lineage]
        ext = ".fastq.gz" if self._fmt == "fastq" else ".fasta.gz"
        tag = f"_R{mate}" if (mate is not None and self._paired) else ""
        return base.with_name(base.name + tag + ext)

    def write(self, lineage: str, mate: int | None, data: bytes) -> None:
        path = self.path_for(lineage, mate)
        handle = self._handles.get(path)
        if handle is None:
            if len(self._handles) >= _MAX_OPEN_LINEAGE_FILES:
                oldest = next(iter(self._handles))
                self._handles.pop(oldest).close()
            handle = _open_out(path, append=path in self._seen)
            self._seen.add(path)
            self._handles[path] = handle
        handle.write(data)

    def touch_all(self) -> None:
        """Create (empty) files for every lineage, mates included."""
        for lineage in self.paths:
            for mate in ((1, 2) if self._paired else (None,)):
                path = self.path_for(lineage, mate)
                if path not in self._seen:
                    _open_out(path).close()
                    self._seen.add(path)

    def close(self) -> None:
        for h in self._handles.values():
            h.close()
        self._handles.clear()


@dataclass
class ExtractionReport:
    """Counts from :func:`extract_reads`."""

    requested: int = 0
    found_1: int = 0
    found_2: int | None = None
    written_1: int = 0
    written_2: int | None = None
    missing: list[str] = field(default_factory=list)  # IDs absent from reads_1
    missing_2: list[str] = field(default_factory=list)  # IDs absent from reads_2
    outputs: list[Path] = field(default_factory=list)
    lineage_files: list[Path] = field(default_factory=list)


def _check_output_paths(
    reads_1: Path,
    reads_2: Path | None,
    fastq_out: Path | None,
    fasta_out: Path | None,
) -> None:
    """Refuse outputs that would overwrite an input or each other.

    Output files are opened for writing (truncated) before the inputs are
    read, so an output equal to an input would destroy the reads.
    """
    paired = reads_2 is not None
    inputs = {p.resolve() for p in (reads_1, reads_2) if p is not None}
    targets: list[Path] = []
    for out_path in (fastq_out, fasta_out):
        if out_path is not None:
            for mate in ((1, 2) if paired else (None,)):
                targets.append(mate_output_path(out_path, mate).resolve())
    for target in targets:
        if target in inputs:
            raise ConfigurationError(
                f"Output {target} is also an input read file; it would be overwritten.",
                suggestion="Choose a different --fastq-out/--fasta-out path.",
            )
    if len(set(targets)) != len(targets):
        raise ConfigurationError("--fastq-out and --fasta-out resolve to the same file.")


def extract_reads(
    *,
    reads_1: Path,
    ids: set[str],
    reads_2: Path | None = None,
    fastq_out: Path | None = None,
    fasta_out: Path | None = None,
    per_lineage_dir: Path | None = None,
    lineage_of: Mapping[str, str] | None = None,
) -> ExtractionReport:
    """Stream ``reads_1`` (and ``reads_2``) and write records whose ID is in ``ids``.

    Args:
        reads_1: FASTQ/FASTA(.gz) with the (first) mates / single-end reads.
        ids: Normalised read IDs to keep (see :func:`normalize_read_id`).
        reads_2: Optional second-mate file. When given, outputs are written to
            ``{stem}_R1``/``{stem}_R2`` files.
        fastq_out: FASTQ destination (gzip when the name ends in ``.gz``).
            Requires FASTQ input.
        fasta_out: FASTA destination; FASTQ input is converted.
        per_lineage_dir: Directory receiving one gzipped file per lineage in
            ``lineage_of`` (input format; ``_R1``/``_R2`` when paired).
        lineage_of: Mapping normalised read ID -> lineage name.

    Raises:
        ConfigurationError: Unusable output/input combination or malformed input.
    """
    if fastq_out is None and fasta_out is None and per_lineage_dir is None:
        raise ConfigurationError("No output requested for read extraction.")
    if per_lineage_dir is not None and not lineage_of and ids:
        raise ConfigurationError("per-lineage output needs a read-to-lineage mapping.")

    in_fmt = detect_format(reads_1)
    if reads_2 is not None and detect_format(reads_2) != in_fmt:
        raise ConfigurationError("--reads-1 and --reads-2 must both be FASTQ or both FASTA.")
    if fastq_out is not None and in_fmt != "fastq":
        raise ConfigurationError(
            "--fastq-out needs FASTQ input; the reads provided are FASTA.",
            suggestion="Use --fasta-out instead.",
        )

    paired = reads_2 is not None
    _check_output_paths(reads_1, reads_2, fastq_out, fasta_out)
    report = ExtractionReport(requested=len(ids))
    mates = [(1, reads_1)] + ([(2, reads_2)] if reads_2 is not None else [])
    found_sets: dict[int, set[str]] = {1: set(), 2: set()}
    written = {1: 0, 2: 0}

    lineage_writers: _LineageWriters | None = None
    if per_lineage_dir is not None:
        assert lineage_of is not None or not ids
        lineage_writers = _LineageWriters(
            per_lineage_dir, (lineage_of or {}).values(), in_fmt, paired
        )
        inputs = {p.resolve() for p in (reads_1, reads_2) if p is not None}
        for lin in lineage_writers.paths:
            for m in ((1, 2) if paired else (None,)):
                if lineage_writers.path_for(lin, m).resolve() in inputs:
                    raise ConfigurationError(
                        f"Per-lineage file {lineage_writers.path_for(lin, m)} would "
                        "overwrite an input read file.",
                        suggestion="Use a separate --per-lineage-dir.",
                    )

    with ExitStack() as stack:
        if lineage_writers is not None:
            stack.callback(lineage_writers.close)
        for idx, (mate, path) in enumerate(mates, start=1):
            sinks: list[tuple[IO[bytes], str]] = []
            for out_path, fmt in ((fastq_out, "fastq"), (fasta_out, "fasta")):
                if out_path is None:
                    continue
                target = mate_output_path(out_path, mate if paired else None)
                sinks.append((stack.enter_context(_open_out(target)), fmt))
                report.outputs.append(target)
            for rec in iter_records(path):
                rid = rec.read_id
                if rid not in ids:
                    continue
                found_sets[idx].add(rid)
                written[idx] += 1
                for handle, fmt in sinks:
                    handle.write(format_record(rec, fmt))
                if lineage_writers is not None and lineage_of is not None:
                    lineage = lineage_of.get(rid)
                    if lineage is not None:
                        lineage_writers.write(
                            lineage, mate if paired else None, format_record(rec, in_fmt)
                        )
        if lineage_writers is not None:
            lineage_writers.touch_all()
            report.lineage_files = sorted(
                {
                    lineage_writers.path_for(lin, m)
                    for lin in lineage_writers.paths
                    for m in ((1, 2) if paired else (None,))
                }
            )

    report.found_1 = len(found_sets[1])
    report.written_1 = written[1]
    report.missing = sorted(ids - found_sets[1])
    if paired:
        report.found_2 = len(found_sets[2])
        report.written_2 = written[2]
        report.missing_2 = sorted(ids - found_sets[2])
    return report
