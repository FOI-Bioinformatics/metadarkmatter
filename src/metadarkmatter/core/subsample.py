"""Reproducible, mate-preserving subsampling of FASTQ read files.

Records are parsed strictly as four-line units, so quality lines that begin
with "@" are not mistaken for headers. Files are streamed; only the set of
selected record indices is held in memory.
"""

from __future__ import annotations

import gzip
import io
import os
import random
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import IO

from metadarkmatter.core.paired_reads import open_reads

# Fast gzip level: favours speed over compression ratio for temporary reads.
_GZIP_LEVEL = 1


@dataclass(frozen=True)
class SubsampleResult:
    """Outcome of a subsampling call."""

    total: int
    kept: int
    seed: int


def _count_records(path: Path) -> int:
    """Count FASTQ records, raising ValueError if the file is truncated."""
    n_lines = 0
    with open_reads(path) as fh:
        for _ in fh:
            n_lines += 1
    if n_lines % 4 != 0:
        raise ValueError(
            f"Truncated FASTQ record in {path}: {n_lines} lines is not a multiple of 4"
        )
    return n_lines // 4


def _open_output(path: Path, gzipped: bool) -> IO[str]:
    if gzipped:
        # Fixed header (no file name, mtime 0) keeps output byte-reproducible.
        handle = path.open("wb")
        raw = gzip.GzipFile(
            filename="",
            mode="wb",
            compresslevel=_GZIP_LEVEL,
            mtime=0,
            fileobj=handle,
        )
        raw.myfileobj = handle  # type: ignore[assignment]  # closed together with the gzip stream
        return io.TextIOWrapper(raw)
    return path.open("w")


def _write_subset(src: Path, dst: Path, indices: Sequence[int], gzipped: bool) -> None:
    """Write the records of ``src`` at the sorted ``indices`` to ``dst``."""
    wanted = iter(indices)
    next_idx = next(wanted, None)
    record = 0
    buf: list[str] = []
    with open_reads(src) as fin, _open_output(dst, gzipped) as fout:
        for line in fin:
            if next_idx is None:
                break
            buf.append(line)
            if len(buf) == 4:
                if record == next_idx:
                    fout.writelines(buf)
                    next_idx = next(wanted, None)
                buf.clear()
                record += 1
    if next_idx is not None:
        raise ValueError(f"Fewer records than expected while reading {src}")


def subsample_reads(
    reads_1: Path, reads_2: Path | None, max_pairs: int, seed: int
) -> SubsampleResult:
    """Subsample FASTQ files in place to at most ``max_pairs`` records.

    The same record indices are kept in both mates, input order is preserved,
    and the selection depends only on the record count, ``max_pairs`` and
    ``seed``. Files are left untouched when they hold ``max_pairs`` records or
    fewer. Output is gzip-compressed when the path ends in ".gz".

    Raises:
        ValueError: If the mates differ in record count or a file ends in a
            truncated record.
    """
    total = _count_records(reads_1)
    if reads_2 is not None:
        total_2 = _count_records(reads_2)
        if total_2 != total:
            raise ValueError(
                f"Read count mismatch between mates: {reads_1} has {total} "
                f"records, {reads_2} has {total_2}"
            )

    if total <= max_pairs:
        return SubsampleResult(total=total, kept=total, seed=seed)

    indices = sorted(random.Random(seed).sample(range(total), max_pairs))
    targets = [reads_1] if reads_2 is None else [reads_1, reads_2]

    temps: list[Path] = []
    try:
        for target in targets:
            fd, name = tempfile.mkstemp(prefix=f".{target.name}.", suffix=".tmp", dir=target.parent)
            os.close(fd)
            tmp = Path(name)
            temps.append(tmp)
            _write_subset(target, tmp, indices, target.name.endswith(".gz"))
        for tmp, target in zip(temps, targets, strict=True):
            tmp.replace(target)
    finally:
        for tmp in temps:
            tmp.unlink(missing_ok=True)

    return SubsampleResult(total=total, kept=max_pairs, seed=seed)
