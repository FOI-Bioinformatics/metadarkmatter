"""Provenance sidecars for classification runs.

``score classify`` (and ``score batch``, once per sample) writes
``<output>.provenance.json`` next to its classification table. The sidecar
records the tool version, the command line, content hashes of the input files,
the fully resolved scoring configuration and a compact result summary, so a
report can state exactly how a table was produced.

Schema ``metadarkmatter.provenance/1`` (consumed by the report generator)::

    {
      "schema": "metadarkmatter.provenance/1",
      "metadarkmatter_version": str,
      "created_utc": ISO-8601 str,
      "command": [argv strings],
      "python": str, "platform": str, "seed": int,
      "inputs": {"alignment", "ani", "metadata", "aai", "id_mapping"},
      "config": ScoringConfig.model_dump(mode="json") | null,
      "outputs": {"classifications": str, "summary": str | null},
      "results": {"total_reads": int, "taxonomic_call_counts": {...}}
    }

Each input is ``{"path": absolute, "size": int, "sha256": str}`` or ``null``.

Everything except ``created_utc`` is a pure function of the inputs and the
command line, so two identical runs produce identical sidecars apart from the
timestamp.
"""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from collections.abc import Mapping, Sequence
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import polars as pl

from metadarkmatter import __version__
from metadarkmatter.core.random import get_seed
from metadarkmatter.models.config import ScoringConfig

PROVENANCE_SCHEMA = "metadarkmatter.provenance/1"
PROVENANCE_SUFFIX = ".provenance.json"

# Input slots that are always present in the sidecar (null when unused).
INPUT_SLOTS: tuple[str, ...] = ("alignment", "ani", "metadata", "aai", "id_mapping", "genomes")

_HASH_CHUNK = 1024 * 1024


def sha256_file(path: Path) -> str:
    """Return the hex SHA-256 of a file, streaming it in 1 MiB chunks."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(_HASH_CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fingerprint_file(path: Path) -> dict[str, Any]:
    """Return ``{"path": absolute path, "size": bytes, "sha256": hex}``."""
    resolved = Path(path).resolve()
    return {
        "path": str(resolved),
        "size": resolved.stat().st_size,
        "sha256": sha256_file(resolved),
    }


def fingerprint_dir(path: Path) -> dict[str, Any]:
    """Fingerprint a directory without reading file contents.

    Returns ``{"path", "files", "size", "listing_sha256"}``, where
    ``listing_sha256`` hashes the sorted (relative path, size) listing of all
    regular files. Genome directories can hold gigabytes, so contents are not
    hashed; a renamed, added, removed or resized file changes the value.
    """
    resolved = Path(path).resolve()
    entries = sorted(
        (str(f.relative_to(resolved)), f.stat().st_size)
        for f in resolved.rglob("*") if f.is_file()
    )
    digest = hashlib.sha256()
    for name, size in entries:
        digest.update(f"{name}\t{size}\n".encode())
    return {
        "path": str(resolved),
        "files": len(entries),
        "size": sum(size for _, size in entries),
        "listing_sha256": digest.hexdigest(),
    }


def provenance_path_for(output: Path) -> Path:
    """Sidecar location for a classification output (``<output>.provenance.json``)."""
    return output.with_name(output.name + PROVENANCE_SUFFIX)


def summarize_classifications(table: pl.DataFrame | pl.LazyFrame) -> dict[str, Any]:
    """Compact result summary: read count and per-category counts (sorted by name).

    Accepts a LazyFrame so streaming runs can summarise the on-disk output
    without loading it into memory.
    """
    lazy = table.lazy()
    if "taxonomic_call" not in lazy.collect_schema().names():
        return {"total_reads": int(lazy.select(pl.len()).collect().item()),
                "taxonomic_call_counts": {}}
    grouped = lazy.group_by("taxonomic_call").len().collect()
    counts = {
        str(call): int(n)
        for call, n in zip(grouped["taxonomic_call"], grouped["len"], strict=True)
    }
    return {
        "total_reads": sum(counts.values()),
        "taxonomic_call_counts": dict(sorted(counts.items())),
    }


def scan_classifications(path: Path) -> pl.LazyFrame:
    """Lazily scan a classification table written as CSV or Parquet."""
    if path.suffix.lower() == ".parquet":
        return pl.scan_parquet(path)
    return pl.scan_csv(path)


def build_provenance(
    command_argv: Sequence[str],
    inputs: Mapping[str, Path | None],
    config: ScoringConfig | None,
    extra: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Assemble the provenance dictionary.

    Args:
        command_argv: The command line (typically ``sys.argv``).
        inputs: Input name to path. ``None`` entries and the slots in
            :data:`INPUT_SLOTS` that are not supplied are recorded as ``null``.
            Each existing file is hashed by streaming.
        config: The scoring configuration actually used, or ``None``.
        extra: Additional top-level keys, merged last (for example
            ``outputs`` and ``results``).
    """
    input_records: dict[str, Any] = dict.fromkeys(INPUT_SLOTS)
    for name, path in inputs.items():
        if path is None:
            input_records[name] = None
        elif Path(path).is_dir():
            input_records[name] = fingerprint_dir(path)
        else:
            input_records[name] = fingerprint_file(path)

    data: dict[str, Any] = {
        "schema": PROVENANCE_SCHEMA,
        "metadarkmatter_version": __version__,
        "created_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "command": [str(a) for a in command_argv],
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "seed": get_seed(),
        "inputs": input_records,
        "config": config.model_dump(mode="json") if config is not None else None,
    }
    if extra:
        data.update(extra)
    return data


def write_provenance(path: Path, data: Mapping[str, Any]) -> None:
    """Write the provenance dictionary as indented JSON (insertion-ordered keys)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
