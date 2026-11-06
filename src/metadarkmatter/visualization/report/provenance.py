"""Loading and presenting the ``score classify`` provenance sidecar.

``mdm score classify`` writes ``<classifications>.provenance.json`` next to
its output (schema ``metadarkmatter.provenance/1``). The report shows a
"Run information" panel built from it; when the sidecar is missing or
unreadable the panel falls back to what can be derived from the data.

The sidecar's ``created_utc`` is displayed only inside that panel. It comes
from the sidecar (not the clock), so a report built twice from the same
inputs stays byte-identical apart from the generation timestamp in the
page header.
"""

from __future__ import annotations

import html
import json
import logging
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

PROVENANCE_SCHEMA = "metadarkmatter.provenance/1"
PROVENANCE_SUFFIX = ".provenance.json"

# Input roles in display order.
INPUT_ROLES: tuple[tuple[str, str], ...] = (
    ("alignment", "Alignment"),
    ("ani", "ANI matrix"),
    ("metadata", "Genome metadata"),
    ("aai", "AAI matrix"),
    ("id_mapping", "ID mapping"),
    ("genomes", "Genome directory (ID mapping)"),
)


def default_provenance_path(classifications: Path | str) -> Path:
    """Sidecar location for a classifications file: ``<file>.provenance.json``."""
    path = Path(classifications)
    return path.with_name(path.name + PROVENANCE_SUFFIX)


def load_provenance(path: Path | str) -> dict[str, Any] | None:
    """Read and validate a provenance sidecar.

    Returns ``None`` (after logging a warning) when the file cannot be
    read, is not a JSON object, or declares an unsupported schema.
    """
    p = Path(path)
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        logger.warning("Could not read provenance file %s: %s", p, e)
        return None
    if not isinstance(data, dict):
        logger.warning("Provenance file %s is not a JSON object", p)
        return None
    schema = data.get("schema")
    if not isinstance(schema, str) or not schema.startswith("metadarkmatter.provenance/"):
        logger.warning("Provenance file %s has unsupported schema %r", p, schema)
        return None
    if schema != PROVENANCE_SCHEMA:
        logger.warning(
            "Provenance file %s uses schema %s (expected %s); showing known fields only",
            p, schema, PROVENANCE_SCHEMA,
        )
    return data


def human_size(num_bytes: object) -> str:
    """Format a byte count (``1536`` -> ``1.5 KB``); non-numbers give ``-``."""
    if isinstance(num_bytes, bool) or not isinstance(num_bytes, (int, float)):
        return "-"
    size = float(num_bytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024 or unit == "TB":
            return f"{int(size)} B" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024
    return "-"  # pragma: no cover


def _e(value: object) -> str:
    return html.escape("" if value is None else str(value))


def render_rows(rows: list[tuple[str, str]]) -> str:
    """Render ``(label, already-escaped HTML value)`` pairs as a ``<dl>``."""
    if not rows:
        return ""
    body = "".join(f"<dt>{_e(k)}</dt><dd>{v}</dd>" for k, v in rows)
    return f"<dl>{body}</dl>"


def file_row_value(info: object, fallback_name: str | None = None) -> str:
    """HTML for an input-file entry (name, size, short sha256)."""
    if not isinstance(info, dict):
        return _e(fallback_name) if fallback_name else "not used"
    path = str(info.get("path", ""))
    name = Path(path).name if path else (fallback_name or "")
    parts = [f'<code title="{_e(path)}">{_e(name)}</code>']
    if "files" in info:
        parts.append(f"{info.get('files')} files")
    if "size" in info:
        parts.append(human_size(info.get("size")))
    sha = info.get("sha256")
    if isinstance(sha, str) and sha:
        parts.append(f'sha256 <code title="{_e(sha)}">{_e(sha[:12])}</code>')
    return " &middot; ".join(parts)


def format_setting(value: object) -> str:
    """Compact display of a ScoringConfig value."""
    if value is None:
        return "none"
    if isinstance(value, bool):
        return "on" if value else "off"
    if isinstance(value, float):
        return f"{value:g}"
    return str(value)
