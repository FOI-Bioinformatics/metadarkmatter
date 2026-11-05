"""Helpers for embedding data (JSON, CSV exports) inside the HTML report."""

from __future__ import annotations

import html
import json
from collections.abc import Sequence
from typing import Any


def json_for_script(data: object, *, compact: bool = False) -> str:
    """Serialize ``data`` as JSON safe to embed inside a ``<script>`` element.

    ``<``, ``>`` and ``&`` are emitted as unicode escapes so values such as
    ``</script>`` in metadata species names cannot terminate the script
    block; the result is still a valid JS/JSON literal.

    Args:
        data: JSON-serialisable object.
        compact: Use the shortest separators (no spaces) to minimise size.
    """
    text = (
        json.dumps(data, separators=(",", ":"), allow_nan=False)
        if compact
        else json.dumps(data)
    )
    return (
        text.replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("&", "\\u0026")
        .replace(chr(0x2028), "\\u2028")
        .replace(chr(0x2029), "\\u2029")
    )


def json_script_tag(element_id: str, data: object) -> str:
    """Return a ``<script type="application/json">`` element holding ``data``."""
    return (
        f'<script type="application/json" id="{html.escape(element_id, quote=True)}">'
        f"{json_for_script(data, compact=True)}</script>"
    )


def csv_export_payload(
    filename: str,
    header: Sequence[str],
    rows: Sequence[Sequence[Any]],
) -> dict[str, Any]:
    """Build the JSON payload consumed by the generic ``mdmDownloadCsv`` JS."""
    return {
        "filename": filename,
        "header": list(header),
        "rows": [list(r) for r in rows],
    }


def csv_export_button(export_id: str, payload: dict[str, Any], label: str = "Download CSV") -> str:
    """Return a download button plus its embedded data for a client-side CSV.

    The data is stored in a JSON ``<script>`` element so the button works
    offline and the CSV is generated entirely in the browser.
    """
    safe_id = html.escape(export_id, quote=True)
    return (
        f'<button type="button" class="csv-btn" data-export="{safe_id}" '
        f"onclick=\"mdmDownloadCsv('{safe_id}')\">{html.escape(label)}</button>"
        + json_script_tag(f"mdm-export-{export_id}", payload)
    )


CSV_EXPORT_JS: str = r"""
// Generic client-side CSV download for tables embedded as JSON
function mdmCsvCell(v) {
    if (v === null || v === undefined) return '';
    var s = String(v);
    if (/[",\r\n]/.test(s)) s = '"' + s.replace(/"/g, '""') + '"';
    return s;
}

function mdmDownloadBlob(text, filename, mime) {
    var blob = new Blob([text], { type: mime });
    var url = URL.createObjectURL(blob);
    var a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    setTimeout(function() { URL.revokeObjectURL(url); }, 1000);
}

function mdmDownloadCsv(exportId) {
    var el = document.getElementById('mdm-export-' + exportId);
    if (!el) return;
    var d = JSON.parse(el.textContent);
    var lines = [d.header.map(mdmCsvCell).join(',')];
    d.rows.forEach(function(r) { lines.push(r.map(mdmCsvCell).join(',')); });
    mdmDownloadBlob(lines.join('\n') + '\n', d.filename, 'text/csv;charset=utf-8');
}
"""
