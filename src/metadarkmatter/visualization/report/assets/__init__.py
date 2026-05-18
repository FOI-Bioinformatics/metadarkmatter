"""
Bundled and on-demand JavaScript assets for the HTML report.

Provides script-tag builders for Plotly and D3. Plotly is shipped with
the Python `plotly` package and is always embeddable inline. D3 is
not bundled by default; if a vendored copy is present in this
directory (``d3.v7.min.js``) it is inlined, otherwise the CDN script tag
is emitted. Use the ``vendor`` subcommand to download and cache D3
locally for offline use.
"""

from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger(__name__)

_ASSETS_DIR = Path(__file__).parent
D3_FILENAME = "d3.v7.min.js"
D3_CDN_URL = "https://d3js.org/d3.v7.min.js"


def get_plotly_script_tag(mode: str) -> str:
    """Return a ``<script>`` block providing Plotly.

    Args:
        mode: ``"offline"`` to inline the JS shipped with the plotly
            Python package, ``"cdn"`` to emit a remote ``<script src>``.

    Returns:
        HTML ``<script>`` block (single tag) suitable for inclusion in the
        report ``<head>``.
    """
    if mode == "cdn":
        return (
            '    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>'
        )

    try:
        from plotly.offline import get_plotlyjs
    except ImportError:  # pragma: no cover - plotly is a hard dep
        logger.warning(
            "plotly.offline.get_plotlyjs unavailable; falling back to CDN.",
        )
        return (
            '    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>'
        )

    return f"<script>\n{get_plotlyjs()}\n</script>"


def get_d3_script_tag(mode: str) -> str:
    """Return a ``<script>`` block providing D3 v7.

    Args:
        mode: ``"offline"`` to inline a vendored copy of D3 if present
            (with CDN fallback when not vendored), ``"cdn"`` to always emit
            the CDN tag.

    Returns:
        HTML ``<script>`` block.
    """
    if mode == "cdn":
        return f'<script src="{D3_CDN_URL}"></script>'

    vendored = _ASSETS_DIR / D3_FILENAME
    if vendored.is_file():
        return f"<script>\n{vendored.read_text(encoding='utf-8')}\n</script>"

    logger.info(
        "D3 not vendored at %s; falling back to CDN. "
        "Run 'metadarkmatter report vendor-d3' to enable offline phylogeny.",
        vendored,
    )
    return f'<script src="{D3_CDN_URL}"></script>'
