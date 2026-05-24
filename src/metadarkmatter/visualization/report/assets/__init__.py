"""
Bundled and on-demand JavaScript assets for the HTML report.

Provides script-tag builders for Plotly and D3. Plotly is shipped with
the Python `plotly` package and is always embeddable inline. D3 is
not bundled by default; if a vendored copy is present in this
directory (``d3.v7.min.js``) it is inlined, otherwise the CDN script
tag is emitted. When CDN script tags are used, ``Subresource
Integrity`` attributes (``integrity`` + ``crossorigin="anonymous"``)
are attached when a precomputed SHA-384 has been cached via
``metadarkmatter report vendor-cdn-sri``. Use the ``vendor-d3``
subcommand to download and cache D3 locally for offline use.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

_ASSETS_DIR = Path(__file__).parent
D3_FILENAME = "d3.v7.min.js"
D3_CDN_URL = "https://d3js.org/d3.v7.min.js"
PLOTLY_CDN_URL = "https://cdn.plot.ly/plotly-2.27.0.min.js"
SRI_CACHE_FILENAME = "cdn_sri.json"


def _load_sri_cache() -> dict[str, str]:
    """Return the SHA-384 cache for CDN URLs, or an empty dict.

    Read on demand so a freshly vendored file is picked up without a
    process restart, and so a malformed JSON file is silently ignored
    rather than blocking report generation.
    """
    path = _ASSETS_DIR / SRI_CACHE_FILENAME
    if not path.is_file():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        logger.debug("Could not read SRI cache at %s: %s", path, exc)
        return {}


def _cdn_script_tag(url: str, sri_key: str) -> str:
    """Build a CDN ``<script src>`` tag with optional integrity attribute."""
    cache = _load_sri_cache()
    digest = cache.get(sri_key)
    if digest:
        return (
            f'<script src="{url}" integrity="sha384-{digest}" '
            'crossorigin="anonymous"></script>'
        )
    return f'<script src="{url}"></script>'


def get_plotly_script_tag(mode: str) -> str:
    """Return a ``<script>`` block providing Plotly.

    Args:
        mode: ``"offline"`` to inline the JS shipped with the plotly
            Python package, ``"cdn"`` to emit a remote ``<script src>``.
            CDN mode attaches an ``integrity`` hash when one has been
            cached via ``metadarkmatter report vendor-cdn-sri``.

    Returns:
        HTML ``<script>`` block (single tag) suitable for inclusion in the
        report ``<head>``.
    """
    if mode == "cdn":
        return "    " + _cdn_script_tag(PLOTLY_CDN_URL, "plotly")

    try:
        from plotly.offline import get_plotlyjs
    except ImportError:  # pragma: no cover - plotly is a hard dep
        logger.warning(
            "plotly.offline.get_plotlyjs unavailable; falling back to CDN.",
        )
        return "    " + _cdn_script_tag(PLOTLY_CDN_URL, "plotly")

    return f"<script>\n{get_plotlyjs()}\n</script>"


def get_d3_script_tag(mode: str) -> str:
    """Return a ``<script>`` block providing D3 v7.

    Args:
        mode: ``"offline"`` to inline a vendored copy of D3 if present
            (with CDN fallback when not vendored), ``"cdn"`` to always emit
            the CDN tag. CDN mode attaches an ``integrity`` hash when one
            has been cached via ``metadarkmatter report vendor-cdn-sri``.

    Returns:
        HTML ``<script>`` block.
    """
    if mode == "cdn":
        return _cdn_script_tag(D3_CDN_URL, "d3")

    vendored = _ASSETS_DIR / D3_FILENAME
    if vendored.is_file():
        return f"<script>\n{vendored.read_text(encoding='utf-8')}\n</script>"

    logger.info(
        "D3 not vendored at %s; falling back to CDN. "
        "Run 'metadarkmatter report vendor-d3' to enable offline phylogeny.",
        vendored,
    )
    return _cdn_script_tag(D3_CDN_URL, "d3")
