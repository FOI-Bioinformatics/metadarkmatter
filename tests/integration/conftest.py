"""Shared conftest for integration tests."""

from __future__ import annotations

import pytest

from metadarkmatter.external.mmseqs2 import MMseqs2


def pytest_collection_modifyitems(config, items):
    """Skip tests marked requires_mmseqs2 if MMseqs2 not available."""
    if not MMseqs2.check_available():
        skip_mmseqs2 = pytest.mark.skip(reason="MMseqs2 not installed")
        for item in items:
            if "requires_mmseqs2" in item.keywords:
                item.add_marker(skip_mmseqs2)
