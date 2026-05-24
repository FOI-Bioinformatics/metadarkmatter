"""Tests for SRI integrity attribute injection on CDN script tags."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import metadarkmatter.visualization.report.assets as assets
from metadarkmatter.visualization.report.assets import (
    SRI_CACHE_FILENAME,
    get_d3_script_tag,
    get_plotly_script_tag,
)


@pytest.fixture
def sri_cache(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Point the assets module at a temp directory and return the cache path."""
    monkeypatch.setattr(assets, "_ASSETS_DIR", tmp_path)
    return tmp_path / SRI_CACHE_FILENAME


def test_cdn_emits_bare_tag_without_cache(sri_cache: Path) -> None:
    # No cache file exists yet.
    tag = get_plotly_script_tag("cdn")
    assert "integrity=" not in tag
    assert "https://cdn.plot.ly" in tag


def test_cdn_emits_integrity_when_cache_present(sri_cache: Path) -> None:
    sri_cache.write_text(json.dumps({"plotly": "AAA", "d3": "BBB"}))
    plotly_tag = get_plotly_script_tag("cdn")
    d3_tag = get_d3_script_tag("cdn")
    assert 'integrity="sha384-AAA"' in plotly_tag
    assert 'crossorigin="anonymous"' in plotly_tag
    assert 'integrity="sha384-BBB"' in d3_tag
    assert 'crossorigin="anonymous"' in d3_tag


def test_malformed_cache_falls_back_to_bare_tag(sri_cache: Path) -> None:
    sri_cache.write_text("not json {{")
    tag = get_d3_script_tag("cdn")
    assert "integrity=" not in tag


def test_partial_cache_still_works(sri_cache: Path) -> None:
    # Cache only knows about Plotly; D3 should emit bare.
    sri_cache.write_text(json.dumps({"plotly": "AAA"}))
    plotly_tag = get_plotly_script_tag("cdn")
    d3_tag = get_d3_script_tag("cdn")
    assert 'integrity="sha384-AAA"' in plotly_tag
    assert "integrity=" not in d3_tag


def test_offline_mode_ignores_sri_cache(sri_cache: Path) -> None:
    sri_cache.write_text(json.dumps({"plotly": "AAA"}))
    tag = get_plotly_script_tag("offline")
    # Offline embeds the JS body inline; no remote script tag at all.
    assert "<script>" in tag
    assert "integrity=" not in tag
