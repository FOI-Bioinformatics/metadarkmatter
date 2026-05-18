"""Tests for the on-disk response cache in GTDBClient."""

from __future__ import annotations

import json
import time
from pathlib import Path
from unittest.mock import MagicMock

import httpx
import pytest

from metadarkmatter.clients.gtdb import (
    USE_DEFAULT_CACHE,
    GTDBClient,
    _default_cache_dir,
)


@pytest.fixture
def cached_client(tmp_path: Path) -> GTDBClient:
    """Return a GTDBClient pointed at an isolated cache directory."""
    return GTDBClient(
        max_retries=1,
        retry_delay=0.01,
        cache_dir=tmp_path / "cache",
        cache_ttl=60.0,
    )


def _mock_http_client_returning(payload: dict) -> MagicMock:
    response = MagicMock()
    response.json.return_value = payload
    response.raise_for_status = MagicMock()

    client = MagicMock(spec=httpx.Client)
    client.get.return_value = response
    return client


def test_cache_hit_skips_http(cached_client: GTDBClient) -> None:
    cached_client._client = _mock_http_client_returning({"rows": [{"gid": "GCF_001"}]})

    first = cached_client._get("/taxon/X/genomes-detail", {"sp_reps_only": "true"})
    second = cached_client._get("/taxon/X/genomes-detail", {"sp_reps_only": "true"})

    assert first == {"rows": [{"gid": "GCF_001"}]}
    assert first == second
    # Second call should have hit the cache, not the network.
    assert cached_client._client.get.call_count == 1


def test_cache_disabled_by_default() -> None:
    client = GTDBClient(max_retries=1, retry_delay=0.01)
    assert client.cache_dir is None


def test_explicit_cache_dir_is_used(tmp_path: Path) -> None:
    client = GTDBClient(cache_dir=tmp_path / "x")
    assert client.cache_dir == tmp_path / "x"


def test_use_default_cache_sentinel_resolves(monkeypatch: pytest.MonkeyPatch) -> None:
    # Force the env-var path: empty string disables, an explicit path overrides.
    monkeypatch.setenv("METADARKMATTER_GTDB_CACHE_DIR", "/tmp/mdm_cache_test")
    client = GTDBClient(cache_dir=USE_DEFAULT_CACHE)
    assert client.cache_dir == Path("/tmp/mdm_cache_test")


def test_env_var_empty_string_disables_default_cache(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("METADARKMATTER_GTDB_CACHE_DIR", "")
    assert _default_cache_dir() is None


def test_expired_cache_is_refreshed(tmp_path: Path) -> None:
    client = GTDBClient(
        max_retries=1,
        retry_delay=0.01,
        cache_dir=tmp_path / "cache",
        cache_ttl=0.01,  # essentially immediate expiry
    )
    client._client = _mock_http_client_returning({"rows": []})

    client._get("/taxon/Y/genomes-detail", {"sp_reps_only": "true"})
    time.sleep(0.05)
    client._get("/taxon/Y/genomes-detail", {"sp_reps_only": "true"})

    # Two HTTP calls because the first cache entry expired.
    assert client._client.get.call_count == 2


def test_cache_file_is_written_with_metadata(cached_client: GTDBClient) -> None:
    cached_client._client = _mock_http_client_returning({"rows": [{"gid": "G1"}]})

    cached_client._get("/taxon/Z/genomes-detail", {"sp_reps_only": "true"})

    cache_files = list(cached_client.cache_dir.glob("*.json"))
    assert len(cache_files) == 1
    payload = json.loads(cache_files[0].read_text())
    assert payload["endpoint"] == "/taxon/Z/genomes-detail"
    assert payload["params"] == {"sp_reps_only": "true"}
    assert payload["data"] == {"rows": [{"gid": "G1"}]}
    assert "timestamp" in payload


def test_corrupt_cache_file_is_ignored(cached_client: GTDBClient) -> None:
    cached_client._client = _mock_http_client_returning({"rows": []})

    # Pre-populate the cache file with garbage so the read returns None.
    cache_path = cached_client._cache_path("/taxon/Q/genomes-detail", {"a": "b"})
    assert cache_path is not None
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text("not json {{")

    result = cached_client._get("/taxon/Q/genomes-detail", {"a": "b"})
    assert result == {"rows": []}
    assert cached_client._client.get.call_count == 1
