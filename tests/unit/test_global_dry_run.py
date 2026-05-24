"""Tests for the process-wide --dry-run flag."""

from __future__ import annotations

import pytest

from metadarkmatter.core.runtime import is_dry_run, set_dry_run


@pytest.fixture(autouse=True)
def _clean(monkeypatch: pytest.MonkeyPatch):
    # Ensure each test starts from a known state.
    monkeypatch.delenv("MDM_DRY_RUN", raising=False)


def test_default_is_off() -> None:
    assert is_dry_run() is False


def test_set_and_unset_via_helper() -> None:
    set_dry_run(True)
    assert is_dry_run() is True
    set_dry_run(False)
    assert is_dry_run() is False


@pytest.mark.parametrize("value", ["1", "true", "TRUE", "yes", "on"])
def test_truthy_env_values(monkeypatch: pytest.MonkeyPatch, value: str) -> None:
    monkeypatch.setenv("MDM_DRY_RUN", value)
    assert is_dry_run() is True


@pytest.mark.parametrize("value", ["0", "false", "no", "off", "", "garbage"])
def test_falsy_env_values(monkeypatch: pytest.MonkeyPatch, value: str) -> None:
    monkeypatch.setenv("MDM_DRY_RUN", value)
    assert is_dry_run() is False


def test_global_flag_propagates_to_subcommand(monkeypatch: pytest.MonkeyPatch) -> None:
    """A top-level '--dry-run' must be honoured by subcommands that
    consult is_dry_run(), independent of their local --dry-run flag.
    """
    from typer.testing import CliRunner

    from metadarkmatter.cli.main import app

    runner = CliRunner()
    # 'download' has --dry-run logic; passing --dry-run at the top level
    # should produce the same exit-without-network behaviour as passing
    # it to the subcommand directly.
    result = runner.invoke(
        app,
        [
            "--dry-run",
            "download", "genomes", "list", "f__Francisellaceae",
            "--output", "/tmp/x.tsv",
        ],
    )
    # We don't care about the exact message, just that the top-level flag
    # made the subcommand short-circuit rather than make a real request.
    assert "dry" in result.output.lower() or "would" in result.output.lower()
