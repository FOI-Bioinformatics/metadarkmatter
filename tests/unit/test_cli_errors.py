"""Tests for the centralized CLI error handler (cli/errors.py)."""

from __future__ import annotations

import typer
from typer.testing import CliRunner

from metadarkmatter.cli.errors import handle_cli_errors, wrap_app_commands
from metadarkmatter.core.exceptions import GenomeCoverageWarning, MetadarkmatterError

runner = CliRunner()


def _make_app() -> typer.Typer:
    app = typer.Typer()

    @app.command()
    def boom_mdm() -> None:
        raise MetadarkmatterError("something broke", suggestion="try this instead")

    @app.command()
    def boom_value() -> None:
        raise ValueError("not square")

    @app.command()
    def boom_unexpected() -> None:
        raise RuntimeError("kaboom")

    @app.command()
    def boom_warning() -> None:
        raise GenomeCoverageWarning(1, 10, 10.0, ["GCF_1"])

    @app.command()
    def ok() -> None:
        typer.echo("done")

    wrap_app_commands(app)
    return app


def test_metadarkmatter_error_shows_message_and_suggestion():
    result = runner.invoke(_make_app(), ["boom-mdm"])
    assert result.exit_code == 1
    assert "something broke" in result.output
    assert "try this instead" in result.output
    # No raw traceback without --debug.
    assert "Traceback" not in result.output


def test_value_error_is_friendly_invalid_input():
    result = runner.invoke(_make_app(), ["boom-value"])
    assert result.exit_code == 1
    assert "Invalid input" in result.output
    assert "not square" in result.output


def test_unexpected_error_is_caught():
    result = runner.invoke(_make_app(), ["boom-unexpected"])
    assert result.exit_code == 1
    assert "Unexpected error" in result.output
    assert "kaboom" in result.output


def test_warning_error_exits_zero_via_exit_code():
    result = runner.invoke(_make_app(), ["boom-warning"])
    # GenomeCoverageWarning carries exit_code = 0 (advisory, not a hard failure).
    assert result.exit_code == 0


def test_successful_command_passes_through():
    result = runner.invoke(_make_app(), ["ok"])
    assert result.exit_code == 0
    assert "done" in result.output


def test_debug_env_emits_traceback(monkeypatch):
    monkeypatch.setenv("MDM_DEBUG", "1")
    result = runner.invoke(_make_app(), ["boom-unexpected"])
    assert result.exit_code == 1
    assert "Traceback" in result.output


def test_decorator_preserves_signature():
    @handle_cli_errors
    def example(x: int, y: str = "z") -> None:
        """Docstring kept."""

    assert example.__name__ == "example"
    assert example.__doc__ == "Docstring kept."
