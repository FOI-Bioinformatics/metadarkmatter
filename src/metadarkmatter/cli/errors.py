"""
Centralized CLI error handling.

Provides a single, signature-preserving decorator (:func:`handle_cli_errors`)
that wraps every Typer command so that:

- :class:`~metadarkmatter.core.exceptions.MetadarkmatterError` subclasses
  (including the ``external`` tool errors, which already inherit from it)
  print their ``message`` and ``suggestion`` cleanly and exit with the
  error's ``exit_code``.
- Common library/runtime errors (Polars, file system, memory) are
  translated into a friendly message plus an actionable suggestion.
- Full tracebacks are printed only when the global ``--debug`` flag is set
  (``MDM_DEBUG=1``); otherwise users see a one-line error, not a stack trace.

The decorator is applied automatically to every registered command by
:func:`wrap_app_commands`, called once per sub-app in ``cli/main.py``,
so individual command modules need no per-command changes.
"""

from __future__ import annotations

import functools
from collections.abc import Callable
from typing import TypeVar

import typer
from rich.console import Console

from metadarkmatter.core.exceptions import MetadarkmatterError
from metadarkmatter.core.runtime import is_debug

_err_console = Console(stderr=True)

F = TypeVar("F", bound=Callable[..., object])


def _print_error(message: str, suggestion: str | None = None) -> None:
    """Render an error message (and optional suggestion) to stderr."""
    _err_console.print(f"[bold red]Error:[/bold red] {message}")
    if suggestion:
        _err_console.print(f"\n[yellow]Suggestion:[/yellow] {suggestion}")


def _maybe_traceback() -> None:
    """Print a full traceback only when --debug / MDM_DEBUG is enabled."""
    if is_debug():
        _err_console.print_exception()


def _exit(code: int) -> typer.Exit:
    return typer.Exit(code=code)


def handle_cli_errors(func: F) -> F:
    """Wrap a Typer command callback with friendly error translation.

    Signature-preserving (via :func:`functools.wraps`) so Typer can still
    introspect the wrapped function to build its options. ``typer.Exit`` and
    keyboard interrupts pass through untouched.
    """

    @functools.wraps(func)
    def wrapper(*args: object, **kwargs: object) -> object:
        try:
            return func(*args, **kwargs)
        except typer.Exit:
            # Intended control-flow exits (including --version, dry-run) pass through.
            raise
        except KeyboardInterrupt:
            _err_console.print("\n[yellow]Interrupted.[/yellow]")
            raise _exit(130) from None
        except MetadarkmatterError as exc:
            # Covers ANIMatrixError, BlastFileError, Tool*Error, GenomeCoverageWarning, etc.
            _print_error(exc.message, exc.suggestion)
            _maybe_traceback()
            raise _exit(exc.exit_code) from None
        except FileNotFoundError as exc:
            target = exc.filename or str(exc)
            _print_error(
                f"File not found: {target}",
                "Check the path is correct and the file exists. "
                "Relative paths are resolved from the current working directory.",
            )
            _maybe_traceback()
            raise _exit(1) from None
        except PermissionError as exc:
            target = exc.filename or str(exc)
            _print_error(
                f"Permission denied: {target}",
                "Check file/directory permissions and that the output location is writable.",
            )
            _maybe_traceback()
            raise _exit(1) from None
        except MemoryError:
            _print_error(
                "Ran out of memory while processing.",
                "For large alignment files, use '--streaming --chunk-size 500000' "
                "to keep memory bounded.",
            )
            _maybe_traceback()
            raise _exit(1) from None
        except ValueError as exc:
            # At the CLI boundary a ValueError is almost always invalid input
            # or configuration (e.g. a malformed matrix), not an internal crash.
            _print_error(
                f"Invalid input: {exc}",
                "Check the input files and options. 'mdm validate' can pre-check "
                "alignment, ANI matrix, and classification files.",
            )
            _maybe_traceback()
            raise _exit(1) from None
        except Exception as exc:  # last-resort CLI boundary
            if _is_polars_error(exc):
                _print_error(
                    f"Data processing error: {exc}",
                    "This usually means a malformed alignment/matrix file or an "
                    "unexpected column layout. Validate inputs with 'mdm validate'.",
                )
                _maybe_traceback()
                raise _exit(1) from None
            if isinstance(exc, OSError):
                _print_error(
                    f"I/O error: {exc}",
                    "Check disk space, the path, and file permissions.",
                )
                _maybe_traceback()
                raise _exit(1) from None
            _print_error(
                f"Unexpected error: {exc}",
                "Re-run with --debug for a full traceback. If this persists, "
                "please report it with the traceback and the command you ran.",
            )
            _maybe_traceback()
            raise _exit(1) from None

    return wrapper  # type: ignore[return-value]


def _is_polars_error(exc: BaseException) -> bool:
    """Return True if ``exc`` is a Polars error, without importing polars eagerly."""
    try:
        import polars as pl
    except ImportError:  # pragma: no cover - polars is a hard dependency
        return False
    return isinstance(exc, pl.exceptions.PolarsError)


def wrap_app_commands(app: typer.Typer) -> None:
    """Wrap every command registered on ``app`` with :func:`handle_cli_errors`.

    Iterates ``app.registered_commands`` and replaces each command's callback
    with the wrapped version. Idempotent: already-wrapped callbacks (marked
    with ``_mdm_error_wrapped``) are skipped so repeated calls are safe.
    """
    for command in app.registered_commands:
        callback = command.callback
        if callback is None or getattr(callback, "_mdm_error_wrapped", False):
            continue
        wrapped = handle_cli_errors(callback)
        wrapped._mdm_error_wrapped = True
        command.callback = wrapped
