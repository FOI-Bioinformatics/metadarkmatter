"""
Shared CLI option definitions.

Centralizes the option types that recur across the ~18 subcommands so that
their flag spelling, help text, and validation live in one place. Commands
consume these via ``Annotated`` parameters; the *default value* still lives
at each command's call site (a Typer requirement), but the flag definition
is shared::

    from metadarkmatter.cli.options import Quiet, Threads

    def compute(..., threads: Threads = 4, quiet: Quiet = False) -> None:
        ...

This keeps ``--quiet``/``--verbose``/``--dry-run``/``--threads`` consistent
across commands and makes a help-text change a single edit.
"""

from __future__ import annotations

from typing import Annotated

import typer

Quiet = Annotated[
    bool,
    typer.Option("--quiet", "-q", help="Suppress progress output (useful for scripting)."),
]

Verbose = Annotated[
    bool,
    typer.Option("--verbose", help="Print additional informational output."),
]

DryRun = Annotated[
    bool,
    typer.Option(
        "--dry-run",
        help="Show what would run without executing. Honours the global --dry-run too.",
    ),
]

Threads = Annotated[
    int,
    typer.Option("--threads", "-t", min=1, help="Number of threads for parallel execution."),
]
