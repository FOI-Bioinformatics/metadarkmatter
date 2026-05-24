"""
Process-wide runtime flags.

A small home for cross-cutting boolean state that needs to flow from
the top-level CLI down into individual subcommands without threading
explicit arguments through every call site. Currently exposes a
process-wide ``--dry-run`` flag backed by the ``MDM_DRY_RUN``
environment variable.
"""

from __future__ import annotations

import os

_DRY_RUN_ENV_VAR = "MDM_DRY_RUN"


def set_dry_run(enabled: bool) -> None:
    """Toggle the process-wide dry-run flag.

    Implemented as an env var so child processes (subprocess.run calls
    in the external-tool wrappers) inherit the setting automatically.
    """
    if enabled:
        os.environ[_DRY_RUN_ENV_VAR] = "1"
    else:
        os.environ.pop(_DRY_RUN_ENV_VAR, None)


def is_dry_run() -> bool:
    """Return True when the global dry-run flag is set.

    Honours ``MDM_DRY_RUN=1``. Per-command ``--dry-run`` flags should
    OR with this value so the global flag is always honoured even
    when a subcommand sets its own boolean.
    """
    raw = os.environ.get(_DRY_RUN_ENV_VAR)
    if raw is None:
        return False
    return raw.strip().lower() in {"1", "true", "yes", "on"}
