"""
Central source of truth for random seeds.

Every randomness-using call site in metadarkmatter resolves its seed
through :func:`get_seed`, so a single environment variable
(``METADARKMATTER_SEED``) controls every stochastic step in a run.
This makes whole-pipeline determinism testable and lets users probe
sensitivity to seed choice without editing source.
"""

from __future__ import annotations

import os

DEFAULT_SEED = 42
_ENV_VAR = "METADARKMATTER_SEED"


def get_seed() -> int:
    """Return the active random seed.

    Reads ``METADARKMATTER_SEED`` from the environment (parsed as an
    integer) and falls back to :data:`DEFAULT_SEED` when the variable is
    unset, empty, or unparseable. Pure function with no side effects.
    """
    raw = os.environ.get(_ENV_VAR)
    if raw is None or raw.strip() == "":
        return DEFAULT_SEED
    try:
        return int(raw)
    except ValueError:
        return DEFAULT_SEED
