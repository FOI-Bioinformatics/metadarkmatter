"""
Central source of truth for random seeds.

Every randomness-using call site in metadarkmatter resolves its seed
through :func:`get_seed`, so a single environment variable
(``METADARKMATTER_SEED``) controls every stochastic step in a run.
This makes whole-pipeline determinism testable and lets users probe
sensitivity to seed choice without editing source.
"""

from __future__ import annotations

import logging
import os

DEFAULT_SEED = 42
_ENV_VAR = "METADARKMATTER_SEED"
# numpy, scikit-learn and polars accept seeds in [0, 2**32 - 1].
_MAX_SEED = 2**32 - 1

logger = logging.getLogger(__name__)


def get_seed() -> int:
    """Return the active random seed.

    Reads ``METADARKMATTER_SEED`` from the environment (parsed as an
    integer) and falls back to :data:`DEFAULT_SEED` when the variable is
    unset, empty, unparseable, or outside 0..2**32-1 (the range numpy,
    scikit-learn and polars accept; a warning is logged in that case).
    """
    raw = os.environ.get(_ENV_VAR)
    if raw is None or raw.strip() == "":
        return DEFAULT_SEED
    try:
        seed = int(raw)
    except ValueError:
        return DEFAULT_SEED
    if not 0 <= seed <= _MAX_SEED:
        logger.warning(
            "%s=%s is outside 0..%d; using the default seed %d.",
            _ENV_VAR, raw, _MAX_SEED, DEFAULT_SEED,
        )
        return DEFAULT_SEED
    return seed
