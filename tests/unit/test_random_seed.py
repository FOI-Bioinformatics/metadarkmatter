"""Tests for the centralised random-seed helper and pipeline determinism."""

from __future__ import annotations

import polars as pl
import pytest

from metadarkmatter.core.random import DEFAULT_SEED, get_seed
from metadarkmatter.visualization.plots.base import subsample_dataframe

# ---------------------------------------------------------------------------
# get_seed()
# ---------------------------------------------------------------------------


def test_default_seed(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("METADARKMATTER_SEED", raising=False)
    assert get_seed() == DEFAULT_SEED


def test_env_override(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("METADARKMATTER_SEED", "12345")
    assert get_seed() == 12345


def test_empty_env_uses_default(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("METADARKMATTER_SEED", "   ")
    assert get_seed() == DEFAULT_SEED


def test_unparseable_env_uses_default(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("METADARKMATTER_SEED", "not-an-int")
    assert get_seed() == DEFAULT_SEED


# ---------------------------------------------------------------------------
# Determinism: identical seed produces identical output
# ---------------------------------------------------------------------------


@pytest.fixture
def big_df() -> pl.DataFrame:
    # Big enough to trigger subsampling and to have enough rows for the
    # seed to produce visibly different orderings.
    return pl.DataFrame(
        {
            "x": list(range(1000)),
            "y": [i * 0.1 for i in range(1000)],
        }
    )


def test_subsample_is_deterministic_with_same_seed(
    big_df: pl.DataFrame, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("METADARKMATTER_SEED", "7")
    out1 = subsample_dataframe(big_df, max_points=50, seed=get_seed())
    out2 = subsample_dataframe(big_df, max_points=50, seed=get_seed())
    assert out1.equals(out2)


def test_subsample_differs_with_different_seed(
    big_df: pl.DataFrame, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("METADARKMATTER_SEED", "7")
    out_a = subsample_dataframe(big_df, max_points=50, seed=get_seed())
    monkeypatch.setenv("METADARKMATTER_SEED", "99")
    out_b = subsample_dataframe(big_df, max_points=50, seed=get_seed())
    # With 1000 rows and 50 samples, two distinct seeds should yield
    # distinguishable subsets virtually always.
    assert not out_a.equals(out_b)


# ---------------------------------------------------------------------------
# Sanity check that the call sites have actually been threaded
# ---------------------------------------------------------------------------


def test_recruitment_plots_import_routes_through_get_seed() -> None:
    import metadarkmatter.visualization.recruitment_plots as mod

    src = mod.__file__
    with open(src, encoding="utf-8") as fh:
        text = fh.read()
    # Hardcoded literal seed=42 should not appear anywhere here.
    assert "seed=42" not in text
    assert "get_seed()" in text


def test_adaptive_module_uses_get_seed() -> None:
    import metadarkmatter.core.classification.adaptive as mod

    with open(mod.__file__, encoding="utf-8") as fh:
        text = fh.read()
    assert "random_state=42" not in text
    assert "random_state=get_seed()" in text
