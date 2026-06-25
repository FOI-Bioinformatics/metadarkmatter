"""Tests for build_classify_config (config precedence: --config > flags > preset > defaults)."""

from __future__ import annotations

from pathlib import Path

import pytest

from metadarkmatter.cli.score import (
    THRESHOLD_PRESETS,
    _ClassifyOptions,
    build_classify_config,
)
from metadarkmatter.core.exceptions import ConfigurationError


def _opts(**overrides: object) -> _ClassifyOptions:
    """An options bundle matching the CLI defaults, with selective overrides."""
    defaults: dict[str, object] = {
        "alignment_mode": "nucleotide",
        "bitscore_threshold": 95.0,
        "min_alignment_length": 100,
        "min_alignment_fraction": 0.3,
        "max_evalue": 0.0,
        "min_percent_identity": 0.0,
        "min_bitscore": 0.0,
        "min_read_length": 0,
        "min_query_coverage": 0.0,
        "coverage_weight_mode": "linear",
        "coverage_weight_strength": 0.5,
        "uncertainty_mode": "second",
        "single_hit_uncertainty_threshold": 10.0,
        "target_family": None,
        "family_ratio_threshold": 0.8,
        "include_legacy_scores": False,
        "bayesian": False,
        "config_file": None,
        "preset": None,
    }
    defaults.update(overrides)
    return _ClassifyOptions(**defaults)  # type: ignore[arg-type]


def test_defaults_produce_default_config():
    config, messages = build_classify_config(_opts())
    assert config.alignment_mode == "nucleotide"
    assert config.bitscore_threshold_pct == 95.0
    assert messages == []


def test_flag_overrides_without_preset():
    config, _ = build_classify_config(_opts(min_alignment_length=150, max_evalue=1e-5))
    assert config.min_alignment_length == 150
    assert config.max_evalue == 1e-5


def test_preset_applies_and_messages():
    config, messages = build_classify_config(_opts(preset="gtdb-strict"))
    assert config == THRESHOLD_PRESETS["gtdb-strict"]
    assert any("Using preset: gtdb-strict" in m for m in messages)


def test_preset_case_insensitive():
    config, _ = build_classify_config(_opts(preset="GTDB-STRICT"))
    assert config == THRESHOLD_PRESETS["gtdb-strict"]


def test_unknown_preset_raises_configuration_error():
    with pytest.raises(ConfigurationError, match="Unknown preset"):
        build_classify_config(_opts(preset="does-not-exist"))


def test_preset_with_explicit_overrides_rebuilds():
    # An explicit non-default mode triggers the preset-override rebuild path.
    config, _ = build_classify_config(_opts(preset="gtdb-strict", alignment_mode="protein"))
    assert config.alignment_mode == "protein"
    # preset thresholds are preserved
    assert config.novelty_known_max == THRESHOLD_PRESETS["gtdb-strict"].novelty_known_max


def test_config_file_loads_and_flags_override(tmp_path: Path):
    from metadarkmatter.models.config import ScoringConfig

    yaml_path = tmp_path / "cfg.yaml"
    ScoringConfig(min_alignment_length=200).to_yaml(yaml_path)

    # No override -> YAML value wins
    config, messages = build_classify_config(_opts(config_file=yaml_path))
    assert config.min_alignment_length == 200
    assert any("Loaded config from" in m for m in messages)

    # CLI flag overrides YAML, and a deprecation warning is emitted
    config2, messages2 = build_classify_config(
        _opts(config_file=yaml_path, min_alignment_length=150)
    )
    assert config2.min_alignment_length == 150
    assert any("deprecated when using --config" in m for m in messages2)


def test_bayesian_flag_emits_deprecation():
    _, messages = build_classify_config(_opts(bayesian=True))
    assert any("--bayesian is deprecated" in m for m in messages)
