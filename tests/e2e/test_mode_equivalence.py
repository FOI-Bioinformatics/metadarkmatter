"""
E2E mode equivalence tests for metadarkmatter CLI.

Tests that the default (vectorized) and streaming processing modes
produce identical or functionally equivalent results.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable

import polars as pl
import pytest

from tests.utils.assertions import CLIAssertions


pytestmark = pytest.mark.e2e


class TestModeEquivalence:
    """
    Tests verifying that default and streaming modes produce identical results.

    Both modes should produce the same classification output for the same input,
    with only performance characteristics differing.
    """

    def test_default_vs_streaming_equivalence(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """Default and streaming modes should produce identical output."""
        blast_path, ani_path = standard_dataset

        # Run default mode (VectorizedClassifier)
        default_result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="default.csv",
        )
        CLIAssertions.assert_success(default_result)

        # Run streaming mode
        streaming_result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="streaming.csv",
            streaming=True,
        )
        CLIAssertions.assert_success(streaming_result)

        # Compare outputs
        default_df = pl.read_csv(e2e_temp_dir / "default.csv")
        streaming_df = pl.read_csv(e2e_temp_dir / "streaming.csv")

        assert len(default_df) == len(streaming_df), (
            f"Row count mismatch: default={len(default_df)}, streaming={len(streaming_df)}"
        )

        # Sort and compare
        default_sorted = default_df.sort("read_id")
        streaming_sorted = streaming_df.sort("read_id")

        for col in ["read_id", "taxonomic_call", "is_novel"]:
            assert default_sorted[col].to_list() == streaming_sorted[col].to_list(), (
                f"Column {col} differs between default and streaming modes"
            )

    def test_both_modes_same_summary_stats(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """Both modes should produce the same summary statistics."""
        import json

        blast_path, ani_path = standard_dataset
        modes = [
            ("default", {}),
            ("streaming", {"streaming": True}),
        ]

        summaries = {}

        for mode_name, mode_kwargs in modes:
            result = run_classify(
                blast=blast_path,
                ani=ani_path,
                output_name=f"{mode_name}_summary.csv",
                summary_name=f"{mode_name}_summary.json",
                **mode_kwargs,
            )
            CLIAssertions.assert_success(result)

            summary_path = e2e_temp_dir / f"{mode_name}_summary.json"
            with open(summary_path) as f:
                summaries[mode_name] = json.load(f)

        # Compare streaming to default
        default = summaries["default"]
        streaming = summaries["streaming"]

        assert streaming["total_reads"] == default["total_reads"], (
            "streaming total_reads differs from default"
        )
        assert streaming["known_species"] == default["known_species"], (
            "streaming known_species differs from default"
        )
        assert streaming["novel_species"] == default["novel_species"], (
            "streaming novel_species differs from default"
        )
        assert streaming["novel_genus"] == default["novel_genus"], (
            "streaming novel_genus differs from default"
        )
        assert streaming["conserved_regions"] == default["conserved_regions"], (
            "streaming conserved_regions differs from default"
        )
