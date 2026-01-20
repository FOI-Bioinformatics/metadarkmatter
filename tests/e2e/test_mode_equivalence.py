"""
E2E mode equivalence tests for metadarkmatter CLI.

Tests that all processing modes (standard, fast, parallel, streaming)
produce identical or functionally equivalent results.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable

import polars as pl
import pytest

from tests.utils.assertions import CLIAssertions, ClassificationAssertions


pytestmark = pytest.mark.e2e


class TestModeEquivalence:
    """
    Tests verifying that different processing modes produce identical results.

    All modes should produce the same classification output for the same input,
    with only performance characteristics differing.
    """

    def test_standard_vs_fast_equivalence(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """Standard and fast modes should produce identical output."""
        blast_path, ani_path = standard_dataset

        # Run standard mode
        standard_result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="standard.csv",
        )
        CLIAssertions.assert_success(standard_result)

        # Run fast mode
        fast_result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="fast.csv",
            fast=True,
        )
        CLIAssertions.assert_success(fast_result)

        # Compare outputs
        standard_df = pl.read_csv(e2e_temp_dir / "standard.csv")
        fast_df = pl.read_csv(e2e_temp_dir / "fast.csv")

        # Should have same number of rows
        assert len(standard_df) == len(fast_df), (
            f"Row count mismatch: standard={len(standard_df)}, fast={len(fast_df)}"
        )

        # Sort by read_id for comparison
        standard_sorted = standard_df.sort("read_id")
        fast_sorted = fast_df.sort("read_id")

        # Compare key columns
        for col in ["read_id", "taxonomic_call", "is_novel"]:
            assert standard_sorted[col].to_list() == fast_sorted[col].to_list(), (
                f"Column {col} differs between standard and fast modes"
            )

        # Numeric columns should be very close (floating point tolerance)
        for col in ["novelty_index", "placement_uncertainty"]:
            standard_vals = standard_sorted[col].to_list()
            fast_vals = fast_sorted[col].to_list()
            for i, (s, f) in enumerate(zip(standard_vals, fast_vals, strict=False)):
                assert abs(s - f) < 0.01, (
                    f"Column {col} row {i} differs: standard={s}, fast={f}"
                )

    def test_standard_vs_parallel_equivalence(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """Standard and parallel modes should produce identical output."""
        blast_path, ani_path = standard_dataset

        # Run standard mode
        standard_result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="std_for_parallel.csv",
        )
        CLIAssertions.assert_success(standard_result)

        # Run parallel mode
        parallel_result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="parallel.csv",
            parallel=True,
        )
        CLIAssertions.assert_success(parallel_result)

        # Compare outputs
        standard_df = pl.read_csv(e2e_temp_dir / "std_for_parallel.csv")
        parallel_df = pl.read_csv(e2e_temp_dir / "parallel.csv")

        assert len(standard_df) == len(parallel_df), (
            f"Row count mismatch: standard={len(standard_df)}, parallel={len(parallel_df)}"
        )

        # Sort and compare
        standard_sorted = standard_df.sort("read_id")
        parallel_sorted = parallel_df.sort("read_id")

        for col in ["read_id", "taxonomic_call", "is_novel"]:
            assert standard_sorted[col].to_list() == parallel_sorted[col].to_list(), (
                f"Column {col} differs between standard and parallel modes"
            )

    def test_standard_vs_streaming_equivalence(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """Standard and streaming modes should produce identical output."""
        blast_path, ani_path = standard_dataset

        # Run standard mode
        standard_result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="std_for_streaming.csv",
        )
        CLIAssertions.assert_success(standard_result)

        # Run streaming mode
        streaming_result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="streaming.csv",
            streaming=True,
        )
        CLIAssertions.assert_success(streaming_result)

        # Compare outputs
        standard_df = pl.read_csv(e2e_temp_dir / "std_for_streaming.csv")
        streaming_df = pl.read_csv(e2e_temp_dir / "streaming.csv")

        assert len(standard_df) == len(streaming_df), (
            f"Row count mismatch: standard={len(standard_df)}, streaming={len(streaming_df)}"
        )

        # Sort and compare
        standard_sorted = standard_df.sort("read_id")
        streaming_sorted = streaming_df.sort("read_id")

        for col in ["read_id", "taxonomic_call", "is_novel"]:
            assert standard_sorted[col].to_list() == streaming_sorted[col].to_list(), (
                f"Column {col} differs between standard and streaming modes"
            )

    def test_fast_vs_parallel_equivalence(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """Fast and parallel modes should produce identical output."""
        blast_path, ani_path = standard_dataset

        # Run fast mode
        fast_result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="fast_for_parallel.csv",
            fast=True,
        )
        CLIAssertions.assert_success(fast_result)

        # Run parallel mode
        parallel_result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="parallel_for_fast.csv",
            parallel=True,
        )
        CLIAssertions.assert_success(parallel_result)

        # Compare outputs
        fast_df = pl.read_csv(e2e_temp_dir / "fast_for_parallel.csv")
        parallel_df = pl.read_csv(e2e_temp_dir / "parallel_for_fast.csv")

        assert len(fast_df) == len(parallel_df)

        # Sort and compare taxonomic calls
        fast_sorted = fast_df.sort("read_id")
        parallel_sorted = parallel_df.sort("read_id")

        assert fast_sorted["taxonomic_call"].to_list() == parallel_sorted["taxonomic_call"].to_list()

    def test_all_modes_same_summary_stats(
        self,
        standard_dataset: tuple[Path, Path],
        run_classify: Callable,
        e2e_temp_dir: Path,
    ):
        """All modes should produce the same summary statistics."""
        blast_path, ani_path = standard_dataset
        modes = [
            ("standard", {}),
            ("fast", {"fast": True}),
            ("parallel", {"parallel": True}),
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

            # Load summary
            import json
            summary_path = e2e_temp_dir / f"{mode_name}_summary.json"
            with open(summary_path) as f:
                summaries[mode_name] = json.load(f)

        # Compare all summaries to standard
        standard = summaries["standard"]
        for mode_name, summary in summaries.items():
            if mode_name == "standard":
                continue

            assert summary["total_reads"] == standard["total_reads"], (
                f"{mode_name} total_reads differs from standard"
            )
            assert summary["known_species"] == standard["known_species"], (
                f"{mode_name} known_species differs from standard"
            )
            assert summary["novel_species"] == standard["novel_species"], (
                f"{mode_name} novel_species differs from standard"
            )
            assert summary["novel_genus"] == standard["novel_genus"], (
                f"{mode_name} novel_genus differs from standard"
            )
            assert summary["conserved_regions"] == standard["conserved_regions"], (
                f"{mode_name} conserved_regions differs from standard"
            )


class TestBatchModeEquivalence:
    """Tests for batch processing mode equivalence."""

    def test_batch_standard_vs_fast(
        self,
        batch_dataset: tuple[Path, Path],
        run_batch: Callable,
        e2e_temp_dir: Path,
    ):
        """Batch standard and fast modes should produce same total reads."""
        blast_dir, ani_path = batch_dataset

        # Run standard batch
        standard_result = run_batch(
            blast_dir=blast_dir,
            ani=ani_path,
            output_dir_name="batch_standard",
        )
        CLIAssertions.assert_success(standard_result)

        # Run fast batch
        fast_result = run_batch(
            blast_dir=blast_dir,
            ani=ani_path,
            output_dir_name="batch_fast",
            fast=True,
        )
        CLIAssertions.assert_success(fast_result)

        # Compare total reads
        standard_dir = e2e_temp_dir / "batch_standard"
        fast_dir = e2e_temp_dir / "batch_fast"

        standard_total = sum(
            len(pl.read_csv(f)) for f in standard_dir.glob("*_classifications.csv")
        )
        fast_total = sum(
            len(pl.read_csv(f)) for f in fast_dir.glob("*_classifications.csv")
        )

        assert standard_total == fast_total, (
            f"Total reads mismatch: standard={standard_total}, fast={fast_total}"
        )

    def test_batch_standard_vs_parallel(
        self,
        batch_dataset: tuple[Path, Path],
        run_batch: Callable,
        e2e_temp_dir: Path,
    ):
        """Batch standard and parallel modes should produce same total reads."""
        blast_dir, ani_path = batch_dataset

        # Run standard batch
        standard_result = run_batch(
            blast_dir=blast_dir,
            ani=ani_path,
            output_dir_name="batch_std_par",
        )
        CLIAssertions.assert_success(standard_result)

        # Run parallel batch
        parallel_result = run_batch(
            blast_dir=blast_dir,
            ani=ani_path,
            output_dir_name="batch_parallel",
            parallel=True,
        )
        CLIAssertions.assert_success(parallel_result)

        # Compare
        standard_dir = e2e_temp_dir / "batch_std_par"
        parallel_dir = e2e_temp_dir / "batch_parallel"

        standard_total = sum(
            len(pl.read_csv(f)) for f in standard_dir.glob("*_classifications.csv")
        )
        parallel_total = sum(
            len(pl.read_csv(f)) for f in parallel_dir.glob("*_classifications.csv")
        )

        assert standard_total == parallel_total
