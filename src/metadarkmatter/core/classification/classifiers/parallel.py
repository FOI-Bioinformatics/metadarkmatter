"""
Parallel classifier for large-scale BLAST file processing.

This module provides ParallelClassifier, which uses multiprocessing to
distribute classification workload across CPU cores for files with
millions of reads.
"""

from __future__ import annotations

import multiprocessing as mp
from collections.abc import Iterator
from functools import partial
from pathlib import Path
from typing import Any, TypeAlias

import numpy as np
import polars as pl

from metadarkmatter.core.classification.ani_matrix import ANIMatrix
from metadarkmatter.core.classification.classifiers.base import _classify_chunk_worker
from metadarkmatter.core.constants import calculate_confidence_score
from metadarkmatter.core.io_utils import write_dataframe
from metadarkmatter.core.parsers import StreamingBlastParser
from metadarkmatter.models.classification import TaxonomicCall, TAXONOMIC_TO_DIVERSITY
from metadarkmatter.models.config import ScoringConfig

# Type aliases for complex data structures used in parallel classification
# HitData = (qseqid, sseqid, pident, bitscore, genome_name)
HitData: TypeAlias = tuple[str, str, float, float, str]
# ReadChunk = (read_id, list of hits)
ReadChunk: TypeAlias = tuple[str, list[HitData]]
# ClassificationResult dict with standard keys
ClassificationDict: TypeAlias = dict[str, Any]

class ParallelClassifier:
    """
    Parallel classifier for very large BLAST files.

    Uses multiprocessing to distribute classification across CPU cores.
    Ideal for files with 10M+ reads where single-threaded processing
    becomes a bottleneck.

    The ANI matrix is shared across processes using shared memory,
    minimizing memory overhead.

    Example:
        ani_matrix = ANIMatrix.from_file(ani_path)
        parallel = ParallelClassifier(ani_matrix, config, num_workers=8)
        df = parallel.classify_file(blast_path)
    """

    # Bounds for chunk_size parameter
    MIN_CHUNK_SIZE = 100
    MAX_CHUNK_SIZE = 1_000_000

    def __init__(
        self,
        ani_matrix: ANIMatrix,
        config: ScoringConfig | None = None,
        num_workers: int | None = None,
        chunk_size: int = 50_000,
    ) -> None:
        """
        Initialize parallel classifier.

        Args:
            ani_matrix: Precomputed ANI matrix.
            config: Scoring configuration (uses defaults if None).
            num_workers: Number of worker processes (default: CPU count - 1).
            chunk_size: Reads per worker chunk (default: 50K).
                Valid range: 100 to 1,000,000.
                Memory usage scales with chunk_size * num_workers.
                50K reads/chunk * 8 workers ~= 400K reads in memory at once.
                Each read with 10 hits ~= 1KB, so ~400MB peak memory.
                Reduce for memory-constrained systems; increase for faster I/O.

        Raises:
            ValueError: If chunk_size is outside valid range.
        """
        if not self.MIN_CHUNK_SIZE <= chunk_size <= self.MAX_CHUNK_SIZE:
            msg = (
                f"chunk_size must be between {self.MIN_CHUNK_SIZE:,} and "
                f"{self.MAX_CHUNK_SIZE:,}, got {chunk_size:,}"
            )
            raise ValueError(msg)

        self.ani_matrix = ani_matrix
        self.config = config or ScoringConfig()
        self.num_workers = num_workers or max(1, mp.cpu_count() - 1)
        self.chunk_size = chunk_size

        # Get effective thresholds based on alignment mode (protein vs nucleotide)
        eff = self.config.get_effective_thresholds()

        # Pre-serialize config for workers using effective thresholds
        # Includes confidence score scaling parameters for mode-specific scoring
        self._config_dict = {
            # Classification thresholds
            "bitscore_threshold_pct": self.config.bitscore_threshold_pct,
            "genus_bitscore_threshold_pct": self.config.genus_bitscore_threshold_pct,
            "novelty_known_max": eff["novelty_known_max"],
            "novelty_novel_species_min": eff["novelty_novel_species_min"],
            "novelty_novel_species_max": eff["novelty_novel_species_max"],
            "novelty_novel_genus_min": eff["novelty_novel_genus_min"],
            "novelty_novel_genus_max": eff["novelty_novel_genus_max"],
            "uncertainty_known_max": eff["uncertainty_known_max"],
            "uncertainty_novel_species_max": eff["uncertainty_novel_species_max"],
            "uncertainty_novel_genus_max": eff["uncertainty_novel_genus_max"],
            "uncertainty_conserved_min": eff["uncertainty_conserved_min"],
            "genus_uncertainty_ambiguous_min": self.config.genus_uncertainty_ambiguous_min,
            "identity_gap_ambiguous_max": self.config.identity_gap_ambiguous_max,
            "default_ani": self.ani_matrix._default_ani,
            # Confidence score scaling parameters (mode-specific)
            "margin_divisor_known": eff["margin_divisor_known"],
            "margin_divisor_novel_species": eff["margin_divisor_novel_species"],
            "margin_divisor_novel_genus": eff["margin_divisor_novel_genus"],
            "identity_gap_thresholds": eff["identity_gap_thresholds"],
            "identity_score_base": eff["identity_score_base"],
            "identity_score_range": eff["identity_score_range"],
            # Coverage weighting parameters
            "coverage_weight_mode": self.config.coverage_weight_mode,
            "coverage_weight_strength": self.config.coverage_weight_strength,
            # Uncertainty calculation mode
            "uncertainty_mode": self.config.uncertainty_mode,
            # Enhanced scoring options
            "enhanced_scoring": self.config.enhanced_scoring,
            "infer_single_hit_uncertainty": self.config.infer_single_hit_uncertainty,
        }

    def _iter_chunks(
        self,
        blast_path: Path,
    ) -> Iterator[list[tuple[str, list[tuple[str, str, float, float, str]]]]]:
        """Yield chunks of reads lazily to avoid memory accumulation.

        This generator yields chunks on demand rather than materializing
        all chunks upfront, enabling memory-efficient parallel processing.

        Args:
            blast_path: Path to BLAST tabular output

        Yields:
            Lists of (read_id, hits_data) tuples, one chunk at a time
        """
        parser = StreamingBlastParser(blast_path)
        current_chunk: list[tuple[str, list[tuple[str, str, float, float, str]]]] = []

        for result in parser.iter_reads_fast():
            # Convert to serializable format
            hits_data = [
                (h.qseqid, h.sseqid, h.pident, h.bitscore, h.genome_name)
                for h in result.hits
            ]
            current_chunk.append((result.read_id, hits_data))

            if len(current_chunk) >= self.chunk_size:
                yield current_chunk
                current_chunk = []

        if current_chunk:
            yield current_chunk

    def classify_file(
        self,
        blast_path: Path,
        progress_callback: Callable[[int], None] | None = None,
    ) -> pl.DataFrame:
        """
        Classify BLAST file using parallel processing.

        Distributes work across multiple CPU cores for faster processing.
        Uses lazy chunk iteration to minimize memory usage.
        Results are collected and returned as a single DataFrame.

        Args:
            blast_path: Path to BLAST tabular output
            progress_callback: Optional callback(chunks_completed) for progress

        Returns:
            DataFrame with classification results
        """
        # Process chunks in parallel using lazy iteration
        all_results: list[dict] = []

        # Get shared data for workers
        ani_array = self.ani_matrix._ani_array
        genome_to_idx = self.ani_matrix._genome_to_idx

        # Create worker function with fixed arguments
        worker_fn = partial(
            _classify_chunk_worker,
            ani_array=ani_array,
            genome_to_idx=genome_to_idx,
            config_dict=self._config_dict,
        )

        # Use lazy chunk iteration to avoid memory accumulation
        # pool.imap() processes chunks as they're yielded by the generator
        chunk_iter = self._iter_chunks(blast_path)

        with mp.Pool(processes=self.num_workers) as pool:
            for i, chunk_results in enumerate(pool.imap(worker_fn, chunk_iter)):
                all_results.extend(chunk_results)
                if progress_callback:
                    progress_callback(i + 1)

        if not all_results:
            return self._empty_dataframe()

        return pl.DataFrame(all_results)

    def stream_to_file(
        self,
        blast_path: Path,
        output_path: Path,
        output_format: str = "parquet",
        progress_callback: Callable[[int], None] | None = None,
    ) -> int:
        """
        Classify and stream results to file using parallel processing.

        Combines parallel classification with chunked output for
        maximum scalability on very large files.

        Args:
            blast_path: Path to BLAST tabular output
            output_path: Output file path
            output_format: 'parquet' or 'csv'
            progress_callback: Optional callback(chunks_completed)

        Returns:
            Total reads classified
        """
        df = self.classify_file(blast_path, progress_callback)
        write_dataframe(df, output_path, output_format)
        return len(df)

    def _empty_dataframe(self) -> pl.DataFrame:
        """Return empty DataFrame with correct schema."""
        return pl.DataFrame(
            schema={
                "read_id": pl.Utf8,
                "best_match_genome": pl.Utf8,
                "top_hit_identity": pl.Float64,
                "novelty_index": pl.Float64,
                "placement_uncertainty": pl.Float64,
                "num_ambiguous_hits": pl.Int64,
                "second_hit_identity": pl.Float64,
                "identity_gap": pl.Float64,
                "confidence_score": pl.Float64,
                "taxonomic_call": pl.Utf8,
                "diversity_status": pl.Utf8,
                "is_novel": pl.Boolean,
            }
        )


