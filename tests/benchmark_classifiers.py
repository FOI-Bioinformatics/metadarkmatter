#!/usr/bin/env python
"""
Benchmark suite for metadarkmatter classifiers.

Compares performance of different classification modes:
- Vectorized (VectorizedClassifier.classify_file)
- Streaming (VectorizedClassifier.stream_to_file)

Usage:
    python tests/benchmark_classifiers.py [--reads N] [--genomes N]

Example:
    python tests/benchmark_classifiers.py --reads 100000 --genomes 50
"""

from __future__ import annotations

import argparse
import gc
import tempfile
import time
from pathlib import Path

import numpy as np
import polars as pl


def generate_test_data(
    num_reads: int,
    num_genomes: int,
    output_dir: Path,
) -> tuple[Path, Path]:
    """Generate synthetic BLAST output and ANI matrix for benchmarking."""
    np.random.seed(42)

    # Generate genome names
    genomes = [f"GCF_{i:09d}.1" for i in range(1, num_genomes + 1)]

    # Generate BLAST output
    # Each read has 1-5 hits with realistic identity distribution
    blast_data = []
    for i in range(num_reads):
        qseqid = f"read_{i:08d}"
        n_hits = np.random.randint(1, 6)
        for _ in range(n_hits):
            sseqid = np.random.choice(genomes)
            # Identity distribution biased toward higher values (metagenomic pattern)
            pident = np.clip(np.random.beta(5, 1) * 30 + 70, 70, 100)
            length = np.random.randint(75, 250)
            bitscore = pident * length / 100 * 1.5  # Rough BLASTN approximation
            blast_data.append({
                "qseqid": qseqid,
                "sseqid": sseqid,
                "pident": round(pident, 2),
                "length": length,
                "mismatch": int((100 - pident) * length / 100),
                "gapopen": np.random.randint(0, 3),
                "qstart": 1,
                "qend": length,
                "sstart": np.random.randint(1, 100000),
                "send": np.random.randint(1, 100000),
                "evalue": float(f"1e-{np.random.randint(5, 50)}"),
                "bitscore": round(bitscore, 1),
            })

    blast_df = pl.DataFrame(blast_data)

    # Generate ANI matrix (symmetric, realistic values)
    ani_data = {"genome": genomes}
    for i, g in enumerate(genomes):
        ani_values = []
        for j, g2 in enumerate(genomes):
            if i == j:
                ani_values.append(100.0)
            elif abs(i - j) == 1:
                # Adjacent genomes: high ANI (same species cluster)
                ani_values.append(round(np.random.uniform(95, 99), 2))
            elif abs(i - j) <= 3:
                # Close genomes: moderate ANI (same genus)
                ani_values.append(round(np.random.uniform(85, 95), 2))
            else:
                # Distant genomes: low ANI
                ani_values.append(round(np.random.uniform(75, 85), 2))
        ani_data[g] = ani_values
    ani_df = pl.DataFrame(ani_data)

    # Write files
    blast_path = output_dir / "benchmark.blast.tsv"
    ani_path = output_dir / "benchmark.ani.csv"

    blast_df.write_csv(blast_path, separator="\t", include_header=False)
    ani_df.write_csv(ani_path)

    return blast_path, ani_path, len(blast_df)


def benchmark_vectorized(ani_matrix, config, blast_path: Path, output_path: Path) -> dict:
    """Benchmark vectorized classifier."""
    from metadarkmatter.core.ani_placement import VectorizedClassifier

    gc.collect()
    vectorized = VectorizedClassifier(ani_matrix=ani_matrix, config=config)

    start = time.perf_counter()
    df = vectorized.classify_file(blast_path)
    elapsed = time.perf_counter() - start

    df.write_parquet(output_path, compression="zstd")

    return {
        "mode": "vectorized (--parallel)",
        "reads": len(df),
        "elapsed": elapsed,
        "rate": len(df) / elapsed,
    }


def benchmark_streaming(ani_matrix, config, blast_path: Path, output_path: Path) -> dict:
    """Benchmark streaming classifier."""
    from metadarkmatter.core.ani_placement import VectorizedClassifier

    gc.collect()
    vectorized = VectorizedClassifier(ani_matrix=ani_matrix, config=config)

    start = time.perf_counter()
    num_reads = vectorized.stream_to_file(
        blast_path,
        output_path,
        output_format="parquet",
        partition_size=500_000,
    )
    elapsed = time.perf_counter() - start

    return {
        "mode": "streaming (--streaming)",
        "reads": num_reads,
        "elapsed": elapsed,
        "rate": num_reads / elapsed,
    }


def main():
    parser = argparse.ArgumentParser(description="Benchmark metadarkmatter classifiers")
    parser.add_argument("--reads", type=int, default=100_000, help="Number of reads")
    parser.add_argument("--genomes", type=int, default=50, help="Number of genomes")
    args = parser.parse_args()

    print("=" * 70)
    print("METADARKMATTER CLASSIFIER BENCHMARK")
    print("=" * 70)
    print(f"\nConfiguration:")
    print(f"  Reads: {args.reads:,}")
    print(f"  Genomes: {args.genomes}")
    print()

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Generate test data
        print("Generating test data...")
        blast_path, ani_path, num_alignments = generate_test_data(
            args.reads, args.genomes, tmpdir
        )
        print(f"  BLAST file: {blast_path.stat().st_size / 1024 / 1024:.1f} MB")
        print(f"  Total alignments: {num_alignments:,}")
        print(f"  ANI matrix: {args.genomes}x{args.genomes}")
        print()

        # Load ANI matrix once
        from metadarkmatter.core.ani_placement import ANIMatrix
        from metadarkmatter.models.config import ScoringConfig

        ani_matrix = ANIMatrix.from_file(ani_path)
        config = ScoringConfig()

        results = []

        # Run benchmarks
        print("Running benchmarks...")
        print("-" * 70)

        output = tmpdir / "vectorized.parquet"
        result = benchmark_vectorized(ani_matrix, config, blast_path, output)
        results.append(result)
        print(f"  Vectorized: {result['elapsed']:6.2f}s  ({result['rate']:10,.0f} reads/s)")

        output = tmpdir / "streaming.parquet"
        result = benchmark_streaming(ani_matrix, config, blast_path, output)
        results.append(result)
        print(f"  Streaming:  {result['elapsed']:6.2f}s  ({result['rate']:10,.0f} reads/s)")

        print("-" * 70)
        print()

        # Summary
        print("SUMMARY")
        print("-" * 70)
        baseline = results[0]
        print(f"{'Mode':<30} {'Time (s)':>10} {'Rate':>15} {'Speedup':>10}")
        print("-" * 70)
        for r in results:
            speedup = baseline["elapsed"] / r["elapsed"]
            print(
                f"{r['mode']:<30} {r['elapsed']:10.2f} {r['rate']:15,.0f} {speedup:9.1f}x"
            )
        print()

        # Recommendations
        print("RECOMMENDATIONS")
        print("-" * 70)
        print("  Small files (<1M reads):     Use --parallel (fastest)")
        print("  Medium files (1-50M reads):  Use --parallel (good balance)")
        print("  Large files (50-100M reads): Use --parallel (may need 32GB+ RAM)")
        print("  Huge files (100M+ reads):    Use --streaming (bounded memory)")
        print()


if __name__ == "__main__":
    main()
