"""
FastANI wrapper for computing Average Nucleotide Identity.

Provides Python interface for:
- FastANI: All-vs-all ANI computation between genome sets
"""

from __future__ import annotations

import logging
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import ClassVar

from metadarkmatter.external.base import ExternalTool

logger = logging.getLogger(__name__)


class FastANI(ExternalTool):
    """Wrapper for fastANI genome comparison tool.

    FastANI computes Average Nucleotide Identity (ANI) between bacterial
    genomes using a fast k-mer based algorithm. Outputs include pairwise
    ANI values and an optional lower triangular PHYLIP matrix.

    Note:
        FastANI does not report ANI values below approximately 80%.
        Such pairs will be absent from the output.

    Example:
        >>> fastani = FastANI()
        >>> result = fastani.run(
        ...     query_list=Path("genomes_list.txt"),
        ...     reference_list=Path("genomes_list.txt"),
        ...     output=Path("ani_results.txt"),
        ...     threads=16,
        ...     matrix=True,
        ... )
    """

    TOOL_NAME: ClassVar[str] = "fastANI"
    TOOL_ALIASES: ClassVar[tuple[str, ...]] = ("fastani",)
    INSTALL_HINT: ClassVar[str] = "conda install -c bioconda fastani"

    def build_command(
        self,
        *,
        query_list: Path,
        reference_list: Path,
        output: Path,
        threads: int = 4,
        kmer_size: int = 16,
        frag_len: int = 3000,
        min_fraction: float = 0.2,
        matrix: bool = True,
    ) -> list[str]:
        """Build fastANI all-vs-all command.

        Args:
            query_list: File containing paths to query genomes (one per line).
            reference_list: File containing paths to reference genomes.
                For all-vs-all, use the same file as query_list.
            output: Output file path for ANI results.
            threads: Number of threads for parallel execution.
            kmer_size: K-mer size for comparison (default: 16).
            frag_len: Fragment length (default: 3000).
            min_fraction: Minimum fraction of genome fragments required
                for ANI calculation.
            matrix: If True, output PHYLIP-format distance matrix as
                additional file with .matrix extension.

        Returns:
            Command as list of strings.
        """
        exe = str(self.get_executable())
        cmd = [exe]

        # Query and reference lists
        cmd.extend(["--ql", str(query_list)])
        cmd.extend(["--rl", str(reference_list)])

        # Output
        cmd.extend(["-o", str(output)])

        # Threading
        cmd.extend(["-t", str(threads)])

        # Algorithm parameters
        cmd.extend(["-k", str(kmer_size)])
        cmd.extend(["--fragLen", str(frag_len)])
        cmd.extend(["--minFraction", str(min_fraction)])

        # Matrix output
        if matrix:
            cmd.append("--matrix")

        return cmd

    def compute_parallel_batches(
        self,
        *,
        query_list: Path,
        reference_list: Path,
        output: Path,
        tmp_dir: Path,
        threads: int = 4,
        batches: int = 1,
        kmer_size: int = 16,
        frag_len: int = 3000,
        min_fraction: float = 0.2,
    ) -> None:
        """Run fastANI across multiple parallel process batches.

        The query list is partitioned into ``batches`` equal-sized
        chunks. Each batch runs fastANI in its own subprocess with
        ``threads // batches`` internal threads, and the per-batch
        outputs are concatenated into ``output``. The reference list
        is shared across all batches so the union of outputs is a
        complete all-vs-all comparison (modulo fastANI's own >=80%
        ANI cutoff).

        With ``batches=1`` this is equivalent to a single fastANI
        invocation. For large genome sets (>100s) batching often
        beats fastANI's internal threading because each batch process
        owns its own memory and avoids contention on the shared
        global state.

        Args:
            query_list: File of query genome paths, one per line.
            reference_list: File of reference genome paths.
            output: Final merged output path.
            tmp_dir: Directory for per-batch intermediate files;
                must exist and be writable.
            threads: Total CPU budget shared across batches.
            batches: Number of parallel batches (default 1).
            kmer_size: fastANI k-mer size.
            frag_len: fastANI fragment length.
            min_fraction: fastANI minimum aligned fraction.

        Raises:
            ToolExecutionError: If any batch fails.
        """
        if batches < 1:
            msg = "batches must be >= 1"
            raise ValueError(msg)

        queries = [
            line.strip()
            for line in query_list.read_text().splitlines()
            if line.strip()
        ]
        if not queries:
            msg = f"Query list is empty: {query_list}"
            raise ValueError(msg)

        # Cap batches at the query count so we never make empty sub-lists.
        effective_batches = min(batches, len(queries))
        threads_per_batch = max(1, threads // effective_batches)

        # Partition queries into roughly equal contiguous chunks. Keep the
        # original order so the merged output remains deterministic.
        per_batch = -(-len(queries) // effective_batches)  # ceil division
        chunks: list[list[str]] = [
            queries[i : i + per_batch] for i in range(0, len(queries), per_batch)
        ]

        # Single-batch fast path: avoid the subprocess pool overhead.
        if effective_batches == 1:
            self.run_or_raise(
                query_list=query_list,
                reference_list=reference_list,
                output=output,
                threads=threads,
                kmer_size=kmer_size,
                frag_len=frag_len,
                min_fraction=min_fraction,
                matrix=False,
            )
            return

        # Write per-batch query list files.
        batch_inputs: list[tuple[Path, Path]] = []
        for i, chunk in enumerate(chunks):
            chunk_list = tmp_dir / f"query_batch_{i}.txt"
            chunk_list.write_text("\n".join(chunk) + "\n")
            chunk_out = tmp_dir / f"fastani_batch_{i}.tsv"
            batch_inputs.append((chunk_list, chunk_out))

        exe_path = str(self.get_executable())

        # Run each batch as its own fastANI subprocess. We use a
        # ThreadPoolExecutor because subprocess.run releases the GIL for
        # the duration of the external tool call, so threads give the
        # same throughput as processes without the pickling overhead.
        logger.info(
            "fastANI: %d query genomes across %d parallel batches "
            "(%d threads per batch).",
            len(queries), effective_batches, threads_per_batch,
        )
        futures = {}
        with ThreadPoolExecutor(max_workers=effective_batches) as pool:
            for chunk_list, chunk_out in batch_inputs:
                fut = pool.submit(
                    _run_fastani_batch,
                    exe_path,
                    chunk_list,
                    reference_list,
                    chunk_out,
                    threads_per_batch,
                    kmer_size,
                    frag_len,
                    min_fraction,
                )
                futures[fut] = chunk_out

            for fut in as_completed(futures):
                fut.result()  # propagate exceptions

        # Concatenate per-batch outputs in the order we submitted them so
        # the result is byte-identical across runs.
        with output.open("w") as merged:
            for _, chunk_out in batch_inputs:
                if chunk_out.exists():
                    merged.write(chunk_out.read_text())


def _run_fastani_batch(
    exe_path: str,
    query_list: Path,
    reference_list: Path,
    output: Path,
    threads: int,
    kmer_size: int,
    frag_len: int,
    min_fraction: float,
) -> None:
    """Module-level worker: run a single fastANI invocation.

    Raises a plain RuntimeError on failure; the caller decides how to
    surface it.
    """
    cmd = [
        exe_path,
        "--ql", str(query_list),
        "--rl", str(reference_list),
        "-o", str(output),
        "-t", str(threads),
        "-k", str(kmer_size),
        "--fragLen", str(frag_len),
        "--minFraction", str(min_fraction),
    ]
    result = subprocess.run(
        cmd, capture_output=True, text=True, check=False
    )
    if result.returncode != 0:
        msg = (
            f"fastANI batch failed (exit {result.returncode}): "
            f"{result.stderr.strip()[:500]}"
        )
        raise RuntimeError(msg)


def create_genome_list_file(
    genome_dir: Path,
    output_path: Path,
    patterns: tuple[str, ...] = ("*.fna", "*.fa", "*.fasta", "*.fna.gz", "*.fa.gz"),
) -> list[Path]:
    """Create a genome list file for fastANI from a directory of genomes.

    Args:
        genome_dir: Directory containing genome FASTA files.
        output_path: Path to write the genome list file.
        patterns: Glob patterns for genome files.

    Returns:
        List of genome paths found.

    Raises:
        FileNotFoundError: If no genomes match the patterns.
    """
    genome_paths: list[Path] = []

    for pattern in patterns:
        genome_paths.extend(genome_dir.glob(pattern))

    # Remove duplicates and sort
    genome_paths = sorted(set(genome_paths))

    if not genome_paths:
        msg = f"No genome files found in {genome_dir} matching patterns: {patterns}"
        raise FileNotFoundError(msg)

    # Write list file with absolute paths
    with output_path.open("w") as f:
        for path in genome_paths:
            f.write(f"{path.resolve()}\n")

    return genome_paths
