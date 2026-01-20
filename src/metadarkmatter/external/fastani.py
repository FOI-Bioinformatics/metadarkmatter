"""
FastANI wrapper for computing Average Nucleotide Identity.

Provides Python interface for:
- FastANI: All-vs-all ANI computation between genome sets
"""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

from metadarkmatter.external.base import ExternalTool


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
