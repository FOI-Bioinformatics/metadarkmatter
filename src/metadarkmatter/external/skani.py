"""
Skani wrapper for computing Average Nucleotide Identity.

Provides Python interface for:
- Skani: Fast ANI computation with aligned fraction metrics
"""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

from metadarkmatter.external.base import ExternalTool


class Skani(ExternalTool):
    """Wrapper for skani genome comparison tool.

    Skani uses sparse chaining for fast ANI estimation, providing
    ANI values along with aligned fraction (AF) metrics. The triangle
    subcommand performs efficient all-vs-all comparison.

    Skani is generally faster and more memory-efficient than fastANI,
    especially for large genome sets.

    Example:
        >>> skani = Skani()
        >>> result = skani.run(
        ...     genomes=[Path("genome1.fa"), Path("genome2.fa")],
        ...     output=Path("ani_results.txt"),
        ...     threads=16,
        ...     full_matrix=True,
        ... )
    """

    TOOL_NAME: ClassVar[str] = "skani"
    INSTALL_HINT: ClassVar[str] = "conda install -c bioconda skani"

    MAX_ARGV_GENOMES: ClassVar[int] = 1000
    GENOME_SUFFIXES: ClassVar[tuple[str, ...]] = (
        ".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz",
    )

    @classmethod
    def list_genome_files(cls, directory: Path) -> list[Path]:
        """Sorted FASTA files directly inside ``directory``."""
        return sorted(
            p for p in directory.iterdir()
            if p.is_file() and p.name.lower().endswith(cls.GENOME_SUFFIXES)
        )

    def build_command(
        self,
        *,
        genomes: list[Path] | Path,
        output: Path,
        threads: int = 4,
        compression_factor: int | None = None,
        marker_compression_factor: int | None = None,
        min_af: float | None = None,
        screen_threshold: float | None = None,
        full_matrix: bool = True,
        sparse: bool = False,
        no_learned_ani: bool = False,
        mode: str | None = None,
    ) -> list[str]:
        """Build skani triangle (all-vs-all) command.

        Args:
            genomes: Directory containing genome files, or list of genome paths.
            output: Output file path.
            threads: Number of threads for parallel execution.
            compression_factor: Compression factor for sketching.
                Lower values increase sensitivity but use more memory.
            marker_compression_factor: Marker compression factor.
            min_af: Minimum aligned fraction to report (0-100).
                Default is 15. Lower values (e.g., 5) capture more distant
                comparisons but with lower confidence.
            screen_threshold: Pre-filter pairs below this ANI threshold (0-100).
                Default is 80. Lower values (e.g., 50) attempt to compute ANI
                for more distant genome pairs.
            full_matrix: If True, output full square matrix instead of
                lower triangular format.
            sparse: If True, output sparse edge-list format instead of matrix.
                Useful for very large datasets.
            no_learned_ani: If True, disable the learned ANI regression model.
                By default, skani uses learned ANI when c >= 70 and >= 150,000
                bases are aligned.
            mode: skani sensitivity mode: None/"default", "medium" (-c 70)
                or "slow" (-c 30). Mutually exclusive with compression_factor
                in skani, so ``compression_factor`` wins when both are given.

        Returns:
            Command as list of strings.
        """
        exe = str(self.get_executable())
        cmd = [exe, "triangle"]

        # Input genomes
        if isinstance(genomes, Path):
            # Directory: expand to the sorted FASTA file list (skani does not
            # glob; a literal "dir/*" would be passed as one filename).
            files = self.list_genome_files(genomes)
            if len(files) > self.MAX_ARGV_GENOMES:
                list_path = output.with_name(output.name + ".genomes.list")
                list_path.parent.mkdir(parents=True, exist_ok=True)
                list_path.write_text("".join(f"{p}\n" for p in files))
                cmd.extend(["-l", str(list_path)])
            else:
                cmd.extend(str(p) for p in files)
        else:
            # List of paths
            cmd.extend(str(p) for p in genomes)

        # Output
        cmd.extend(["-o", str(output)])

        # Threading
        cmd.extend(["-t", str(threads)])

        # Algorithm parameters
        if compression_factor is not None:
            cmd.extend(["-c", str(compression_factor)])

        if mode in ("medium", "slow") and compression_factor is None:
            cmd.append(f"--{mode}")
        elif mode not in (None, "default", "medium", "slow"):
            msg = f"Unknown skani mode: {mode!r} (expected default, medium or slow)"
            raise ValueError(msg)

        if marker_compression_factor is not None:
            cmd.extend(["-m", str(marker_compression_factor)])

        if min_af is not None:
            cmd.extend(["--min-af", str(min_af)])

        if screen_threshold is not None:
            cmd.extend(["-s", str(screen_threshold)])

        if no_learned_ani:
            cmd.append("--no-learned-ani")

        # Output format
        if full_matrix:
            cmd.append("--full-matrix")

        if sparse:
            cmd.append("--sparse")

        return cmd

    def build_dist_command(
        self,
        *,
        query: Path | list[Path],
        reference: Path | list[Path],
        output: Path,
        threads: int = 4,
    ) -> list[str]:
        """Build skani dist (query vs reference) command.

        Use this for comparing query genomes against a reference set,
        rather than all-vs-all comparison.

        Args:
            query: Query genome(s) or directory.
            reference: Reference genome(s) or directory.
            output: Output file path.
            threads: Number of threads.

        Returns:
            Command as list of strings.
        """
        exe = str(self.get_executable())
        cmd = [exe, "dist"]

        # Query
        if isinstance(query, Path):
            cmd.extend(["-q", str(query)])
        else:
            for q in query:
                cmd.extend(["-q", str(q)])

        # Reference
        if isinstance(reference, Path):
            cmd.extend(["-r", str(reference)])
        else:
            for r in reference:
                cmd.extend(["-r", str(r)])

        # Output and threading
        cmd.extend(["-o", str(output)])
        cmd.extend(["-t", str(threads)])

        return cmd
