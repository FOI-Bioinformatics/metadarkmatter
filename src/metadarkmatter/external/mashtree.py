"""
Mashtree wrapper for building phylogenetic trees from genome assemblies.

Mashtree uses Mash (MinHash) distances between genomes to construct a
neighbor-joining tree. It is fast and suitable for quick phylogenetic
overviews of genome collections.
"""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

from metadarkmatter.external.base import ExternalTool


class Mashtree(ExternalTool):
    """Wrapper for mashtree phylogenetic tree builder.

    Mashtree builds a tree from genome assemblies using Mash distances.
    It takes FASTA files as input and writes a Newick tree to stdout.

    Example:
        >>> mt = Mashtree()
        >>> result = mt.run(genomes=[Path("g1.fna"), Path("g2.fna")], threads=4)
        >>> newick = result.stdout.strip()
    """

    TOOL_NAME: ClassVar[str] = "mashtree"
    INSTALL_HINT: ClassVar[str] = "conda install -c bioconda mashtree"

    def build_command(
        self,
        *,
        genomes: list[Path],
        threads: int = 4,
    ) -> list[str]:
        """Build mashtree command.

        Mashtree reads genome FASTA files and writes a Newick tree to stdout.

        Args:
            genomes: List of paths to genome FASTA files.
            threads: Number of CPU threads for parallel execution.

        Returns:
            Command as list of strings.
        """
        exe = str(self.get_executable())
        cmd = [exe, "--numcpus", str(threads)]
        cmd.extend(str(g) for g in genomes)
        return cmd
