"""
Prodigal gene prediction wrapper.

Prodigal (PROkaryotic DYnamic programming Gene-finding ALgorithm) is a
gene prediction tool for prokaryotic genomes. It can predict protein-coding
genes in finished or draft genomes.

Reference:
    Hyatt D, et al. (2010). Prodigal: prokaryotic gene recognition and
    translation initiation site identification. BMC Bioinformatics, 11:119.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import ClassVar

from metadarkmatter.external.base import ExternalTool, ToolResult

logger = logging.getLogger(__name__)


class Prodigal(ExternalTool):
    """Wrapper for Prodigal gene prediction tool.

    Prodigal predicts protein-coding genes in prokaryotic genomes.
    It supports both single genome mode (for isolate genomes) and
    meta mode (for metagenomic contigs).

    Install via conda:
        conda install -c bioconda prodigal
    """

    TOOL_NAME: ClassVar[str] = "prodigal"
    TOOL_ALIASES: ClassVar[tuple[str, ...]] = ()
    INSTALL_HINT: ClassVar[str] = "conda install -c bioconda prodigal"

    def build_command(self, **kwargs: object) -> list[str]:
        """Build Prodigal command.

        Keyword Args:
            input_file: Input genome FASTA file
            output_proteins: Output protein FASTA file (.faa)
            output_genes: Output gene nucleotide FASTA file (.fna) [optional]
            output_gff: Output GFF3 annotation file [optional]
            procedure: Prediction procedure - 'single' or 'meta' (default: single)
            quiet: Suppress output (default: True)
            closed_ends: Treat edges as closed (default: False)

        Returns:
            Command as list of strings
        """
        input_file = kwargs.get("input_file")
        output_proteins = kwargs.get("output_proteins")
        output_genes = kwargs.get("output_genes")
        output_gff = kwargs.get("output_gff")
        procedure = kwargs.get("procedure", "single")
        quiet = kwargs.get("quiet", True)
        closed_ends = kwargs.get("closed_ends", False)

        if not input_file:
            msg = "input_file parameter is required"
            raise ValueError(msg)
        if not output_proteins:
            msg = "output_proteins parameter is required"
            raise ValueError(msg)

        cmd = [
            str(self.get_executable()),
            "-i", str(input_file),
            "-a", str(output_proteins),
        ]

        # Add optional outputs
        if output_genes:
            cmd.extend(["-d", str(output_genes)])
        if output_gff:
            cmd.extend(["-f", "gff", "-o", str(output_gff)])

        # Procedure mode
        if procedure == "meta":
            cmd.append("-p")
            cmd.append("meta")
        else:
            # Single genome mode (default)
            cmd.append("-p")
            cmd.append("single")

        # Quiet mode
        if quiet:
            cmd.append("-q")

        # Closed ends (don't allow genes to run off edges)
        if closed_ends:
            cmd.append("-c")

        return cmd

    def predict(
        self,
        input_file: Path,
        output_proteins: Path,
        *,
        output_genes: Path | None = None,
        output_gff: Path | None = None,
        procedure: str = "single",
        quiet: bool = True,
        closed_ends: bool = False,
        timeout: float | None = 600.0,
    ) -> ToolResult:
        """Run Prodigal gene prediction on a genome.

        Args:
            input_file: Input genome FASTA file (.fna)
            output_proteins: Output protein FASTA file (.faa)
            output_genes: Output gene nucleotide FASTA file (optional)
            output_gff: Output GFF3 annotation file (optional)
            procedure: 'single' for isolate genomes, 'meta' for metagenomes
            quiet: Suppress Prodigal output
            closed_ends: Don't allow genes to run off contig edges
            timeout: Timeout in seconds (default 10 minutes)

        Returns:
            ToolResult with command output
        """
        # Create output directory if needed
        output_proteins.parent.mkdir(parents=True, exist_ok=True)

        return self.run(
            input_file=input_file,
            output_proteins=output_proteins,
            output_genes=output_genes,
            output_gff=output_gff,
            procedure=procedure,
            quiet=quiet,
            closed_ends=closed_ends,
            timeout=timeout,
        )

    def predict_with_prefix(
        self,
        input_file: Path,
        output_proteins: Path,
        genome_accession: str,
        *,
        procedure: str = "single",
        quiet: bool = True,
        timeout: float | None = 600.0,
    ) -> ToolResult:
        """Run Prodigal and reformat headers with genome accession prefix.

        This produces protein files compatible with AAI computation,
        where headers are formatted as: >{accession}|{protein_id}

        Args:
            input_file: Input genome FASTA file (.fna)
            output_proteins: Output protein FASTA file (.faa)
            genome_accession: Genome accession to use as prefix
            procedure: 'single' for isolate genomes, 'meta' for metagenomes
            quiet: Suppress Prodigal output
            timeout: Timeout in seconds

        Returns:
            ToolResult with command output
        """
        import tempfile

        # Run Prodigal to temp file
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".faa", delete=False
        ) as tmp:
            tmp_path = Path(tmp.name)

        try:
            result = self.predict(
                input_file=input_file,
                output_proteins=tmp_path,
                procedure=procedure,
                quiet=quiet,
                timeout=timeout,
            )

            if not result.success:
                return result

            # Reformat headers with genome accession prefix
            self._reformat_headers(tmp_path, output_proteins, genome_accession)

            return result

        finally:
            # Clean up temp file
            if tmp_path.exists():
                tmp_path.unlink()

    def _reformat_headers(
        self,
        input_file: Path,
        output_file: Path,
        genome_accession: str,
    ) -> int:
        """Reformat Prodigal protein headers with genome accession prefix.

        Prodigal headers look like:
            >contig_1 # 1 # 1374 # 1 # ID=1_1;partial=10;...

        This reformats to:
            >GCF_000123456.1|contig_1_1

        Args:
            input_file: Input Prodigal protein FASTA
            output_file: Output reformatted FASTA
            genome_accession: Genome accession to use as prefix

        Returns:
            Number of proteins reformatted
        """
        # Pattern to extract contig and protein ID from Prodigal header
        # Example: >contig_1 # 1 # 1374 # 1 # ID=1_1;partial=10;...
        header_pattern = re.compile(r"^>(\S+)\s+#.*ID=(\d+_\d+);")

        protein_count = 0
        output_file.parent.mkdir(parents=True, exist_ok=True)

        with input_file.open("r") as fin, output_file.open("w") as fout:
            for line in fin:
                if line.startswith(">"):
                    match = header_pattern.match(line)
                    if match:
                        contig = match.group(1)
                        protein_id = match.group(2)
                        # Format: >accession|contig_proteinid
                        new_header = f">{genome_accession}|{contig}_{protein_id}\n"
                        fout.write(new_header)
                        protein_count += 1
                    else:
                        # Fallback: just prefix with accession
                        old_id = line[1:].split()[0]
                        fout.write(f">{genome_accession}|{old_id}\n")
                        protein_count += 1
                else:
                    fout.write(line)

        return protein_count
