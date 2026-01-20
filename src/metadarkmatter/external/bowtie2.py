"""
Bowtie2 wrapper classes.

Provides Python interfaces for:
- Bowtie2Build: Building Bowtie2 indices from reference genomes
- Bowtie2: Aligning reads to indexed reference genomes
"""

from __future__ import annotations

from pathlib import Path

from metadarkmatter.external.base import ExternalTool


class Bowtie2Build(ExternalTool):
    """Wrapper for bowtie2-build index construction.

    Builds a Bowtie2 index from one or more reference sequences.
    The index is required before read alignment with Bowtie2.

    Example:
        >>> builder = Bowtie2Build()
        >>> result = builder.run(
        ...     reference=Path("pangenome.fasta"),
        ...     index_prefix=Path("pangenome_db"),
        ...     threads=8,
        ... )
    """

    TOOL_NAME = "bowtie2-build"
    INSTALL_HINT = "conda install -c bioconda bowtie2"

    def build_command(
        self,
        *,
        reference: Path,
        index_prefix: Path,
        threads: int = 4,
        large_index: bool = False,
        seed: int | None = None,
        quiet: bool = False,
    ) -> list[str]:
        """Build bowtie2-build command.

        Args:
            reference: FASTA file with reference sequences.
            index_prefix: Prefix for output index files.
            threads: Number of threads for parallel construction.
            large_index: Force large index (for references > 4 billion bases).
            seed: Random seed for reproducibility.
            quiet: Suppress progress messages.

        Returns:
            Command as list of strings.
        """
        exe = str(self.get_executable())
        cmd = [exe]

        # Threading
        cmd.extend(["--threads", str(threads)])

        # Options
        if large_index:
            cmd.append("--large-index")

        if seed is not None:
            cmd.extend(["--seed", str(seed)])

        if quiet:
            cmd.append("--quiet")

        # Required positional arguments
        cmd.append(str(reference))
        cmd.append(str(index_prefix))

        return cmd


class Bowtie2(ExternalTool):
    """Wrapper for Bowtie2 read aligner.

    Aligns short reads to a reference genome using the Bowtie2 algorithm.
    Supports both single-end and paired-end reads.

    Example:
        >>> aligner = Bowtie2()
        >>> result = aligner.run(
        ...     index_prefix=Path("pangenome_db"),
        ...     reads_1=Path("sample_R1.fastq.gz"),
        ...     reads_2=Path("sample_R2.fastq.gz"),
        ...     output_sam=Path("aligned.sam"),
        ...     mode="local",
        ...     max_alignments=100,
        ... )
    """

    TOOL_NAME = "bowtie2"
    INSTALL_HINT = "conda install -c bioconda bowtie2"

    def build_command(
        self,
        *,
        index_prefix: Path,
        reads_1: Path,
        output_sam: Path,
        reads_2: Path | None = None,
        threads: int = 4,
        mode: str = "local",
        max_alignments: int = 1,
        no_unal: bool = True,
        no_discordant: bool = False,
        no_mixed: bool = False,
        very_sensitive: bool = False,
        seed: int | None = None,
        min_insert: int | None = None,
        max_insert: int | None = None,
        interleaved: bool = False,
        quiet: bool = False,
    ) -> list[str]:
        """Build bowtie2 command.

        Args:
            index_prefix: Prefix of Bowtie2 index files.
            reads_1: Forward reads file (or interleaved).
            output_sam: Output SAM file path.
            reads_2: Reverse reads file for paired-end.
            threads: Number of threads.
            mode: Alignment mode - "local" or "end-to-end".
            max_alignments: Report up to this many alignments per read (-k).
            no_unal: Do not include unaligned reads in output.
            no_discordant: Do not report discordant paired-end alignments.
            no_mixed: Do not report unpaired alignments for paired reads.
            very_sensitive: Use very-sensitive preset (slower but more accurate).
            seed: Random seed for reproducibility.
            min_insert: Minimum insert size for paired-end.
            max_insert: Maximum insert size for paired-end.
            interleaved: Reads are interleaved in single file.
            quiet: Suppress progress messages.

        Returns:
            Command as list of strings.
        """
        exe = str(self.get_executable())
        cmd = [exe]

        # Index
        cmd.extend(["-x", str(index_prefix)])

        # Input reads
        if reads_2 is not None:
            cmd.extend(["-1", str(reads_1)])
            cmd.extend(["-2", str(reads_2)])
        elif interleaved:
            cmd.append("--interleaved")
            cmd.extend(["-U", str(reads_1)])
        else:
            cmd.extend(["-U", str(reads_1)])

        # Output
        cmd.extend(["-S", str(output_sam)])

        # Threading
        cmd.extend(["--threads", str(threads)])

        # Alignment mode
        if mode == "local":
            cmd.append("--local")
        elif mode == "end-to-end":
            cmd.append("--end-to-end")

        # Sensitivity preset
        if very_sensitive:
            if mode == "local":
                cmd.append("--very-sensitive-local")
            else:
                cmd.append("--very-sensitive")

        # Multiple alignments
        if max_alignments > 1:
            cmd.extend(["-k", str(max_alignments)])

        # Filtering options
        if no_unal:
            cmd.append("--no-unal")

        if no_discordant:
            cmd.append("--no-discordant")

        if no_mixed:
            cmd.append("--no-mixed")

        # Insert size constraints
        if min_insert is not None:
            cmd.extend(["--minins", str(min_insert)])

        if max_insert is not None:
            cmd.extend(["--maxins", str(max_insert)])

        # Misc options
        if seed is not None:
            cmd.extend(["--seed", str(seed)])

        if quiet:
            cmd.append("--quiet")

        return cmd


def concatenate_genomes(
    genome_dir: Path,
    output_path: Path,
    contig_mapping_path: Path | None = None,
    pattern: str = "*.fna",
) -> list[str]:
    """Concatenate multiple genome FASTA files into a pangenome.

    Creates a single FASTA file from all genomes in a directory with
    standardized headers in format: {accession}|{original_contig_id}

    This ensures reliable genome identification from mapping hits.

    Args:
        genome_dir: Directory containing genome FASTA files.
        output_path: Output path for concatenated pangenome.
        contig_mapping_path: Output path for contig mapping TSV (optional).
        pattern: Glob pattern for genome files.

    Returns:
        List of genome accessions included in the pangenome.

    Example:
        >>> genomes = concatenate_genomes(
        ...     genome_dir=Path("./references/"),
        ...     output_path=Path("./pangenome.fasta"),
        ...     contig_mapping_path=Path("./contig_mapping.tsv"),
        ... )
        >>> print(f"Concatenated {len(genomes)} genomes")
    """
    from metadarkmatter.core.genome_utils import (
        concatenate_genomes_with_mapping,
        extract_accession_from_filename,
    )

    # If no contig mapping path provided, use a default location
    if contig_mapping_path is None:
        contig_mapping_path = output_path.with_suffix(".contig_mapping.tsv")

    genome_count, _ = concatenate_genomes_with_mapping(
        genome_dir=genome_dir,
        output_fasta=output_path,
        contig_mapping_path=contig_mapping_path,
        pattern=pattern,
    )

    if genome_count == 0:
        raise FileNotFoundError(
            f"No genome files found matching '{pattern}' in {genome_dir}"
        )

    # Return list of genome accessions for backwards compatibility
    genome_files = sorted(genome_dir.glob(pattern))
    return [extract_accession_from_filename(f.name) for f in genome_files]
