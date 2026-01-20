"""
Utilities for genome concatenation and header standardization.

Provides functions for creating pangenome FASTA files with standardized
contig identifiers that enable reliable mapping back to source genomes.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import polars as pl


def extract_accession_from_filename(filename: str) -> str:
    """Extract genome accession from NCBI filename.

    Handles patterns like:
    - GCF_000005845.2_ASM584v2_genomic.fna -> GCF_000005845.2
    - GCA_000001405.15_GRCh38_genomic.fna -> GCA_000001405.15
    - GCF_000005845.2.fna -> GCF_000005845.2

    Args:
        filename: Genome filename (with or without extension)

    Returns:
        Extracted accession string
    """
    # Remove common suffixes first
    stem = filename
    for suffix in [".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta"]:
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break

    # Match RefSeq/GenBank accession pattern
    match = re.match(r"(GCF_\d+\.\d+|GCA_\d+\.\d+)", stem)
    if match:
        return match.group(1)

    # Fallback: use stem before first underscore if it looks like an accession
    parts = stem.split("_")
    if len(parts) >= 2 and parts[0] in ("GCF", "GCA"):
        return f"{parts[0]}_{parts[1]}"

    # Last resort: return the whole stem
    return stem


def concatenate_genomes_with_mapping(
    genome_dir: Path,
    output_fasta: Path,
    contig_mapping_path: Path,
    pattern: str = "*.fna",
) -> tuple[int, int]:
    """Concatenate genomes with standardized headers and create contig mapping.

    Rewrites FASTA headers to a standardized format that enables reliable
    genome identification from BLAST or mapping hits:

        >{accession}|{original_contig_id}

    Creates a contig mapping TSV file with columns:
        - contig_id: Standardized contig identifier
        - genome_accession: Parent genome accession
        - original_header: Original FASTA header (for reference)

    Args:
        genome_dir: Directory containing genome FASTA files
        output_fasta: Output path for concatenated FASTA
        contig_mapping_path: Output path for contig mapping TSV
        pattern: Glob pattern for genome files

    Returns:
        Tuple of (genome_count, contig_count)
    """
    genome_files = sorted(genome_dir.glob(pattern))

    # Try alternative patterns if none found
    if not genome_files:
        for alt_pattern in ["*.fa", "*.fasta", "*.fna.gz", "*.fa.gz", "*.fasta.gz"]:
            genome_files = sorted(genome_dir.glob(alt_pattern))
            if genome_files:
                break

    if not genome_files:
        return 0, 0

    mappings: list[dict[str, str]] = []
    contig_count = 0

    with output_fasta.open("w") as out_f:
        for genome_file in genome_files:
            accession = extract_accession_from_filename(genome_file.name)

            # Handle gzipped files
            open_func = gzip.open if genome_file.suffix == ".gz" else Path.open
            with open_func(genome_file, "rt") as in_f:
                for line in in_f:
                    if line.startswith(">"):
                        original_header = line[1:].strip()
                        # Extract just the contig ID (first word)
                        original_contig = original_header.split()[0]

                        # Create standardized ID: accession|contig
                        standardized_id = f"{accession}|{original_contig}"
                        out_f.write(f">{standardized_id}\n")

                        # Record mapping
                        mappings.append(
                            {
                                "contig_id": standardized_id,
                                "genome_accession": accession,
                                "original_header": original_header,
                            }
                        )
                        contig_count += 1
                    else:
                        out_f.write(line)

    # Write contig mapping
    if mappings:
        df = pl.DataFrame(mappings)
        df.write_csv(contig_mapping_path, separator="\t")

    return len(genome_files), contig_count


def extract_accession_from_contig_id(contig_id: str) -> str:
    """Extract genome accession from standardized contig ID.

    Args:
        contig_id: Standardized ID in format {accession}|{original_contig}

    Returns:
        Genome accession
    """
    if "|" in contig_id:
        return contig_id.split("|")[0]

    # Fallback: try to extract accession pattern
    match = re.match(r"(GCF_\d+\.\d+|GCA_\d+\.\d+)", contig_id)
    if match:
        return match.group(1)

    return contig_id
