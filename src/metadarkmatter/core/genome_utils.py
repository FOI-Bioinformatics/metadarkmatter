"""
Utilities for genome concatenation and header standardization.

Provides functions for creating pangenome FASTA files with standardized
contig identifiers that enable reliable mapping back to source genomes.
"""

from __future__ import annotations

import gzip
import logging
import re
import zlib
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import polars as pl

logger = logging.getLogger(__name__)


def extract_accession_from_filename(filename: str) -> str:
    """Extract genome accession from NCBI filename.

    Handles patterns like:
    - GCF_000005845.2_ASM584v2_genomic.fna -> GCF_000005845.2
    - GCA_000001405.15_GRCh38_genomic.fna -> GCA_000001405.15
    - GCF_000005845.2.fna -> GCF_000005845.2
    - Francisellaceae_Genus_species_GCA_000710735.1.fasta -> GCA_000710735.1

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
    match = re.search(r"(GCF_\d+\.\d+|GCA_\d+\.\d+)", stem)
    if match:
        return match.group(1)

    # Fallback: use stem before first underscore if it looks like an accession
    parts = stem.split("_")
    if len(parts) >= 2 and parts[0] in ("GCF", "GCA"):
        return f"{parts[0]}_{parts[1]}"

    # Last resort: the stem without the NCBI "_genomic" suffix, matching
    # ani_matrix_builder.extract_genome_name_from_path so the BLAST database
    # and ANI matrix labels agree for files without an accession.
    return stem.removesuffix("_genomic")


# Genome FASTA extensions recognised throughout the package.
GENOME_PATTERNS: tuple[str, ...] = (
    "*.fna", "*.fa", "*.fasta", "*.fna.gz", "*.fa.gz", "*.fasta.gz",
)
DEFAULT_GENOME_PATTERN = "*.fna"


def find_genome_files(genome_dir: Path, pattern: str = DEFAULT_GENOME_PATTERN) -> list[Path]:
    """Genome FASTA files in ``genome_dir`` (sorted, regular files only).

    With the default pattern every supported extension is included, so a
    directory that mixes plain and gzipped genomes keeps all of them (the
    alignment database and the ANI matrix must see the same genome set).
    Any other pattern is used as given; when it matches nothing, all
    supported extensions are tried.
    """
    found = set() if pattern == DEFAULT_GENOME_PATTERN else set(genome_dir.glob(pattern))
    if not found:
        # Default pattern, or an explicit pattern that matched nothing.
        found = {f for pat in GENOME_PATTERNS for f in genome_dir.glob(pat)}
    return sorted(f for f in found if f.is_file())


def _read_lines(handle: Any, path: Path) -> Iterator[str]:
    """Yield lines; an unreadable/truncated/corrupt file is a ConfigurationError."""
    from metadarkmatter.core.exceptions import ConfigurationError

    try:
        yield from handle
    except (OSError, EOFError, UnicodeDecodeError, zlib.error) as e:
        raise ConfigurationError(
            f"Could not read genome file {path}: {e}",
            suggestion="Check that the file is complete (e.g. 'gzip -t').",
        ) from e


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
    from metadarkmatter.core.exceptions import ConfigurationError

    genome_files = find_genome_files(genome_dir, pattern)
    if not genome_files:
        return 0, 0

    # Two files with one accession (e.g. GCF_x.1.fna and GCF_x.1_ASM_genomic.fna)
    # would give duplicate sequence IDs: makeblastdb fails obscurely and
    # MMseqs2 silently builds an ambiguous database.
    by_accession: dict[str, Path] = {}
    for genome_file in genome_files:
        acc = extract_accession_from_filename(genome_file.name)
        if acc in by_accession:
            raise ConfigurationError(
                f"Two genome files map to accession {acc}: "
                f"{by_accession[acc].name} and {genome_file.name}.",
                suggestion="Keep one file per genome in the directory.",
            )
        by_accession[acc] = genome_file

    mappings: list[dict[str, str]] = []
    contig_count = 0

    empty_genomes: list[str] = []
    with output_fasta.open("w") as out_f:
        for genome_file in genome_files:
            accession = extract_accession_from_filename(genome_file.name)
            contigs_before = contig_count

            # Handle gzipped files
            open_func = gzip.open if genome_file.suffix == ".gz" else Path.open
            with open_func(genome_file, "rt") as in_f:
                for line_no, line in enumerate(_read_lines(in_f, genome_file), start=1):
                    if contig_count == contigs_before and line.strip() and not line.startswith(">"):
                        # Text before the first header would be appended to
                        # the previous genome's last contig.
                        raise ConfigurationError(
                            f"{genome_file} line {line_no}: data before the first FASTA "
                            "header; the file is not FASTA (or is truncated at the start).",
                        )
                    if line.startswith(">"):
                        original_header = line[1:].strip()
                        if not original_header:
                            raise ConfigurationError(
                                f"Empty FASTA header at line {line_no} of {genome_file}.",
                                suggestion="Give every sequence an ID after '>'.",
                            )
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
                        # Guarantee a newline so a file lacking a trailing
                        # newline does not glue the next genome's header onto
                        # its last sequence line.
                        out_f.write(line if line.endswith("\n") else line + "\n")
            if contig_count == contigs_before:
                empty_genomes.append(genome_file.name)

    if contig_count == 0:
        output_fasta.unlink(missing_ok=True)
        raise ConfigurationError(
            f"No FASTA records in the {len(genome_files)} genome file(s) in {genome_dir}.",
            suggestion="Check that the genome files are complete FASTA files.",
        )
    if empty_genomes:
        logger.warning(
            "%d genome file(s) contain no sequences and were skipped: %s",
            len(empty_genomes), ", ".join(empty_genomes[:5]),
        )

    # Write contig mapping
    if mappings:
        df = pl.DataFrame(mappings)
        df.write_csv(contig_mapping_path, separator="\t")

    return len(genome_files) - len(empty_genomes), contig_count


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
