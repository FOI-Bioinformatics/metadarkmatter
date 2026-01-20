"""
ANI matrix building and conversion utilities.

Provides functions to:
- Parse fastANI output to nested dictionary format
- Parse skani output to nested dictionary format
- Convert parsed output to metadarkmatter CSV matrix format
- Validate ANI matrix coverage against BLAST results
"""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

import polars as pl


def parse_fastani_output(output_path: Path) -> dict[str, dict[str, float]]:
    """Parse fastANI output into nested dictionary.

    FastANI outputs tab-delimited format:
    query_path, ref_path, ANI, matched_fragments, total_fragments

    Note that fastANI does not output pairs with ANI below approximately 80%.
    Missing pairs are assigned a default value of 0.0 to indicate no
    reliable ANI estimate could be obtained.

    Args:
        output_path: Path to fastANI tabular output file.

    Returns:
        Nested dict mapping genome1 -> genome2 -> ANI value.
        All values are symmetric and diagonal entries are 100.0.
    """
    ani_dict: dict[str, dict[str, float]] = defaultdict(dict)
    all_genomes: set[str] = set()

    with output_path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                continue

            query_path, ref_path, ani_str = parts[:3]

            # Extract genome name from path
            query_name = extract_genome_name_from_path(query_path)
            ref_name = extract_genome_name_from_path(ref_path)

            ani_value = float(ani_str)

            all_genomes.add(query_name)
            all_genomes.add(ref_name)

            # Store both directions for symmetry
            ani_dict[query_name][ref_name] = ani_value
            ani_dict[ref_name][query_name] = ani_value

    # Fill diagonal with 100.0 and missing pairs with 0.0
    for genome in all_genomes:
        ani_dict[genome][genome] = 100.0
        for other in all_genomes:
            if other not in ani_dict[genome]:
                ani_dict[genome][other] = 0.0

    return dict(ani_dict)


def parse_skani_output(
    output_path: Path,
    full_matrix: bool = True,
) -> dict[str, dict[str, float]]:
    """Parse skani triangle output into nested dictionary.

    Args:
        output_path: Path to skani output file.
        full_matrix: If True, parse as full matrix format (from --full-matrix).
            If False, parse as pairwise/sparse format.

    Returns:
        Nested dict mapping genome1 -> genome2 -> ANI value.
    """
    if full_matrix:
        return _parse_skani_full_matrix(output_path)
    return _parse_skani_pairwise(output_path)


def _parse_skani_full_matrix(output_path: Path) -> dict[str, dict[str, float]]:
    """Parse skani full matrix format (PHYLIP-like).

    Format:
        N (number of genomes)
        genome1 val11 val12 val13 ...
        genome2 val21 val22 val23 ...
        ...
    """
    ani_dict: dict[str, dict[str, float]] = {}

    with output_path.open() as f:
        lines = [line.strip() for line in f if line.strip()]

    if not lines:
        return {}

    # First line is number of genomes
    n_genomes = int(lines[0])

    # Following lines contain genome name and values
    genome_names: list[str] = []
    values: list[list[float]] = []

    for line in lines[1 : n_genomes + 1]:
        parts = line.split()
        genome_name = extract_genome_name_from_path(parts[0])
        genome_names.append(genome_name)
        row_values = [float(v) for v in parts[1:]]
        values.append(row_values)

    # Build nested dict
    for i, genome1 in enumerate(genome_names):
        ani_dict[genome1] = {}
        for j, genome2 in enumerate(genome_names):
            ani_dict[genome1][genome2] = values[i][j]

    return ani_dict


def _parse_skani_pairwise(output_path: Path) -> dict[str, dict[str, float]]:
    """Parse skani pairwise (sparse/default) output format.

    Format (tab-delimited):
        Ref_file Query_file ANI Align_fraction_ref Align_fraction_query ...
    """
    ani_dict: dict[str, dict[str, float]] = defaultdict(dict)
    all_genomes: set[str] = set()

    with output_path.open() as f:
        for line in f:
            line = line.strip()
            # Skip empty lines and header
            if not line or line.startswith("Ref"):
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                continue

            ref_path, query_path, ani_str = parts[:3]

            ref_name = extract_genome_name_from_path(ref_path)
            query_name = extract_genome_name_from_path(query_path)
            ani_value = float(ani_str)

            all_genomes.add(ref_name)
            all_genomes.add(query_name)

            # Store both directions for symmetry
            ani_dict[ref_name][query_name] = ani_value
            ani_dict[query_name][ref_name] = ani_value

    # Fill diagonal with 100.0
    for genome in all_genomes:
        ani_dict[genome][genome] = 100.0

    return dict(ani_dict)


def extract_genome_name_from_path(path_str: str) -> str:
    """Extract genome accession from file path or BLAST sseqid.

    Handles common formats:
    - File paths: /path/to/GCF_000123456.1_ASM123v1_genomic.fna -> GCF_000123456.1
    - BLAST sseqid with pipe format: GCF_000123456.1|NZ_CP007557.1 -> GCF_000123456.1
    - Custom names: custom_genome|contig1 -> custom_genome
    - Simple filenames: genome.fasta -> genome

    Args:
        path_str: Path string or BLAST sseqid to extract genome name from.

    Returns:
        Extracted genome accession or cleaned name.
    """
    # Get filename without path
    filename = Path(path_str).name

    # Handle standardized pipe format (from BLAST sseqid): accession|contig
    # This should be checked early as it's the format from BLAST output
    if "|" in filename:
        filename = filename.split("|")[0]

    # Remove common extensions (order matters for compound extensions)
    for ext in [".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta"]:
        if filename.endswith(ext):
            filename = filename[: -len(ext)]
            break

    # Remove _genomic suffix common in NCBI downloads
    filename = filename.removesuffix("_genomic")

    # Try to extract RefSeq/GenBank accession pattern
    match = re.match(r"(GCF_\d+\.\d+|GCA_\d+\.\d+)", filename)
    if match:
        return match.group(1)

    # Fallback: return cleaned filename
    return filename


def ani_dict_to_csv(
    ani_dict: dict[str, dict[str, float]],
    output_path: Path,
    compress: bool = False,
) -> int:
    """Convert ANI dictionary to metadarkmatter CSV matrix format.

    Output format:
    - First column named 'genome' contains genome names
    - Header row contains genome names
    - Values are ANI percentages (0-100)
    - Square, symmetric matrix

    Args:
        ani_dict: Nested dict from parse functions.
        output_path: Output CSV file path (use .csv.gz for compressed).
        compress: If True, gzip the output file.

    Returns:
        Number of genomes in matrix.

    Raises:
        ValueError: If ANI dictionary is empty.
    """
    import gzip

    if not ani_dict:
        msg = "ANI dictionary is empty"
        raise ValueError(msg)

    # Sort genome names for consistent output
    genomes = sorted(ani_dict.keys())
    n_genomes = len(genomes)

    # Build data dict for Polars DataFrame
    data: dict[str, list[str] | list[float]] = {"genome": genomes}
    for genome in genomes:
        data[genome] = [ani_dict[g].get(genome, 0.0) for g in genomes]

    # Create DataFrame
    df = pl.DataFrame(data)

    # Write to file (compressed or not)
    if compress or str(output_path).endswith(".gz"):
        # Write to gzipped file
        csv_bytes = df.write_csv().encode("utf-8")
        with gzip.open(output_path, "wb") as f:
            f.write(csv_bytes)
    else:
        df.write_csv(output_path)

    return n_genomes


def validate_ani_coverage(
    ani_genomes: set[str],
    blast_genomes: set[str],
) -> tuple[int, int, float, set[str]]:
    """Validate ANI matrix coverage against genomes found in BLAST results.

    Args:
        ani_genomes: Set of genome names in the ANI matrix.
        blast_genomes: Set of genome names found in BLAST results.

    Returns:
        Tuple of:
        - matched_count: Number of BLAST genomes found in ANI matrix
        - total_blast: Total number of unique genomes in BLAST results
        - coverage_pct: Percentage of BLAST genomes covered by ANI matrix
        - missing_genomes: Set of BLAST genomes not in ANI matrix
    """
    matched = blast_genomes & ani_genomes
    missing = blast_genomes - ani_genomes

    total = len(blast_genomes)
    matched_count = len(matched)
    coverage_pct = (100.0 * matched_count / total) if total > 0 else 0.0

    return matched_count, total, coverage_pct, missing
