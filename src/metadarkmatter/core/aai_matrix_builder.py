"""
AAI matrix building and computation utilities.

Average Amino Acid Identity (AAI) provides a more reliable metric than ANI
for genus-level classification. While ANI becomes unreliable below approximately
80% identity, AAI maintains accuracy at the genus level where AAI values
typically range from 45-65%.

AAI genus boundaries (Riesco & Trujillo 2024):
- Same genus: AAI > 65%
- Genus boundary zone: AAI 58-65%
- Different genus: AAI < 58%

This module provides:
- AAIMatrix: High-performance AAI matrix with NumPy backend
- compute_aai_matrix: Compute AAI from genome protein files
- Reciprocal best hit (RBH) filtering for ortholog identification
"""

from __future__ import annotations

import gzip
import logging
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import polars as pl

from metadarkmatter.core.constants import AAI_DEFAULT_UNRELATED

logger = logging.getLogger(__name__)


class AAIMatrix:
    """
    High-performance AAI matrix using NumPy arrays with integer indexing.

    Mirrors the ANIMatrix interface but stores Average Amino Acid Identity
    values for genus-level classification decisions.

    Performance characteristics match ANIMatrix:
    - Memory: 8MB for 1000 genomes
    - Lookup: O(1) array access after genome name lookup

    Missing AAI values (0.0 in the matrix) are treated as distant organisms
    and return the default_aai value (typically 50%) instead of 0.0.
    """

    __slots__ = ("_aai_array", "_default_aai", "_genome_to_idx", "_genomes", "_num_genomes")

    def __init__(
        self,
        aai_dict: dict[str, dict[str, float]],
        default_aai: float = AAI_DEFAULT_UNRELATED,
    ) -> None:
        """
        Initialize AAI matrix from nested dictionary.

        Args:
            aai_dict: Nested dict {genome1: {genome2: aai_value}}
            default_aai: Default AAI value for missing pairs (default: 50.0).
                When Diamond cannot compute AAI between distant genomes,
                the value is stored as 0.0. This parameter specifies what value
                to return instead.
        """
        self._default_aai = default_aai

        # Sort genomes for consistent indexing
        self._genomes: tuple[str, ...] = tuple(sorted(aai_dict.keys()))
        self._num_genomes: int = len(self._genomes)

        # Create genome -> index mapping for O(1) lookups
        self._genome_to_idx: dict[str, int] = {
            genome: idx for idx, genome in enumerate(self._genomes)
        }

        # Build NumPy array: float32 is sufficient for AAI values (0-100)
        self._aai_array: np.ndarray = np.zeros(
            (self._num_genomes, self._num_genomes),
            dtype=np.float32,
        )

        # Fill diagonal with 100.0 (self-AAI)
        np.fill_diagonal(self._aai_array, 100.0)

        # Populate matrix from dictionary
        for genome1, inner_dict in aai_dict.items():
            i = self._genome_to_idx[genome1]
            for genome2, aai_value in inner_dict.items():
                if genome2 in self._genome_to_idx:
                    j = self._genome_to_idx[genome2]
                    self._aai_array[i, j] = aai_value

    @classmethod
    def from_file(cls, path: Path, default_aai: float = AAI_DEFAULT_UNRELATED) -> AAIMatrix:
        """
        Load AAI matrix from CSV file.

        Expected format (same as ANI matrix):
        - First column named 'genome' contains genome names
        - Header row contains genome names
        - Values are AAI percentages (0-100)
        - Square, symmetric matrix

        Args:
            path: Path to AAI matrix CSV
            default_aai: Default AAI value for missing pairs (default: 50.0)

        Returns:
            AAIMatrix instance
        """
        aai_dict = parse_aai_csv(path)
        return cls(aai_dict, default_aai=default_aai)

    @property
    def genomes(self) -> set[str]:
        """Set of genome names in the matrix."""
        return set(self._genomes)

    def __len__(self) -> int:
        """Number of genomes in the matrix."""
        return self._num_genomes

    def get_aai(self, genome1: str, genome2: str) -> float:
        """
        Get AAI value between two genomes.

        Args:
            genome1: First genome identifier
            genome2: Second genome identifier

        Returns:
            AAI value (0-100). Returns default_aai for missing pairs
            where both genomes exist in the matrix. Returns 0.0 only
            if a genome is not in the matrix at all.
        """
        # Fast path: same genome
        if genome1 == genome2:
            return 100.0

        # Get indices
        i = self._genome_to_idx.get(genome1)
        if i is None:
            return 0.0

        j = self._genome_to_idx.get(genome2)
        if j is None:
            return 0.0

        # O(1) array access
        aai_value = float(self._aai_array[i, j])

        # If AAI is 0.0 (not computed), return default for distant organisms
        if aai_value == 0.0:
            return self._default_aai

        return aai_value

    def get_aai_by_idx(self, idx1: int, idx2: int) -> float:
        """
        Get AAI value by pre-computed integer indices.

        Fastest lookup method when genome indices are already known.

        Args:
            idx1: First genome index
            idx2: Second genome index

        Returns:
            AAI value (0-100)
        """
        return float(self._aai_array[idx1, idx2])

    def get_genome_idx(self, genome: str) -> int | None:
        """
        Get integer index for a genome name.

        Args:
            genome: Genome identifier

        Returns:
            Integer index, or None if genome not found
        """
        return self._genome_to_idx.get(genome)

    def has_genome(self, genome: str) -> bool:
        """Check if genome is present in AAI matrix."""
        return genome in self._genome_to_idx

    def memory_usage_bytes(self) -> int:
        """Estimate memory usage in bytes."""
        array_bytes = self._aai_array.nbytes
        dict_bytes = self._num_genomes * 100  # ~100 bytes per entry
        return array_bytes + dict_bytes


def parse_aai_csv(path: Path) -> dict[str, dict[str, float]]:
    """
    Parse AAI matrix from CSV file into nested dictionary.

    Expected format:
    - First column named 'genome' contains genome names
    - Header row contains genome names as column headers
    - Values are AAI percentages (0-100)

    Args:
        path: Path to AAI matrix CSV file

    Returns:
        Nested dict mapping genome1 -> genome2 -> AAI value
    """
    # Read CSV with Polars
    df = pl.read_csv(path)

    # Get genome names from first column
    if "genome" not in df.columns:
        msg = "AAI matrix CSV must have 'genome' column"
        raise ValueError(msg)

    genome_names = df["genome"].to_list()

    # Build nested dict
    aai_dict: dict[str, dict[str, float]] = {}

    for row_genome in genome_names:
        aai_dict[row_genome] = {}
        row = df.filter(pl.col("genome") == row_genome).drop("genome")

        for col_genome in genome_names:
            if col_genome in row.columns:
                value = row[col_genome].item()
                aai_dict[row_genome][col_genome] = float(value) if value is not None else 0.0

    return aai_dict


def extract_genome_from_protein_header(header: str) -> str:
    """
    Extract genome accession from protein FASTA header.

    Handles common formats:
    - NCBI format: >WP_123456789.1 MULTISPECIES: protein [genus species]
    - Prodigal format: >GCF_000123456.1_1 # start # end # strand # ID=...
    - RefSeq format: >NP_123456.1 protein [Genus species]
    - Custom: >GCF_000123456.1|protein_id description

    The genome name is extracted from the prefix before the first space,
    underscore after accession, or pipe character.

    Args:
        header: FASTA header line (with or without leading '>')

    Returns:
        Genome accession (e.g., 'GCF_000123456.1')
    """
    import re

    # Remove leading '>' if present
    if header.startswith(">"):
        header = header[1:]

    # Get the ID portion (before first space)
    protein_id = header.split()[0] if header else ""

    # Handle pipe-separated format: genome|protein_id
    if "|" in protein_id:
        parts = protein_id.split("|")
        # Check if first part looks like a genome accession
        if re.match(r"^GC[AF]_\d+\.\d+", parts[0]):
            return parts[0]
        # Otherwise return second part's prefix
        protein_id = parts[-1]

    # Handle Prodigal format: genome_contig_proteinnum
    # e.g., GCF_000123456.1_1 or GCF_000123456.1_NZ_CP007557.1_1
    match = re.match(r"^(GC[AF]_\d+\.\d+)", protein_id)
    if match:
        return match.group(1)

    # Handle underscore-separated with numeric suffix: genome_123
    parts = protein_id.rsplit("_", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]

    # Fallback: return as-is (might be a genome name already)
    return protein_id


def concatenate_protein_files(
    protein_dir: Path,
    output_path: Path,
    patterns: list[str] | None = None,
) -> dict[str, str]:
    """
    Concatenate genome protein files into a single FASTA for Diamond.

    Rewrites headers to include genome accession prefix for tracking.
    Format: >genome_accession|original_protein_id

    Args:
        protein_dir: Directory containing protein FASTA files (.faa)
        output_path: Output concatenated FASTA path
        patterns: Glob patterns for protein files (default: *.faa, *.faa.gz)

    Returns:
        Dict mapping genome_accession -> protein_file_path
    """
    if patterns is None:
        patterns = ["*.faa", "*.faa.gz", "*.fasta", "*.fasta.gz"]

    genome_to_file: dict[str, str] = {}

    # Find all protein files
    protein_files: list[Path] = []
    for pattern in patterns:
        protein_files.extend(protein_dir.glob(pattern))

    if not protein_files:
        msg = f"No protein files found in {protein_dir} with patterns {patterns}"
        raise FileNotFoundError(msg)

    logger.info("Found %d protein files in %s", len(protein_files), protein_dir)

    with output_path.open("w") as out_f:
        for prot_file in sorted(protein_files):
            # Extract genome accession from filename
            genome_name = extract_genome_name_from_filename(prot_file)
            genome_to_file[genome_name] = str(prot_file)

            # Open file (handle gzipped)
            opener = gzip.open if str(prot_file).endswith(".gz") else open
            with opener(prot_file, "rt") as in_f:
                for line in in_f:
                    if line.startswith(">"):
                        # Rewrite header: >genome|original_id
                        original_id = line[1:].split()[0]
                        out_f.write(f">{genome_name}|{original_id}\n")
                    else:
                        out_f.write(line)

    logger.info("Concatenated %d genomes to %s", len(genome_to_file), output_path)
    return genome_to_file


def extract_genome_name_from_filename(path: Path) -> str:
    """
    Extract genome accession from protein filename.

    Handles common NCBI naming patterns:
    - GCF_000123456.1_ASM123v1_protein.faa -> GCF_000123456.1
    - GCF_000123456.1.faa -> GCF_000123456.1
    - custom_genome.faa -> custom_genome

    Args:
        path: Path to protein file

    Returns:
        Genome accession or name
    """
    import re

    filename = path.name

    # Remove extensions
    for ext in [".faa.gz", ".fasta.gz", ".faa", ".fasta"]:
        if filename.endswith(ext):
            filename = filename[: -len(ext)]
            break

    # Remove common suffixes
    for suffix in ["_protein", "_genomic", "_cds_from_genomic"]:
        if filename.endswith(suffix):
            filename = filename[: -len(suffix)]
            break

    # Try to extract RefSeq/GenBank accession
    match = re.match(r"(GC[AF]_\d+\.\d+)", filename)
    if match:
        return match.group(1)

    # Fallback: return cleaned filename
    return filename


def compute_rbh_pairs(
    diamond_output: Path,
    min_identity: float = 30.0,
    min_coverage: float = 0.5,
) -> pl.DataFrame:
    """
    Compute reciprocal best hits (RBH) from Diamond blastp output.

    RBH pairs are the basis for AAI calculation - they represent
    putative orthologous proteins between two genomes.

    Args:
        diamond_output: Path to Diamond blastp output (format 6)
        min_identity: Minimum percent identity to consider (default: 30%)
        min_coverage: Minimum query coverage to consider (default: 50%)

    Returns:
        DataFrame with columns: genome1, genome2, protein1, protein2, pident
        containing reciprocal best hit pairs
    """
    # Read Diamond output
    # Format: qseqid sseqid pident length qlen slen
    df = pl.read_csv(
        diamond_output,
        separator="\t",
        has_header=False,
        new_columns=["qseqid", "sseqid", "pident", "length", "qlen", "slen"],
    )

    # Extract genome names from protein IDs (format: genome|protein_id)
    df = df.with_columns([
        pl.col("qseqid").str.split("|").list.first().alias("query_genome"),
        pl.col("qseqid").str.split("|").list.last().alias("query_protein"),
        pl.col("sseqid").str.split("|").list.first().alias("subject_genome"),
        pl.col("sseqid").str.split("|").list.last().alias("subject_protein"),
    ])

    # Filter by identity and coverage
    df = df.filter(
        (pl.col("pident") >= min_identity)
        & ((pl.col("length") / pl.col("qlen")) >= min_coverage)
    )

    # Remove self-hits (same genome)
    df = df.filter(pl.col("query_genome") != pl.col("subject_genome"))

    # Find best hit for each query protein to each target genome
    best_hits = (
        df.group_by(["query_genome", "query_protein", "subject_genome"])
        .agg([
            pl.col("subject_protein").sort_by("pident", descending=True).first(),
            pl.col("pident").max(),
        ])
        .rename({"pident": "best_pident"})
    )

    # Find reciprocal best hits
    # Join forward and reverse directions
    forward = best_hits.select([
        pl.col("query_genome").alias("genome1"),
        pl.col("query_protein").alias("protein1"),
        pl.col("subject_genome").alias("genome2"),
        pl.col("subject_protein").alias("protein2"),
        pl.col("best_pident").alias("pident_forward"),
    ])

    reverse = best_hits.select([
        pl.col("query_genome").alias("genome2"),
        pl.col("query_protein").alias("protein2"),
        pl.col("subject_genome").alias("genome1"),
        pl.col("subject_protein").alias("protein1"),
        pl.col("best_pident").alias("pident_reverse"),
    ])

    # Join to find pairs where A->B best hit and B->A best hit match
    rbh = forward.join(
        reverse,
        on=["genome1", "protein1", "genome2", "protein2"],
        how="inner",
    )

    # Average the two identity values for the pair
    rbh = rbh.with_columns([
        ((pl.col("pident_forward") + pl.col("pident_reverse")) / 2.0).alias("pident"),
    ]).select(["genome1", "genome2", "protein1", "protein2", "pident"])

    return rbh


def compute_aai_from_rbh(rbh_df: pl.DataFrame) -> dict[str, dict[str, float]]:
    """
    Compute AAI matrix from reciprocal best hit pairs.

    AAI is calculated as the mean percent identity of all RBH pairs
    between two genomes.

    Args:
        rbh_df: DataFrame with RBH pairs (from compute_rbh_pairs)

    Returns:
        Nested dict mapping genome1 -> genome2 -> AAI value
    """
    # Compute mean identity for each genome pair
    aai_values = (
        rbh_df.group_by(["genome1", "genome2"])
        .agg([
            pl.col("pident").mean().alias("aai"),
            pl.col("pident").len().alias("num_rbh"),
        ])
    )

    # Build nested dict
    aai_dict: dict[str, dict[str, float]] = defaultdict(dict)
    all_genomes: set[str] = set()

    for row in aai_values.iter_rows(named=True):
        g1, g2 = row["genome1"], row["genome2"]
        aai = row["aai"]

        all_genomes.add(g1)
        all_genomes.add(g2)

        # Store both directions for symmetry
        aai_dict[g1][g2] = aai
        aai_dict[g2][g1] = aai

    # Fill diagonal with 100.0
    for genome in all_genomes:
        aai_dict[genome][genome] = 100.0

    return dict(aai_dict)


def aai_dict_to_csv(
    aai_dict: dict[str, dict[str, float]],
    output_path: Path,
    compress: bool = False,
) -> int:
    """
    Convert AAI dictionary to CSV matrix format.

    Output format matches ANI matrix:
    - First column named 'genome' contains genome names
    - Header row contains genome names
    - Values are AAI percentages (0-100)
    - Square, symmetric matrix

    Args:
        aai_dict: Nested dict from compute functions
        output_path: Output CSV file path
        compress: If True, gzip the output file

    Returns:
        Number of genomes in matrix

    Raises:
        ValueError: If AAI dictionary is empty
    """
    if not aai_dict:
        msg = "AAI dictionary is empty"
        raise ValueError(msg)

    # Sort genome names for consistent output
    genomes = sorted(aai_dict.keys())
    n_genomes = len(genomes)

    # Build data dict for Polars DataFrame
    data: dict[str, list[str] | list[float]] = {"genome": genomes}
    for genome in genomes:
        data[genome] = [aai_dict[g].get(genome, 0.0) for g in genomes]

    # Create DataFrame
    df = pl.DataFrame(data)

    # Write to file
    if compress or str(output_path).endswith(".gz"):
        csv_bytes = df.write_csv().encode("utf-8")
        with gzip.open(output_path, "wb") as f:
            f.write(csv_bytes)
    else:
        df.write_csv(output_path)

    return n_genomes


def extract_genome_proteins(
    concat_fasta: Path,
    genome: str,
    output_path: Path,
) -> int:
    """
    Extract proteins for a single genome from concatenated FASTA.

    Args:
        concat_fasta: Path to concatenated protein FASTA
        genome: Genome accession to extract (prefix before |)
        output_path: Output FASTA path for extracted proteins

    Returns:
        Number of proteins extracted
    """
    prefix = f">{genome}|"
    count = 0
    writing = False

    with concat_fasta.open() as in_f, output_path.open("w") as out_f:
        for line in in_f:
            if line.startswith(">"):
                writing = line.startswith(prefix)
                if writing:
                    count += 1
            if writing:
                out_f.write(line)

    return count


def compute_aai_matrix(
    protein_dir: Path,
    output_path: Path,
    threads: int = 4,
    min_identity: float = 30.0,
    min_coverage: float = 0.5,
    evalue: float = 1e-5,
    keep_intermediates: bool = False,
    batch_mode: bool = True,
    sensitive: bool = False,
) -> dict[str, dict[str, float]]:
    """
    Compute AAI matrix from genome protein files.

    This is the main entry point for AAI computation. It:
    1. Concatenates all protein files into a single FASTA
    2. Creates a Diamond database
    3. Runs Diamond blastp (batch mode queries each genome separately)
    4. Computes reciprocal best hits (RBH)
    5. Calculates AAI from mean RBH identity per genome pair
    6. Writes the AAI matrix to CSV

    Args:
        protein_dir: Directory containing protein FASTA files (.faa)
        output_path: Output AAI matrix CSV path
        threads: Number of CPU threads for Diamond
        min_identity: Minimum percent identity for hits (default: 30%)
        min_coverage: Minimum query coverage for hits (default: 50%)
        evalue: E-value threshold for Diamond (default: 1e-5)
        keep_intermediates: Keep intermediate files (default: False)
        batch_mode: Process each genome separately to avoid Diamond
            dropping queries in large all-vs-all runs (default: True)
        sensitive: Use Diamond sensitive mode for divergent sequences.
            Slower but more accurate for distant homologs (default: False)

    Returns:
        Nested dict mapping genome1 -> genome2 -> AAI value

    Raises:
        FileNotFoundError: If no protein files found
        ToolNotFoundError: If Diamond is not installed
    """
    from metadarkmatter.external.diamond import Diamond

    diamond = Diamond()

    # Create temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)

        # Step 1: Concatenate protein files
        concat_fasta = tmpdir_path / "all_proteins.faa"
        logger.info("Concatenating protein files...")
        genome_map = concatenate_protein_files(protein_dir, concat_fasta)
        logger.info("Found %d genomes", len(genome_map))

        # Step 2: Create Diamond database
        db_path = tmpdir_path / "proteins"
        logger.info("Creating Diamond database...")
        diamond.make_database(concat_fasta, db_path)

        # Step 3: Run Diamond blastp
        blastp_output = tmpdir_path / "all_vs_all.tsv"

        if batch_mode:
            # Batch mode: query each genome separately to avoid Diamond
            # silently dropping queries when processing large files
            logger.info(
                "Running Diamond blastp in batch mode (%d genomes, threads=%d)...",
                len(genome_map),
                threads,
            )

            genome_query_dir = tmpdir_path / "genome_queries"
            genome_query_dir.mkdir()
            genome_output_dir = tmpdir_path / "genome_outputs"
            genome_output_dir.mkdir()

            genomes = sorted(genome_map.keys())

            # Helper function to process a single genome
            def process_genome(genome: str, genome_index: int) -> tuple[str, int, bool]:
                """Extract proteins and run Diamond for one genome."""
                query_file = genome_query_dir / f"{genome}.faa"
                n_proteins = extract_genome_proteins(concat_fasta, genome, query_file)

                if n_proteins == 0:
                    logger.warning("No proteins found for genome %s", genome)
                    return (genome, 0, False)

                output_file = genome_output_dir / f"{genome}.tsv"
                logger.info(
                    "  [%d/%d] %s (%d proteins)",
                    genome_index,
                    len(genomes),
                    genome,
                    n_proteins,
                )

                # Reduce threads per job to allow parallel processing
                # If threads=8, use threads_per_job=2 with max_workers=4
                threads_per_job = max(1, threads // 4)

                diamond.blastp(
                    query=query_file,
                    database=db_path,
                    output=output_file,
                    threads=threads_per_job,
                    evalue=evalue,
                    min_identity=min_identity,
                    sensitive=sensitive,
                    capture_output=False,  # Prevent pipe buffer blocking
                )
                return (genome, n_proteins, True)

            # Process genomes in parallel (max_workers=4 for threads=8)
            max_workers = max(1, threads // 2)
            logger.info("Processing %d genomes with %d parallel workers...", len(genomes), max_workers)

            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                # Submit all jobs
                future_to_genome = {
                    executor.submit(process_genome, genome, i): genome
                    for i, genome in enumerate(genomes, 1)
                }

                # Wait for completion
                completed = 0
                for future in as_completed(future_to_genome):
                    genome, n_proteins, success = future.result()
                    completed += 1
                    if success:
                        logger.info("Completed %d/%d genomes", completed, len(genomes))

            # Merge all output files
            logger.info("Merging batch results...")
            with blastp_output.open("w") as out_f:
                for genome in genomes:
                    output_file = genome_output_dir / f"{genome}.tsv"
                    if output_file.exists():
                        with output_file.open() as in_f:
                            for line in in_f:
                                out_f.write(line)

        else:
            # Original all-vs-all mode (may drop queries for large datasets)
            logger.info("Running Diamond blastp (threads=%d)...", threads)
            diamond.blastp(
                query=concat_fasta,
                database=db_path,
                output=blastp_output,
                threads=threads,
                evalue=evalue,
                min_identity=min_identity,
                sensitive=sensitive,
            )

        # Step 4: Compute reciprocal best hits
        logger.info("Computing reciprocal best hits...")
        rbh_df = compute_rbh_pairs(
            blastp_output,
            min_identity=min_identity,
            min_coverage=min_coverage,
        )
        logger.info("Found %d RBH pairs", len(rbh_df))

        # Step 5: Calculate AAI from RBH
        logger.info("Computing AAI matrix...")
        aai_dict = compute_aai_from_rbh(rbh_df)

        # Step 6: Write to CSV
        n_genomes = aai_dict_to_csv(aai_dict, output_path)
        logger.info("Wrote AAI matrix for %d genomes to %s", n_genomes, output_path)

        # Optionally keep intermediates
        if keep_intermediates:
            import shutil

            intermediate_dir = output_path.parent / "aai_intermediates"
            intermediate_dir.mkdir(exist_ok=True)
            shutil.copy(blastp_output, intermediate_dir / "diamond_hits.tsv")
            rbh_df.write_csv(intermediate_dir / "rbh_pairs.csv")
            logger.info("Intermediate files saved to %s", intermediate_dir)

    return aai_dict
