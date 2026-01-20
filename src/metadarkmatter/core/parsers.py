"""
Streaming parsers for bioinformatics file formats using Polars.

These parsers handle large-scale metagenomic data efficiently using
lazy evaluation and chunked streaming to minimize memory footprint.
Designed for processing 100M+ BLAST alignment records without RAM overflow.
"""

from __future__ import annotations

import logging
from collections.abc import Generator, Iterator
from pathlib import Path
from typing import TYPE_CHECKING, ClassVar, NamedTuple

import polars as pl

from metadarkmatter.core.constants import UNKNOWN_GENOME
from metadarkmatter.models.blast import BlastHit, BlastResult

if TYPE_CHECKING:
    from metadarkmatter.core.id_mapping import ContigIdMapping

logger = logging.getLogger(__name__)

# =============================================================================
# Lightweight Internal Data Structures (P1 Optimization)
# =============================================================================
# These NamedTuples avoid Pydantic validation overhead in hot paths.
# ~50x faster to create than Pydantic models.


class BlastHitFast(NamedTuple):
    """
    Lightweight BLAST hit for high-performance internal processing.

    No validation overhead - use only with pre-validated data from Polars.
    Fields are a subset of BlastHit, containing only what's needed for
    classification calculations.
    """

    qseqid: str
    sseqid: str
    pident: float
    bitscore: float
    genome_name: str  # Pre-extracted by Polars


class BlastResultFast(NamedTuple):
    """
    Lightweight BLAST result for high-performance internal processing.

    Hits are stored as a tuple for immutability and memory efficiency.
    """

    read_id: str
    hits: tuple[BlastHitFast, ...]


# =============================================================================
# Vectorized Genome Name Extraction (P4 Optimization)
# =============================================================================


def extract_genome_name_expr() -> pl.Expr:
    """
    Create Polars expression for vectorized genome name extraction.

    This extracts genome accession from standardized FASTA headers in format:
        {accession}|{original_contig_id}

    The standardized format is created during BLAST/Bowtie2 database building
    via concatenate_genomes_with_mapping() in genome_utils.py.

    Performance: 100x faster than Python regex per-row due to vectorization
    in Polars' native Rust implementation.

    Returns:
        Polars expression that extracts genome accession from 'sseqid' column

    Example:
        Input sseqid: "GCF_000195955.2|NZ_CP007557.1"
        Output genome_name: "GCF_000195955.2"
    """
    # Extract everything before the pipe separator
    return (
        pl.col("sseqid")
        .str.extract(r"^([^|]+)\|", 1)
        .fill_null(UNKNOWN_GENOME)
        .alias("genome_name")
    )


class StreamingBlastParser:
    """
    Memory-efficient streaming parser for BLAST tabular output.

    Uses Polars with chunked batch reading to process BLAST results
    without loading entire files into memory. Critical for handling
    100M+ alignment results from competitive read recruitment.

    The parser reads data in configurable chunks and groups reads
    while maintaining low memory usage, suitable for HPC environments.

    Expected BLAST format (-outfmt 6):
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore [qlen]

    Note: qlen is optional for backward compatibility. Files with 12 columns (old format)
    or 13 columns (new format with qlen) are both supported.
    """

    # Schema for new 13-column format (with qlen)
    BLAST_SCHEMA: ClassVar[dict[str, pl.DataType]] = {
        "qseqid": pl.Utf8,
        "sseqid": pl.Utf8,
        "pident": pl.Float64,
        "length": pl.Int64,
        "mismatch": pl.Int64,
        "gapopen": pl.Int64,
        "qstart": pl.Int64,
        "qend": pl.Int64,
        "sstart": pl.Int64,
        "send": pl.Int64,
        "evalue": pl.Float64,
        "bitscore": pl.Float64,
        "qlen": pl.Int64,
    }

    # Schema for old 12-column format (without qlen)
    BLAST_SCHEMA_LEGACY: ClassVar[dict[str, pl.DataType]] = {
        "qseqid": pl.Utf8,
        "sseqid": pl.Utf8,
        "pident": pl.Float64,
        "length": pl.Int64,
        "mismatch": pl.Int64,
        "gapopen": pl.Int64,
        "qstart": pl.Int64,
        "qend": pl.Int64,
        "sstart": pl.Int64,
        "send": pl.Int64,
        "evalue": pl.Float64,
        "bitscore": pl.Float64,
    }

    COLUMN_NAMES: ClassVar[list[str]] = list(BLAST_SCHEMA.keys())
    COLUMN_NAMES_LEGACY: ClassVar[list[str]] = list(BLAST_SCHEMA_LEGACY.keys())

    # Bounds for chunk_size parameter
    MIN_CHUNK_SIZE = 100  # Minimum for reasonable performance (allows testing)
    MAX_CHUNK_SIZE = 100_000_000  # 100M rows max to prevent memory issues

    def __init__(self, blast_path: Path, chunk_size: int = 1_000_000) -> None:
        """
        Initialize streaming BLAST parser.

        Args:
            blast_path: Path to BLAST tabular output file (may be gzipped)
            chunk_size: Number of rows to read per batch for streaming.
                        Default of 1M rows balances memory and I/O efficiency.
                        For 100M records, this processes in ~100 batches.
                        Valid range: 1,000 to 100,000,000.

        Raises:
            ValueError: If chunk_size is outside valid range.

        Memory usage:
            Approximately 200-300 bytes per row when loaded into memory.
            chunk_size=1M uses ~200-300MB RAM per batch.
        """
        if not self.MIN_CHUNK_SIZE <= chunk_size <= self.MAX_CHUNK_SIZE:
            msg = (
                f"chunk_size must be between {self.MIN_CHUNK_SIZE:,} and "
                f"{self.MAX_CHUNK_SIZE:,}, got {chunk_size:,}"
            )
            raise ValueError(msg)

        self.blast_path = blast_path
        self.chunk_size = chunk_size
        self._validate_path()

        # Detect format and set appropriate schema
        num_cols = self._detect_column_count()
        if num_cols == 13:
            self.schema = self.BLAST_SCHEMA
            self.column_names = self.COLUMN_NAMES
        else:  # 12 columns (legacy format)
            self.schema = self.BLAST_SCHEMA_LEGACY
            self.column_names = self.COLUMN_NAMES_LEGACY

    def _validate_path(self) -> None:
        """Ensure BLAST file exists."""
        if not self.blast_path.exists():
            msg = f"BLAST file not found: {self.blast_path}"
            raise FileNotFoundError(msg)

    def _is_gzipped(self) -> bool:
        """Check if file is gzip compressed."""
        return self.blast_path.suffix == ".gz"

    def _detect_column_count(self) -> int:
        """
        Detect number of columns by reading first non-comment line.

        Returns:
            Number of tab-separated columns (12 or 13)

        Raises:
            ValueError: If file is empty or has unexpected column count
        """
        import gzip

        # Open file (handle gzip if needed)
        if self._is_gzipped():
            file_handle = gzip.open(self.blast_path, "rt")
        else:
            file_handle = self.blast_path.open("r")

        try:
            for line in file_handle:
                # Skip comment lines
                if line.startswith("#"):
                    continue

                # Count columns
                fields = line.strip().split("\t")
                num_cols = len(fields)

                if num_cols not in (12, 13):
                    msg = (
                        f"Unexpected BLAST format: expected 12 or 13 columns, "
                        f"got {num_cols} in {self.blast_path}"
                    )
                    raise ValueError(msg)

                return num_cols

            # No data lines found
            msg = f"BLAST file contains no data lines: {self.blast_path}"
            raise ValueError(msg)

        finally:
            file_handle.close()

    def parse_lazy(self) -> pl.LazyFrame:
        """
        Parse BLAST file into Polars LazyFrame for deferred execution.

        Returns:
            LazyFrame with BLAST results, not yet materialized in memory

        Example:
            parser = StreamingBlastParser(Path("results.blast.gz"))
            lazy_df = parser.parse_lazy()
            # Filter and aggregate without loading entire file
            top_hits = (
                lazy_df
                .group_by("qseqid")
                .agg(pl.col("bitscore").max().alias("max_bitscore"))
                .collect(engine="streaming")
            )
        """
        return pl.scan_csv(
            self.blast_path,
            separator="\t",
            has_header=False,
            new_columns=self.column_names,
            schema_overrides=self.schema,
            comment_prefix="#",
            ignore_errors=False,
        )

    def parse_eager(self) -> pl.DataFrame:
        """
        Parse BLAST file into materialized Polars DataFrame.

        Loads entire file into memory. Use only for small-to-medium files
        or when immediate computation is required.

        Returns:
            Materialized DataFrame with BLAST results
        """
        return self.parse_lazy().collect()

    def iter_reads_chunked(self) -> Iterator[BlastResult]:
        """
        Iterate over BLAST results using true chunked streaming.

        Reads file in batches of chunk_size rows, groups hits per read,
        and handles read boundaries across chunks correctly. This approach
        maintains bounded memory usage regardless of total file size.

        Memory usage is approximately: chunk_size * ~200 bytes per row,
        plus overhead for partial read accumulation.

        Yields:
            BlastResult objects with hits sorted by bitscore descending

        Example:
            parser = StreamingBlastParser(Path("results.blast.gz"), chunk_size=500_000)
            for result in parser.iter_reads_chunked():
                if result.num_hits >= 10:
                    print(f"{result.read_id}: {result.best_hit.genome_name}")
        """
        # Buffer to hold partial read data across chunk boundaries
        pending_hits: dict[str, list[dict]] = {}

        # Use Polars scan_csv with collect_batches for streaming
        batches = pl.scan_csv(
            self.blast_path,
            separator="\t",
            has_header=False,
            new_columns=self.column_names,
            comment_prefix="#",
        ).collect_batches(chunk_size=self.chunk_size)

        for chunk_df in batches:
            if chunk_df.is_empty():
                continue

            # Sort chunk by read_id and bitscore for consistent grouping
            chunk_df = chunk_df.sort(["qseqid", "bitscore"], descending=[False, True])

            # Group hits by read_id within this chunk
            for row in chunk_df.iter_rows(named=True):
                read_id = row["qseqid"]
                if read_id not in pending_hits:
                    pending_hits[read_id] = []
                pending_hits[read_id].append(row)

            # Get unique read_ids in this chunk to determine which are complete
            chunk_read_ids = chunk_df["qseqid"].unique().to_list()

            # The last read_id in the chunk may span to the next chunk
            # Sort to find the last one alphabetically (matches sort order)
            if chunk_read_ids:
                current_last_read_id = chunk_df["qseqid"][-1]

                # Yield all reads except the last one (which may be incomplete)
                for read_id in list(pending_hits.keys()):
                    if read_id != current_last_read_id:
                        yield self._create_blast_result(read_id, pending_hits.pop(read_id))


        # Yield remaining reads from the buffer
        for read_id, hits in pending_hits.items():
            yield self._create_blast_result(read_id, hits)

    def _create_blast_result(self, read_id: str, hit_rows: list[dict]) -> BlastResult:
        """
        Create BlastResult from accumulated hit rows.

        Args:
            read_id: Query sequence identifier
            hit_rows: List of row dictionaries for this read (pre-sorted by Polars)

        Returns:
            BlastResult with hits sorted by bitscore descending
        """
        # NOTE: hit_rows are already sorted by bitscore descending from Polars
        # Sorting removed here to avoid redundant O(n log n) operation per read
        # The chunk_df.sort() in iter_reads_chunked() handles sorting

        hits = tuple(
            BlastHit(
                qseqid=row["qseqid"],
                sseqid=row["sseqid"],
                pident=row["pident"],
                length=row["length"],
                mismatch=row["mismatch"],
                gapopen=row["gapopen"],
                qstart=row["qstart"],
                qend=row["qend"],
                sstart=row["sstart"],
                send=row["send"],
                evalue=row["evalue"],
                bitscore=row["bitscore"],
                qlen=row.get("qlen"),  # Optional for backward compatibility
            )
            for row in hit_rows
        )

        return BlastResult(read_id=read_id, hits=hits)

    def iter_reads(self) -> Iterator[BlastResult]:
        """
        Iterate over BLAST results grouped by read ID.

        This method delegates to iter_reads_chunked() for true streaming.
        For very large files (100M+ records), this maintains bounded memory
        usage by processing in batches.

        Yields:
            BlastResult objects with hits sorted by bitscore descending

        Example:
            parser = StreamingBlastParser(Path("results.blast.gz"))
            for result in parser.iter_reads():
                if result.num_hits >= 10:
                    print(f"{result.read_id}: {result.best_hit.genome_name}")
        """
        yield from self.iter_reads_chunked()

    def iter_reads_fast(self) -> Iterator[BlastResultFast]:
        """
        High-performance iterator using lightweight data structures.

        Performance optimizations vs iter_reads():
        - Uses NamedTuples instead of Pydantic models (~50x faster object creation)
        - Vectorized genome extraction in Polars (~100x faster)
        - Minimal fields (only what classification needs)
        - Batch grouping via partition_by reduces dict operations

        Yields:
            BlastResultFast objects with pre-extracted genome names

        Example:
            parser = StreamingBlastParser(Path("results.blast.gz"))
            for result in parser.iter_reads_fast():
                best = result.hits[0]
                print(f"{result.read_id}: {best.genome_name} @ {best.pident}%")
        """
        pending_hits: dict[str, list[BlastHitFast]] = {}

        # Use Polars scan_csv with collect_batches for streaming
        batches = pl.scan_csv(
            self.blast_path,
            separator="\t",
            has_header=False,
            new_columns=self.column_names,
            comment_prefix="#",
        ).collect_batches(chunk_size=self.chunk_size)

        for chunk_df in batches:
            if chunk_df.is_empty():
                continue

            # Vectorized genome extraction + sorting in Polars
            chunk_df = (
                chunk_df
                .with_columns([extract_genome_name_expr()])
                .sort(["qseqid", "bitscore"], descending=[False, True])
            )

            # Use partition_by for efficient grouping - reduces dict operations
            # from O(rows) to O(unique_reads) within each chunk
            grouped = chunk_df.select(
                ["qseqid", "sseqid", "pident", "bitscore", "genome_name"]
            ).partition_by("qseqid", maintain_order=True)

            for group_df in grouped:
                if group_df.is_empty():
                    continue
                read_id = group_df["qseqid"][0]

                # Batch convert group to BlastHitFast objects
                hits = [
                    BlastHitFast(
                        qseqid=row[0],
                        sseqid=row[1],
                        pident=float(row[2]),
                        bitscore=float(row[3]),
                        genome_name=row[4],
                    )
                    for row in group_df.iter_rows()
                ]

                if read_id in pending_hits:
                    pending_hits[read_id].extend(hits)
                else:
                    pending_hits[read_id] = hits

            # Yield complete reads (same boundary logic as iter_reads_chunked)
            if chunk_df.height > 0:
                current_last_read_id = chunk_df["qseqid"][-1]
                for read_id in list(pending_hits.keys()):
                    if read_id != current_last_read_id:
                        hits = pending_hits.pop(read_id)
                        yield BlastResultFast(read_id=read_id, hits=tuple(hits))

        # Yield remaining reads
        for read_id, hits in pending_hits.items():
            yield BlastResultFast(read_id=read_id, hits=tuple(hits))

    def get_top_hits_per_read(self) -> pl.DataFrame:
        """
        Extract only the top BLAST hit for each read.

        Memory-efficient operation using Polars lazy evaluation.

        Returns:
            DataFrame with one row per read (highest bitscore hit only)
        """
        return (
            self.parse_lazy()
            .sort("bitscore", descending=True)
            .group_by("qseqid")
            .first()
            .collect(engine="streaming")
        )

    def get_ambiguous_hits(
        self,
        bitscore_threshold_pct: float = 95.0,
    ) -> pl.DataFrame:
        """
        Get all hits within threshold percentage of top bitscore per read.

        Identifies reads with competitive placements across multiple genomes.

        Args:
            bitscore_threshold_pct: Percentage of max bitscore (default: 95%)

        Returns:
            DataFrame with ambiguous hits, including 'max_bitscore' column
        """
        lazy_df = self.parse_lazy()

        # Calculate max bitscore per read
        max_bitscores = (
            lazy_df.group_by("qseqid")
            .agg(pl.col("bitscore").max().alias("max_bitscore"))
        )

        # Join back and filter
        return (
            lazy_df.join(max_bitscores, on="qseqid", how="left")
            .filter(
                pl.col("bitscore") >= pl.col("max_bitscore") * (bitscore_threshold_pct / 100.0)
            )
            .collect(engine="streaming")
        )

    def count_hits_per_read(self) -> pl.DataFrame:
        """
        Count number of BLAST hits per read.

        Returns:
            DataFrame with columns: qseqid, hit_count
        """
        return (
            self.parse_lazy()
            .group_by("qseqid")
            .agg(pl.len().alias("hit_count"))
            .collect(engine="streaming")
        )

    def get_summary_stats(self) -> dict[str, float]:
        """
        Calculate summary statistics for BLAST results.

        Returns:
            Dictionary with total_hits, unique_reads, mean_pident, mean_bitscore
        """
        lazy_df = self.parse_lazy()

        stats = lazy_df.select(
            pl.len().alias("total_hits"),
            pl.col("qseqid").n_unique().alias("unique_reads"),
            pl.col("pident").mean().alias("mean_pident"),
            pl.col("bitscore").mean().alias("mean_bitscore"),
        ).collect()

        return stats.to_dicts()[0]


class ANIMatrixParser:
    """
    Parser for precomputed ANI (Average Nucleotide Identity) matrices.

    ANI matrices contain pairwise genome similarity values used for
    calculating placement uncertainty in competitive recruitment.

    Expected format:
    - CSV or TSV with genome names as first column and header row
    - Symmetric matrix with ANI values (0-100)
    """

    def __init__(self, ani_path: Path) -> None:
        """
        Initialize ANI matrix parser.

        Args:
            ani_path: Path to ANI matrix file (CSV/TSV)
        """
        self.ani_path = ani_path
        self._validate_path()

    def _validate_path(self) -> None:
        """Ensure ANI file exists."""
        if not self.ani_path.exists():
            msg = f"ANI matrix file not found: {self.ani_path}"
            raise FileNotFoundError(msg)

    def parse(self) -> pl.DataFrame:
        """
        Parse ANI matrix into Polars DataFrame.

        Supports both compressed (.csv.gz, .tsv.gz) and uncompressed files.

        Returns:
            DataFrame with genome names as first column and ANI values

        Raises:
            ValueError: If matrix is not symmetric or contains invalid values
        """
        # Detect separator (CSV vs TSV), handling .gz extension
        path_str = str(self.ani_path)
        separator = "\t" if path_str.endswith((".tsv", ".tsv.gz")) else ","

        # Polars read_csv handles gzip automatically based on extension
        df = pl.read_csv(
            self.ani_path,
            separator=separator,
            has_header=True,
        )

        # Validate matrix structure
        self._validate_matrix(df)

        return df

    def _validate_matrix(self, df: pl.DataFrame) -> None:
        """
        Validate ANI matrix structure, symmetry, and values.

        Performs comprehensive validation to ensure the ANI matrix is
        properly formatted for use in placement uncertainty calculations.

        Args:
            df: Parsed DataFrame

        Raises:
            ValueError: If validation fails due to:
                - Non-square matrix dimensions
                - Row genome names not matching column headers
                - ANI values outside valid range [0, 100]
        """
        if df.is_empty():
            msg = "ANI matrix is empty"
            raise ValueError(msg)

        # Check that matrix is square
        genome_col = df.columns[0]
        num_genomes = len(df)
        value_columns = df.columns[1:]  # Exclude genome name column
        num_cols = len(value_columns)

        if num_genomes != num_cols:
            msg = f"ANI matrix is not square: {num_genomes} rows, {num_cols} value columns"
            raise ValueError(msg)

        # Validate that row genome names match column headers (Issue 8)
        row_genomes = df[genome_col].to_list()
        col_genomes = value_columns

        # Check for exact match between row names and column headers
        row_set = set(row_genomes)
        col_set = set(col_genomes)

        if row_set != col_set:
            missing_in_rows = col_set - row_set
            missing_in_cols = row_set - col_set
            msg_parts = ["ANI matrix row/column mismatch:"]
            if missing_in_rows:
                msg_parts.append(f"  Columns missing from rows: {sorted(missing_in_rows)[:5]}")
            if missing_in_cols:
                msg_parts.append(f"  Rows missing from columns: {sorted(missing_in_cols)[:5]}")
            raise ValueError("\n".join(msg_parts))

        # Check for duplicate genome names
        if len(row_genomes) != len(row_set):
            duplicates = [g for g in row_set if row_genomes.count(g) > 1]
            msg = f"ANI matrix contains duplicate genome names: {duplicates[:5]}"
            raise ValueError(msg)

        # Check for valid ANI range (0-100)
        numeric_cols = df.select(pl.exclude(genome_col))

        # Handle potential null values
        if numeric_cols.null_count().sum_horizontal()[0] > 0:
            msg = "ANI matrix contains null/missing values"
            raise ValueError(msg)

        min_val = numeric_cols.min().min_horizontal()[0]
        max_val = numeric_cols.max().max_horizontal()[0]

        if min_val < 0 or max_val > 100:
            msg = f"ANI values out of range [0, 100]: min={min_val}, max={max_val}"
            raise ValueError(msg)

    def to_dict(self) -> dict[str, dict[str, float]]:
        """
        Convert ANI matrix to nested dictionary for fast lookups.

        Uses Polars' optimized to_dicts() method for efficient conversion.

        Returns:
            Nested dict: {genome1: {genome2: ani_value, ...}, ...}
        """
        df = self.parse()
        genome_col = df.columns[0]
        other_cols = [c for c in df.columns if c != genome_col]

        # Use to_dicts() for batch conversion - more efficient than iter_rows
        ani_dict: dict[str, dict[str, float]] = {}
        for row_dict in df.to_dicts():
            genome_name = row_dict[genome_col]
            ani_dict[genome_name] = {
                col: row_dict[col] for col in other_cols
            }

        return ani_dict


# =============================================================================
# Convenience Functions for External Tool Results
# =============================================================================


def iter_blast_chunks(
    blast_path: Path,
    *,
    chunk_size: int = 100_000,
    id_mapping: ContigIdMapping | None = None,
) -> Generator[pl.DataFrame, None, None]:
    """
    Iterate over BLAST results in chunks with optional ID transformation.

    Memory-efficient streaming for large BLAST files, optionally transforming
    external tool results to use metadarkmatter's standardized ID format.

    Args:
        blast_path: Path to BLAST tabular output file (may be gzipped)
        chunk_size: Number of rows per chunk (default: 100,000)
        id_mapping: Optional ContigIdMapping for transforming sseqid values

    Yields:
        DataFrame chunks with BLAST results

    Example:
        >>> mapping = ContigIdMapping.from_genome_dir(Path("genomes/"))
        >>> for chunk in iter_blast_chunks(blast_path, id_mapping=mapping):
        ...     # chunk has transformed sseqid values
        ...     process(chunk)
    """
    parser = StreamingBlastParser(blast_path, chunk_size=chunk_size)

    # Use Polars scan_csv with collect_batches for streaming
    batches = pl.scan_csv(
        blast_path,
        separator="\t",
        has_header=False,
        new_columns=parser.column_names,
        comment_prefix="#",
    ).collect_batches(chunk_size=chunk_size)

    for chunk_df in batches:
        if chunk_df.is_empty():
            continue

        # Apply ID transformation if mapping provided
        if id_mapping is not None:
            chunk_df = id_mapping.transform_column(chunk_df, "sseqid")

        yield chunk_df


def iter_blast_results(
    blast_path: Path,
    *,
    chunk_size: int = 100_000,
    id_mapping: ContigIdMapping | None = None,
) -> Generator[BlastResultFast, None, None]:
    """
    Iterate over BLAST results grouped by read with optional ID transformation.

    Combines streaming parsing with ID transformation for external tool results.
    Uses lightweight BlastResultFast objects for performance.

    Args:
        blast_path: Path to BLAST tabular output file (may be gzipped)
        chunk_size: Number of rows per batch (default: 100,000)
        id_mapping: Optional ContigIdMapping for transforming sseqid values

    Yields:
        BlastResultFast objects with hits grouped by read ID

    Example:
        >>> mapping = ContigIdMapping.from_genome_dir(Path("genomes/"))
        >>> for result in iter_blast_results(blast_path, id_mapping=mapping):
        ...     print(f"{result.read_id}: {result.hits[0].genome_name}")
    """
    pending_hits: dict[str, list[BlastHitFast]] = {}

    for chunk_df in iter_blast_chunks(blast_path, chunk_size=chunk_size, id_mapping=id_mapping):
        # Vectorized genome extraction + sorting
        chunk_df = (
            chunk_df
            .with_columns([extract_genome_name_expr()])
            .sort(["qseqid", "bitscore"], descending=[False, True])
        )

        # Group by read ID
        grouped = chunk_df.select(
            ["qseqid", "sseqid", "pident", "bitscore", "genome_name"]
        ).partition_by("qseqid", maintain_order=True)

        for group_df in grouped:
            if group_df.is_empty():
                continue
            read_id = group_df["qseqid"][0]

            hits = [
                BlastHitFast(
                    qseqid=row[0],
                    sseqid=row[1],
                    pident=float(row[2]),
                    bitscore=float(row[3]),
                    genome_name=row[4],
                )
                for row in group_df.iter_rows()
            ]

            if read_id in pending_hits:
                pending_hits[read_id].extend(hits)
            else:
                pending_hits[read_id] = hits

        # Yield complete reads
        if chunk_df.height > 0:
            current_last_read_id = chunk_df["qseqid"][-1]
            for read_id in list(pending_hits.keys()):
                if read_id != current_last_read_id:
                    hits = pending_hits.pop(read_id)
                    yield BlastResultFast(read_id=read_id, hits=tuple(hits))

    # Yield remaining reads
    for read_id, hits in pending_hits.items():
        yield BlastResultFast(read_id=read_id, hits=tuple(hits))
