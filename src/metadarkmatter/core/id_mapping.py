"""
ID mapping utilities for external tool results.

Provides functionality to transform contig IDs in external BLAST, Bowtie2,
and Kraken2 results to the standardized format expected by metadarkmatter.
"""

from __future__ import annotations

import gzip
import logging
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

from metadarkmatter.core.genome_utils import extract_accession_from_filename

logger = logging.getLogger(__name__)


@dataclass
class ContigIdMapping:
    """Mapping from original contig IDs to genome accessions.

    This class enables transformation of external tool results to use
    metadarkmatter's standardized ID format: {accession}|{original_contig_id}

    Example:
        >>> mapping = ContigIdMapping.from_genome_dir(Path("genomes/"))
        >>> mapping.transform_sseqid("NZ_CP007557.1")
        'GCF_000195955.2|NZ_CP007557.1'
    """

    contig_to_accession: dict[str, str] = field(default_factory=dict)

    @classmethod
    def from_genome_dir(
        cls,
        genome_dir: Path,
        pattern: str = "*.fna",
    ) -> ContigIdMapping:
        """Generate mapping from genome FASTA directory.

        Parses genome FASTA files to build a mapping from original contig IDs
        to their parent genome accession. Accessions are extracted from filenames.

        Args:
            genome_dir: Directory containing genome FASTA files
            pattern: Glob pattern for genome files

        Returns:
            ContigIdMapping instance

        Raises:
            FileNotFoundError: If genome_dir does not exist
            ValueError: If no genome files found
        """
        if not genome_dir.exists():
            raise FileNotFoundError(f"Genome directory not found: {genome_dir}")

        genome_files = sorted(genome_dir.glob(pattern))

        # Try alternative patterns
        if not genome_files:
            for alt_pattern in [
                "*.fa",
                "*.fasta",
                "*.fna.gz",
                "*.fa.gz",
                "*.fasta.gz",
            ]:
                genome_files = sorted(genome_dir.glob(alt_pattern))
                if genome_files:
                    break

        if not genome_files:
            raise ValueError(
                f"No genome files found in {genome_dir} with pattern {pattern}"
            )

        mapping: dict[str, str] = {}
        duplicates: list[str] = []

        for genome_file in genome_files:
            accession = extract_accession_from_filename(genome_file.name)

            if not accession:
                logger.warning(
                    "Could not extract accession from filename: %s", genome_file.name
                )
                continue

            # Handle gzipped files
            open_func = gzip.open if genome_file.suffix == ".gz" else Path.open

            try:
                with open_func(genome_file, "rt") as f:  # type: ignore[operator]
                    for line in f:
                        if line.startswith(">"):
                            # Extract contig ID (first word after >)
                            header = line[1:].strip()
                            contig_id = header.split()[0]

                            if contig_id in mapping:
                                duplicates.append(contig_id)
                            mapping[contig_id] = accession
            except Exception as e:
                logger.warning("Error reading %s: %s", genome_file, e)
                continue

        if duplicates:
            logger.warning(
                "Found %d duplicate contig IDs (last genome wins): %s...",
                len(duplicates),
                duplicates[:3],
            )

        logger.info(
            "Built ID mapping with %d contigs from %d genomes",
            len(mapping),
            len(genome_files),
        )

        return cls(contig_to_accession=mapping)

    @classmethod
    def from_tsv(cls, path: Path) -> ContigIdMapping:
        """Load mapping from TSV file.

        Expected format (tab-separated):
            original_contig_id	genome_accession
            NZ_CP007557.1	GCF_000195955.2
            NC_006570.2	GCF_000008985.1

        Args:
            path: Path to TSV file

        Returns:
            ContigIdMapping instance
        """
        df = pl.read_csv(path, separator="\t")

        # Validate columns
        required = {"original_contig_id", "genome_accession"}
        if not required.issubset(df.columns):
            raise ValueError(
                f"TSV must have columns: {required}, found: {df.columns}"
            )

        mapping = dict(
            zip(
                df["original_contig_id"].to_list(),
                df["genome_accession"].to_list(),
                strict=True,
            )
        )

        logger.info("Loaded ID mapping with %d entries from %s", len(mapping), path)

        return cls(contig_to_accession=mapping)

    def to_tsv(self, path: Path) -> None:
        """Save mapping to TSV file.

        Args:
            path: Output path for TSV file
        """
        if not self.contig_to_accession:
            raise ValueError("Cannot save empty mapping")

        df = pl.DataFrame(
            {
                "original_contig_id": list(self.contig_to_accession.keys()),
                "genome_accession": list(self.contig_to_accession.values()),
            }
        )
        df.write_csv(path, separator="\t")
        logger.info("Saved ID mapping with %d entries to %s", len(df), path)

    def transform_sseqid(self, original_id: str) -> str:
        """Transform single ID to standardized format.

        Args:
            original_id: Original contig ID (e.g., NZ_CP007557.1)

        Returns:
            Standardized ID (e.g., GCF_000195955.2|NZ_CP007557.1)
            or original ID if not found in mapping
        """
        # Already in standardized format
        if "|" in original_id:
            return original_id

        accession = self.contig_to_accession.get(original_id)
        if accession:
            return f"{accession}|{original_id}"

        return original_id

    def transform_column(
        self,
        df: pl.DataFrame,
        col: str = "sseqid",
    ) -> pl.DataFrame:
        """Transform ID column in DataFrame (vectorized).

        Efficiently transforms all IDs in a column using Polars join operations.
        IDs not found in the mapping are passed through unchanged.

        Args:
            df: DataFrame with ID column to transform
            col: Name of column containing IDs

        Returns:
            DataFrame with transformed IDs in the specified column
        """
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in DataFrame")

        if not self.contig_to_accession:
            logger.warning("Empty ID mapping, returning DataFrame unchanged")
            return df

        # Build mapping DataFrame for join
        mapping_df = pl.DataFrame(
            {
                "_original_id": list(self.contig_to_accession.keys()),
                "_accession": list(self.contig_to_accession.values()),
            }
        )

        # Join and transform
        result = (
            df.join(mapping_df, left_on=col, right_on="_original_id", how="left")
            .with_columns(
                pl.when(
                    pl.col("_accession").is_not_null()
                    & ~pl.col(col).str.contains(r"\|", literal=False)
                )
                .then(pl.col("_accession") + "|" + pl.col(col))
                .otherwise(pl.col(col))
                .alias(col)
            )
            .drop("_accession")
        )

        return result

    def __len__(self) -> int:
        """Return number of entries in mapping."""
        return len(self.contig_to_accession)

    def __contains__(self, item: str) -> bool:
        """Check if contig ID is in mapping."""
        return item in self.contig_to_accession
