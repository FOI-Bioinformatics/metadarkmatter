"""
Data models for genome accession handling.

Provides models for storing and managing genome accessions
retrieved from GTDB and downloaded from NCBI.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Self

import polars as pl
from pydantic import BaseModel, Field


class GenomeAccession(BaseModel):
    """Single genome accession with metadata from GTDB.

    Attributes:
        accession: NCBI accession (e.g., GCF_000005845.2)
        gtdb_taxonomy: Full GTDB taxonomy string
        species: Species name
        genome_size: Genome size in base pairs (if available)
    """

    accession: str = Field(description="NCBI accession (e.g., GCF_000005845.2)")
    gtdb_taxonomy: str = Field(description="Full GTDB taxonomy string")
    species: str = Field(description="Species name")
    genome_size: int | None = Field(default=None, description="Genome size in bp")

    model_config = {"frozen": True}

    @property
    def genus(self) -> str:
        """Extract genus from GTDB taxonomy string."""
        for part in self.gtdb_taxonomy.split(";"):
            part = part.strip()
            if part.startswith("g__"):
                return part[3:]  # Remove "g__" prefix
        return ""

    @property
    def family(self) -> str:
        """Extract family from GTDB taxonomy string."""
        for part in self.gtdb_taxonomy.split(";"):
            part = part.strip()
            if part.startswith("f__"):
                return part[3:]  # Remove "f__" prefix
        return ""


class AccessionList(BaseModel):
    """Collection of genome accessions with summary statistics.

    Provides methods for reading/writing TSV files and computing
    summary statistics about the accessions.

    Attributes:
        taxon: Query taxon string (e.g., f__Enterobacteriaceae)
        accessions: List of genome accessions
        genus_counts: Count of genomes per genus
    """

    taxon: str = Field(description="Query taxon")
    accessions: list[GenomeAccession] = Field(default_factory=list)
    genus_counts: dict[str, int] = Field(default_factory=dict)

    @classmethod
    def from_tsv(cls, path: Path) -> Self:
        """Load accession list from TSV file.

        Args:
            path: Path to TSV file with accession data

        Returns:
            AccessionList populated from file
        """
        df = pl.read_csv(path, separator="\t")

        accessions = [
            GenomeAccession(
                accession=row["accession"],
                gtdb_taxonomy=row.get("gtdb_taxonomy", ""),
                species=row.get("species", ""),
                genome_size=row.get("genome_size"),
            )
            for row in df.iter_rows(named=True)
        ]

        # Compute genus counts
        genus_counts = Counter(a.genus for a in accessions if a.genus)

        # Try to get taxon from header or infer from data
        taxon = ""
        if "taxon" in df.columns and len(df) > 0:
            taxon = df["taxon"][0]

        return cls(
            taxon=taxon,
            accessions=accessions,
            genus_counts=dict(genus_counts),
        )

    def to_tsv(self, path: Path, *, include_metadata: bool = False) -> None:
        """Save accession list to TSV file.

        Args:
            path: Output path for TSV file
            include_metadata: Include additional columns (genome_size, etc.)
        """
        columns = {
            "accession": [a.accession for a in self.accessions],
            "gtdb_taxonomy": [a.gtdb_taxonomy for a in self.accessions],
            "species": [a.species for a in self.accessions],
        }

        if include_metadata:
            columns["genome_size"] = [a.genome_size for a in self.accessions]
            columns["genus"] = [a.genus for a in self.accessions]
            columns["family"] = [a.family for a in self.accessions]

        df = pl.DataFrame(columns)
        df.write_csv(path, separator="\t")

    @property
    def total_count(self) -> int:
        """Total number of accessions."""
        return len(self.accessions)

    @property
    def species_count(self) -> int:
        """Number of unique species."""
        return len({a.species for a in self.accessions if a.species})

    @property
    def genus_count(self) -> int:
        """Number of unique genera."""
        return len(self.genus_counts)

    def get_accession_strings(self) -> list[str]:
        """Get list of accession strings for NCBI download.

        Returns:
            List of accession strings
        """
        return [a.accession for a in self.accessions]

    def to_metadata_tsv(self, path: Path) -> None:
        """Save genome metadata for use throughout the pipeline.

        Creates a standardized metadata file with taxonomy information
        that can be joined with classification results for species-level
        aggregation.

        Columns: accession, species, genus, family, gtdb_taxonomy

        Args:
            path: Output path for metadata TSV file
        """
        columns = {
            "accession": [a.accession for a in self.accessions],
            "species": [a.species for a in self.accessions],
            "genus": [a.genus for a in self.accessions],
            "family": [a.family for a in self.accessions],
            "gtdb_taxonomy": [a.gtdb_taxonomy for a in self.accessions],
        }

        df = pl.DataFrame(columns)
        df.write_csv(path, separator="\t")
