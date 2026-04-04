"""
Data models for genome accession handling.

Provides models for storing and managing genome accessions
retrieved from GTDB and downloaded from NCBI.
"""

from __future__ import annotations

import logging
import re
from collections import Counter
from pathlib import Path
from typing import Self

import polars as pl
from pydantic import BaseModel, Field, field_validator

logger = logging.getLogger(__name__)

_ACCESSION_PATTERN = re.compile(r"^GC[AF]_\d{9}\.\d+$")
_GTDB_PREFIXES = ("d__", "p__", "c__", "o__", "f__", "g__", "s__")


class GenomeAccession(BaseModel):
    """Single genome accession with metadata from GTDB.

    Attributes:
        accession: NCBI accession (e.g., GCF_000005845.2)
        gtdb_taxonomy: Full GTDB taxonomy string
        species: Species name
        genome_size: Genome size in base pairs (if available)
        is_representative: Whether this genome is a GTDB species representative
    """

    accession: str = Field(description="NCBI accession (e.g., GCF_000005845.2)")
    gtdb_taxonomy: str = Field(description="Full GTDB taxonomy string")
    species: str = Field(description="Species name")
    genome_size: int | None = Field(default=None, description="Genome size in bp")
    is_representative: bool = Field(default=False, description="GTDB species representative")

    model_config = {"frozen": True}

    @field_validator("accession")
    @classmethod
    def validate_accession(cls, v: str) -> str:
        """Validate NCBI accession matches expected pattern."""
        if not _ACCESSION_PATTERN.match(v):
            msg = (
                f"Invalid NCBI accession '{v}': "
                "expected format GCF_XXXXXXXXX.N or GCA_XXXXXXXXX.N "
                "(e.g., GCF_000005845.2)"
            )
            raise ValueError(msg)
        return v

    @field_validator("genome_size")
    @classmethod
    def validate_genome_size(cls, v: int | None) -> int | None:
        """Validate genome size is within plausible bounds."""
        if v is not None:
            if v <= 0:
                msg = f"Genome size must be positive, got {v}"
                raise ValueError(msg)
            if v >= 50_000_000_000:
                msg = f"Genome size {v} bp exceeds 50 Gbp upper bound"
                raise ValueError(msg)
        return v

    @field_validator("gtdb_taxonomy")
    @classmethod
    def validate_gtdb_taxonomy(cls, v: str) -> str:
        """Warn if GTDB taxonomy string is missing expected rank prefixes."""
        if v:
            missing = [p for p in _GTDB_PREFIXES if p not in v]
            if missing:
                logger.warning(
                    "GTDB taxonomy '%s' missing expected prefixes: %s",
                    v[:80],
                    ", ".join(missing),
                )
        return v

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
        representative_map: Maps species name to representative accession.
            Used to populate the representative column in metadata TSV.
    """

    taxon: str = Field(description="Query taxon")
    accessions: list[GenomeAccession] = Field(default_factory=list)
    genus_counts: dict[str, int] = Field(default_factory=dict)
    representative_map: dict[str, str] = Field(
        default_factory=dict,
        description="Maps species name to representative genome accession",
    )

    @classmethod
    def from_tsv(cls, path: Path) -> Self:
        """Load accession list from TSV file.

        Args:
            path: Path to TSV file with accession data

        Returns:
            AccessionList populated from file
        """
        df = pl.read_csv(path, separator="\t")

        has_representative = "representative" in df.columns

        accessions = [
            GenomeAccession(
                accession=row["accession"],
                gtdb_taxonomy=row.get("gtdb_taxonomy", ""),
                species=row.get("species", ""),
                genome_size=row.get("genome_size"),
                is_representative=(
                    row.get("representative", "") == row["accession"]
                    if has_representative else False
                ),
            )
            for row in df.iter_rows(named=True)
        ]

        # Compute genus counts
        genus_counts = Counter(a.genus for a in accessions if a.genus)

        # Build representative map from file if present
        representative_map: dict[str, str] = {}
        if has_representative:
            for row in df.iter_rows(named=True):
                species = row.get("species", "")
                rep = row.get("representative", "")
                if species and rep and rep == row["accession"]:
                    # This row IS the representative for its species
                    representative_map[species] = rep
                elif species and rep and species not in representative_map:
                    # Store mapping but don't overwrite if we already found the actual rep
                    representative_map[species] = rep

        # Try to get taxon from header or infer from data
        taxon = ""
        if "taxon" in df.columns and len(df) > 0:
            taxon = df["taxon"][0]

        return cls(
            taxon=taxon,
            accessions=accessions,
            genus_counts=dict(genus_counts),
            representative_map=representative_map,
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

        Columns: accession, species, genus, family, representative, gtdb_taxonomy

        The representative column contains the accession of the GTDB species
        representative for each genome's species. When representative_map is
        empty, each genome maps to itself.

        Args:
            path: Output path for metadata TSV file
        """
        # Build representative column: look up each genome's species in the map
        representatives = []
        for a in self.accessions:
            if self.representative_map:
                rep = self.representative_map.get(a.species, a.accession)
            else:
                # No mapping available: each genome is its own representative
                rep = a.accession
            representatives.append(rep)

        columns = {
            "accession": [a.accession for a in self.accessions],
            "species": [a.species for a in self.accessions],
            "genus": [a.genus for a in self.accessions],
            "family": [a.family for a in self.accessions],
            "representative": representatives,
            "gtdb_taxonomy": [a.gtdb_taxonomy for a in self.accessions],
        }

        df = pl.DataFrame(columns)
        df.write_csv(path, separator="\t")
