"""
Genome metadata handling for species-level classification.

Provides classes for loading and joining genome metadata with classification
results, enabling species-level aggregation in reports.

Includes GTDB-style naming for novel taxa:
- Novel species: "{Genus} sp. MDM{id}" (e.g., "Francisella sp. MDM001")
- Novel genus: "{Family} gen. nov. MDM{id}" (e.g., "Francisellaceae gen. nov. MDM001")
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import polars as pl


class GenomeMetadata:
    """Genome metadata for species-level aggregation.

    Loads metadata from TSV file and provides efficient lookups and joining
    with classification results. The metadata file is created during the
    'download genomes list' command with columns:
        accession, species, genus, family, gtdb_taxonomy

    Example:
        >>> metadata = GenomeMetadata.from_file(Path("genome_metadata.tsv"))
        >>> classifications = pl.read_csv("classifications.csv")
        >>> enriched = metadata.join_classifications(classifications)
        >>> # enriched now has 'species' and 'genus' columns
    """

    def __init__(self, metadata_df: pl.DataFrame) -> None:
        """Initialize from DataFrame.

        Args:
            metadata_df: DataFrame with columns: accession, species, genus, family, gtdb_taxonomy
                Optional column: representative (species representative accession)
        """
        self._df = metadata_df

        # Build lookup dictionaries for fast access
        self._accession_to_species: dict[str, str] = dict(
            zip(
                metadata_df["accession"].to_list(),
                metadata_df["species"].to_list(),
                strict=False,
            )
        )
        self._accession_to_genus: dict[str, str] = dict(
            zip(
                metadata_df["accession"].to_list(),
                metadata_df["genus"].to_list(),
                strict=False,
            )
        )
        # Family lookup (optional column)
        if "family" in metadata_df.columns:
            self._accession_to_family: dict[str, str] = dict(
                zip(
                    metadata_df["accession"].to_list(),
                    metadata_df["family"].to_list(),
                    strict=False,
                )
            )
        else:
            self._accession_to_family = {}

        # Representative lookup (optional column)
        if "representative" in metadata_df.columns:
            self._accession_to_representative: dict[str, str] = dict(
                zip(
                    metadata_df["accession"].to_list(),
                    metadata_df["representative"].to_list(),
                    strict=False,
                )
            )
        else:
            self._accession_to_representative = {}

    @classmethod
    def from_file(cls, path: Path) -> GenomeMetadata:
        """Load metadata from TSV file.

        Args:
            path: Path to genome_metadata.tsv file

        Returns:
            GenomeMetadata instance

        Raises:
            FileNotFoundError: If file does not exist
            ValueError: If file is missing required columns
        """
        if not path.exists():
            msg = f"Metadata file not found: {path}"
            raise FileNotFoundError(msg)

        df = pl.read_csv(path, separator="\t")

        # Validate required columns
        required = {"accession", "species", "genus"}
        missing = required - set(df.columns)
        if missing:
            msg = f"Metadata file missing required columns: {missing}"
            raise ValueError(msg)

        return cls(df)

    @property
    def dataframe(self) -> pl.DataFrame:
        """Get underlying DataFrame for custom operations."""
        return self._df

    @property
    def genome_count(self) -> int:
        """Number of genomes in metadata."""
        return len(self._df)

    @property
    def species_count(self) -> int:
        """Number of unique species in metadata."""
        return self._df["species"].n_unique()

    @property
    def genus_count(self) -> int:
        """Number of unique genera in metadata."""
        return self._df["genus"].n_unique()

    def get_species(self, accession: str) -> str | None:
        """Get species name for a genome accession.

        Args:
            accession: Genome accession (e.g., GCF_000195955.2)

        Returns:
            Species name or None if not found
        """
        return self._accession_to_species.get(accession)

    def get_genus(self, accession: str) -> str | None:
        """Get genus name for a genome accession.

        Args:
            accession: Genome accession (e.g., GCF_000195955.2)

        Returns:
            Genus name or None if not found
        """
        return self._accession_to_genus.get(accession)

    def get_family(self, accession: str) -> str | None:
        """Get family name for a genome accession.

        Args:
            accession: Genome accession (e.g., GCF_000195955.2)

        Returns:
            Family name or None if not found
        """
        return self._accession_to_family.get(accession)

    def infer_target_family(self) -> str | None:
        """Infer the target family from genome metadata.

        Returns the most common family among genomes in this metadata set.
        Used as fallback when --target-family is not explicitly provided.

        Returns:
            Most common family name, or None if family column is absent or empty.
        """
        if "family" not in self._df.columns:
            return None
        if self._df.is_empty():
            return None
        family_counts = (
            self._df.group_by("family")
            .len()
            .sort("len", descending=True)
        )
        if family_counts.is_empty():
            return None
        return family_counts["family"][0]

    def get_representative(self, accession: str) -> str:
        """Get the species representative accession for a genome.

        If no representative column exists in metadata, returns the
        accession itself (identity mapping for backwards compatibility).

        Args:
            accession: Genome accession (e.g., GCF_000195955.2)

        Returns:
            Representative genome accession
        """
        if not self._accession_to_representative:
            return accession
        return self._accession_to_representative.get(accession, accession)

    def build_representative_mapping(self) -> dict[str, str]:
        """Build a complete accession-to-representative mapping dict.

        Returns a dictionary mapping every genome accession to its species
        representative accession. If no representative column exists, returns
        an identity mapping (each genome maps to itself).

        Returns:
            Dictionary of {accession: representative_accession}
        """
        if not self._accession_to_representative:
            # Identity mapping: each genome is its own representative
            return {acc: acc for acc in self._accession_to_species}
        return dict(self._accession_to_representative)

    @property
    def has_representatives(self) -> bool:
        """True if representative column exists and has non-self mappings."""
        if not self._accession_to_representative:
            return False
        return any(
            acc != rep
            for acc, rep in self._accession_to_representative.items()
        )

    @property
    def representative_count(self) -> int:
        """Number of unique representative genomes."""
        if not self._accession_to_representative:
            return self.genome_count
        return len(set(self._accession_to_representative.values()))

    @staticmethod
    def generate_novel_id(read_id: str, accession: str) -> str:
        """Generate a deterministic short ID for a novel taxon.

        Creates a reproducible 3-character alphanumeric ID based on
        the read ID and closest reference genome. This ensures the same
        read always gets the same novel taxon ID.

        Args:
            read_id: The read identifier
            accession: Closest reference genome accession

        Returns:
            3-character alphanumeric ID (e.g., "A7X")
        """
        # Create deterministic hash from read_id + accession
        hash_input = f"{read_id}:{accession}"
        hash_bytes = hashlib.md5(hash_input.encode()).digest()  # noqa: S324

        # Convert first 2 bytes to alphanumeric (A-Z, 0-9)
        chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
        id_chars = [
            chars[hash_bytes[0] % len(chars)],
            chars[hash_bytes[1] % len(chars)],
            chars[hash_bytes[2] % len(chars)],
        ]
        return "".join(id_chars)

    def suggest_novel_species_name(
        self,
        accession: str,
        read_id: str | None = None,
    ) -> str:
        """Suggest a GTDB-style name for a novel species.

        Format: "{Genus} sp. MDM{id}"

        Args:
            accession: Closest reference genome accession
            read_id: Optional read ID for deterministic ID generation

        Returns:
            Suggested novel species name (e.g., "Francisella sp. MDM-A7X")
        """
        genus = self.get_genus(accession) or "Unknown"

        if read_id:
            novel_id = self.generate_novel_id(read_id, accession)
            return f"{genus} sp. MDM-{novel_id}"
        else:
            return f"{genus} sp. nov."

    def suggest_novel_genus_name(
        self,
        accession: str,
        read_id: str | None = None,
    ) -> str:
        """Suggest a GTDB-style name for a novel genus.

        Format: "{Family} gen. nov. MDM{id}"

        Args:
            accession: Closest reference genome accession
            read_id: Optional read ID for deterministic ID generation

        Returns:
            Suggested novel genus name (e.g., "Francisellaceae gen. nov. MDM-A7X")
        """
        family = self.get_family(accession)
        if not family:
            # Fall back to genus-based naming
            genus = self.get_genus(accession) or "Unknown"
            family = f"{genus}-related"

        if read_id:
            novel_id = self.generate_novel_id(read_id, accession)
            return f"{family} gen. nov. MDM-{novel_id}"
        else:
            return f"{family} gen. nov."

    def join_classifications(self, df: pl.DataFrame) -> pl.DataFrame:
        """Join metadata to classification DataFrame.

        Adds species and genus columns based on the best_match_genome column.
        Unknown genomes are filled with "Unknown species" / "Unknown genus".

        Args:
            df: Classification DataFrame with 'best_match_genome' column

        Returns:
            DataFrame with added 'species' and 'genus' columns
        """
        if "best_match_genome" not in df.columns:
            msg = "Classification DataFrame must have 'best_match_genome' column"
            raise ValueError(msg)

        # Select columns to join (include family if available)
        join_cols = ["accession", "species", "genus"]
        if "family" in self._df.columns:
            join_cols.append("family")

        result = df.join(
            self._df.select(join_cols),
            left_on="best_match_genome",
            right_on="accession",
            how="left",
        ).with_columns(
            [
                pl.col("species").fill_null("Unknown species"),
                pl.col("genus").fill_null("Unknown genus"),
            ]
        )

        if "family" in self._df.columns:
            result = result.with_columns(
                pl.col("family").fill_null("Unknown family")
            )

        return result

    def join_classifications_with_novel_names(
        self,
        df: pl.DataFrame,
    ) -> pl.DataFrame:
        """Join metadata and add suggested names for novel taxa.

        Extends join_classifications by adding a 'suggested_name' column
        that provides GTDB-style names for novel species and genera.

        For known species: uses the actual species name
        For novel species: "{Genus} sp. MDM-{id}"
        For novel genus: "{Family} gen. nov. MDM-{id}"
        For other: uses closest reference species

        Args:
            df: Classification DataFrame with 'best_match_genome', 'taxonomic_call',
                and 'read_id' columns

        Returns:
            DataFrame with added 'species', 'genus', 'family', and 'suggested_name' columns
        """
        required = {"best_match_genome", "taxonomic_call", "read_id"}
        missing = required - set(df.columns)
        if missing:
            msg = f"Classification DataFrame missing required columns: {missing}"
            raise ValueError(msg)

        # First do standard join
        result = self.join_classifications(df)

        # Generate deterministic novel IDs using struct to pass both columns
        result = result.with_columns(
            pl.struct(["read_id", "best_match_genome"])
            .map_elements(
                lambda x: self.generate_novel_id(x["read_id"], x["best_match_genome"]),
                return_dtype=pl.Utf8,
            )
            .alias("_novel_id")
        )

        # Add family column if not present
        if "family" not in result.columns:
            result = result.with_columns(pl.lit("Unknown family").alias("family"))

        # Generate suggested names based on taxonomic call
        result = result.with_columns(
            pl.when(pl.col("taxonomic_call") == "Novel Species")
            .then(pl.col("genus") + pl.lit(" sp. MDM-") + pl.col("_novel_id"))
            .when(pl.col("taxonomic_call") == "Novel Genus")
            .then(pl.col("family") + pl.lit(" gen. nov. MDM-") + pl.col("_novel_id"))
            .when(pl.col("taxonomic_call") == "Known Species")
            .then(pl.col("species"))
            .otherwise(pl.col("species") + pl.lit(" (") + pl.col("taxonomic_call") + pl.lit(")"))
            .alias("suggested_name")
        ).drop("_novel_id")

        return result

    def aggregate_by_species(self, df: pl.DataFrame) -> pl.DataFrame:
        """Aggregate classification counts by species.

        Groups reads by species and calculates summary statistics.
        Requires 'species' column (use join_classifications first).

        Args:
            df: Classification DataFrame with 'species' column

        Returns:
            DataFrame with species-level statistics:
                - species: Species name
                - read_count: Number of reads
                - mean_novelty: Mean novelty index
                - mean_identity: Mean percent identity
                - mean_uncertainty: Mean placement uncertainty
                - genome_count: Number of genomes for this species
        """
        if "species" not in df.columns:
            msg = "DataFrame must have 'species' column. Use join_classifications first."
            raise ValueError(msg)

        return (
            df.group_by("species")
            .agg(
                [
                    pl.len().alias("read_count"),
                    pl.col("novelty_index").mean().alias("mean_novelty"),
                    pl.col("top_hit_identity").mean().alias("mean_identity"),
                    pl.col("placement_uncertainty").mean().alias("mean_uncertainty"),
                    pl.col("best_match_genome").n_unique().alias("genome_count"),
                ]
            )
            .sort("read_count", descending=True)
        )

    def aggregate_by_genus(self, df: pl.DataFrame) -> pl.DataFrame:
        """Aggregate classification counts by genus.

        Groups reads by genus and calculates summary statistics.
        Requires 'genus' column (use join_classifications first).

        Args:
            df: Classification DataFrame with 'genus' column

        Returns:
            DataFrame with genus-level statistics
        """
        if "genus" not in df.columns:
            msg = "DataFrame must have 'genus' column. Use join_classifications first."
            raise ValueError(msg)

        return (
            df.group_by("genus")
            .agg(
                [
                    pl.len().alias("read_count"),
                    pl.col("novelty_index").mean().alias("mean_novelty"),
                    pl.col("top_hit_identity").mean().alias("mean_identity"),
                    pl.col("placement_uncertainty").mean().alias("mean_uncertainty"),
                    pl.col("best_match_genome").n_unique().alias("genome_count"),
                    pl.col("species").n_unique().alias("species_count"),
                ]
            )
            .sort("read_count", descending=True)
        )
