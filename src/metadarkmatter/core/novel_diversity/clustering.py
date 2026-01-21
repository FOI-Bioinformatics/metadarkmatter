"""
Novel diversity clustering for identifying putative novel taxa.

Groups novel reads by their nearest reference genome and novelty level
to form clusters representing putative novel species or genera.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Literal

import polars as pl

from metadarkmatter.core.novel_diversity.models import (
    NovelCluster,
    NovelDiversitySummary,
)

if TYPE_CHECKING:
    from metadarkmatter.core.metadata import GenomeMetadata

logger = logging.getLogger(__name__)


class NovelDiversityAnalyzer:
    """
    Analyzes novel diversity by clustering reads into putative novel taxa.

    Clusters are formed based on:
    - Best match genome (closest reference)
    - Novelty band (5% intervals)
    - Taxonomic call (Novel Species or Novel Genus)

    This enables researchers to:
    1. Identify distinct novel organisms in their sample
    2. Prioritize clusters for further investigation
    3. Understand phylogenetic context of novel diversity

    Example:
        >>> analyzer = NovelDiversityAnalyzer(
        ...     classifications=df,
        ...     metadata=genome_metadata,
        ... )
        >>> clusters = analyzer.cluster_novel_reads()
        >>> high_priority = [c for c in clusters if c.confidence == "High"]
    """

    def __init__(
        self,
        classifications: pl.DataFrame,
        metadata: GenomeMetadata | None = None,
        novelty_band_size: float = 5.0,
        min_cluster_size: int = 3,
    ) -> None:
        """
        Initialize the analyzer.

        Args:
            classifications: DataFrame with classification results.
                Required columns: read_id, best_match_genome, novelty_index,
                placement_uncertainty, taxonomic_call
                Optional: discovery_score (for enhanced prioritization)
            metadata: Genome metadata for species/genus/family lookups
            novelty_band_size: Size of novelty bands in percent (default 5.0)
            min_cluster_size: Minimum reads to form a cluster (default 3)
        """
        self._df = classifications
        self._metadata = metadata
        self._novelty_band_size = novelty_band_size
        self._min_cluster_size = min_cluster_size
        self._clusters: list[NovelCluster] | None = None
        self._summary: NovelDiversitySummary | None = None

        # Validate required columns
        required = {"read_id", "best_match_genome", "novelty_index",
                    "placement_uncertainty", "taxonomic_call"}
        missing = required - set(classifications.columns)
        if missing:
            msg = f"Classification DataFrame missing required columns: {missing}"
            raise ValueError(msg)

    def cluster_novel_reads(self) -> list[NovelCluster]:
        """
        Cluster novel reads into putative novel taxa.

        Returns:
            List of NovelCluster objects representing putative novel taxa

        Algorithm:
            1. Filter to Novel Species and Novel Genus reads only
            2. Compute novelty band for each read
            3. Group by (best_match_genome, novelty_band, taxonomic_call)
            4. Aggregate statistics for each group
            5. Filter by minimum cluster size
            6. Generate cluster IDs and suggested names
            7. Assign confidence ratings
            8. Infer phylogenetic context
        """
        if self._clusters is not None:
            return self._clusters

        # Step 1: Filter to novel reads only
        novel_df = self._df.filter(
            pl.col("taxonomic_call").is_in(["Novel Species", "Novel Genus"])
        )

        if len(novel_df) == 0:
            logger.info("No novel reads found for clustering")
            self._clusters = []
            return self._clusters

        # Step 2: Compute novelty band
        novel_df = novel_df.with_columns(
            (pl.col("novelty_index") / self._novelty_band_size)
            .floor()
            .cast(pl.Int32)
            .mul(int(self._novelty_band_size))
            .alias("novelty_band")
        )

        # Step 2b: Compute effective uncertainty
        # For single-hit reads (uncertainty=0), infer from novelty index
        # Higher novelty = more uncertainty about exact phylogenetic placement
        has_ambiguous_hits = "num_ambiguous_hits" in novel_df.columns

        if has_ambiguous_hits:
            # Use novelty-based inference for single hits (where uncertainty is 0)
            # Inferred uncertainty = novelty_index * 0.4 (capped at 15%)
            # This reflects that higher divergence = less certain placement
            novel_df = novel_df.with_columns(
                pl.when(pl.col("num_ambiguous_hits") <= 1)
                .then(
                    pl.min_horizontal(
                        pl.col("novelty_index") * 0.4,
                        pl.lit(15.0)
                    )
                )
                .otherwise(pl.col("placement_uncertainty"))
                .alias("effective_uncertainty")
            )
        else:
            # If num_ambiguous_hits not available, infer for zero-uncertainty reads
            novel_df = novel_df.with_columns(
                pl.when(pl.col("placement_uncertainty") < 0.1)
                .then(
                    pl.min_horizontal(
                        pl.col("novelty_index") * 0.4,
                        pl.lit(15.0)
                    )
                )
                .otherwise(pl.col("placement_uncertainty"))
                .alias("effective_uncertainty")
            )

        # Step 3: Group by clustering key
        has_discovery = "discovery_score" in novel_df.columns

        agg_exprs = [
            pl.len().alias("read_count"),
            pl.col("novelty_index").mean().alias("mean_novelty"),
            pl.col("novelty_index").min().alias("novelty_min"),
            pl.col("novelty_index").max().alias("novelty_max"),
            pl.col("effective_uncertainty").mean().alias("mean_uncertainty"),
            pl.col("read_id").alias("read_ids"),
        ]

        if has_discovery:
            agg_exprs.append(
                pl.col("discovery_score").mean().alias("mean_discovery")
            )

        grouped = (
            novel_df
            .group_by(["best_match_genome", "novelty_band", "taxonomic_call"])
            .agg(agg_exprs)
        )

        # Step 4: Filter by minimum cluster size
        grouped = grouped.filter(pl.col("read_count") >= self._min_cluster_size)

        if len(grouped) == 0:
            logger.info(
                f"No clusters with >= {self._min_cluster_size} reads found"
            )
            self._clusters = []
            return self._clusters

        # Sort by read count descending
        grouped = grouped.sort("read_count", descending=True)

        # Step 5-8: Build NovelCluster objects
        clusters: list[NovelCluster] = []
        nsp_counter = 0
        ngn_counter = 0

        for row in grouped.iter_rows(named=True):
            genome = row["best_match_genome"]
            tax_call = row["taxonomic_call"]
            novelty_band = row["novelty_band"]
            read_count = row["read_count"]
            mean_novelty = row["mean_novelty"]
            novelty_min = row["novelty_min"]
            novelty_max = row["novelty_max"]
            mean_uncertainty = row["mean_uncertainty"]
            read_ids = row["read_ids"]
            mean_discovery = row.get("mean_discovery")

            # Generate cluster ID
            if tax_call == "Novel Species":
                nsp_counter += 1
                cluster_id = f"NSP_{nsp_counter:03d}"
            else:
                ngn_counter += 1
                cluster_id = f"NGN_{ngn_counter:03d}"

            # Get taxonomy from metadata
            species, genus, family = self._get_taxonomy(genome)

            # Generate suggested name
            suggested_name = self._generate_suggested_name(
                tax_call, genus, family, cluster_id
            )

            # Assign confidence rating
            confidence = self._assign_confidence(
                read_count, mean_uncertainty, mean_discovery
            )

            # Infer phylogenetic context
            phylo_context = self._infer_phylogenetic_context(
                tax_call, genus, family
            )

            cluster = NovelCluster(
                cluster_id=cluster_id,
                taxonomic_call=tax_call,
                nearest_genome=genome,
                nearest_species=species,
                nearest_genus=genus,
                nearest_family=family,
                novelty_band=novelty_band,
                read_count=read_count,
                mean_novelty_index=mean_novelty,
                novelty_min=novelty_min,
                novelty_max=novelty_max,
                mean_placement_uncertainty=mean_uncertainty,
                mean_discovery_score=mean_discovery,
                suggested_name=suggested_name,
                confidence=confidence,
                phylogenetic_context=phylo_context,
                read_ids=read_ids,
            )
            clusters.append(cluster)

        logger.info(f"Created {len(clusters)} novel diversity clusters")
        self._clusters = clusters
        return clusters

    def get_summary(self) -> NovelDiversitySummary:
        """
        Get summary statistics for novel diversity.

        Returns:
            NovelDiversitySummary with aggregate statistics
        """
        if self._summary is not None:
            return self._summary

        clusters = self.cluster_novel_reads()

        if not clusters:
            self._summary = NovelDiversitySummary()
            return self._summary

        # Count reads by type
        novel_df = self._df.filter(
            pl.col("taxonomic_call").is_in(["Novel Species", "Novel Genus"])
        )
        nsp_reads = len(novel_df.filter(pl.col("taxonomic_call") == "Novel Species"))
        ngn_reads = len(novel_df.filter(pl.col("taxonomic_call") == "Novel Genus"))

        # Count clusters by type and confidence
        nsp_clusters = sum(1 for c in clusters if c.taxonomic_call == "Novel Species")
        ngn_clusters = sum(1 for c in clusters if c.taxonomic_call == "Novel Genus")
        high_conf = sum(1 for c in clusters if c.confidence == "High")
        med_conf = sum(1 for c in clusters if c.confidence == "Medium")
        low_conf = sum(1 for c in clusters if c.confidence == "Low")

        # Unique genera/families with novel diversity
        genera = set(c.nearest_genus for c in clusters if c.taxonomic_call == "Novel Species")
        families = set(c.nearest_family for c in clusters if c.taxonomic_call == "Novel Genus")

        # Cluster size statistics
        sizes = [c.read_count for c in clusters]
        mean_size = sum(sizes) / len(sizes) if sizes else 0.0
        max_size = max(sizes) if sizes else 0

        self._summary = NovelDiversitySummary(
            total_novel_reads=nsp_reads + ngn_reads,
            novel_species_reads=nsp_reads,
            novel_genus_reads=ngn_reads,
            total_clusters=len(clusters),
            novel_species_clusters=nsp_clusters,
            novel_genus_clusters=ngn_clusters,
            high_confidence_clusters=high_conf,
            medium_confidence_clusters=med_conf,
            low_confidence_clusters=low_conf,
            genera_with_novel_species=len(genera),
            families_with_novel_genera=len(families),
            mean_cluster_size=round(mean_size, 1),
            largest_cluster_size=max_size,
        )
        return self._summary

    def to_dataframe(self) -> pl.DataFrame:
        """
        Export clusters as a Polars DataFrame.

        Returns:
            DataFrame with cluster summary information
        """
        clusters = self.cluster_novel_reads()

        if not clusters:
            return pl.DataFrame()

        data = [c.to_summary_dict() for c in clusters]
        return pl.DataFrame(data)

    def to_dict(self) -> dict:
        """
        Export as JSON-serializable dictionary.

        Returns:
            Dictionary with 'summary' and 'clusters' keys
        """
        clusters = self.cluster_novel_reads()
        summary = self.get_summary()

        return {
            "summary": summary.model_dump(),
            "clusters": [c.model_dump() for c in clusters],
        }

    def _get_taxonomy(self, accession: str) -> tuple[str, str, str]:
        """
        Get species, genus, and family for a genome accession.

        Args:
            accession: Genome accession

        Returns:
            Tuple of (species, genus, family)
        """
        if self._metadata is None:
            return ("Unknown species", "Unknown", "Unknown family")

        species = self._metadata.get_species(accession) or "Unknown species"
        genus = self._metadata.get_genus(accession) or "Unknown"
        family = self._metadata.get_family(accession) or "Unknown family"

        return (species, genus, family)

    def _generate_suggested_name(
        self,
        taxonomic_call: str,
        genus: str,
        family: str,
        cluster_id: str,
    ) -> str:
        """
        Generate a GTDB-style provisional name for a novel taxon.

        Args:
            taxonomic_call: "Novel Species" or "Novel Genus"
            genus: Nearest genus name
            family: Nearest family name
            cluster_id: Unique cluster identifier

        Returns:
            Suggested provisional name
        """
        # Extract numeric suffix from cluster_id (e.g., "001" from "NSP_001")
        suffix = cluster_id.split("_")[1]

        if taxonomic_call == "Novel Species":
            return f"{genus} sp. MDM-{suffix}"
        else:  # Novel Genus
            if family and family != "Unknown family":
                return f"{family} gen. nov. MDM-{suffix}"
            else:
                return f"{genus}-related gen. nov. MDM-{suffix}"

    def _assign_confidence(
        self,
        read_count: int,
        mean_uncertainty: float,
        mean_discovery: float | None,
    ) -> Literal["High", "Medium", "Low"]:
        """
        Assign a confidence rating to a cluster.

        Criteria:
            High: read_count >= 10, mean_uncertainty < 5%, discovery_score >= 75
            Medium: read_count >= 5, mean_uncertainty < 10%, discovery_score >= 50
            Low: All other clusters

        Args:
            read_count: Number of reads in cluster
            mean_uncertainty: Mean placement uncertainty
            mean_discovery: Mean discovery score (optional)

        Returns:
            Confidence rating string
        """
        # Use discovery score if available, otherwise use simplified criteria
        if mean_discovery is not None:
            if read_count >= 10 and mean_uncertainty < 5 and mean_discovery >= 75:
                return "High"
            elif read_count >= 5 and mean_uncertainty < 10 and mean_discovery >= 50:
                return "Medium"
        else:
            # Simplified criteria without discovery score
            if read_count >= 10 and mean_uncertainty < 5:
                return "High"
            elif read_count >= 5 and mean_uncertainty < 10:
                return "Medium"

        return "Low"

    def _infer_phylogenetic_context(
        self,
        taxonomic_call: str,
        genus: str,
        family: str,
    ) -> str:
        """
        Infer phylogenetic placement description for a cluster.

        Args:
            taxonomic_call: "Novel Species" or "Novel Genus"
            genus: Nearest genus name
            family: Nearest family name

        Returns:
            Phylogenetic context description
        """
        if taxonomic_call == "Novel Species":
            if genus and genus != "Unknown":
                return f"Novel species within {genus}"
            else:
                return "Novel species (genus unknown)"
        else:  # Novel Genus
            if family and family != "Unknown family":
                return f"Novel genus within {family}"
            elif genus and genus != "Unknown":
                return f"Novel genus sister to {genus}"
            else:
                return "Novel genus (family unknown)"
