"""
Novel diversity clustering for identifying putative novel taxa.

Groups novel reads by their nearest reference genome and novelty level
to form clusters representing putative novel species or genera.  When an
ANI matrix is provided, genomes within a shared ANI neighborhood are
merged before clustering, reducing inflation caused by reads scattering
across equidistant references.
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
    from metadarkmatter.core.classification.ani_matrix import ANIMatrix
    from metadarkmatter.core.metadata import GenomeMetadata

logger = logging.getLogger(__name__)


class NovelDiversityAnalyzer:
    """
    Analyzes novel diversity by clustering reads into putative novel taxa.

    Clusters are formed based on:
    - Best match genome (closest reference), optionally merged by ANI neighborhood
    - Novelty band (5% intervals)
    - Taxonomic call (Novel Species or Novel Genus)

    Adjacent bands sharing the same (representative genome, taxonomic call)
    are merged when their novelty ranges overlap or nearly overlap, reducing
    boundary artifacts at hard band edges.

    This enables researchers to:
    1. Identify distinct novel organisms in their sample
    2. Prioritize clusters for further investigation
    3. Understand phylogenetic context of novel diversity

    Example:
        >>> analyzer = NovelDiversityAnalyzer(
        ...     classifications=df,
        ...     metadata=genome_metadata,
        ...     ani_matrix=ani_matrix_obj,
        ... )
        >>> clusters = analyzer.cluster_novel_reads()
        >>> high_priority = [c for c in clusters if c.confidence == "High"]
    """

    def __init__(
        self,
        classifications: pl.DataFrame,
        metadata: GenomeMetadata | None = None,
        ani_matrix: ANIMatrix | None = None,
        novelty_band_size: float = 5.0,
        min_cluster_size: int = 3,
        genome_neighborhood_threshold: float = 80.0,
    ) -> None:
        """
        Initialize the analyzer.

        Args:
            classifications: DataFrame with classification results.
                Required columns: read_id, best_match_genome, novelty_index,
                placement_uncertainty, taxonomic_call
                Optional: discovery_score, confidence_score
            metadata: Genome metadata for species/genus/family lookups
            ani_matrix: ANI matrix for genome neighborhood merging
            novelty_band_size: Size of novelty bands in percent (default 5.0)
            min_cluster_size: Minimum reads to form a cluster (default 3)
            genome_neighborhood_threshold: ANI threshold for grouping genomes
                into the same neighborhood (default 80.0)
        """
        self._df = classifications
        self._metadata = metadata
        self._ani_matrix = ani_matrix
        self._novelty_band_size = novelty_band_size
        self._min_cluster_size = min_cluster_size
        self._genome_neighborhood_threshold = genome_neighborhood_threshold
        self._clusters: list[NovelCluster] | None = None
        self._summary: NovelDiversitySummary | None = None

        # Validate required columns
        required = {"read_id", "best_match_genome", "novelty_index",
                    "placement_uncertainty", "taxonomic_call"}
        missing = required - set(classifications.columns)
        if missing:
            msg = f"Classification DataFrame missing required columns: {missing}"
            raise ValueError(msg)

    # ------------------------------------------------------------------
    # Genome neighborhood merging
    # ------------------------------------------------------------------

    def _build_genome_neighborhoods(
        self,
        novel_genomes: list[str],
    ) -> dict[str, str]:
        """
        Group genomes into neighborhoods using union-find on the ANI matrix.

        Genomes with pairwise ANI >= ``genome_neighborhood_threshold`` are
        placed in the same connected component.  The representative of each
        component is the genome with the most novel reads.

        Args:
            novel_genomes: Unique best_match_genome values from novel reads.

        Returns:
            Mapping from each genome to its neighborhood representative.
        """
        if self._ani_matrix is None or len(novel_genomes) <= 1:
            return {g: g for g in novel_genomes}

        # Count reads per genome for representative selection
        novel_df = self._df.filter(
            pl.col("taxonomic_call").is_in(["Novel Species", "Novel Genus"])
        )
        genome_counts: dict[str, int] = {}
        for g, cnt in (
            novel_df.group_by("best_match_genome")
            .agg(pl.len().alias("n"))
            .iter_rows()
        ):
            genome_counts[g] = cnt

        # Union-find
        parent: dict[str, str] = {g: g for g in novel_genomes}

        def find(x: str) -> str:
            while parent[x] != x:
                parent[x] = parent[parent[x]]  # path compression
                x = parent[x]
            return x

        def union(x: str, y: str) -> None:
            px, py = find(x), find(y)
            if px != py:
                parent[px] = py

        threshold = self._genome_neighborhood_threshold

        # Check all pairs in the ANI matrix
        for i, g1 in enumerate(novel_genomes):
            if not self._ani_matrix.has_genome(g1):
                continue
            for g2 in novel_genomes[i + 1:]:
                if not self._ani_matrix.has_genome(g2):
                    continue
                ani = self._ani_matrix.get_ani(g1, g2)
                if ani >= threshold:
                    union(g1, g2)

        # Group by root
        components: dict[str, list[str]] = {}
        for g in novel_genomes:
            root = find(g)
            components.setdefault(root, []).append(g)

        # Pick representative = genome with most novel reads in each component
        mapping: dict[str, str] = {}
        for members in components.values():
            representative = max(members, key=lambda g: genome_counts.get(g, 0))
            for g in members:
                mapping[g] = representative

        n_merged = sum(1 for v in components.values() if len(v) > 1)
        if n_merged > 0:
            logger.info(
                f"Merged {sum(len(v) for v in components.values())} genomes "
                f"into {len(components)} neighborhoods ({n_merged} multi-genome)"
            )

        return mapping

    # ------------------------------------------------------------------
    # Adjacent-band merging
    # ------------------------------------------------------------------

    def _merge_adjacent_bands(
        self,
        rows: list[dict],
    ) -> list[dict]:
        """
        Merge clusters that share (representative_genome, taxonomic_call)
        and occupy adjacent novelty bands with overlapping or near-overlapping
        novelty ranges.

        Args:
            rows: List of row dicts from the initial grouping, each containing
                  representative_genome, taxonomic_call, novelty_band,
                  read_count, mean_novelty, novelty_min, novelty_max,
                  mean_uncertainty, read_ids, and optionally mean_discovery,
                  mean_bayesian_conf, contributing_genomes.

        Returns:
            Merged list of row dicts.
        """
        if not rows:
            return rows

        band_size = self._novelty_band_size

        # Sort by (representative_genome, taxonomic_call, novelty_band)
        rows_sorted = sorted(
            rows,
            key=lambda r: (r["representative_genome"], r["taxonomic_call"], r["novelty_band"]),
        )

        merged: list[dict] = []
        current = rows_sorted[0].copy()

        for nxt in rows_sorted[1:]:
            same_group = (
                nxt["representative_genome"] == current["representative_genome"]
                and nxt["taxonomic_call"] == current["taxonomic_call"]
            )
            adjacent_band = (nxt["novelty_band"] - current["novelty_band"]) <= band_size
            # Overlapping or near-overlapping novelty ranges (2% tolerance)
            ranges_close = (current["novelty_max"] >= nxt["novelty_min"] - 2.0)

            if same_group and adjacent_band and ranges_close:
                # Merge nxt into current
                total_reads = current["read_count"] + nxt["read_count"]
                w_curr = current["read_count"] / total_reads
                w_nxt = nxt["read_count"] / total_reads

                current["mean_novelty"] = (
                    current["mean_novelty"] * w_curr + nxt["mean_novelty"] * w_nxt
                )
                current["mean_uncertainty"] = (
                    current["mean_uncertainty"] * w_curr + nxt["mean_uncertainty"] * w_nxt
                )
                current["novelty_min"] = min(current["novelty_min"], nxt["novelty_min"])
                current["novelty_max"] = max(current["novelty_max"], nxt["novelty_max"])
                current["read_count"] = total_reads
                current["read_ids"] = current["read_ids"] + nxt["read_ids"]

                # Update novelty_band to the band of the lower edge
                current["novelty_band"] = min(current["novelty_band"], nxt["novelty_band"])

                # Merge optional fields
                if "mean_discovery" in current and current["mean_discovery"] is not None:
                    nxt_disc = nxt.get("mean_discovery")
                    if nxt_disc is not None:
                        current["mean_discovery"] = (
                            current["mean_discovery"] * w_curr + nxt_disc * w_nxt
                        )

                if "mean_bayesian_conf" in current and current["mean_bayesian_conf"] is not None:
                    nxt_bc = nxt.get("mean_bayesian_conf")
                    if nxt_bc is not None:
                        current["mean_bayesian_conf"] = (
                            current["mean_bayesian_conf"] * w_curr + nxt_bc * w_nxt
                        )

                # Union contributing genomes
                curr_cg = set(current.get("contributing_genomes") or [])
                nxt_cg = set(nxt.get("contributing_genomes") or [])
                current["contributing_genomes"] = sorted(curr_cg | nxt_cg)

            else:
                merged.append(current)
                current = nxt.copy()

        merged.append(current)
        return merged

    # ------------------------------------------------------------------
    # Main clustering
    # ------------------------------------------------------------------

    def cluster_novel_reads(self) -> list[NovelCluster]:
        """
        Cluster novel reads into putative novel taxa.

        Returns:
            List of NovelCluster objects representing putative novel taxa

        Algorithm:
            1. Filter to Novel Species and Novel Genus reads only
            2. Compute novelty band for each read
            3. Map genomes to neighborhood representatives (ANI-aware)
            4. Group by (representative_genome, novelty_band, taxonomic_call)
            5. Filter by minimum cluster size
            6. Merge adjacent bands with overlapping novelty ranges
            7. Generate cluster IDs and suggested names
            8. Assign confidence ratings (Bayesian-informed when available)
            9. Infer phylogenetic context
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
        has_ambiguous_hits = "num_ambiguous_hits" in novel_df.columns

        if has_ambiguous_hits:
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

        # Step 3: Build genome neighborhoods and map to representatives
        unique_genomes = novel_df["best_match_genome"].unique().to_list()
        genome_neighborhood = self._build_genome_neighborhoods(unique_genomes)

        novel_df = novel_df.with_columns(
            pl.col("best_match_genome")
            .replace_strict(genome_neighborhood, default=pl.col("best_match_genome"))
            .alias("representative_genome")
        )

        # Step 4: Group by clustering key (now using representative_genome)
        has_discovery = "discovery_score" in novel_df.columns
        has_bayesian_conf = "confidence_score" in novel_df.columns

        agg_exprs = [
            pl.len().alias("read_count"),
            pl.col("novelty_index").mean().alias("mean_novelty"),
            pl.col("novelty_index").min().alias("novelty_min"),
            pl.col("novelty_index").max().alias("novelty_max"),
            pl.col("effective_uncertainty").mean().alias("mean_uncertainty"),
            pl.col("read_id").alias("read_ids"),
            pl.col("best_match_genome").unique().alias("contributing_genomes"),
        ]

        if has_discovery:
            agg_exprs.append(
                pl.col("discovery_score").mean().alias("mean_discovery")
            )

        if has_bayesian_conf:
            agg_exprs.append(
                pl.col("confidence_score").mean().alias("mean_bayesian_conf")
            )

        grouped = (
            novel_df
            .group_by(["representative_genome", "novelty_band", "taxonomic_call"])
            .agg(agg_exprs)
        )

        # Step 5: Filter by minimum cluster size
        grouped = grouped.filter(pl.col("read_count") >= self._min_cluster_size)

        if len(grouped) == 0:
            logger.info(
                f"No clusters with >= {self._min_cluster_size} reads found"
            )
            self._clusters = []
            return self._clusters

        # Convert to list of dicts for band merging
        rows = grouped.to_dicts()

        # Step 6: Merge adjacent bands
        rows = self._merge_adjacent_bands(rows)

        # Re-filter after merging (shouldn't shrink, but be safe)
        rows = [r for r in rows if r["read_count"] >= self._min_cluster_size]

        # Sort by read count descending
        rows.sort(key=lambda r: r["read_count"], reverse=True)

        # Step 7-9: Build NovelCluster objects
        clusters: list[NovelCluster] = []
        nsp_counter = 0
        ngn_counter = 0

        for row in rows:
            genome = row["representative_genome"]
            tax_call = row["taxonomic_call"]
            novelty_band = row["novelty_band"]
            read_count = row["read_count"]
            mean_novelty = row["mean_novelty"]
            novelty_min = row["novelty_min"]
            novelty_max = row["novelty_max"]
            mean_uncertainty = row["mean_uncertainty"]
            read_ids = row["read_ids"]
            mean_discovery = row.get("mean_discovery")
            mean_bayesian_conf = row.get("mean_bayesian_conf")
            contributing_genomes = row.get("contributing_genomes") or []

            # Ensure contributing_genomes is a flat list of strings
            if contributing_genomes and isinstance(contributing_genomes[0], list):
                flat: list[str] = []
                for item in contributing_genomes:
                    if isinstance(item, list):
                        flat.extend(item)
                    else:
                        flat.append(item)
                contributing_genomes = sorted(set(flat))
            else:
                contributing_genomes = sorted(set(contributing_genomes))

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
                read_count, mean_uncertainty, mean_discovery,
                mean_bayesian_confidence=mean_bayesian_conf,
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
                mean_bayesian_confidence=mean_bayesian_conf,
                contributing_genomes=contributing_genomes,
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
        *,
        mean_bayesian_confidence: float | None = None,
    ) -> Literal["High", "Medium", "Low"]:
        """
        Assign a confidence rating to a cluster.

        When Bayesian per-read confidence scores are available, ratings are
        based on ``mean_bayesian_confidence`` combined with read count:
            High:   mean_bayesian_confidence >= 70 AND read_count >= 10
            Medium: mean_bayesian_confidence >= 50 AND read_count >= 5
            Low:    everything else

        Without Bayesian scores, the legacy criteria apply:
            High:   read_count >= 10, mean_uncertainty < 5%, discovery >= 75
            Medium: read_count >= 5, mean_uncertainty < 10%, discovery >= 50
            Low:    everything else

        Args:
            read_count: Number of reads in cluster
            mean_uncertainty: Mean placement uncertainty
            mean_discovery: Mean discovery score (optional)
            mean_bayesian_confidence: Mean per-read Bayesian confidence (0-100)

        Returns:
            Confidence rating string
        """
        # Prefer Bayesian confidence when available
        if mean_bayesian_confidence is not None:
            if mean_bayesian_confidence >= 70 and read_count >= 10:
                return "High"
            elif mean_bayesian_confidence >= 50 and read_count >= 5:
                return "Medium"
            return "Low"

        # Legacy: use discovery score if available, otherwise simplified criteria
        if mean_discovery is not None:
            if read_count >= 10 and mean_uncertainty < 5 and mean_discovery >= 75:
                return "High"
            elif read_count >= 5 and mean_uncertainty < 10 and mean_discovery >= 50:
                return "Medium"
        else:
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
