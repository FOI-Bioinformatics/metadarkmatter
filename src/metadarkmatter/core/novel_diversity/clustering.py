"""
Novel diversity clustering for identifying putative novel taxa.

Novel reads (Novel Species / Novel Genus) are grouped into "novel lineage
groups". Three grouping methods are available:

- ``genus`` (default when metadata gives genus names): reads are grouped by
  the GTDB genus of their nearest reference; within a genus, references are
  linked when their ANI is >= ``genome_neighborhood_threshold`` (80%), so
  distant lineages of one genus stay apart. Reads of one novel organism
  usually spread over several references of the same genus, so this keeps
  them together without chaining a whole family, as unrestricted
  single-linkage ANI merging does in dense clades. Groups are not split by
  call or novelty band; the group's call is the majority call of its reads.
- ``reference``: one group per nearest reference (species representative).
  Used when no genus information is available.
- ``legacy``: the previous method, (ANI neighbourhood of the nearest
  reference) x 5% novelty band x call, with adjacent-band merging.

A group is a candidate lineage, not a confirmed single organism: a novel
genus whose reads hit one existing genus is pooled with that genus.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any, Literal

import polars as pl

from metadarkmatter.core.novel_diversity.models import (
    NovelCluster,
    NovelDiversitySummary,
)

if TYPE_CHECKING:
    from metadarkmatter.core.classification.ani_matrix import ANIMatrix
    from metadarkmatter.core.metadata import GenomeMetadata

logger = logging.getLogger(__name__)

Grouping = Literal["auto", "genus", "reference", "legacy"]
GROUPING_CHOICES: tuple[str, ...] = ("auto", "genus", "reference", "legacy")


class NovelDiversityAnalyzer:
    """
    Analyzes novel diversity by grouping novel reads into novel lineage groups.

    The grouping method is set by ``grouping`` (see the module docstring):
    ``genus`` (default with genus metadata), ``reference`` (fallback without
    it) or ``legacy`` (neighbourhood x novelty band x call).

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
        merge_band_tolerance: float = 2.0,
        grouping: Grouping | str = "auto",
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
            genome_neighborhood_threshold: ANI threshold for linking
                references: within a genus (genus grouping) or across all
                references (legacy grouping). Default 80.0
            merge_band_tolerance: Adjacent-band merge tolerance (legacy only)
            grouping: "genus", "reference", "legacy", or "auto" (genus when
                metadata provides genus names for the novel reads' references,
                otherwise reference). See the module docstring.
        """
        if novelty_band_size <= 0 or float(novelty_band_size) != int(novelty_band_size):
            msg = f"novelty_band_size must be a positive whole number, got {novelty_band_size}"
            raise ValueError(msg)
        if grouping not in GROUPING_CHOICES:
            msg = f"grouping must be one of {GROUPING_CHOICES}, got {grouping!r}"
            raise ValueError(msg)
        self._df = classifications
        self._metadata = metadata
        self._ani_matrix = ani_matrix
        self._novelty_band_size = novelty_band_size
        self._min_cluster_size = min_cluster_size
        self._genome_neighborhood_threshold = genome_neighborhood_threshold
        self._merge_band_tolerance = merge_band_tolerance
        self._grouping_requested: str = grouping
        self.grouping: str | None = None  # resolved method, set by cluster_novel_reads
        self._clusters: list[NovelCluster] | None = None
        self._summary: NovelDiversitySummary | None = None

        # Validate required columns
        required = {"read_id", "best_match_genome", "novelty_index",
                    "placement_uncertainty", "taxonomic_call"}
        missing = required - set(classifications.columns)
        if missing:
            msg = f"Classification DataFrame missing required columns: {missing}"
            raise ValueError(msg)

    def _representative_mapping(self) -> dict[str, str]:
        """Return {accession: species representative} from metadata, if any.

        Empty when no metadata is available or no genome maps to a different
        representative (identity mapping), so the default path is unchanged.
        """
        if self._metadata is None:
            return {}
        try:
            mapping = self._metadata.build_representative_mapping()
        except Exception:  # pragma: no cover - defensive for stub metadata
            return {}
        if not isinstance(mapping, dict):
            return {}
        return {k: v for k, v in mapping.items() if k != v}

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
        rep_map = self._representative_mapping()
        for g, cnt in (
            novel_df.group_by("best_match_genome")
            .agg(pl.len().alias("n"))
            .iter_rows()
        ):
            g = rep_map.get(g, g)
            genome_counts[g] = genome_counts.get(g, 0) + cnt

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

        # Pick representative = genome with most novel reads in each component.
        # Break count ties by accession so the choice is deterministic across
        # processes (polars .unique() does not guarantee order).
        mapping: dict[str, str] = {}
        for members in components.values():
            representative = max(members, key=lambda g: (genome_counts.get(g, 0), g))
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
        rows: list[dict[str, Any]],
    ) -> list[dict[str, Any]]:
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

        merged: list[dict[str, Any]] = []
        current = rows_sorted[0].copy()
        # Highest band absorbed into ``current``. Adjacency is judged against
        # it (not the reported lower-edge band) so chains of consecutive bands
        # merge transitively.
        last_band = current["novelty_band"]

        for nxt in rows_sorted[1:]:
            same_group = (
                nxt["representative_genome"] == current["representative_genome"]
                and nxt["taxonomic_call"] == current["taxonomic_call"]
            )
            adjacent_band = (nxt["novelty_band"] - last_band) <= band_size
            # Overlapping or near-overlapping novelty ranges
            ranges_close = (
                current["novelty_max"] >= nxt["novelty_min"] - self._merge_band_tolerance
            )

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
                last_band = max(last_band, nxt["novelty_band"])

                # Merge optional fields
                if current.get("mean_inferred") is not None and nxt.get("mean_inferred") is not None:
                    current["mean_inferred"] = (
                        current["mean_inferred"] * w_curr + nxt["mean_inferred"] * w_nxt
                    )
                elif nxt.get("mean_inferred") is None:
                    pass
                else:
                    current["mean_inferred"] = nxt["mean_inferred"]
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
                last_band = current["novelty_band"]

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
            2. Map hit genomes to their species representatives
            3. Group reads by the resolved grouping method (genus, reference
               or legacy neighbourhood x band x call)
            4. Filter by minimum group size
            5. Generate IDs (NSP_/NGN_ by the group's majority call), names,
               confidence ratings and phylogenetic context
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

        # Step 2b: Report the real placement uncertainty.  Single-hit reads
        # have placement_uncertainty == 0 (no competing genome); that value is
        # reported as-is rather than replaced by a novelty-derived guess.
        novel_df = novel_df.with_columns(
            pl.col("placement_uncertainty").fill_null(0.0).alias("effective_uncertainty")
        )

        # Step 2c: Map each hit genome to its species representative so that
        # clustering, neighborhood analysis and tree placement operate on
        # genomes present in the (representative-only) ANI matrix.  The raw
        # genome is kept for display in ``contributing_genomes``.
        rep_map = self._representative_mapping()
        if rep_map:
            novel_df = novel_df.with_columns(
                pl.col("best_match_genome")
                .replace_strict(rep_map, default=pl.col("best_match_genome"))
                .alias("cluster_genome")
            )
        else:
            novel_df = novel_df.with_columns(
                pl.col("best_match_genome").alias("cluster_genome")
            )

        grouping = self._resolve_grouping(novel_df)
        self.grouping = grouping
        if grouping == "legacy":
            rows = self._legacy_rows(novel_df)
        else:
            rows = self._lineage_rows(novel_df, grouping)

        # Filter by minimum cluster size
        rows = [r for r in rows if r["read_count"] >= self._min_cluster_size]
        return self._build_clusters(rows, grouping)

    # ------------------------------------------------------------------
    # Grouping methods
    # ------------------------------------------------------------------

    def _genus_of(self, accession: str) -> str | None:
        if self._metadata is None:
            return None
        genus = self._metadata.get_genus(accession)
        if not genus or genus == "Unknown":
            return None
        return str(genus)

    def _resolve_grouping(self, novel_df: pl.DataFrame) -> str:
        """Resolve ``auto`` to ``genus`` or ``reference``."""
        requested = self._grouping_requested
        if requested in ("reference", "legacy"):
            return requested
        genomes = novel_df["cluster_genome"].unique().to_list()
        has_genus = any(self._genus_of(g) is not None for g in genomes)
        if has_genus:
            return "genus"
        logger.warning(
            "No genus information for the novel reads' references (no --metadata, "
            "or no genus column); grouping novel reads by nearest reference "
            "instead of by genus. Reads of one novel organism may be split "
            "across several groups."
        )
        return "reference"

    def _genus_components(self, genomes: list[str]) -> dict[str, str]:
        """Map each genome to a group key: its genus, split by ANI within the genus.

        References of one genus are linked (single linkage) when their ANI is
        >= ``genome_neighborhood_threshold``; links never cross genera, so a
        dense family cannot chain into one group. Without an ANI matrix every
        genus is one group; with one, a reference absent from it is linked to
        nothing. Genomes without a genus form their own group.
        """
        parent = {g: g for g in genomes}

        def find(x: str) -> str:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        genus_of = {g: self._genus_of(g) for g in genomes}
        by_genus: dict[str, list[str]] = {}
        for g in genomes:
            if genus_of[g] is not None:
                by_genus.setdefault(str(genus_of[g]), []).append(g)
        ani = self._ani_matrix
        for members in by_genus.values():
            for i, a in enumerate(members):
                for b in members[i + 1:]:
                    # Without an ANI matrix the whole genus is one group. With
                    # one, link only on a measured ANI; a reference missing
                    # from the matrix must not chain distant lineages together.
                    linked = ani is None or (
                        ani.has_genome(a)
                        and ani.has_genome(b)
                        and ani.get_ani(a, b) >= self._genome_neighborhood_threshold
                    )
                    if linked:
                        ra, rb = find(a), find(b)
                        if ra != rb:
                            parent[max(ra, rb)] = min(ra, rb)
        return {
            g: (f"genus:{genus_of[g]}:{find(g)}" if genus_of[g] is not None else f"genome:{g}")
            for g in genomes
        }

    def _lineage_rows(self, novel_df: pl.DataFrame, grouping: str) -> list[dict[str, Any]]:
        """One group per genus (or per reference) of the nearest reference.

        Genomes without a genus in the metadata form their own group.
        """
        genomes = sorted(novel_df["cluster_genome"].unique().to_list())
        if grouping == "genus":
            key_map = self._genus_components(genomes)
        else:
            key_map = {g: f"genome:{g}" for g in genomes}
        df = novel_df.with_columns(
            pl.col("cluster_genome").replace_strict(key_map).alias("_group_key")
        )

        agg_exprs = [
            pl.len().alias("read_count"),
            pl.col("novelty_index").mean().alias("mean_novelty"),
            pl.col("novelty_index").min().alias("novelty_min"),
            pl.col("novelty_index").max().alias("novelty_max"),
            pl.col("effective_uncertainty").mean().alias("mean_uncertainty"),
            pl.col("read_id").sort().alias("read_ids"),
            pl.col("best_match_genome").unique().sort().alias("contributing_genomes"),
            (pl.col("taxonomic_call") == "Novel Species").sum().alias("n_species"),
            (pl.col("taxonomic_call") == "Novel Genus").sum().alias("n_genus"),
        ]
        if "inferred_uncertainty" in df.columns:
            agg_exprs.append(pl.col("inferred_uncertainty").mean().alias("mean_inferred"))
        if "discovery_score" in df.columns:
            agg_exprs.append(pl.col("discovery_score").mean().alias("mean_discovery"))
        if "confidence_score" in df.columns:
            agg_exprs.append(pl.col("confidence_score").mean().alias("mean_bayesian_conf"))
        grouped = df.group_by("_group_key").agg(agg_exprs).sort("_group_key")

        # Nearest reference of each group: the genome with the most reads
        # (ties broken by accession, so the choice is reproducible).
        counts = (
            df.group_by(["_group_key", "cluster_genome"]).agg(pl.len().alias("n"))
            .sort(["_group_key", "n", "cluster_genome"], descending=[False, True, False])
        )
        top_genome: dict[str, str] = {}
        for key, genome, _n in counts.iter_rows():
            top_genome.setdefault(key, genome)

        rows: list[dict[str, Any]] = []
        band = self._novelty_band_size
        for row in grouped.to_dicts():
            n_sp, n_gn = int(row.pop("n_species")), int(row.pop("n_genus"))
            row["taxonomic_call"] = "Novel Genus" if n_gn > n_sp else "Novel Species"
            row["novel_species_reads"] = n_sp
            row["novel_genus_reads"] = n_gn
            row["representative_genome"] = top_genome[row.pop("_group_key")]
            row["novelty_band"] = int(row["novelty_min"] // band) * int(band)
            rows.append(row)
        return rows

    def _legacy_rows(self, novel_df: pl.DataFrame) -> list[dict[str, Any]]:
        """Previous method: neighbourhood x novelty band x call."""
        # Step 3: Build genome neighborhoods and map to representatives.
        # Sort for a deterministic genome order (polars .unique() is unordered),
        # so neighborhood construction and tie-breaks are reproducible.
        unique_genomes = sorted(novel_df["cluster_genome"].unique().to_list())
        genome_neighborhood = self._build_genome_neighborhoods(unique_genomes)

        novel_df = novel_df.with_columns(
            pl.col("cluster_genome")
            .replace_strict(genome_neighborhood, default=pl.col("cluster_genome"))
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

        if "inferred_uncertainty" in novel_df.columns:
            agg_exprs.append(
                pl.col("inferred_uncertainty").mean().alias("mean_inferred")
            )

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

        # Step 5: Merge adjacent bands first, so that sub-threshold groups in
        # neighbouring bands can combine into a cluster that passes the
        # minimum size (filtering first would silently drop those reads).
        return self._merge_adjacent_bands(grouped.to_dicts())

    def _build_clusters(self, rows: list[dict[str, Any]], grouping: str) -> list[NovelCluster]:
        if not rows:
            logger.info(
                f"No clusters with >= {self._min_cluster_size} reads found"
            )
            self._clusters = []
            return self._clusters

        # Sort by read count descending, with deterministic tiebreakers
        # (representative genome, then call) so sequential cluster IDs are
        # assigned reproducibly when read counts tie.
        rows.sort(
            key=lambda r: (-r["read_count"], r["representative_genome"], r["taxonomic_call"])
        )

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
            mean_inferred = row.get("mean_inferred")
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
                mean_inferred_uncertainty=mean_inferred,
                mean_discovery_score=mean_discovery,
                suggested_name=suggested_name,
                confidence=confidence,
                phylogenetic_context=phylo_context,
                read_ids=read_ids,
                mean_bayesian_confidence=mean_bayesian_conf,
                contributing_genomes=contributing_genomes,
                grouping=grouping,
                novel_species_reads=row.get("novel_species_reads"),
                novel_genus_reads=row.get("novel_genus_reads"),
            )
            clusters.append(cluster)

        logger.info(f"Created {len(clusters)} novel lineage groups ({grouping} grouping)")
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

        # Count reads by call (whole table, including reads in groups below
        # the minimum size)
        novel_df = self._df.filter(
            pl.col("taxonomic_call").is_in(["Novel Species", "Novel Genus"])
        )
        nsp_reads = len(novel_df.filter(pl.col("taxonomic_call") == "Novel Species"))
        ngn_reads = len(novel_df.filter(pl.col("taxonomic_call") == "Novel Genus"))

        if not clusters:
            self._summary = NovelDiversitySummary(
                total_novel_reads=nsp_reads + ngn_reads,
                novel_species_reads=nsp_reads,
                novel_genus_reads=ngn_reads,
            )
            return self._summary

        # Count clusters by type and confidence
        nsp_clusters = sum(1 for c in clusters if c.taxonomic_call == "Novel Species")
        ngn_clusters = sum(1 for c in clusters if c.taxonomic_call == "Novel Genus")
        high_conf = sum(1 for c in clusters if c.confidence == "High")
        med_conf = sum(1 for c in clusters if c.confidence == "Medium")
        low_conf = sum(1 for c in clusters if c.confidence == "Low")

        # Unique genera/families with novel diversity
        genera = {c.nearest_genus for c in clusters if c.taxonomic_call == "Novel Species"}
        families = {c.nearest_family for c in clusters if c.taxonomic_call == "Novel Genus"}

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
            grouped_reads=sum(sizes),
            novel_species_cluster_reads=sum(
                c.read_count for c in clusters if c.taxonomic_call == "Novel Species"
            ),
            novel_genus_cluster_reads=sum(
                c.read_count for c in clusters if c.taxonomic_call == "Novel Genus"
            ),
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

    def to_dict(self) -> dict[str, Any]:
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


def build_report_analyzer(
    classifications: pl.DataFrame,
    metadata: GenomeMetadata | None,
    ani_matrix_df: pl.DataFrame | None,
    grouping: Grouping | str = "auto",
) -> NovelDiversityAnalyzer:
    """Create the analyzer used by every report section.

    The Novel Diversity tab and the phylogeny tree both call this, so they
    share one definition of clusters (IDs, read counts, confidence).
    """
    ani_obj = None
    if ani_matrix_df is not None and len(ani_matrix_df) > 0:
        from metadarkmatter.core.classification.ani_matrix import ANIMatrix

        try:
            ani_obj = ANIMatrix.from_dataframe(ani_matrix_df)
        except Exception as e:
            logger.warning(f"Could not convert ANI matrix for clustering: {e}")
    return NovelDiversityAnalyzer(
        classifications=classifications,
        metadata=metadata,
        ani_matrix=ani_obj,
        novelty_band_size=5.0,
        min_cluster_size=3,
        grouping=grouping,
    )
