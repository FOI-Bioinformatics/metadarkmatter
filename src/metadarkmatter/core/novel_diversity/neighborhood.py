"""
Phylogenetic neighborhood analysis for novel clusters.

Estimates the evolutionary context of each novel cluster by triangulating
ANI distances to reference genomes, grouping references by genus, and
computing placement support metrics. The result is a compact neighborhood
profile attached to each cluster.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from metadarkmatter.core.novel_diversity.models import (
    GenusDistance,
    PhylogeneticNeighborhood,
)

if TYPE_CHECKING:
    from metadarkmatter.core.classification.ani_matrix import ANIMatrix
    from metadarkmatter.core.novel_diversity.models import NovelCluster

logger = logging.getLogger(__name__)

_MAX_NEAREST_GENERA = 5
_ANI_FLOOR = 70.0


class PhylogeneticNeighborhoodAnalyzer:
    """Compute neighborhood profiles for novel clusters.

    For each cluster the analyzer estimates ANI to all reference genomes
    in the matrix, groups them by genus, and derives a placement support
    score together with a human-readable context string.

    Args:
        ani_matrix: Reference genome ANI matrix.
        genus_map: Optional mapping from genome accession to genus name.
            When absent, genomes are grouped into anonymous neighborhoods
            using a union-find at the genus boundary threshold.
        genus_boundary: ANI threshold separating genera (default 80.0).
    """

    def __init__(
        self,
        ani_matrix: ANIMatrix,
        genus_map: dict[str, str] | None = None,
        genus_boundary: float = 80.0,
    ) -> None:
        self._ani_matrix = ani_matrix
        self._genus_map = genus_map
        self._genus_boundary = genus_boundary

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def analyze(self, clusters: list[NovelCluster]) -> list[NovelCluster]:
        """Attach a phylogenetic neighborhood profile to each cluster.

        Args:
            clusters: Novel clusters produced by the clustering stage.

        Returns:
            New cluster objects with the ``neighborhood`` field populated.
        """
        genus_map = self._genus_map or self._infer_genus_groups()

        results: list[NovelCluster] = []
        for cluster in clusters:
            nbr = self._build_neighborhood(cluster, genus_map)
            results.append(cluster.model_copy(update={"neighborhood": nbr}))

        logger.info(
            "Computed neighborhood profiles for %d cluster(s)", len(results)
        )
        return results

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _build_neighborhood(
        self,
        cluster: NovelCluster,
        genus_map: dict[str, str],
    ) -> PhylogeneticNeighborhood:
        """Build a single neighborhood profile for *cluster*."""
        nearest_genome = cluster.nearest_genome
        nearest_ani = 100.0 - cluster.mean_novelty_index

        # Step 1 -- estimate ANI from cluster to every reference genome
        estimated: dict[str, float] = {}
        for genome in self._ani_matrix.genomes:
            if genome == nearest_genome:
                estimated[genome] = nearest_ani
            else:
                nearest_to_ref = self._ani_matrix.get_ani(
                    nearest_genome, genome
                )
                divergence_novel = 100.0 - nearest_ani
                divergence_ref = 100.0 - nearest_to_ref
                est = 100.0 - (divergence_novel + divergence_ref * 0.7)
                estimated[genome] = max(_ANI_FLOOR, est)

        # Step 2 -- group by genus and take the best ANI per genus
        genus_genomes: dict[str, list[str]] = {}
        for genome in self._ani_matrix.genomes:
            genus = genus_map.get(genome, "Unknown")
            genus_genomes.setdefault(genus, []).append(genome)

        genus_distances: list[GenusDistance] = []
        for genus, members in genus_genomes.items():
            best_genome = max(members, key=lambda g: estimated[g])
            best_ani = estimated[best_genome]
            genus_distances.append(
                GenusDistance(
                    genus=genus,
                    representative_genome=best_genome,
                    estimated_ani=round(best_ani, 2),
                    num_genomes_in_genus=len(members),
                )
            )

        genus_distances.sort(key=lambda gd: gd.estimated_ani, reverse=True)
        nearest_genera = genus_distances[:_MAX_NEAREST_GENERA]

        # Step 3 -- isolation score
        if len(nearest_genera) >= 2:
            isolation_score = (
                nearest_genera[0].estimated_ani - nearest_genera[1].estimated_ani
            )
        else:
            isolation_score = 0.0

        # Step 4 -- neighborhood density
        density_threshold = self._genus_boundary - 5.0
        neighborhood_density = sum(
            1 for gd in genus_distances if gd.estimated_ani >= density_threshold
        )

        # Step 5 -- placement support (0-100)
        placement_support = self._compute_placement_support(
            cluster, nearest_genera, isolation_score
        )

        # Step 6 -- context text
        phylogenetic_context = self._generate_context_text(
            cluster, nearest_genera, isolation_score, placement_support
        )

        return PhylogeneticNeighborhood(
            cluster_id=cluster.cluster_id,
            nearest_genera=nearest_genera,
            placement_support=round(placement_support, 1),
            isolation_score=round(isolation_score, 2),
            neighborhood_density=neighborhood_density,
            phylogenetic_context=phylogenetic_context,
            genus_boundary_ani=self._genus_boundary,
        )

    # ------------------------------------------------------------------

    def _compute_placement_support(
        self,
        cluster: NovelCluster,
        nearest_genera: list[GenusDistance],
        isolation_score: float,
    ) -> float:
        """Derive a composite placement support score (0-100)."""
        confidence_raw = cluster.mean_bayesian_confidence or 50.0

        if cluster.taxonomic_call == "Novel Genus":
            isolation_component = min(isolation_score / 10.0, 1.0) * 40
            top_ani = nearest_genera[0].estimated_ani if nearest_genera else 0.0
            gap = self._genus_boundary - top_ani
            boundary_component = min(max(gap, 0.0) / 10.0, 1.0) * 30
            read_component = min(cluster.read_count / 20, 1.0) * 15
            confidence_component = (confidence_raw / 100.0) * 15
        else:
            # Novel Species -- boundary component is not applicable
            isolation_component = min(isolation_score / 10.0, 1.0) * 55
            boundary_component = 0.0
            read_component = min(cluster.read_count / 20, 1.0) * 22.5
            confidence_component = (confidence_raw / 100.0) * 22.5

        total = (
            isolation_component
            + boundary_component
            + read_component
            + confidence_component
        )
        return min(total, 100.0)

    # ------------------------------------------------------------------

    @staticmethod
    def _generate_context_text(
        cluster: NovelCluster,
        nearest_genera: list[GenusDistance],
        isolation_score: float,
        placement_support: float,
    ) -> str:
        """Produce a one-line human-readable placement summary."""
        if not nearest_genera:
            return f"No reference genomes available. Support: {placement_support:.0f}/100."

        top = nearest_genera[0]

        if cluster.taxonomic_call == "Novel Genus":
            if len(nearest_genera) >= 2:
                second = nearest_genera[1]
                return (
                    f"Sister to {top.genus} ({top.estimated_ani:.0f}% ANI), "
                    f"{isolation_score:.0f}% isolated from {second.genus}. "
                    f"Support: {placement_support:.0f}/100."
                )
            return (
                f"Sister to {top.genus} ({top.estimated_ani:.0f}% ANI), "
                f"no close neighboring genera. "
                f"Support: {placement_support:.0f}/100."
            )

        # Novel Species
        return (
            f"Within {top.genus}, closest to {top.representative_genome} "
            f"({top.estimated_ani:.0f}% ANI). "
            f"Support: {placement_support:.0f}/100."
        )

    # ------------------------------------------------------------------

    def _infer_genus_groups(self) -> dict[str, str]:
        """Group genomes into anonymous genera using union-find at genus boundary.

        When no genus_map is provided, genomes sharing ANI above the genus
        boundary threshold are placed in the same group.  Groups are labeled
        as "Group 1", "Group 2", etc.

        Returns:
            Mapping from genome accession to group label.
        """
        genomes = sorted(self._ani_matrix.genomes)
        if not genomes:
            return {}

        parent: dict[str, str] = {g: g for g in genomes}

        def find(x: str) -> str:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(x: str, y: str) -> None:
            px, py = find(x), find(y)
            if px != py:
                parent[px] = py

        threshold = self._genus_boundary
        for i, g1 in enumerate(genomes):
            for g2 in genomes[i + 1 :]:
                ani = self._ani_matrix.get_ani(g1, g2)
                if ani >= threshold:
                    union(g1, g2)

        # Assign sequential group labels
        components: dict[str, list[str]] = {}
        for g in genomes:
            root = find(g)
            components.setdefault(root, []).append(g)

        group_map: dict[str, str] = {}
        for idx, members in enumerate(sorted(components.values(), key=lambda m: m[0]), start=1):
            label = f"Group {idx}"
            for g in members:
                group_map[g] = label

        return group_map
