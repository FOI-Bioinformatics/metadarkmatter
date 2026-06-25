"""ReportGenerator novel section."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import polars as pl

from metadarkmatter.visualization.report.components import (
    build_phylogenetic_context_heatmap,
)
from metadarkmatter.visualization.report.novel_section import (
    NOVEL_EMPTY_TEMPLATE,
    PHYLOGENETIC_HEATMAP_INTRO_TEMPLATE,
    PHYLOGENETIC_HEATMAP_LEGEND_TEMPLATE,
    build_cluster_scatter_figure,
    build_cluster_table_html,
    build_novel_summary_html,
    build_sunburst_figure,
)
from metadarkmatter.visualization.report.sections._base import _ReportBase
from metadarkmatter.visualization.report.templates import (
    COLLAPSIBLE_PANEL_TEMPLATE,
    EMPTY_SECTION_TEMPLATE,
    KPI_CARD_TEMPLATE,
    KPI_STRIP_TEMPLATE,
    PLOT_CONTAINER_TEMPLATE,
    TAB_SECTION_TEMPLATE,
)

if TYPE_CHECKING:

    from metadarkmatter.core.classification.adaptive import AdaptiveGenusThreshold
    from metadarkmatter.core.novel_diversity.models import NovelCluster


logger = logging.getLogger(__name__)


class NovelMixin(_ReportBase):
    """Novel report-section builders."""

    def _build_novel_diversity_section(self) -> str:
        """Build the novel diversity tab section with cluster analysis."""
        from metadarkmatter.core.novel_diversity import NovelDiversityAnalyzer

        content_parts = []

        # Track genus boundary result across nested scopes so it can
        # be used in KPI cards and the scatter plot later.
        genus_boundary_ani: float = 80.0
        genus_boundary_result = None  # AdaptiveGenusThreshold or None

        # Create analyzer and cluster novel reads
        try:
            from metadarkmatter.core.classification.ani_matrix import ANIMatrix

            ani_matrix_obj = None
            if self.ani_matrix is not None and len(self.ani_matrix) > 0:
                try:
                    ani_matrix_obj = ANIMatrix.from_dataframe(self.ani_matrix)
                except Exception as e:
                    logger.warning(f"Could not convert ANI matrix for clustering: {e}")

            analyzer = NovelDiversityAnalyzer(
                classifications=self.df,
                metadata=self.genome_metadata,
                ani_matrix=ani_matrix_obj,
                novelty_band_size=5.0,
                min_cluster_size=3,
            )
            clusters = analyzer.cluster_novel_reads()
            summary = analyzer.get_summary()

            # Compute phylogenetic neighborhoods if ANI matrix available
            if ani_matrix_obj is not None:
                from metadarkmatter.core.novel_diversity.neighborhood import (
                    PhylogeneticNeighborhoodAnalyzer,
                )

                genus_map = None
                if self.genome_metadata is not None:
                    genus_map = {}
                    for genome in ani_matrix_obj.genomes:
                        genus = self.genome_metadata.get_genus(genome)
                        if genus:
                            genus_map[genome] = genus

                # Detect adaptive genus boundary (updates method-level variables)
                try:
                    from metadarkmatter.core.classification.adaptive import (
                        detect_genus_boundary,
                    )

                    genus_boundary_result = detect_genus_boundary(
                        ani_matrix_obj,
                        genus_map=genus_map,
                    )
                    genus_boundary_ani = genus_boundary_result.genus_boundary
                except Exception as e_gb:
                    logger.debug("Genus boundary detection skipped: %s", e_gb)

                try:
                    nbr_analyzer = PhylogeneticNeighborhoodAnalyzer(
                        ani_matrix=ani_matrix_obj,
                        genus_map=genus_map,
                        genus_boundary=genus_boundary_ani,
                    )
                    clusters = nbr_analyzer.analyze(clusters)
                except Exception as e_nbr:
                    logger.warning(
                        "Neighborhood analysis failed: %s", e_nbr
                    )

        except Exception as e:
            logger.warning(f"Could not analyze novel diversity: {e}")
            content = EMPTY_SECTION_TEMPLATE.format(
                message="Could not analyze novel diversity. Check classification data."
            )
            return TAB_SECTION_TEMPLATE.format(
                tab_id="novel-diversity",
                active_class="",
                section_title="Novel Diversity",
                content=content,
            )

        if not clusters:
            content = NOVEL_EMPTY_TEMPLATE
            return TAB_SECTION_TEMPLATE.format(
                tab_id="novel-diversity",
                active_class="",
                section_title="Novel Diversity",
                content=content,
            )

        content_parts.append(
            self._build_novel_kpi_strip(clusters, genus_boundary_result)
        )

        # Build summary section
        content_parts.append(build_novel_summary_html(summary))

        # Build scatter plot of cluster quality
        scatter_fig = build_cluster_scatter_figure(
            clusters, width=1000, height=500, genus_boundary=genus_boundary_ani,
        )
        scatter_id = "plot-novel-scatter"
        self._register_plot(scatter_id, scatter_fig)

        content_parts.append(PLOT_CONTAINER_TEMPLATE.format(
            extra_class="full-width",
            title="Cluster Quality Overview",
            description=(
                "Scatter plot showing cluster novelty vs quality metrics. "
                "Marker size indicates read count. Colors indicate confidence rating."
            ),
            plot_id=scatter_id,
        ))

        # Build sunburst chart for phylogenetic context
        sunburst_fig = build_sunburst_figure(clusters, width=600, height=600)
        sunburst_id = "plot-novel-sunburst"
        self._register_plot(sunburst_id, sunburst_fig)

        content_parts.append(PLOT_CONTAINER_TEMPLATE.format(
            extra_class="full-width",
            title="Phylogenetic Context",
            description=(
                "Sunburst chart showing novel clusters within taxonomic hierarchy. "
                "Click to explore family/genus/cluster relationships."
            ),
            plot_id=sunburst_id,
        ))

        heatmap_panel = self._build_novel_heatmap_panel(clusters)
        if heatmap_panel:
            content_parts.append(heatmap_panel)

        # Build cluster table (includes phylogenetic placement)
        content_parts.append(build_cluster_table_html(clusters))

        content = "\n".join(content_parts)

        return TAB_SECTION_TEMPLATE.format(
            tab_id="novel-diversity",
            active_class="",
            section_title="Novel Diversity",
            content=content,
        )

    def _build_novel_kpi_strip(
        self, clusters: list[NovelCluster], genus_boundary_result: AdaptiveGenusThreshold | None
    ) -> str:
        """Build the KPI metric strip for the novel-diversity tab."""
        # Build KPI metric strip for novel diversity
        novel_kpi_cards = []
        cluster_count = len(clusters) if clusters else 0
        novel_kpi_cards.append(KPI_CARD_TEMPLATE.format(
            color_class="accent", value=str(cluster_count), label="Clusters",
        ))
        novel_kpi_cards.append(KPI_CARD_TEMPLATE.format(
            color_class="novel", value=str(self.summary.novel_species), label="Novel Species",
        ))
        novel_kpi_cards.append(KPI_CARD_TEMPLATE.format(
            color_class="novel", value=str(self.summary.novel_genus), label="Novel Genus",
        ))
        if self.summary.has_bayesian:
            high_conf_novel = len(self.df.filter(
                (pl.col("taxonomic_call").is_in(["Novel Species", "Novel Genus"])) &
                (pl.col("posterior_entropy") < 1.0)
            ))
            novel_kpi_cards.append(KPI_CARD_TEMPLATE.format(
                color_class="known", value=str(high_conf_novel), label="High Confidence",
            ))
        if genus_boundary_result is not None and genus_boundary_result.method != "fallback":
            novel_kpi_cards.append(KPI_CARD_TEMPLATE.format(
                color_class="accent",
                value=f"{genus_boundary_result.genus_boundary:.0f}%",
                label="Genus Boundary",
            ))
        return KPI_STRIP_TEMPLATE.format(cards="\n".join(novel_kpi_cards))

    def _build_novel_heatmap_panel(self, clusters: list[NovelCluster]) -> str:
        """Build the collapsible phylogenetic-context heatmap panel HTML."""
        # Build phylogenetic context heatmap(s) based on available matrices
        # Prefer AAI for protein mode, ANI for nucleotide mode
        # Show both when both are available (ANI for species, AAI for genus context)
        is_protein_mode = self.config.alignment_mode == "protein"

        # Determine which matrices to use
        heatmaps_to_build: list[tuple[pl.DataFrame, str]] = []

        if is_protein_mode:
            # Protein mode: prefer AAI, fallback to ANI
            if self.aai_matrix is not None and len(self.aai_matrix) > 0:
                heatmaps_to_build.append((self.aai_matrix, "AAI"))
            elif self.ani_matrix is not None and len(self.ani_matrix) > 0:
                heatmaps_to_build.append((self.ani_matrix, "ANI"))
        else:
            # Nucleotide mode: use ANI primarily
            if self.ani_matrix is not None and len(self.ani_matrix) > 0:
                heatmaps_to_build.append((self.ani_matrix, "ANI"))
            # Also show AAI if available (better for Novel Genus clusters)
            if self.aai_matrix is not None and len(self.aai_matrix) > 0:
                heatmaps_to_build.append((self.aai_matrix, "AAI"))

        heatmap_parts: list[str] = []
        for matrix, sim_type in heatmaps_to_build:
            try:
                # Get genome accessions from matrix
                matrix_cols = matrix.columns
                if "genome" in matrix_cols:
                    genome_accessions = matrix["genome"].to_list()
                else:
                    genome_accessions = list(matrix_cols)

                genome_labels_map = self._get_genome_labels_map(genome_accessions)

                # Build the phylogenetic context heatmap
                phylo_fig, phylo_metadata, _clustering_ok = build_phylogenetic_context_heatmap(
                    similarity_matrix=matrix,
                    novel_clusters=clusters,
                    genome_labels_map=genome_labels_map,
                    similarity_type=sim_type,
                    max_references=self.config.max_phylo_references,
                    max_clusters=self.config.max_phylo_clusters,
                )

                if phylo_metadata.get("n_total", 0) > 0:
                    # Add intro section with note if clusters were truncated
                    n_clusters_shown = phylo_metadata.get("n_clusters", 0)
                    total_clusters = len(clusters)
                    clusters_note = ""
                    if n_clusters_shown < total_clusters:
                        clusters_note = f" (top {n_clusters_shown} of {total_clusters} by read count)"

                    # Customize description based on similarity type
                    if sim_type == "AAI":
                        pass
                    else:
                        pass

                    heatmap_parts.append(PHYLOGENETIC_HEATMAP_INTRO_TEMPLATE.format(
                        n_references=phylo_metadata.get("n_references", 0),
                        n_clusters=n_clusters_shown,
                        clusters_note=clusters_note,
                    ))

                    # Register and add the heatmap
                    phylo_heatmap_id = f"plot-novel-phylo-heatmap-{sim_type.lower()}"
                    self._register_plot(phylo_heatmap_id, phylo_fig)

                    heatmap_parts.append(PLOT_CONTAINER_TEMPLATE.format(
                        extra_class="full-width",
                        title=f"Phylogenetic Context Heatmap ({sim_type})",
                        description=(
                            f"Extended {sim_type} heatmap showing novel clusters alongside reference genomes. "
                            f"Entries marked with [*] are novel clusters. Hover for details."
                        ),
                        plot_id=phylo_heatmap_id,
                    ))

                    # Add legend (only for first heatmap to avoid repetition)
                    if matrix is heatmaps_to_build[0][0]:
                        heatmap_parts.append(PHYLOGENETIC_HEATMAP_LEGEND_TEMPLATE)

            except Exception as e:
                logger.warning(f"Could not build {sim_type} phylogenetic context heatmap: {e}")

        # Wrap heatmaps in a collapsible panel if any were built
        if not heatmap_parts:
            return ""
        heatmap_html = "\n".join(heatmap_parts)
        return COLLAPSIBLE_PANEL_TEMPLATE.format(
            default_state="",
            panel_id="novel-heatmap",
            title="Phylogenetic Context Heatmap",
            content=heatmap_html,
        )
