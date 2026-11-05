"""ReportGenerator reference section."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import polars as pl

from metadarkmatter.visualization.plots.base import TAXONOMY_COLORS, place_legend_below
from metadarkmatter.visualization.report.breakdowns import (
    MAX_EXPORT_MATRIX_GENOMES,
    build_species_table_html,
    compute_species_breakdown,
    matrix_csv,
    species_breakdown_csv,
)
from metadarkmatter.visualization.report.components import (
    build_aai_heatmap,
    build_aai_stats_cards,
    build_ani_heatmap,
    build_ani_stats_cards,
)
from metadarkmatter.visualization.report.embed import csv_export_button, csv_export_payload
from metadarkmatter.visualization.report.sections._base import _ReportBase
from metadarkmatter.visualization.report.templates import (
    COLLAPSIBLE_PANEL_TEMPLATE,
    EMPTY_SECTION_TEMPLATE,
    PLOT_CONTAINER_SIMPLE_TEMPLATE,
    PLOT_CONTAINER_TEMPLATE,
    TAB_SECTION_TEMPLATE,
    TWO_COLUMN_ROW_TEMPLATE,
)

if TYPE_CHECKING:

    pass


logger = logging.getLogger(__name__)

# Species (or genomes) drawn in the stacked bar chart.
TOP_SPECIES_BAR_COUNT = 20


class ReferenceMixin(_ReportBase):
    """Reference report-section builders."""

    def _build_reference_section(self) -> str:
        """Build the Reference tab with species bar, heatmaps, phylogeny, and recruitment.

        Merges content from the former Species & Genomes, ANI, AAI,
        Phylogeny, and Recruitment tabs into a single unified section.

        Returns:
            HTML string for the Reference tab section.
        """
        content_parts: list[str] = []

        # 1. Species bar chart (compact, full-width)
        species_bar = self._build_species_bar_chart()
        if species_bar:
            content_parts.append(species_bar)
            content_parts.append(self._build_species_breakdown_table())

        # 2. ANI/AAI heatmaps
        ani_content = self._build_ani_content()
        aai_content = self._build_aai_content()

        if ani_content and aai_content:
            # Side-by-side layout
            heatmap_row = TWO_COLUMN_ROW_TEMPLATE.format(
                left_flex="1",
                right_flex="1",
                left_content=ani_content,
                right_content=aai_content,
            )
            content_parts.append(heatmap_row)
        elif ani_content:
            content_parts.append(ani_content)
        elif aai_content:
            content_parts.append(aai_content)

        # 3. Phylogenetic tree (collapsible, default collapsed)
        phylo_html = self._build_phylogeny_content()
        if phylo_html:
            content_parts.append(COLLAPSIBLE_PANEL_TEMPLATE.format(
                default_state="",
                panel_id="phylogeny",
                title="Phylogenetic Tree",
                content=phylo_html,
            ))

        # 4. Recruitment plots (collapsible, default collapsed, only if data)
        recruit_html = self._build_recruitment_content()
        if recruit_html:
            content_parts.append(COLLAPSIBLE_PANEL_TEMPLATE.format(
                default_state="",
                panel_id="recruitment",
                title="Recruitment Plots",
                content=recruit_html,
            ))

        if not content_parts:
            content_parts.append(EMPTY_SECTION_TEMPLATE.format(
                message="No reference data provided. Supply ANI matrix and/or genome metadata."
            ))

        return TAB_SECTION_TEMPLATE.format(
            tab_id="reference",
            active_class="",
            section_title="Reference",
            content="\n".join(content_parts),
        )

    def _build_ani_content(self) -> str | None:
        """Build ANI heatmap HTML content without tab wrapper.

        Returns:
            HTML string with ANI heatmap and statistics cards, or None if
            no ANI matrix is available.
        """
        if self.ani_matrix is None or len(self.ani_matrix) == 0:
            return None

        # Create genome labels map
        matrix_cols = self.ani_matrix.columns
        if "genome" in matrix_cols:
            genome_accessions = self.ani_matrix["genome"].to_list()
        else:
            genome_accessions = list(matrix_cols)

        genome_labels_map = self._get_genome_labels_map(genome_accessions)

        # Build full heatmap
        heatmap_fig, stats, _clustering_succeeded = build_ani_heatmap(
            self.ani_matrix,
            genome_labels_map,
        )

        ani_id = "plot-ani-heatmap"
        self._register_plot(ani_id, heatmap_fig)
        stats_html = build_ani_stats_cards(stats) + self._matrix_export_html(
            "ani", self.ani_matrix
        )

        # Check for representative genome toggle
        rep_matrix = self._filter_matrix_to_representatives(self.ani_matrix)

        if rep_matrix is not None:
            # Build representative-only heatmap
            rep_cols = rep_matrix.columns
            if "genome" in rep_cols:
                rep_accessions = rep_matrix["genome"].to_list()
            else:
                rep_accessions = list(rep_cols)
            rep_labels = self._get_genome_labels_map(rep_accessions)
            rep_fig, rep_stats, _ = build_ani_heatmap(rep_matrix, rep_labels)
            rep_id = "plot-ani-heatmap-reps"
            self._register_plot(rep_id, rep_fig)
            rep_stats_html = build_ani_stats_cards(rep_stats)

            n_all = len(genome_accessions)
            n_reps = len(rep_accessions)

            return self._build_toggle_heatmap_content(
                metric="ani",
                title="Average Nucleotide Identity Matrix",
                description=(
                    "Heatmap showing pairwise ANI values between reference genomes, "
                    "hierarchically clustered by similarity. "
                    "Red = high ANI (similar), blue = low ANI (distant). "
                    "ANI ~95% indicates species boundary (Jain et al. 2018)."
                ),
                all_plot_id=ani_id,
                all_stats_html=stats_html,
                rep_plot_id=rep_id,
                rep_stats_html=rep_stats_html,
                n_all=n_all,
                n_reps=n_reps,
            )

        return stats_html + PLOT_CONTAINER_TEMPLATE.format(
            extra_class="full-width",
            title="Average Nucleotide Identity Matrix",
            description=(
                "Heatmap showing pairwise ANI values between reference genomes, "
                "hierarchically clustered by similarity. "
                "Red = high ANI (similar), blue = low ANI (distant). "
                "ANI ~95% indicates species boundary (Jain et al. 2018)."
            ),
            plot_id=ani_id,
        )

    def _build_aai_content(self) -> str | None:
        """Build AAI heatmap HTML content without tab wrapper.

        Returns:
            HTML string with AAI heatmap and statistics cards, or None if
            no AAI matrix is available.
        """
        if self.aai_matrix is None or len(self.aai_matrix) == 0:
            return None

        # Create genome labels map
        matrix_cols = self.aai_matrix.columns
        if "genome" in matrix_cols:
            genome_accessions = self.aai_matrix["genome"].to_list()
        else:
            genome_accessions = list(matrix_cols)

        genome_labels_map = self._get_genome_labels_map(genome_accessions)

        # Build full heatmap
        heatmap_fig, stats, _clustering_succeeded = build_aai_heatmap(
            self.aai_matrix,
            genome_labels_map,
        )

        aai_id = "plot-aai-heatmap"
        self._register_plot(aai_id, heatmap_fig)
        stats_html = build_aai_stats_cards(stats) + self._matrix_export_html(
            "aai", self.aai_matrix
        )

        # Check for representative genome toggle
        rep_matrix = self._filter_matrix_to_representatives(self.aai_matrix)

        if rep_matrix is not None:
            rep_cols = rep_matrix.columns
            if "genome" in rep_cols:
                rep_accessions = rep_matrix["genome"].to_list()
            else:
                rep_accessions = list(rep_cols)
            rep_labels = self._get_genome_labels_map(rep_accessions)
            rep_fig, rep_stats, _ = build_aai_heatmap(rep_matrix, rep_labels)
            rep_id = "plot-aai-heatmap-reps"
            self._register_plot(rep_id, rep_fig)
            rep_stats_html = build_aai_stats_cards(rep_stats)

            n_all = len(genome_accessions)
            n_reps = len(rep_accessions)

            return self._build_toggle_heatmap_content(
                metric="aai",
                title="Average Amino Acid Identity Matrix",
                description=(
                    "Heatmap showing pairwise AAI values between reference genomes, "
                    "hierarchically clustered by similarity. "
                    "Red = high AAI (similar), blue = low AAI (distant). "
                    "AAI ~65% indicates genus boundary (Riesco & Trujillo 2024)."
                ),
                all_plot_id=aai_id,
                all_stats_html=stats_html,
                rep_plot_id=rep_id,
                rep_stats_html=rep_stats_html,
                n_all=n_all,
                n_reps=n_reps,
            )

        return stats_html + PLOT_CONTAINER_TEMPLATE.format(
            extra_class="full-width",
            title="Average Amino Acid Identity Matrix",
            description=(
                "Heatmap showing pairwise AAI values between reference genomes, "
                "hierarchically clustered by similarity. "
                "Red = high AAI (similar), blue = low AAI (distant). "
                "AAI ~65% indicates genus boundary (Riesco & Trujillo 2024)."
            ),
            plot_id=aai_id,
        )

    def _build_toggle_heatmap_content(
        self,
        *,
        metric: str,
        title: str,
        description: str,
        all_plot_id: str,
        all_stats_html: str,
        rep_plot_id: str,
        rep_stats_html: str,
        n_all: int,
        n_reps: int,
    ) -> str:
        """Build HTML content with toggle between all genomes and representatives.

        Args:
            metric: Short identifier ('ani' or 'aai') used for element IDs
            title: Plot title text
            description: Plot description text
            all_plot_id: Plotly div ID for the all-genomes heatmap
            all_stats_html: Statistics cards HTML for all genomes
            rep_plot_id: Plotly div ID for the representative-only heatmap
            rep_stats_html: Statistics cards HTML for representatives
            n_all: Total number of genomes
            n_reps: Number of representative genomes

        Returns:
            HTML string with toggle controls and both heatmaps
        """
        return f"""
        <div class="heatmap-toggle-controls">
            <button class="heatmap-toggle-btn active"
                    id="{metric}-toggle-all"
                    onclick="toggleHeatmapView('{metric}', 'all')">
                All Genomes ({n_all})
            </button>
            <button class="heatmap-toggle-btn"
                    id="{metric}-toggle-reps"
                    onclick="toggleHeatmapView('{metric}', 'reps')">
                Representatives Only ({n_reps})
            </button>
        </div>
        <div id="{metric}-view-all">
            {all_stats_html}
            <div class="plot-container full-width">
                <div class="plot-title">{title}</div>
                <div class="plot-description">{description}</div>
                <div id="{all_plot_id}" class="plotly-chart"></div>
            </div>
        </div>
        <div id="{metric}-view-reps" style="display:none;">
            {rep_stats_html}
            <div class="plot-container full-width">
                <div class="plot-title">{title} (Representatives Only)</div>
                <div class="plot-description">
                    Showing {n_reps} species representative genomes
                    (one per species from GTDB). {description}
                </div>
                <div id="{rep_plot_id}" class="plotly-chart"></div>
            </div>
        </div>
        """

    def _filter_matrix_to_representatives(
        self, matrix: pl.DataFrame
    ) -> pl.DataFrame | None:
        """Filter a similarity matrix to representative genomes only.

        Args:
            matrix: Similarity matrix DataFrame (wide format, optional 'genome' column)

        Returns:
            Filtered matrix with only representative genomes, or None if
            metadata lacks representative information or filtering would
            remove all genomes.
        """
        if self.genome_metadata is None or not self.genome_metadata.has_representatives:
            return None

        rep_set = set(self.genome_metadata.build_representative_mapping().values())

        # Identify accessions in the matrix
        cols = matrix.columns
        has_genome_col = "genome" in cols
        accessions = matrix["genome"].to_list() if has_genome_col else list(cols)

        # Filter to representatives present in the matrix
        keep = [acc for acc in accessions if acc in rep_set]
        if len(keep) < 2:
            return None  # Need at least 2 genomes for a meaningful heatmap

        if has_genome_col:
            filtered = matrix.filter(pl.col("genome").is_in(keep))
            filtered = filtered.select(["genome", *keep])
        else:
            # Columns are accessions; filter both rows and columns
            indices = [i for i, acc in enumerate(accessions) if acc in keep]
            filtered = matrix.select([accessions[i] for i in indices])
            filtered = filtered[indices]

        return filtered

    def _build_species_bar_chart(self) -> str:
        """Stacked horizontal bars: top species (or genomes) by category.

        Each bar is one species (the genome accession when no species
        metadata is available), split by Bayesian ``taxonomic_call`` so the
        chart shows both how many reads a species received and how novel
        those reads were. Off-target reads are excluded.

        Returns:
            HTML for the chart container, or "" when there is nothing to plot.
        """
        import plotly.graph_objects as go

        table, categories, uses_species = compute_species_breakdown(
            self.df, self.genome_metadata
        )
        if len(table) == 0:
            return ""
        top = table.head(TOP_SPECIES_BAR_COUNT)
        labels = top["label"].to_list()[::-1]

        fig = go.Figure()
        for cat in categories:
            fig.add_trace(go.Bar(
                y=labels,
                x=top[cat].to_list()[::-1],
                name=cat,
                orientation="h",
                marker={"color": TAXONOMY_COLORS.get(cat, "#999999"),
                        "line": {"color": "rgba(0,0,0,0.35)", "width": 0.5}},
                hovertemplate=f"<b>%{{y}}</b><br>{cat}: %{{x:,}} reads<extra></extra>",
            ))
        fig.add_trace(go.Scatter(
            y=labels,
            x=top["reads"].to_list()[::-1],
            mode="text",
            text=[f"{n:,}" for n in top["reads"].to_list()[::-1]],
            textposition="middle right",
            showlegend=False,
            hoverinfo="skip",
            cliponaxis=False,
        ))
        noun = "Species" if uses_species else "Genomes"
        fig.update_layout(
            title=f"Top {len(top)} {noun} by Read Count, Stacked by Classification",
            xaxis_title="Number of Reads",
            yaxis_title="",
            barmode="stack",
            template="plotly_white",
            height=max(400, 30 * len(top) + 140),
            margin={"l": 40, "r": 60, "t": 60, "b": 60},
        )
        fig.update_yaxes(automargin=True)
        place_legend_below(fig)
        species_bar_id = "plot-species-bar"
        self._register_plot(species_bar_id, fig)

        return PLOT_CONTAINER_SIMPLE_TEMPLATE.format(
            extra_class="full-width",
            plot_id=species_bar_id,
        )

    def _build_species_breakdown_table(self) -> str:
        """Sortable per-species (or per-genome) table with a CSV download."""
        table, categories, uses_species = compute_species_breakdown(
            self.df, self.genome_metadata
        )
        if len(table) == 0:
            return ""
        noun = "species" if uses_species else "genome"
        header, rows = species_breakdown_csv(table, categories, uses_species=uses_species)
        export = csv_export_button(
            "species-breakdown",
            csv_export_payload(f"{noun}_breakdown.csv", header, rows),
            label="Download CSV",
        )
        return (
            '<div class="plot-container full-width" id="species-breakdown">'
            f'<div class="plot-title">Breakdown by {noun}</div>'
            '<div class="plot-description">On-target reads per '
            f"{noun} (best-match genome&rsquo;s species when metadata is available, "
            "otherwise the genome accession), split by classification. "
            "Novel % is Novel Species plus Novel Genus reads over all reads of the row. "
            "Click a column header to sort.</div>"
            f'<div class="export-bar">{export}</div>'
            + build_species_table_html(table, categories, uses_species=uses_species)
            + "</div>"
        )

    def _matrix_export_html(self, name: str, matrix: pl.DataFrame | None) -> str:
        """CSV download button for a similarity matrix (or a size note)."""
        if matrix is None or len(matrix) == 0:
            return ""
        exported = matrix_csv(matrix)
        if exported is None:
            return (
                '<p class="mdm-note">The matrix has more than '
                f"{MAX_EXPORT_MATRIX_GENOMES} genomes and is not embedded for download; "
                "use the matrix file passed to the report.</p>"
            )
        header, rows = exported
        return '<div class="export-bar">' + csv_export_button(
            f"{name}-matrix",
            csv_export_payload(f"{name}_matrix.csv", header, rows),
            label=f"Download {name.upper()} matrix (CSV)",
        ) + "</div>"

    def _build_recruitment_content(self) -> str | None:
        """Build recruitment plot HTML content without tab wrapper.

        Returns:
            HTML string with recruitment plots, or None if no recruitment
            data is available.
        """
        if self.recruitment_data is None or len(self.recruitment_data) == 0:
            return None

        # Import here to avoid circular imports
        from metadarkmatter.visualization.recruitment_plots import (
            RecruitmentPlotGenerator,
        )

        # Create recruitment plot
        gen = RecruitmentPlotGenerator(self.recruitment_data)
        fig = gen.create_figure(
            max_points=self.config.max_scatter_points,
            width=1100,
            height=600,
        )

        recruit_id = "plot-recruitment"
        self._register_plot(recruit_id, fig)

        # Multi-genome view
        multi_fig = gen.create_multi_genome_figure(
            max_points_per_genome=20000,
            width=1100,
        )
        multi_id = "plot-recruitment-multi"
        self._register_plot(multi_id, multi_fig)

        return (
            PLOT_CONTAINER_TEMPLATE.format(
                extra_class="full-width",
                title="Read Recruitment Overview",
                description=(
                    "Scatter plot showing percent identity vs genome position "
                    "for all mapped reads."
                ),
                plot_id=recruit_id,
            )
            + PLOT_CONTAINER_TEMPLATE.format(
                extra_class="full-width",
                title="Per-Genome Recruitment",
                description=(
                    "Individual recruitment plots for top reference genomes "
                    "by read count."
                ),
                plot_id=multi_id,
            )
        )
