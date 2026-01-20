"""
Report generator for unified HTML reports.

Combines multiple plot components into a single, self-contained HTML report
with tabbed navigation, interactive data tables, and professional styling.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

import polars as pl

logger = logging.getLogger(__name__)

from metadarkmatter.core.io_utils import read_dataframe
from metadarkmatter.core.metadata import GenomeMetadata
from metadarkmatter.visualization.plots.base import (
    PlotConfig,
    ThresholdConfig,
    format_count,
)
from metadarkmatter.visualization.report.components import (
    build_aai_heatmap,
    build_aai_stats_cards,
    build_ani_heatmap,
    build_ani_stats_cards,
)
from metadarkmatter.visualization.plots.classification_charts import (
    ClassificationBarChart,
    ClassificationDonutChart,
    DiversityDonutChart,
    DiversitySunburstChart,
)
from metadarkmatter.visualization.plots.distributions import (
    NoveltyHistogram,
    UncertaintyHistogram,
)
from metadarkmatter.visualization.plots.scatter_2d import (
    NoveltyUncertaintyScatter,
)
from metadarkmatter.visualization.report.styles import get_css_styles
from metadarkmatter.visualization.report.templates import (
    ANI_NOT_PROVIDED_MESSAGE,
    CATEGORY_BREAKDOWN_TEMPLATE,
    DATA_COLUMN_GUIDE_TEMPLATE,
    DATA_QUICK_FILTERS_TEMPLATE,
    DATA_SUMMARY_TEMPLATE,
    DATA_TABLE_JS,
    DATA_TABLE_TEMPLATE,
    DISTRIBUTIONS_SUMMARY_TEMPLATE,
    DIVERSITY_SUMMARY_TEMPLATE,
    EMPTY_SECTION_TEMPLATE,
    GENOME_HIGHLIGHTS_TEMPLATE,
    GENOME_INTERPRETATION_TEMPLATE,
    GENOMES_SUMMARY_TEMPLATE,
    METRIC_CARD_TEMPLATE,
    METRIC_CARDS_CONTAINER,
    PLOT_CONTAINER_SIMPLE_TEMPLATE,
    PLOT_CONTAINER_TEMPLATE,
    PLOT_ROW_TEMPLATE,
    RECRUITMENT_NOT_PROVIDED_MESSAGE,
    REPORT_BASE_TEMPLATE,
    SCATTER_INTERPRETATION_TEMPLATE,
    TAB_NAVIGATION_JS,
    TAB_SECTION_TEMPLATE,
    TABLE_ROW_TEMPLATE,
    get_cell_class,
)

if TYPE_CHECKING:
    import plotly.graph_objects as go


# Version for report footer
try:
    from metadarkmatter import __version__
except ImportError:
    __version__ = "0.1.0"


@dataclass
class TaxonomicSummary:
    """Summary statistics from classification results."""

    total_reads: int = 0
    known_species: int = 0
    novel_species: int = 0
    novel_genus: int = 0
    species_boundary: int = 0
    ambiguous: int = 0
    ambiguous_within_genus: int = 0
    conserved_regions: int = 0
    unclassified: int = 0
    # High-level diversity grouping
    diversity_known: int = 0
    diversity_novel: int = 0
    diversity_uncertain: int = 0
    mean_novelty_index: float = 0.0
    mean_placement_uncertainty: float = 0.0
    mean_top_hit_identity: float = 0.0

    @property
    def novel_percentage(self) -> float:
        """Percentage of reads classified as novel (species or genus)."""
        if self.total_reads == 0:
            return 0.0
        return (self.novel_species + self.novel_genus) / self.total_reads * 100

    @property
    def diversity_known_pct(self) -> float:
        """Percentage of reads with confident known diversity."""
        if self.total_reads == 0:
            return 0.0
        return self.diversity_known / self.total_reads * 100

    @property
    def diversity_novel_pct(self) -> float:
        """Percentage of reads with confident novel diversity."""
        if self.total_reads == 0:
            return 0.0
        return self.diversity_novel / self.total_reads * 100

    @property
    def diversity_uncertain_pct(self) -> float:
        """Percentage of reads with uncertain diversity status."""
        if self.total_reads == 0:
            return 0.0
        return self.diversity_uncertain / self.total_reads * 100

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for template rendering."""
        return {
            "total_reads": self.total_reads,
            "known_species": self.known_species,
            "novel_species": self.novel_species,
            "novel_genus": self.novel_genus,
            "species_boundary": self.species_boundary,
            "ambiguous": self.ambiguous,
            "ambiguous_within_genus": self.ambiguous_within_genus,
            "conserved_regions": self.conserved_regions,
            "unclassified": self.unclassified,
            "diversity_known": self.diversity_known,
            "diversity_novel": self.diversity_novel,
            "diversity_uncertain": self.diversity_uncertain,
            "mean_novelty_index": self.mean_novelty_index,
            "mean_placement_uncertainty": self.mean_placement_uncertainty,
            "mean_top_hit_identity": self.mean_top_hit_identity,
            "novel_percentage": self.novel_percentage,
            "diversity_known_pct": self.diversity_known_pct,
            "diversity_novel_pct": self.diversity_novel_pct,
            "diversity_uncertain_pct": self.diversity_uncertain_pct,
        }


@dataclass
class ReportConfig:
    """Configuration for report generation."""

    sample_name: str = "Sample"
    title: str = "Metadarkmatter Classification Report"
    theme: str = "light"
    page_size: int = 100
    max_table_rows: int = 10000
    max_scatter_points: int = 50000
    include_plotlyjs: str = "cdn"  # 'cdn', 'embed', or path

    plot_config: PlotConfig = field(default_factory=PlotConfig)
    thresholds: ThresholdConfig = field(default_factory=ThresholdConfig)


class ReportGenerator:
    """
    Generates unified HTML reports from classification results.

    Creates a self-contained HTML file with:
    - Overview metrics and charts
    - Distribution histograms
    - 2D scatter plot (novelty vs uncertainty)
    - Recruitment plots (if BAM provided)
    - Genome breakdown charts
    - ANI matrix heatmap (if provided)
    - Interactive data table
    """

    def __init__(
        self,
        classifications: pl.DataFrame,
        config: ReportConfig | None = None,
        recruitment_data: pl.DataFrame | None = None,
        ani_matrix: pl.DataFrame | None = None,
        aai_matrix: pl.DataFrame | None = None,
        genome_metadata: GenomeMetadata | None = None,
    ) -> None:
        """
        Initialize report generator.

        Args:
            classifications: DataFrame with classification results. Must include:
                - read_id: Read identifier
                - best_match_genome: Best matching genome
                - top_hit_identity: Percent identity to best hit
                - novelty_index: 100 - top_hit_identity
                - placement_uncertainty: ANI-based uncertainty
                - taxonomic_call: Classification category
            config: Report configuration
            recruitment_data: Optional recruitment plot data
            ani_matrix: Optional ANI matrix for heatmap
            aai_matrix: Optional AAI matrix for genus-level heatmap
            genome_metadata: Optional genome metadata for species-level aggregation
        """
        self.df = classifications
        self.config = config or ReportConfig()
        self.recruitment_data = recruitment_data
        self.ani_matrix = ani_matrix
        self.aai_matrix = aai_matrix
        self.genome_metadata = genome_metadata

        # Compute summary statistics
        self.summary = self._compute_summary()

        # Store generated plot HTML/JS
        self._plot_data: dict[str, tuple[str, str]] = {}  # {plot_id: (div_html, js)}

    def _compute_summary(self) -> TaxonomicSummary:
        """Compute summary statistics from classifications."""
        # Count by classification
        counts = (
            self.df.group_by("taxonomic_call")
            .len()
            .to_dict(as_series=False)
        )

        count_map = dict(zip(counts.get("taxonomic_call", []),
                            counts.get("len", []), strict=True))

        # Mean metrics
        mean_novelty = self.df["novelty_index"].mean() or 0.0
        mean_uncertainty = self.df["placement_uncertainty"].mean() or 0.0
        mean_identity = self.df["top_hit_identity"].mean() or 0.0

        # Compute diversity groupings
        known_species = count_map.get("Known Species", 0)
        novel_species = count_map.get("Novel Species", 0)
        novel_genus = count_map.get("Novel Genus", 0)
        species_boundary = count_map.get("Species Boundary", 0)
        ambiguous = count_map.get("Ambiguous", 0)
        ambiguous_within_genus = count_map.get("Ambiguous Within Genus", 0)
        conserved_regions = count_map.get("Conserved Region", 0)
        unclassified = count_map.get("Unclassified", 0)

        # High-level diversity grouping
        diversity_known = known_species
        diversity_novel = novel_species + novel_genus
        diversity_uncertain = (
            species_boundary + ambiguous + ambiguous_within_genus +
            conserved_regions + unclassified
        )

        return TaxonomicSummary(
            total_reads=len(self.df),
            known_species=known_species,
            novel_species=novel_species,
            novel_genus=novel_genus,
            species_boundary=species_boundary,
            ambiguous=ambiguous,
            ambiguous_within_genus=ambiguous_within_genus,
            conserved_regions=conserved_regions,
            unclassified=unclassified,
            diversity_known=diversity_known,
            diversity_novel=diversity_novel,
            diversity_uncertain=diversity_uncertain,
            mean_novelty_index=mean_novelty,
            mean_placement_uncertainty=mean_uncertainty,
            mean_top_hit_identity=mean_identity,
        )

    def _get_genome_label(self, accession: str, max_species_len: int = 25) -> str:
        """
        Get a display label for a genome accession.

        If metadata is available, returns "Accession (Species)" format.
        Otherwise returns just the accession.

        Args:
            accession: Genome accession (e.g., GCF_019852205.1)
            max_species_len: Maximum length for species name before truncation

        Returns:
            Formatted label string
        """
        if self.genome_metadata is None:
            return accession

        # Look up species from metadata
        species = self.genome_metadata.get_species(accession)
        if species and species != "Unknown":
            # Truncate species name if too long
            if len(species) > max_species_len:
                species = species[:max_species_len-3] + "..."
            return f"{accession} ({species})"
        return accession

    def _get_genome_labels_map(self, accessions: list[str]) -> dict[str, str]:
        """
        Create a mapping from accessions to display labels.

        Args:
            accessions: List of genome accessions

        Returns:
            Dictionary mapping accession -> display label
        """
        return {acc: self._get_genome_label(acc) for acc in accessions}

    def generate(self, output_path: Path | str) -> None:
        """
        Generate and save the HTML report.

        Args:
            output_path: Path for output HTML file
        """
        html = self._build_html()

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(html, encoding="utf-8")

    def _build_html(self) -> str:
        """Build the complete HTML report."""
        # Generate all sections
        content_sections = []

        # Overview tab (always visible first)
        content_sections.append(self._build_overview_section())

        # Distributions tab
        content_sections.append(self._build_distributions_section())

        # Recruitment tab
        content_sections.append(self._build_recruitment_section())

        # Species tab (if metadata provided)
        content_sections.append(self._build_species_section())

        # Genomes tab
        content_sections.append(self._build_genomes_section())

        # ANI tab
        content_sections.append(self._build_ani_section())

        # AAI tab
        content_sections.append(self._build_aai_section())

        # Data table tab
        content_sections.append(self._build_data_section())

        content = "\n".join(content_sections)

        # Collect all Plotly JS initialization code
        plotly_js = self._build_plotly_js()

        # Combine JavaScript
        js_scripts = TAB_NAVIGATION_JS + "\n" + DATA_TABLE_JS.format(
            page_size=self.config.page_size
        )

        # Build final HTML
        return REPORT_BASE_TEMPLATE.format(
            title=self.config.title,
            sample_name=self.config.sample_name,
            timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            total_reads=format_count(self.summary.total_reads),
            version=__version__,
            css_styles=get_css_styles(self.config.theme),
            content=content,
            plotly_js=plotly_js,
            js_scripts=js_scripts,
        )

    def _build_overview_section(self) -> str:
        """Build the overview tab section with hero diversity summary."""
        s = self.summary

        # Calculate percentages for all categories
        total = s.total_reads if s.total_reads > 0 else 1  # Avoid division by zero

        # Diversity summary - the hero section
        diversity_html = DIVERSITY_SUMMARY_TEMPLATE.format(
            known_pct=s.diversity_known_pct,
            known_count=s.diversity_known,
            novel_pct=s.diversity_novel_pct,
            novel_count=s.diversity_novel,
            uncertain_pct=s.diversity_uncertain_pct,
            uncertain_count=s.diversity_uncertain,
        )

        # Category breakdown with detailed counts
        category_html = CATEGORY_BREAKDOWN_TEMPLATE.format(
            known_species=s.known_species,
            known_species_pct=s.known_species / total * 100,
            novel_species=s.novel_species,
            novel_species_pct=s.novel_species / total * 100,
            novel_genus=s.novel_genus,
            novel_genus_pct=s.novel_genus / total * 100,
            species_boundary=s.species_boundary,
            species_boundary_pct=s.species_boundary / total * 100,
            ambiguous=s.ambiguous,
            ambiguous_pct=s.ambiguous / total * 100,
            conserved_regions=s.conserved_regions,
            conserved_pct=s.conserved_regions / total * 100,
        )

        # Sunburst chart for interactive exploration - larger size
        sunburst = DiversitySunburstChart(
            s.to_dict(),
            config=PlotConfig(width=800, height=650),
            title="Diversity Classification",
        )
        sunburst_id = "plot-diversity-sunburst"
        self._register_plot(sunburst_id, sunburst.create_figure())

        # Build plot row with sunburst chart
        plots_row = PLOT_ROW_TEMPLATE.format(
            plots=(
                PLOT_CONTAINER_TEMPLATE.format(
                    extra_class="full-width",
                    title="Interactive Diversity Chart",
                    description="Click segments to explore. Inner ring shows diversity status (Known/Novel/Uncertain), outer ring shows detailed categories.",
                    plot_id=sunburst_id,
                )
            )
        )

        content = diversity_html + category_html + plots_row

        return TAB_SECTION_TEMPLATE.format(
            tab_id="overview",
            active_class="active",
            section_title="Overview",
            content=content,
        )

    def _build_metric_cards(self) -> str:
        """Build the metric cards HTML."""
        s = self.summary

        cards_html = ""

        # Total reads
        cards_html += METRIC_CARD_TEMPLATE.format(
            value=format_count(s.total_reads),
            label="Total Reads",
            subtext="Classified reads",
            color_class="",
        )

        # Known species
        known_pct = (s.known_species / s.total_reads * 100) if s.total_reads else 0
        cards_html += METRIC_CARD_TEMPLATE.format(
            value=format_count(s.known_species),
            label="Known Species",
            subtext=f"{known_pct:.1f}% of reads",
            color_class="success",
        )

        # Novel species
        novel_sp_pct = (s.novel_species / s.total_reads * 100) if s.total_reads else 0
        cards_html += METRIC_CARD_TEMPLATE.format(
            value=format_count(s.novel_species),
            label="Novel Species",
            subtext=f"{novel_sp_pct:.1f}% of reads",
            color_class="warning",
        )

        # Novel genus
        novel_gen_pct = (s.novel_genus / s.total_reads * 100) if s.total_reads else 0
        cards_html += METRIC_CARD_TEMPLATE.format(
            value=format_count(s.novel_genus),
            label="Novel Genus",
            subtext=f"{novel_gen_pct:.1f}% of reads",
            color_class="danger",
        )

        # Ambiguous
        ambig_pct = (s.ambiguous / s.total_reads * 100) if s.total_reads else 0
        cards_html += METRIC_CARD_TEMPLATE.format(
            value=format_count(s.ambiguous),
            label="Ambiguous",
            subtext=f"{ambig_pct:.1f}% of reads",
            color_class="muted",
        )

        # Novel diversity (combined)
        cards_html += METRIC_CARD_TEMPLATE.format(
            value=f"{s.novel_percentage:.1f}%",
            label="Novel Diversity",
            subtext="Novel species + genus",
            color_class="info",
        )

        # Mean novelty index
        cards_html += METRIC_CARD_TEMPLATE.format(
            value=f"{s.mean_novelty_index:.2f}",
            label="Mean Novelty",
            subtext="Average novelty index",
            color_class="",
        )

        return METRIC_CARDS_CONTAINER.format(cards=cards_html)

    def _build_distributions_section(self) -> str:
        """Build the distributions tab section with improved UX."""
        s = self.summary

        # Calculate interpretations based on mean values
        if s.mean_novelty_index < 5:
            novelty_interp = "Mostly known taxa"
        elif s.mean_novelty_index < 15:
            novelty_interp = "Some novel diversity"
        else:
            novelty_interp = "High novel diversity"

        if s.mean_placement_uncertainty < 2:
            uncertainty_interp = "Confident placements"
        elif s.mean_placement_uncertainty < 5:
            uncertainty_interp = "Moderate confidence"
        else:
            uncertainty_interp = "Many ambiguous hits"

        # Count confident classifications (low uncertainty)
        confident_count = len(self.df.filter(
            self.df["placement_uncertainty"] < 2
        ))
        confident_pct = (confident_count / s.total_reads * 100) if s.total_reads else 0

        # Summary section
        summary_html = DISTRIBUTIONS_SUMMARY_TEMPLATE.format(
            mean_identity=s.mean_top_hit_identity,
            mean_novelty=s.mean_novelty_index,
            novelty_interpretation=novelty_interp,
            mean_uncertainty=s.mean_placement_uncertainty,
            uncertainty_interpretation=uncertainty_interp,
            confident_pct=confident_pct,
        )

        # 2D scatter (the key diagnostic plot) - show first as it's the primary view
        scatter = NoveltyUncertaintyScatter(
            self.df,
            config=PlotConfig(width=1100, height=650),
            thresholds=self.config.thresholds,
            max_points=self.config.max_scatter_points,
        )
        scatter_id = "plot-scatter-2d"
        self._register_plot(scatter_id, scatter.create_figure())

        # Interpretation guide
        interpretation_html = SCATTER_INTERPRETATION_TEMPLATE

        scatter_plot = PLOT_CONTAINER_TEMPLATE.format(
            extra_class="full-width",
            title="Novelty vs Uncertainty Landscape",
            description=(
                "Each point represents a read. Position indicates divergence from references (x-axis) "
                "and confidence in placement (y-axis). Colors show classification categories."
            ),
            plot_id=scatter_id,
        )

        # Novelty histogram
        novelty_hist = NoveltyHistogram(
            self.df,
            config=PlotConfig(width=550, height=400),
            thresholds=self.config.thresholds,
        )
        novelty_id = "plot-novelty-hist"
        self._register_plot(novelty_id, novelty_hist.create_figure())

        # Uncertainty histogram
        uncertainty_hist = UncertaintyHistogram(
            self.df,
            config=PlotConfig(width=550, height=400),
            thresholds=self.config.thresholds,
        )
        uncertainty_id = "plot-uncertainty-hist"
        self._register_plot(uncertainty_id, uncertainty_hist.create_figure())

        # Build histogram row
        hist_row = PLOT_ROW_TEMPLATE.format(
            plots=(
                PLOT_CONTAINER_TEMPLATE.format(
                    extra_class="half-width",
                    title="Novelty Index Distribution",
                    description=(
                        "Low values (0-5) = known species. "
                        "Medium (5-20) = novel species. "
                        "High (>20) = novel genus or more divergent."
                    ),
                    plot_id=novelty_id,
                )
                + PLOT_CONTAINER_TEMPLATE.format(
                    extra_class="half-width",
                    title="Placement Uncertainty Distribution",
                    description=(
                        "Low values (<2) = confident placement. "
                        "Medium (2-5) = species boundary zone. "
                        "High (>5) = conserved regions or ambiguous."
                    ),
                    plot_id=uncertainty_id,
                )
            )
        )

        content = summary_html + interpretation_html + scatter_plot + hist_row

        return TAB_SECTION_TEMPLATE.format(
            tab_id="distributions",
            active_class="",
            section_title="Distributions",
            content=content,
        )

    def _build_recruitment_section(self) -> str:
        """Build the recruitment plots tab section."""
        if self.recruitment_data is None or len(self.recruitment_data) == 0:
            content = EMPTY_SECTION_TEMPLATE.format(
                message=RECRUITMENT_NOT_PROVIDED_MESSAGE
            )
        else:
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

            content = (
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

        return TAB_SECTION_TEMPLATE.format(
            tab_id="recruitment",
            active_class="",
            section_title="Recruitment Plots",
            content=content,
        )

    def _build_species_section(self) -> str:
        """Build the species breakdown tab section."""
        if self.genome_metadata is None:
            content = EMPTY_SECTION_TEMPLATE.format(
                message=(
                    "Species breakdown not available. "
                    "Use --metadata option with genome_metadata.tsv "
                    "to enable species-level aggregation."
                )
            )
        else:
            import plotly.graph_objects as go

            # Join metadata to get species information
            if "species" not in self.df.columns:
                # Need to join species information
                enriched_df = self.genome_metadata.join_classifications(self.df)
            else:
                enriched_df = self.df

            # Aggregate by species
            species_counts = (
                enriched_df.group_by("species")
                .agg([
                    pl.len().alias("count"),
                    pl.col("novelty_index").mean().alias("mean_novelty"),
                    pl.col("top_hit_identity").mean().alias("mean_identity"),
                    pl.col("best_match_genome").n_unique().alias("genome_count"),
                ])
                .sort("count", descending=True)
                .head(20)
            )

            if len(species_counts) == 0:
                content = EMPTY_SECTION_TEMPLATE.format(
                    message="No species data available."
                )
            else:
                # Horizontal bar chart of species read counts
                bar_fig = go.Figure()
                bar_fig.add_trace(
                    go.Bar(
                        y=species_counts["species"].to_list()[::-1],
                        x=species_counts["count"].to_list()[::-1],
                        orientation="h",
                        marker_color="#22c55e",
                        text=[f"{c:,}" for c in species_counts["count"].to_list()[::-1]],
                        textposition="outside",
                        hovertemplate=(
                            "<b>%{y}</b><br>"
                            "Reads: %{x:,}<extra></extra>"
                        ),
                    )
                )
                bar_fig.update_layout(
                    title="Top 20 Species by Read Count",
                    xaxis_title="Number of Reads",
                    yaxis_title="",
                    template="plotly_white",
                    height=max(400, 30 * len(species_counts)),
                    margin={"l": 300, "r": 40, "t": 60, "b": 60},
                )
                species_bar_id = "plot-species-bar"
                self._register_plot(species_bar_id, bar_fig)

                # Pie chart showing species proportion
                pie_fig = go.Figure()
                pie_fig.add_trace(
                    go.Pie(
                        labels=species_counts["species"].head(10).to_list(),
                        values=species_counts["count"].head(10).to_list(),
                        textinfo="percent+label",
                        textposition="inside",
                        hole=0.4,
                        hovertemplate=(
                            "<b>%{label}</b><br>"
                            "Reads: %{value:,}<br>"
                            "Proportion: %{percent}<extra></extra>"
                        ),
                    )
                )
                pie_fig.update_layout(
                    title="Top 10 Species Composition",
                    template="plotly_white",
                    height=500,
                    showlegend=True,
                    legend=dict(
                        orientation="h",
                        yanchor="bottom",
                        y=-0.3,
                        xanchor="center",
                        x=0.5,
                    ),
                )
                species_pie_id = "plot-species-pie"
                self._register_plot(species_pie_id, pie_fig)

                # Summary statistics
                total_species = species_counts["species"].n_unique()

                content = (
                    f'<div class="metric-cards">'
                    f'<div class="metric-card">'
                    f'<div class="metric-value">{total_species}</div>'
                    f'<div class="metric-label">Unique Species</div>'
                    f'<div class="metric-subtext">Detected in sample</div>'
                    f'</div>'
                    f'<div class="metric-card">'
                    f'<div class="metric-value">{self.genome_metadata.species_count}</div>'
                    f'<div class="metric-label">Reference Species</div>'
                    f'<div class="metric-subtext">In metadata database</div>'
                    f'</div>'
                    f'</div>'
                    + PLOT_ROW_TEMPLATE.format(
                        plots=(
                            PLOT_CONTAINER_TEMPLATE.format(
                                extra_class="half-width",
                                title="Species Composition",
                                description=(
                                    "Pie chart showing the proportion of reads "
                                    "assigned to each species."
                                ),
                                plot_id=species_pie_id,
                            )
                            + PLOT_CONTAINER_SIMPLE_TEMPLATE.format(
                                extra_class="half-width",
                                plot_id=species_bar_id,
                            )
                        )
                    )
                )

        return TAB_SECTION_TEMPLATE.format(
            tab_id="species",
            active_class="",
            section_title="Species Breakdown",
            content=content,
        )

    def _build_genomes_section(self) -> str:
        """Build the genomes breakdown tab section with improved UX."""
        # Per-genome read counts with statistics
        genome_counts = (
            self.df.group_by("best_match_genome")
            .agg([
                pl.len().alias("count"),
                pl.col("novelty_index").mean().alias("mean_novelty"),
                pl.col("top_hit_identity").mean().alias("mean_identity"),
            ])
            .sort("count", descending=True)
        )

        if len(genome_counts) == 0:
            content = EMPTY_SECTION_TEMPLATE.format(
                message="No genome data available."
            )
        else:
            import plotly.graph_objects as go

            # Calculate summary statistics
            total_genomes = len(genome_counts)
            ref_genomes = (
                self.genome_metadata.genome_count
                if self.genome_metadata else total_genomes
            )

            # Top genome stats
            top_genome_acc = genome_counts["best_match_genome"][0]
            top_genome_count = genome_counts["count"][0]
            top_genome_pct = (
                (top_genome_count / self.summary.total_reads * 100)
                if self.summary.total_reads else 0
            )
            top_genome_name = self._get_genome_label(top_genome_acc, max_species_len=25)

            # Top 3 genomes concentration
            top3_count = genome_counts["count"].head(3).sum()
            top3_pct = (
                (top3_count / self.summary.total_reads * 100)
                if self.summary.total_reads else 0
            )
            concentration_hint = (
                "Concentrated" if top3_pct > 70 else
                "Moderate spread" if top3_pct > 40 else
                "Widely distributed"
            )

            # Median reads per genome
            median_reads = int(genome_counts["count"].median() or 0)

            # Build summary section
            summary_html = GENOMES_SUMMARY_TEMPLATE.format(
                total_genomes=total_genomes,
                ref_genomes=ref_genomes,
                top_genome_pct=top_genome_pct,
                top_genome_name=top_genome_name[:30],
                top3_pct=top3_pct,
                concentration_hint=concentration_hint,
                median_reads=median_reads,
            )

            # Find key genomes for highlights
            # Helper to get species name
            def get_species(acc: str) -> str:
                if self.genome_metadata:
                    species = self.genome_metadata.get_species(acc)
                    if species and species != "Unknown":
                        return species
                return "Unknown species"

            # Most abundant (already have it)
            top_identity = genome_counts["mean_identity"][0]
            top_species = get_species(top_genome_acc)

            # Most novel (highest mean novelty with at least 5 reads)
            novel_df = genome_counts.filter(pl.col("count") >= 5)
            if len(novel_df) > 0:
                most_novel_idx = novel_df["mean_novelty"].arg_max()
                novel_genome_acc = novel_df["best_match_genome"][most_novel_idx]
                novel_reads = novel_df["count"][most_novel_idx]
                novel_novelty = novel_df["mean_novelty"][most_novel_idx]
            else:
                novel_genome_acc = top_genome_acc
                novel_reads = top_genome_count
                novel_novelty = genome_counts["mean_novelty"][0]
            novel_species = get_species(novel_genome_acc)

            # Most confident (highest mean identity with at least 5 reads)
            confident_df = genome_counts.filter(pl.col("count") >= 5)
            if len(confident_df) > 0:
                most_confident_idx = confident_df["mean_identity"].arg_max()
                confident_genome_acc = confident_df["best_match_genome"][
                    most_confident_idx
                ]
                confident_reads = confident_df["count"][most_confident_idx]
                confident_identity = confident_df["mean_identity"][most_confident_idx]
            else:
                confident_genome_acc = top_genome_acc
                confident_reads = top_genome_count
                confident_identity = top_identity
            confident_species = get_species(confident_genome_acc)

            # Build highlights section
            highlights_html = GENOME_HIGHLIGHTS_TEMPLATE.format(
                top_species=top_species,
                top_accession=top_genome_acc,
                top_reads=top_genome_count,
                top_identity=top_identity,
                novel_species=novel_species,
                novel_accession=novel_genome_acc,
                novel_reads=novel_reads,
                novel_novelty=novel_novelty,
                confident_species=confident_species,
                confident_accession=confident_genome_acc,
                confident_reads=confident_reads,
                confident_identity=confident_identity,
            )

            # Interpretation guide
            interpretation_html = GENOME_INTERPRETATION_TEMPLATE

            # Get top 20 genomes for plotting
            genome_counts_top20 = genome_counts.head(20)
            genome_accessions = genome_counts_top20["best_match_genome"].to_list()
            genome_labels = [
                self._get_genome_label(acc, max_species_len=30)
                for acc in genome_accessions
            ]

            # Horizontal bar chart of genome read counts
            bar_fig = go.Figure()
            bar_fig.add_trace(
                go.Bar(
                    y=genome_labels[::-1],
                    x=genome_counts_top20["count"].to_list()[::-1],
                    orientation="h",
                    marker_color="#667eea",
                    text=[
                        f"{c:,}"
                        for c in genome_counts_top20["count"].to_list()[::-1]
                    ],
                    textposition="outside",
                    hovertemplate=(
                        "<b>%{y}</b><br>"
                        "Reads: %{x:,}<extra></extra>"
                    ),
                )
            )
            bar_fig.update_layout(
                title="Top 20 Reference Genomes by Read Count",
                xaxis_title="Number of Reads",
                yaxis_title="",
                template="plotly_white",
                height=max(400, 35 * len(genome_counts_top20)),
                margin={"l": 350, "r": 40, "t": 60, "b": 60},
            )
            genome_bar_id = "plot-genome-bar"
            self._register_plot(genome_bar_id, bar_fig)

            # Box plot of identity distribution per genome
            # Get identity data for top 10 genomes
            top_genomes = genome_counts["best_match_genome"].head(10).to_list()
            box_data = self.df.filter(pl.col("best_match_genome").is_in(top_genomes))

            box_fig = go.Figure()
            for genome in top_genomes:
                gdata = box_data.filter(pl.col("best_match_genome") == genome)
                # Use species-annotated label for box plot
                label = self._get_genome_label(genome, max_species_len=20)
                box_fig.add_trace(
                    go.Box(
                        y=gdata["top_hit_identity"].to_list(),
                        name=label[:40] + ("..." if len(label) > 40 else ""),
                        boxpoints=False,
                        hoverinfo="y+name",
                    )
                )

            box_fig.update_layout(
                title="Identity Distribution by Genome (Top 10)",
                yaxis_title="Percent Identity (%)",
                template="plotly_white",
                showlegend=False,
                height=500,
            )
            genome_box_id = "plot-genome-box"
            self._register_plot(genome_box_id, box_fig)

            # Build content with new templates
            plots_html = (
                PLOT_CONTAINER_TEMPLATE.format(
                    extra_class="full-width",
                    title="Read Counts by Reference Genome",
                    description=(
                        "Bar chart showing the number of reads assigned "
                        "to each reference genome."
                    ),
                    plot_id=genome_bar_id,
                )
                + PLOT_CONTAINER_TEMPLATE.format(
                    extra_class="full-width",
                    title="Identity Distribution by Genome",
                    description=(
                        "Box plots showing the distribution of alignment identity "
                        "for top genomes."
                    ),
                    plot_id=genome_box_id,
                )
            )

            content = summary_html + highlights_html + interpretation_html + plots_html

        return TAB_SECTION_TEMPLATE.format(
            tab_id="genomes",
            active_class="",
            section_title="Genome Breakdown",
            content=content,
        )

    def _build_ani_section(self) -> str:
        """Build the ANI matrix heatmap tab section."""
        if self.ani_matrix is None or len(self.ani_matrix) == 0:
            content = EMPTY_SECTION_TEMPLATE.format(
                message=ANI_NOT_PROVIDED_MESSAGE
            )
        else:
            # Create genome labels map
            matrix_cols = self.ani_matrix.columns
            if "genome" in matrix_cols:
                genome_accessions = self.ani_matrix["genome"].to_list()
            else:
                genome_accessions = list(matrix_cols)

            genome_labels_map = self._get_genome_labels_map(genome_accessions)

            # Build heatmap and compute statistics using extracted components
            heatmap_fig, stats, clustering_succeeded = build_ani_heatmap(
                self.ani_matrix,
                genome_labels_map,
            )

            # Register plot
            ani_id = "plot-ani-heatmap"
            self._register_plot(ani_id, heatmap_fig)

            # Generate statistics cards
            stats_html = build_ani_stats_cards(stats)

            content = stats_html + PLOT_CONTAINER_TEMPLATE.format(
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

        return TAB_SECTION_TEMPLATE.format(
            tab_id="ani",
            active_class="",
            section_title="ANI Matrix",
            content=content,
        )

    def _build_aai_section(self) -> str:
        """Build the AAI matrix heatmap tab section for genus-level comparisons."""
        if self.aai_matrix is None or len(self.aai_matrix) == 0:
            content = EMPTY_SECTION_TEMPLATE.format(
                message=(
                    "AAI matrix not provided. "
                    "Use --aai option with AAI matrix CSV to enable genus-level heatmap. "
                    "Generate with: metadarkmatter aai compute"
                )
            )
        else:
            # Create genome labels map
            matrix_cols = self.aai_matrix.columns
            if "genome" in matrix_cols:
                genome_accessions = self.aai_matrix["genome"].to_list()
            else:
                genome_accessions = list(matrix_cols)

            genome_labels_map = self._get_genome_labels_map(genome_accessions)

            # Build heatmap and compute statistics using extracted components
            heatmap_fig, stats, clustering_succeeded = build_aai_heatmap(
                self.aai_matrix,
                genome_labels_map,
            )

            # Register plot
            aai_id = "plot-aai-heatmap"
            self._register_plot(aai_id, heatmap_fig)

            # Generate statistics cards
            stats_html = build_aai_stats_cards(stats)

            content = stats_html + PLOT_CONTAINER_TEMPLATE.format(
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

        return TAB_SECTION_TEMPLATE.format(
            tab_id="aai",
            active_class="",
            section_title="AAI Matrix",
            content=content,
        )

    def _build_data_section(self) -> str:
        """Build the interactive data table section with improved UX."""
        s = self.summary

        # Limit rows for performance
        table_df = self.df.head(self.config.max_table_rows)

        # Build summary section
        summary_html = DATA_SUMMARY_TEMPLATE.format(
            total_rows=len(table_df),
            known_count=s.diversity_known,
            novel_count=s.diversity_novel,
            uncertain_count=s.diversity_uncertain,
        )

        # Build quick filters section
        quick_filters_html = DATA_QUICK_FILTERS_TEMPLATE.format(
            known_count=s.known_species,
            novel_species_count=s.novel_species,
            novel_genus_count=s.novel_genus,
            ambiguous_count=s.ambiguous,
            conserved_count=s.conserved_regions,
        )

        # Column guide
        column_guide_html = DATA_COLUMN_GUIDE_TEMPLATE

        # Build table rows with additional columns
        rows_html = ""
        for row in table_df.iter_rows(named=True):
            genome = row.get("best_match_genome", "")

            # Get species and genus from metadata if available
            species = ""
            genus = ""
            if self.genome_metadata is not None:
                species = self.genome_metadata.get_species(genome) or ""
                genus = self.genome_metadata.get_genus(genome) or ""

            # Get additional columns from dataframe
            ambiguous_hits = row.get("num_ambiguous_hits", 0)
            identity_gap = row.get("identity_gap")
            # Format identity_gap: show "-" if not available, otherwise format as number
            identity_gap_str = "-" if identity_gap is None else f"{identity_gap:.2f}"
            is_novel = row.get("is_novel", False)
            if is_novel is None:
                # Compute from taxonomic_call if not in dataframe
                call = row.get("taxonomic_call", "")
                is_novel = call in ("Novel Species", "Novel Genus")

            rows_html += TABLE_ROW_TEMPLATE.format(
                read_id=row.get("read_id", ""),
                best_match=genome[:40] if len(genome) > 40 else genome,
                species=species[:30] if len(species) > 30 else species,
                genus=genus[:20] if len(genus) > 20 else genus,
                identity=row.get("top_hit_identity", 0),
                novelty=row.get("novelty_index", 0),
                uncertainty=row.get("placement_uncertainty", 0),
                ambiguous_hits=ambiguous_hits,
                identity_gap=identity_gap_str,
                is_novel="Yes" if is_novel else "No",
                classification=row.get("taxonomic_call", ""),
                cell_class=get_cell_class(row.get("taxonomic_call", "")),
            )

        table_html = DATA_TABLE_TEMPLATE.format(
            table_rows=rows_html,
            page_size=self.config.page_size,
            total_rows=len(table_df),
        )

        # Combine all sections
        content = summary_html + quick_filters_html + column_guide_html + table_html

        return TAB_SECTION_TEMPLATE.format(
            tab_id="data",
            active_class="",
            section_title="Classification Data",
            content=content,
        )

    def _register_plot(self, plot_id: str, fig: go.Figure) -> None:
        """Register a plot for later JS initialization."""
        plot_json = fig.to_json()
        self._plot_data[plot_id] = plot_json

    def _build_plotly_js(self) -> str:
        """Build Plotly.js initialization code for all plots."""
        js_lines = ["<script>"]
        js_lines.append("document.addEventListener('DOMContentLoaded', function() {")

        for plot_id, plot_json in self._plot_data.items():
            js_lines.append(f"  var data_{plot_id.replace('-', '_')} = {plot_json};")
            js_lines.append(
                f"  Plotly.newPlot('{plot_id}', "
                f"data_{plot_id.replace('-', '_')}.data, "
                f"data_{plot_id.replace('-', '_')}.layout, "
                f"{{responsive: true}});"
            )

        js_lines.append("});")
        js_lines.append("</script>")

        return "\n".join(js_lines)


def generate_report(
    classifications_path: Path | str,
    output_path: Path | str,
    sample_name: str = "Sample",
    ani_matrix_path: Path | str | None = None,
    recruitment_data_path: Path | str | None = None,
    theme: str = "light",
) -> None:
    """
    Convenience function to generate a report from file paths.

    Args:
        classifications_path: Path to classifications CSV file
        output_path: Output path for HTML report
        sample_name: Sample name for report header
        ani_matrix_path: Optional path to ANI matrix CSV
        recruitment_data_path: Optional path to recruitment data CSV
        theme: Color theme ('light' or 'dark')
    """
    # Load classifications
    classifications = read_dataframe(classifications_path)

    # Load optional data
    ani_matrix = None
    if ani_matrix_path:
        ani_matrix = read_dataframe(ani_matrix_path)

    recruitment_data = None
    if recruitment_data_path:
        recruitment_data = read_dataframe(recruitment_data_path)

    # Create config
    config = ReportConfig(
        sample_name=sample_name,
        theme=theme,
    )

    # Generate report
    generator = ReportGenerator(
        classifications=classifications,
        config=config,
        recruitment_data=recruitment_data,
        ani_matrix=ani_matrix,
    )
    generator.generate(output_path)
