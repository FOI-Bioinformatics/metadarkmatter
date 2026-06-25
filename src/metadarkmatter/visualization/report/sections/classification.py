"""ReportGenerator classification section."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import polars as pl

from metadarkmatter.visualization.plots.base import (
    BAYESIAN_CATEGORY_COLORS,
    PlotConfig,
)
from metadarkmatter.visualization.plots.classification_charts import (
    DiversitySunburstChart,
)
from metadarkmatter.visualization.plots.distributions import (
    NoveltyHistogram,
    UncertaintyHistogram,
)
from metadarkmatter.visualization.plots.scatter_2d import (
    NoveltyUncertaintyScatter,
)
from metadarkmatter.visualization.report.models import (
    _safe_float,
)
from metadarkmatter.visualization.report.sections._base import _ReportBase
from metadarkmatter.visualization.report.templates import (
    BAYESIAN_INTERPRETATION_TEMPLATE,
    COLLAPSIBLE_PANEL_TEMPLATE,
    CONFIDENCE_SUMMARY_TEMPLATE,
    ENHANCED_SCORING_SUMMARY_TEMPLATE,
    ENHANCED_SCORING_UNCERTAINTY_TYPES_TEMPLATE,
    KPI_CARD_TEMPLATE,
    KPI_STRIP_TEMPLATE,
    PLOT_CONTAINER_SIMPLE_TEMPLATE,
    PLOT_CONTAINER_TEMPLATE,
    PLOT_ROW_TEMPLATE,
    TAB_SECTION_TEMPLATE,
    TWO_COLUMN_ROW_TEMPLATE,
)

if TYPE_CHECKING:

    pass


logger = logging.getLogger(__name__)


class ClassificationMixin(_ReportBase):
    """Classification report-section builders."""

    def _build_classification_section(self) -> str:
        """Build the Classification tab with scatter, sunburst, histograms, and collapsible confidence."""
        s = self.summary

        # 1. Metric Strip (KPI-style)
        kpi_cards = []
        kpi_cards.append(KPI_CARD_TEMPLATE.format(
            color_class="accent",
            value=f"{s.mean_novelty_index:.1f}",
            label="Mean Novelty",
        ))
        kpi_cards.append(KPI_CARD_TEMPLATE.format(
            color_class="accent",
            value=f"{s.mean_placement_uncertainty:.1f}",
            label="Mean Uncertainty",
        ))
        if s.has_bayesian:
            kpi_cards.append(KPI_CARD_TEMPLATE.format(
                color_class="accent",
                value=f"{s.mean_posterior_entropy:.2f}",
                label="Mean Entropy",
            ))
            kpi_cards.append(KPI_CARD_TEMPLATE.format(
                color_class="known",
                value=f"{s.high_confidence_pct:.0f}%",
                label="High Confidence",
            ))
        kpi_html = KPI_STRIP_TEMPLATE.format(cards="\n".join(kpi_cards))

        # 2. Scatter plot (primary diagnostic)
        scatter = NoveltyUncertaintyScatter(
            self.df,
            config=PlotConfig(width=800, height=550),
            thresholds=self.config.thresholds,
            max_points=self.config.max_scatter_points,
        )
        scatter_id = "plot-scatter-2d"
        self._register_plot(scatter_id, scatter.create_figure())

        scatter_plot = PLOT_CONTAINER_SIMPLE_TEMPLATE.format(
            extra_class="",
            plot_id=scatter_id,
        )

        # 3. Sunburst
        sunburst = DiversitySunburstChart(
            s.to_dict(),
            config=PlotConfig(width=500, height=550),
            title="Diversity Classification",
        )
        sunburst_id = "plot-diversity-sunburst"
        self._register_plot(sunburst_id, sunburst.create_figure())

        sunburst_plot = PLOT_CONTAINER_SIMPLE_TEMPLATE.format(
            extra_class="",
            plot_id=sunburst_id,
        )

        # Two-column: scatter (60%) + sunburst (40%)
        scatter_row = TWO_COLUMN_ROW_TEMPLATE.format(
            left_flex="3",
            right_flex="2",
            left_content=scatter_plot,
            right_content=sunburst_plot,
        )

        # 4. Histograms side-by-side
        novelty_hist = NoveltyHistogram(
            self.df,
            config=PlotConfig(width=550, height=400),
            thresholds=self.config.thresholds,
            bin_size=self.config.novelty_bin_size,
        )
        novelty_id = "plot-novelty-hist"
        self._register_plot(novelty_id, novelty_hist.create_figure())

        uncertainty_hist = UncertaintyHistogram(
            self.df,
            config=PlotConfig(width=550, height=400),
            thresholds=self.config.thresholds,
            bin_size=self.config.uncertainty_bin_size,
        )
        uncertainty_id = "plot-uncertainty-hist"
        self._register_plot(uncertainty_id, uncertainty_hist.create_figure())

        hist_row = PLOT_ROW_TEMPLATE.format(
            plots=(
                PLOT_CONTAINER_SIMPLE_TEMPLATE.format(
                    extra_class="half-width",
                    plot_id=novelty_id,
                )
                + PLOT_CONTAINER_SIMPLE_TEMPLATE.format(
                    extra_class="half-width",
                    plot_id=uncertainty_id,
                )
            )
        )

        # 5. Confidence Analysis (collapsible, default collapsed)
        confidence_content = self._build_confidence_panel_content()
        if confidence_content:
            confidence_panel = COLLAPSIBLE_PANEL_TEMPLATE.format(
                default_state="",
                panel_id="confidence-analysis",
                title="Confidence Analysis",
                content=confidence_content,
            )
        else:
            confidence_panel = ""

        content = kpi_html + scatter_row + hist_row + confidence_panel

        return TAB_SECTION_TEMPLATE.format(
            tab_id="classification",
            active_class="",
            section_title="Classification",
            content=content,
        )

    def _build_confidence_panel_content(self) -> str:
        """Build inner HTML content for the collapsible Confidence Analysis panel.

        Extracts the body content from what was previously the standalone
        Classification Confidence tab. Returns empty string if no Bayesian
        or enhanced scoring data is available.
        """
        s = self.summary

        if not (s.has_bayesian or s.has_enhanced_scoring or s.has_inferred_uncertainty):
            return ""

        parts = [
            self._confidence_summary_cards(),
            self._confidence_entropy_distribution(),
            self._confidence_posterior_bar(),
            self._confidence_discovery_distribution(),
            self._confidence_uncertainty_types(),
            self._confidence_landscape_scatter(),
        ]
        return "\n".join(p for p in parts if p)

    def _confidence_summary_cards(self) -> str:
        """Build the merged Bayesian / enhanced-scoring summary cards."""
        s = self.summary
        parts: list[str] = []
        # --- Merged summary cards ---
        if s.has_bayesian:
            summary_html = CONFIDENCE_SUMMARY_TEMPLATE.format(
                mean_entropy=s.mean_posterior_entropy,
                high_confidence_pct=s.high_confidence_pct,
                high_confidence_count=s.high_confidence_count,
                map_agreement_pct=s.map_agreement_pct,
                map_agreement_count=s.map_agreement_count,
                total_reads=s.total_reads,
                single_hit_pct=s.single_hit_pct,
                single_hit_count=s.single_hit_count,
                mean_discovery_score=s.mean_discovery_score or 0.0,
                novel_count=s.novel_with_discovery_score,
                boundary_pct=s.boundary_pct,
                boundary_count=s.boundary_count,
            )
            parts.append(summary_html)
            parts.append(BAYESIAN_INTERPRETATION_TEMPLATE)
        elif s.has_enhanced_scoring or s.has_inferred_uncertainty:
            # Fallback: show discovery-only summary when no Bayesian data
            summary_html = ENHANCED_SCORING_SUMMARY_TEMPLATE.format(
                single_hit_pct=s.single_hit_pct,
                single_hit_count=s.single_hit_count,
                total_reads=s.total_reads,
                mean_inferred_uncertainty=s.mean_inferred_uncertainty or 0.0,
                high_priority_count=s.high_priority_discoveries,
                mean_discovery_score=s.mean_discovery_score or 0.0,
                novel_count=s.novel_with_discovery_score,
            )
            parts.append(summary_html)

        return "\n".join(parts)

    def _confidence_entropy_distribution(self) -> str:
        """Build the posterior-entropy distribution histogram block."""
        s = self.summary
        parts: list[str] = []
        # --- Posterior Entropy Distribution ---
        if s.has_bayesian and "posterior_entropy" in self.df.columns:
            entropy_id = "plot-bayesian-entropy"
            self._build_histogram(
                data=self.df["posterior_entropy"].to_list(),
                plot_id=entropy_id,
                title="Posterior Entropy Distribution",
                x_label="Shannon Entropy (bits)",
                nbins=40,
                thresholds=[
                    (0.5, "dash", "#22c55e", "Very confident"),
                    (1.0, "dot", "#f59e0b", "Confident"),
                    (1.5, "dot", "#ef4444", "Uncertain"),
                ],
                x_range=[0, 2.7],
            )

            parts.append(PLOT_CONTAINER_TEMPLATE.format(
                extra_class="full-width",
                title="Posterior Entropy Distribution",
                description=(
                    "Shannon entropy of the posterior distribution per read. "
                    "Low entropy (left) indicates confident classification; "
                    "high entropy (right) indicates reads near boundaries "
                    "where multiple categories are plausible."
                ),
                plot_id=entropy_id,
            ))

        return "\n".join(parts)

    def _confidence_posterior_bar(self) -> str:
        """Build the stacked mean-posterior bar chart by MAP category."""
        import plotly.graph_objects as go

        s = self.summary
        parts: list[str] = []
        # --- Stacked Posterior Bar by Category ---
        if s.has_bayesian and "taxonomic_call" in self.df.columns:
            posterior_cols = [
                ("p_known_species", "Known Species"),
                ("p_novel_species", "Novel Species"),
                ("p_novel_genus", "Novel Genus"),
                ("p_species_boundary", "Species Boundary"),
                ("p_ambiguous", "Ambiguous"),
                ("p_unclassified", "Unclassified"),
            ]

            bar_fig = go.Figure()
            all_map_cats = [
                "Known Species", "Novel Species", "Novel Genus",
                "Species Boundary", "Ambiguous", "Unclassified",
            ]
            for col_name, cat_label in posterior_cols:
                if col_name not in self.df.columns:
                    continue
                means = []
                cats = []
                for map_cat in all_map_cats:
                    cat_df = self.df.filter(pl.col("taxonomic_call") == map_cat)
                    if len(cat_df) > 0:
                        cats.append(map_cat)
                        means.append(_safe_float(cat_df[col_name].mean()) or 0.0)
                if means:
                    bar_fig.add_trace(go.Bar(
                        x=cats,
                        y=means,
                        name=cat_label,
                        marker_color=BAYESIAN_CATEGORY_COLORS.get(cat_label, "#94a3b8"),
                    ))

            bar_fig.update_layout(
                title="Mean Posterior Composition by MAP Category",
                xaxis_title="Bayesian MAP Category",
                yaxis_title="Mean Posterior Probability",
                template="plotly_white",
                height=450,
                barmode="stack",
                legend={"orientation": "h", "y": -0.2},
                yaxis={"range": [0, 1.05]},
            )
            bar_id = "plot-bayesian-posterior-bar"
            self._register_plot(bar_id, bar_fig)

            parts.append(PLOT_CONTAINER_TEMPLATE.format(
                extra_class="full-width",
                title="Posterior Composition by Category",
                description=(
                    "Mean posterior probabilities grouped by the Bayesian MAP classification. "
                    "A well-separated classifier shows tall bars for the dominant category "
                    "in each group. Mixed bars indicate regions of classification uncertainty."
                ),
                plot_id=bar_id,
            ))

        return "\n".join(parts)

    def _confidence_discovery_distribution(self) -> str:
        """Build the discovery-score distribution histogram block."""
        s = self.summary
        parts: list[str] = []
        # --- Discovery Score Distribution ---
        if s.has_enhanced_scoring and "discovery_score" in self.df.columns:
            discovery_df = self.df.filter(pl.col("discovery_score").is_not_null())
            if len(discovery_df) > 0:
                discovery_hist_id = "plot-discovery-hist"
                self._build_histogram(
                    data=discovery_df["discovery_score"].to_list(),
                    plot_id=discovery_hist_id,
                    title="Discovery Score Distribution (Novel Reads)",
                    x_label="Discovery Score",
                    nbins=25,
                    thresholds=[
                        (75, "dash", "#22c55e", "High Priority (75+)"),
                        (50, "dot", "#f59e0b", "Moderate (50+)"),
                        (25, "dot", "#ef4444", "Low (25+)"),
                    ],
                )

                parts.append(PLOT_CONTAINER_TEMPLATE.format(
                    extra_class="full-width",
                    title="Discovery Score Distribution",
                    description=(
                        "Distribution of discovery scores for novel reads. Higher scores indicate "
                        "more reliable discoveries. Green line marks high-priority threshold (75+)."
                    ),
                    plot_id=discovery_hist_id,
                ))

        return "\n".join(parts)

    def _confidence_uncertainty_types(self) -> str:
        """Build the measured-vs-inferred uncertainty breakdown block."""
        s = self.summary
        parts: list[str] = []
        # --- Uncertainty types breakdown ---
        if s.has_inferred_uncertainty and "uncertainty_type" in self.df.columns:
            measured_df = self.df.filter(pl.col("uncertainty_type") == "measured")
            inferred_df = self.df.filter(pl.col("uncertainty_type") == "inferred")
            measured_count = len(measured_df)
            inferred_count = len(inferred_df)
            total = s.total_reads if s.total_reads > 0 else 1

            uncertainty_types_html = ENHANCED_SCORING_UNCERTAINTY_TYPES_TEMPLATE.format(
                measured_count=measured_count,
                measured_pct=measured_count / total * 100,
                inferred_count=inferred_count,
                inferred_pct=inferred_count / total * 100,
                single_hit_pct=s.single_hit_pct,
            )
            parts.append(uncertainty_types_html)

        return "\n".join(parts)

    def _confidence_landscape_scatter(self) -> str:
        """Build the novelty-vs-entropy confidence-landscape scatter block."""
        s = self.summary
        parts: list[str] = []
        # --- Confidence Landscape Scatter ---
        if s.has_bayesian and "posterior_entropy" in self.df.columns and "novelty_index" in self.df.columns:
            scatter_id = "plot-bayesian-landscape"
            self._build_category_scatter(
                x_col="novelty_index",
                y_col="posterior_entropy",
                plot_id=scatter_id,
                title="Classification Confidence Landscape",
                x_label="Novelty Index (%)",
                y_label="Posterior Entropy (bits)",
                y_range=[0, 2.7],
            )

            parts.append(PLOT_CONTAINER_TEMPLATE.format(
                extra_class="full-width",
                title="Confidence Landscape",
                description=(
                    "Novelty index vs posterior entropy colored by Bayesian MAP category. "
                    "Reads with high entropy (top) lie near classification boundaries. "
                    "Vertical bands of high entropy correspond to threshold boundaries."
                ),
                plot_id=scatter_id,
            ))

        return "\n".join(parts)
