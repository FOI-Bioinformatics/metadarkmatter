"""
2D scatter plot: Novelty Index vs Placement Uncertainty.

The central diagnostic visualization for classification interpretation.
Shows the distribution of reads in novelty-uncertainty space with
classification regions overlaid.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import plotly.graph_objects as go

from metadarkmatter.visualization.plots.base import (
    TAXONOMY_COLORS,
    BasePlot,
    PlotConfig,
    ThresholdConfig,
    subsample_dataframe,
)

if TYPE_CHECKING:
    import polars as pl


class NoveltyUncertaintyScatter(BasePlot):
    """
    2D scatter plot showing novelty vs uncertainty with points
    colored by taxonomic classification.

    This is the key diagnostic plot that visualizes:
    - X-axis: Novelty Index (100 - top hit identity)
    - Y-axis: Placement Uncertainty (100 - ANI between top hits)
    - Color: Taxonomic classification category

    Includes classification region boxes showing thresholds for
    Known Species, Novel Species, Novel Genus, and Conserved Region.

    Supports toggling visibility of single-hit reads (reads with only one
    BLAST hit) via dropdown menu.
    """

    def __init__(
        self,
        data: pl.DataFrame,
        config: PlotConfig | None = None,
        thresholds: ThresholdConfig | None = None,
        max_points: int = 50000,
        show_regions: bool = True,
        show_legend: bool = True,
        title: str = "Novelty Index vs Placement Uncertainty",
        enable_single_hit_toggle: bool = True,
    ) -> None:
        """
        Initialize 2D scatter plot.

        Args:
            data: DataFrame with columns: novelty_index, placement_uncertainty, taxonomic_call
            config: Plot configuration
            thresholds: Classification thresholds for region boxes
            max_points: Maximum points to display (subsampling for performance)
            show_regions: Whether to show classification region boxes
            show_legend: Whether to show legend
            title: Plot title
            enable_single_hit_toggle: Whether to add toggle for single-hit reads
        """
        super().__init__(config, thresholds)
        self.data = data
        self.max_points = max_points
        self.show_regions = show_regions
        self.show_legend = show_legend
        self.title = title
        self.enable_single_hit_toggle = enable_single_hit_toggle

    def create_figure(self) -> go.Figure:
        """Create 2D scatter with classification regions."""
        import polars as pl

        # Subsample if needed for performance
        plot_data = subsample_dataframe(self.data, self.max_points)

        fig = go.Figure()

        # Add classification regions first (behind points)
        if self.show_regions:
            self._add_classification_regions(fig)

        # Check if we can enable single-hit toggle
        has_ambiguous_hits_col = "num_ambiguous_hits" in plot_data.columns
        use_toggle = self.enable_single_hit_toggle and has_ambiguous_hits_col

        # Add scatter traces for each classification category
        categories = [
            ("Known Species", TAXONOMY_COLORS["Known Species"]),
            ("Novel Species", TAXONOMY_COLORS["Novel Species"]),
            ("Novel Genus", TAXONOMY_COLORS["Novel Genus"]),
            ("Conserved Region", TAXONOMY_COLORS["Conserved Region"]),
        ]

        # Track trace names for toggle functionality
        trace_names = []

        for category, color in categories:
            category_data = plot_data.filter(pl.col("taxonomic_call") == category)

            if len(category_data) == 0:
                continue

            if use_toggle:
                # Split into single-hit and multi-hit traces
                single_hit = category_data.filter(pl.col("num_ambiguous_hits") <= 1)
                multi_hit = category_data.filter(pl.col("num_ambiguous_hits") > 1)

                # Add multi-hit trace (always visible)
                if len(multi_hit) > 0:
                    fig.add_trace(
                        go.Scattergl(
                            x=multi_hit["novelty_index"].to_list(),
                            y=multi_hit["placement_uncertainty"].to_list(),
                            mode="markers",
                            name=f"{category} multi-hit ({len(multi_hit):,})",
                            legendgroup=category,
                            marker={
                                "color": color,
                                "size": 5,
                                "opacity": 0.6,
                            },
                            hovertemplate=(
                                f"<b>{category}</b> (multi-hit)<br>"
                                "Novelty: %{x:.2f}<br>"
                                "Uncertainty: %{y:.2f}<br>"
                                "<extra></extra>"
                            ),
                        )
                    )
                    trace_names.append(("multi", category))

                # Add single-hit trace (can be toggled)
                if len(single_hit) > 0:
                    fig.add_trace(
                        go.Scattergl(
                            x=single_hit["novelty_index"].to_list(),
                            y=single_hit["placement_uncertainty"].to_list(),
                            mode="markers",
                            name=f"{category} single-hit ({len(single_hit):,})",
                            legendgroup=category,
                            marker={
                                "color": color,
                                "size": 5,
                                "opacity": 0.4,
                                "symbol": "circle-open",
                            },
                            hovertemplate=(
                                f"<b>{category}</b> (single-hit)<br>"
                                "Novelty: %{x:.2f}<br>"
                                "Uncertainty: %{y:.2f}<br>"
                                "<extra></extra>"
                            ),
                        )
                    )
                    trace_names.append(("single", category))
            else:
                # Original behavior: single trace per category
                fig.add_trace(
                    go.Scattergl(
                        x=category_data["novelty_index"].to_list(),
                        y=category_data["placement_uncertainty"].to_list(),
                        mode="markers",
                        name=f"{category} ({len(category_data):,})",
                        marker={
                            "color": color,
                            "size": 5,
                            "opacity": 0.6,
                        },
                        hovertemplate=(
                            f"<b>{category}</b><br>"
                            "Novelty: %{x:.2f}<br>"
                            "Uncertainty: %{y:.2f}<br>"
                            "<extra></extra>"
                        ),
                    )
                )
                trace_names.append(("all", category))

        # Add toggle buttons if enabled
        if use_toggle and len(trace_names) > 0:
            # Build visibility lists for each mode
            all_visible = [True] * len(trace_names)
            multi_only = [t[0] == "multi" for t in trace_names]
            single_only = [t[0] == "single" for t in trace_names]

            fig.update_layout(
                updatemenus=[
                    {
                        "type": "dropdown",
                        "direction": "down",
                        "x": 0.01,
                        "y": 0.99,
                        "xanchor": "left",
                        "yanchor": "top",
                        "showactive": True,
                        "buttons": [
                            {
                                "label": "All Reads",
                                "method": "update",
                                "args": [{"visible": all_visible}],
                            },
                            {
                                "label": "Multi-hit Only",
                                "method": "update",
                                "args": [{"visible": multi_only}],
                            },
                            {
                                "label": "Single-hit Only",
                                "method": "update",
                                "args": [{"visible": single_only}],
                            },
                        ],
                    }
                ]
            )

        # Layout
        fig.update_layout(
            xaxis_title="Novelty Index (100 - % Identity)",
            yaxis_title="Placement Uncertainty (100 - ANI)",
            legend={
                "yanchor": "top",
                "y": 0.99,
                "xanchor": "right",
                "x": 0.99,
                "bgcolor": "rgba(255, 255, 255, 0.8)",
                "bordercolor": "#ddd",
                "borderwidth": 1,
            },
            showlegend=self.show_legend,
            **self.config.to_layout_dict(self.title),
        )

        # Set axis ranges
        max_novelty = max(30, plot_data["novelty_index"].max() * 1.1)
        max_uncertainty = max(20, plot_data["placement_uncertainty"].max() * 1.1)

        fig.update_xaxes(range=[0, max_novelty])
        fig.update_yaxes(range=[0, max_uncertainty])

        return self._apply_config(fig)

    def _add_classification_regions(self, fig: go.Figure) -> None:
        """Add shaded regions showing classification boundaries."""
        t = self.thresholds

        # Region definitions: (x0, x1, y0, y1, color, label)
        regions = [
            # Known Species region (low novelty, low uncertainty)
            (0, t.novelty_known_max, 0, t.uncertainty_known_max,
             TAXONOMY_COLORS["Known Species"], "Known Species Region"),

            # Novel Species region (medium novelty, low uncertainty)
            (t.novelty_novel_species_min, t.novelty_novel_species_max,
             0, t.uncertainty_novel_species_max,
             TAXONOMY_COLORS["Novel Species"], "Novel Species Region"),

            # Novel Genus region (high novelty, low uncertainty)
            (t.novelty_novel_genus_min, t.novelty_novel_genus_max,
             0, t.uncertainty_novel_genus_max,
             TAXONOMY_COLORS["Novel Genus"], "Novel Genus Region"),

            # Conserved Region (high uncertainty, any novelty)
            (0, 30, t.uncertainty_conserved_min, 25,
             TAXONOMY_COLORS["Conserved Region"], "Conserved Region"),
        ]

        for x0, x1, y0, y1, color, label in regions:
            fig.add_shape(
                type="rect",
                x0=x0,
                x1=x1,
                y0=y0,
                y1=y1,
                fillcolor=color,
                opacity=0.1,
                line={"width": 1, "color": color, "dash": "dot"},
                layer="below",
            )

            # Add region label
            fig.add_annotation(
                x=(x0 + x1) / 2,
                y=y1 - 0.5,
                text=label,
                showarrow=False,
                font={"size": 9, "color": color},
                opacity=0.7,
            )


class ConfidenceNoveltyScatter(BasePlot):
    """
    2D scatter plot showing confidence score vs novelty index.

    This plot helps visualize the relationship between how confident
    the classification is (confidence_score) and how divergent the
    read is from known references (novelty_index).

    Confidence score is calculated as: 100 - novelty_index - placement_uncertainty
    Higher values indicate more confident classifications.
    """

    def __init__(
        self,
        data: pl.DataFrame,
        config: PlotConfig | None = None,
        thresholds: ThresholdConfig | None = None,
        max_points: int = 50000,
        show_legend: bool = True,
        title: str = "Confidence Score vs Novelty Index",
        enable_single_hit_toggle: bool = True,
    ) -> None:
        """
        Initialize confidence vs novelty scatter plot.

        Args:
            data: DataFrame with columns: novelty_index, confidence_score, taxonomic_call
            config: Plot configuration
            thresholds: Classification thresholds
            max_points: Maximum points to display (subsampling for performance)
            show_legend: Whether to show legend
            title: Plot title
            enable_single_hit_toggle: Whether to add toggle for single-hit reads
        """
        super().__init__(config, thresholds)
        self.data = data
        self.max_points = max_points
        self.show_legend = show_legend
        self.title = title
        self.enable_single_hit_toggle = enable_single_hit_toggle

    def create_figure(self) -> go.Figure:
        """Create scatter plot of confidence score vs novelty index."""
        import polars as pl

        # Check if confidence_score column exists
        if "confidence_score" not in self.data.columns:
            # Return empty figure with message
            fig = go.Figure()
            fig.add_annotation(
                x=0.5,
                y=0.5,
                text="confidence_score column not available",
                showarrow=False,
                xref="paper",
                yref="paper",
                font={"size": 14},
            )
            return fig

        # Subsample if needed for performance
        plot_data = subsample_dataframe(self.data, self.max_points)

        fig = go.Figure()

        # Check if we can enable single-hit toggle
        has_ambiguous_hits_col = "num_ambiguous_hits" in plot_data.columns
        use_toggle = self.enable_single_hit_toggle and has_ambiguous_hits_col

        # Add scatter traces for each classification category
        categories = [
            ("Known Species", TAXONOMY_COLORS["Known Species"]),
            ("Novel Species", TAXONOMY_COLORS["Novel Species"]),
            ("Novel Genus", TAXONOMY_COLORS["Novel Genus"]),
            ("Conserved Region", TAXONOMY_COLORS["Conserved Region"]),
        ]

        # Track trace names for toggle functionality
        trace_names = []

        for category, color in categories:
            category_data = plot_data.filter(pl.col("taxonomic_call") == category)

            if len(category_data) == 0:
                continue

            if use_toggle:
                # Split into single-hit and multi-hit traces
                single_hit = category_data.filter(pl.col("num_ambiguous_hits") <= 1)
                multi_hit = category_data.filter(pl.col("num_ambiguous_hits") > 1)

                # Add multi-hit trace (always visible)
                if len(multi_hit) > 0:
                    fig.add_trace(
                        go.Scattergl(
                            x=multi_hit["novelty_index"].to_list(),
                            y=multi_hit["confidence_score"].to_list(),
                            mode="markers",
                            name=f"{category} multi-hit ({len(multi_hit):,})",
                            legendgroup=category,
                            marker={
                                "color": color,
                                "size": 5,
                                "opacity": 0.6,
                            },
                            hovertemplate=(
                                f"<b>{category}</b> (multi-hit)<br>"
                                "Novelty: %{x:.2f}<br>"
                                "Confidence: %{y:.2f}<br>"
                                "<extra></extra>"
                            ),
                        )
                    )
                    trace_names.append(("multi", category))

                # Add single-hit trace (can be toggled)
                if len(single_hit) > 0:
                    fig.add_trace(
                        go.Scattergl(
                            x=single_hit["novelty_index"].to_list(),
                            y=single_hit["confidence_score"].to_list(),
                            mode="markers",
                            name=f"{category} single-hit ({len(single_hit):,})",
                            legendgroup=category,
                            marker={
                                "color": color,
                                "size": 5,
                                "opacity": 0.4,
                                "symbol": "circle-open",
                            },
                            hovertemplate=(
                                f"<b>{category}</b> (single-hit)<br>"
                                "Novelty: %{x:.2f}<br>"
                                "Confidence: %{y:.2f}<br>"
                                "<extra></extra>"
                            ),
                        )
                    )
                    trace_names.append(("single", category))
            else:
                # Original behavior: single trace per category
                fig.add_trace(
                    go.Scattergl(
                        x=category_data["novelty_index"].to_list(),
                        y=category_data["confidence_score"].to_list(),
                        mode="markers",
                        name=f"{category} ({len(category_data):,})",
                        marker={
                            "color": color,
                            "size": 5,
                            "opacity": 0.6,
                        },
                        hovertemplate=(
                            f"<b>{category}</b><br>"
                            "Novelty: %{x:.2f}<br>"
                            "Confidence: %{y:.2f}<br>"
                            "<extra></extra>"
                        ),
                    )
                )
                trace_names.append(("all", category))

        # Add toggle buttons if enabled
        if use_toggle and len(trace_names) > 0:
            # Build visibility lists for each mode
            all_visible = [True] * len(trace_names)
            multi_only = [t[0] == "multi" for t in trace_names]
            single_only = [t[0] == "single" for t in trace_names]

            fig.update_layout(
                updatemenus=[
                    {
                        "type": "dropdown",
                        "direction": "down",
                        "x": 0.01,
                        "y": 0.99,
                        "xanchor": "left",
                        "yanchor": "top",
                        "showactive": True,
                        "buttons": [
                            {
                                "label": "All Reads",
                                "method": "update",
                                "args": [{"visible": all_visible}],
                            },
                            {
                                "label": "Multi-hit Only",
                                "method": "update",
                                "args": [{"visible": multi_only}],
                            },
                            {
                                "label": "Single-hit Only",
                                "method": "update",
                                "args": [{"visible": single_only}],
                            },
                        ],
                    }
                ]
            )

        # Add reference lines
        # Confidence threshold line (if available from thresholds)
        # confidence_score = 100 - novelty - uncertainty
        # For known species: novelty < 5%, uncertainty < 2% => confidence > 93%
        fig.add_hline(
            y=50.0,
            line_dash="dash",
            line_color="gray",
            opacity=0.5,
            annotation_text="50% confidence threshold",
            annotation_position="bottom right",
        )

        # Layout
        fig.update_layout(
            xaxis_title="Novelty Index (100 - % Identity)",
            yaxis_title="Confidence Score (%)",
            legend={
                "yanchor": "top",
                "y": 0.99,
                "xanchor": "right",
                "x": 0.99,
                "bgcolor": "rgba(255, 255, 255, 0.8)",
                "bordercolor": "#ddd",
                "borderwidth": 1,
            },
            showlegend=self.show_legend,
            **self.config.to_layout_dict(self.title),
        )

        # Set axis ranges
        max_novelty = max(30, plot_data["novelty_index"].max() * 1.1)

        fig.update_xaxes(range=[0, max_novelty])
        fig.update_yaxes(range=[0, 105])

        return self._apply_config(fig)


class NoveltyUncertaintyDensity(BasePlot):
    """
    2D density/heatmap plot showing concentration of reads
    in novelty-uncertainty space.

    Alternative to scatter plot for very large datasets where
    individual points would be overwhelming.
    """

    def __init__(
        self,
        data: pl.DataFrame,
        config: PlotConfig | None = None,
        thresholds: ThresholdConfig | None = None,
        nbinsx: int = 50,
        nbinsy: int = 50,
        title: str = "Read Density: Novelty vs Uncertainty",
    ) -> None:
        """
        Initialize density plot.

        Args:
            data: DataFrame with novelty_index and placement_uncertainty columns
            config: Plot configuration
            thresholds: Classification thresholds
            nbinsx: Number of bins on x-axis
            nbinsy: Number of bins on y-axis
            title: Plot title
        """
        super().__init__(config, thresholds)
        self.data = data
        self.nbinsx = nbinsx
        self.nbinsy = nbinsy
        self.title = title

    def create_figure(self) -> go.Figure:
        """Create 2D density heatmap."""
        novelty = self.data["novelty_index"].to_list()
        uncertainty = self.data["placement_uncertainty"].to_list()

        fig = go.Figure()

        fig.add_trace(
            go.Histogram2d(
                x=novelty,
                y=uncertainty,
                nbinsx=self.nbinsx,
                nbinsy=self.nbinsy,
                colorscale="Viridis",
                colorbar={"title": "Read Count"},
                hovertemplate=(
                    "Novelty: %{x:.1f}<br>"
                    "Uncertainty: %{y:.1f}<br>"
                    "Count: %{z}<extra></extra>"
                ),
            )
        )

        # Add threshold lines
        t = self.thresholds

        # Vertical lines for novelty thresholds
        for value in [t.novelty_known_max, t.novelty_novel_species_min,
                      t.novelty_novel_species_max, t.novelty_novel_genus_max]:
            fig.add_vline(x=value, line_dash="dash", line_color="white", opacity=0.5)

        # Horizontal lines for uncertainty thresholds
        for value in [t.uncertainty_known_max, t.uncertainty_novel_genus_max,
                      t.uncertainty_conserved_min]:
            fig.add_hline(y=value, line_dash="dash", line_color="white", opacity=0.5)

        fig.update_layout(
            xaxis_title="Novelty Index",
            yaxis_title="Placement Uncertainty",
            **self.config.to_layout_dict(self.title),
        )

        return self._apply_config(fig)


class ClassificationScatterMatrix(BasePlot):
    """
    Scatter matrix showing relationships between all numeric
    classification metrics: novelty, uncertainty, identity.
    """

    def __init__(
        self,
        data: pl.DataFrame,
        config: PlotConfig | None = None,
        thresholds: ThresholdConfig | None = None,
        max_points: int = 10000,
    ) -> None:
        """
        Initialize scatter matrix.

        Args:
            data: DataFrame with numeric classification columns
            config: Plot configuration
            thresholds: Classification thresholds
            max_points: Maximum points per subplot
        """
        super().__init__(config, thresholds)
        self.data = data
        self.max_points = max_points

    def create_figure(self) -> go.Figure:
        """Create scatter matrix of classification metrics."""
        import plotly.express as px

        # Subsample
        plot_data = subsample_dataframe(self.data, self.max_points)

        # Convert to pandas for plotly express
        df = plot_data.select([
            "novelty_index",
            "placement_uncertainty",
            "top_hit_identity",
            "taxonomic_call",
        ]).to_pandas()

        fig = px.scatter_matrix(
            df,
            dimensions=["novelty_index", "placement_uncertainty", "top_hit_identity"],
            color="taxonomic_call",
            color_discrete_map=TAXONOMY_COLORS,
            opacity=0.5,
            title="Classification Metrics Scatter Matrix",
        )

        fig.update_traces(diagonal_visible=False, showupperhalf=False)
        fig.update_layout(**self.config.to_layout_dict())

        return self._apply_config(fig)
