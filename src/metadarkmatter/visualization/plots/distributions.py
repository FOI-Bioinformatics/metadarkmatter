"""
Distribution plots: histograms for novelty index, placement uncertainty,
and percent identity distributions.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import plotly.graph_objects as go

from metadarkmatter.visualization.plots.base import (
    TAXONOMY_COLORS,
    BasePlot,
    PlotConfig,
    ThresholdConfig,
)

if TYPE_CHECKING:
    import polars as pl


class NoveltyHistogram(BasePlot):
    """
    Histogram showing distribution of novelty index values.

    Includes vertical threshold lines indicating classification boundaries:
    - 2.0: Known species boundary
    - 5.0: Novel species minimum
    - 15.0: Novel species/genus boundary
    - 25.0: Novel genus maximum
    """

    def __init__(
        self,
        data: pl.DataFrame,
        config: PlotConfig | None = None,
        thresholds: ThresholdConfig | None = None,
        nbins: int | None = None,
        bin_size: float | None = None,
        show_thresholds: bool = True,
        title: str = "Novelty Index Distribution",
    ) -> None:
        """
        Initialize novelty histogram.

        Args:
            data: DataFrame with 'novelty_index' column
            config: Plot configuration
            thresholds: Classification thresholds
            nbins: Number of histogram bins (ignored if bin_size is set)
            bin_size: Explicit bin width (e.g., 1.0 for 1% bins, 0.5 for 0.5% bins)
            show_thresholds: Whether to show threshold lines
            title: Plot title
        """
        super().__init__(config, thresholds)
        self.data = data
        self.nbins = nbins if nbins is not None else 50
        self.bin_size = bin_size
        self.show_thresholds = show_thresholds
        self.title = title

    def create_figure(self) -> go.Figure:
        """Create novelty index histogram with classification thresholds."""
        novelty_values = self.data["novelty_index"].to_list()

        fig = go.Figure()

        # Build histogram kwargs
        hist_kwargs = {
            "x": novelty_values,
            "name": "Novelty Index",
            "marker_color": "#667eea",
            "opacity": 0.75,
            "hovertemplate": "Novelty: %{x:.1f}<br>Count: %{y}<extra></extra>",
        }

        # Use explicit bin_size if provided, otherwise use nbins
        if self.bin_size is not None:
            hist_kwargs["xbins"] = {"size": self.bin_size}
        else:
            hist_kwargs["nbinsx"] = self.nbins

        # Main histogram
        fig.add_trace(go.Histogram(**hist_kwargs))

        # Add threshold lines
        if self.show_thresholds:
            thresholds = [
                (
                    self.thresholds.novelty_known_max,
                    "Known",
                    TAXONOMY_COLORS["Known Species"],
                ),
                (
                    self.thresholds.novelty_novel_species_min,
                    "Novel Sp. Min",
                    TAXONOMY_COLORS["Novel Species"],
                ),
                (
                    self.thresholds.novelty_novel_species_max,
                    "Novel Sp. Max",
                    TAXONOMY_COLORS["Novel Species"],
                ),
                (
                    self.thresholds.novelty_novel_genus_max,
                    "Novel Gen. Max",
                    TAXONOMY_COLORS["Novel Genus"],
                ),
            ]

            for value, label, color in thresholds:
                fig.add_vline(
                    x=value,
                    line_dash="dash",
                    line_color=color,
                    opacity=0.8,
                    annotation_text=label,
                    annotation_position="top",
                    annotation_font_size=9,
                )

        # Layout
        fig.update_layout(
            xaxis_title="Novelty Index (100 - % Identity)",
            yaxis_title="Number of Reads",
            bargap=0.05,
            **self.config.to_layout_dict(self.title),
        )

        fig.update_xaxes(range=[0, max(30, max(novelty_values) * 1.1)])

        return self._apply_config(fig)


class UncertaintyHistogram(BasePlot):
    """
    Histogram showing distribution of placement uncertainty values.

    Includes vertical threshold lines indicating classification boundaries:
    - 0.5: Known/novel species max uncertainty
    - 2.0: Novel genus max uncertainty
    - 5.0: Conserved region min uncertainty
    """

    def __init__(
        self,
        data: pl.DataFrame,
        config: PlotConfig | None = None,
        thresholds: ThresholdConfig | None = None,
        nbins: int | None = None,
        bin_size: float | None = None,
        show_thresholds: bool = True,
        title: str = "Placement Uncertainty Distribution",
    ) -> None:
        """
        Initialize uncertainty histogram.

        Args:
            data: DataFrame with 'placement_uncertainty' column
            config: Plot configuration
            thresholds: Classification thresholds
            nbins: Number of histogram bins (ignored if bin_size is set)
            bin_size: Explicit bin width (e.g., 1.0 for 1% bins, 0.5 for 0.5% bins)
            show_thresholds: Whether to show threshold lines
            title: Plot title
        """
        super().__init__(config, thresholds)
        self.data = data
        self.nbins = nbins if nbins is not None else 50
        self.bin_size = bin_size
        self.show_thresholds = show_thresholds
        self.title = title

    def create_figure(self) -> go.Figure:
        """Create placement uncertainty histogram with thresholds."""
        uncertainty_values = self.data["placement_uncertainty"].to_list()

        fig = go.Figure()

        # Build histogram kwargs
        hist_kwargs = {
            "x": uncertainty_values,
            "name": "Placement Uncertainty",
            "marker_color": "#764ba2",
            "opacity": 0.75,
            "hovertemplate": "Uncertainty: %{x:.1f}<br>Count: %{y}<extra></extra>",
        }

        # Use explicit bin_size if provided, otherwise use nbins
        if self.bin_size is not None:
            hist_kwargs["xbins"] = {"size": self.bin_size}
        else:
            hist_kwargs["nbinsx"] = self.nbins

        # Main histogram
        fig.add_trace(go.Histogram(**hist_kwargs))

        # Add threshold lines
        if self.show_thresholds:
            thresholds = [
                (
                    self.thresholds.uncertainty_known_max,
                    "Known Max",
                    TAXONOMY_COLORS["Known Species"],
                ),
                (
                    self.thresholds.uncertainty_novel_genus_max,
                    "Novel Max",
                    TAXONOMY_COLORS["Novel Genus"],
                ),
                (
                    self.thresholds.uncertainty_conserved_min,
                    "Conserved Min",
                    TAXONOMY_COLORS["Conserved Region"],
                ),
            ]

            for value, label, color in thresholds:
                fig.add_vline(
                    x=value,
                    line_dash="dash",
                    line_color=color,
                    opacity=0.8,
                    annotation_text=label,
                    annotation_position="top",
                    annotation_font_size=9,
                )

        # Layout
        fig.update_layout(
            xaxis_title="Placement Uncertainty (100 - ANI)",
            yaxis_title="Number of Reads",
            bargap=0.05,
            **self.config.to_layout_dict(self.title),
        )

        fig.update_xaxes(range=[0, max(20, max(uncertainty_values) * 1.1)])

        return self._apply_config(fig)


class IdentityHistogram(BasePlot):
    """
    Histogram showing distribution of top hit identity values.

    Shows the percent identity of the best BLAST hit for each read,
    with threshold bands indicating classification regions.
    """

    def __init__(
        self,
        data: pl.DataFrame,
        config: PlotConfig | None = None,
        thresholds: ThresholdConfig | None = None,
        nbins: int = 50,
        show_thresholds: bool = True,
        title: str = "Top Hit Identity Distribution",
    ) -> None:
        """
        Initialize identity histogram.

        Args:
            data: DataFrame with 'top_hit_identity' column
            config: Plot configuration
            thresholds: Classification thresholds
            nbins: Number of histogram bins
            show_thresholds: Whether to show threshold lines
            title: Plot title
        """
        super().__init__(config, thresholds)
        self.data = data
        self.nbins = nbins
        self.show_thresholds = show_thresholds
        self.title = title

    def create_figure(self) -> go.Figure:
        """Create identity histogram with classification regions."""
        identity_values = self.data["top_hit_identity"].to_list()

        fig = go.Figure()

        # Main histogram
        fig.add_trace(
            go.Histogram(
                x=identity_values,
                nbinsx=self.nbins,
                name="Top Hit Identity",
                marker_color="#2ecc71",
                opacity=0.75,
                hovertemplate="Identity: %{x:.1f}%<br>Count: %{y}<extra></extra>",
            )
        )

        # Add threshold lines
        if self.show_thresholds:
            thresholds = [
                (
                    self.thresholds.identity_known_min,
                    "Known Species",
                    TAXONOMY_COLORS["Known Species"],
                ),
                (
                    self.thresholds.identity_novel_species_min,
                    "Novel Species",
                    TAXONOMY_COLORS["Novel Species"],
                ),
                (
                    self.thresholds.identity_novel_genus_min,
                    "Novel Genus",
                    TAXONOMY_COLORS["Novel Genus"],
                ),
            ]

            for value, label, color in thresholds:
                fig.add_vline(
                    x=value,
                    line_dash="dash",
                    line_color=color,
                    opacity=0.8,
                    annotation_text=label,
                    annotation_position="top",
                    annotation_font_size=9,
                )

        # Layout
        fig.update_layout(
            xaxis_title="Top Hit Identity (%)",
            yaxis_title="Number of Reads",
            bargap=0.05,
            **self.config.to_layout_dict(self.title),
        )

        fig.update_xaxes(range=[min(70, min(identity_values) - 5), 100])

        return self._apply_config(fig)


class CombinedDistributionPlot(BasePlot):
    """
    Combined subplot showing novelty and uncertainty distributions side by side.

    Provides a compact two-panel view of both key metrics.
    """

    def __init__(
        self,
        data: pl.DataFrame,
        config: PlotConfig | None = None,
        thresholds: ThresholdConfig | None = None,
        nbins: int = 40,
        show_thresholds: bool = True,
    ) -> None:
        """
        Initialize combined distribution plot.

        Args:
            data: DataFrame with 'novelty_index' and 'placement_uncertainty' columns
            config: Plot configuration
            thresholds: Classification thresholds
            nbins: Number of histogram bins
            show_thresholds: Whether to show threshold lines
        """
        super().__init__(config, thresholds)
        self.data = data
        self.nbins = nbins
        self.show_thresholds = show_thresholds

    def create_figure(self) -> go.Figure:
        """Create combined novelty and uncertainty histograms."""
        from plotly.subplots import make_subplots

        novelty_values = self.data["novelty_index"].to_list()
        uncertainty_values = self.data["placement_uncertainty"].to_list()

        fig = make_subplots(
            rows=1,
            cols=2,
            subplot_titles=("Novelty Index Distribution", "Placement Uncertainty Distribution"),
            horizontal_spacing=0.1,
        )

        # Novelty histogram
        fig.add_trace(
            go.Histogram(
                x=novelty_values,
                nbinsx=self.nbins,
                name="Novelty",
                marker_color="#667eea",
                opacity=0.75,
                showlegend=False,
            ),
            row=1,
            col=1,
        )

        # Uncertainty histogram
        fig.add_trace(
            go.Histogram(
                x=uncertainty_values,
                nbinsx=self.nbins,
                name="Uncertainty",
                marker_color="#764ba2",
                opacity=0.75,
                showlegend=False,
            ),
            row=1,
            col=2,
        )

        # Add threshold lines
        if self.show_thresholds:
            # Novelty thresholds
            for value in [2.0, 5.0, 15.0, 25.0]:
                fig.add_vline(
                    x=value,
                    line_dash="dash",
                    line_color="#888",
                    opacity=0.5,
                    row=1,
                    col=1,
                )

            # Uncertainty thresholds
            for value in [0.5, 2.0, 5.0]:
                fig.add_vline(
                    x=value,
                    line_dash="dash",
                    line_color="#888",
                    opacity=0.5,
                    row=1,
                    col=2,
                )

        # Layout
        fig.update_layout(
            **self.config.to_layout_dict("Distribution of Classification Metrics"),
        )

        fig.update_xaxes(title_text="Novelty Index", row=1, col=1)
        fig.update_xaxes(title_text="Placement Uncertainty", row=1, col=2)
        fig.update_yaxes(title_text="Count", row=1, col=1)

        return self._apply_config(fig)
