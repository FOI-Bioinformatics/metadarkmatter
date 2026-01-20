"""
Multi-sample comparison plots.

Visualizations for comparing classification results across multiple samples,
enabling cross-sample analysis of microbial diversity patterns.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import plotly.graph_objects as go

from metadarkmatter.visualization.plots.base import (
    SEQUENTIAL_PALETTE,
    TAXONOMY_COLORS,
    BasePlot,
    PlotConfig,
)

if TYPE_CHECKING:
    import polars as pl


class MultiSampleBarChart(BasePlot):
    """
    Grouped or stacked bar chart comparing classification proportions across samples.

    Shows the distribution of taxonomic classifications for each sample,
    enabling visual comparison of diversity patterns.
    """

    def __init__(
        self,
        summaries: dict[str, dict[str, Any]],
        config: PlotConfig | None = None,
        stacked: bool = False,
        normalize: bool = True,
        title: str = "Classification Comparison",
    ) -> None:
        """
        Initialize multi-sample bar chart.

        Args:
            summaries: Dictionary mapping sample names to their summary dicts.
                Each summary should contain: known_species, novel_species,
                novel_genus, conserved_regions, total_reads
            config: Plot configuration
            stacked: If True, create stacked bars; if False, grouped bars
            normalize: If True, show percentages; if False, show counts
            title: Chart title
        """
        super().__init__(config)
        self.summaries = summaries
        self.stacked = stacked
        self.normalize = normalize
        self.title = title

    def create_figure(self) -> go.Figure:
        """Create multi-sample comparison bar chart."""
        sample_names = list(self.summaries.keys())
        categories = ["Known Species", "Novel Species", "Novel Genus", "Conserved Region"]

        fig = go.Figure()

        for category in categories:
            key = category.lower().replace(" ", "_")
            if key == "conserved_region":
                key = "conserved_regions"

            values = []
            for sample in sample_names:
                summary = self.summaries[sample]
                count = summary.get(key, 0)
                total = summary.get("total_reads", 1)

                if self.normalize:
                    values.append(count / total * 100 if total > 0 else 0)
                else:
                    values.append(count)

            fig.add_trace(
                go.Bar(
                    name=category,
                    x=sample_names,
                    y=values,
                    marker_color=TAXONOMY_COLORS[category],
                    text=[f"{v:.1f}%" if self.normalize else f"{v:,}" for v in values],
                    textposition="inside" if self.stacked else "outside",
                    hovertemplate=(
                        f"<b>{category}</b><br>"
                        "Sample: %{x}<br>"
                        f"{'Percentage' if self.normalize else 'Count'}: "
                        f"%{{y:.1f}}{'%' if self.normalize else ''}<extra></extra>"
                    ),
                )
            )

        barmode = "stack" if self.stacked else "group"
        yaxis_title = "Percentage of Reads" if self.normalize else "Number of Reads"

        fig.update_layout(
            barmode=barmode,
            xaxis_title="Sample",
            yaxis_title=yaxis_title,
            legend={
                "orientation": "h",
                "yanchor": "bottom",
                "y": 1.02,
                "xanchor": "center",
                "x": 0.5,
            },
            **self.config.to_layout_dict(self.title),
        )

        return self._apply_config(fig)


class MultiSampleNoveltyComparison(BasePlot):
    """
    Box plot comparing novelty index distributions across samples.

    Shows the spread and central tendency of novelty for each sample,
    highlighting differences in divergence from reference genomes.
    """

    def __init__(
        self,
        sample_data: dict[str, pl.DataFrame],
        config: PlotConfig | None = None,
        title: str = "Novelty Index Comparison",
    ) -> None:
        """
        Initialize novelty comparison plot.

        Args:
            sample_data: Dictionary mapping sample names to DataFrames.
                Each DataFrame must have 'novelty_index' column.
            config: Plot configuration
            title: Chart title
        """
        super().__init__(config)
        self.sample_data = sample_data
        self.title = title

    def create_figure(self) -> go.Figure:
        """Create novelty index comparison box plot."""
        fig = go.Figure()

        for i, (sample_name, df) in enumerate(self.sample_data.items()):
            novelty_values = df["novelty_index"].to_list()

            fig.add_trace(
                go.Box(
                    y=novelty_values,
                    name=sample_name,
                    marker_color=SEQUENTIAL_PALETTE[i % len(SEQUENTIAL_PALETTE)],
                    boxpoints="outliers",
                    hoverinfo="y+name",
                )
            )

        fig.update_layout(
            yaxis_title="Novelty Index (100 - % Identity)",
            showlegend=False,
            **self.config.to_layout_dict(self.title),
        )

        return self._apply_config(fig)


class MultiSampleScatterMatrix(BasePlot):
    """
    Scatter plot matrix showing pairwise sample comparisons.

    Compares mean novelty vs mean uncertainty for each sample,
    positioned as points in the 2D diversity space.
    """

    def __init__(
        self,
        summaries: dict[str, dict[str, Any]],
        config: PlotConfig | None = None,
        title: str = "Sample Diversity Overview",
    ) -> None:
        """
        Initialize sample scatter plot.

        Args:
            summaries: Dictionary mapping sample names to summary dicts.
                Each summary should contain: mean_novelty_index,
                mean_placement_uncertainty, total_reads
            config: Plot configuration
            title: Chart title
        """
        super().__init__(config)
        self.summaries = summaries
        self.title = title

    def create_figure(self) -> go.Figure:
        """Create sample diversity scatter plot."""
        sample_names = list(self.summaries.keys())

        x_values = [self.summaries[s].get("mean_novelty_index", 0) for s in sample_names]
        y_values = [self.summaries[s].get("mean_placement_uncertainty", 0) for s in sample_names]
        sizes = [self.summaries[s].get("total_reads", 1) for s in sample_names]

        # Normalize sizes for display
        max_size = max(sizes) if sizes else 1
        marker_sizes = [max(10, 50 * (s / max_size)) for s in sizes]

        fig = go.Figure()

        fig.add_trace(
            go.Scatter(
                x=x_values,
                y=y_values,
                mode="markers+text",
                text=sample_names,
                textposition="top center",
                marker={
                    "size": marker_sizes,
                    "color": [SEQUENTIAL_PALETTE[i % len(SEQUENTIAL_PALETTE)]
                              for i in range(len(sample_names))],
                    "opacity": 0.7,
                    "line": {"width": 1, "color": "#333"},
                },
                hovertemplate=(
                    "<b>%{text}</b><br>"
                    "Mean Novelty: %{x:.2f}<br>"
                    "Mean Uncertainty: %{y:.2f}<br>"
                    "<extra></extra>"
                ),
            )
        )

        fig.update_layout(
            xaxis_title="Mean Novelty Index",
            yaxis_title="Mean Placement Uncertainty",
            **self.config.to_layout_dict(self.title),
        )

        return self._apply_config(fig)


class MultiSampleHeatmap(BasePlot):
    """
    Heatmap showing classification proportions across samples.

    Provides a compact overview of diversity patterns, with samples
    as rows and classification categories as columns.
    """

    def __init__(
        self,
        summaries: dict[str, dict[str, Any]],
        config: PlotConfig | None = None,
        title: str = "Sample Classification Heatmap",
    ) -> None:
        """
        Initialize classification heatmap.

        Args:
            summaries: Dictionary mapping sample names to summary dicts.
            config: Plot configuration
            title: Chart title
        """
        super().__init__(config)
        self.summaries = summaries
        self.title = title

    def create_figure(self) -> go.Figure:
        """Create sample classification heatmap."""
        sample_names = list(self.summaries.keys())
        categories = ["Known Species", "Novel Species", "Novel Genus", "Conserved Region"]

        # Build matrix of percentages
        z_data = []
        for sample in sample_names:
            summary = self.summaries[sample]
            total = summary.get("total_reads", 1)
            row = []
            for category in categories:
                key = category.lower().replace(" ", "_")
                if key == "conserved_region":
                    key = "conserved_regions"
                count = summary.get(key, 0)
                pct = count / total * 100 if total > 0 else 0
                row.append(pct)
            z_data.append(row)

        fig = go.Figure()

        fig.add_trace(
            go.Heatmap(
                z=z_data,
                x=categories,
                y=sample_names,
                colorscale="RdYlGn_r",
                zmin=0,
                zmax=100,
                text=[[f"{v:.1f}%" for v in row] for row in z_data],
                texttemplate="%{text}",
                textfont={"size": 10},
                hovertemplate=(
                    "Sample: %{y}<br>"
                    "Category: %{x}<br>"
                    "Percentage: %{z:.1f}%<extra></extra>"
                ),
                colorbar={"title": "% of Reads"},
            )
        )

        fig.update_layout(
            xaxis_title="Classification",
            yaxis_title="Sample",
            **self.config.to_layout_dict(self.title),
        )

        return self._apply_config(fig)


class MultiSampleTimeSeries(BasePlot):
    """
    Line plot showing diversity metrics over time or sample order.

    Useful for tracking changes in microbial diversity across
    a temporal series or sampling gradient.
    """

    def __init__(
        self,
        summaries: dict[str, dict[str, Any]],
        config: PlotConfig | None = None,
        metric: str = "novel_percentage",
        title: str = "Diversity Trend",
    ) -> None:
        """
        Initialize time series plot.

        Args:
            summaries: Ordered dictionary mapping sample names to summaries.
            config: Plot configuration
            metric: Which metric to plot. Options:
                - "novel_percentage": Combined novel species + genus
                - "mean_novelty": Mean novelty index
                - "mean_uncertainty": Mean placement uncertainty
            title: Chart title
        """
        super().__init__(config)
        self.summaries = summaries
        self.metric = metric
        self.title = title

    def create_figure(self) -> go.Figure:
        """Create diversity time series plot."""
        sample_names = list(self.summaries.keys())

        if self.metric == "novel_percentage":
            values = []
            for s in sample_names:
                summary = self.summaries[s]
                total = summary.get("total_reads", 1)
                novel = summary.get("novel_species", 0) + summary.get("novel_genus", 0)
                values.append(novel / total * 100 if total > 0 else 0)
            y_title = "Novel Diversity (%)"
        elif self.metric == "mean_novelty":
            values = [self.summaries[s].get("mean_novelty_index", 0) for s in sample_names]
            y_title = "Mean Novelty Index"
        else:  # mean_uncertainty
            values = [self.summaries[s].get("mean_placement_uncertainty", 0) for s in sample_names]
            y_title = "Mean Placement Uncertainty"

        fig = go.Figure()

        fig.add_trace(
            go.Scatter(
                x=sample_names,
                y=values,
                mode="lines+markers",
                line={"color": "#667eea", "width": 2},
                marker={"size": 8, "color": "#667eea"},
                hovertemplate=(
                    "<b>%{x}</b><br>"
                    f"{y_title}: %{{y:.2f}}<extra></extra>"
                ),
            )
        )

        fig.update_layout(
            xaxis_title="Sample",
            yaxis_title=y_title,
            **self.config.to_layout_dict(self.title),
        )

        return self._apply_config(fig)
