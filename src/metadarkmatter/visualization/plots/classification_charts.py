"""
Classification summary charts: donut charts, bar charts, and summary visualizations.
"""

from __future__ import annotations

from typing import Any

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from metadarkmatter.visualization.plots.base import (
    DIVERSITY_COLORS,
    TAXONOMY_COLORS,
    BasePlot,
    PlotConfig,
    format_count,
)


class ClassificationDonutChart(BasePlot):
    """
    Donut chart showing classification category distribution.

    Displays the proportion of reads in each taxonomic classification
    category with counts and percentages.
    """

    def __init__(
        self,
        summary: dict[str, Any],
        config: PlotConfig | None = None,
        title: str = "Classification Distribution",
        hole_size: float = 0.5,
    ) -> None:
        """
        Initialize donut chart.

        Args:
            summary: TaxonomicSummary dictionary with counts
            config: Plot configuration
            title: Chart title
            hole_size: Size of center hole (0-1)
        """
        super().__init__(config)
        self.summary = summary
        self.title = title
        self.hole_size = hole_size

    def create_figure(self) -> go.Figure:
        """Create donut chart with classification proportions."""
        # Extract values in consistent order
        categories = [
            "Known Species",
            "Novel Species",
            "Novel Genus",
            "Species Boundary",
            "Ambiguous",
            "Ambiguous Within Genus",
            "Conserved Region",
            "Unclassified",
        ]
        values = [
            self.summary.get("known_species", 0),
            self.summary.get("novel_species", 0),
            self.summary.get("novel_genus", 0),
            self.summary.get("species_boundary", 0),
            self.summary.get("ambiguous", 0),
            self.summary.get("ambiguous_within_genus", 0),
            self.summary.get("conserved_regions", 0),
            self.summary.get("unclassified", 0),
        ]
        colors = [TAXONOMY_COLORS[cat] for cat in categories]

        total = sum(values)

        # Create custom labels with counts and percentages
        labels = []
        for cat, val in zip(categories, values, strict=True):
            pct = (val / total * 100) if total > 0 else 0
            labels.append(f"{cat}<br>{format_count(val)} ({pct:.1f}%)")

        fig = go.Figure()

        # Pull out novel categories slightly for emphasis
        pull = [0, 0.02, 0.02, 0, 0, 0, 0, 0]

        fig.add_trace(
            go.Pie(
                values=values,
                labels=categories,
                hole=self.hole_size,
                marker={"colors": colors},
                textinfo="percent",
                textposition="outside",
                hovertemplate=(
                    "<b>%{label}</b><br>"
                    "Count: %{value:,}<br>"
                    "Percentage: %{percent}<extra></extra>"
                ),
                pull=pull,
            )
        )

        # Add center annotation with total
        fig.add_annotation(
            text=f"<b>{format_count(total)}</b><br>Total Reads",
            x=0.5,
            y=0.5,
            font={"size": 14},
            showarrow=False,
        )

        fig.update_layout(
            showlegend=True,
            legend={
                "orientation": "h",
                "yanchor": "bottom",
                "y": -0.15,
                "xanchor": "center",
                "x": 0.5,
            },
            **self.config.to_layout_dict(self.title),
        )

        return self._apply_config(fig)


class ClassificationBarChart(BasePlot):
    """
    Horizontal bar chart showing classification counts.

    Provides a clear view of absolute numbers in each category
    with percentage labels.
    """

    def __init__(
        self,
        summary: dict[str, Any],
        config: PlotConfig | None = None,
        title: str = "Classification Counts",
        orientation: str = "h",
    ) -> None:
        """
        Initialize bar chart.

        Args:
            summary: TaxonomicSummary dictionary with counts
            config: Plot configuration
            title: Chart title
            orientation: 'h' for horizontal, 'v' for vertical
        """
        super().__init__(config)
        self.summary = summary
        self.title = title
        self.orientation = orientation

    def create_figure(self) -> go.Figure:
        """Create bar chart with classification counts."""
        categories = [
            "Known Species",
            "Novel Species",
            "Novel Genus",
            "Species Boundary",
            "Ambiguous",
            "Ambiguous Within Genus",
            "Conserved Region",
            "Unclassified",
        ]
        values = [
            self.summary.get("known_species", 0),
            self.summary.get("novel_species", 0),
            self.summary.get("novel_genus", 0),
            self.summary.get("species_boundary", 0),
            self.summary.get("ambiguous", 0),
            self.summary.get("ambiguous_within_genus", 0),
            self.summary.get("conserved_regions", 0),
            self.summary.get("unclassified", 0),
        ]
        colors = [TAXONOMY_COLORS[cat] for cat in categories]

        total = sum(values)

        fig = go.Figure()

        if self.orientation == "h":
            # Reverse for horizontal (so Known Species is at top)
            fig.add_trace(
                go.Bar(
                    y=categories[::-1],
                    x=values[::-1],
                    orientation="h",
                    marker_color=colors[::-1],
                    text=[
                        f"{v:,} ({v/total*100:.1f}%)" if total > 0 else "0"
                        for v in values[::-1]
                    ],
                    textposition="outside",
                    hovertemplate=(
                        "<b>%{y}</b><br>"
                        "Count: %{x:,}<extra></extra>"
                    ),
                )
            )
            fig.update_layout(xaxis_title="Number of Reads")
        else:
            fig.add_trace(
                go.Bar(
                    x=categories,
                    y=values,
                    marker_color=colors,
                    text=[
                        f"{v:,}<br>({v/total*100:.1f}%)" if total > 0 else "0"
                        for v in values
                    ],
                    textposition="outside",
                    hovertemplate=(
                        "<b>%{x}</b><br>"
                        "Count: %{y:,}<extra></extra>"
                    ),
                )
            )
            fig.update_layout(yaxis_title="Number of Reads")

        fig.update_layout(
            showlegend=False,
            **self.config.to_layout_dict(self.title),
        )

        return self._apply_config(fig)


class ClassificationSummaryPlot(BasePlot):
    """
    Combined summary visualization with donut chart and statistics.

    Provides an overview of classification results in a single figure.
    """

    def __init__(
        self,
        summary: dict[str, Any],
        config: PlotConfig | None = None,
        title: str = "Classification Summary",
    ) -> None:
        """
        Initialize summary plot.

        Args:
            summary: TaxonomicSummary dictionary
            config: Plot configuration
            title: Chart title
        """
        super().__init__(config)
        self.summary = summary
        self.title = title

    def create_figure(self) -> go.Figure:
        """Create combined summary with donut and bar charts."""
        fig = make_subplots(
            rows=1,
            cols=2,
            specs=[[{"type": "pie"}, {"type": "bar"}]],
            subplot_titles=("Distribution", "Counts"),
            horizontal_spacing=0.15,
        )

        categories = [
            "Known Species",
            "Novel Species",
            "Novel Genus",
            "Species Boundary",
            "Ambiguous",
            "Ambiguous Within Genus",
            "Conserved Region",
            "Unclassified",
        ]
        values = [
            self.summary.get("known_species", 0),
            self.summary.get("novel_species", 0),
            self.summary.get("novel_genus", 0),
            self.summary.get("species_boundary", 0),
            self.summary.get("ambiguous", 0),
            self.summary.get("ambiguous_within_genus", 0),
            self.summary.get("conserved_regions", 0),
            self.summary.get("unclassified", 0),
        ]
        colors = [TAXONOMY_COLORS[cat] for cat in categories]

        # Donut chart
        fig.add_trace(
            go.Pie(
                values=values,
                labels=categories,
                hole=0.4,
                marker={"colors": colors},
                textinfo="percent",
                showlegend=True,
                domain={"x": [0, 0.45]},
            ),
            row=1,
            col=1,
        )

        # Bar chart
        fig.add_trace(
            go.Bar(
                y=categories[::-1],
                x=values[::-1],
                orientation="h",
                marker_color=colors[::-1],
                showlegend=False,
            ),
            row=1,
            col=2,
        )

        fig.update_layout(
            **self.config.to_layout_dict(self.title),
        )

        return self._apply_config(fig)


class NovelDiversityGauge(BasePlot):
    """
    Gauge chart showing novel diversity percentage.

    Highlights the percentage of reads classified as novel
    (Novel Species + Novel Genus).
    """

    def __init__(
        self,
        summary: dict[str, Any],
        config: PlotConfig | None = None,
        title: str = "Novel Diversity",
    ) -> None:
        """
        Initialize gauge chart.

        Args:
            summary: TaxonomicSummary dictionary
            config: Plot configuration
            title: Chart title
        """
        super().__init__(config)
        self.summary = summary
        self.title = title

    def create_figure(self) -> go.Figure:
        """Create gauge showing novel diversity percentage."""
        novel_species = self.summary.get("novel_species", 0)
        novel_genus = self.summary.get("novel_genus", 0)
        total = self.summary.get("total_reads", 1)

        novel_pct = (novel_species + novel_genus) / total * 100 if total > 0 else 0

        fig = go.Figure()

        fig.add_trace(
            go.Indicator(
                mode="gauge+number+delta",
                value=novel_pct,
                number={"suffix": "%", "font": {"size": 32}},
                title={"text": self.title, "font": {"size": 16}},
                gauge={
                    "axis": {"range": [0, 100], "ticksuffix": "%"},
                    "bar": {"color": TAXONOMY_COLORS["Novel Species"]},
                    "steps": [
                        {"range": [0, 25], "color": "#e8f5e9"},
                        {"range": [25, 50], "color": "#fff3e0"},
                        {"range": [50, 75], "color": "#ffecb3"},
                        {"range": [75, 100], "color": "#ffcdd2"},
                    ],
                    "threshold": {
                        "line": {"color": "red", "width": 2},
                        "thickness": 0.75,
                        "value": 50,
                    },
                },
            )
        )

        fig.update_layout(
            **self.config.to_layout_dict(),
        )

        return self._apply_config(fig)


class ClassificationMetricsCards(BasePlot):
    """
    Visual representation of key metrics as card-style indicators.
    """

    def __init__(
        self,
        summary: dict[str, Any],
        config: PlotConfig | None = None,
    ) -> None:
        """
        Initialize metrics cards.

        Args:
            summary: TaxonomicSummary dictionary
            config: Plot configuration
        """
        super().__init__(config)
        self.summary = summary

    def create_figure(self) -> go.Figure:
        """Create indicator cards for key metrics."""
        fig = make_subplots(
            rows=1,
            cols=4,
            specs=[[{"type": "indicator"}] * 4],
            subplot_titles=(
                "Total Reads",
                "Novel Diversity",
                "Mean Novelty",
                "Mean Uncertainty",
            ),
        )

        total = self.summary.get("total_reads", 0)
        novel = self.summary.get("novel_species", 0) + self.summary.get("novel_genus", 0)
        novel_pct = (novel / total * 100) if total > 0 else 0

        # Total reads
        fig.add_trace(
            go.Indicator(
                mode="number",
                value=total,
                number={"font": {"size": 36}},
            ),
            row=1,
            col=1,
        )

        # Novel diversity
        fig.add_trace(
            go.Indicator(
                mode="number",
                value=novel_pct,
                number={"suffix": "%", "font": {"size": 36}},
            ),
            row=1,
            col=2,
        )

        # Mean novelty
        fig.add_trace(
            go.Indicator(
                mode="number",
                value=self.summary.get("mean_novelty_index", 0),
                number={"font": {"size": 36}},
            ),
            row=1,
            col=3,
        )

        # Mean uncertainty
        fig.add_trace(
            go.Indicator(
                mode="number",
                value=self.summary.get("mean_placement_uncertainty", 0),
                number={"font": {"size": 36}},
            ),
            row=1,
            col=4,
        )

        fig.update_layout(
            height=200,
            **self.config.to_layout_dict(),
        )

        return self._apply_config(fig)


class DiversityDonutChart(BasePlot):
    """
    Donut chart showing high-level diversity distribution.

    Groups reads into three categories:
    - Known: Confident match to known diversity (Known Species)
    - Novel: Confident novel diversity (Novel Species, Novel Genus)
    - Uncertain: Cannot confidently classify (all other categories)
    """

    def __init__(
        self,
        summary: dict[str, Any],
        config: PlotConfig | None = None,
        title: str = "Diversity Summary",
        hole_size: float = 0.5,
    ) -> None:
        """
        Initialize diversity donut chart.

        Args:
            summary: TaxonomicSummary dictionary with diversity counts
            config: Plot configuration
            title: Chart title
            hole_size: Size of center hole (0-1)
        """
        super().__init__(config)
        self.summary = summary
        self.title = title
        self.hole_size = hole_size

    def create_figure(self) -> go.Figure:
        """Create donut chart with diversity proportions."""
        categories = ["Known", "Novel", "Uncertain"]
        values = [
            self.summary.get("diversity_known", 0),
            self.summary.get("diversity_novel", 0),
            self.summary.get("diversity_uncertain", 0),
        ]
        colors = [DIVERSITY_COLORS[cat] for cat in categories]

        total = sum(values)

        fig = go.Figure()

        # Pull out novel category slightly for emphasis
        pull = [0, 0.03, 0]

        fig.add_trace(
            go.Pie(
                values=values,
                labels=categories,
                hole=self.hole_size,
                marker={"colors": colors},
                textinfo="percent+label",
                textposition="outside",
                hovertemplate=(
                    "<b>%{label}</b><br>"
                    "Count: %{value:,}<br>"
                    "Percentage: %{percent}<extra></extra>"
                ),
                pull=pull,
            )
        )

        # Add center annotation with total
        fig.add_annotation(
            text=f"<b>{format_count(total)}</b><br>Total Reads",
            x=0.5,
            y=0.5,
            font={"size": 14},
            showarrow=False,
        )

        fig.update_layout(
            showlegend=True,
            legend={
                "orientation": "h",
                "yanchor": "bottom",
                "y": -0.15,
                "xanchor": "center",
                "x": 0.5,
            },
            **self.config.to_layout_dict(self.title),
        )

        return self._apply_config(fig)


class DiversitySunburstChart(BasePlot):
    """
    Sunburst chart showing diversity hierarchy.

    Inner ring shows high-level diversity status (Known/Novel/Uncertain).
    Outer ring shows detailed taxonomic categories grouped under their parent.
    This makes it easy to compare the two levels visually.
    """

    # Define the hierarchy: diversity_status -> taxonomic categories
    DIVERSITY_HIERARCHY: dict[str, list[str]] = {
        "Known": ["Known Species"],
        "Novel": ["Novel Species", "Novel Genus"],
        "Uncertain": [
            "Species Boundary",
            "Ambiguous",
            "Ambiguous Within Genus",
            "Conserved Region",
            "Unclassified",
        ],
    }

    # Map taxonomic categories to summary dict keys
    CATEGORY_TO_KEY: dict[str, str] = {
        "Known Species": "known_species",
        "Novel Species": "novel_species",
        "Novel Genus": "novel_genus",
        "Species Boundary": "species_boundary",
        "Ambiguous": "ambiguous",
        "Ambiguous Within Genus": "ambiguous_within_genus",
        "Conserved Region": "conserved_regions",
        "Unclassified": "unclassified",
    }

    def __init__(
        self,
        summary: dict[str, Any],
        config: PlotConfig | None = None,
        title: str = "Diversity Classification",
    ) -> None:
        """
        Initialize sunburst chart.

        Args:
            summary: TaxonomicSummary dictionary with counts
            config: Plot configuration
            title: Chart title
        """
        super().__init__(config)
        self.summary = summary
        self.title = title

    def create_figure(self) -> go.Figure:
        """Create sunburst chart showing diversity hierarchy."""
        # Build hierarchical data for sunburst
        ids = []
        labels = []
        parents = []
        values = []
        colors = []

        total = self.summary.get("total_reads", 0)

        # Add root (empty string parent for top level)
        # We don't add an explicit root - sunburst handles this

        # Add diversity status level (inner ring)
        for diversity_status in ["Known", "Novel", "Uncertain"]:
            ids.append(diversity_status)
            labels.append(diversity_status)
            parents.append("")  # Top level has no parent
            # Value is sum of child categories
            child_sum = sum(
                self.summary.get(self.CATEGORY_TO_KEY[cat], 0)
                for cat in self.DIVERSITY_HIERARCHY[diversity_status]
            )
            values.append(child_sum)
            colors.append(DIVERSITY_COLORS[diversity_status])

        # Add taxonomic category level (outer ring)
        for diversity_status, categories in self.DIVERSITY_HIERARCHY.items():
            base_color = DIVERSITY_COLORS[diversity_status]
            for i, category in enumerate(categories):
                cat_id = f"{diversity_status}-{category}"
                ids.append(cat_id)
                labels.append(category)
                parents.append(diversity_status)
                values.append(self.summary.get(self.CATEGORY_TO_KEY[category], 0))
                # Use taxonomy color for detailed view
                colors.append(TAXONOMY_COLORS.get(category, base_color))

        fig = go.Figure()

        fig.add_trace(
            go.Sunburst(
                ids=ids,
                labels=labels,
                parents=parents,
                values=values,
                branchvalues="total",
                marker={"colors": colors},
                hovertemplate=(
                    "<b>%{label}</b><br>"
                    "Count: %{value:,}<br>"
                    "Percentage: %{percentRoot:.1%}<extra></extra>"
                ),
                textinfo="label+percent root",
                insidetextorientation="horizontal",
            )
        )

        # Add center annotation with total
        fig.add_annotation(
            text=f"<b>{format_count(total)}</b><br>Total",
            x=0.5,
            y=0.5,
            font={"size": 12},
            showarrow=False,
        )

        fig.update_layout(
            **self.config.to_layout_dict(self.title),
        )

        return self._apply_config(fig)
