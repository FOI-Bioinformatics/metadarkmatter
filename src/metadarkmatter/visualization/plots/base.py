"""
Base classes and utilities for plot generation.

Defines color palettes, common styling, and abstract interfaces
for all visualization components.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import plotly.graph_objects as go

if TYPE_CHECKING:
    import polars as pl


# =============================================================================
# Color Palettes
# =============================================================================

# Taxonomic classification colors (consistent with methodology)
TAXONOMY_COLORS: dict[str, str] = {
    "Known Species": "#2ecc71",      # Green
    "Novel Species": "#f39c12",      # Orange
    "Novel Genus": "#e74c3c",        # Red
    "Species Boundary": "#3498db",   # Blue - at species boundary zone
    "Ambiguous": "#1abc9c",          # Teal - high U within genus
    "Ambiguous Within Genus": "#9b59b6",  # Purple
    "Conserved Region": "#95a5a6",   # Gray
    "Unclassified": "#7f8c8d",       # Dark gray
}

# High-level diversity status colors (for summary views)
DIVERSITY_COLORS: dict[str, str] = {
    "Known": "#2ecc71",       # Green - confident known diversity
    "Novel": "#e74c3c",       # Red - confident novel diversity
    "Uncertain": "#95a5a6",   # Gray - uncertain classification
}

# Bayesian 6-category colors used in report generator charts
BAYESIAN_CATEGORY_COLORS: dict[str, str] = {
    "Known Species": "#22c55e",
    "Novel Species": "#f59e0b",
    "Novel Genus": "#ef4444",
    "Species Boundary": "#a855f7",
    "Ambiguous": "#94a3b8",
    "Unclassified": "#64748b",
}

# Sequential palette for categorical data
SEQUENTIAL_PALETTE: list[str] = [
    "#1f77b4",  # Blue
    "#ff7f0e",  # Orange
    "#2ca02c",  # Green
    "#d62728",  # Red
    "#9467bd",  # Purple
    "#8c564b",  # Brown
    "#e377c2",  # Pink
    "#7f7f7f",  # Gray
    "#bcbd22",  # Yellow-green
    "#17becf",  # Cyan
]

# ANI heatmap colorscale (75-100% range)
ANI_COLORSCALE: list[list[Any]] = [
    [0.0, "#d73027"],    # 75% - Red (distant)
    [0.25, "#fc8d59"],   # 81.25%
    [0.5, "#fee090"],    # 87.5% - Yellow (species boundary)
    [0.75, "#91bfdb"],   # 93.75%
    [1.0, "#4575b4"],    # 100% - Blue (identical)
]

# Identity threshold colors for recruitment plots
IDENTITY_BAND_COLORS: dict[str, str] = {
    "known": "#2ecc71",    # Green (98-100%)
    "novel_species": "#f39c12",  # Orange (85-98%)
    "novel_genus": "#e74c3c",    # Red (75-85%)
}


# =============================================================================
# Configuration
# =============================================================================

@dataclass
class PlotConfig:
    """Common configuration for all plot types."""

    width: int = 800
    height: int = 500
    template: str = "plotly_white"
    font_family: str = "Arial, Helvetica, sans-serif"
    title_font_size: int = 16
    axis_font_size: int = 12
    legend_font_size: int = 11
    margin: dict[str, int] = field(
        default_factory=lambda: {"l": 60, "r": 40, "t": 60, "b": 60}
    )
    show_legend: bool = True

    def to_layout_dict(
        self,
        title: str | None = None,
        include_legend: bool = False,
    ) -> dict[str, Any]:
        """Convert config to Plotly layout dictionary.

        Args:
            title: Optional title text to include in layout.
            include_legend: Whether to include legend settings (default False
                to avoid conflicts when plots set legend explicitly).
        """
        layout: dict[str, Any] = {
            "template": self.template,
            "font": {
                "family": self.font_family,
                "size": self.axis_font_size,
            },
            "margin": self.margin,
        }
        if title:
            layout["title"] = {"text": title, "font": {"size": self.title_font_size}}
        if include_legend:
            layout["legend"] = {"font": {"size": self.legend_font_size}}
            layout["showlegend"] = self.show_legend
        return layout


@dataclass
class ThresholdConfig:
    """Classification threshold configuration."""

    # Novelty thresholds (must match ScoringConfig defaults)
    novelty_known_max: float = 4.0
    novelty_novel_species_min: float = 4.0
    novelty_novel_species_max: float = 20.0
    novelty_novel_genus_min: float = 20.0
    novelty_novel_genus_max: float = 25.0

    # Uncertainty thresholds (must match ScoringConfig defaults)
    uncertainty_known_max: float = 1.5
    uncertainty_novel_species_max: float = 1.5
    uncertainty_novel_genus_max: float = 1.5
    uncertainty_conserved_min: float = 5.0

    # Identity thresholds (for recruitment plots)
    identity_known_min: float = 98.0
    identity_novel_species_min: float = 80.0
    identity_novel_genus_min: float = 75.0


# =============================================================================
# Base Plot Class
# =============================================================================

class BasePlot(ABC):
    """Abstract base class for all plot generators."""

    def __init__(
        self,
        config: PlotConfig | None = None,
        thresholds: ThresholdConfig | None = None,
    ) -> None:
        """
        Initialize plot generator.

        Args:
            config: Plot configuration (dimensions, styling)
            thresholds: Classification thresholds for annotations
        """
        self.config = config or PlotConfig()
        self.thresholds = thresholds or ThresholdConfig()

    @abstractmethod
    def create_figure(self) -> go.Figure:
        """Create and return the Plotly figure."""
        ...

    def to_html_div(self, include_plotlyjs: bool = False) -> str:
        """
        Export plot as HTML div for embedding in reports.

        Args:
            include_plotlyjs: Whether to include Plotly.js library

        Returns:
            HTML string containing the plot div
        """
        fig = self.create_figure()
        return fig.to_html(
            full_html=False,
            include_plotlyjs="cdn" if include_plotlyjs else False,
            div_id=self._get_div_id(),
        )

    def to_html(self, include_plotlyjs: bool = True) -> str:
        """
        Export plot as standalone HTML file content.

        Args:
            include_plotlyjs: Whether to embed Plotly.js library

        Returns:
            Complete HTML document string
        """
        fig = self.create_figure()
        return fig.to_html(
            full_html=True,
            include_plotlyjs=True if include_plotlyjs else "cdn",
        )

    def to_json(self) -> str:
        """Export plot as JSON for data interchange."""
        fig = self.create_figure()
        return fig.to_json()

    def save(self, path: str, **kwargs: Any) -> None:
        """
        Save plot to file.

        Args:
            path: Output file path (.html, .png, .json)
            **kwargs: Additional arguments passed to write method
        """
        fig = self.create_figure()

        if path.endswith(".html"):
            fig.write_html(path, include_plotlyjs=True, **kwargs)
        elif path.endswith(".json"):
            fig.write_json(path, **kwargs)
        elif path.endswith((".png", ".jpg", ".jpeg", ".svg", ".pdf")):
            fig.write_image(path, **kwargs)
        else:
            msg = f"Unsupported file format: {path}"
            raise ValueError(msg)

    def _get_div_id(self) -> str:
        """Generate unique div ID for the plot."""
        return f"plot-{self.__class__.__name__.lower()}"

    def _apply_config(self, fig: go.Figure) -> go.Figure:
        """Apply common configuration to figure."""
        fig.update_layout(
            **self.config.to_layout_dict(),
            width=self.config.width,
            height=self.config.height,
        )
        return fig

    def _add_threshold_line(
        self,
        fig: go.Figure,
        value: float,
        axis: str = "x",
        color: str = "#888888",
        dash: str = "dash",
        annotation: str | None = None,
    ) -> None:
        """
        Add vertical or horizontal threshold line to figure.

        Args:
            fig: Plotly figure to modify
            value: Threshold value
            axis: 'x' for vertical line, 'y' for horizontal line
            color: Line color
            dash: Line style ('solid', 'dash', 'dot', 'dashdot')
            annotation: Optional text annotation for the line
        """
        if axis == "x":
            fig.add_vline(
                x=value,
                line_color=color,
                line_dash=dash,
                opacity=0.7,
            )
            if annotation:
                fig.add_annotation(
                    x=value,
                    y=1.02,
                    yref="paper",
                    text=annotation,
                    showarrow=False,
                    font={"size": 10, "color": color},
                )
        else:
            fig.add_hline(
                y=value,
                line_color=color,
                line_dash=dash,
                opacity=0.7,
            )
            if annotation:
                fig.add_annotation(
                    x=1.02,
                    xref="paper",
                    y=value,
                    text=annotation,
                    showarrow=False,
                    font={"size": 10, "color": color},
                )


# =============================================================================
# Utility Functions
# =============================================================================

def get_taxonomy_color(classification: str) -> str:
    """Get color for a taxonomic classification category."""
    return TAXONOMY_COLORS.get(classification, "#808080")


def subsample_dataframe(
    df: pl.DataFrame,
    max_points: int,
    seed: int = 42,
) -> pl.DataFrame:
    """
    Subsample DataFrame if it exceeds max_points.

    Args:
        df: Input DataFrame
        max_points: Maximum number of rows to keep
        seed: Random seed for reproducibility

    Returns:
        Subsampled DataFrame (or original if under limit)
    """
    if len(df) <= max_points:
        return df
    return df.sample(n=max_points, seed=seed)


def format_count(count: int) -> str:
    """Format large counts with K/M suffixes."""
    if count >= 1_000_000:
        return f"{count / 1_000_000:.1f}M"
    if count >= 1_000:
        return f"{count / 1_000:.1f}K"
    return str(count)


def format_percentage(value: float, decimals: int = 1) -> str:
    """Format value as percentage string."""
    return f"{value:.{decimals}f}%"
