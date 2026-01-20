"""
Recruitment plot visualization using Plotly.

Creates interactive scatter plots showing read recruitment patterns
across reference genomes, with horizontal bands indicating identity
thresholds for taxonomic classification.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import polars as pl

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    go = None


@dataclass(frozen=True)
class IdentityBand:
    """Horizontal band indicating identity threshold region."""

    name: str
    min_identity: float
    max_identity: float
    color: str
    opacity: float = 0.15


# Default classification bands matching the scoring thresholds
DEFAULT_BANDS = (
    IdentityBand(
        name="Known Species",
        min_identity=98.0,
        max_identity=100.0,
        color="#2ecc71",  # Green
    ),
    IdentityBand(
        name="Novel Species",
        min_identity=85.0,
        max_identity=98.0,
        color="#f39c12",  # Orange
    ),
    IdentityBand(
        name="Novel Genus",
        min_identity=75.0,
        max_identity=85.0,
        color="#e74c3c",  # Red
    ),
)


class RecruitmentPlotGenerator:
    """Generator for recruitment plot visualizations.

    Creates interactive Plotly scatter plots showing alignment
    percent identity vs genome position for read recruitment analysis.
    """

    def __init__(
        self,
        data: pl.DataFrame,
        bands: tuple[IdentityBand, ...] = DEFAULT_BANDS,
    ):
        """Initialize with recruitment data.

        Args:
            data: Polars DataFrame with columns:
                - genome_name: Reference genome
                - position: Alignment position
                - percent_identity: Alignment identity
            bands: Identity threshold bands to display.
        """
        if not PLOTLY_AVAILABLE:
            raise ImportError(
                "Plotly is required for visualization. "
                "Install with: pip install plotly"
            )

        self.data = data
        self.bands = bands

    def create_figure(
        self,
        *,
        title: str = "Read Recruitment Plot",
        width: int = 1200,
        height: int = 800,
        show_bands: bool = True,
        genome: str | None = None,
        max_points: int = 100000,
        point_size: int = 3,
        point_opacity: float = 0.5,
    ) -> go.Figure:
        """Create recruitment plot figure.

        Args:
            title: Plot title.
            width: Figure width in pixels.
            height: Figure height in pixels.
            show_bands: Show identity threshold bands.
            genome: Filter to specific genome (None for all).
            max_points: Maximum points to plot (subsample if exceeded).
            point_size: Size of scatter points.
            point_opacity: Opacity of scatter points.

        Returns:
            Plotly Figure object.
        """
        # Filter to genome if specified
        plot_data = self.data
        if genome is not None:
            plot_data = plot_data.filter(pl.col("genome_name") == genome)
            title = f"{title} - {genome}"

        # Subsample if too many points
        if len(plot_data) > max_points:
            plot_data = plot_data.sample(n=max_points, seed=42)

        # Create figure
        fig = go.Figure()

        # Add identity bands as background shapes
        if show_bands:
            x_max = plot_data["position"].max() or 1e6

            for band in self.bands:
                fig.add_shape(
                    type="rect",
                    x0=0,
                    x1=x_max,
                    y0=band.min_identity,
                    y1=band.max_identity,
                    fillcolor=band.color,
                    opacity=band.opacity,
                    line=dict(width=0),
                    layer="below",
                )

                # Add band label
                fig.add_annotation(
                    x=x_max * 0.02,
                    y=(band.min_identity + band.max_identity) / 2,
                    text=band.name,
                    showarrow=False,
                    font=dict(size=10, color=band.color),
                    xanchor="left",
                )

        # Get unique genomes for coloring
        genomes = plot_data["genome_name"].unique().sort().to_list()

        # Color palette for genomes
        colors = [
            "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
            "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
        ]

        # Add scatter trace for each genome
        for i, genome_name in enumerate(genomes):
            genome_data = plot_data.filter(pl.col("genome_name") == genome_name)

            fig.add_trace(
                go.Scattergl(  # Use Scattergl for performance with large datasets
                    x=genome_data["position"].to_list(),
                    y=genome_data["percent_identity"].to_list(),
                    mode="markers",
                    name=genome_name,
                    marker=dict(
                        size=point_size,
                        color=colors[i % len(colors)],
                        opacity=point_opacity,
                    ),
                    hovertemplate=(
                        f"<b>{genome_name}</b><br>"
                        "Position: %{x:,.0f}<br>"
                        "Identity: %{y:.1f}%<br>"
                        "<extra></extra>"
                    ),
                )
            )

        # Layout
        fig.update_layout(
            title=dict(text=title, x=0.5, xanchor="center"),
            xaxis=dict(
                title="Genome Position (bp)",
                showgrid=True,
                gridwidth=1,
                gridcolor="rgba(128,128,128,0.2)",
            ),
            yaxis=dict(
                title="Percent Identity (%)",
                range=[65, 102],
                showgrid=True,
                gridwidth=1,
                gridcolor="rgba(128,128,128,0.2)",
            ),
            width=width,
            height=height,
            template="plotly_white",
            legend=dict(
                yanchor="top",
                y=0.99,
                xanchor="right",
                x=0.99,
                bgcolor="rgba(255,255,255,0.8)",
            ),
            hovermode="closest",
        )

        return fig

    def create_multi_genome_figure(
        self,
        *,
        genomes: list[str] | None = None,
        title: str = "Read Recruitment by Genome",
        width: int = 1400,
        height_per_genome: int = 300,
        show_bands: bool = True,
        max_points_per_genome: int = 50000,
    ) -> go.Figure:
        """Create multi-panel figure with one subplot per genome.

        Args:
            genomes: List of genomes to include (None for top 5 by reads).
            title: Overall figure title.
            width: Figure width in pixels.
            height_per_genome: Height per genome subplot.
            show_bands: Show identity threshold bands.
            max_points_per_genome: Maximum points per genome subplot.

        Returns:
            Plotly Figure with subplots.
        """
        # Determine genomes to plot
        if genomes is None:
            genome_counts = (
                self.data.group_by("genome_name")
                .len()
                .sort("len", descending=True)
                .head(5)
            )
            genomes = genome_counts["genome_name"].to_list()

        n_genomes = len(genomes)

        # Create subplots
        fig = make_subplots(
            rows=n_genomes,
            cols=1,
            shared_xaxes=False,
            subplot_titles=genomes,
            vertical_spacing=0.05,
        )

        colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]

        for i, genome_name in enumerate(genomes, 1):
            genome_data = self.data.filter(pl.col("genome_name") == genome_name)

            if len(genome_data) > max_points_per_genome:
                genome_data = genome_data.sample(n=max_points_per_genome, seed=42)

            # Add identity bands
            if show_bands and len(genome_data) > 0:
                x_max = genome_data["position"].max() or 1e6

                for band in self.bands:
                    fig.add_shape(
                        type="rect",
                        x0=0,
                        x1=x_max,
                        y0=band.min_identity,
                        y1=band.max_identity,
                        fillcolor=band.color,
                        opacity=band.opacity,
                        line=dict(width=0),
                        row=i,
                        col=1,
                        layer="below",
                    )

            # Add scatter
            fig.add_trace(
                go.Scattergl(
                    x=genome_data["position"].to_list(),
                    y=genome_data["percent_identity"].to_list(),
                    mode="markers",
                    name=genome_name,
                    marker=dict(
                        size=2,
                        color=colors[(i - 1) % len(colors)],
                        opacity=0.5,
                    ),
                    showlegend=False,
                    hovertemplate=(
                        f"<b>{genome_name}</b><br>"
                        "Position: %{x:,.0f}<br>"
                        "Identity: %{y:.1f}%<br>"
                        "<extra></extra>"
                    ),
                ),
                row=i,
                col=1,
            )

            # Update y-axis range
            fig.update_yaxes(range=[65, 102], row=i, col=1)

        # Update layout
        total_height = height_per_genome * n_genomes + 100
        fig.update_layout(
            title=dict(text=title, x=0.5, xanchor="center"),
            width=width,
            height=total_height,
            template="plotly_white",
        )

        return fig

    def save(
        self,
        output_path: Path,
        output_format: Literal["html", "png", "json"] = "html",
        **kwargs,
    ) -> None:
        """Save plot to file.

        Args:
            output_path: Output file path.
            output_format: Format - 'html', 'png', or 'json'.
            **kwargs: Additional arguments passed to create_figure().
        """
        fig = self.create_figure(**kwargs)

        if output_format == "html":
            fig.write_html(output_path, include_plotlyjs=True)
        elif output_format == "png":
            fig.write_image(output_path)
        elif output_format == "json":
            fig.write_json(output_path)
        else:
            raise ValueError(f"Unknown format: {output_format}")

    def export_for_anvio(self, output_path: Path) -> None:
        """Export data in format compatible with anvi'o.

        Creates a TSV file that can be imported into anvi'o for
        visualization with anvi-interactive.

        Args:
            output_path: Output TSV file path.
        """
        # Export columns needed for anvi'o visualization
        export_df = self.data.select([
            "genome_name",
            "contig",
            "position",
            "percent_identity",
            "alignment_length",
        ]).rename({
            "genome_name": "split_name",
            "contig": "contig_name",
            "position": "pos",
            "percent_identity": "percent_id",
            "alignment_length": "alignment_len",
        })

        export_df.write_csv(output_path, separator="\t")


def create_recruitment_plot(
    data: pl.DataFrame,
    output_path: Path,
    output_format: Literal["html", "png", "json"] = "html",
    **kwargs,
) -> None:
    """Convenience function to create and save a recruitment plot.

    Args:
        data: Polars DataFrame with recruitment data.
        output_path: Output file path.
        output_format: Format - 'html', 'png', or 'json'.
        **kwargs: Additional arguments passed to create_figure().
    """
    generator = RecruitmentPlotGenerator(data)
    generator.save(output_path, output_format=output_format, **kwargs)
