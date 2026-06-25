"""Shared state and plotting primitives for ReportGenerator section mixins."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

import polars as pl

from metadarkmatter.core.metadata import GenomeMetadata
from metadarkmatter.visualization.plots.base import (
    BAYESIAN_CATEGORY_COLORS,
)
from metadarkmatter.visualization.report.models import (
    ReportConfig,
    TaxonomicSummary,
)

if TYPE_CHECKING:
    import plotly.graph_objects as go



logger = logging.getLogger(__name__)


class _ReportBase:
    """Instance attributes and low-level plot helpers shared by all
    report-section mixins. ReportGenerator supplies the real values via
    its __init__; the mixins only read them."""

    df: pl.DataFrame
    config: ReportConfig
    recruitment_data: pl.DataFrame | None
    ani_matrix: pl.DataFrame | None
    aai_matrix: pl.DataFrame | None
    genome_metadata: GenomeMetadata | None
    summary: TaxonomicSummary
    _plot_data: dict[str, tuple[str, str]]

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

    def _register_plot(self, plot_id: str, fig: go.Figure) -> None:
        """Register a plot for later JS initialization."""
        plot_json = fig.to_json()
        self._plot_data[plot_id] = plot_json

    def _build_histogram(
        self,
        data: list[float],
        plot_id: str,
        title: str,
        x_label: str,
        y_label: str = "Number of Reads",
        nbins: int = 30,
        color: str = "#667eea",
        thresholds: list[tuple[float, str, str, str]] | None = None,
        height: int = 400,
        x_range: list[float] | None = None,
    ) -> str:
        """Build a histogram plot and return PLOT_CONTAINER_TEMPLATE HTML.

        Args:
            data: Values to plot.
            plot_id: Unique DOM id for the plot div.
            title: Plot title.
            x_label: X-axis label.
            y_label: Y-axis label.
            nbins: Number of histogram bins.
            color: Bar fill color.
            thresholds: Optional list of (x_value, line_dash, line_color, annotation_text).
            height: Plot height in pixels.
            x_range: Optional [min, max] for x-axis.

        Returns:
            HTML string (empty if data is empty).
        """
        import plotly.graph_objects as go

        if not data:
            return ""

        fig = go.Figure()
        fig.add_trace(go.Histogram(
            x=data,
            nbinsx=nbins,
            marker_color=color,
            hovertemplate=f"{x_label}: %{{x:.2f}}<br>Count: %{{y}}<extra></extra>",
        ))

        if thresholds:
            for x_val, dash, lcolor, annotation in thresholds:
                fig.add_vline(x=x_val, line_dash=dash, line_color=lcolor,
                              annotation_text=annotation)

        layout_kwargs: dict[str, Any] = {
            "title": title,
            "xaxis_title": x_label,
            "yaxis_title": y_label,
            "template": "plotly_white",
            "height": height,
            "showlegend": False,
        }
        if x_range is not None:
            layout_kwargs["xaxis"] = {"range": x_range}
        fig.update_layout(**layout_kwargs)

        self._register_plot(plot_id, fig)
        return ""

    def _build_category_scatter(
        self,
        x_col: str,
        y_col: str,
        plot_id: str,
        title: str,
        x_label: str,
        y_label: str,
        categories: dict[str, str] | None = None,
        height: int = 500,
        max_points: int | None = None,
        y_range: list[float] | None = None,
    ) -> str:
        """Build a per-category scatter plot and register it.

        Args:
            x_col: DataFrame column for x-axis values.
            y_col: DataFrame column for y-axis values.
            plot_id: Unique DOM id.
            title: Plot title.
            x_label: X-axis label.
            y_label: Y-axis label.
            categories: Dict mapping category name to color. Defaults to BAYESIAN_CATEGORY_COLORS.
            height: Plot height in pixels.
            max_points: If set, subsample the DataFrame.
            y_range: Optional [min, max] for y-axis.

        Returns:
            Empty string (plot is registered internally).
        """
        import plotly.graph_objects as go

        if x_col not in self.df.columns or y_col not in self.df.columns:
            return ""

        if categories is None:
            categories = BAYESIAN_CATEGORY_COLORS

        plot_df = self.df
        effective_max = max_points or self.config.max_scatter_points
        if len(plot_df) > effective_max:
            from metadarkmatter.core.random import get_seed
            plot_df = plot_df.sample(n=effective_max, seed=get_seed())

        fig = go.Figure()
        for call, color in categories.items():
            call_df = plot_df.filter(pl.col("taxonomic_call") == call)
            if len(call_df) > 0:
                fig.add_trace(go.Scattergl(
                    x=call_df[x_col].to_list(),
                    y=call_df[y_col].to_list(),
                    mode="markers",
                    name=call,
                    marker={"color": color, "size": 4, "opacity": 0.5},
                    hovertemplate=(
                        f"<b>{call}</b><br>"
                        f"{x_label}: %{{x:.1f}}<br>"
                        f"{y_label}: %{{y:.2f}}<extra></extra>"
                    ),
                ))

        layout_kwargs: dict[str, Any] = {
            "title": title,
            "xaxis_title": x_label,
            "yaxis_title": y_label,
            "template": "plotly_white",
            "height": height,
            "legend": {"orientation": "h", "y": -0.15},
        }
        if y_range is not None:
            layout_kwargs["yaxis"] = {"range": y_range}
        fig.update_layout(**layout_kwargs)

        self._register_plot(plot_id, fig)
        return ""

    def _build_plotly_cdn_script(self) -> str:
        """Build the Plotly.js script tag based on config.

        Uses the report-mode setting to decide between an inline copy of
        plotly.js (from the installed plotly Python package) and the public
        CDN. The 'include_plotlyjs' value 'embed' is treated as a request to
        keep this tag empty because the JS is already inlined elsewhere.
        """
        from metadarkmatter.visualization.report.assets import get_plotly_script_tag

        if self.config.include_plotlyjs == "embed":
            return ""
        return get_plotly_script_tag(self.config.report_mode)

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

        # Toggle function for ANI/AAI heatmap views
        js_lines.append("""
function toggleHeatmapView(metric, view) {
    var allView = document.getElementById(metric + '-view-all');
    var repsView = document.getElementById(metric + '-view-reps');
    var allBtn = document.getElementById(metric + '-toggle-all');
    var repsBtn = document.getElementById(metric + '-toggle-reps');
    if (!allView || !repsView) return;
    if (view === 'all') {
        allView.style.display = '';
        repsView.style.display = 'none';
        allBtn.classList.add('active');
        repsBtn.classList.remove('active');
    } else {
        allView.style.display = 'none';
        repsView.style.display = '';
        allBtn.classList.remove('active');
        repsBtn.classList.add('active');
        // Resize Plotly chart in the newly visible container
        var plotDiv = repsView.querySelector('.plotly-chart');
        if (plotDiv && typeof Plotly !== 'undefined') {
            Plotly.Plots.resize(plotDiv);
        }
    }
}
""")

        js_lines.append("</script>")

        return "\n".join(js_lines)

    def _build_phylogeny_content(self) -> str:
        """Provided by PhylogenyMixin; declared here for cross-mixin typing."""
        raise NotImplementedError
