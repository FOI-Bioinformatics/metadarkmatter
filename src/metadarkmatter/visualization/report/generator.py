"""
Report generator for unified HTML reports.

Combines multiple plot components into a single, self-contained HTML report
with tabbed navigation, interactive data tables, and professional styling.
"""

from __future__ import annotations

import html
import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl

from metadarkmatter.core.io_utils import read_dataframe
from metadarkmatter.core.metadata import GenomeMetadata
from metadarkmatter.visualization.plots.base import (
    format_count,
)
from metadarkmatter.visualization.report.models import (
    ReportConfig,
    TaxonomicSummary,
)
from metadarkmatter.visualization.report.sections._base import _ReportBase
from metadarkmatter.visualization.report.sections.classification import ClassificationMixin
from metadarkmatter.visualization.report.sections.data import DataMixin
from metadarkmatter.visualization.report.sections.novel import NovelMixin
from metadarkmatter.visualization.report.sections.phylogeny import PhylogenyMixin
from metadarkmatter.visualization.report.sections.reference import ReferenceMixin
from metadarkmatter.visualization.report.sections.summary import SummaryMixin
from metadarkmatter.visualization.report.styles import get_css_styles
from metadarkmatter.visualization.report.templates import (
    COLLAPSIBLE_PANEL_TEMPLATE,
    DATA_TABLE_JS,
    METHODS_SECTION_TEMPLATE,
    REPORT_BASE_TEMPLATE,
    TAB_NAVIGATION_JS,
)

if TYPE_CHECKING:

    pass


logger = logging.getLogger(__name__)

# Version for report footer
try:
    from metadarkmatter import __version__
except ImportError:
    __version__ = "0.1.0"

__all__ = [
    "ReportConfig",
    "ReportGenerator",
    "TaxonomicSummary",
    "generate_report",
]


class ReportGenerator(
    SummaryMixin,
    ClassificationMixin,
    NovelMixin,
    ReferenceMixin,
    PhylogenyMixin,
    DataMixin,
    _ReportBase,
):
    """Generate a unified, self-contained HTML report from classifications.

    Section-building methods live in the mixins under ``report.sections``;
    this class owns construction, the top-level ``generate`` entry point,
    and the HTML/navigation assembly."""

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
        content_sections = []

        # Summary tab (always, active)
        content_sections.append(self._build_summary_section())

        # Classification tab (always)
        content_sections.append(self._build_classification_section())

        # Novel Diversity tab (conditional)
        if self.summary.diversity_novel > 0:
            content_sections.append(self._build_novel_diversity_section())

        # Reference tab (always)
        content_sections.append(self._build_reference_section())

        # Data table tab (always)
        content_sections.append(self._build_data_section())

        content = "\n".join(content_sections)

        # Build dynamic navigation
        navigation = self._build_navigation()

        # Methods footer (outside tab container)
        methods_footer = self._build_methods_footer()

        # Collect all Plotly JS initialization code
        plotly_js = self._build_plotly_js()

        # Combine JavaScript
        js_scripts = TAB_NAVIGATION_JS + "\n" + DATA_TABLE_JS.format(
            page_size=self.config.page_size
        )

        # Build Plotly CDN script tag (empty when embedding JS inline)
        plotly_cdn_script = self._build_plotly_cdn_script()

        # Build final HTML
        return REPORT_BASE_TEMPLATE.format(
            title=html.escape(self.config.title),
            sample_name=html.escape(self.config.sample_name),
            timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            total_reads=format_count(self.summary.total_reads),
            version=__version__,
            css_styles=get_css_styles(self.config.theme),
            navigation=navigation,
            content=content,
            methods_footer=methods_footer,
            plotly_cdn_script=plotly_cdn_script,
            plotly_js=plotly_js,
            js_scripts=js_scripts,
        )

    def _build_navigation(self) -> str:
        """Build dynamic navigation based on available data."""
        tabs = [
            ("summary", "Summary", True),
            ("classification", "Classification", False),
        ]

        if self.summary.diversity_novel > 0:
            tabs.append(("novel-diversity", "Novel Diversity", False))

        tabs.extend([
            ("reference", "Reference", False),
            ("data", "Data", False),
        ])

        nav_items = []
        for tab_id, label, is_active in tabs:
            active_class = " active" if is_active else ""
            safe_tab_id = html.escape(tab_id, quote=True)
            safe_label = html.escape(label)
            nav_items.append(
                f'        <button class="tab-btn{active_class}" '
                f"onclick=\"showTab('{safe_tab_id}')\">{safe_label}</button>"
            )

        return "\n".join(nav_items)

    def _build_methods_footer(self) -> str:
        """Build Methods as a collapsible footer panel below all tabs."""
        return (
            '<div class="methods-footer">'
            + COLLAPSIBLE_PANEL_TEMPLATE.format(
                default_state="",
                panel_id="methods",
                title="Methods",
                content=METHODS_SECTION_TEMPLATE,
            )
            + "</div>"
        )



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
    classifications = read_dataframe(Path(classifications_path))

    # Load optional data
    ani_matrix = None
    if ani_matrix_path:
        ani_matrix = read_dataframe(Path(ani_matrix_path))

    recruitment_data = None
    if recruitment_data_path:
        recruitment_data = read_dataframe(Path(recruitment_data_path))

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
