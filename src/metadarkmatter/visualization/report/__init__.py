"""
Report generation module for metadarkmatter.

Provides unified HTML report generation combining multiple visualizations
into a single, self-contained, interactive document.
"""

from metadarkmatter.visualization.report.generator import (
    ReportConfig,
    ReportGenerator,
    TaxonomicSummary,
    generate_report,
)
from metadarkmatter.visualization.report.multi_generator import (
    MultiSampleConfig,
    MultiSampleReportGenerator,
    generate_multi_sample_report,
)
from metadarkmatter.visualization.report.styles import (
    DARK_THEME,
    LIGHT_THEME,
    get_css_styles,
)

__all__ = [
    "DARK_THEME",
    "LIGHT_THEME",
    "MultiSampleConfig",
    "MultiSampleReportGenerator",
    "ReportConfig",
    "ReportGenerator",
    "TaxonomicSummary",
    "generate_multi_sample_report",
    "generate_report",
    "get_css_styles",
]
