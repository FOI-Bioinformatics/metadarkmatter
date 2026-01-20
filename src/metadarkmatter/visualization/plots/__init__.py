"""
Plot components for metadarkmatter reporting.

Provides modular, reusable plot generators for various visualization types
including distributions, scatter plots, heatmaps, and genome breakdowns.
"""

from metadarkmatter.visualization.plots.base import (
    ANI_COLORSCALE,
    SEQUENTIAL_PALETTE,
    TAXONOMY_COLORS,
    BasePlot,
    PlotConfig,
    ThresholdConfig,
    format_count,
    format_percentage,
    subsample_dataframe,
)
from metadarkmatter.visualization.plots.classification_charts import (
    ClassificationBarChart,
    ClassificationDonutChart,
    ClassificationMetricsCards,
    ClassificationSummaryPlot,
    NovelDiversityGauge,
)
from metadarkmatter.visualization.plots.distributions import (
    CombinedDistributionPlot,
    IdentityHistogram,
    NoveltyHistogram,
    UncertaintyHistogram,
)
from metadarkmatter.visualization.plots.multi_sample import (
    MultiSampleBarChart,
    MultiSampleHeatmap,
    MultiSampleNoveltyComparison,
    MultiSampleScatterMatrix,
    MultiSampleTimeSeries,
)
from metadarkmatter.visualization.plots.scatter_2d import (
    ClassificationScatterMatrix,
    NoveltyUncertaintyDensity,
    NoveltyUncertaintyScatter,
)

__all__ = [
    "ANI_COLORSCALE",
    "SEQUENTIAL_PALETTE",
    "TAXONOMY_COLORS",
    "BasePlot",
    "ClassificationBarChart",
    "ClassificationDonutChart",
    "ClassificationMetricsCards",
    "ClassificationScatterMatrix",
    "ClassificationSummaryPlot",
    "CombinedDistributionPlot",
    "IdentityHistogram",
    "MultiSampleBarChart",
    "MultiSampleHeatmap",
    "MultiSampleNoveltyComparison",
    "MultiSampleScatterMatrix",
    "MultiSampleTimeSeries",
    "NovelDiversityGauge",
    "NoveltyHistogram",
    "NoveltyUncertaintyDensity",
    "NoveltyUncertaintyScatter",
    "PlotConfig",
    "ThresholdConfig",
    "UncertaintyHistogram",
    "format_count",
    "format_percentage",
    "subsample_dataframe",
]
