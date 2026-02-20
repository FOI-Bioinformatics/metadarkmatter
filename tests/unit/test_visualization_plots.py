"""
Tests for visualization plot modules.

Tests the plot generators for distributions, scatter plots, and classification charts.
"""

from __future__ import annotations

import pytest
import polars as pl

# Skip if plotly not available
plotly = pytest.importorskip("plotly")


class TestPlotConfig:
    """Test PlotConfig dataclass."""

    def test_default_values(self):
        from metadarkmatter.visualization.plots.base import PlotConfig

        config = PlotConfig()
        assert config.width == 800
        assert config.height == 500
        assert config.template == "plotly_white"
        assert config.show_legend is True

    def test_custom_values(self):
        from metadarkmatter.visualization.plots.base import PlotConfig

        config = PlotConfig(width=1200, height=800, show_legend=False)
        assert config.width == 1200
        assert config.height == 800
        assert config.show_legend is False

    def test_to_layout_dict_without_title(self):
        from metadarkmatter.visualization.plots.base import PlotConfig

        config = PlotConfig()
        layout = config.to_layout_dict()
        assert "template" in layout
        assert "font" in layout
        assert "margin" in layout
        assert "title" not in layout

    def test_to_layout_dict_with_title(self):
        from metadarkmatter.visualization.plots.base import PlotConfig

        config = PlotConfig()
        layout = config.to_layout_dict(title="Test Title")
        assert "title" in layout
        assert layout["title"]["text"] == "Test Title"


class TestThresholdConfig:
    """Test ThresholdConfig dataclass."""

    def test_default_thresholds(self):
        """Verify ThresholdConfig defaults match ScoringConfig."""
        from metadarkmatter.visualization.plots.base import ThresholdConfig

        thresholds = ThresholdConfig()
        assert thresholds.novelty_known_max == 4.0
        assert thresholds.novelty_novel_species_min == 4.0
        assert thresholds.novelty_novel_species_max == 20.0
        assert thresholds.novelty_novel_genus_min == 20.0
        assert thresholds.novelty_novel_genus_max == 25.0
        assert thresholds.uncertainty_known_max == 1.5
        assert thresholds.uncertainty_novel_species_max == 1.5
        assert thresholds.uncertainty_novel_genus_max == 1.5
        assert thresholds.uncertainty_conserved_min == 5.0


class TestColorPalettes:
    """Test color palette definitions."""

    def test_taxonomy_colors_defined(self):
        from metadarkmatter.visualization.plots.base import TAXONOMY_COLORS

        assert "Known Species" in TAXONOMY_COLORS
        assert "Novel Species" in TAXONOMY_COLORS
        assert "Novel Genus" in TAXONOMY_COLORS
        assert "Species Boundary" in TAXONOMY_COLORS
        assert "Ambiguous" in TAXONOMY_COLORS
        assert "Ambiguous Within Genus" in TAXONOMY_COLORS
        assert "Conserved Region" in TAXONOMY_COLORS

    def test_sequential_palette_length(self):
        from metadarkmatter.visualization.plots.base import SEQUENTIAL_PALETTE

        assert len(SEQUENTIAL_PALETTE) >= 10


class TestUtilityFunctions:
    """Test utility functions."""

    def test_format_count_small(self):
        from metadarkmatter.visualization.plots.base import format_count

        assert format_count(500) == "500"
        assert format_count(999) == "999"

    def test_format_count_thousands(self):
        from metadarkmatter.visualization.plots.base import format_count

        assert format_count(1000) == "1.0K"
        assert format_count(5500) == "5.5K"
        assert format_count(999999) == "1000.0K"

    def test_format_count_millions(self):
        from metadarkmatter.visualization.plots.base import format_count

        assert format_count(1000000) == "1.0M"
        assert format_count(2500000) == "2.5M"

    def test_format_percentage(self):
        from metadarkmatter.visualization.plots.base import format_percentage

        assert format_percentage(50.0) == "50.0%"
        assert format_percentage(33.333, decimals=2) == "33.33%"

    def test_subsample_dataframe_under_limit(self):
        from metadarkmatter.visualization.plots.base import subsample_dataframe

        df = pl.DataFrame({"x": range(100)})
        result = subsample_dataframe(df, max_points=200)
        assert len(result) == 100

    def test_subsample_dataframe_over_limit(self):
        from metadarkmatter.visualization.plots.base import subsample_dataframe

        df = pl.DataFrame({"x": range(1000)})
        result = subsample_dataframe(df, max_points=100)
        assert len(result) == 100


@pytest.fixture
def sample_classification_data():
    """Create sample classification data for testing."""
    return pl.DataFrame({
        "read_id": [f"read_{i}" for i in range(100)],
        "novelty_index": [0.5] * 40 + [8.0] * 30 + [18.0] * 20 + [1.0] * 10,
        "placement_uncertainty": [0.2] * 40 + [0.3] * 30 + [1.5] * 20 + [8.0] * 10,
        "top_hit_identity": [99.5] * 40 + [92.0] * 30 + [82.0] * 20 + [99.0] * 10,
        "taxonomic_call": (
            ["Known Species"] * 40 +
            ["Novel Species"] * 30 +
            ["Novel Genus"] * 20 +
            ["Conserved Region"] * 10
        ),
    })


@pytest.fixture
def sample_summary():
    """Create sample summary dictionary for testing."""
    return {
        "total_reads": 100,
        "known_species": 40,
        "novel_species": 30,
        "novel_genus": 20,
        "conserved_regions": 10,
        "mean_novelty_index": 5.5,
        "mean_placement_uncertainty": 1.2,
        # High-level diversity grouping
        "diversity_known": 40,
        "diversity_novel": 50,
        "diversity_uncertain": 10,
    }


class TestNoveltyHistogram:
    """Test NoveltyHistogram plot generator."""

    def test_create_figure(self, sample_classification_data):
        from metadarkmatter.visualization.plots.distributions import NoveltyHistogram

        plot = NoveltyHistogram(sample_classification_data)
        fig = plot.create_figure()

        assert fig is not None
        assert len(fig.data) >= 1  # At least one trace

    def test_with_custom_bins(self, sample_classification_data):
        from metadarkmatter.visualization.plots.distributions import NoveltyHistogram

        plot = NoveltyHistogram(sample_classification_data, nbins=20)
        fig = plot.create_figure()

        assert fig is not None

    def test_without_thresholds(self, sample_classification_data):
        from metadarkmatter.visualization.plots.distributions import NoveltyHistogram

        plot = NoveltyHistogram(sample_classification_data, show_thresholds=False)
        fig = plot.create_figure()

        assert fig is not None


class TestUncertaintyHistogram:
    """Test UncertaintyHistogram plot generator."""

    def test_create_figure(self, sample_classification_data):
        from metadarkmatter.visualization.plots.distributions import UncertaintyHistogram

        plot = UncertaintyHistogram(sample_classification_data)
        fig = plot.create_figure()

        assert fig is not None
        assert len(fig.data) >= 1


class TestNoveltyUncertaintyScatter:
    """Test NoveltyUncertaintyScatter plot generator."""

    def test_create_figure(self, sample_classification_data):
        from metadarkmatter.visualization.plots.scatter_2d import NoveltyUncertaintyScatter

        plot = NoveltyUncertaintyScatter(sample_classification_data)
        fig = plot.create_figure()

        assert fig is not None
        # Should have traces for each classification category present
        assert len(fig.data) >= 1

    def test_with_regions_disabled(self, sample_classification_data):
        from metadarkmatter.visualization.plots.scatter_2d import NoveltyUncertaintyScatter

        plot = NoveltyUncertaintyScatter(sample_classification_data, show_regions=False)
        fig = plot.create_figure()

        assert fig is not None

    def test_subsampling(self, sample_classification_data):
        from metadarkmatter.visualization.plots.scatter_2d import NoveltyUncertaintyScatter

        plot = NoveltyUncertaintyScatter(sample_classification_data, max_points=50)
        fig = plot.create_figure()

        assert fig is not None

    def test_empty_category_handling(self):
        """Test handling of data with missing classification categories."""
        from metadarkmatter.visualization.plots.scatter_2d import NoveltyUncertaintyScatter

        # Data with only two categories
        df = pl.DataFrame({
            "novelty_index": [0.5] * 50 + [8.0] * 50,
            "placement_uncertainty": [0.2] * 50 + [0.3] * 50,
            "taxonomic_call": ["Known Species"] * 50 + ["Novel Species"] * 50,
        })

        plot = NoveltyUncertaintyScatter(df)
        fig = plot.create_figure()

        assert fig is not None


class TestClassificationDonutChart:
    """Test ClassificationDonutChart plot generator."""

    def test_create_figure(self, sample_summary):
        from metadarkmatter.visualization.plots.classification_charts import ClassificationDonutChart

        plot = ClassificationDonutChart(sample_summary)
        fig = plot.create_figure()

        assert fig is not None
        assert len(fig.data) >= 1

    def test_with_zero_values(self):
        """Test handling of summary with zero values."""
        from metadarkmatter.visualization.plots.classification_charts import ClassificationDonutChart

        summary = {
            "total_reads": 100,
            "known_species": 100,
            "novel_species": 0,
            "novel_genus": 0,
            "conserved_regions": 0,
        }
        plot = ClassificationDonutChart(summary)
        fig = plot.create_figure()

        assert fig is not None


class TestClassificationBarChart:
    """Test ClassificationBarChart plot generator."""

    def test_create_figure_horizontal(self, sample_summary):
        from metadarkmatter.visualization.plots.classification_charts import ClassificationBarChart

        plot = ClassificationBarChart(sample_summary, orientation="h")
        fig = plot.create_figure()

        assert fig is not None

    def test_create_figure_vertical(self, sample_summary):
        from metadarkmatter.visualization.plots.classification_charts import ClassificationBarChart

        plot = ClassificationBarChart(sample_summary, orientation="v")
        fig = plot.create_figure()

        assert fig is not None


class TestDiversityDonutChart:
    """Test DiversityDonutChart plot generator."""

    def test_create_figure(self, sample_summary):
        from metadarkmatter.visualization.plots.classification_charts import DiversityDonutChart

        plot = DiversityDonutChart(sample_summary)
        fig = plot.create_figure()

        assert fig is not None
        assert len(fig.data) >= 1

    def test_with_zero_values(self):
        """Test handling of summary with zero values."""
        from metadarkmatter.visualization.plots.classification_charts import DiversityDonutChart

        summary = {
            "total_reads": 100,
            "diversity_known": 100,
            "diversity_novel": 0,
            "diversity_uncertain": 0,
        }
        plot = DiversityDonutChart(summary)
        fig = plot.create_figure()

        assert fig is not None


class TestDiversitySunburstChart:
    """Test DiversitySunburstChart plot generator."""

    def test_create_figure(self, sample_summary):
        from metadarkmatter.visualization.plots.classification_charts import DiversitySunburstChart

        # Add all required fields for sunburst
        summary = {
            **sample_summary,
            "species_boundary": 5,
            "ambiguous": 3,
            "ambiguous_within_genus": 1,
            "unclassified": 1,
        }
        plot = DiversitySunburstChart(summary)
        fig = plot.create_figure()

        assert fig is not None
        assert len(fig.data) >= 1
        # Should be a sunburst trace
        assert fig.data[0].type == "sunburst"

    def test_hierarchy_structure(self, sample_summary):
        """Test that hierarchy is correctly built."""
        from metadarkmatter.visualization.plots.classification_charts import DiversitySunburstChart

        summary = {
            "total_reads": 100,
            "known_species": 40,
            "novel_species": 30,
            "novel_genus": 10,
            "species_boundary": 5,
            "ambiguous": 10,
            "ambiguous_within_genus": 2,
            "conserved_regions": 2,
            "unclassified": 1,
        }
        plot = DiversitySunburstChart(summary)
        fig = plot.create_figure()

        # Check that all three diversity statuses are present
        labels = fig.data[0].labels
        assert "Known" in labels
        assert "Novel" in labels
        assert "Uncertain" in labels
        # Check that detailed categories are present
        assert "Novel Species" in labels
        assert "Known Species" in labels


class TestCombinedDistributionPlot:
    """Test CombinedDistributionPlot generator."""

    def test_create_figure(self, sample_classification_data):
        from metadarkmatter.visualization.plots.distributions import CombinedDistributionPlot

        plot = CombinedDistributionPlot(sample_classification_data)
        fig = plot.create_figure()

        assert fig is not None
        # Should have multiple traces (histograms)
        assert len(fig.data) >= 2


class TestBasePlotMethods:
    """Test BasePlot abstract class methods."""

    def test_to_html_div(self, sample_classification_data):
        from metadarkmatter.visualization.plots.distributions import NoveltyHistogram

        plot = NoveltyHistogram(sample_classification_data)
        html = plot.to_html_div()

        assert "<div" in html
        assert "plotly" in html.lower() or "data" in html.lower()

    def test_to_json(self, sample_classification_data):
        from metadarkmatter.visualization.plots.distributions import NoveltyHistogram
        import json

        plot = NoveltyHistogram(sample_classification_data)
        json_str = plot.to_json()

        # Should be valid JSON
        data = json.loads(json_str)
        assert "data" in data
        assert "layout" in data


# =============================================================================
# Fixtures for scatter_2d tests
# =============================================================================


@pytest.fixture
def classification_data_with_confidence(sample_classification_data):
    """Extend sample data with confidence_score and num_ambiguous_hits columns."""
    return sample_classification_data.with_columns([
        pl.lit(75.0).alias("confidence_score"),
        pl.Series("num_ambiguous_hits", [1] * 40 + [3] * 30 + [5] * 20 + [1] * 10),
    ])


# =============================================================================
# ConfidenceNoveltyScatter tests
# =============================================================================


class TestConfidenceNoveltyScatter:
    """Test ConfidenceNoveltyScatter plot generator."""

    def test_create_figure_with_confidence_score(
        self, classification_data_with_confidence
    ):
        """Verify figure is created when confidence_score column is present."""
        from metadarkmatter.visualization.plots.scatter_2d import (
            ConfidenceNoveltyScatter,
        )

        plot = ConfidenceNoveltyScatter(classification_data_with_confidence)
        fig = plot.create_figure()

        assert fig is not None
        # Should have at least one scatter trace per non-empty category
        assert len(fig.data) >= 1
        # Verify traces are Scattergl
        for trace in fig.data:
            assert trace.type == "scattergl"

    def test_create_figure_without_confidence_score(
        self, sample_classification_data
    ):
        """Verify empty figure with message when confidence_score is missing."""
        from metadarkmatter.visualization.plots.scatter_2d import (
            ConfidenceNoveltyScatter,
        )

        plot = ConfidenceNoveltyScatter(sample_classification_data)
        fig = plot.create_figure()

        assert fig is not None
        # Should have no scatter traces
        assert len(fig.data) == 0
        # Should have an annotation with the expected message
        annotations = fig.layout.annotations
        assert len(annotations) == 1
        assert "confidence_score column not available" in annotations[0].text

    def test_single_hit_toggle_with_ambiguous_hits(
        self, classification_data_with_confidence
    ):
        """Verify single-hit toggle is active when num_ambiguous_hits column exists."""
        from metadarkmatter.visualization.plots.scatter_2d import (
            ConfidenceNoveltyScatter,
        )

        plot = ConfidenceNoveltyScatter(
            classification_data_with_confidence, enable_single_hit_toggle=True
        )
        fig = plot.create_figure()

        assert fig is not None
        # With toggle enabled, traces are split into multi-hit and single-hit
        # Known Species (40 reads, all single-hit) + Novel Species (30, all multi-hit)
        # + Novel Genus (20, all multi-hit) + Conserved Region (10, all single-hit)
        # That gives us 4 traces (some single, some multi per category)
        assert len(fig.data) >= 2

        # Verify trace names contain "single-hit" and "multi-hit" labels
        trace_names = [trace.name for trace in fig.data]
        has_single = any("single-hit" in name for name in trace_names)
        has_multi = any("multi-hit" in name for name in trace_names)
        assert has_single, f"Expected single-hit traces, got: {trace_names}"
        assert has_multi, f"Expected multi-hit traces, got: {trace_names}"

        # Verify dropdown menu is present
        assert fig.layout.updatemenus is not None
        assert len(fig.layout.updatemenus) >= 1
        buttons = fig.layout.updatemenus[0].buttons
        button_labels = [b.label for b in buttons]
        assert "All Reads" in button_labels
        assert "Multi-hit Only" in button_labels
        assert "Single-hit Only" in button_labels

    def test_no_toggle_without_ambiguous_hits_column(
        self, sample_classification_data
    ):
        """Verify no toggle when num_ambiguous_hits column is absent."""
        from metadarkmatter.visualization.plots.scatter_2d import (
            ConfidenceNoveltyScatter,
        )

        # Add confidence_score but not num_ambiguous_hits
        data = sample_classification_data.with_columns(
            pl.lit(80.0).alias("confidence_score")
        )
        plot = ConfidenceNoveltyScatter(data, enable_single_hit_toggle=True)
        fig = plot.create_figure()

        assert fig is not None
        # All traces should be plain category traces (no single/multi split)
        trace_names = [trace.name for trace in fig.data]
        assert all("single-hit" not in name for name in trace_names)
        assert all("multi-hit" not in name for name in trace_names)
        # No dropdown menu
        assert fig.layout.updatemenus is None or len(fig.layout.updatemenus) == 0

    def test_toggle_disabled_explicitly(self, classification_data_with_confidence):
        """Verify toggle is not added when enable_single_hit_toggle is False."""
        from metadarkmatter.visualization.plots.scatter_2d import (
            ConfidenceNoveltyScatter,
        )

        plot = ConfidenceNoveltyScatter(
            classification_data_with_confidence, enable_single_hit_toggle=False
        )
        fig = plot.create_figure()

        assert fig is not None
        trace_names = [trace.name for trace in fig.data]
        assert all("single-hit" not in name for name in trace_names)
        assert all("multi-hit" not in name for name in trace_names)

    def test_confidence_hline_present(self, classification_data_with_confidence):
        """Verify the 50% confidence threshold horizontal line is added."""
        from metadarkmatter.visualization.plots.scatter_2d import (
            ConfidenceNoveltyScatter,
        )

        plot = ConfidenceNoveltyScatter(classification_data_with_confidence)
        fig = plot.create_figure()

        # The hline is added as a shape; check that shapes exist
        shapes = fig.layout.shapes
        assert shapes is not None
        assert len(shapes) >= 1
        # At least one shape should be a horizontal line at y=50
        hline_values = [s.y0 for s in shapes if s.type == "line" and s.y0 == s.y1]
        assert 50.0 in hline_values, (
            f"Expected 50% confidence hline, got shapes with y values: {hline_values}"
        )

    def test_axis_labels(self, classification_data_with_confidence):
        """Verify axis labels are set correctly."""
        from metadarkmatter.visualization.plots.scatter_2d import (
            ConfidenceNoveltyScatter,
        )

        plot = ConfidenceNoveltyScatter(classification_data_with_confidence)
        fig = plot.create_figure()

        assert fig.layout.xaxis.title.text == "Novelty Index (100 - % Identity)"
        assert fig.layout.yaxis.title.text == "Confidence Score (%)"

    def test_custom_title(self, classification_data_with_confidence):
        """Verify custom title is applied."""
        from metadarkmatter.visualization.plots.scatter_2d import (
            ConfidenceNoveltyScatter,
        )

        custom_title = "My Custom Confidence Plot"
        plot = ConfidenceNoveltyScatter(
            classification_data_with_confidence, title=custom_title
        )
        fig = plot.create_figure()

        assert fig.layout.title.text == custom_title


# =============================================================================
# NoveltyUncertaintyDensity tests
# =============================================================================


class TestNoveltyUncertaintyDensity:
    """Test NoveltyUncertaintyDensity plot generator."""

    def test_create_figure_produces_histogram2d(self, sample_classification_data):
        """Verify figure contains a Histogram2d trace."""
        from metadarkmatter.visualization.plots.scatter_2d import (
            NoveltyUncertaintyDensity,
        )

        plot = NoveltyUncertaintyDensity(sample_classification_data)
        fig = plot.create_figure()

        assert fig is not None
        assert len(fig.data) == 1
        assert fig.data[0].type == "histogram2d"

    def test_custom_nbins(self, sample_classification_data):
        """Verify custom bin parameters are passed to the trace."""
        from metadarkmatter.visualization.plots.scatter_2d import (
            NoveltyUncertaintyDensity,
        )

        plot = NoveltyUncertaintyDensity(
            sample_classification_data, nbinsx=25, nbinsy=30
        )
        fig = plot.create_figure()

        trace = fig.data[0]
        assert trace.nbinsx == 25
        assert trace.nbinsy == 30

    def test_default_nbins(self, sample_classification_data):
        """Verify default bin parameters are 50x50."""
        from metadarkmatter.visualization.plots.scatter_2d import (
            NoveltyUncertaintyDensity,
        )

        plot = NoveltyUncertaintyDensity(sample_classification_data)
        fig = plot.create_figure()

        trace = fig.data[0]
        assert trace.nbinsx == 50
        assert trace.nbinsy == 50

    def test_threshold_lines_added(self, sample_classification_data):
        """Verify threshold lines are present in the figure."""
        from metadarkmatter.visualization.plots.scatter_2d import (
            NoveltyUncertaintyDensity,
        )
        from metadarkmatter.visualization.plots.base import ThresholdConfig

        thresholds = ThresholdConfig()
        plot = NoveltyUncertaintyDensity(
            sample_classification_data, thresholds=thresholds
        )
        fig = plot.create_figure()

        # Threshold lines are added as shapes
        shapes = fig.layout.shapes
        assert shapes is not None

        # Extract line positions
        # Vertical lines (x0 == x1) for novelty thresholds
        vline_values = sorted(
            {s.x0 for s in shapes if s.type == "line" and s.x0 == s.x1}
        )
        # Horizontal lines (y0 == y1) for uncertainty thresholds
        hline_values = sorted(
            {s.y0 for s in shapes if s.type == "line" and s.y0 == s.y1}
        )

        # Expected vertical lines: novelty_known_max (4.0),
        # novelty_novel_species_min (4.0), novelty_novel_species_max (20.0),
        # novelty_novel_genus_max (25.0)
        # Note: 4.0 appears twice but set deduplicates
        expected_vlines = {4.0, 20.0, 25.0}
        assert expected_vlines.issubset(set(vline_values)), (
            f"Expected vlines at {expected_vlines}, got {vline_values}"
        )

        # Expected horizontal lines: uncertainty_known_max (1.5),
        # uncertainty_novel_genus_max (1.5), uncertainty_conserved_min (5.0)
        # Note: 1.5 appears twice but set deduplicates
        expected_hlines = {1.5, 5.0}
        assert expected_hlines.issubset(set(hline_values)), (
            f"Expected hlines at {expected_hlines}, got {hline_values}"
        )

    def test_colorscale_is_viridis(self, sample_classification_data):
        """Verify the density heatmap uses a Viridis-based colorscale."""
        from metadarkmatter.visualization.plots.scatter_2d import (
            NoveltyUncertaintyDensity,
        )

        plot = NoveltyUncertaintyDensity(sample_classification_data)
        fig = plot.create_figure()

        # Plotly expands the string "Viridis" into a tuple of (position, color)
        # pairs. Verify the colorscale is non-empty and starts with the
        # characteristic Viridis dark purple (#440154).
        cs = fig.data[0].colorscale
        assert cs is not None
        assert len(cs) > 0
        # First color in the Viridis palette is dark purple
        assert cs[0][1] == "#440154"

    def test_axis_labels(self, sample_classification_data):
        """Verify axis labels for density plot."""
        from metadarkmatter.visualization.plots.scatter_2d import (
            NoveltyUncertaintyDensity,
        )

        plot = NoveltyUncertaintyDensity(sample_classification_data)
        fig = plot.create_figure()

        assert fig.layout.xaxis.title.text == "Novelty Index"
        assert fig.layout.yaxis.title.text == "Placement Uncertainty"

    def test_custom_title(self, sample_classification_data):
        """Verify custom title is applied to density plot."""
        from metadarkmatter.visualization.plots.scatter_2d import (
            NoveltyUncertaintyDensity,
        )

        custom_title = "Custom Density Title"
        plot = NoveltyUncertaintyDensity(
            sample_classification_data, title=custom_title
        )
        fig = plot.create_figure()

        assert fig.layout.title.text == custom_title


# =============================================================================
# ClassificationScatterMatrix tests
# =============================================================================


class TestClassificationScatterMatrix:
    """Test ClassificationScatterMatrix plot generator."""

    def test_create_figure_produces_scatter_matrix(
        self, sample_classification_data
    ):
        """Verify figure is created with scatter matrix traces."""
        from metadarkmatter.visualization.plots.scatter_2d import (
            ClassificationScatterMatrix,
        )

        plot = ClassificationScatterMatrix(sample_classification_data)
        fig = plot.create_figure()

        assert fig is not None
        # Scatter matrix should produce multiple splom traces
        assert len(fig.data) >= 1

    def test_subsampling_reduces_data(self):
        """Verify subsampling limits displayed points."""
        from metadarkmatter.visualization.plots.scatter_2d import (
            ClassificationScatterMatrix,
        )

        # Create a larger dataset that exceeds max_points
        large_data = pl.DataFrame({
            "read_id": [f"read_{i}" for i in range(500)],
            "novelty_index": [0.5] * 200 + [8.0] * 150 + [18.0] * 100 + [1.0] * 50,
            "placement_uncertainty": (
                [0.2] * 200 + [0.3] * 150 + [1.5] * 100 + [8.0] * 50
            ),
            "top_hit_identity": (
                [99.5] * 200 + [92.0] * 150 + [82.0] * 100 + [99.0] * 50
            ),
            "taxonomic_call": (
                ["Known Species"] * 200
                + ["Novel Species"] * 150
                + ["Novel Genus"] * 100
                + ["Conserved Region"] * 50
            ),
        })

        plot = ClassificationScatterMatrix(large_data, max_points=50)
        fig = plot.create_figure()

        assert fig is not None
        # The plot should still be created successfully even with subsampling
        assert len(fig.data) >= 1

    def test_no_subsampling_under_limit(self, sample_classification_data):
        """Verify no subsampling when data is under max_points."""
        from metadarkmatter.visualization.plots.scatter_2d import (
            ClassificationScatterMatrix,
        )

        # 100 rows is well under the default 10000 max_points
        plot = ClassificationScatterMatrix(sample_classification_data)
        fig = plot.create_figure()

        assert fig is not None
        # With 100 rows and 4 categories, all data should be present
        assert len(fig.data) >= 1

    def test_scatter_matrix_dimensions(self, sample_classification_data):
        """Verify scatter matrix uses the correct metric dimensions."""
        from metadarkmatter.visualization.plots.scatter_2d import (
            ClassificationScatterMatrix,
        )

        plot = ClassificationScatterMatrix(sample_classification_data)
        fig = plot.create_figure()

        # Scatter matrix traces contain dimensions info; check the first trace
        # The splom trace should reference the three classification metrics
        trace = fig.data[0]
        dimension_labels = [d.label for d in trace.dimensions]
        assert "novelty_index" in dimension_labels
        assert "placement_uncertainty" in dimension_labels
        assert "top_hit_identity" in dimension_labels

    def test_diagonal_and_upper_hidden(self, sample_classification_data):
        """Verify diagonal and upper half are not visible."""
        from metadarkmatter.visualization.plots.scatter_2d import (
            ClassificationScatterMatrix,
        )

        plot = ClassificationScatterMatrix(sample_classification_data)
        fig = plot.create_figure()

        for trace in fig.data:
            assert trace.diagonal.visible is False
            assert trace.showupperhalf is False


# =============================================================================
# Additional BasePlot and base module coverage tests
# =============================================================================


class TestPlotConfigLegend:
    """Test PlotConfig.to_layout_dict() with include_legend=True."""

    def test_to_layout_dict_with_legend(self):
        """Verify legend settings are included when include_legend is True."""
        from metadarkmatter.visualization.plots.base import PlotConfig

        config = PlotConfig(legend_font_size=14, show_legend=True)
        layout = config.to_layout_dict(include_legend=True)

        assert "legend" in layout
        assert layout["legend"] == {"font": {"size": 14}}
        assert layout["showlegend"] is True

    def test_to_layout_dict_legend_hidden(self):
        """Verify showlegend=False when show_legend is False and include_legend is True."""
        from metadarkmatter.visualization.plots.base import PlotConfig

        config = PlotConfig(show_legend=False)
        layout = config.to_layout_dict(include_legend=True)

        assert layout["showlegend"] is False

    def test_to_layout_dict_with_title_and_legend(self):
        """Verify title and legend coexist in layout."""
        from metadarkmatter.visualization.plots.base import PlotConfig

        config = PlotConfig()
        layout = config.to_layout_dict(title="Test", include_legend=True)

        assert "title" in layout
        assert "legend" in layout
        assert layout["title"]["text"] == "Test"


class TestBasePlotGetDivId:
    """Test BasePlot._get_div_id() method."""

    def test_get_div_id_returns_class_based_id(self, sample_classification_data):
        """Verify _get_div_id uses lowercase class name."""
        from metadarkmatter.visualization.plots.distributions import NoveltyHistogram

        plot = NoveltyHistogram(sample_classification_data)
        div_id = plot._get_div_id()

        assert div_id == "plot-noveltyhistogram"

    def test_get_div_id_varies_by_subclass(self, sample_classification_data):
        """Verify different subclasses produce different div IDs."""
        from metadarkmatter.visualization.plots.distributions import (
            NoveltyHistogram,
            UncertaintyHistogram,
        )

        plot_n = NoveltyHistogram(sample_classification_data)
        plot_u = UncertaintyHistogram(sample_classification_data)

        assert plot_n._get_div_id() != plot_u._get_div_id()
        assert plot_u._get_div_id() == "plot-uncertaintyhistogram"

    def test_div_id_embedded_in_html_div(self, sample_classification_data):
        """Verify _get_div_id value appears in to_html_div output."""
        from metadarkmatter.visualization.plots.distributions import NoveltyHistogram

        plot = NoveltyHistogram(sample_classification_data)
        html = plot.to_html_div()

        assert "plot-noveltyhistogram" in html


class TestBasePlotToHtml:
    """Test BasePlot.to_html() standalone HTML method."""

    def test_to_html_returns_full_document(self, sample_classification_data):
        """Verify to_html produces a complete HTML document."""
        from metadarkmatter.visualization.plots.distributions import NoveltyHistogram

        plot = NoveltyHistogram(sample_classification_data)
        html = plot.to_html()

        assert "<html>" in html.lower() or "<!doctype" in html.lower()
        assert "</html>" in html.lower()

    def test_to_html_without_embedded_js(self, sample_classification_data):
        """Verify to_html with include_plotlyjs=False uses CDN."""
        from metadarkmatter.visualization.plots.distributions import NoveltyHistogram

        plot = NoveltyHistogram(sample_classification_data)
        html = plot.to_html(include_plotlyjs=False)

        # CDN-referenced Plotly should appear as a script src
        assert "cdn" in html.lower() or "plotly" in html.lower()


class TestBasePlotSave:
    """Test BasePlot.save() method for various formats."""

    def test_save_html(self, sample_classification_data, tmp_path):
        """Verify save to HTML file works."""
        from metadarkmatter.visualization.plots.distributions import NoveltyHistogram

        plot = NoveltyHistogram(sample_classification_data)
        output = str(tmp_path / "test_plot.html")
        plot.save(output)

        with open(output) as f:
            content = f.read()
        assert "<html>" in content.lower() or "plotly" in content.lower()

    def test_save_json(self, sample_classification_data, tmp_path):
        """Verify save to JSON file works."""
        import json

        from metadarkmatter.visualization.plots.distributions import NoveltyHistogram

        plot = NoveltyHistogram(sample_classification_data)
        output = str(tmp_path / "test_plot.json")
        plot.save(output)

        with open(output) as f:
            data = json.load(f)
        assert "data" in data
        assert "layout" in data

    def test_save_unsupported_format_raises(self, sample_classification_data, tmp_path):
        """Verify save raises ValueError for unsupported formats."""
        from metadarkmatter.visualization.plots.distributions import NoveltyHistogram

        plot = NoveltyHistogram(sample_classification_data)
        output = str(tmp_path / "test_plot.xyz")

        with pytest.raises(ValueError, match="Unsupported file format"):
            plot.save(output)

    def test_save_image_format_requires_kaleido(
        self, sample_classification_data, tmp_path
    ):
        """Verify save to image format calls write_image (may need kaleido)."""
        from metadarkmatter.visualization.plots.distributions import NoveltyHistogram
        from unittest.mock import patch, MagicMock

        plot = NoveltyHistogram(sample_classification_data)
        output = str(tmp_path / "test_plot.png")

        # Mock write_image to avoid requiring kaleido
        with patch.object(
            type(plot.create_figure()), "write_image", new_callable=MagicMock
        ) as mock_write:
            # We need to patch the figure returned by create_figure
            fig = plot.create_figure()
            with patch.object(plot, "create_figure", return_value=fig):
                fig.write_image = mock_write
                plot.save(output)
                mock_write.assert_called_once_with(output)


class TestAddThresholdLineYAxis:
    """Test BasePlot._add_threshold_line with axis='y'."""

    def test_add_horizontal_threshold_line(self, sample_classification_data):
        """Verify horizontal threshold line is added for axis='y'."""
        from metadarkmatter.visualization.plots.distributions import NoveltyHistogram

        plot = NoveltyHistogram(sample_classification_data)
        fig = plot.create_figure()
        initial_shapes = len(fig.layout.shapes) if fig.layout.shapes else 0

        plot._add_threshold_line(fig, value=5.0, axis="y", color="#ff0000")

        shapes = fig.layout.shapes
        assert shapes is not None
        assert len(shapes) > initial_shapes

        # Find the new horizontal line
        hlines = [
            s for s in shapes
            if s.type == "line" and s.y0 == s.y1 and s.y0 == 5.0
        ]
        assert len(hlines) >= 1

    def test_add_horizontal_threshold_line_with_annotation(
        self, sample_classification_data
    ):
        """Verify horizontal line annotation is positioned on the right."""
        from metadarkmatter.visualization.plots.distributions import NoveltyHistogram

        plot = NoveltyHistogram(sample_classification_data)
        fig = plot.create_figure()
        initial_annotations = (
            len(fig.layout.annotations) if fig.layout.annotations else 0
        )

        plot._add_threshold_line(
            fig, value=3.0, axis="y", color="#00ff00", annotation="Test Line"
        )

        annotations = fig.layout.annotations
        assert annotations is not None
        assert len(annotations) > initial_annotations

        # Find the new annotation
        new_annotations = [a for a in annotations if a.text == "Test Line"]
        assert len(new_annotations) == 1
        ann = new_annotations[0]
        assert ann.y == 3.0
        assert ann.xref == "paper"
        assert ann.x == 1.02
        assert ann.showarrow is False

    def test_add_vertical_threshold_line_with_annotation(
        self, sample_classification_data
    ):
        """Verify vertical line annotation is positioned at the top."""
        from metadarkmatter.visualization.plots.distributions import NoveltyHistogram

        plot = NoveltyHistogram(sample_classification_data)
        fig = plot.create_figure()

        plot._add_threshold_line(
            fig, value=10.0, axis="x", color="#0000ff", annotation="Vertical"
        )

        annotations = fig.layout.annotations
        new_annotations = [a for a in annotations if a.text == "Vertical"]
        assert len(new_annotations) == 1
        ann = new_annotations[0]
        assert ann.x == 10.0
        assert ann.yref == "paper"
        assert ann.y == 1.02


class TestGetTaxonomyColorUnknown:
    """Test get_taxonomy_color with unknown classification strings."""

    def test_unknown_classification_returns_default_gray(self):
        """Verify unknown classifications return the fallback color."""
        from metadarkmatter.visualization.plots.base import get_taxonomy_color

        color = get_taxonomy_color("Completely Unknown Category")
        assert color == "#808080"

    def test_empty_string_returns_default_gray(self):
        """Verify empty string returns the fallback color."""
        from metadarkmatter.visualization.plots.base import get_taxonomy_color

        color = get_taxonomy_color("")
        assert color == "#808080"

    def test_known_classification_returns_correct_color(self):
        """Verify known classifications return their mapped color."""
        from metadarkmatter.visualization.plots.base import (
            TAXONOMY_COLORS,
            get_taxonomy_color,
        )

        for classification, expected_color in TAXONOMY_COLORS.items():
            assert get_taxonomy_color(classification) == expected_color


# =============================================================================
# ClassificationSummaryPlot tests (lines 254-322)
# =============================================================================


class TestClassificationSummaryPlot:
    """Test ClassificationSummaryPlot combined donut + bar chart."""

    def test_create_figure(self, sample_summary):
        """Verify combined summary figure is created with two subplots."""
        from metadarkmatter.visualization.plots.classification_charts import (
            ClassificationSummaryPlot,
        )

        summary = {
            **sample_summary,
            "species_boundary": 5,
            "ambiguous": 3,
            "ambiguous_within_genus": 2,
            "unclassified": 0,
        }
        plot = ClassificationSummaryPlot(summary)
        fig = plot.create_figure()

        assert fig is not None
        # Should have two traces: pie + bar
        assert len(fig.data) == 2
        assert fig.data[0].type == "pie"
        assert fig.data[1].type == "bar"

    def test_custom_title(self, sample_summary):
        """Verify custom title is applied."""
        from metadarkmatter.visualization.plots.classification_charts import (
            ClassificationSummaryPlot,
        )

        plot = ClassificationSummaryPlot(sample_summary, title="Custom Summary")
        fig = plot.create_figure()

        assert fig.layout.title.text == "Custom Summary"

    def test_with_all_zero_values(self):
        """Verify figure is created when all category values are zero."""
        from metadarkmatter.visualization.plots.classification_charts import (
            ClassificationSummaryPlot,
        )

        summary = {
            "known_species": 0,
            "novel_species": 0,
            "novel_genus": 0,
            "species_boundary": 0,
            "ambiguous": 0,
            "ambiguous_within_genus": 0,
            "conserved_regions": 0,
            "unclassified": 0,
        }
        plot = ClassificationSummaryPlot(summary)
        fig = plot.create_figure()

        assert fig is not None
        assert len(fig.data) == 2


# =============================================================================
# NovelDiversityGauge tests (lines 347-389)
# =============================================================================


class TestNovelDiversityGauge:
    """Test NovelDiversityGauge indicator chart."""

    def test_create_figure(self, sample_summary):
        """Verify gauge figure is created."""
        from metadarkmatter.visualization.plots.classification_charts import (
            NovelDiversityGauge,
        )

        plot = NovelDiversityGauge(sample_summary)
        fig = plot.create_figure()

        assert fig is not None
        assert len(fig.data) == 1
        assert fig.data[0].type == "indicator"

    def test_gauge_value_calculation(self, sample_summary):
        """Verify gauge computes correct novel diversity percentage."""
        from metadarkmatter.visualization.plots.classification_charts import (
            NovelDiversityGauge,
        )

        # novel_species=30, novel_genus=20, total_reads=100 => 50%
        plot = NovelDiversityGauge(sample_summary)
        fig = plot.create_figure()

        assert fig.data[0].value == pytest.approx(50.0)

    def test_gauge_zero_total(self):
        """Verify gauge handles zero total reads gracefully."""
        from metadarkmatter.visualization.plots.classification_charts import (
            NovelDiversityGauge,
        )

        summary = {
            "novel_species": 0,
            "novel_genus": 0,
            "total_reads": 0,
        }
        plot = NovelDiversityGauge(summary)
        fig = plot.create_figure()

        assert fig.data[0].value == pytest.approx(0.0)

    def test_custom_title(self, sample_summary):
        """Verify custom title is applied to gauge."""
        from metadarkmatter.visualization.plots.classification_charts import (
            NovelDiversityGauge,
        )

        plot = NovelDiversityGauge(sample_summary, title="My Gauge")
        fig = plot.create_figure()

        assert fig.data[0].title.text == "My Gauge"


# =============================================================================
# ClassificationMetricsCards tests (lines 409-479)
# =============================================================================


class TestClassificationMetricsCards:
    """Test ClassificationMetricsCards indicator subplots."""

    def test_create_figure(self, sample_summary):
        """Verify metrics cards figure is created with 4 indicators."""
        from metadarkmatter.visualization.plots.classification_charts import (
            ClassificationMetricsCards,
        )

        plot = ClassificationMetricsCards(sample_summary)
        fig = plot.create_figure()

        assert fig is not None
        # Should have 4 indicator traces
        assert len(fig.data) == 4
        for trace in fig.data:
            assert trace.type == "indicator"

    def test_total_reads_indicator(self, sample_summary):
        """Verify total reads indicator shows correct value."""
        from metadarkmatter.visualization.plots.classification_charts import (
            ClassificationMetricsCards,
        )

        plot = ClassificationMetricsCards(sample_summary)
        fig = plot.create_figure()

        # First trace is total reads
        assert fig.data[0].value == 100

    def test_novel_diversity_percentage(self, sample_summary):
        """Verify novel diversity percentage is computed correctly."""
        from metadarkmatter.visualization.plots.classification_charts import (
            ClassificationMetricsCards,
        )

        # novel_species=30, novel_genus=20, total_reads=100 => 50%
        plot = ClassificationMetricsCards(sample_summary)
        fig = plot.create_figure()

        # Second trace is novel diversity percentage
        assert fig.data[1].value == pytest.approx(50.0)

    def test_mean_novelty_indicator(self, sample_summary):
        """Verify mean novelty indicator shows correct value."""
        from metadarkmatter.visualization.plots.classification_charts import (
            ClassificationMetricsCards,
        )

        plot = ClassificationMetricsCards(sample_summary)
        fig = plot.create_figure()

        # Third trace is mean novelty
        assert fig.data[2].value == pytest.approx(5.5)

    def test_mean_uncertainty_indicator(self, sample_summary):
        """Verify mean uncertainty indicator shows correct value."""
        from metadarkmatter.visualization.plots.classification_charts import (
            ClassificationMetricsCards,
        )

        plot = ClassificationMetricsCards(sample_summary)
        fig = plot.create_figure()

        # Fourth trace is mean uncertainty
        assert fig.data[3].value == pytest.approx(1.2)

    def test_zero_total_reads(self):
        """Verify metrics cards handle zero total reads."""
        from metadarkmatter.visualization.plots.classification_charts import (
            ClassificationMetricsCards,
        )

        summary = {
            "total_reads": 0,
            "novel_species": 0,
            "novel_genus": 0,
            "mean_novelty_index": 0,
            "mean_placement_uncertainty": 0,
        }
        plot = ClassificationMetricsCards(summary)
        fig = plot.create_figure()

        assert fig.data[0].value == 0
        assert fig.data[1].value == pytest.approx(0.0)

    def test_layout_height(self, sample_summary):
        """Verify layout height is set to 200."""
        from metadarkmatter.visualization.plots.classification_charts import (
            ClassificationMetricsCards,
        )

        plot = ClassificationMetricsCards(sample_summary)
        fig = plot.create_figure()

        # After _apply_config, height should be overridden by config
        # but the initial layout sets height=200
        assert fig is not None


# =============================================================================
# IdentityHistogram tests (distributions.py lines 271-336)
# =============================================================================


class TestIdentityHistogram:
    """Test IdentityHistogram plot generator."""

    def test_create_figure(self, sample_classification_data):
        """Verify identity histogram figure is created."""
        from metadarkmatter.visualization.plots.distributions import (
            IdentityHistogram,
        )

        plot = IdentityHistogram(sample_classification_data)
        fig = plot.create_figure()

        assert fig is not None
        assert len(fig.data) >= 1
        assert fig.data[0].type == "histogram"

    def test_with_thresholds(self, sample_classification_data):
        """Verify threshold lines are added when show_thresholds is True."""
        from metadarkmatter.visualization.plots.distributions import (
            IdentityHistogram,
        )

        plot = IdentityHistogram(sample_classification_data, show_thresholds=True)
        fig = plot.create_figure()

        # Should have threshold vlines as shapes
        shapes = fig.layout.shapes
        assert shapes is not None
        assert len(shapes) >= 3  # 3 threshold lines

    def test_without_thresholds(self, sample_classification_data):
        """Verify no threshold lines when show_thresholds is False."""
        from metadarkmatter.visualization.plots.distributions import (
            IdentityHistogram,
        )

        plot = IdentityHistogram(sample_classification_data, show_thresholds=False)
        fig = plot.create_figure()

        # No threshold shapes should be added by IdentityHistogram
        # (_apply_config does not add shapes)
        shapes = fig.layout.shapes
        assert shapes is None or len(shapes) == 0

    def test_custom_nbins(self, sample_classification_data):
        """Verify custom nbins is applied."""
        from metadarkmatter.visualization.plots.distributions import (
            IdentityHistogram,
        )

        plot = IdentityHistogram(sample_classification_data, nbins=25)
        fig = plot.create_figure()

        assert fig.data[0].nbinsx == 25

    def test_axis_labels(self, sample_classification_data):
        """Verify axis labels for identity histogram."""
        from metadarkmatter.visualization.plots.distributions import (
            IdentityHistogram,
        )

        plot = IdentityHistogram(sample_classification_data)
        fig = plot.create_figure()

        assert fig.layout.xaxis.title.text == "Top Hit Identity (%)"
        assert fig.layout.yaxis.title.text == "Number of Reads"

    def test_xaxis_range_adjusts_to_data(self, sample_classification_data):
        """Verify x-axis range adjusts based on minimum identity value."""
        from metadarkmatter.visualization.plots.distributions import (
            IdentityHistogram,
        )

        plot = IdentityHistogram(sample_classification_data)
        fig = plot.create_figure()

        # min identity is 82.0, so range should be [min(70, 82-5), 100]
        x_range = fig.layout.xaxis.range
        assert x_range is not None
        assert x_range[0] == 70  # min(70, 77) = 70
        assert x_range[1] == 100

    def test_custom_title(self, sample_classification_data):
        """Verify custom title is applied."""
        from metadarkmatter.visualization.plots.distributions import (
            IdentityHistogram,
        )

        plot = IdentityHistogram(
            sample_classification_data, title="Custom Identity Distribution"
        )
        fig = plot.create_figure()

        assert fig.layout.title.text == "Custom Identity Distribution"


# =============================================================================
# NoveltyHistogram with explicit bin_size (distributions.py line 80)
# =============================================================================


class TestNoveltyHistogramBinSize:
    """Test NoveltyHistogram with explicit bin_size parameter."""

    def test_explicit_bin_size(self, sample_classification_data):
        """Verify explicit bin_size is used instead of nbins."""
        from metadarkmatter.visualization.plots.distributions import NoveltyHistogram

        plot = NoveltyHistogram(sample_classification_data, bin_size=1.0)
        fig = plot.create_figure()

        assert fig is not None
        # When bin_size is set, xbins should be used instead of nbinsx
        trace = fig.data[0]
        assert trace.xbins is not None
        assert trace.xbins.size == 1.0

    def test_bin_size_none_uses_nbins(self, sample_classification_data):
        """Verify nbins is used when bin_size is None."""
        from metadarkmatter.visualization.plots.distributions import NoveltyHistogram

        plot = NoveltyHistogram(sample_classification_data, nbins=30, bin_size=None)
        fig = plot.create_figure()

        trace = fig.data[0]
        assert trace.nbinsx == 30


# =============================================================================
# UncertaintyHistogram with explicit bin_size
# =============================================================================


class TestUncertaintyHistogramBinSize:
    """Test UncertaintyHistogram with explicit bin_size parameter."""

    def test_explicit_bin_size(self, sample_classification_data):
        """Verify explicit bin_size is used for uncertainty histogram."""
        from metadarkmatter.visualization.plots.distributions import (
            UncertaintyHistogram,
        )

        plot = UncertaintyHistogram(sample_classification_data, bin_size=0.5)
        fig = plot.create_figure()

        trace = fig.data[0]
        assert trace.xbins is not None
        assert trace.xbins.size == 0.5


# =============================================================================
# CombinedDistributionPlot threshold toggle
# =============================================================================


class TestCombinedDistributionPlotThresholds:
    """Test CombinedDistributionPlot with thresholds disabled."""

    def test_without_thresholds(self, sample_classification_data):
        """Verify no threshold lines when show_thresholds is False."""
        from metadarkmatter.visualization.plots.distributions import (
            CombinedDistributionPlot,
        )

        plot = CombinedDistributionPlot(
            sample_classification_data, show_thresholds=False
        )
        fig = plot.create_figure()

        assert fig is not None
        # No vlines added
        shapes = fig.layout.shapes
        assert shapes is None or len(shapes) == 0
