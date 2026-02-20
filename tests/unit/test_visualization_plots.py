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
