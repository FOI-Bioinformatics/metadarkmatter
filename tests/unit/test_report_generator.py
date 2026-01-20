"""
Tests for the report generator module.

Tests the ReportGenerator class and related functions.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import polars as pl

# Skip if plotly not available
plotly = pytest.importorskip("plotly")


@pytest.fixture
def sample_classifications():
    """Create sample classification data for testing."""
    return pl.DataFrame({
        "read_id": [f"read_{i}" for i in range(200)],
        "best_match_genome": ["GCF_000001.1"] * 80 + ["GCF_000002.1"] * 70 + ["GCF_000003.1"] * 50,
        "top_hit_identity": [99.5] * 80 + [92.0] * 70 + [82.0] * 50,
        "novelty_index": [0.5] * 80 + [8.0] * 70 + [18.0] * 50,
        "placement_uncertainty": [0.2] * 80 + [0.3] * 70 + [1.5] * 50,
        "num_ambiguous_hits": [1] * 80 + [2] * 70 + [3] * 50,
        "taxonomic_call": ["Known Species"] * 80 + ["Novel Species"] * 70 + ["Novel Genus"] * 50,
        "is_novel": [False] * 80 + [True] * 120,
    })


@pytest.fixture
def sample_ani_matrix():
    """Create sample ANI matrix for testing."""
    return pl.DataFrame({
        "genome": ["GCF_000001.1", "GCF_000002.1", "GCF_000003.1"],
        "GCF_000001.1": [100.0, 85.5, 78.2],
        "GCF_000002.1": [85.5, 100.0, 82.1],
        "GCF_000003.1": [78.2, 82.1, 100.0],
    })


class TestTaxonomicSummary:
    """Test TaxonomicSummary dataclass."""

    def test_novel_percentage_calculation(self):
        from metadarkmatter.visualization.report.generator import TaxonomicSummary

        summary = TaxonomicSummary(
            total_reads=100,
            known_species=50,
            novel_species=30,
            novel_genus=10,
            conserved_regions=10,
        )
        assert summary.novel_percentage == 40.0  # (30 + 10) / 100 * 100

    def test_novel_percentage_zero_reads(self):
        from metadarkmatter.visualization.report.generator import TaxonomicSummary

        summary = TaxonomicSummary(total_reads=0)
        assert summary.novel_percentage == 0.0

    def test_to_dict(self):
        from metadarkmatter.visualization.report.generator import TaxonomicSummary

        summary = TaxonomicSummary(
            total_reads=100,
            known_species=60,
            novel_species=25,
            novel_genus=10,
            conserved_regions=5,
        )
        d = summary.to_dict()

        assert d["total_reads"] == 100
        assert d["known_species"] == 60
        assert d["novel_percentage"] == 35.0


class TestReportConfig:
    """Test ReportConfig dataclass."""

    def test_default_values(self):
        from metadarkmatter.visualization.report.generator import ReportConfig

        config = ReportConfig()
        assert config.sample_name == "Sample"
        assert config.theme == "light"
        assert config.page_size == 100
        assert config.max_scatter_points == 50000

    def test_custom_values(self):
        from metadarkmatter.visualization.report.generator import ReportConfig

        config = ReportConfig(
            sample_name="Test Sample",
            theme="dark",
            max_scatter_points=10000,
        )
        assert config.sample_name == "Test Sample"
        assert config.theme == "dark"
        assert config.max_scatter_points == 10000


class TestReportGenerator:
    """Test ReportGenerator class."""

    def test_compute_summary(self, sample_classifications):
        from metadarkmatter.visualization.report.generator import ReportGenerator

        generator = ReportGenerator(sample_classifications)

        assert generator.summary.total_reads == 200
        assert generator.summary.known_species == 80
        assert generator.summary.novel_species == 70
        assert generator.summary.novel_genus == 50

    def test_generate_report(self, sample_classifications, tmp_path):
        from metadarkmatter.visualization.report.generator import ReportGenerator

        output_path = tmp_path / "test_report.html"
        generator = ReportGenerator(sample_classifications)
        generator.generate(output_path)

        assert output_path.exists()
        content = output_path.read_text()

        # Check key elements are present
        assert "<!DOCTYPE html>" in content
        assert "Metadarkmatter" in content
        assert "Overview" in content
        assert "Distributions" in content

    def test_generate_with_ani_matrix(self, sample_classifications, sample_ani_matrix, tmp_path):
        from metadarkmatter.visualization.report.generator import ReportGenerator

        output_path = tmp_path / "test_report_ani.html"
        generator = ReportGenerator(
            sample_classifications,
            ani_matrix=sample_ani_matrix,
        )
        generator.generate(output_path)

        assert output_path.exists()
        content = output_path.read_text()

        # ANI section should not show "not provided" message
        assert "ANI matrix not provided" not in content

    def test_generate_dark_theme(self, sample_classifications, tmp_path):
        from metadarkmatter.visualization.report.generator import (
            ReportGenerator,
            ReportConfig,
        )

        output_path = tmp_path / "test_report_dark.html"
        config = ReportConfig(theme="dark")
        generator = ReportGenerator(sample_classifications, config=config)
        generator.generate(output_path)

        assert output_path.exists()
        content = output_path.read_text()

        # Check dark theme variables are present
        assert "--bg-primary: #1a1a2e" in content or "--bg-primary" in content

    def test_build_metric_cards(self, sample_classifications):
        from metadarkmatter.visualization.report.generator import ReportGenerator

        generator = ReportGenerator(sample_classifications)
        cards_html = generator._build_metric_cards()

        # Check metric cards are generated
        assert "Total Reads" in cards_html
        assert "Known Species" in cards_html
        assert "Novel Species" in cards_html

    def test_build_data_section(self, sample_classifications):
        from metadarkmatter.visualization.report.generator import ReportGenerator

        generator = ReportGenerator(sample_classifications)
        section = generator._build_data_section()

        # Check data table elements
        assert "dataTable" in section
        assert "tableSearch" in section


class TestGenerateReportFunction:
    """Test convenience generate_report function."""

    def test_generate_from_file(self, sample_classifications, tmp_path):
        from metadarkmatter.visualization.report.generator import generate_report

        # Write sample data to file
        input_path = tmp_path / "classifications.csv"
        sample_classifications.write_csv(input_path)

        output_path = tmp_path / "report.html"

        generate_report(
            classifications_path=input_path,
            output_path=output_path,
            sample_name="Test Sample",
        )

        assert output_path.exists()


class TestCSSStyles:
    """Test CSS style generation."""

    def test_light_theme(self):
        from metadarkmatter.visualization.report.styles import get_css_styles

        css = get_css_styles("light")
        assert "--bg-primary: #ffffff" in css

    def test_dark_theme(self):
        from metadarkmatter.visualization.report.styles import get_css_styles

        css = get_css_styles("dark")
        assert "--bg-primary: #1a1a2e" in css


class TestHTMLTemplates:
    """Test HTML template strings."""

    def test_report_base_template_placeholders(self):
        from metadarkmatter.visualization.report.templates import REPORT_BASE_TEMPLATE

        # Check required placeholders exist
        assert "{title}" in REPORT_BASE_TEMPLATE
        assert "{sample_name}" in REPORT_BASE_TEMPLATE
        assert "{content}" in REPORT_BASE_TEMPLATE
        assert "{css_styles}" in REPORT_BASE_TEMPLATE

    def test_data_table_js_format(self):
        from metadarkmatter.visualization.report.templates import DATA_TABLE_JS

        # Should be able to format with page_size
        formatted = DATA_TABLE_JS.format(page_size=100)
        assert "const pageSize = 100" in formatted

    def test_get_cell_class(self):
        from metadarkmatter.visualization.report.templates import get_cell_class

        assert get_cell_class("Known Species") == "cell-known"
        assert get_cell_class("Novel Species") == "cell-novel-species"
        assert get_cell_class("Unknown") == ""
