"""
Tests for the recruitment_plots visualization module.

Covers IdentityBand dataclass, DEFAULT_BANDS constants,
RecruitmentPlotGenerator (single and multi-genome figures, save, anvi'o export),
and the create_recruitment_plot convenience function.
"""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl
import pytest

# Skip entire module if plotly is not installed
plotly = pytest.importorskip("plotly")

from metadarkmatter.visualization.recruitment_plots import (
    DEFAULT_BANDS,
    IdentityBand,
    RecruitmentPlotGenerator,
    create_recruitment_plot,
)


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def recruitment_data() -> pl.DataFrame:
    """Synthetic recruitment data with three genomes spanning different identity ranges."""
    return pl.DataFrame({
        "position": list(range(0, 1000, 10)) * 3,
        "percent_identity": [98.5] * 100 + [92.0] * 100 + [78.0] * 100,
        "genome_name": ["GCF_001"] * 100 + ["GCF_002"] * 100 + ["GCF_003"] * 100,
    })


@pytest.fixture
def recruitment_data_with_anvio_columns() -> pl.DataFrame:
    """Recruitment data including the extra columns required for anvi'o export."""
    n = 50
    return pl.DataFrame({
        "position": list(range(0, n * 10, 10)) * 2,
        "percent_identity": [97.0] * n + [85.0] * n,
        "genome_name": ["GCF_A"] * n + ["GCF_B"] * n,
        "contig": [f"contig_{i % 5}" for i in range(n)] * 2,
        "alignment_length": [150] * (n * 2),
    })


@pytest.fixture
def single_genome_data() -> pl.DataFrame:
    """Recruitment data containing a single genome."""
    return pl.DataFrame({
        "position": list(range(0, 500, 5)),
        "percent_identity": [99.0] * 100,
        "genome_name": ["GCF_SOLO"] * 100,
    })


@pytest.fixture
def temp_dir(tmp_path: Path) -> Path:
    """Provide a temporary directory via pytest's built-in tmp_path."""
    return tmp_path


# =============================================================================
# IdentityBand dataclass tests
# =============================================================================


class TestIdentityBand:
    """Test the IdentityBand frozen dataclass."""

    def test_creation_with_all_fields(self):
        band = IdentityBand(
            name="Test Band",
            min_identity=90.0,
            max_identity=100.0,
            color="#abcdef",
            opacity=0.25,
        )
        assert band.name == "Test Band"
        assert band.min_identity == 90.0
        assert band.max_identity == 100.0
        assert band.color == "#abcdef"
        assert band.opacity == 0.25

    def test_default_opacity(self):
        band = IdentityBand(
            name="Default Opacity",
            min_identity=80.0,
            max_identity=90.0,
            color="#000000",
        )
        assert band.opacity == 0.15

    def test_frozen_immutability(self):
        band = IdentityBand(
            name="Frozen",
            min_identity=70.0,
            max_identity=80.0,
            color="#111111",
        )
        with pytest.raises(AttributeError):
            band.name = "Modified"  # type: ignore[misc]


# =============================================================================
# DEFAULT_BANDS tests
# =============================================================================


class TestDefaultBands:
    """Verify the module-level DEFAULT_BANDS constant."""

    def test_contains_three_bands(self):
        assert len(DEFAULT_BANDS) == 3

    def test_band_names(self):
        names = [b.name for b in DEFAULT_BANDS]
        assert "Known Species" in names
        assert "Novel Species" in names
        assert "Novel Genus" in names

    def test_known_species_band_range(self):
        known = [b for b in DEFAULT_BANDS if b.name == "Known Species"][0]
        assert known.min_identity == 98.0
        assert known.max_identity == 100.0

    def test_novel_species_band_range(self):
        novel_sp = [b for b in DEFAULT_BANDS if b.name == "Novel Species"][0]
        assert novel_sp.min_identity == 85.0
        assert novel_sp.max_identity == 98.0

    def test_novel_genus_band_range(self):
        novel_gen = [b for b in DEFAULT_BANDS if b.name == "Novel Genus"][0]
        assert novel_gen.min_identity == 75.0
        assert novel_gen.max_identity == 85.0

    def test_bands_are_contiguous(self):
        """The top of Novel Genus should meet the bottom of Novel Species, etc."""
        genus = [b for b in DEFAULT_BANDS if b.name == "Novel Genus"][0]
        species = [b for b in DEFAULT_BANDS if b.name == "Novel Species"][0]
        known = [b for b in DEFAULT_BANDS if b.name == "Known Species"][0]

        assert genus.max_identity == species.min_identity
        assert species.max_identity == known.min_identity


# =============================================================================
# RecruitmentPlotGenerator tests
# =============================================================================


class TestRecruitmentPlotGeneratorInit:
    """Test constructor behaviour of RecruitmentPlotGenerator."""

    def test_stores_data_and_default_bands(self, recruitment_data):
        gen = RecruitmentPlotGenerator(recruitment_data)
        assert gen.data is recruitment_data
        assert gen.bands is DEFAULT_BANDS

    def test_custom_bands(self, recruitment_data):
        custom = (
            IdentityBand("Only Band", 90.0, 100.0, "#ffffff"),
        )
        gen = RecruitmentPlotGenerator(recruitment_data, bands=custom)
        assert len(gen.bands) == 1
        assert gen.bands[0].name == "Only Band"


class TestCreateFigure:
    """Test create_figure() for single-genome and multi-genome data."""

    def test_returns_plotly_figure(self, recruitment_data):
        gen = RecruitmentPlotGenerator(recruitment_data)
        fig = gen.create_figure()
        assert isinstance(fig, plotly.graph_objects.Figure)

    def test_figure_has_traces_per_genome(self, recruitment_data):
        gen = RecruitmentPlotGenerator(recruitment_data)
        fig = gen.create_figure()
        # One ScatterGL trace per unique genome
        unique_genomes = recruitment_data["genome_name"].unique().len()
        assert len(fig.data) == unique_genomes

    def test_figure_with_bands_has_shapes(self, recruitment_data):
        gen = RecruitmentPlotGenerator(recruitment_data)
        fig = gen.create_figure(show_bands=True)
        # 3 default bands -> 3 shapes
        assert len(fig.layout.shapes) == len(DEFAULT_BANDS)

    def test_figure_without_bands_has_no_shapes(self, recruitment_data):
        gen = RecruitmentPlotGenerator(recruitment_data)
        fig = gen.create_figure(show_bands=False)
        assert len(fig.layout.shapes) == 0

    def test_figure_filter_single_genome(self, recruitment_data):
        gen = RecruitmentPlotGenerator(recruitment_data)
        fig = gen.create_figure(genome="GCF_001")
        # Only one trace when filtered to a single genome
        assert len(fig.data) == 1
        assert fig.data[0].name == "GCF_001"

    def test_title_appended_when_genome_specified(self, recruitment_data):
        gen = RecruitmentPlotGenerator(recruitment_data)
        fig = gen.create_figure(title="My Plot", genome="GCF_002")
        assert "GCF_002" in fig.layout.title.text

    def test_custom_dimensions(self, recruitment_data):
        gen = RecruitmentPlotGenerator(recruitment_data)
        fig = gen.create_figure(width=800, height=600)
        assert fig.layout.width == 800
        assert fig.layout.height == 600

    def test_subsampling(self, recruitment_data):
        """When max_points < row count, the traces should contain fewer points."""
        gen = RecruitmentPlotGenerator(recruitment_data)
        fig = gen.create_figure(max_points=50)
        total_points = sum(len(trace.x) for trace in fig.data)
        assert total_points <= 50

    def test_single_genome_data(self, single_genome_data):
        gen = RecruitmentPlotGenerator(single_genome_data)
        fig = gen.create_figure()
        assert len(fig.data) == 1
        assert fig.data[0].name == "GCF_SOLO"

    def test_point_size_and_opacity_propagated(self, single_genome_data):
        gen = RecruitmentPlotGenerator(single_genome_data)
        fig = gen.create_figure(point_size=7, point_opacity=0.9)
        marker = fig.data[0].marker
        assert marker.size == 7
        assert marker.opacity == 0.9


class TestCreateMultiGenomeFigure:
    """Test create_multi_genome_figure() subplot generation."""

    def test_returns_plotly_figure(self, recruitment_data):
        gen = RecruitmentPlotGenerator(recruitment_data)
        fig = gen.create_multi_genome_figure()
        assert isinstance(fig, plotly.graph_objects.Figure)

    def test_default_top5_genomes(self, recruitment_data):
        gen = RecruitmentPlotGenerator(recruitment_data)
        fig = gen.create_multi_genome_figure()
        # Data has 3 genomes, all should appear (fewer than the top-5 cap)
        assert len(fig.data) == 3

    def test_explicit_genome_list(self, recruitment_data):
        gen = RecruitmentPlotGenerator(recruitment_data)
        fig = gen.create_multi_genome_figure(genomes=["GCF_001", "GCF_003"])
        assert len(fig.data) == 2

    def test_single_genome_subplot(self, recruitment_data):
        gen = RecruitmentPlotGenerator(recruitment_data)
        fig = gen.create_multi_genome_figure(genomes=["GCF_002"])
        assert len(fig.data) == 1

    def test_figure_height_scales_with_genomes(self, recruitment_data):
        gen = RecruitmentPlotGenerator(recruitment_data)
        fig2 = gen.create_multi_genome_figure(
            genomes=["GCF_001", "GCF_002"], height_per_genome=300
        )
        fig3 = gen.create_multi_genome_figure(
            genomes=["GCF_001", "GCF_002", "GCF_003"], height_per_genome=300
        )
        assert fig3.layout.height > fig2.layout.height

    def test_bands_shown_by_default(self, recruitment_data):
        gen = RecruitmentPlotGenerator(recruitment_data)
        fig = gen.create_multi_genome_figure()
        # Each genome subplot gets band shapes
        assert len(fig.layout.shapes) == 3 * len(DEFAULT_BANDS)

    def test_bands_hidden(self, recruitment_data):
        gen = RecruitmentPlotGenerator(recruitment_data)
        fig = gen.create_multi_genome_figure(show_bands=False)
        assert len(fig.layout.shapes) == 0

    def test_custom_title(self, recruitment_data):
        gen = RecruitmentPlotGenerator(recruitment_data)
        fig = gen.create_multi_genome_figure(title="Custom Multi Title")
        assert fig.layout.title.text == "Custom Multi Title"

    def test_subsampling_per_genome(self, recruitment_data):
        """Each genome subplot should be subsampled independently."""
        gen = RecruitmentPlotGenerator(recruitment_data)
        fig = gen.create_multi_genome_figure(max_points_per_genome=20)
        for trace in fig.data:
            assert len(trace.x) <= 20


# =============================================================================
# Save tests
# =============================================================================


class TestSave:
    """Test save() for supported output formats."""

    def test_save_html(self, recruitment_data, temp_dir):
        gen = RecruitmentPlotGenerator(recruitment_data)
        out = temp_dir / "plot.html"
        gen.save(out, output_format="html")
        assert out.exists()
        content = out.read_text()
        assert "<html" in content.lower() or "<div" in content.lower()

    def test_save_json(self, recruitment_data, temp_dir):
        gen = RecruitmentPlotGenerator(recruitment_data)
        out = temp_dir / "plot.json"
        gen.save(out, output_format="json")
        assert out.exists()
        data = json.loads(out.read_text())
        assert "data" in data
        assert "layout" in data

    def test_save_unknown_format_raises(self, recruitment_data, temp_dir):
        gen = RecruitmentPlotGenerator(recruitment_data)
        out = temp_dir / "plot.xyz"
        with pytest.raises(ValueError, match="Unknown format"):
            gen.save(out, output_format="xyz")  # type: ignore[arg-type]

    def test_save_passes_kwargs_to_create_figure(self, recruitment_data, temp_dir):
        """Verify that extra kwargs (e.g. title) are forwarded to create_figure."""
        gen = RecruitmentPlotGenerator(recruitment_data)
        out = temp_dir / "titled.html"
        gen.save(out, output_format="html", title="Custom Title")
        assert out.exists()


# =============================================================================
# export_for_anvio tests
# =============================================================================


class TestExportForAnvio:
    """Test the anvi'o-compatible TSV export."""

    def test_creates_tsv_file(self, recruitment_data_with_anvio_columns, temp_dir):
        gen = RecruitmentPlotGenerator(recruitment_data_with_anvio_columns)
        out = temp_dir / "anvio.tsv"
        gen.export_for_anvio(out)
        assert out.exists()

    def test_tsv_column_names(self, recruitment_data_with_anvio_columns, temp_dir):
        gen = RecruitmentPlotGenerator(recruitment_data_with_anvio_columns)
        out = temp_dir / "anvio.tsv"
        gen.export_for_anvio(out)

        result = pl.read_csv(out, separator="\t")
        expected_cols = {"split_name", "contig_name", "pos", "percent_id", "alignment_len"}
        assert set(result.columns) == expected_cols

    def test_tsv_row_count_matches(self, recruitment_data_with_anvio_columns, temp_dir):
        gen = RecruitmentPlotGenerator(recruitment_data_with_anvio_columns)
        out = temp_dir / "anvio.tsv"
        gen.export_for_anvio(out)

        result = pl.read_csv(out, separator="\t")
        assert len(result) == len(recruitment_data_with_anvio_columns)

    def test_tsv_values_match_source(self, recruitment_data_with_anvio_columns, temp_dir):
        gen = RecruitmentPlotGenerator(recruitment_data_with_anvio_columns)
        out = temp_dir / "anvio.tsv"
        gen.export_for_anvio(out)

        result = pl.read_csv(out, separator="\t")
        # Verify data integrity for a few key columns
        assert result["split_name"].to_list() == recruitment_data_with_anvio_columns["genome_name"].to_list()
        assert result["pos"].to_list() == recruitment_data_with_anvio_columns["position"].to_list()


# =============================================================================
# Convenience function test
# =============================================================================


class TestCreateRecruitmentPlotFunction:
    """Test the module-level create_recruitment_plot convenience function."""

    def test_creates_html_file(self, recruitment_data, temp_dir):
        out = temp_dir / "convenience.html"
        create_recruitment_plot(recruitment_data, out, output_format="html")
        assert out.exists()
        assert out.stat().st_size > 0

    def test_creates_json_file(self, recruitment_data, temp_dir):
        out = temp_dir / "convenience.json"
        create_recruitment_plot(recruitment_data, out, output_format="json")
        assert out.exists()
        data = json.loads(out.read_text())
        assert "data" in data

    def test_forwards_kwargs(self, recruitment_data, temp_dir):
        out = temp_dir / "kwargs.html"
        create_recruitment_plot(
            recruitment_data, out, output_format="html", title="Via Convenience"
        )
        assert out.exists()


# =============================================================================
# Edge-case tests
# =============================================================================


class TestEdgeCases:
    """Edge cases and boundary conditions."""

    def test_empty_dataframe(self):
        """An empty DataFrame should still produce a figure without errors."""
        df = pl.DataFrame({
            "position": pl.Series([], dtype=pl.Int64),
            "percent_identity": pl.Series([], dtype=pl.Float64),
            "genome_name": pl.Series([], dtype=pl.Utf8),
        })
        gen = RecruitmentPlotGenerator(df)
        fig = gen.create_figure()
        assert isinstance(fig, plotly.graph_objects.Figure)
        assert len(fig.data) == 0

    def test_single_row(self):
        """A single-row DataFrame should work."""
        df = pl.DataFrame({
            "position": [500],
            "percent_identity": [99.9],
            "genome_name": ["GCF_SINGLE"],
        })
        gen = RecruitmentPlotGenerator(df)
        fig = gen.create_figure()
        assert len(fig.data) == 1
        assert len(fig.data[0].x) == 1

    def test_genome_filter_no_match(self, recruitment_data):
        """Filtering to a nonexistent genome should produce an empty figure."""
        gen = RecruitmentPlotGenerator(recruitment_data)
        fig = gen.create_figure(genome="NONEXISTENT")
        assert len(fig.data) == 0

    def test_many_genomes_color_cycling(self):
        """When there are more genomes than palette entries, colors should cycle."""
        n_genomes = 15
        rows_per = 10
        df = pl.DataFrame({
            "position": list(range(rows_per)) * n_genomes,
            "percent_identity": [95.0] * (rows_per * n_genomes),
            "genome_name": [f"GCF_{i:03d}" for i in range(n_genomes) for _ in range(rows_per)],
        })
        gen = RecruitmentPlotGenerator(df)
        fig = gen.create_figure()
        assert len(fig.data) == n_genomes

    def test_custom_bands_in_figure(self, single_genome_data):
        """Custom bands should appear as shapes in the figure."""
        custom_bands = (
            IdentityBand("High", 95.0, 100.0, "#00ff00"),
            IdentityBand("Low", 70.0, 95.0, "#ff0000"),
        )
        gen = RecruitmentPlotGenerator(single_genome_data, bands=custom_bands)
        fig = gen.create_figure(show_bands=True)
        assert len(fig.layout.shapes) == 2

    def test_multi_genome_empty_genome_list_raises(self, recruitment_data):
        """Passing an empty genome list raises because Plotly requires rows >= 1."""
        gen = RecruitmentPlotGenerator(recruitment_data)
        with pytest.raises(ValueError):
            gen.create_multi_genome_figure(genomes=[])
