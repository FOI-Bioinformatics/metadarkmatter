# tests/integration/test_report_phylogeny.py
"""Integration tests for phylogeny in report generation."""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.report import app

runner = CliRunner()


@pytest.fixture
def sample_data(tmp_path):
    """Create sample data files for testing."""
    # Classifications with novel reads
    classifications = pl.DataFrame({
        "read_id": [f"read_{i}" for i in range(100)],
        "best_match_genome": (
            ["GCF_000001.1"] * 40 +
            ["GCF_000002.1"] * 30 +
            ["GCF_000003.1"] * 30
        ),
        "top_hit_identity": [99.5] * 40 + [88.0] * 30 + [78.0] * 30,
        "novelty_index": [0.5] * 40 + [12.0] * 30 + [22.0] * 30,
        "placement_uncertainty": [0.2] * 40 + [1.5] * 30 + [1.8] * 30,
        "num_ambiguous_hits": [1] * 40 + [2] * 30 + [3] * 30,
        "taxonomic_call": (
            ["Known Species"] * 40 +
            ["Novel Species"] * 30 +
            ["Novel Genus"] * 30
        ),
        "is_novel": [False] * 40 + [True] * 60,
    })

    # ANI matrix
    ani_matrix = pl.DataFrame({
        "genome": ["GCF_000001.1", "GCF_000002.1", "GCF_000003.1"],
        "GCF_000001.1": [100.0, 85.5, 78.2],
        "GCF_000002.1": [85.5, 100.0, 82.1],
        "GCF_000003.1": [78.2, 82.1, 100.0],
    })

    # Save files
    class Data:
        pass

    data = Data()
    data.classifications = tmp_path / "classifications.csv"
    data.ani = tmp_path / "ani_matrix.csv"
    data.output = tmp_path / "report.html"

    classifications.write_csv(data.classifications)
    ani_matrix.write_csv(data.ani)

    return data


class TestReportPhylogenyIntegration:
    """Integration tests for phylogeny tab in reports."""

    def test_report_with_phylogeny_tab(self, sample_data):
        """Report includes phylogeny tab when ANI matrix provided."""
        result = runner.invoke(app, [
            "generate",
            "--classifications", str(sample_data.classifications),
            "--ani", str(sample_data.ani),
            "--output", str(sample_data.output),
        ])

        assert result.exit_code == 0
        html = sample_data.output.read_text()

        # Check phylogeny tab exists
        assert "Phylogeny" in html
        assert "phylogeny-container" in html
        assert "TREE_DATA" in html

    def test_report_phylogeny_includes_novel_clusters(self, sample_data):
        """Phylogeny tab shows novel clusters."""
        result = runner.invoke(app, [
            "generate",
            "--classifications", str(sample_data.classifications),
            "--ani", str(sample_data.ani),
            "--output", str(sample_data.output),
        ])

        assert result.exit_code == 0
        html = sample_data.output.read_text()

        # Check novel clusters are included (clusters are named with [NOVEL] prefix)
        assert "[NOVEL]" in html or "annotations" in html

    def test_report_no_phylogeny_flag(self, sample_data):
        """--no-phylogeny skips phylogeny tab."""
        result = runner.invoke(app, [
            "generate",
            "--classifications", str(sample_data.classifications),
            "--ani", str(sample_data.ani),
            "--no-phylogeny",
            "--output", str(sample_data.output),
        ])

        assert result.exit_code == 0
        html = sample_data.output.read_text()

        # Check phylogeny tab is NOT present
        assert "phylogeny-container" not in html

    def test_report_without_ani(self, sample_data, tmp_path):
        """Report without ANI matrix has no phylogeny tab."""
        output = tmp_path / "report_no_ani.html"

        result = runner.invoke(app, [
            "generate",
            "--classifications", str(sample_data.classifications),
            "--output", str(output),
        ])

        assert result.exit_code == 0
        html = output.read_text()

        # Check phylogeny tab is NOT present (no ANI = no tree)
        assert "phylogeny-container" not in html

    def test_report_with_user_tree(self, sample_data, tmp_path):
        """Report uses user-provided tree file."""
        tree_file = tmp_path / "tree.nwk"
        tree_file.write_text(
            "((GCF_000001.1:0.1,GCF_000002.1:0.2):0.3,GCF_000003.1:0.4);"
        )

        result = runner.invoke(app, [
            "generate",
            "--classifications", str(sample_data.classifications),
            "--ani", str(sample_data.ani),
            "--tree", str(tree_file),
            "--output", str(sample_data.output),
        ])

        assert result.exit_code == 0
        html = sample_data.output.read_text()

        # Check tree was used - the tree source note should indicate user-provided
        assert "phylogeny-container" in html
        assert "user-provided" in html.lower()

    def test_phylogeny_tab_has_d3_script(self, sample_data):
        """Phylogeny tab includes D3.js library."""
        result = runner.invoke(app, [
            "generate",
            "--classifications", str(sample_data.classifications),
            "--ani", str(sample_data.ani),
            "--output", str(sample_data.output),
        ])

        assert result.exit_code == 0
        html = sample_data.output.read_text()

        # Check D3.js is included
        assert "d3.v7.min.js" in html or "d3js.org" in html

    def test_phylogeny_tab_has_legend(self, sample_data):
        """Phylogeny tab includes legend for node types."""
        result = runner.invoke(app, [
            "generate",
            "--classifications", str(sample_data.classifications),
            "--ani", str(sample_data.ani),
            "--output", str(sample_data.output),
        ])

        assert result.exit_code == 0
        html = sample_data.output.read_text()

        # Check legend elements are present
        assert "Reference Genome" in html
        assert "Novel Species" in html
        assert "Novel Genus" in html
