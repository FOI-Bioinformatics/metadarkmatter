"""Tests for the tree CLI command."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.main import app

runner = CliRunner()


@pytest.fixture()
def ani_csv(tmp_path: Path) -> Path:
    """Create a minimal ANI matrix CSV."""
    ani = pd.DataFrame(
        {
            "A": [100.0, 95.0, 80.0],
            "B": [95.0, 100.0, 82.0],
            "C": [80.0, 82.0, 100.0],
        },
        index=["A", "B", "C"],
    )
    path = tmp_path / "ani_matrix.csv"
    ani.to_csv(path)
    return path


class TestTreeBuildNJ:
    """Test tree build with NJ method."""

    def test_nj_produces_newick(self, ani_csv: Path, tmp_path: Path) -> None:
        """NJ method produces valid Newick output file."""
        output = tmp_path / "tree.nwk"
        result = runner.invoke(app, [
            "tree", "build",
            "--method", "nj",
            "--ani", str(ani_csv),
            "--output", str(output),
        ])
        assert result.exit_code == 0, result.output
        assert output.exists()
        content = output.read_text().strip()
        assert content.endswith(";")
        assert "A" in content

    def test_nj_without_ani_fails(self, tmp_path: Path) -> None:
        """NJ without --ani should fail."""
        output = tmp_path / "tree.nwk"
        result = runner.invoke(app, [
            "tree", "build",
            "--method", "nj",
            "--output", str(output),
        ])
        assert result.exit_code != 0


class TestTreeBuildUPGMA:
    """Test tree build with UPGMA method."""

    def test_upgma_produces_newick(self, ani_csv: Path, tmp_path: Path) -> None:
        """UPGMA method produces valid Newick output file."""
        output = tmp_path / "tree.nwk"
        result = runner.invoke(app, [
            "tree", "build",
            "--method", "upgma",
            "--ani", str(ani_csv),
            "--output", str(output),
        ])
        assert result.exit_code == 0, result.output
        assert output.exists()
        content = output.read_text().strip()
        assert content.endswith(";")


class TestTreeBuildMashtree:
    """Test tree build with Mashtree method."""

    def test_mashtree_without_genomes_fails(self, tmp_path: Path) -> None:
        """Mashtree without --genomes should fail."""
        output = tmp_path / "tree.nwk"
        result = runner.invoke(app, [
            "tree", "build",
            "--method", "mashtree",
            "--output", str(output),
        ])
        assert result.exit_code != 0


class TestInputValidation:
    """Test CLI input validation."""

    def test_nonexistent_ani_fails(self, tmp_path: Path) -> None:
        """Non-existent ANI file should fail."""
        output = tmp_path / "tree.nwk"
        result = runner.invoke(app, [
            "tree", "build",
            "--method", "nj",
            "--ani", str(tmp_path / "missing.csv"),
            "--output", str(output),
        ])
        assert result.exit_code != 0

    def test_too_few_genomes_in_ani(self, tmp_path: Path) -> None:
        """ANI matrix with < 3 genomes should fail gracefully."""
        ani = pd.DataFrame(
            {"A": [100.0, 95.0], "B": [95.0, 100.0]},
            index=["A", "B"],
        )
        ani_path = tmp_path / "small_ani.csv"
        ani.to_csv(ani_path)

        output = tmp_path / "tree.nwk"
        result = runner.invoke(app, [
            "tree", "build",
            "--method", "nj",
            "--ani", str(ani_path),
            "--output", str(output),
        ])
        assert result.exit_code != 0
