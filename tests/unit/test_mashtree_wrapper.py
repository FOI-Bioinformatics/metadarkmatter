"""Tests for Mashtree external tool wrapper."""

from __future__ import annotations

from pathlib import Path

import pytest

from metadarkmatter.external.mashtree import Mashtree


class TestMashtreeCommand:
    """Test Mashtree command building."""

    @pytest.fixture(autouse=True)
    def mock_executable(self) -> None:
        """Mock the executable resolver so tests don't need mashtree installed."""
        Mashtree.set_executable_resolver(lambda name: f"/usr/bin/{name}")
        Mashtree.clear_cache()
        yield
        Mashtree.reset_executable_resolver()

    def test_build_command_basic(self, tmp_path: Path) -> None:
        """Basic command builds correctly."""
        g1 = tmp_path / "genome1.fna"
        g2 = tmp_path / "genome2.fna"
        g1.touch()
        g2.touch()

        mt = Mashtree()
        cmd = mt.build_command(genomes=[g1, g2], threads=4)

        assert cmd[0] == "/usr/bin/mashtree"
        assert "--numcpus" in cmd
        assert "4" in cmd
        assert str(g1) in cmd
        assert str(g2) in cmd

    def test_build_command_custom_threads(self, tmp_path: Path) -> None:
        """Threads parameter is passed correctly."""
        g1 = tmp_path / "genome1.fna"
        g1.touch()

        mt = Mashtree()
        cmd = mt.build_command(genomes=[g1], threads=16)

        idx = cmd.index("--numcpus")
        assert cmd[idx + 1] == "16"

    def test_tool_name(self) -> None:
        """Tool name is 'mashtree'."""
        assert Mashtree.TOOL_NAME == "mashtree"

    def test_install_hint(self) -> None:
        """Install hint includes conda."""
        assert "conda" in Mashtree.INSTALL_HINT

    def test_check_available_with_mock(self) -> None:
        """check_available returns True when mock resolver finds it."""
        assert Mashtree.check_available() is True

    def test_check_available_not_installed(self) -> None:
        """check_available returns False when not installed."""
        Mashtree.set_executable_resolver(lambda name: None)
        Mashtree.clear_cache()
        assert Mashtree.check_available() is False
