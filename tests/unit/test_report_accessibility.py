"""Accessibility and palette tests for the HTML report."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("plotly")

from metadarkmatter.visualization.plots.base import (
    BAYESIAN_CATEGORY_COLORS,
    DIVERSITY_COLORS,
    OKABE_ITO,
    TAXONOMY_COLORS,
)
from tests.unit.report_helpers import render


class TestPalette:
    def test_taxonomy_colors_are_okabe_ito_and_distinct(self):
        values = list(TAXONOMY_COLORS.values())
        assert len(values) == len(set(values))
        assert set(values) <= set(OKABE_ITO.values())
        assert len(TAXONOMY_COLORS) >= 9

    def test_novel_categories_use_warm_saturated_hues(self):
        assert TAXONOMY_COLORS["Novel Species"] == OKABE_ITO["orange"]
        assert TAXONOMY_COLORS["Novel Genus"] == OKABE_ITO["vermillion"]

    def test_derived_palettes_stay_consistent(self):
        for name, color in BAYESIAN_CATEGORY_COLORS.items():
            assert TAXONOMY_COLORS[name] == color
        assert DIVERSITY_COLORS["Novel"] == TAXONOMY_COLORS["Novel Genus"]

    def test_text_variants_meet_contrast_on_white(self):
        from metadarkmatter.visualization.report.styles import PALETTE_LIGHT

        def lum(hex_color: str) -> float:
            rgb = [int(hex_color[i : i + 2], 16) / 255 for i in (1, 3, 5)]
            lin = [c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4 for c in rgb]
            return 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]

        found = re.findall(r"--(cat-[a-z-]+-text):\s*(#[0-9a-fA-F]{6})", PALETTE_LIGHT)
        assert len(found) >= 9
        for name, value in found:
            ratio = 1.05 / (lum(value) + 0.05)
            assert ratio >= 4.5, f"{name} {value} has contrast {ratio:.2f}"


@pytest.fixture(scope="module")
def html_text(tmp_path_factory):
    return render(tmp_path_factory.mktemp("a11y"))


class TestTabsAndPanels:
    def test_tablist_roles(self, html_text):
        assert 'role="tablist"' in html_text
        tabs = re.findall(r'<button[^>]*role="tab"[^>]*>', html_text)
        assert len(tabs) >= 4
        selected = [t for t in tabs if 'aria-selected="true"' in t]
        assert len(selected) == 1
        for t in tabs:
            controls = re.search(r'aria-controls="([^"]+)"', t).group(1)
            assert re.search(
                rf'<section id="{controls}" class="tab-section[^"]*" role="tabpanel" '
                rf'aria-labelledby="tab-{controls}"',
                html_text,
            )

    def test_inactive_tabs_not_in_tab_order(self, html_text):
        tabs = re.findall(r'<button[^>]*role="tab"[^>]*>', html_text)
        for t in tabs:
            expected = "0" if 'aria-selected="true"' in t else "-1"
            assert f'tabindex="{expected}"' in t

    def test_keyboard_and_hash_handlers_present(self, html_text):
        for token in ("ArrowRight", "ArrowLeft", "Home", "End", "hashchange", "history.replaceState"):
            assert token in html_text

    def test_collapsible_headers_are_buttons(self, html_text):
        headers = re.findall(r'<button[^>]*class="collapsible-header"[^>]*>', html_text)
        assert headers
        for h in headers:
            assert "aria-expanded=" in h and "aria-controls=" in h
        assert '<div class="collapsible-header"' not in html_text

    def test_finding_links_are_anchors_with_href(self, html_text):
        for m in re.finditer(r'<a class="finding-link"[^>]*>', html_text):
            assert "href=" in m.group(0)
