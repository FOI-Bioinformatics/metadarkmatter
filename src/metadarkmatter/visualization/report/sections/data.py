"""ReportGenerator data section."""

from __future__ import annotations

import html
import logging
from typing import TYPE_CHECKING

import polars as pl

from metadarkmatter.visualization.report.models import (
    _safe_float,
)
from metadarkmatter.visualization.report.sections._base import _ReportBase
from metadarkmatter.visualization.report.templates import (
    DATA_COLUMN_GUIDE_TEMPLATE,
    DATA_QUICK_FILTERS_TEMPLATE,
    DATA_SUMMARY_TEMPLATE,
    DATA_TABLE_TEMPLATE,
    TAB_SECTION_TEMPLATE,
    TABLE_ROW_TEMPLATE,
    get_cell_class,
)

if TYPE_CHECKING:

    pass


logger = logging.getLogger(__name__)


class DataMixin(_ReportBase):
    """Data report-section builders."""

    def _build_data_section(self) -> str:
        """Build the interactive data table section with improved UX."""
        s = self.summary
        total_reads = len(self.df)

        # Limit rows for performance
        table_df = self.df.head(self.config.max_table_rows)
        rows_truncated = total_reads > self.config.max_table_rows

        # Build summary section with truncation notice if applicable
        truncation_notice = ""
        if rows_truncated:
            truncation_notice = (
                f'<div class="truncation-notice">'
                f'<strong>Note:</strong> Showing first {self.config.max_table_rows:,} of '
                f'{total_reads:,} reads. Export full data for complete analysis.'
                f'</div>'
            )

        # Summary counts reflect the FULL dataset (not truncated table)
        off_target_stat = ""
        if s.off_target > 0:
            off_target_stat = (
                '<div class="data-stat off-target">'
                f'<span class="stat-value">{s.off_target:,}</span>'
                '<span class="stat-label">Off-target</span>'
                '</div>'
            )
        summary_html = DATA_SUMMARY_TEMPLATE.format(
            total_rows=total_reads,
            known_count=s.diversity_known,
            novel_count=s.diversity_novel,
            uncertain_count=s.diversity_uncertain,
            off_target_stat=off_target_stat,
        )

        # Build quick filters section with all classification categories
        off_target_chip = ""
        if s.off_target > 0:
            off_target_chip = (
                f'<button class="filter-chip off-target" data-filter="Off-target" '
                f'onclick="quickFilter(this, \'Off-target\')">Off-target ({s.off_target})</button>'
            )
        quick_filters_html = DATA_QUICK_FILTERS_TEMPLATE.format(
            known_count=s.known_species,
            novel_species_count=s.novel_species,
            novel_genus_count=s.novel_genus,
            species_boundary_count=s.species_boundary,
            ambiguous_count=s.ambiguous,
            ambiguous_wg_count=s.ambiguous_within_genus,
            conserved_count=s.conserved_regions,
            unclassified_count=s.unclassified,
            off_target_chip=off_target_chip,
        )

        # Column guide with similarity type (ANI or AAI)
        column_guide_html = DATA_COLUMN_GUIDE_TEMPLATE.format(
            similarity_type=self.config.similarity_type,
        )

        rows_html = self._build_data_table_rows(table_df)

        table_html = DATA_TABLE_TEMPLATE.format(
            table_rows=rows_html,
            page_size=self.config.page_size,
            total_rows=len(table_df),
        )

        # Add JavaScript to show enhanced scoring columns and guide if data available
        enhanced_js = ""
        if self.summary.has_enhanced_scoring or self.summary.has_inferred_uncertainty:
            enhanced_js = """
<script>
document.addEventListener('DOMContentLoaded', function() {
    // Show enhanced scoring column options
    var divider = document.getElementById('enhancedColsDivider');
    if (divider) divider.style.display = 'block';
    for (var i = 11; i <= 15; i++) {
        var label = document.getElementById('colLabel' + i);
        if (label) label.style.display = 'block';
    }
    // Show enhanced guide section
    var enhancedGuide = document.getElementById('enhancedGuide');
    if (enhancedGuide) enhancedGuide.style.display = 'block';
});
</script>
"""

        # Combine all sections
        content = (
            truncation_notice + summary_html + quick_filters_html +
            column_guide_html + table_html + enhanced_js
        )

        return TAB_SECTION_TEMPLATE.format(
            tab_id="data",
            active_class="",
            section_title="Classification Data",
            content=content,
        )

    def _build_data_table_rows(self, table_df: pl.DataFrame) -> str:
        """Render the data-table body rows for the given (truncated) frame."""
        # Build table rows with additional columns
        rows_html = ""
        for row in table_df.iter_rows(named=True):
            genome = row.get("best_match_genome", "")

            # Get species and genus from metadata if available
            species = ""
            genus = ""
            if self.genome_metadata is not None:
                species = self.genome_metadata.get_species(genome) or ""
                genus = self.genome_metadata.get_genus(genome) or ""

            # Get additional columns from dataframe
            ambiguous_hits = row.get("num_ambiguous_hits", 0)
            identity_gap = _safe_float(row.get("identity_gap"))
            # Format identity_gap: use absolute value, show "-" if not available
            identity_gap_str = "-" if identity_gap is None else f"{abs(identity_gap):.2f}"
            is_novel = row.get("is_novel", False)
            if is_novel is None:
                # Compute from taxonomic_call if not in dataframe
                call = row.get("taxonomic_call", "")
                is_novel = call in ("Novel Species", "Novel Genus")

            # Enhanced scoring fields (with defaults for when not available)
            # Values may be None, float, or string (when Polars infers
            # mostly-empty columns as Utf8), so coerce before formatting.
            uncertainty_type = row.get("uncertainty_type", "-")
            inferred_unc = _safe_float(row.get("inferred_uncertainty"))
            inferred_unc_str = "-" if inferred_unc is None else f"{inferred_unc:.1f}"
            discovery = _safe_float(row.get("discovery_score"))
            discovery_str = "-" if discovery is None else f"{discovery:.1f}"
            id_conf = _safe_float(row.get("identity_confidence"))
            id_conf_str = "-" if id_conf is None else f"{id_conf:.1f}"
            pl_conf = _safe_float(row.get("placement_confidence"))
            pl_conf_str = "-" if pl_conf is None else f"{pl_conf:.1f}"

            rows_html += TABLE_ROW_TEMPLATE.format(
                read_id=html.escape(str(row.get("read_id", ""))),
                best_match=html.escape(
                    genome[:40] if len(genome) > 40 else genome
                ),
                species=html.escape(
                    species[:30] if len(species) > 30 else species
                ),
                genus=html.escape(
                    genus[:20] if len(genus) > 20 else genus
                ),
                identity=row.get("top_hit_identity", 0),
                novelty=row.get("novelty_index", 0),
                uncertainty=row.get("placement_uncertainty", 0),
                ambiguous_hits=ambiguous_hits,
                identity_gap=identity_gap_str,
                is_novel="Yes" if is_novel else "No",
                classification=html.escape(
                    str(row.get("taxonomic_call", ""))
                ),
                cell_class=get_cell_class(row.get("taxonomic_call", "")),
                uncertainty_type=html.escape(
                    str(uncertainty_type or "-")
                ),
                inferred_uncertainty=inferred_unc_str,
                discovery_score=discovery_str,
                identity_confidence=id_conf_str,
                placement_confidence=pl_conf_str,
            )

        return rows_html
