"""ReportGenerator data section.

The table itself is rendered client-side from a compact embedded JSON
payload (see :mod:`metadarkmatter.visualization.report.data_table`).
"""

from __future__ import annotations

import html
import logging

from metadarkmatter.visualization.report.data_table import (
    SAMPLING_DESCRIPTION,
    active_columns,
    build_table_payload,
    select_table_rows,
)
from metadarkmatter.visualization.report.embed import json_script_tag
from metadarkmatter.visualization.report.sections._base import _ReportBase
from metadarkmatter.visualization.report.templates import (
    DATA_COLUMN_GUIDE_TEMPLATE,
    DATA_QUICK_FILTERS_TEMPLATE,
    DATA_SUMMARY_TEMPLATE,
    DATA_TABLE_TEMPLATE,
    TAB_SECTION_TEMPLATE,
)

logger = logging.getLogger(__name__)


class DataMixin(_ReportBase):
    """Data report-section builders."""

    def _build_truncation_notice(self, shown: int, total: int) -> str:
        """Notice for a sampled table (empty when every read is shown)."""
        if shown >= total:
            return ""
        source = self.config.classifications_name or "the classifications file"
        return (
            '<div class="truncation-notice" role="note">'
            f"<strong>Table truncated:</strong> showing a deterministic sample of "
            f"{shown:,} of {total:,} reads ({html.escape(SAMPLING_DESCRIPTION)}; "
            f"original order preserved). Summary counts and plots use all {total:,} reads. "
            f"The complete per-read results are in the classification output "
            f"(<code>{html.escape(source)}</code>); raise or remove "
            "<code>--max-table-rows</code> to embed more rows."
            "</div>"
        )

    def _build_data_section(self) -> str:
        """Build the interactive data table section."""
        s = self.summary
        total_reads = len(self.df)

        table_df = select_table_rows(self.df, self.config.max_table_rows)
        truncation_notice = self._build_truncation_notice(len(table_df), total_reads)

        # Summary counts reflect the FULL dataset (not a truncated table)
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

        off_target_chip = ""
        if s.off_target > 0:
            off_target_chip = (
                '<button class="filter-chip off-target" aria-pressed="false" '
                'data-filter="Off-target" '
                f"onclick=\"quickFilter(this, 'Off-target')\">Off-target ({s.off_target})</button>"
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

        column_guide_html = DATA_COLUMN_GUIDE_TEMPLATE.format(
            similarity_type=self.config.similarity_type,
        )

        columns = active_columns(table_df)
        column_menu = "\n".join(
            f'                        <label><input type="checkbox" data-col="{c.key}"'
            f'{" checked" if c.default_visible else ""} '
            f"onchange=\"toggleColumn('{c.key}')\"> {html.escape(c.label)}</label>"
            for c in columns
        )
        hidden_style = ' style="display:none;"'
        table_head = "\n".join(
            f'                    <th scope="col" data-col="{c.key}" aria-sort="none"'
            f'{"" if c.default_visible else hidden_style}>'
            f'<button type="button" class="th-sort" onclick="sortTable(\'{c.key}\')">'
            f'{html.escape(c.label)} <span class="sort-icon" aria-hidden="true">-</span>'
            "</button></th>"
            for c in columns
        )

        payload = build_table_payload(
            table_df,
            self.genome_metadata,
            page_size=self.config.page_size,
            total_reads=total_reads,
        )
        table_html = DATA_TABLE_TEMPLATE.format(
            column_menu=column_menu,
            table_head=table_head,
            page_size=self.config.page_size,
            total_rows=len(table_df),
            table_data=json_script_tag("mdm-table-data", payload),
        )

        content = (
            truncation_notice + summary_html + quick_filters_html
            + column_guide_html + table_html
        )

        return TAB_SECTION_TEMPLATE.format(
            tab_id="data",
            active_class="",
            section_title="Classification Data",
            content=content,
        )
