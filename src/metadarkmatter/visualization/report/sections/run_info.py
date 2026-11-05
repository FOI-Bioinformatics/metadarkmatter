"""ReportGenerator "Run information" panel (Summary tab)."""

from __future__ import annotations

import logging
from typing import Any

import polars as pl

from metadarkmatter.visualization.report.components import count_missing_pairs
from metadarkmatter.visualization.report.provenance import (
    INPUT_ROLES,
    _e,
    file_row_value,
    format_setting,
    render_rows,
)
from metadarkmatter.visualization.report.sections._base import _ReportBase
from metadarkmatter.visualization.report.templates import COLLAPSIBLE_PANEL_TEMPLATE

logger = logging.getLogger(__name__)

try:
    from metadarkmatter import __version__ as _REPORT_VERSION
except ImportError:  # pragma: no cover
    _REPORT_VERSION = "unknown"

# ScoringConfig fields shown as "filters"; value 0 / None means disabled.
_FILTER_FIELDS: tuple[tuple[str, str], ...] = (
    ("min_percent_identity", "Min. percent identity"),
    ("min_bitscore", "Min. bitscore"),
    ("max_evalue", "Max. E-value"),
    ("min_alignment_length", "Min. alignment length"),
    ("min_alignment_fraction", "Min. alignment fraction"),
    ("min_read_length", "Min. read length"),
    ("min_query_coverage", "Min. query coverage"),
    ("target_family", "Target family"),
)


class RunInfoMixin(_ReportBase):
    """Builds the collapsible Run information panel."""

    def _run_info_settings_rows(self, cfg: dict[str, Any]) -> list[tuple[str, str]]:
        rows: list[tuple[str, str]] = []

        def add(label: str, key: str) -> None:
            if key in cfg:
                rows.append((label, _e(format_setting(cfg[key]))))

        add("Alignment mode", "alignment_mode")
        add("Uncertainty mode", "uncertainty_mode")
        if "coverage_weight_mode" in cfg:
            mode = cfg["coverage_weight_mode"]
            text = format_setting(mode)
            if mode != "none" and "coverage_weight_strength" in cfg:
                text += f" (strength {format_setting(cfg['coverage_weight_strength'])})"
            rows.append(("Coverage weighting", _e(text)))
        if "clade_novel_species" in cfg:
            rows.append((
                "Clade-novelty rule",
                "on" if cfg["clade_novel_species"] else "off",
            ))
        add("Bitscore window (% of top hit)", "bitscore_threshold_pct")
        active = [
            f"{label} = {format_setting(cfg[key])}"
            for key, label in _FILTER_FIELDS
            if key in cfg and cfg[key] not in (0, 0.0, None, "")
        ]
        if any(key in cfg for key, _ in _FILTER_FIELDS):
            rows.append((
                "Hit filters",
                _e("; ".join(active)) if active else "none (all thresholds disabled)",
            ))
        return rows

    def _derived_settings_rows(self) -> list[tuple[str, str]]:
        """What can be said about the run from the data alone."""
        rows: list[tuple[str, str]] = [(
            "Alignment mode (from --alignment-mode)",
            _e(self.config.alignment_mode),
        )]
        if "clade_novel_species" in self.df.columns:
            flagged = int(self.df["clade_novel_species"].cast(pl.Boolean, strict=False).sum() or 0)
            rows.append(("Clade-novelty rule", f"column present; {flagged:,} reads flagged"))
        if self.summary.has_family_validation:
            rows.append(("Family validation", "active (family bitscore columns present)"))
        return rows

    def _run_info_read_rows(self, prov: dict[str, Any] | None) -> list[tuple[str, str]]:
        s = self.summary
        rows: list[tuple[str, str]] = [("Reads in classification table", f"{s.total_reads:,}")]
        if prov:
            results = prov.get("results")
            declared = results.get("total_reads") if isinstance(results, dict) else None
            if isinstance(declared, int) and declared != s.total_reads:
                rows.append((
                    "Reads declared by provenance",
                    f"{declared:,} (differs from the table: the table may have been "
                    "edited or filtered after classification)",
                ))
        rows.append((
            "Known / Novel / Uncertain",
            f"{s.diversity_known:,} / {s.diversity_novel:,} / {s.diversity_uncertain:,}",
        ))
        if s.off_target:
            rows.append(("Off-target", f"{s.off_target:,}"))
        if "uncertainty_type" in self.df.columns:
            types = self.df["uncertainty_type"].drop_nulls()
            if len(types):
                measured = int((types == "measured").sum())
                inferred = int((types == "inferred").sum())
                n = len(self.df) or 1
                rows.append((
                    "Uncertainty: measured / inferred",
                    f"{measured:,} ({measured / n * 100:.1f}%) / "
                    f"{inferred:,} ({inferred / n * 100:.1f}%)",
                ))
                rows.append((
                    "",
                    "Measured = from ANI between competing genomes; inferred = "
                    "estimated for single-hit reads. Only measured values enter taxonomic_call.",
                ))
        return rows

    def _run_info_matrix_rows(self) -> list[tuple[str, str]]:
        rows: list[tuple[str, str]] = []
        for label, matrix in (("ANI matrix", self.ani_matrix), ("AAI matrix", self.aai_matrix)):
            if matrix is None or len(matrix) == 0:
                continue
            missing, total = count_missing_pairs(matrix)
            n_genomes = len(matrix)
            if total == 0:
                rows.append((f"{label} completeness", f"{n_genomes} genome(s); no pairs to evaluate"))
                continue
            computed = total - missing
            rows.append((
                f"{label} completeness",
                f"{computed:,} of {total:,} genome pairs computed "
                f"({computed / total * 100:.1f}%) across {n_genomes:,} genomes",
            ))
        return rows

    def _build_run_info_panel(self) -> str:
        """Collapsible "Run information" panel for the Summary tab."""
        prov = self.provenance
        parts: list[str] = ['<div class="run-info">']

        if prov:
            parts.append(
                '<p class="mdm-note">Recorded by <code>mdm score classify</code> in the '
                "provenance sidecar next to the classification table.</p>"
            )
        else:
            parts.append(
                '<p class="mdm-note"><strong>Provenance was not available</strong> '
                "(no <code>&lt;classifications&gt;.provenance.json</code> sidecar was found "
                "or passed with <code>--provenance</code>). Only values derived from the "
                "classification table and the files given to the report are shown; "
                "classification settings cannot be verified.</p>"
            )

        # Software
        software: list[tuple[str, str]] = []
        report_version = _e(_REPORT_VERSION)
        if prov and prov.get("metadarkmatter_version"):
            run_version = _e(prov["metadarkmatter_version"])
            software.append(("metadarkmatter (classification)", run_version))
            if run_version != report_version:
                software.append(("metadarkmatter (this report)", report_version))
        else:
            software.append(("metadarkmatter (this report)", report_version))
        if prov:
            for label, key in (("Python", "python"), ("Platform", "platform"), ("Random seed", "seed"),
                               ("Classification created (UTC)", "created_utc")):
                if prov.get(key) is not None:
                    software.append((label, _e(prov[key])))
            command = prov.get("command")
            if isinstance(command, list) and command:
                software.append(("Command", f"<code>{_e(' '.join(map(str, command)))}</code>"))
        parts.append("<h4>Software</h4>" + render_rows(software))

        # Inputs
        inputs: list[tuple[str, str]] = []
        names = self.config.input_names
        prov_inputs = prov.get("inputs") if prov else None
        if isinstance(prov_inputs, dict):
            for key, label in INPUT_ROLES:
                if prov_inputs.get(key):
                    inputs.append((label, file_row_value(prov_inputs[key])))
        if names.get("classifications"):
            inputs.append(("Classification table", f"<code>{_e(names['classifications'])}</code>"))
        recorded = {label for label, _ in inputs}
        for key, label in (("ani", "ANI matrix"), ("aai", "AAI matrix"),
                           ("metadata", "Genome metadata"), ("recruitment", "Recruitment data")):
            if names.get(key) and label not in recorded:
                inputs.append((label, f"<code>{_e(names[key])}</code>"))
        if inputs:
            parts.append("<h4>Input files</h4>" + render_rows(inputs))

        # Settings
        cfg = prov.get("config") if prov else None
        if isinstance(cfg, dict):
            settings = self._run_info_settings_rows(cfg)
            parts.append("<h4>Classification settings</h4>" + render_rows(settings))
        else:
            parts.append("<h4>Classification settings</h4>" + render_rows(self._derived_settings_rows()))

        # Reads and matrix completeness
        data_rows = self._run_info_read_rows(prov) + self._run_info_matrix_rows()
        rows_html = "".join(
            f"<dt>{_e(k)}</dt><dd>{v}</dd>" if k else f"<dt></dt><dd><small>{_e(v)}</small></dd>"
            for k, v in data_rows
        )
        parts.append(f"<h4>Reads and reference data</h4><dl>{rows_html}</dl>")
        parts.append("</div>")

        return COLLAPSIBLE_PANEL_TEMPLATE.format(
            default_state="",
            panel_id="run-info",
            title="Run information",
            content="\n".join(parts),
        )
