"""ReportGenerator summary section."""

from __future__ import annotations

import contextlib
import logging
from typing import TYPE_CHECKING

import polars as pl

from metadarkmatter.visualization.plots.base import (
    format_count,
)
from metadarkmatter.visualization.report.models import (
    TaxonomicSummary,
    _BayesianMetrics,
    _EnhancedScoringMetrics,
    _safe_float,
)
from metadarkmatter.visualization.report.sections._base import _ReportBase
from metadarkmatter.visualization.report.templates import (
    CATEGORY_GRID_CELL,
    CATEGORY_GRID_TEMPLATE,
    KPI_CARD_TEMPLATE,
    KPI_STRIP_TEMPLATE,
    METRIC_CARD_TEMPLATE,
    METRIC_CARDS_CONTAINER,
    OVERVIEW_FINDING_CARD_TEMPLATE,
    OVERVIEW_KEY_FINDINGS_TEMPLATE,
    TAB_SECTION_TEMPLATE,
    TOP_SPECIES_ROW_TEMPLATE,
    TOP_SPECIES_TABLE_TEMPLATE,
    TWO_COLUMN_ROW_TEMPLATE,
)

if TYPE_CHECKING:

    pass


logger = logging.getLogger(__name__)


class SummaryMixin(_ReportBase):
    """Summary report-section builders."""

    def _compute_summary(self) -> TaxonomicSummary:
        """Compute summary statistics from classifications."""
        # Count by classification
        counts = (
            self.df.group_by("taxonomic_call")
            .len()
            .to_dict(as_series=False)
        )

        count_map = dict(zip(counts.get("taxonomic_call", []),
                            counts.get("len", []), strict=True))

        # Mean metrics
        mean_novelty = _safe_float(self.df["novelty_index"].mean()) or 0.0
        mean_uncertainty = _safe_float(self.df["placement_uncertainty"].mean()) or 0.0
        mean_identity = _safe_float(self.df["top_hit_identity"].mean()) or 0.0

        # Compute diversity groupings
        known_species = count_map.get("Known Species", 0)
        novel_species = count_map.get("Novel Species", 0)
        novel_genus = count_map.get("Novel Genus", 0)
        species_boundary = count_map.get("Species Boundary", 0)
        ambiguous = count_map.get("Ambiguous", 0)
        ambiguous_within_genus = count_map.get("Ambiguous Within Genus", 0)
        conserved_regions = count_map.get("Conserved Region", 0)
        unclassified = count_map.get("Unclassified", 0)
        off_target = count_map.get("Off-target", 0)

        # Family validation detection
        has_family_validation = "family_bitscore_ratio" in self.df.columns
        target_family = ""

        # High-level diversity grouping
        diversity_known = known_species
        diversity_novel = novel_species + novel_genus
        diversity_uncertain = (
            species_boundary + ambiguous + ambiguous_within_genus +
            conserved_regions + unclassified
        )
        diversity_off_target = off_target

        # Check for enhanced scoring columns.  The column must exist AND
        # have non-null numeric data for classified (non-Off-target) reads;
        # Off-target rows always carry placeholder zeros which should not
        # trigger the "Discovery Scores" tab.
        has_enhanced_scoring = False
        if "alignment_quality" in self.df.columns:
            classified_df = self.df.filter(pl.col("taxonomic_call") != "Off-target")
            aq_col = classified_df["alignment_quality"]
            # Cast to Float64 in case Polars inferred String from CSV
            with contextlib.suppress(Exception):
                aq_col = aq_col.cast(pl.Float64, strict=False)
            has_enhanced_scoring = aq_col.drop_nulls().len() > 0
        has_inferred_uncertainty = "inferred_uncertainty" in self.df.columns

        # Calculate single-hit statistics
        single_hit_count = 0
        single_hit_pct = 0.0
        if "num_ambiguous_hits" in self.df.columns:
            single_hit_df = self.df.filter(pl.col("num_ambiguous_hits") <= 1)
            single_hit_count = len(single_hit_df)
            single_hit_pct = (single_hit_count / len(self.df) * 100) if len(self.df) > 0 else 0.0

        enh = self._compute_enhanced_scoring_metrics(
            has_enhanced_scoring, has_inferred_uncertainty
        )
        mean_inferred_uncertainty = enh.mean_inferred_uncertainty
        mean_alignment_quality = enh.mean_alignment_quality
        mean_identity_confidence = enh.mean_identity_confidence
        mean_placement_confidence = enh.mean_placement_confidence
        mean_discovery_score = enh.mean_discovery_score
        novel_with_discovery_score = enh.novel_with_discovery_score
        high_priority_discoveries = enh.high_priority_discoveries

        # Check for Bayesian posterior columns
        has_bayesian = "posterior_entropy" in self.df.columns
        bay = self._compute_bayesian_metrics(has_bayesian)
        mean_posterior_entropy = bay.mean_posterior_entropy
        high_confidence_count = bay.high_confidence_count
        high_confidence_pct = bay.high_confidence_pct
        boundary_count = bay.boundary_count
        boundary_pct = bay.boundary_pct
        map_agreement_count = bay.map_agreement_count
        map_agreement_pct = bay.map_agreement_pct

        return TaxonomicSummary(
            total_reads=len(self.df),
            known_species=known_species,
            novel_species=novel_species,
            novel_genus=novel_genus,
            species_boundary=species_boundary,
            ambiguous=ambiguous,
            ambiguous_within_genus=ambiguous_within_genus,
            conserved_regions=conserved_regions,
            unclassified=unclassified,
            diversity_known=diversity_known,
            diversity_novel=diversity_novel,
            diversity_uncertain=diversity_uncertain,
            diversity_off_target=diversity_off_target,
            mean_novelty_index=mean_novelty,
            mean_placement_uncertainty=mean_uncertainty,
            mean_top_hit_identity=mean_identity,
            # Family validation
            off_target=off_target,
            has_family_validation=has_family_validation,
            target_family=target_family,
            # Enhanced scoring
            has_enhanced_scoring=has_enhanced_scoring,
            has_inferred_uncertainty=has_inferred_uncertainty,
            single_hit_count=single_hit_count,
            single_hit_pct=single_hit_pct,
            mean_inferred_uncertainty=mean_inferred_uncertainty,
            mean_alignment_quality=mean_alignment_quality,
            mean_identity_confidence=mean_identity_confidence,
            mean_placement_confidence=mean_placement_confidence,
            mean_discovery_score=mean_discovery_score,
            novel_with_discovery_score=novel_with_discovery_score,
            high_priority_discoveries=high_priority_discoveries,
            # Bayesian
            has_bayesian=has_bayesian,
            mean_posterior_entropy=mean_posterior_entropy,
            high_confidence_count=high_confidence_count,
            high_confidence_pct=high_confidence_pct,
            boundary_count=boundary_count,
            boundary_pct=boundary_pct,
            map_agreement_count=map_agreement_count,
            map_agreement_pct=map_agreement_pct,
        )

    def _compute_enhanced_scoring_metrics(
        self, has_enhanced_scoring: bool, has_inferred_uncertainty: bool
    ) -> _EnhancedScoringMetrics:
        """Compute mean enhanced-scoring / discovery metrics."""
        # Calculate enhanced scoring metrics if available
        mean_inferred_uncertainty = None
        mean_alignment_quality = None
        mean_identity_confidence = None
        mean_placement_confidence = None
        mean_discovery_score = None
        novel_with_discovery_score = 0
        high_priority_discoveries = 0

        if has_inferred_uncertainty:
            inferred_df = self.df.filter(pl.col("inferred_uncertainty").is_not_null())
            if len(inferred_df) > 0:
                mean_inferred_uncertainty = _safe_float(inferred_df["inferred_uncertainty"].mean()) or 0.0

        if has_enhanced_scoring:
            if "alignment_quality" in self.df.columns:
                mean_alignment_quality = _safe_float(self.df["alignment_quality"].mean()) or 0.0
            if "identity_confidence" in self.df.columns:
                mean_identity_confidence = _safe_float(self.df["identity_confidence"].mean()) or 0.0
            if "placement_confidence" in self.df.columns:
                mean_placement_confidence = _safe_float(self.df["placement_confidence"].mean()) or 0.0
            if "discovery_score" in self.df.columns:
                discovery_df = self.df.filter(pl.col("discovery_score").is_not_null())
                novel_with_discovery_score = len(discovery_df)
                if novel_with_discovery_score > 0:
                    mean_discovery_score = _safe_float(discovery_df["discovery_score"].mean()) or 0.0
                    high_priority_discoveries = len(
                        discovery_df.filter(pl.col("discovery_score") >= 75)
                    )

        return _EnhancedScoringMetrics(
            mean_inferred_uncertainty=mean_inferred_uncertainty,
            mean_alignment_quality=mean_alignment_quality,
            mean_identity_confidence=mean_identity_confidence,
            mean_placement_confidence=mean_placement_confidence,
            mean_discovery_score=mean_discovery_score,
            novel_with_discovery_score=novel_with_discovery_score,
            high_priority_discoveries=high_priority_discoveries,
        )

    def _compute_bayesian_metrics(self, has_bayesian: bool) -> _BayesianMetrics:
        """Compute mean posterior entropy and confidence/agreement counts."""
        mean_posterior_entropy = 0.0
        high_confidence_count = 0
        high_confidence_pct = 0.0
        boundary_count = 0
        boundary_pct = 0.0
        map_agreement_count = 0
        map_agreement_pct = 0.0

        if has_bayesian:
            total_n = len(self.df) if len(self.df) > 0 else 1
            mean_posterior_entropy = _safe_float(self.df["posterior_entropy"].mean()) or 0.0
            high_confidence_count = len(
                self.df.filter(pl.col("posterior_entropy") < 1.0)
            )
            high_confidence_pct = high_confidence_count / total_n * 100
            boundary_count = len(
                self.df.filter(pl.col("posterior_entropy") > 1.5)
            )
            boundary_pct = boundary_count / total_n * 100
            # In the Bayesian-primary workflow, taxonomic_call IS the Bayesian
            # MAP + Stage 2 result, so agreement with legacy_call shows how
            # often the two classifiers concur.
            if "legacy_call" in self.df.columns:
                map_agreement_count = len(
                    self.df.filter(
                        pl.col("legacy_call") == pl.col("taxonomic_call")
                    )
                )
                map_agreement_pct = map_agreement_count / total_n * 100

        return _BayesianMetrics(
            mean_posterior_entropy=mean_posterior_entropy,
            high_confidence_count=high_confidence_count,
            high_confidence_pct=high_confidence_pct,
            boundary_count=boundary_count,
            boundary_pct=boundary_pct,
            map_agreement_count=map_agreement_count,
            map_agreement_pct=map_agreement_pct,
        )

    def _build_summary_section(self) -> str:
        """Build the Summary tab with KPI strip, diversity bar, findings, and category grid."""
        s = self.summary
        total = s.total_reads if s.total_reads > 0 else 1

        # 1. KPI Strip
        kpi_cards = []
        kpi_cards.append(KPI_CARD_TEMPLATE.format(
            color_class="accent", value=format_count(s.total_reads), label="Total Reads"
        ))
        kpi_cards.append(KPI_CARD_TEMPLATE.format(
            color_class="known", value=f"{s.diversity_known_pct:.1f}%", label="Known"
        ))
        kpi_cards.append(KPI_CARD_TEMPLATE.format(
            color_class="novel", value=f"{s.diversity_novel_pct:.1f}%", label="Novel"
        ))
        kpi_cards.append(KPI_CARD_TEMPLATE.format(
            color_class="uncertain", value=f"{s.diversity_uncertain_pct:.1f}%", label="Uncertain"
        ))
        if s.diversity_off_target > 0:
            kpi_cards.append(KPI_CARD_TEMPLATE.format(
                color_class="off-target",
                value=f"{s.diversity_off_target_pct:.1f}%",
                label="Off-target",
            ))
        kpi_html = KPI_STRIP_TEMPLATE.format(cards="\n".join(kpi_cards))

        # 2. Diversity Bar (prominent stacked bar with legend)
        if s.diversity_off_target > 0:
            ot_pct = s.diversity_off_target_pct
            off_target_bar = (
                f'        <div class="bar-segment off-target" style="width: {ot_pct}%;" '
                f'title="Off-target: {ot_pct:.1f}%"></div>'
            )
            off_target_legend = (
                '        <span class="legend-item">'
                '<span class="legend-dot off-target"></span> Off-target</span>'
            )
        else:
            off_target_bar = ""
            off_target_legend = ""

        bar_html = f'''
    <div class="diversity-bar">
        <div class="bar-segment known" style="width: {s.diversity_known_pct}%;" title="Known: {s.diversity_known_pct:.1f}%"></div>
        <div class="bar-segment novel" style="width: {s.diversity_novel_pct}%;" title="Novel: {s.diversity_novel_pct:.1f}%"></div>
        <div class="bar-segment uncertain" style="width: {s.diversity_uncertain_pct}%;" title="Uncertain: {s.diversity_uncertain_pct:.1f}%"></div>
{off_target_bar}
    </div>
    <div class="diversity-bar-legend">
        <span class="legend-item"><span class="legend-dot known"></span> Known ({s.diversity_known_pct:.1f}%)</span>
        <span class="legend-item"><span class="legend-dot novel"></span> Novel ({s.diversity_novel_pct:.1f}%)</span>
        <span class="legend-item"><span class="legend-dot uncertain"></span> Uncertain ({s.diversity_uncertain_pct:.1f}%)</span>
{off_target_legend}
    </div>
'''

        # 3. Two-column: Key Findings + Top Species
        key_findings_html = self._build_key_findings()
        top_species_html = self._build_top_species_table()

        two_col_html = TWO_COLUMN_ROW_TEMPLATE.format(
            left_flex="3",
            right_flex="2",
            left_content=key_findings_html,
            right_content=top_species_html,
        )

        # 4. Category Grid (compact 2-row layout replacing verbose breakdown)
        categories = [
            ("Known Species", s.known_species, "known"),
            ("Novel Species", s.novel_species, "novel-species"),
            ("Novel Genus", s.novel_genus, "novel-genus"),
            ("Species Boundary", s.species_boundary, "species-boundary"),
            ("Ambiguous", s.ambiguous, "ambiguous"),
            ("Conserved", s.conserved_regions, "conserved"),
            ("Unclassified", s.unclassified, "unclassified"),
        ]
        if s.off_target > 0:
            categories.append(("Off-target", s.off_target, "off-target"))

        grid_cells = []
        for name, count, color_class in categories:
            grid_cells.append(CATEGORY_GRID_CELL.format(
                color_class=color_class,
                name=name,
                count=count,
                pct=count / total * 100,
            ))
        category_html = CATEGORY_GRID_TEMPLATE.format(cells="\n".join(grid_cells))

        content = kpi_html + bar_html + two_col_html + category_html

        return TAB_SECTION_TEMPLATE.format(
            tab_id="summary",
            active_class="active",
            section_title="Summary",
            content=content,
        )

    def _build_top_species_table(self) -> str:
        """Build a compact top species mini-table for the summary tab."""
        if self.genome_metadata is None:
            return ""

        # Enrich with species info
        if "species" not in self.df.columns:
            enriched_df = self.genome_metadata.join_classifications(self.df)
        else:
            enriched_df = self.df

        # Exclude off-target reads
        if "taxonomic_call" in enriched_df.columns:
            base_df = enriched_df.filter(pl.col("taxonomic_call") != "Off-target")
        else:
            base_df = enriched_df

        if "species" not in base_df.columns:
            return ""

        total = len(base_df) if len(base_df) > 0 else 1
        species_counts = (
            base_df.group_by("species")
            .len()
            .sort("len", descending=True)
            .head(8)
        )

        if len(species_counts) == 0:
            return ""

        rows = []
        for row in species_counts.iter_rows(named=True):
            rows.append(TOP_SPECIES_ROW_TEMPLATE.format(
                species=row["species"],
                count=f'{row["len"]:,}',
                pct=row["len"] / total * 100,
            ))

        return TOP_SPECIES_TABLE_TEMPLATE.format(rows="\n".join(rows))

    def _build_key_findings(self) -> str:
        """Build key findings cards for the overview tab."""
        s = self.summary
        cards = []

        # Top species card (if metadata available)
        if self.genome_metadata is not None and "species" not in self.df.columns:
            enriched_df = self.genome_metadata.join_classifications(self.df)
        elif "species" in self.df.columns:
            enriched_df = self.df
        else:
            enriched_df = None

        if enriched_df is not None and "species" in enriched_df.columns:
            # Exclude off-target reads from top species computation
            if "taxonomic_call" in enriched_df.columns:
                enriched_df = enriched_df.filter(
                    pl.col("taxonomic_call") != "Off-target"
                )
            species_counts = (
                enriched_df.group_by("species")
                .len()
                .sort("len", descending=True)
            )
            if len(species_counts) > 0:
                top_sp = species_counts["species"][0]
                top_count = species_counts["len"][0]
                top_pct = top_count / s.total_reads * 100 if s.total_reads > 0 else 0
                cards.append(OVERVIEW_FINDING_CARD_TEMPLATE.format(
                    card_class="species",
                    icon="S",
                    headline=f"Dominant species: {top_sp}",
                    detail=f"{top_pct:.0f}% of reads ({top_count:,})",
                    link="",
                ))

        # Novel diversity headline
        if s.diversity_novel > 0:
            novel_detail = (
                f"{s.novel_species:,} novel species + "
                f"{s.novel_genus:,} novel genus reads detected"
            )
            link_html = (
                '<a class="finding-link" onclick="showTab(\'novel-diversity\')">'
                'See Novel Diversity tab</a>'
            )
            cards.append(OVERVIEW_FINDING_CARD_TEMPLATE.format(
                card_class="novel",
                icon="N",
                headline=f"{s.diversity_novel:,} novel diversity reads ({s.diversity_novel_pct:.1f}%)",
                detail=novel_detail,
                link=link_html,
            ))

        # Classification confidence
        if s.has_bayesian:
            cards.append(OVERVIEW_FINDING_CARD_TEMPLATE.format(
                card_class="confidence",
                icon="Q",
                headline=f"{s.high_confidence_pct:.0f}% high-confidence classifications",
                detail=f"{s.high_confidence_count:,} reads with entropy < 1.0",
                link="",
            ))

        action_card = self._build_overview_action_items_card()
        if action_card:
            cards.append(action_card)

        if not cards:
            return ""

        return OVERVIEW_KEY_FINDINGS_TEMPLATE.format(
            findings_cards="\n".join(cards),
        )

    def _build_overview_action_items_card(self) -> str:
        """Build the conditional action-items finding card ("" if none)."""
        s = self.summary
        # Conditional action items
        action_items = []
        if s.high_priority_discoveries > 0:
            action_items.append(
                f"{s.high_priority_discoveries} high-priority novel reads "
                f"(discovery score >= 75) warrant experimental validation"
            )
        if s.has_family_validation and s.off_target > 0:
            off_target_pct = s.off_target / s.total_reads * 100 if s.total_reads > 0 else 0
            if off_target_pct > 5:
                action_items.append(
                    f"{off_target_pct:.0f}% off-target reads detected"
                )
            # Count off-target reads with species/genus-level in-family divergence
            if "in_family_novelty_index" in self.df.columns:
                th = self.config.thresholds
                off_target_novel = self.df.filter(
                    (pl.col("taxonomic_call") == "Off-target")
                    & (pl.col("in_family_novelty_index") >= th.novelty_novel_species_min)
                    & (pl.col("in_family_novelty_index") <= th.novelty_novel_genus_max)
                )
                if len(off_target_novel) > 0:
                    action_items.append(
                        f"{len(off_target_novel)} off-target reads show "
                        f"species/genus-level divergence from in-family references"
                    )

        if action_items:
            link_html = ""
            if s.has_family_validation and s.off_target > 0:
                link_html = (
                    '<a class="finding-link" onclick="showTab(\'data\')">'
                    "Review off-target reads in Data tab</a>"
                )
            return OVERVIEW_FINDING_CARD_TEMPLATE.format(
                card_class="action",
                icon="!",
                headline="Action items",
                detail="; ".join(action_items),
                link=link_html,
            )
        return ""

    def _build_metric_cards(self) -> str:
        """Build the metric cards HTML."""
        s = self.summary

        cards_html = ""

        # Total reads
        cards_html += METRIC_CARD_TEMPLATE.format(
            value=format_count(s.total_reads),
            label="Total Reads",
            subtext="Classified reads",
            color_class="",
        )

        # Known species
        known_pct = (s.known_species / s.total_reads * 100) if s.total_reads else 0
        cards_html += METRIC_CARD_TEMPLATE.format(
            value=format_count(s.known_species),
            label="Known Species",
            subtext=f"{known_pct:.1f}% of reads",
            color_class="success",
        )

        # Novel species
        novel_sp_pct = (s.novel_species / s.total_reads * 100) if s.total_reads else 0
        cards_html += METRIC_CARD_TEMPLATE.format(
            value=format_count(s.novel_species),
            label="Novel Species",
            subtext=f"{novel_sp_pct:.1f}% of reads",
            color_class="warning",
        )

        # Novel genus
        novel_gen_pct = (s.novel_genus / s.total_reads * 100) if s.total_reads else 0
        cards_html += METRIC_CARD_TEMPLATE.format(
            value=format_count(s.novel_genus),
            label="Novel Genus",
            subtext=f"{novel_gen_pct:.1f}% of reads",
            color_class="danger",
        )

        # Ambiguous
        ambig_pct = (s.ambiguous / s.total_reads * 100) if s.total_reads else 0
        cards_html += METRIC_CARD_TEMPLATE.format(
            value=format_count(s.ambiguous),
            label="Ambiguous",
            subtext=f"{ambig_pct:.1f}% of reads",
            color_class="muted",
        )

        # Novel diversity (combined)
        cards_html += METRIC_CARD_TEMPLATE.format(
            value=f"{s.novel_percentage:.1f}%",
            label="Novel Diversity",
            subtext="Novel species + genus",
            color_class="info",
        )

        # Mean novelty index
        cards_html += METRIC_CARD_TEMPLATE.format(
            value=f"{s.mean_novelty_index:.2f}",
            label="Mean Novelty",
            subtext="Average novelty index",
            color_class="",
        )

        return METRIC_CARDS_CONTAINER.format(cards=cards_html)
