"""
Report generator for unified HTML reports.

Combines multiple plot components into a single, self-contained HTML report
with tabbed navigation, interactive data tables, and professional styling.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

import polars as pl

logger = logging.getLogger(__name__)

from metadarkmatter.core.io_utils import read_dataframe
from metadarkmatter.core.metadata import GenomeMetadata
from metadarkmatter.visualization.plots.base import (
    BAYESIAN_CATEGORY_COLORS,
    PlotConfig,
    ThresholdConfig,
    format_count,
)
from metadarkmatter.visualization.report.components import (
    build_aai_heatmap,
    build_aai_stats_cards,
    build_ani_heatmap,
    build_ani_stats_cards,
    build_phylogenetic_context_heatmap,
)
from metadarkmatter.visualization.plots.classification_charts import (
    ClassificationBarChart,
    ClassificationDonutChart,
    DiversityDonutChart,
    DiversitySunburstChart,
)
from metadarkmatter.visualization.plots.distributions import (
    NoveltyHistogram,
    UncertaintyHistogram,
)
from metadarkmatter.visualization.plots.scatter_2d import (
    ConfidenceNoveltyScatter,
    NoveltyUncertaintyScatter,
)
from metadarkmatter.visualization.report.styles import get_css_styles
from metadarkmatter.visualization.report.templates import (
    CATEGORY_GRID_CELL,
    CATEGORY_GRID_TEMPLATE,
    COLLAPSIBLE_PANEL_TEMPLATE,
    DATA_COLUMN_GUIDE_TEMPLATE,
    DATA_QUICK_FILTERS_TEMPLATE,
    DATA_SUMMARY_TEMPLATE,
    DATA_TABLE_JS,
    DATA_TABLE_TEMPLATE,
    EMPTY_SECTION_TEMPLATE,
    ENHANCED_SCORING_SUMMARY_TEMPLATE,
    ENHANCED_SCORING_UNCERTAINTY_TYPES_TEMPLATE,
    KPI_CARD_TEMPLATE,
    KPI_STRIP_TEMPLATE,
    METRIC_CARD_TEMPLATE,
    METRIC_CARDS_CONTAINER,
    METHODS_SECTION_TEMPLATE,
    OVERVIEW_FINDING_CARD_TEMPLATE,
    OVERVIEW_KEY_FINDINGS_TEMPLATE,
    PLOT_CONTAINER_SIMPLE_TEMPLATE,
    PLOT_CONTAINER_TEMPLATE,
    PLOT_ROW_TEMPLATE,
    REPORT_BASE_TEMPLATE,
    TAB_NAVIGATION_JS,
    TAB_SECTION_TEMPLATE,
    TABLE_ROW_TEMPLATE,
    TOP_SPECIES_ROW_TEMPLATE,
    TOP_SPECIES_TABLE_TEMPLATE,
    TWO_COLUMN_ROW_TEMPLATE,
    BAYESIAN_INTERPRETATION_TEMPLATE,
    CONFIDENCE_SUMMARY_TEMPLATE,
    BORDERLINE_ANALYSIS_TEMPLATE,
    FAMILY_VALIDATION_SECTION_TEMPLATE,
    get_cell_class,
)
from metadarkmatter.visualization.report.novel_section import (
    NOVEL_EMPTY_TEMPLATE,
    PHYLOGENETIC_HEATMAP_INTRO_TEMPLATE,
    PHYLOGENETIC_HEATMAP_LEGEND_TEMPLATE,
    build_cluster_scatter_figure,
    build_cluster_table_html,
    build_novel_summary_html,
    build_sunburst_figure,
)

if TYPE_CHECKING:
    import pandas as pd
    import plotly.graph_objects as go


# Version for report footer
try:
    from metadarkmatter import __version__
except ImportError:
    __version__ = "0.1.0"


def _safe_float(value: object) -> float | None:
    """Convert a value to float, returning None for missing/unconvertible values."""
    if value is None:
        return None
    try:
        return float(value)
    except (ValueError, TypeError):
        return None


@dataclass
class TaxonomicSummary:
    """Summary statistics from classification results."""

    total_reads: int = 0
    known_species: int = 0
    novel_species: int = 0
    novel_genus: int = 0
    species_boundary: int = 0
    ambiguous: int = 0
    ambiguous_within_genus: int = 0
    conserved_regions: int = 0
    unclassified: int = 0
    # High-level diversity grouping
    diversity_known: int = 0
    diversity_novel: int = 0
    diversity_uncertain: int = 0
    diversity_off_target: int = 0
    mean_novelty_index: float = 0.0
    mean_placement_uncertainty: float = 0.0
    mean_top_hit_identity: float = 0.0
    # Enhanced scoring metrics (optional, only when --enhanced-scoring used)
    has_enhanced_scoring: bool = False
    has_inferred_uncertainty: bool = False
    single_hit_count: int = 0
    single_hit_pct: float = 0.0
    mean_inferred_uncertainty: float | None = None
    mean_alignment_quality: float | None = None
    mean_identity_confidence: float | None = None
    mean_placement_confidence: float | None = None
    mean_discovery_score: float | None = None
    novel_with_discovery_score: int = 0
    high_priority_discoveries: int = 0  # discovery_score >= 75
    # Family validation metrics (optional, only when family validation active)
    off_target: int = 0
    has_family_validation: bool = False
    target_family: str = ""
    # Bayesian posterior metrics (always present in Bayesian-primary workflow)
    has_bayesian: bool = False
    mean_posterior_entropy: float = 0.0
    high_confidence_count: int = 0  # entropy < 1.0
    high_confidence_pct: float = 0.0
    boundary_count: int = 0  # entropy > 1.5
    boundary_pct: float = 0.0
    map_agreement_count: int = 0
    map_agreement_pct: float = 0.0

    @property
    def novel_percentage(self) -> float:
        """Percentage of reads classified as novel (species or genus)."""
        if self.total_reads == 0:
            return 0.0
        return (self.novel_species + self.novel_genus) / self.total_reads * 100

    @property
    def diversity_known_pct(self) -> float:
        """Percentage of reads with confident known diversity."""
        if self.total_reads == 0:
            return 0.0
        return self.diversity_known / self.total_reads * 100

    @property
    def diversity_novel_pct(self) -> float:
        """Percentage of reads with confident novel diversity."""
        if self.total_reads == 0:
            return 0.0
        return self.diversity_novel / self.total_reads * 100

    @property
    def diversity_uncertain_pct(self) -> float:
        """Percentage of reads with uncertain diversity status."""
        if self.total_reads == 0:
            return 0.0
        return self.diversity_uncertain / self.total_reads * 100

    @property
    def diversity_off_target_pct(self) -> float:
        """Percentage of reads classified as off-target."""
        if self.total_reads == 0:
            return 0.0
        return self.diversity_off_target / self.total_reads * 100

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for template rendering."""
        result = {
            "total_reads": self.total_reads,
            "known_species": self.known_species,
            "novel_species": self.novel_species,
            "novel_genus": self.novel_genus,
            "species_boundary": self.species_boundary,
            "ambiguous": self.ambiguous,
            "ambiguous_within_genus": self.ambiguous_within_genus,
            "conserved_regions": self.conserved_regions,
            "unclassified": self.unclassified,
            "diversity_known": self.diversity_known,
            "diversity_novel": self.diversity_novel,
            "diversity_uncertain": self.diversity_uncertain,
            "diversity_off_target": self.diversity_off_target,
            "mean_novelty_index": self.mean_novelty_index,
            "mean_placement_uncertainty": self.mean_placement_uncertainty,
            "mean_top_hit_identity": self.mean_top_hit_identity,
            "novel_percentage": self.novel_percentage,
            "diversity_known_pct": self.diversity_known_pct,
            "diversity_novel_pct": self.diversity_novel_pct,
            "diversity_uncertain_pct": self.diversity_uncertain_pct,
            "diversity_off_target_pct": self.diversity_off_target_pct,
            # Family validation
            "off_target": self.off_target,
            "has_family_validation": self.has_family_validation,
            "target_family": self.target_family,
            # Enhanced scoring
            "has_enhanced_scoring": self.has_enhanced_scoring,
            "has_inferred_uncertainty": self.has_inferred_uncertainty,
            "single_hit_count": self.single_hit_count,
            "single_hit_pct": self.single_hit_pct,
        }
        # Add enhanced scoring metrics if available
        if self.mean_inferred_uncertainty is not None:
            result["mean_inferred_uncertainty"] = self.mean_inferred_uncertainty
        if self.mean_alignment_quality is not None:
            result["mean_alignment_quality"] = self.mean_alignment_quality
        if self.mean_identity_confidence is not None:
            result["mean_identity_confidence"] = self.mean_identity_confidence
        if self.mean_placement_confidence is not None:
            result["mean_placement_confidence"] = self.mean_placement_confidence
        if self.mean_discovery_score is not None:
            result["mean_discovery_score"] = self.mean_discovery_score
            result["novel_with_discovery_score"] = self.novel_with_discovery_score
            result["high_priority_discoveries"] = self.high_priority_discoveries
        return result


@dataclass
class ReportConfig:
    """Configuration for report generation."""

    sample_name: str = "Sample"
    title: str = "Metadarkmatter Classification Report"
    theme: str = "light"
    page_size: int = 100
    max_table_rows: int = 10000  # Configurable via CLI --max-table-rows
    max_scatter_points: int = 50000
    include_plotlyjs: str = "cdn"  # 'cdn', 'embed', or path

    # Histogram bin sizes (None = auto-calculate based on nbins)
    novelty_bin_size: float | None = None  # e.g., 1.0 for 1% bins
    uncertainty_bin_size: float | None = 1.0  # Default to 1% bins for finer resolution

    # Alignment mode indicator for display
    alignment_mode: str = "nucleotide"  # "nucleotide" (ANI) or "protein" (AAI)

    # Phylogenetic context heatmap limits
    max_phylo_clusters: int = 20  # Maximum novel clusters in heatmap
    max_phylo_references: int = 50  # Maximum reference genomes in heatmap

    # Phylogeny tab options
    skip_phylogeny: bool = False  # Skip phylogeny tab generation
    user_tree_path: Path | None = None  # Path to user-provided Newick tree

    plot_config: PlotConfig = field(default_factory=PlotConfig)
    thresholds: ThresholdConfig = field(default_factory=ThresholdConfig)

    @property
    def similarity_type(self) -> str:
        """Return ANI or AAI based on alignment mode."""
        return "AAI" if self.alignment_mode == "protein" else "ANI"


class ReportGenerator:
    """
    Generates unified HTML reports from classification results.

    Creates a self-contained HTML file with:
    - Overview metrics and charts
    - Distribution histograms
    - 2D scatter plot (novelty vs uncertainty)
    - Recruitment plots (if BAM provided)
    - Genome breakdown charts
    - ANI matrix heatmap (if provided)
    - Interactive data table
    """

    def __init__(
        self,
        classifications: pl.DataFrame,
        config: ReportConfig | None = None,
        recruitment_data: pl.DataFrame | None = None,
        ani_matrix: pl.DataFrame | None = None,
        aai_matrix: pl.DataFrame | None = None,
        genome_metadata: GenomeMetadata | None = None,
    ) -> None:
        """
        Initialize report generator.

        Args:
            classifications: DataFrame with classification results. Must include:
                - read_id: Read identifier
                - best_match_genome: Best matching genome
                - top_hit_identity: Percent identity to best hit
                - novelty_index: 100 - top_hit_identity
                - placement_uncertainty: ANI-based uncertainty
                - taxonomic_call: Classification category
            config: Report configuration
            recruitment_data: Optional recruitment plot data
            ani_matrix: Optional ANI matrix for heatmap
            aai_matrix: Optional AAI matrix for genus-level heatmap
            genome_metadata: Optional genome metadata for species-level aggregation
        """
        self.df = classifications
        self.config = config or ReportConfig()
        self.recruitment_data = recruitment_data
        self.ani_matrix = ani_matrix
        self.aai_matrix = aai_matrix
        self.genome_metadata = genome_metadata

        # Compute summary statistics
        self.summary = self._compute_summary()

        # Store generated plot HTML/JS
        self._plot_data: dict[str, tuple[str, str]] = {}  # {plot_id: (div_html, js)}

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
        mean_novelty = self.df["novelty_index"].mean() or 0.0
        mean_uncertainty = self.df["placement_uncertainty"].mean() or 0.0
        mean_identity = self.df["top_hit_identity"].mean() or 0.0

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
            try:
                aq_col = aq_col.cast(pl.Float64, strict=False)
            except Exception:
                pass
            has_enhanced_scoring = aq_col.drop_nulls().len() > 0
        has_inferred_uncertainty = "inferred_uncertainty" in self.df.columns

        # Calculate single-hit statistics
        single_hit_count = 0
        single_hit_pct = 0.0
        if "num_ambiguous_hits" in self.df.columns:
            single_hit_df = self.df.filter(pl.col("num_ambiguous_hits") <= 1)
            single_hit_count = len(single_hit_df)
            single_hit_pct = (single_hit_count / len(self.df) * 100) if len(self.df) > 0 else 0.0

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
                mean_inferred_uncertainty = inferred_df["inferred_uncertainty"].mean() or 0.0

        if has_enhanced_scoring:
            if "alignment_quality" in self.df.columns:
                mean_alignment_quality = self.df["alignment_quality"].mean() or 0.0
            if "identity_confidence" in self.df.columns:
                mean_identity_confidence = self.df["identity_confidence"].mean() or 0.0
            if "placement_confidence" in self.df.columns:
                mean_placement_confidence = self.df["placement_confidence"].mean() or 0.0
            if "discovery_score" in self.df.columns:
                discovery_df = self.df.filter(pl.col("discovery_score").is_not_null())
                novel_with_discovery_score = len(discovery_df)
                if novel_with_discovery_score > 0:
                    mean_discovery_score = discovery_df["discovery_score"].mean() or 0.0
                    high_priority_discoveries = len(
                        discovery_df.filter(pl.col("discovery_score") >= 75)
                    )

        # Check for Bayesian posterior columns
        has_bayesian = "posterior_entropy" in self.df.columns
        mean_posterior_entropy = 0.0
        high_confidence_count = 0
        high_confidence_pct = 0.0
        boundary_count = 0
        boundary_pct = 0.0
        map_agreement_count = 0
        map_agreement_pct = 0.0

        if has_bayesian:
            total_n = len(self.df) if len(self.df) > 0 else 1
            mean_posterior_entropy = self.df["posterior_entropy"].mean() or 0.0
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

    def _get_genome_label(self, accession: str, max_species_len: int = 25) -> str:
        """
        Get a display label for a genome accession.

        If metadata is available, returns "Accession (Species)" format.
        Otherwise returns just the accession.

        Args:
            accession: Genome accession (e.g., GCF_019852205.1)
            max_species_len: Maximum length for species name before truncation

        Returns:
            Formatted label string
        """
        if self.genome_metadata is None:
            return accession

        # Look up species from metadata
        species = self.genome_metadata.get_species(accession)
        if species and species != "Unknown":
            # Truncate species name if too long
            if len(species) > max_species_len:
                species = species[:max_species_len-3] + "..."
            return f"{accession} ({species})"
        return accession

    def _get_genome_labels_map(self, accessions: list[str]) -> dict[str, str]:
        """
        Create a mapping from accessions to display labels.

        Args:
            accessions: List of genome accessions

        Returns:
            Dictionary mapping accession -> display label
        """
        return {acc: self._get_genome_label(acc) for acc in accessions}

    def generate(self, output_path: Path | str) -> None:
        """
        Generate and save the HTML report.

        Args:
            output_path: Path for output HTML file
        """
        html = self._build_html()

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(html, encoding="utf-8")

    def _build_html(self) -> str:
        """Build the complete HTML report."""
        content_sections = []

        # Summary tab (always, active)
        content_sections.append(self._build_summary_section())

        # Classification tab (always)
        content_sections.append(self._build_classification_section())

        # Novel Diversity tab (conditional)
        if self.summary.diversity_novel > 0:
            content_sections.append(self._build_novel_diversity_section())

        # Reference tab (always)
        content_sections.append(self._build_reference_section())

        # Data table tab (always)
        content_sections.append(self._build_data_section())

        content = "\n".join(content_sections)

        # Build dynamic navigation
        navigation = self._build_navigation()

        # Methods footer (outside tab container)
        methods_footer = self._build_methods_footer()

        # Collect all Plotly JS initialization code
        plotly_js = self._build_plotly_js()

        # Combine JavaScript
        js_scripts = TAB_NAVIGATION_JS + "\n" + DATA_TABLE_JS.format(
            page_size=self.config.page_size
        )

        # Build final HTML
        return REPORT_BASE_TEMPLATE.format(
            title=self.config.title,
            sample_name=self.config.sample_name,
            timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            total_reads=format_count(self.summary.total_reads),
            version=__version__,
            css_styles=get_css_styles(self.config.theme),
            navigation=navigation,
            content=content,
            methods_footer=methods_footer,
            plotly_js=plotly_js,
            js_scripts=js_scripts,
        )

    def _build_navigation(self) -> str:
        """Build dynamic navigation based on available data."""
        tabs = [
            ("summary", "Summary", True),
            ("classification", "Classification", False),
        ]

        if self.summary.diversity_novel > 0:
            tabs.append(("novel-diversity", "Novel Diversity", False))

        tabs.extend([
            ("reference", "Reference", False),
            ("data", "Data", False),
        ])

        nav_items = []
        for tab_id, label, is_active in tabs:
            active_class = " active" if is_active else ""
            nav_items.append(
                f'        <button class="tab-btn{active_class}" '
                f"onclick=\"showTab('{tab_id}')\">{label}</button>"
            )

        return "\n".join(nav_items)

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
            cards.append(OVERVIEW_FINDING_CARD_TEMPLATE.format(
                card_class="action",
                icon="!",
                headline="Action items",
                detail="; ".join(action_items),
                link=link_html,
            ))

        if not cards:
            return ""

        return OVERVIEW_KEY_FINDINGS_TEMPLATE.format(
            findings_cards="\n".join(cards),
        )

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

    def _build_classification_section(self) -> str:
        """Build the Classification tab with scatter, sunburst, histograms, and collapsible confidence."""
        s = self.summary

        # 1. Metric Strip (KPI-style)
        kpi_cards = []
        kpi_cards.append(KPI_CARD_TEMPLATE.format(
            color_class="accent",
            value=f"{s.mean_novelty_index:.1f}",
            label="Mean Novelty",
        ))
        kpi_cards.append(KPI_CARD_TEMPLATE.format(
            color_class="accent",
            value=f"{s.mean_placement_uncertainty:.1f}",
            label="Mean Uncertainty",
        ))
        if s.has_bayesian:
            kpi_cards.append(KPI_CARD_TEMPLATE.format(
                color_class="accent",
                value=f"{s.mean_posterior_entropy:.2f}",
                label="Mean Entropy",
            ))
            kpi_cards.append(KPI_CARD_TEMPLATE.format(
                color_class="known",
                value=f"{s.high_confidence_pct:.0f}%",
                label="High Confidence",
            ))
        kpi_html = KPI_STRIP_TEMPLATE.format(cards="\n".join(kpi_cards))

        # 2. Scatter plot (primary diagnostic)
        scatter = NoveltyUncertaintyScatter(
            self.df,
            config=PlotConfig(width=800, height=550),
            thresholds=self.config.thresholds,
            max_points=self.config.max_scatter_points,
        )
        scatter_id = "plot-scatter-2d"
        self._register_plot(scatter_id, scatter.create_figure())

        scatter_plot = PLOT_CONTAINER_SIMPLE_TEMPLATE.format(
            extra_class="",
            plot_id=scatter_id,
        )

        # 3. Sunburst
        sunburst = DiversitySunburstChart(
            s.to_dict(),
            config=PlotConfig(width=500, height=550),
            title="Diversity Classification",
        )
        sunburst_id = "plot-diversity-sunburst"
        self._register_plot(sunburst_id, sunburst.create_figure())

        sunburst_plot = PLOT_CONTAINER_SIMPLE_TEMPLATE.format(
            extra_class="",
            plot_id=sunburst_id,
        )

        # Two-column: scatter (60%) + sunburst (40%)
        scatter_row = TWO_COLUMN_ROW_TEMPLATE.format(
            left_flex="3",
            right_flex="2",
            left_content=scatter_plot,
            right_content=sunburst_plot,
        )

        # 4. Histograms side-by-side
        novelty_hist = NoveltyHistogram(
            self.df,
            config=PlotConfig(width=550, height=400),
            thresholds=self.config.thresholds,
            bin_size=self.config.novelty_bin_size,
        )
        novelty_id = "plot-novelty-hist"
        self._register_plot(novelty_id, novelty_hist.create_figure())

        uncertainty_hist = UncertaintyHistogram(
            self.df,
            config=PlotConfig(width=550, height=400),
            thresholds=self.config.thresholds,
            bin_size=self.config.uncertainty_bin_size,
        )
        uncertainty_id = "plot-uncertainty-hist"
        self._register_plot(uncertainty_id, uncertainty_hist.create_figure())

        hist_row = PLOT_ROW_TEMPLATE.format(
            plots=(
                PLOT_CONTAINER_SIMPLE_TEMPLATE.format(
                    extra_class="half-width",
                    plot_id=novelty_id,
                )
                + PLOT_CONTAINER_SIMPLE_TEMPLATE.format(
                    extra_class="half-width",
                    plot_id=uncertainty_id,
                )
            )
        )

        # 5. Confidence Analysis (collapsible, default collapsed)
        confidence_content = self._build_confidence_panel_content()
        if confidence_content:
            confidence_panel = COLLAPSIBLE_PANEL_TEMPLATE.format(
                default_state="",
                panel_id="confidence-analysis",
                title="Confidence Analysis",
                content=confidence_content,
            )
        else:
            confidence_panel = ""

        content = kpi_html + scatter_row + hist_row + confidence_panel

        return TAB_SECTION_TEMPLATE.format(
            tab_id="classification",
            active_class="",
            section_title="Classification",
            content=content,
        )

    def _build_confidence_panel_content(self) -> str:
        """Build inner HTML content for the collapsible Confidence Analysis panel.

        Extracts the body content from what was previously the standalone
        Classification Confidence tab. Returns empty string if no Bayesian
        or enhanced scoring data is available.
        """
        import plotly.graph_objects as go

        s = self.summary

        if not (s.has_bayesian or s.has_enhanced_scoring or s.has_inferred_uncertainty):
            return ""

        content_parts: list[str] = []

        # --- Merged summary cards ---
        if s.has_bayesian:
            summary_html = CONFIDENCE_SUMMARY_TEMPLATE.format(
                mean_entropy=s.mean_posterior_entropy,
                high_confidence_pct=s.high_confidence_pct,
                high_confidence_count=s.high_confidence_count,
                map_agreement_pct=s.map_agreement_pct,
                map_agreement_count=s.map_agreement_count,
                total_reads=s.total_reads,
                single_hit_pct=s.single_hit_pct,
                single_hit_count=s.single_hit_count,
                mean_discovery_score=s.mean_discovery_score or 0.0,
                novel_count=s.novel_with_discovery_score,
                boundary_pct=s.boundary_pct,
                boundary_count=s.boundary_count,
            )
            content_parts.append(summary_html)
            content_parts.append(BAYESIAN_INTERPRETATION_TEMPLATE)
        elif s.has_enhanced_scoring or s.has_inferred_uncertainty:
            # Fallback: show discovery-only summary when no Bayesian data
            summary_html = ENHANCED_SCORING_SUMMARY_TEMPLATE.format(
                single_hit_pct=s.single_hit_pct,
                single_hit_count=s.single_hit_count,
                total_reads=s.total_reads,
                mean_inferred_uncertainty=s.mean_inferred_uncertainty or 0.0,
                high_priority_count=s.high_priority_discoveries,
                mean_discovery_score=s.mean_discovery_score or 0.0,
                novel_count=s.novel_with_discovery_score,
            )
            content_parts.append(summary_html)

        # --- Posterior Entropy Distribution ---
        if s.has_bayesian and "posterior_entropy" in self.df.columns:
            entropy_id = "plot-bayesian-entropy"
            self._build_histogram(
                data=self.df["posterior_entropy"].to_list(),
                plot_id=entropy_id,
                title="Posterior Entropy Distribution",
                x_label="Shannon Entropy (bits)",
                nbins=40,
                thresholds=[
                    (0.5, "dash", "#22c55e", "Very confident"),
                    (1.0, "dot", "#f59e0b", "Confident"),
                    (1.5, "dot", "#ef4444", "Uncertain"),
                ],
                x_range=[0, 2.7],
            )

            content_parts.append(PLOT_CONTAINER_TEMPLATE.format(
                extra_class="full-width",
                title="Posterior Entropy Distribution",
                description=(
                    "Shannon entropy of the posterior distribution per read. "
                    "Low entropy (left) indicates confident classification; "
                    "high entropy (right) indicates reads near boundaries "
                    "where multiple categories are plausible."
                ),
                plot_id=entropy_id,
            ))

        # --- Stacked Posterior Bar by Category ---
        if s.has_bayesian and "taxonomic_call" in self.df.columns:
            posterior_cols = [
                ("p_known_species", "Known Species"),
                ("p_novel_species", "Novel Species"),
                ("p_novel_genus", "Novel Genus"),
                ("p_species_boundary", "Species Boundary"),
                ("p_ambiguous", "Ambiguous"),
                ("p_unclassified", "Unclassified"),
            ]

            bar_fig = go.Figure()
            all_map_cats = [
                "Known Species", "Novel Species", "Novel Genus",
                "Species Boundary", "Ambiguous", "Unclassified",
            ]
            for col_name, cat_label in posterior_cols:
                if col_name not in self.df.columns:
                    continue
                means = []
                cats = []
                for map_cat in all_map_cats:
                    cat_df = self.df.filter(pl.col("taxonomic_call") == map_cat)
                    if len(cat_df) > 0:
                        cats.append(map_cat)
                        means.append(cat_df[col_name].mean() or 0.0)
                if means:
                    bar_fig.add_trace(go.Bar(
                        x=cats,
                        y=means,
                        name=cat_label,
                        marker_color=BAYESIAN_CATEGORY_COLORS.get(cat_label, "#94a3b8"),
                    ))

            bar_fig.update_layout(
                title="Mean Posterior Composition by MAP Category",
                xaxis_title="Bayesian MAP Category",
                yaxis_title="Mean Posterior Probability",
                template="plotly_white",
                height=450,
                barmode="stack",
                legend={"orientation": "h", "y": -0.2},
                yaxis={"range": [0, 1.05]},
            )
            bar_id = "plot-bayesian-posterior-bar"
            self._register_plot(bar_id, bar_fig)

            content_parts.append(PLOT_CONTAINER_TEMPLATE.format(
                extra_class="full-width",
                title="Posterior Composition by Category",
                description=(
                    "Mean posterior probabilities grouped by the Bayesian MAP classification. "
                    "A well-separated classifier shows tall bars for the dominant category "
                    "in each group. Mixed bars indicate regions of classification uncertainty."
                ),
                plot_id=bar_id,
            ))

        # --- Discovery Score Distribution ---
        if s.has_enhanced_scoring and "discovery_score" in self.df.columns:
            discovery_df = self.df.filter(pl.col("discovery_score").is_not_null())
            if len(discovery_df) > 0:
                discovery_hist_id = "plot-discovery-hist"
                self._build_histogram(
                    data=discovery_df["discovery_score"].to_list(),
                    plot_id=discovery_hist_id,
                    title="Discovery Score Distribution (Novel Reads)",
                    x_label="Discovery Score",
                    nbins=25,
                    thresholds=[
                        (75, "dash", "#22c55e", "High Priority (75+)"),
                        (50, "dot", "#f59e0b", "Moderate (50+)"),
                        (25, "dot", "#ef4444", "Low (25+)"),
                    ],
                )

                content_parts.append(PLOT_CONTAINER_TEMPLATE.format(
                    extra_class="full-width",
                    title="Discovery Score Distribution",
                    description=(
                        "Distribution of discovery scores for novel reads. Higher scores indicate "
                        "more reliable discoveries. Green line marks high-priority threshold (75+)."
                    ),
                    plot_id=discovery_hist_id,
                ))

        # --- Uncertainty types breakdown ---
        if s.has_inferred_uncertainty and "uncertainty_type" in self.df.columns:
            measured_df = self.df.filter(pl.col("uncertainty_type") == "measured")
            inferred_df = self.df.filter(pl.col("uncertainty_type") == "inferred")
            measured_count = len(measured_df)
            inferred_count = len(inferred_df)
            total = s.total_reads if s.total_reads > 0 else 1

            uncertainty_types_html = ENHANCED_SCORING_UNCERTAINTY_TYPES_TEMPLATE.format(
                measured_count=measured_count,
                measured_pct=measured_count / total * 100,
                inferred_count=inferred_count,
                inferred_pct=inferred_count / total * 100,
                single_hit_pct=s.single_hit_pct,
            )
            content_parts.append(uncertainty_types_html)

        # --- Confidence Landscape Scatter ---
        if s.has_bayesian and "posterior_entropy" in self.df.columns and "novelty_index" in self.df.columns:
            scatter_id = "plot-bayesian-landscape"
            self._build_category_scatter(
                x_col="novelty_index",
                y_col="posterior_entropy",
                plot_id=scatter_id,
                title="Classification Confidence Landscape",
                x_label="Novelty Index (%)",
                y_label="Posterior Entropy (bits)",
                y_range=[0, 2.7],
            )

            content_parts.append(PLOT_CONTAINER_TEMPLATE.format(
                extra_class="full-width",
                title="Confidence Landscape",
                description=(
                    "Novelty index vs posterior entropy colored by Bayesian MAP category. "
                    "Reads with high entropy (top) lie near classification boundaries. "
                    "Vertical bands of high entropy correspond to threshold boundaries."
                ),
                plot_id=scatter_id,
            ))

        return "\n".join(content_parts)

    def _build_novel_diversity_section(self) -> str:
        """Build the novel diversity tab section with cluster analysis."""
        from metadarkmatter.core.novel_diversity import NovelDiversityAnalyzer

        content_parts = []

        # Create analyzer and cluster novel reads
        try:
            from metadarkmatter.core.classification.ani_matrix import ANIMatrix

            ani_matrix_obj = None
            if self.ani_matrix is not None and len(self.ani_matrix) > 0:
                try:
                    ani_matrix_obj = ANIMatrix.from_dataframe(self.ani_matrix)
                except Exception as e:
                    logger.warning(f"Could not convert ANI matrix for clustering: {e}")

            analyzer = NovelDiversityAnalyzer(
                classifications=self.df,
                metadata=self.genome_metadata,
                ani_matrix=ani_matrix_obj,
                novelty_band_size=5.0,
                min_cluster_size=3,
            )
            clusters = analyzer.cluster_novel_reads()
            summary = analyzer.get_summary()
        except Exception as e:
            logger.warning(f"Could not analyze novel diversity: {e}")
            content = EMPTY_SECTION_TEMPLATE.format(
                message="Could not analyze novel diversity. Check classification data."
            )
            return TAB_SECTION_TEMPLATE.format(
                tab_id="novel-diversity",
                active_class="",
                section_title="Novel Diversity",
                content=content,
            )

        if not clusters:
            content = NOVEL_EMPTY_TEMPLATE
            return TAB_SECTION_TEMPLATE.format(
                tab_id="novel-diversity",
                active_class="",
                section_title="Novel Diversity",
                content=content,
            )

        # Build KPI metric strip for novel diversity
        novel_kpi_cards = []
        cluster_count = len(clusters) if clusters else 0
        novel_kpi_cards.append(KPI_CARD_TEMPLATE.format(
            color_class="accent", value=str(cluster_count), label="Clusters",
        ))
        novel_kpi_cards.append(KPI_CARD_TEMPLATE.format(
            color_class="novel", value=str(self.summary.novel_species), label="Novel Species",
        ))
        novel_kpi_cards.append(KPI_CARD_TEMPLATE.format(
            color_class="novel", value=str(self.summary.novel_genus), label="Novel Genus",
        ))
        if self.summary.has_bayesian:
            high_conf_novel = len(self.df.filter(
                (pl.col("taxonomic_call").is_in(["Novel Species", "Novel Genus"])) &
                (pl.col("posterior_entropy") < 1.0)
            ))
            novel_kpi_cards.append(KPI_CARD_TEMPLATE.format(
                color_class="known", value=str(high_conf_novel), label="High Confidence",
            ))
        kpi_html = KPI_STRIP_TEMPLATE.format(cards="\n".join(novel_kpi_cards))
        content_parts.append(kpi_html)

        # Build summary section
        content_parts.append(build_novel_summary_html(summary))

        # Build scatter plot of cluster quality
        scatter_fig = build_cluster_scatter_figure(clusters, width=1000, height=500)
        scatter_id = "plot-novel-scatter"
        self._register_plot(scatter_id, scatter_fig)

        content_parts.append(PLOT_CONTAINER_TEMPLATE.format(
            extra_class="full-width",
            title="Cluster Quality Overview",
            description=(
                "Scatter plot showing cluster novelty vs quality metrics. "
                "Marker size indicates read count. Colors indicate confidence rating."
            ),
            plot_id=scatter_id,
        ))

        # Build sunburst chart for phylogenetic context
        sunburst_fig = build_sunburst_figure(clusters, width=600, height=600)
        sunburst_id = "plot-novel-sunburst"
        self._register_plot(sunburst_id, sunburst_fig)

        content_parts.append(PLOT_CONTAINER_TEMPLATE.format(
            extra_class="full-width",
            title="Phylogenetic Context",
            description=(
                "Sunburst chart showing novel clusters within taxonomic hierarchy. "
                "Click to explore family/genus/cluster relationships."
            ),
            plot_id=sunburst_id,
        ))

        # Build phylogenetic context heatmap(s) based on available matrices
        # Prefer AAI for protein mode, ANI for nucleotide mode
        # Show both when both are available (ANI for species, AAI for genus context)
        is_protein_mode = self.config.alignment_mode == "protein"

        # Determine which matrices to use
        heatmaps_to_build: list[tuple[pl.DataFrame, str]] = []

        if is_protein_mode:
            # Protein mode: prefer AAI, fallback to ANI
            if self.aai_matrix is not None and len(self.aai_matrix) > 0:
                heatmaps_to_build.append((self.aai_matrix, "AAI"))
            elif self.ani_matrix is not None and len(self.ani_matrix) > 0:
                heatmaps_to_build.append((self.ani_matrix, "ANI"))
        else:
            # Nucleotide mode: use ANI primarily
            if self.ani_matrix is not None and len(self.ani_matrix) > 0:
                heatmaps_to_build.append((self.ani_matrix, "ANI"))
            # Also show AAI if available (better for Novel Genus clusters)
            if self.aai_matrix is not None and len(self.aai_matrix) > 0:
                heatmaps_to_build.append((self.aai_matrix, "AAI"))

        heatmap_parts: list[str] = []
        for matrix, sim_type in heatmaps_to_build:
            try:
                # Get genome accessions from matrix
                matrix_cols = matrix.columns
                if "genome" in matrix_cols:
                    genome_accessions = matrix["genome"].to_list()
                else:
                    genome_accessions = list(matrix_cols)

                genome_labels_map = self._get_genome_labels_map(genome_accessions)

                # Build the phylogenetic context heatmap
                phylo_fig, phylo_metadata, clustering_ok = build_phylogenetic_context_heatmap(
                    similarity_matrix=matrix,
                    novel_clusters=clusters,
                    genome_labels_map=genome_labels_map,
                    similarity_type=sim_type,
                    max_references=self.config.max_phylo_references,
                    max_clusters=self.config.max_phylo_clusters,
                )

                if phylo_metadata.get("n_total", 0) > 0:
                    # Add intro section with note if clusters were truncated
                    n_clusters_shown = phylo_metadata.get("n_clusters", 0)
                    total_clusters = len(clusters)
                    clusters_note = ""
                    if n_clusters_shown < total_clusters:
                        clusters_note = f" (top {n_clusters_shown} of {total_clusters} by read count)"

                    # Customize description based on similarity type
                    if sim_type == "AAI":
                        intro_description = (
                            "This AAI-based heatmap is particularly useful for visualizing "
                            "Novel Genus candidates, where protein-level similarity provides "
                            "better resolution than nucleotide comparisons."
                        )
                    else:
                        intro_description = (
                            "This ANI-based heatmap shows novel clusters positioned alongside "
                            "reference genomes. Best for visualizing Novel Species placement."
                        )

                    heatmap_parts.append(PHYLOGENETIC_HEATMAP_INTRO_TEMPLATE.format(
                        n_references=phylo_metadata.get("n_references", 0),
                        n_clusters=n_clusters_shown,
                        clusters_note=clusters_note,
                    ))

                    # Register and add the heatmap
                    phylo_heatmap_id = f"plot-novel-phylo-heatmap-{sim_type.lower()}"
                    self._register_plot(phylo_heatmap_id, phylo_fig)

                    heatmap_parts.append(PLOT_CONTAINER_TEMPLATE.format(
                        extra_class="full-width",
                        title=f"Phylogenetic Context Heatmap ({sim_type})",
                        description=(
                            f"Extended {sim_type} heatmap showing novel clusters alongside reference genomes. "
                            f"Entries marked with [*] are novel clusters. Hover for details."
                        ),
                        plot_id=phylo_heatmap_id,
                    ))

                    # Add legend (only for first heatmap to avoid repetition)
                    if matrix is heatmaps_to_build[0][0]:
                        heatmap_parts.append(PHYLOGENETIC_HEATMAP_LEGEND_TEMPLATE)

            except Exception as e:
                logger.warning(f"Could not build {sim_type} phylogenetic context heatmap: {e}")

        # Wrap heatmaps in a collapsible panel if any were built
        if heatmap_parts:
            heatmap_html = "\n".join(heatmap_parts)
            content_parts.append(COLLAPSIBLE_PANEL_TEMPLATE.format(
                default_state="",
                panel_id="novel-heatmap",
                title="Phylogenetic Context Heatmap",
                content=heatmap_html,
            ))

        # Build cluster table (includes phylogenetic placement)
        content_parts.append(build_cluster_table_html(clusters))

        content = "\n".join(content_parts)

        return TAB_SECTION_TEMPLATE.format(
            tab_id="novel-diversity",
            active_class="",
            section_title="Novel Diversity",
            content=content,
        )

    def _build_recruitment_content(self) -> str | None:
        """Build recruitment plot HTML content without tab wrapper.

        Returns:
            HTML string with recruitment plots, or None if no recruitment
            data is available.
        """
        if self.recruitment_data is None or len(self.recruitment_data) == 0:
            return None

        # Import here to avoid circular imports
        from metadarkmatter.visualization.recruitment_plots import (
            RecruitmentPlotGenerator,
        )

        # Create recruitment plot
        gen = RecruitmentPlotGenerator(self.recruitment_data)
        fig = gen.create_figure(
            max_points=self.config.max_scatter_points,
            width=1100,
            height=600,
        )

        recruit_id = "plot-recruitment"
        self._register_plot(recruit_id, fig)

        # Multi-genome view
        multi_fig = gen.create_multi_genome_figure(
            max_points_per_genome=20000,
            width=1100,
        )
        multi_id = "plot-recruitment-multi"
        self._register_plot(multi_id, multi_fig)

        return (
            PLOT_CONTAINER_TEMPLATE.format(
                extra_class="full-width",
                title="Read Recruitment Overview",
                description=(
                    "Scatter plot showing percent identity vs genome position "
                    "for all mapped reads."
                ),
                plot_id=recruit_id,
            )
            + PLOT_CONTAINER_TEMPLATE.format(
                extra_class="full-width",
                title="Per-Genome Recruitment",
                description=(
                    "Individual recruitment plots for top reference genomes "
                    "by read count."
                ),
                plot_id=multi_id,
            )
        )

    def _build_species_bar_chart(self) -> str:
        """Build a species bar chart for the Reference tab.

        Extracts only the horizontal bar chart of top species by read count.
        Requires genome metadata to be available.

        Returns:
            HTML string with the species bar chart, or empty string if no
            metadata is available.
        """
        if self.genome_metadata is None:
            return ""

        import plotly.graph_objects as go

        # Join metadata to get species information
        if "species" not in self.df.columns:
            enriched_df = self.genome_metadata.join_classifications(self.df)
        else:
            enriched_df = self.df

        # Only include in-family reads in species breakdown
        if "taxonomic_call" in enriched_df.columns:
            species_base_df = enriched_df.filter(
                pl.col("taxonomic_call") != "Off-target"
            )
        else:
            species_base_df = enriched_df

        # Aggregate by species
        species_counts = (
            species_base_df.group_by("species")
            .agg([
                pl.len().alias("count"),
                pl.col("novelty_index").mean().alias("mean_novelty"),
                pl.col("top_hit_identity").mean().alias("mean_identity"),
                pl.col("best_match_genome").n_unique().alias("genome_count"),
            ])
            .sort("count", descending=True)
            .head(20)
        )

        if len(species_counts) == 0:
            return ""

        # Horizontal bar chart of species read counts
        bar_fig = go.Figure()
        bar_fig.add_trace(
            go.Bar(
                y=species_counts["species"].to_list()[::-1],
                x=species_counts["count"].to_list()[::-1],
                orientation="h",
                marker_color="#22c55e",
                text=[f"{c:,}" for c in species_counts["count"].to_list()[::-1]],
                textposition="outside",
                hovertemplate=(
                    "<b>%{y}</b><br>"
                    "Reads: %{x:,}<extra></extra>"
                ),
            )
        )
        bar_fig.update_layout(
            title="Top 20 Species by Read Count",
            xaxis_title="Number of Reads",
            yaxis_title="",
            template="plotly_white",
            height=max(400, 30 * len(species_counts)),
            margin={"l": 300, "r": 40, "t": 60, "b": 60},
        )
        species_bar_id = "plot-species-bar"
        self._register_plot(species_bar_id, bar_fig)

        return PLOT_CONTAINER_SIMPLE_TEMPLATE.format(
            extra_class="full-width",
            plot_id=species_bar_id,
        )

    def _build_family_validation_section(self) -> str:
        """Build the Family Validation tab content."""
        import plotly.graph_objects as go

        total = self.summary.total_reads
        off_target = self.summary.off_target
        in_family = total - off_target
        off_target_pct = (off_target / total * 100) if total > 0 else 0.0
        validated_pct = 100.0 - off_target_pct

        # Build external genomes breakdown table from off-target reads
        external_families_table = ""
        if (
            "external_best_genome" in self.df.columns
            and "taxonomic_call" in self.df.columns
        ):
            off_target_df = self.df.filter(
                pl.col("taxonomic_call") == "Off-target"
            )
            if len(off_target_df) > 0:
                ext_col = "external_best_genome"
                ident_col = "external_best_identity"

                # Group by external genome, count reads and mean identity
                agg_cols = [pl.len().alias("read_count")]
                has_identity = ident_col in off_target_df.columns
                if has_identity:
                    agg_cols.append(
                        pl.col(ident_col).mean().alias("mean_identity")
                    )
                has_if_ident = "in_family_identity" in off_target_df.columns
                if has_if_ident:
                    agg_cols.append(
                        pl.col("in_family_identity").mean().alias("mean_in_family_identity")
                    )

                ext_summary = (
                    off_target_df
                    .with_columns(
                        pl.col(ext_col).fill_null("unknown").alias(ext_col)
                    )
                    .group_by(ext_col)
                    .agg(agg_cols)
                    .sort("read_count", descending=True)
                    .head(20)
                )

                if len(ext_summary) > 0:
                    if_ident_header = (
                        "<th>In-Family Identity</th>" if has_if_ident else ""
                    )
                    rows_html = ""
                    for row in ext_summary.iter_rows(named=True):
                        genome = row[ext_col]
                        count = row["read_count"]
                        pct = count / off_target * 100 if off_target > 0 else 0
                        ident_str = (
                            f"{row['mean_identity']:.1f}%"
                            if has_identity and row.get("mean_identity")
                            is not None
                            else "-"
                        )
                        if_ident_cell = ""
                        if has_if_ident:
                            if_val = row.get("mean_in_family_identity")
                            if_ident_cell = (
                                f"<td>{if_val:.1f}%</td>"
                                if if_val is not None
                                else "<td>-</td>"
                            )
                        rows_html += (
                            f"<tr><td>{genome}</td>"
                            f"<td>{count:,}</td>"
                            f"<td>{pct:.1f}%</td>"
                            f"<td>{ident_str}</td>"
                            f"{if_ident_cell}</tr>\n"
                        )

                    external_families_table = (
                        '<h3>Off-target Best External Matches</h3>\n'
                        '<div style="overflow-x:auto;">\n'
                        '<table class="data-table">\n'
                        "<thead><tr>"
                        "<th>External Genome</th>"
                        "<th>Reads</th>"
                        "<th>% of Off-target</th>"
                        "<th>Mean Identity</th>"
                        f"{if_ident_header}"
                        "</tr></thead>\n<tbody>\n"
                        + rows_html
                        + "</tbody></table></div>\n"
                    )

        # Borderline analysis: classify off-target reads by in-family novelty
        borderline_analysis = ""
        if (
            "in_family_novelty_index" in self.df.columns
            and "taxonomic_call" in self.df.columns
        ):
            ot_with_novelty = self.df.filter(
                pl.col("taxonomic_call") == "Off-target"
            )
            if len(ot_with_novelty) > 0:
                th = self.config.thresholds
                n_total_ot = len(ot_with_novelty)
                ni = ot_with_novelty["in_family_novelty_index"]

                known_count = len(ni.filter(ni < th.novelty_known_max))
                novel_sp_count = len(ni.filter(
                    (ni >= th.novelty_novel_species_min)
                    & (ni < th.novelty_novel_species_max)
                ))
                novel_gen_count = len(ni.filter(
                    (ni >= th.novelty_novel_genus_min)
                    & (ni <= th.novelty_novel_genus_max)
                ))
                divergent_count = len(ni.filter(ni > th.novelty_novel_genus_max))

                def _pct(c: int) -> float:
                    return c / n_total_ot * 100 if n_total_ot > 0 else 0.0

                borderline_analysis = BORDERLINE_ANALYSIS_TEMPLATE.format(
                    known_count=known_count,
                    known_pct=_pct(known_count),
                    novel_sp_count=novel_sp_count,
                    novel_sp_pct=_pct(novel_sp_count),
                    novel_gen_count=novel_gen_count,
                    novel_gen_pct=_pct(novel_gen_count),
                    divergent_count=divergent_count,
                    divergent_pct=_pct(divergent_count),
                    novel_sp_highlight=(
                        " highlight-amber" if novel_sp_count > 0 else ""
                    ),
                    novel_gen_highlight=(
                        " highlight-red" if novel_gen_count > 0 else ""
                    ),
                    novelty_known_max=th.novelty_known_max,
                    novelty_novel_sp_max=th.novelty_novel_species_max,
                    novelty_novel_gen_max=th.novelty_novel_genus_max,
                )

                # Scatter plot: off-target reads in novelty-uncertainty space
                novelty_vals = ot_with_novelty["in_family_novelty_index"].to_list()

                # Compute inferred uncertainty using the same piecewise formula
                # as VectorizedClassifier (lines 773-789)
                known_break = 5.0 + th.novelty_known_max * 0.5
                species_break = (
                    known_break
                    + (th.novelty_novel_species_max - th.novelty_known_max)
                )
                uncertainty_vals = []
                for nv in novelty_vals:
                    if nv < th.novelty_known_max:
                        uncertainty_vals.append(5.0 + nv * 0.5)
                    elif nv < th.novelty_novel_species_max:
                        uncertainty_vals.append(
                            known_break
                            + (nv - th.novelty_known_max) * 1.0
                        )
                    elif nv < th.novelty_novel_genus_max:
                        uncertainty_vals.append(
                            species_break
                            + (nv - th.novelty_novel_species_max) * 1.5
                        )
                    else:
                        uncertainty_vals.append(35.0)

                # Assign would-be categories for coloring
                cat_colors = []
                cat_labels = []
                for nv in novelty_vals:
                    if nv < th.novelty_known_max:
                        cat_colors.append("#22c55e")
                        cat_labels.append("Known Species")
                    elif nv < th.novelty_novel_species_max:
                        cat_colors.append("#f59e0b")
                        cat_labels.append("Novel Species")
                    elif nv <= th.novelty_novel_genus_max:
                        cat_colors.append("#ef4444")
                        cat_labels.append("Novel Genus")
                    else:
                        cat_colors.append("#64748b")
                        cat_labels.append("Very Divergent")

                scatter_fig = go.Figure()

                # Background classification zones
                zone_shapes = [
                    dict(
                        type="rect", x0=0, x1=th.novelty_known_max,
                        y0=0, y1=th.uncertainty_known_max,
                        fillcolor="rgba(34,197,94,0.08)",
                        line=dict(width=0), layer="below",
                    ),
                    dict(
                        type="rect",
                        x0=th.novelty_novel_species_min,
                        x1=th.novelty_novel_species_max,
                        y0=0, y1=th.uncertainty_novel_species_max,
                        fillcolor="rgba(245,158,11,0.08)",
                        line=dict(width=0), layer="below",
                    ),
                    dict(
                        type="rect",
                        x0=th.novelty_novel_genus_min,
                        x1=th.novelty_novel_genus_max,
                        y0=0, y1=th.uncertainty_novel_genus_max,
                        fillcolor="rgba(239,68,68,0.08)",
                        line=dict(width=0), layer="below",
                    ),
                ]

                scatter_fig.add_trace(go.Scatter(
                    x=novelty_vals,
                    y=uncertainty_vals,
                    mode="markers",
                    marker=dict(
                        color=cat_colors,
                        size=6,
                        opacity=0.7,
                    ),
                    text=cat_labels,
                    hovertemplate=(
                        "In-Family Novelty: %{x:.1f}%<br>"
                        "Inferred Uncertainty: %{y:.1f}<br>"
                        "Would-Be: %{text}<extra></extra>"
                    ),
                ))
                scatter_fig.update_layout(
                    xaxis_title="In-Family Novelty Index (%)",
                    yaxis_title="Inferred Uncertainty",
                    template="plotly_white",
                    height=450,
                    shapes=zone_shapes,
                    showlegend=False,
                )
                self._register_plot("family-novelty-scatter", scatter_fig)

        # Format summary stats
        content = FAMILY_VALIDATION_SECTION_TEMPLATE.format(
            validated_pct=validated_pct,
            off_target_count=off_target,
            off_target_pct=off_target_pct,
            target_family=(
                self.summary.target_family or "Inferred from ANI matrix"
            ),
            borderline_analysis=borderline_analysis,
            external_families_table=external_families_table,
        )

        # Histogram of family_bitscore_ratio -- exclude ratio=1.0 reads
        # to avoid a dominant spike that obscures the off-target distribution
        if "family_bitscore_ratio" in self.df.columns:
            all_ratios = self.df["family_bitscore_ratio"].drop_nulls()
            n_total = len(all_ratios)
            n_at_one = len(all_ratios.filter(all_ratios == 1.0))
            sub_ratios = all_ratios.filter(all_ratios < 1.0).to_list()

            if sub_ratios:
                fig = go.Figure(data=[go.Histogram(
                    x=sub_ratios,
                    nbinsx=50,
                    marker_color="#4a90d9",
                    name="Reads",
                )])
                fig.update_layout(
                    xaxis_title="Family Bitscore Ratio",
                    yaxis_title="Read Count",
                    template="plotly_white",
                    height=400,
                )
                fig.add_vline(
                    x=0.8,
                    line_dash="dash",
                    line_color="red",
                    annotation_text="Threshold (0.8)",
                )
                if n_at_one > 0:
                    fig.add_annotation(
                        text=(
                            f"{n_at_one:,} reads at ratio=1.0 "
                            f"({n_at_one / n_total * 100:.0f}% "
                            "of total, not shown)"
                        ),
                        xref="paper", yref="paper",
                        x=0.98, y=0.98,
                        showarrow=False,
                        font=dict(size=11, color="#666"),
                        align="right",
                        bgcolor="rgba(255,255,255,0.8)",
                        bordercolor="#ccc",
                        borderwidth=1,
                    )
                self._register_plot("family-ratio-histogram", fig)
            elif n_at_one > 0:
                # All reads at 1.0 -- show a simple note instead of empty plot
                fig = go.Figure()
                fig.update_layout(
                    template="plotly_white",
                    height=200,
                    annotations=[dict(
                        text=(
                            f"All {n_at_one:,} reads have "
                            "family bitscore ratio = 1.0 "
                            "(entirely in-family)"
                        ),
                        xref="paper", yref="paper",
                        x=0.5, y=0.5,
                        showarrow=False,
                        font=dict(size=14),
                    )],
                    xaxis=dict(visible=False),
                    yaxis=dict(visible=False),
                )
                self._register_plot("family-ratio-histogram", fig)

        # Pie chart of in-family vs off-target
        fig_pie = go.Figure(data=[go.Pie(
            labels=["In-Family", "Off-target"],
            values=[in_family, off_target],
            marker_colors=["#4a90d9", "#e74c3c"],
            hole=0.4,
        )])
        fig_pie.update_layout(
            template="plotly_white",
            height=400,
        )
        self._register_plot("family-pie-chart", fig_pie)

        return TAB_SECTION_TEMPLATE.format(
            tab_id="family-validation",
            active_class="",
            section_title="Family Validation",
            content=content,
        )

    def _build_toggle_heatmap_content(
        self,
        *,
        metric: str,
        title: str,
        description: str,
        all_plot_id: str,
        all_stats_html: str,
        rep_plot_id: str,
        rep_stats_html: str,
        n_all: int,
        n_reps: int,
    ) -> str:
        """Build HTML content with toggle between all genomes and representatives.

        Args:
            metric: Short identifier ('ani' or 'aai') used for element IDs
            title: Plot title text
            description: Plot description text
            all_plot_id: Plotly div ID for the all-genomes heatmap
            all_stats_html: Statistics cards HTML for all genomes
            rep_plot_id: Plotly div ID for the representative-only heatmap
            rep_stats_html: Statistics cards HTML for representatives
            n_all: Total number of genomes
            n_reps: Number of representative genomes

        Returns:
            HTML string with toggle controls and both heatmaps
        """
        return f"""
        <div class="heatmap-toggle-controls">
            <button class="heatmap-toggle-btn active"
                    id="{metric}-toggle-all"
                    onclick="toggleHeatmapView('{metric}', 'all')">
                All Genomes ({n_all})
            </button>
            <button class="heatmap-toggle-btn"
                    id="{metric}-toggle-reps"
                    onclick="toggleHeatmapView('{metric}', 'reps')">
                Representatives Only ({n_reps})
            </button>
        </div>
        <div id="{metric}-view-all">
            {all_stats_html}
            <div class="plot-container full-width">
                <div class="plot-title">{title}</div>
                <div class="plot-description">{description}</div>
                <div id="{all_plot_id}" class="plotly-chart"></div>
            </div>
        </div>
        <div id="{metric}-view-reps" style="display:none;">
            {rep_stats_html}
            <div class="plot-container full-width">
                <div class="plot-title">{title} (Representatives Only)</div>
                <div class="plot-description">
                    Showing {n_reps} species representative genomes
                    (one per species from GTDB). {description}
                </div>
                <div id="{rep_plot_id}" class="plotly-chart"></div>
            </div>
        </div>
        """

    def _filter_matrix_to_representatives(
        self, matrix: pl.DataFrame
    ) -> pl.DataFrame | None:
        """Filter a similarity matrix to representative genomes only.

        Args:
            matrix: Similarity matrix DataFrame (wide format, optional 'genome' column)

        Returns:
            Filtered matrix with only representative genomes, or None if
            metadata lacks representative information or filtering would
            remove all genomes.
        """
        if self.genome_metadata is None or not self.genome_metadata.has_representatives:
            return None

        rep_set = set(self.genome_metadata.build_representative_mapping().values())

        # Identify accessions in the matrix
        cols = matrix.columns
        has_genome_col = "genome" in cols
        if has_genome_col:
            accessions = matrix["genome"].to_list()
        else:
            accessions = list(cols)

        # Filter to representatives present in the matrix
        keep = [acc for acc in accessions if acc in rep_set]
        if len(keep) < 2:
            return None  # Need at least 2 genomes for a meaningful heatmap

        if has_genome_col:
            filtered = matrix.filter(pl.col("genome").is_in(keep))
            filtered = filtered.select(["genome"] + keep)
        else:
            # Columns are accessions; filter both rows and columns
            indices = [i for i, acc in enumerate(accessions) if acc in keep]
            filtered = matrix.select([accessions[i] for i in indices])
            filtered = filtered[indices]

        return filtered

    def _build_ani_content(self) -> str | None:
        """Build ANI heatmap HTML content without tab wrapper.

        Returns:
            HTML string with ANI heatmap and statistics cards, or None if
            no ANI matrix is available.
        """
        if self.ani_matrix is None or len(self.ani_matrix) == 0:
            return None

        # Create genome labels map
        matrix_cols = self.ani_matrix.columns
        if "genome" in matrix_cols:
            genome_accessions = self.ani_matrix["genome"].to_list()
        else:
            genome_accessions = list(matrix_cols)

        genome_labels_map = self._get_genome_labels_map(genome_accessions)

        # Build full heatmap
        heatmap_fig, stats, clustering_succeeded = build_ani_heatmap(
            self.ani_matrix,
            genome_labels_map,
        )

        ani_id = "plot-ani-heatmap"
        self._register_plot(ani_id, heatmap_fig)
        stats_html = build_ani_stats_cards(stats)

        # Check for representative genome toggle
        rep_matrix = self._filter_matrix_to_representatives(self.ani_matrix)

        if rep_matrix is not None:
            # Build representative-only heatmap
            rep_cols = rep_matrix.columns
            if "genome" in rep_cols:
                rep_accessions = rep_matrix["genome"].to_list()
            else:
                rep_accessions = list(rep_cols)
            rep_labels = self._get_genome_labels_map(rep_accessions)
            rep_fig, rep_stats, _ = build_ani_heatmap(rep_matrix, rep_labels)
            rep_id = "plot-ani-heatmap-reps"
            self._register_plot(rep_id, rep_fig)
            rep_stats_html = build_ani_stats_cards(rep_stats)

            n_all = len(genome_accessions)
            n_reps = len(rep_accessions)

            return self._build_toggle_heatmap_content(
                metric="ani",
                title="Average Nucleotide Identity Matrix",
                description=(
                    "Heatmap showing pairwise ANI values between reference genomes, "
                    "hierarchically clustered by similarity. "
                    "Red = high ANI (similar), blue = low ANI (distant). "
                    "ANI ~95% indicates species boundary (Jain et al. 2018)."
                ),
                all_plot_id=ani_id,
                all_stats_html=stats_html,
                rep_plot_id=rep_id,
                rep_stats_html=rep_stats_html,
                n_all=n_all,
                n_reps=n_reps,
            )

        return stats_html + PLOT_CONTAINER_TEMPLATE.format(
            extra_class="full-width",
            title="Average Nucleotide Identity Matrix",
            description=(
                "Heatmap showing pairwise ANI values between reference genomes, "
                "hierarchically clustered by similarity. "
                "Red = high ANI (similar), blue = low ANI (distant). "
                "ANI ~95% indicates species boundary (Jain et al. 2018)."
            ),
            plot_id=ani_id,
        )

    def _build_aai_content(self) -> str | None:
        """Build AAI heatmap HTML content without tab wrapper.

        Returns:
            HTML string with AAI heatmap and statistics cards, or None if
            no AAI matrix is available.
        """
        if self.aai_matrix is None or len(self.aai_matrix) == 0:
            return None

        # Create genome labels map
        matrix_cols = self.aai_matrix.columns
        if "genome" in matrix_cols:
            genome_accessions = self.aai_matrix["genome"].to_list()
        else:
            genome_accessions = list(matrix_cols)

        genome_labels_map = self._get_genome_labels_map(genome_accessions)

        # Build full heatmap
        heatmap_fig, stats, clustering_succeeded = build_aai_heatmap(
            self.aai_matrix,
            genome_labels_map,
        )

        aai_id = "plot-aai-heatmap"
        self._register_plot(aai_id, heatmap_fig)
        stats_html = build_aai_stats_cards(stats)

        # Check for representative genome toggle
        rep_matrix = self._filter_matrix_to_representatives(self.aai_matrix)

        if rep_matrix is not None:
            rep_cols = rep_matrix.columns
            if "genome" in rep_cols:
                rep_accessions = rep_matrix["genome"].to_list()
            else:
                rep_accessions = list(rep_cols)
            rep_labels = self._get_genome_labels_map(rep_accessions)
            rep_fig, rep_stats, _ = build_aai_heatmap(rep_matrix, rep_labels)
            rep_id = "plot-aai-heatmap-reps"
            self._register_plot(rep_id, rep_fig)
            rep_stats_html = build_aai_stats_cards(rep_stats)

            n_all = len(genome_accessions)
            n_reps = len(rep_accessions)

            return self._build_toggle_heatmap_content(
                metric="aai",
                title="Average Amino Acid Identity Matrix",
                description=(
                    "Heatmap showing pairwise AAI values between reference genomes, "
                    "hierarchically clustered by similarity. "
                    "Red = high AAI (similar), blue = low AAI (distant). "
                    "AAI ~65% indicates genus boundary (Riesco & Trujillo 2024)."
                ),
                all_plot_id=aai_id,
                all_stats_html=stats_html,
                rep_plot_id=rep_id,
                rep_stats_html=rep_stats_html,
                n_all=n_all,
                n_reps=n_reps,
            )

        return stats_html + PLOT_CONTAINER_TEMPLATE.format(
            extra_class="full-width",
            title="Average Amino Acid Identity Matrix",
            description=(
                "Heatmap showing pairwise AAI values between reference genomes, "
                "hierarchically clustered by similarity. "
                "Red = high AAI (similar), blue = low AAI (distant). "
                "AAI ~65% indicates genus boundary (Riesco & Trujillo 2024)."
            ),
            plot_id=aai_id,
        )

    def _build_phylogeny_content(self) -> str | None:
        """Build phylogeny HTML content without tab wrapper.

        Generates an interactive phylogenetic tree from the ANI matrix (or
        a user-provided tree) with novel cluster placement.

        Returns:
            HTML string with phylogeny visualization, or None if the tree
            cannot be built (e.g., skip requested, no ANI matrix, or fewer
            than 3 genomes).
        """
        if self.config.skip_phylogeny:
            return None
        if self.ani_matrix is None or len(self.ani_matrix) < 3:
            return None

        ani_pd = self.ani_matrix.to_pandas()
        if "genome" in ani_pd.columns:
            ani_pd = ani_pd.set_index("genome")

        return self._build_phylogeny_section(
            ani_pd,
            user_tree_path=self.config.user_tree_path,
        )

    def _build_reference_section(self) -> str:
        """Build the Reference tab with species bar, heatmaps, phylogeny, and recruitment.

        Merges content from the former Species & Genomes, ANI, AAI,
        Phylogeny, and Recruitment tabs into a single unified section.

        Returns:
            HTML string for the Reference tab section.
        """
        content_parts: list[str] = []

        # 1. Species bar chart (compact, full-width)
        species_bar = self._build_species_bar_chart()
        if species_bar:
            content_parts.append(species_bar)

        # 2. ANI/AAI heatmaps
        ani_content = self._build_ani_content()
        aai_content = self._build_aai_content()

        if ani_content and aai_content:
            # Side-by-side layout
            heatmap_row = TWO_COLUMN_ROW_TEMPLATE.format(
                left_flex="1",
                right_flex="1",
                left_content=ani_content,
                right_content=aai_content,
            )
            content_parts.append(heatmap_row)
        elif ani_content:
            content_parts.append(ani_content)
        elif aai_content:
            content_parts.append(aai_content)

        # 3. Phylogenetic tree (collapsible, default collapsed)
        phylo_html = self._build_phylogeny_content()
        if phylo_html:
            content_parts.append(COLLAPSIBLE_PANEL_TEMPLATE.format(
                default_state="",
                panel_id="phylogeny",
                title="Phylogenetic Tree",
                content=phylo_html,
            ))

        # 4. Recruitment plots (collapsible, default collapsed, only if data)
        recruit_html = self._build_recruitment_content()
        if recruit_html:
            content_parts.append(COLLAPSIBLE_PANEL_TEMPLATE.format(
                default_state="",
                panel_id="recruitment",
                title="Recruitment Plots",
                content=recruit_html,
            ))

        if not content_parts:
            content_parts.append(EMPTY_SECTION_TEMPLATE.format(
                message="No reference data provided. Supply ANI matrix and/or genome metadata."
            ))

        return TAB_SECTION_TEMPLATE.format(
            tab_id="reference",
            active_class="",
            section_title="Reference",
            content="\n".join(content_parts),
        )

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

        summary_html = DATA_SUMMARY_TEMPLATE.format(
            total_rows=len(table_df),
            known_count=s.diversity_known,
            novel_count=s.diversity_novel,
            uncertain_count=s.diversity_uncertain,
        )

        # Build quick filters section
        quick_filters_html = DATA_QUICK_FILTERS_TEMPLATE.format(
            known_count=s.known_species,
            novel_species_count=s.novel_species,
            novel_genus_count=s.novel_genus,
            ambiguous_count=s.ambiguous,
            conserved_count=s.conserved_regions,
        )

        # Column guide with similarity type (ANI or AAI)
        column_guide_html = DATA_COLUMN_GUIDE_TEMPLATE.format(
            similarity_type=self.config.similarity_type,
        )

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
            if identity_gap is None:
                identity_gap_str = "-"
            else:
                identity_gap_str = f"{abs(identity_gap):.2f}"
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
                read_id=row.get("read_id", ""),
                best_match=genome[:40] if len(genome) > 40 else genome,
                species=species[:30] if len(species) > 30 else species,
                genus=genus[:20] if len(genus) > 20 else genus,
                identity=row.get("top_hit_identity", 0),
                novelty=row.get("novelty_index", 0),
                uncertainty=row.get("placement_uncertainty", 0),
                ambiguous_hits=ambiguous_hits,
                identity_gap=identity_gap_str,
                is_novel="Yes" if is_novel else "No",
                classification=row.get("taxonomic_call", ""),
                cell_class=get_cell_class(row.get("taxonomic_call", "")),
                uncertainty_type=uncertainty_type or "-",
                inferred_uncertainty=inferred_unc_str,
                discovery_score=discovery_str,
                identity_confidence=id_conf_str,
                placement_confidence=pl_conf_str,
            )

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

    def _build_phylogeny_section(
        self,
        ani_matrix_pd: pd.DataFrame | None,
        user_tree_path: Path | None = None,
    ) -> str | None:
        """Build interactive phylogenetic tree section.

        Generates a phylogenetic tree from the ANI matrix (or loads a user-provided
        tree) and places novel clusters at estimated positions. The tree is rendered
        using D3.js for interactive exploration.

        Args:
            ani_matrix_pd: Pandas DataFrame with ANI values. Rows and columns
                should be genome accessions. Should be pre-indexed (no 'genome'
                column as index).
            user_tree_path: Optional path to user-provided Newick tree file.
                If provided, this tree is used instead of building from ANI.

        Returns:
            HTML string with phylogeny section content and D3.js visualization
            code, or None if the tree cannot be built (e.g., fewer than 3 genomes).
        """
        import json

        from metadarkmatter.visualization.report.templates import (
            PHYLOGENY_SECTION_TEMPLATE,
            PHYLOTREE_JS_TEMPLATE,
        )

        if ani_matrix_pd is None or len(ani_matrix_pd) < 3:
            return None

        try:
            from metadarkmatter.core.phylogeny import (
                ani_to_newick,
                extract_novel_clusters,
                load_user_tree,
                place_novel_clusters,
            )

            # Build or load tree
            if user_tree_path is not None:
                newick = load_user_tree(user_tree_path, ani_matrix_pd)
                tree_source_note = "Tree source: user-provided Newick file."
            else:
                newick = ani_to_newick(ani_matrix_pd)
                tree_source_note = "Tree source: neighbor-joining from ANI matrix."

            if newick is None:
                return None

            # Extract novel clusters from classifications
            novel_clusters = extract_novel_clusters(self.df, min_reads=3)

            # Place novel clusters on the tree
            if novel_clusters:
                final_newick = place_novel_clusters(
                    newick, novel_clusters, ani_matrix_pd
                )
            else:
                final_newick = newick

            # Strip BioPython comment blocks from Newick (JavaScript parser can't handle them)
            # Comments look like: NSP_001:0.1[{"classification": "novel_species", ...}]
            # We already have this data in the annotations dict
            import re
            final_newick = re.sub(r'\[[^\]]*\]', '', final_newick)

            # Build annotations dict for novel clusters
            annotations: dict[str, dict[str, Any]] = {}
            for cluster in novel_clusters:
                # Novel clusters are named by their cluster_id
                node_name = cluster.cluster_id
                annotations[node_name] = {
                    "type": cluster.classification,
                    "read_count": cluster.read_count,
                    "mean_novelty": cluster.mean_novelty,
                    "mean_uncertainty": cluster.mean_uncertainty,
                    "confidence": cluster.confidence_rating,
                    "nearest_ref": cluster.best_match_genome,
                    "est_ani": cluster.mean_identity,
                    "is_novel": True,
                }

            # Build genome labels for tree leaf display
            genome_labels = {
                acc: self._get_genome_label(acc, max_species_len=20)
                for acc in ani_matrix_pd.columns
            }

            # Build tree data JSON for D3.js visualization
            tree_data = {
                "newick": final_newick,
                "annotations": annotations,
                "genome_labels": genome_labels,
                "tip_count": len(ani_matrix_pd) + len(novel_clusters),
                "novel_count": len(novel_clusters),
                "source_note": tree_source_note,
            }

            section_html = PHYLOGENY_SECTION_TEMPLATE.format(
                tree_source_note=tree_source_note,
                tree_data_json=json.dumps(tree_data),
            )

            # PHYLOTREE_JS_TEMPLATE uses {{ and }} for JS braces, call .format() to convert
            return section_html + PHYLOTREE_JS_TEMPLATE.format()

        except ImportError as e:
            logger.warning(f"Phylogeny dependencies not available: {e}")
            return None
        except Exception as e:
            logger.warning(f"Error building phylogeny section: {e}")
            return None

    def _build_methods_footer(self) -> str:
        """Build Methods as a collapsible footer panel below all tabs."""
        return (
            '<div class="methods-footer">'
            + COLLAPSIBLE_PANEL_TEMPLATE.format(
                default_state="",
                panel_id="methods",
                title="Methods",
                content=METHODS_SECTION_TEMPLATE,
            )
            + "</div>"
        )

    def _build_histogram(
        self,
        data: list[float],
        plot_id: str,
        title: str,
        x_label: str,
        y_label: str = "Number of Reads",
        nbins: int = 30,
        color: str = "#667eea",
        thresholds: list[tuple[float, str, str, str]] | None = None,
        height: int = 400,
        x_range: list[float] | None = None,
    ) -> str:
        """Build a histogram plot and return PLOT_CONTAINER_TEMPLATE HTML.

        Args:
            data: Values to plot.
            plot_id: Unique DOM id for the plot div.
            title: Plot title.
            x_label: X-axis label.
            y_label: Y-axis label.
            nbins: Number of histogram bins.
            color: Bar fill color.
            thresholds: Optional list of (x_value, line_dash, line_color, annotation_text).
            height: Plot height in pixels.
            x_range: Optional [min, max] for x-axis.

        Returns:
            HTML string (empty if data is empty).
        """
        import plotly.graph_objects as go

        if not data:
            return ""

        fig = go.Figure()
        fig.add_trace(go.Histogram(
            x=data,
            nbinsx=nbins,
            marker_color=color,
            hovertemplate=f"{x_label}: %{{x:.2f}}<br>Count: %{{y}}<extra></extra>",
        ))

        if thresholds:
            for x_val, dash, lcolor, annotation in thresholds:
                fig.add_vline(x=x_val, line_dash=dash, line_color=lcolor,
                              annotation_text=annotation)

        layout_kwargs: dict[str, Any] = {
            "title": title,
            "xaxis_title": x_label,
            "yaxis_title": y_label,
            "template": "plotly_white",
            "height": height,
            "showlegend": False,
        }
        if x_range is not None:
            layout_kwargs["xaxis"] = {"range": x_range}
        fig.update_layout(**layout_kwargs)

        self._register_plot(plot_id, fig)
        return ""

    def _build_category_scatter(
        self,
        x_col: str,
        y_col: str,
        plot_id: str,
        title: str,
        x_label: str,
        y_label: str,
        categories: dict[str, str] | None = None,
        height: int = 500,
        max_points: int | None = None,
        y_range: list[float] | None = None,
    ) -> str:
        """Build a per-category scatter plot and register it.

        Args:
            x_col: DataFrame column for x-axis values.
            y_col: DataFrame column for y-axis values.
            plot_id: Unique DOM id.
            title: Plot title.
            x_label: X-axis label.
            y_label: Y-axis label.
            categories: Dict mapping category name to color. Defaults to BAYESIAN_CATEGORY_COLORS.
            height: Plot height in pixels.
            max_points: If set, subsample the DataFrame.
            y_range: Optional [min, max] for y-axis.

        Returns:
            Empty string (plot is registered internally).
        """
        import plotly.graph_objects as go

        if x_col not in self.df.columns or y_col not in self.df.columns:
            return ""

        if categories is None:
            categories = BAYESIAN_CATEGORY_COLORS

        plot_df = self.df
        effective_max = max_points or self.config.max_scatter_points
        if len(plot_df) > effective_max:
            plot_df = plot_df.sample(n=effective_max, seed=42)

        fig = go.Figure()
        for call, color in categories.items():
            call_df = plot_df.filter(pl.col("taxonomic_call") == call)
            if len(call_df) > 0:
                fig.add_trace(go.Scattergl(
                    x=call_df[x_col].to_list(),
                    y=call_df[y_col].to_list(),
                    mode="markers",
                    name=call,
                    marker={"color": color, "size": 4, "opacity": 0.5},
                    hovertemplate=(
                        f"<b>{call}</b><br>"
                        f"{x_label}: %{{x:.1f}}<br>"
                        f"{y_label}: %{{y:.2f}}<extra></extra>"
                    ),
                ))

        layout_kwargs: dict[str, Any] = {
            "title": title,
            "xaxis_title": x_label,
            "yaxis_title": y_label,
            "template": "plotly_white",
            "height": height,
            "legend": {"orientation": "h", "y": -0.15},
        }
        if y_range is not None:
            layout_kwargs["yaxis"] = {"range": y_range}
        fig.update_layout(**layout_kwargs)

        self._register_plot(plot_id, fig)
        return ""

    def _register_plot(self, plot_id: str, fig: go.Figure) -> None:
        """Register a plot for later JS initialization."""
        plot_json = fig.to_json()
        self._plot_data[plot_id] = plot_json

    def _build_plotly_js(self) -> str:
        """Build Plotly.js initialization code for all plots."""
        js_lines = ["<script>"]
        js_lines.append("document.addEventListener('DOMContentLoaded', function() {")

        for plot_id, plot_json in self._plot_data.items():
            js_lines.append(f"  var data_{plot_id.replace('-', '_')} = {plot_json};")
            js_lines.append(
                f"  Plotly.newPlot('{plot_id}', "
                f"data_{plot_id.replace('-', '_')}.data, "
                f"data_{plot_id.replace('-', '_')}.layout, "
                f"{{responsive: true}});"
            )

        js_lines.append("});")

        # Toggle function for ANI/AAI heatmap views
        js_lines.append("""
function toggleHeatmapView(metric, view) {
    var allView = document.getElementById(metric + '-view-all');
    var repsView = document.getElementById(metric + '-view-reps');
    var allBtn = document.getElementById(metric + '-toggle-all');
    var repsBtn = document.getElementById(metric + '-toggle-reps');
    if (!allView || !repsView) return;
    if (view === 'all') {
        allView.style.display = '';
        repsView.style.display = 'none';
        allBtn.classList.add('active');
        repsBtn.classList.remove('active');
    } else {
        allView.style.display = 'none';
        repsView.style.display = '';
        allBtn.classList.remove('active');
        repsBtn.classList.add('active');
        // Resize Plotly chart in the newly visible container
        var plotDiv = repsView.querySelector('.plotly-chart');
        if (plotDiv && typeof Plotly !== 'undefined') {
            Plotly.Plots.resize(plotDiv);
        }
    }
}
""")

        js_lines.append("</script>")

        return "\n".join(js_lines)


def generate_report(
    classifications_path: Path | str,
    output_path: Path | str,
    sample_name: str = "Sample",
    ani_matrix_path: Path | str | None = None,
    recruitment_data_path: Path | str | None = None,
    theme: str = "light",
) -> None:
    """
    Convenience function to generate a report from file paths.

    Args:
        classifications_path: Path to classifications CSV file
        output_path: Output path for HTML report
        sample_name: Sample name for report header
        ani_matrix_path: Optional path to ANI matrix CSV
        recruitment_data_path: Optional path to recruitment data CSV
        theme: Color theme ('light' or 'dark')
    """
    # Load classifications
    classifications = read_dataframe(classifications_path)

    # Load optional data
    ani_matrix = None
    if ani_matrix_path:
        ani_matrix = read_dataframe(ani_matrix_path)

    recruitment_data = None
    if recruitment_data_path:
        recruitment_data = read_dataframe(recruitment_data_path)

    # Create config
    config = ReportConfig(
        sample_name=sample_name,
        theme=theme,
    )

    # Generate report
    generator = ReportGenerator(
        classifications=classifications,
        config=config,
        recruitment_data=recruitment_data,
        ani_matrix=ani_matrix,
    )
    generator.generate(output_path)
