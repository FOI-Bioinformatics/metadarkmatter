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
    ANI_NOT_PROVIDED_MESSAGE,
    CATEGORY_BREAKDOWN_TEMPLATE,
    DATA_COLUMN_GUIDE_TEMPLATE,
    DATA_QUICK_FILTERS_TEMPLATE,
    DATA_SUMMARY_TEMPLATE,
    DATA_TABLE_JS,
    DATA_TABLE_TEMPLATE,
    DISTRIBUTIONS_SUMMARY_TEMPLATE,
    DIVERSITY_SUMMARY_TEMPLATE,
    EMPTY_SECTION_TEMPLATE,
    ENHANCED_SCORING_CONFIDENCE_TEMPLATE,
    ENHANCED_SCORING_DISCOVERY_GUIDE_TEMPLATE,
    ENHANCED_SCORING_SUMMARY_TEMPLATE,
    ENHANCED_SCORING_UNCERTAINTY_TYPES_TEMPLATE,
    GENOME_HIGHLIGHTS_TEMPLATE,
    GENOME_INTERPRETATION_TEMPLATE,
    GENOMES_SUMMARY_TEMPLATE,
    METRIC_CARD_TEMPLATE,
    METRIC_CARDS_CONTAINER,
    METHODS_SECTION_TEMPLATE,
    PLOT_CONTAINER_SIMPLE_TEMPLATE,
    PLOT_CONTAINER_TEMPLATE,
    PLOT_ROW_TEMPLATE,
    RECRUITMENT_NOT_PROVIDED_MESSAGE,
    REPORT_BASE_TEMPLATE,
    SCATTER_INTERPRETATION_TEMPLATE,
    TAB_NAVIGATION_JS,
    TAB_SECTION_TEMPLATE,
    TABLE_ROW_TEMPLATE,
    BAYESIAN_SUMMARY_TEMPLATE,
    BAYESIAN_INTERPRETATION_TEMPLATE,
    FAMILY_VALIDATION_SECTION_TEMPLATE,
    get_cell_class,
)
from metadarkmatter.visualization.report.novel_section import (
    NOVEL_CONFIDENCE_GUIDE_TEMPLATE,
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
    # Bayesian posterior metrics (optional, only when --bayesian used)
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
            "mean_novelty_index": self.mean_novelty_index,
            "mean_placement_uncertainty": self.mean_placement_uncertainty,
            "mean_top_hit_identity": self.mean_top_hit_identity,
            "novel_percentage": self.novel_percentage,
            "diversity_known_pct": self.diversity_known_pct,
            "diversity_novel_pct": self.diversity_novel_pct,
            "diversity_uncertain_pct": self.diversity_uncertain_pct,
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
            conserved_regions + unclassified + off_target
        )

        # Check for enhanced scoring columns
        has_enhanced_scoring = "alignment_quality" in self.df.columns
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
            if "bayesian_category" in self.df.columns:
                map_agreement_count = len(
                    self.df.filter(
                        pl.col("bayesian_category") == pl.col("taxonomic_call")
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
        # Generate all sections (order matches tab navigation)
        content_sections = []

        # Overview tab (always visible first)
        content_sections.append(self._build_overview_section())

        # Distributions tab
        content_sections.append(self._build_distributions_section())

        # Species tab (if metadata provided)
        content_sections.append(self._build_species_section())

        # Genomes tab
        content_sections.append(self._build_genomes_section())

        # Novel Diversity tab (only if there are novel reads)
        if self.summary.diversity_novel > 0:
            content_sections.append(self._build_novel_diversity_section())

        # Family Validation tab (only if family validation was active)
        if self.summary.has_family_validation:
            content_sections.append(self._build_family_validation_section())

        # Reference ANI tab
        content_sections.append(self._build_ani_section())

        # Reference AAI tab
        content_sections.append(self._build_aai_section())

        # Phylogeny tab (only if ANI matrix provided with >= 3 genomes and not skipped)
        self._phylogeny_html = None
        if not self.config.skip_phylogeny and self.ani_matrix is not None and len(self.ani_matrix) >= 3:
            ani_pd = self.ani_matrix.to_pandas()
            if "genome" in ani_pd.columns:
                ani_pd = ani_pd.set_index("genome")
            self._phylogeny_html = self._build_phylogeny_section(
                ani_pd,
                user_tree_path=self.config.user_tree_path,
            )
            if self._phylogeny_html:
                content_sections.append(
                    TAB_SECTION_TEMPLATE.format(
                        tab_id="phylogeny",
                        active_class="",
                        section_title="Phylogeny",
                        content=self._phylogeny_html,
                    )
                )

        # Discovery Scores tab (only if enhanced scoring data available)
        if self.summary.has_enhanced_scoring or self.summary.has_inferred_uncertainty:
            content_sections.append(self._build_enhanced_scoring_section())

        # Bayesian Confidence tab (only if posterior columns present)
        if self.summary.has_bayesian:
            content_sections.append(self._build_bayesian_section())

        # Recruitment tab
        content_sections.append(self._build_recruitment_section())

        # Data table tab
        content_sections.append(self._build_data_section())

        # Methods tab (always visible - comprehensive documentation)
        content_sections.append(self._build_methods_section())

        content = "\n".join(content_sections)

        # Build dynamic navigation
        navigation = self._build_navigation()

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
            plotly_js=plotly_js,
            js_scripts=js_scripts,
        )

    def _build_navigation(self) -> str:
        """Build dynamic navigation based on available data."""
        # Tab order: General -> Specific -> Novel -> Reference -> Technical -> Data
        tabs = [
            ("overview", "Overview", True),
            ("distributions", "Distributions", False),
            ("species", "Species", False),
            ("genomes", "Genomes", False),
        ]

        # Add Novel Diversity tab if there are novel reads
        if self.summary.diversity_novel > 0:
            tabs.append(("novel-diversity", "Novel Diversity", False))

        # Add Family Validation tab if family validation was active
        if self.summary.has_family_validation:
            tabs.append(("family-validation", "Family Validation", False))

        # Reference matrices (renamed for clarity)
        tabs.extend([
            ("ani", "Reference ANI", False),
            ("aai", "Reference AAI", False),
        ])

        # Add Phylogeny tab if phylogeny section was built
        if hasattr(self, "_phylogeny_html") and self._phylogeny_html is not None:
            tabs.append(("phylogeny", "Phylogeny", False))

        # Add Discovery Scores tab if data available (renamed from Enhanced Scoring)
        if self.summary.has_enhanced_scoring or self.summary.has_inferred_uncertainty:
            tabs.append(("enhanced-scoring", "Discovery Scores", False))

        # Add Bayesian Confidence tab if posterior data available
        if self.summary.has_bayesian:
            tabs.append(("bayesian", "Bayesian Confidence", False))

        tabs.extend([
            ("recruitment", "Recruitment", False),
            ("data", "Data", False),
            ("methods", "Methods", False),
        ])

        nav_items = []
        for tab_id, label, is_active in tabs:
            active_class = " active" if is_active else ""
            nav_items.append(
                f'        <button class="tab-btn{active_class}" '
                f"onclick=\"showTab('{tab_id}')\">{label}</button>"
            )

        return "\n".join(nav_items)

    def _build_overview_section(self) -> str:
        """Build the overview tab section with hero diversity summary."""
        s = self.summary

        # Calculate percentages for all categories
        total = s.total_reads if s.total_reads > 0 else 1  # Avoid division by zero

        # Diversity summary - the hero section
        diversity_html = DIVERSITY_SUMMARY_TEMPLATE.format(
            known_pct=s.diversity_known_pct,
            known_count=s.diversity_known,
            novel_pct=s.diversity_novel_pct,
            novel_count=s.diversity_novel,
            uncertain_pct=s.diversity_uncertain_pct,
            uncertain_count=s.diversity_uncertain,
        )

        # Category breakdown with detailed counts
        category_html = CATEGORY_BREAKDOWN_TEMPLATE.format(
            known_species=s.known_species,
            known_species_pct=s.known_species / total * 100,
            novel_species=s.novel_species,
            novel_species_pct=s.novel_species / total * 100,
            novel_genus=s.novel_genus,
            novel_genus_pct=s.novel_genus / total * 100,
            species_boundary=s.species_boundary,
            species_boundary_pct=s.species_boundary / total * 100,
            ambiguous=s.ambiguous,
            ambiguous_pct=s.ambiguous / total * 100,
            conserved_regions=s.conserved_regions,
            conserved_pct=s.conserved_regions / total * 100,
        )

        # Sunburst chart for interactive exploration - larger size
        sunburst = DiversitySunburstChart(
            s.to_dict(),
            config=PlotConfig(width=800, height=650),
            title="Diversity Classification",
        )
        sunburst_id = "plot-diversity-sunburst"
        self._register_plot(sunburst_id, sunburst.create_figure())

        # Build plot row with sunburst chart
        plots_row = PLOT_ROW_TEMPLATE.format(
            plots=(
                PLOT_CONTAINER_TEMPLATE.format(
                    extra_class="full-width",
                    title="Interactive Diversity Chart",
                    description="Click segments to explore. Inner ring shows diversity status (Known/Novel/Uncertain), outer ring shows detailed categories.",
                    plot_id=sunburst_id,
                )
            )
        )

        content = diversity_html + category_html + plots_row

        return TAB_SECTION_TEMPLATE.format(
            tab_id="overview",
            active_class="active",
            section_title="Overview",
            content=content,
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

    def _build_distributions_section(self) -> str:
        """Build the distributions tab section with improved UX."""
        s = self.summary

        # Calculate interpretations based on mean values
        if s.mean_novelty_index < 5:
            novelty_interp = "Mostly known taxa"
        elif s.mean_novelty_index < 15:
            novelty_interp = "Some novel diversity"
        else:
            novelty_interp = "High novel diversity"

        if s.mean_placement_uncertainty < 2:
            uncertainty_interp = "Confident placements"
        elif s.mean_placement_uncertainty < 5:
            uncertainty_interp = "Moderate confidence"
        else:
            uncertainty_interp = "Many ambiguous hits"

        # Count confident classifications (low uncertainty)
        confident_count = len(self.df.filter(
            self.df["placement_uncertainty"] < 2
        ))
        confident_pct = (confident_count / s.total_reads * 100) if s.total_reads else 0

        # Summary section
        summary_html = DISTRIBUTIONS_SUMMARY_TEMPLATE.format(
            mean_identity=s.mean_top_hit_identity,
            mean_novelty=s.mean_novelty_index,
            novelty_interpretation=novelty_interp,
            mean_uncertainty=s.mean_placement_uncertainty,
            uncertainty_interpretation=uncertainty_interp,
            confident_pct=confident_pct,
        )

        # 2D scatter (the key diagnostic plot) - show first as it's the primary view
        scatter = NoveltyUncertaintyScatter(
            self.df,
            config=PlotConfig(width=1100, height=650),
            thresholds=self.config.thresholds,
            max_points=self.config.max_scatter_points,
        )
        scatter_id = "plot-scatter-2d"
        self._register_plot(scatter_id, scatter.create_figure())

        # Interpretation guide
        interpretation_html = SCATTER_INTERPRETATION_TEMPLATE

        scatter_plot = PLOT_CONTAINER_TEMPLATE.format(
            extra_class="full-width",
            title="Novelty vs Uncertainty Landscape",
            description=(
                "Each point represents a read. Position indicates divergence from references (x-axis) "
                "and confidence in placement (y-axis). Colors show classification categories. "
                "Use dropdown to toggle single-hit reads."
            ),
            plot_id=scatter_id,
        )

        # Confidence vs Novelty scatter plot
        confidence_scatter = ConfidenceNoveltyScatter(
            self.df,
            config=PlotConfig(width=1100, height=650),
            thresholds=self.config.thresholds,
            max_points=self.config.max_scatter_points,
        )
        confidence_id = "plot-confidence-novelty"
        self._register_plot(confidence_id, confidence_scatter.create_figure())

        confidence_plot = PLOT_CONTAINER_TEMPLATE.format(
            extra_class="full-width",
            title="Confidence Score vs Novelty Index",
            description=(
                "Confidence score (100 - novelty - uncertainty) vs novelty index. "
                "High confidence with high novelty indicates reliable novel species detection. "
                "Use dropdown to toggle single-hit reads."
            ),
            plot_id=confidence_id,
        )

        # Novelty histogram
        novelty_hist = NoveltyHistogram(
            self.df,
            config=PlotConfig(width=550, height=400),
            thresholds=self.config.thresholds,
            bin_size=self.config.novelty_bin_size,
        )
        novelty_id = "plot-novelty-hist"
        self._register_plot(novelty_id, novelty_hist.create_figure())

        # Uncertainty histogram
        uncertainty_hist = UncertaintyHistogram(
            self.df,
            config=PlotConfig(width=550, height=400),
            thresholds=self.config.thresholds,
            bin_size=self.config.uncertainty_bin_size,
        )
        uncertainty_id = "plot-uncertainty-hist"
        self._register_plot(uncertainty_id, uncertainty_hist.create_figure())

        # Build histogram row
        hist_row = PLOT_ROW_TEMPLATE.format(
            plots=(
                PLOT_CONTAINER_TEMPLATE.format(
                    extra_class="half-width",
                    title="Novelty Index Distribution",
                    description=(
                        "Low values (0-5) = known species. "
                        "Medium (5-20) = novel species. "
                        "High (>20) = novel genus or more divergent."
                    ),
                    plot_id=novelty_id,
                )
                + PLOT_CONTAINER_TEMPLATE.format(
                    extra_class="half-width",
                    title="Placement Uncertainty Distribution",
                    description=(
                        "Low values (<2) = confident placement. "
                        "Medium (2-5) = species boundary zone. "
                        "High (>5) = conserved regions or ambiguous."
                    ),
                    plot_id=uncertainty_id,
                )
            )
        )

        content = summary_html + interpretation_html + scatter_plot + confidence_plot + hist_row

        return TAB_SECTION_TEMPLATE.format(
            tab_id="distributions",
            active_class="",
            section_title="Distributions",
            content=content,
        )

    def _build_enhanced_scoring_section(self) -> str:
        """Build the enhanced scoring tab section with confidence metrics and discovery scores."""
        import plotly.graph_objects as go

        s = self.summary
        content_parts = []

        # Summary section with key metrics
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

        # Confidence dimensions section (only if enhanced scoring)
        if s.has_enhanced_scoring:
            confidence_html = ENHANCED_SCORING_CONFIDENCE_TEMPLATE.format(
                mean_identity_confidence=s.mean_identity_confidence or 0.0,
                mean_placement_confidence=s.mean_placement_confidence or 0.0,
                mean_alignment_quality=s.mean_alignment_quality or 0.0,
            )
            content_parts.append(confidence_html)

        # Uncertainty types breakdown
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

        # Discovery score interpretation guide
        if s.has_enhanced_scoring and s.novel_with_discovery_score > 0:
            content_parts.append(ENHANCED_SCORING_DISCOVERY_GUIDE_TEMPLATE)

        # Visualizations

        # Discovery Score Distribution (histogram for novel reads)
        if s.has_enhanced_scoring and "discovery_score" in self.df.columns:
            discovery_df = self.df.filter(pl.col("discovery_score").is_not_null())
            if len(discovery_df) > 0:
                scores = discovery_df["discovery_score"].to_list()

                discovery_fig = go.Figure()
                discovery_fig.add_trace(go.Histogram(
                    x=scores,
                    nbinsx=25,
                    marker_color="#667eea",
                    hovertemplate="Score: %{x:.0f}<br>Count: %{y}<extra></extra>",
                ))

                # Add threshold lines
                discovery_fig.add_vline(x=75, line_dash="dash", line_color="#22c55e",
                                        annotation_text="High Priority (75+)")
                discovery_fig.add_vline(x=50, line_dash="dot", line_color="#f59e0b",
                                        annotation_text="Moderate (50+)")
                discovery_fig.add_vline(x=25, line_dash="dot", line_color="#ef4444",
                                        annotation_text="Low (25+)")

                discovery_fig.update_layout(
                    title="Discovery Score Distribution (Novel Reads)",
                    xaxis_title="Discovery Score",
                    yaxis_title="Number of Reads",
                    template="plotly_white",
                    height=400,
                    showlegend=False,
                )
                discovery_hist_id = "plot-discovery-hist"
                self._register_plot(discovery_hist_id, discovery_fig)

                content_parts.append(PLOT_CONTAINER_TEMPLATE.format(
                    extra_class="full-width",
                    title="Discovery Score Distribution",
                    description=(
                        "Distribution of discovery scores for novel reads. Higher scores indicate "
                        "more reliable discoveries. Green line marks high-priority threshold (75+)."
                    ),
                    plot_id=discovery_hist_id,
                ))

        # Confidence Dimensions Scatter Plot
        if s.has_enhanced_scoring and "identity_confidence" in self.df.columns:
            # Sample for performance
            plot_df = self.df
            if len(plot_df) > self.config.max_scatter_points:
                plot_df = plot_df.sample(n=self.config.max_scatter_points, seed=42)

            conf_fig = go.Figure()

            # Color by taxonomic call
            for call, color in [
                ("Known Species", "#22c55e"),
                ("Novel Species", "#f59e0b"),
                ("Novel Genus", "#ef4444"),
                ("Ambiguous", "#94a3b8"),
            ]:
                call_df = plot_df.filter(pl.col("taxonomic_call") == call)
                if len(call_df) > 0:
                    conf_fig.add_trace(go.Scattergl(
                        x=call_df["identity_confidence"].to_list(),
                        y=call_df["placement_confidence"].to_list(),
                        mode="markers",
                        name=call,
                        marker={"color": color, "size": 4, "opacity": 0.6},
                        hovertemplate=(
                            f"<b>{call}</b><br>"
                            "Identity Conf: %{x:.1f}<br>"
                            "Placement Conf: %{y:.1f}<extra></extra>"
                        ),
                    ))

            conf_fig.update_layout(
                title="Confidence Dimensions",
                xaxis_title="Identity Confidence",
                yaxis_title="Placement Confidence",
                template="plotly_white",
                height=500,
                legend={"orientation": "h", "y": -0.15},
            )
            conf_scatter_id = "plot-confidence-dimensions"
            self._register_plot(conf_scatter_id, conf_fig)

            content_parts.append(PLOT_CONTAINER_TEMPLATE.format(
                extra_class="full-width",
                title="Confidence Dimensions Scatter",
                description=(
                    "Identity confidence (reliability of identity measurement) vs "
                    "placement confidence (confidence in genome assignment). "
                    "High values on both axes indicate reliable classifications."
                ),
                plot_id=conf_scatter_id,
            ))

        # Inferred vs Measured Uncertainty Comparison
        if s.has_inferred_uncertainty and "inferred_uncertainty" in self.df.columns:
            # Create comparison histogram
            inferred_df = self.df.filter(pl.col("inferred_uncertainty").is_not_null())
            measured_df = self.df.filter(
                (pl.col("uncertainty_type") == "measured") &
                (pl.col("placement_uncertainty").is_not_null())
            )

            unc_fig = go.Figure()

            if len(inferred_df) > 0:
                unc_fig.add_trace(go.Histogram(
                    x=inferred_df["inferred_uncertainty"].to_list(),
                    name="Inferred (single-hit)",
                    marker_color="#f59e0b",
                    opacity=0.7,
                    nbinsx=30,
                ))

            if len(measured_df) > 0:
                unc_fig.add_trace(go.Histogram(
                    x=measured_df["placement_uncertainty"].to_list(),
                    name="Measured (multi-hit)",
                    marker_color="#22c55e",
                    opacity=0.7,
                    nbinsx=30,
                ))

            unc_fig.update_layout(
                title="Uncertainty Distribution by Type",
                xaxis_title="Uncertainty (%)",
                yaxis_title="Number of Reads",
                template="plotly_white",
                height=400,
                barmode="overlay",
                legend={"orientation": "h", "y": -0.15},
            )
            unc_hist_id = "plot-uncertainty-comparison"
            self._register_plot(unc_hist_id, unc_fig)

            content_parts.append(PLOT_CONTAINER_TEMPLATE.format(
                extra_class="full-width",
                title="Uncertainty by Type",
                description=(
                    "Comparison of measured uncertainty (from ANI between competing hits) "
                    "vs inferred uncertainty (estimated from novelty for single-hit reads). "
                    "Inferred uncertainty prevents false confidence in single-hit reads."
                ),
                plot_id=unc_hist_id,
            ))

        content = "\n".join(content_parts)

        return TAB_SECTION_TEMPLATE.format(
            tab_id="enhanced-scoring",
            active_class="",
            section_title="Discovery Scores",
            content=content,
        )

    def _build_bayesian_section(self) -> str:
        """Build the Bayesian confidence tab with posterior analysis and entropy."""
        import plotly.graph_objects as go

        s = self.summary
        content_parts = []

        # Summary cards
        summary_html = BAYESIAN_SUMMARY_TEMPLATE.format(
            mean_entropy=s.mean_posterior_entropy,
            high_confidence_pct=s.high_confidence_pct,
            high_confidence_count=s.high_confidence_count,
            map_agreement_pct=s.map_agreement_pct,
            map_agreement_count=s.map_agreement_count,
            total_reads=s.total_reads,
            boundary_pct=s.boundary_pct,
            boundary_count=s.boundary_count,
        )
        content_parts.append(summary_html)

        # Interpretation guide
        content_parts.append(BAYESIAN_INTERPRETATION_TEMPLATE)

        # --- Posterior Entropy Distribution ---
        if "posterior_entropy" in self.df.columns:
            entropy_vals = self.df["posterior_entropy"].to_list()

            entropy_fig = go.Figure()
            entropy_fig.add_trace(go.Histogram(
                x=entropy_vals,
                nbinsx=40,
                marker_color="#667eea",
                hovertemplate="Entropy: %{x:.2f}<br>Count: %{y}<extra></extra>",
            ))

            # Add threshold lines
            entropy_fig.add_vline(
                x=0.5, line_dash="dash", line_color="#22c55e",
                annotation_text="Very confident",
            )
            entropy_fig.add_vline(
                x=1.0, line_dash="dot", line_color="#f59e0b",
                annotation_text="Confident",
            )
            entropy_fig.add_vline(
                x=1.5, line_dash="dot", line_color="#ef4444",
                annotation_text="Uncertain",
            )

            entropy_fig.update_layout(
                title="Posterior Entropy Distribution",
                xaxis_title="Shannon Entropy (bits)",
                yaxis_title="Number of Reads",
                template="plotly_white",
                height=400,
                showlegend=False,
                xaxis={"range": [0, 2.1]},
            )
            entropy_id = "plot-bayesian-entropy"
            self._register_plot(entropy_id, entropy_fig)

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
        if "bayesian_category" in self.df.columns:
            category_colors = {
                "Known Species": "#22c55e",
                "Novel Species": "#f59e0b",
                "Novel Genus": "#ef4444",
                "Ambiguous": "#94a3b8",
            }
            posterior_cols = [
                ("p_known_species", "Known Species"),
                ("p_novel_species", "Novel Species"),
                ("p_novel_genus", "Novel Genus"),
                ("p_ambiguous", "Ambiguous"),
            ]

            # Mean posteriors grouped by Bayesian MAP category
            bar_fig = go.Figure()
            for col_name, cat_label in posterior_cols:
                if col_name not in self.df.columns:
                    continue
                means = []
                cats = []
                for map_cat in ["Known Species", "Novel Species", "Novel Genus", "Ambiguous"]:
                    cat_df = self.df.filter(pl.col("bayesian_category") == map_cat)
                    if len(cat_df) > 0:
                        cats.append(map_cat)
                        means.append(cat_df[col_name].mean() or 0.0)
                if means:
                    bar_fig.add_trace(go.Bar(
                        x=cats,
                        y=means,
                        name=cat_label,
                        marker_color=category_colors.get(cat_label, "#94a3b8"),
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

        # --- MAP Confidence vs Novelty Scatter ---
        if "posterior_entropy" in self.df.columns and "novelty_index" in self.df.columns:
            plot_df = self.df
            if len(plot_df) > self.config.max_scatter_points:
                plot_df = plot_df.sample(n=self.config.max_scatter_points, seed=42)

            scatter_fig = go.Figure()

            category_colors_scatter = {
                "Known Species": "#22c55e",
                "Novel Species": "#f59e0b",
                "Novel Genus": "#ef4444",
                "Ambiguous": "#94a3b8",
            }
            for call, color in category_colors_scatter.items():
                call_df = plot_df.filter(pl.col("bayesian_category") == call)
                if len(call_df) > 0:
                    scatter_fig.add_trace(go.Scattergl(
                        x=call_df["novelty_index"].to_list(),
                        y=call_df["posterior_entropy"].to_list(),
                        mode="markers",
                        name=call,
                        marker={"color": color, "size": 4, "opacity": 0.5},
                        hovertemplate=(
                            f"<b>{call}</b><br>"
                            "Novelty: %{{x:.1f}}%<br>"
                            "Entropy: %{{y:.2f}}<extra></extra>"
                        ),
                    ))

            scatter_fig.update_layout(
                title="Classification Confidence Landscape",
                xaxis_title="Novelty Index (%)",
                yaxis_title="Posterior Entropy (bits)",
                template="plotly_white",
                height=500,
                legend={"orientation": "h", "y": -0.15},
                yaxis={"range": [0, 2.1]},
            )
            scatter_id = "plot-bayesian-landscape"
            self._register_plot(scatter_id, scatter_fig)

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

        content = "\n".join(content_parts)

        return TAB_SECTION_TEMPLATE.format(
            tab_id="bayesian",
            active_class="",
            section_title="Bayesian Confidence",
            content=content,
        )

    def _build_novel_diversity_section(self) -> str:
        """Build the novel diversity tab section with cluster analysis."""
        from metadarkmatter.core.novel_diversity import NovelDiversityAnalyzer

        content_parts = []

        # Create analyzer and cluster novel reads
        try:
            analyzer = NovelDiversityAnalyzer(
                classifications=self.df,
                metadata=self.genome_metadata,
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

        # Build summary section
        content_parts.append(build_novel_summary_html(summary))

        # Add confidence guide
        content_parts.append(NOVEL_CONFIDENCE_GUIDE_TEMPLATE)

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

                    content_parts.append(PHYLOGENETIC_HEATMAP_INTRO_TEMPLATE.format(
                        n_references=phylo_metadata.get("n_references", 0),
                        n_clusters=n_clusters_shown,
                        clusters_note=clusters_note,
                    ))

                    # Register and add the heatmap
                    phylo_heatmap_id = f"plot-novel-phylo-heatmap-{sim_type.lower()}"
                    self._register_plot(phylo_heatmap_id, phylo_fig)

                    content_parts.append(PLOT_CONTAINER_TEMPLATE.format(
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
                        content_parts.append(PHYLOGENETIC_HEATMAP_LEGEND_TEMPLATE)

            except Exception as e:
                logger.warning(f"Could not build {sim_type} phylogenetic context heatmap: {e}")

        # Build cluster table (includes phylogenetic placement)
        content_parts.append(build_cluster_table_html(clusters))

        content = "\n".join(content_parts)

        return TAB_SECTION_TEMPLATE.format(
            tab_id="novel-diversity",
            active_class="",
            section_title="Novel Diversity",
            content=content,
        )

    def _build_recruitment_section(self) -> str:
        """Build the recruitment plots tab section."""
        if self.recruitment_data is None or len(self.recruitment_data) == 0:
            content = EMPTY_SECTION_TEMPLATE.format(
                message=RECRUITMENT_NOT_PROVIDED_MESSAGE
            )
        else:
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

            content = (
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

        return TAB_SECTION_TEMPLATE.format(
            tab_id="recruitment",
            active_class="",
            section_title="Recruitment Plots",
            content=content,
        )

    def _build_species_section(self) -> str:
        """Build the species breakdown tab section."""
        if self.genome_metadata is None:
            content = EMPTY_SECTION_TEMPLATE.format(
                message=(
                    "Species breakdown not available. "
                    "Use --metadata option with genome_metadata.tsv "
                    "to enable species-level aggregation."
                )
            )
        else:
            import plotly.graph_objects as go

            # Join metadata to get species information
            if "species" not in self.df.columns:
                # Need to join species information
                enriched_df = self.genome_metadata.join_classifications(self.df)
            else:
                enriched_df = self.df

            # Aggregate by species
            species_counts = (
                enriched_df.group_by("species")
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
                content = EMPTY_SECTION_TEMPLATE.format(
                    message="No species data available."
                )
            else:
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

                # Pie chart showing species proportion
                pie_fig = go.Figure()
                pie_fig.add_trace(
                    go.Pie(
                        labels=species_counts["species"].head(10).to_list(),
                        values=species_counts["count"].head(10).to_list(),
                        textinfo="percent+label",
                        textposition="inside",
                        hole=0.4,
                        hovertemplate=(
                            "<b>%{label}</b><br>"
                            "Reads: %{value:,}<br>"
                            "Proportion: %{percent}<extra></extra>"
                        ),
                    )
                )
                pie_fig.update_layout(
                    title="Top 10 Species Composition",
                    template="plotly_white",
                    height=500,
                    showlegend=True,
                    legend=dict(
                        orientation="h",
                        yanchor="bottom",
                        y=-0.3,
                        xanchor="center",
                        x=0.5,
                    ),
                )
                species_pie_id = "plot-species-pie"
                self._register_plot(species_pie_id, pie_fig)

                # Summary statistics
                total_species = species_counts["species"].n_unique()

                content = (
                    f'<div class="metric-cards">'
                    f'<div class="metric-card">'
                    f'<div class="metric-value">{total_species}</div>'
                    f'<div class="metric-label">Unique Species</div>'
                    f'<div class="metric-subtext">Detected in sample</div>'
                    f'</div>'
                    f'<div class="metric-card">'
                    f'<div class="metric-value">{self.genome_metadata.species_count}</div>'
                    f'<div class="metric-label">Reference Species</div>'
                    f'<div class="metric-subtext">In metadata database</div>'
                    f'</div>'
                    f'</div>'
                    + PLOT_ROW_TEMPLATE.format(
                        plots=(
                            PLOT_CONTAINER_TEMPLATE.format(
                                extra_class="half-width",
                                title="Species Composition",
                                description=(
                                    "Pie chart showing the proportion of reads "
                                    "assigned to each species."
                                ),
                                plot_id=species_pie_id,
                            )
                            + PLOT_CONTAINER_SIMPLE_TEMPLATE.format(
                                extra_class="half-width",
                                plot_id=species_bar_id,
                            )
                        )
                    )
                )

        return TAB_SECTION_TEMPLATE.format(
            tab_id="species",
            active_class="",
            section_title="Species Breakdown",
            content=content,
        )

    def _build_genomes_section(self) -> str:
        """Build the genomes breakdown tab section with improved UX."""
        # Per-genome read counts with statistics
        genome_counts = (
            self.df.group_by("best_match_genome")
            .agg([
                pl.len().alias("count"),
                pl.col("novelty_index").mean().alias("mean_novelty"),
                pl.col("top_hit_identity").mean().alias("mean_identity"),
            ])
            .sort("count", descending=True)
        )

        if len(genome_counts) == 0:
            content = EMPTY_SECTION_TEMPLATE.format(
                message="No genome data available."
            )
        else:
            import plotly.graph_objects as go

            # Calculate summary statistics
            total_genomes = len(genome_counts)
            ref_genomes = (
                self.genome_metadata.genome_count
                if self.genome_metadata else total_genomes
            )

            # Top genome stats
            top_genome_acc = genome_counts["best_match_genome"][0]
            top_genome_count = genome_counts["count"][0]
            top_genome_pct = (
                (top_genome_count / self.summary.total_reads * 100)
                if self.summary.total_reads else 0
            )
            top_genome_name = self._get_genome_label(top_genome_acc, max_species_len=25)

            # Top 3 genomes concentration
            top3_count = genome_counts["count"].head(3).sum()
            top3_pct = (
                (top3_count / self.summary.total_reads * 100)
                if self.summary.total_reads else 0
            )
            concentration_hint = (
                "Concentrated" if top3_pct > 70 else
                "Moderate spread" if top3_pct > 40 else
                "Widely distributed"
            )

            # Median reads per genome
            median_reads = int(genome_counts["count"].median() or 0)

            # Build summary section
            summary_html = GENOMES_SUMMARY_TEMPLATE.format(
                total_genomes=total_genomes,
                ref_genomes=ref_genomes,
                top_genome_pct=top_genome_pct,
                top_genome_name=top_genome_name[:30],
                top3_pct=top3_pct,
                concentration_hint=concentration_hint,
                median_reads=median_reads,
            )

            # Find key genomes for highlights
            # Helper to get species name
            def get_species(acc: str) -> str:
                if self.genome_metadata:
                    species = self.genome_metadata.get_species(acc)
                    if species and species != "Unknown":
                        return species
                return "Unknown species"

            # Most abundant (already have it)
            top_identity = genome_counts["mean_identity"][0]
            top_species = get_species(top_genome_acc)

            # Most novel (highest mean novelty with at least 5 reads)
            novel_df = genome_counts.filter(pl.col("count") >= 5)
            if len(novel_df) > 0:
                most_novel_idx = novel_df["mean_novelty"].arg_max()
                novel_genome_acc = novel_df["best_match_genome"][most_novel_idx]
                novel_reads = novel_df["count"][most_novel_idx]
                novel_novelty = novel_df["mean_novelty"][most_novel_idx]
            else:
                novel_genome_acc = top_genome_acc
                novel_reads = top_genome_count
                novel_novelty = genome_counts["mean_novelty"][0]
            novel_species = get_species(novel_genome_acc)

            # Most confident (highest mean identity with at least 5 reads)
            confident_df = genome_counts.filter(pl.col("count") >= 5)
            if len(confident_df) > 0:
                most_confident_idx = confident_df["mean_identity"].arg_max()
                confident_genome_acc = confident_df["best_match_genome"][
                    most_confident_idx
                ]
                confident_reads = confident_df["count"][most_confident_idx]
                confident_identity = confident_df["mean_identity"][most_confident_idx]
            else:
                confident_genome_acc = top_genome_acc
                confident_reads = top_genome_count
                confident_identity = top_identity
            confident_species = get_species(confident_genome_acc)

            # Build highlights section
            highlights_html = GENOME_HIGHLIGHTS_TEMPLATE.format(
                top_species=top_species,
                top_accession=top_genome_acc,
                top_reads=top_genome_count,
                top_identity=top_identity,
                novel_species=novel_species,
                novel_accession=novel_genome_acc,
                novel_reads=novel_reads,
                novel_novelty=novel_novelty,
                confident_species=confident_species,
                confident_accession=confident_genome_acc,
                confident_reads=confident_reads,
                confident_identity=confident_identity,
            )

            # Interpretation guide
            interpretation_html = GENOME_INTERPRETATION_TEMPLATE

            # Get top 20 genomes for plotting
            genome_counts_top20 = genome_counts.head(20)
            genome_accessions = genome_counts_top20["best_match_genome"].to_list()
            genome_labels = [
                self._get_genome_label(acc, max_species_len=30)
                for acc in genome_accessions
            ]

            # Horizontal bar chart of genome read counts
            bar_fig = go.Figure()
            bar_fig.add_trace(
                go.Bar(
                    y=genome_labels[::-1],
                    x=genome_counts_top20["count"].to_list()[::-1],
                    orientation="h",
                    marker_color="#667eea",
                    text=[
                        f"{c:,}"
                        for c in genome_counts_top20["count"].to_list()[::-1]
                    ],
                    textposition="outside",
                    hovertemplate=(
                        "<b>%{y}</b><br>"
                        "Reads: %{x:,}<extra></extra>"
                    ),
                )
            )
            bar_fig.update_layout(
                title="Top 20 Reference Genomes by Read Count",
                xaxis_title="Number of Reads",
                yaxis_title="",
                template="plotly_white",
                height=max(400, 35 * len(genome_counts_top20)),
                margin={"l": 350, "r": 40, "t": 60, "b": 60},
            )
            genome_bar_id = "plot-genome-bar"
            self._register_plot(genome_bar_id, bar_fig)

            # Box plot of identity distribution per genome
            # Get identity data for top 10 genomes
            top_genomes = genome_counts["best_match_genome"].head(10).to_list()
            box_data = self.df.filter(pl.col("best_match_genome").is_in(top_genomes))

            box_fig = go.Figure()
            for genome in top_genomes:
                gdata = box_data.filter(pl.col("best_match_genome") == genome)
                # Use species-annotated label for box plot
                label = self._get_genome_label(genome, max_species_len=20)
                box_fig.add_trace(
                    go.Box(
                        y=gdata["top_hit_identity"].to_list(),
                        name=label[:40] + ("..." if len(label) > 40 else ""),
                        boxpoints=False,
                        hoverinfo="y+name",
                    )
                )

            box_fig.update_layout(
                title="Identity Distribution by Genome (Top 10)",
                yaxis_title="Percent Identity (%)",
                template="plotly_white",
                showlegend=False,
                height=500,
            )
            genome_box_id = "plot-genome-box"
            self._register_plot(genome_box_id, box_fig)

            # Build content with new templates
            plots_html = (
                PLOT_CONTAINER_TEMPLATE.format(
                    extra_class="full-width",
                    title="Read Counts by Reference Genome",
                    description=(
                        "Bar chart showing the number of reads assigned "
                        "to each reference genome."
                    ),
                    plot_id=genome_bar_id,
                )
                + PLOT_CONTAINER_TEMPLATE.format(
                    extra_class="full-width",
                    title="Identity Distribution by Genome",
                    description=(
                        "Box plots showing the distribution of alignment identity "
                        "for top genomes."
                    ),
                    plot_id=genome_box_id,
                )
            )

            content = summary_html + highlights_html + interpretation_html + plots_html

        return TAB_SECTION_TEMPLATE.format(
            tab_id="genomes",
            active_class="",
            section_title="Genome Breakdown",
            content=content,
        )

    def _build_family_validation_section(self) -> str:
        """Build the Family Validation tab content."""
        import plotly.graph_objects as go

        total = self.summary.total_reads
        off_target = self.summary.off_target
        in_family = total - off_target
        off_target_pct = (off_target / total * 100) if total > 0 else 0.0
        validated_pct = 100.0 - off_target_pct

        # Format summary stats
        content = FAMILY_VALIDATION_SECTION_TEMPLATE.format(
            validated_pct=validated_pct,
            off_target_count=off_target,
            off_target_pct=off_target_pct,
            target_family=(
                self.summary.target_family or "Inferred from ANI matrix"
            ),
            external_families_table="",
        )

        # Histogram of family_bitscore_ratio
        if "family_bitscore_ratio" in self.df.columns:
            ratios = self.df["family_bitscore_ratio"].drop_nulls().to_list()
            if ratios:
                fig = go.Figure(data=[go.Histogram(
                    x=ratios,
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

    def _build_ani_section(self) -> str:
        """Build the ANI matrix heatmap tab section."""
        if self.ani_matrix is None or len(self.ani_matrix) == 0:
            content = EMPTY_SECTION_TEMPLATE.format(
                message=ANI_NOT_PROVIDED_MESSAGE
            )
        else:
            # Create genome labels map
            matrix_cols = self.ani_matrix.columns
            if "genome" in matrix_cols:
                genome_accessions = self.ani_matrix["genome"].to_list()
            else:
                genome_accessions = list(matrix_cols)

            genome_labels_map = self._get_genome_labels_map(genome_accessions)

            # Build heatmap and compute statistics using extracted components
            heatmap_fig, stats, clustering_succeeded = build_ani_heatmap(
                self.ani_matrix,
                genome_labels_map,
            )

            # Register plot
            ani_id = "plot-ani-heatmap"
            self._register_plot(ani_id, heatmap_fig)

            # Generate statistics cards
            stats_html = build_ani_stats_cards(stats)

            content = stats_html + PLOT_CONTAINER_TEMPLATE.format(
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

        return TAB_SECTION_TEMPLATE.format(
            tab_id="ani",
            active_class="",
            section_title="Reference ANI",
            content=content,
        )

    def _build_aai_section(self) -> str:
        """Build the AAI matrix heatmap tab section for genus-level comparisons."""
        if self.aai_matrix is None or len(self.aai_matrix) == 0:
            content = EMPTY_SECTION_TEMPLATE.format(
                message=(
                    "AAI matrix not provided. "
                    "Use --aai option with AAI matrix CSV to enable genus-level heatmap. "
                    "Generate with: metadarkmatter aai compute"
                )
            )
        else:
            # Create genome labels map
            matrix_cols = self.aai_matrix.columns
            if "genome" in matrix_cols:
                genome_accessions = self.aai_matrix["genome"].to_list()
            else:
                genome_accessions = list(matrix_cols)

            genome_labels_map = self._get_genome_labels_map(genome_accessions)

            # Build heatmap and compute statistics using extracted components
            heatmap_fig, stats, clustering_succeeded = build_aai_heatmap(
                self.aai_matrix,
                genome_labels_map,
            )

            # Register plot
            aai_id = "plot-aai-heatmap"
            self._register_plot(aai_id, heatmap_fig)

            # Generate statistics cards
            stats_html = build_aai_stats_cards(stats)

            content = stats_html + PLOT_CONTAINER_TEMPLATE.format(
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

        return TAB_SECTION_TEMPLATE.format(
            tab_id="aai",
            active_class="",
            section_title="Reference AAI",
            content=content,
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
            identity_gap = row.get("identity_gap")
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
            uncertainty_type = row.get("uncertainty_type", "-")
            inferred_unc = row.get("inferred_uncertainty")
            inferred_unc_str = "-" if inferred_unc is None else f"{inferred_unc:.1f}"
            discovery = row.get("discovery_score")
            discovery_str = "-" if discovery is None else f"{discovery:.1f}"
            id_conf = row.get("identity_confidence")
            id_conf_str = "-" if id_conf is None else f"{id_conf:.1f}"
            pl_conf = row.get("placement_confidence")
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

            # Build tree data JSON for D3.js visualization
            tree_data = {
                "newick": final_newick,
                "annotations": annotations,
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

    def _build_methods_section(self) -> str:
        """Build the Methods tab with comprehensive documentation of calculations."""
        return TAB_SECTION_TEMPLATE.format(
            tab_id="methods",
            active_class="",
            section_title="Methods",
            content=METHODS_SECTION_TEMPLATE,
        )

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
