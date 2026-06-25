"""Data models and config for HTML report generation.

Holds the plain-data pieces shared across the report sections: the
scoring/threshold summary (:class:`TaxonomicSummary`), the per-group
metric bundles, the :class:`ReportConfig`, and the ``_safe_float`` helper.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, cast

from metadarkmatter.visualization.plots.base import PlotConfig, ThresholdConfig


def _safe_float(value: object) -> float | None:
    """Convert a value to float, returning None for missing/unconvertible values."""
    if value is None:
        return None
    try:
        # value is `object` (often a polars Series.mean() with a wide stub
        # type); the try/except guards genuinely non-numeric inputs.
        return float(cast(Any, value))
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
    # Asset mode for embedded JS dependencies (Plotly, D3). 'offline' inlines
    # the JS from the Python package and any vendored assets so the report is
    # fully self-contained; 'cdn' uses remote script tags (smaller file).
    report_mode: str = "offline"  # 'offline' or 'cdn'

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


@dataclass(frozen=True)
class _EnhancedScoringMetrics:
    """Mean enhanced-scoring / discovery metrics for a classification set."""

    mean_inferred_uncertainty: float | None = None
    mean_alignment_quality: float | None = None
    mean_identity_confidence: float | None = None
    mean_placement_confidence: float | None = None
    mean_discovery_score: float | None = None
    novel_with_discovery_score: int = 0
    high_priority_discoveries: int = 0


@dataclass(frozen=True)
class _BayesianMetrics:
    """Mean posterior entropy and confidence / agreement counts."""

    mean_posterior_entropy: float = 0.0
    high_confidence_count: int = 0
    high_confidence_pct: float = 0.0
    boundary_count: int = 0
    boundary_pct: float = 0.0
    map_agreement_count: int = 0
    map_agreement_pct: float = 0.0
