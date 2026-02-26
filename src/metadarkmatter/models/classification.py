"""
Pydantic models for read classification results.

These models represent the output of the ANI-weighted placement uncertainty
algorithm applied to environmental DNA metagenomic samples.
"""

from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field, computed_field


class DiversityStatus(str, Enum):
    """
    High-level diversity classification for summarizing results.

    Groups detailed taxonomic calls into categories for easier
    interpretation and communication of results.

    Categories:
        KNOWN: Confident match to known diversity (Known Species)
        NOVEL: Confident novel diversity (Novel Species, Novel Genus)
        UNCERTAIN: Cannot confidently classify as known or novel
        OFF_TARGET: Reads with better matches outside the target family
    """

    KNOWN = "Known"
    NOVEL = "Novel"
    UNCERTAIN = "Uncertain"
    OFF_TARGET = "Off-target"


# Mapping from TaxonomicCall to DiversityStatus
TAXONOMIC_TO_DIVERSITY: dict[str, str] = {
    "Known Species": "Known",
    "Novel Species": "Novel",
    "Novel Genus": "Novel",
    "Species Boundary": "Uncertain",
    "Ambiguous": "Uncertain",
    "Ambiguous Within Genus": "Uncertain",
    "Conserved Region": "Uncertain",
    "Unclassified": "Uncertain",
    "Off-target": "Off-target",
}


class TaxonomicCall(str, Enum):
    """
    Taxonomic classification categories for environmental DNA reads.

    Based on novelty index and placement uncertainty thresholds from the
    competitive read recruitment methodology.

    Categories:
        KNOWN_SPECIES: High identity (N < 5%), low ambiguity (U < 2%)
        NOVEL_SPECIES: Moderate divergence (5-20%), low ambiguity (U < 2%)
        NOVEL_GENUS: High divergence (20-25%), low ambiguity (U < 2%)
        SPECIES_BOUNDARY: Moderate placement ambiguity (2% <= U < 5%) indicating
            the read matches multiple closely-related species at the species
            boundary zone (95-98% ANI between competing genomes)
        AMBIGUOUS: High placement ambiguity (U >= 5%) within genus, or identity
            gap ambiguity. Indicates conserved region shared by multiple species
            within the same genus.
        AMBIGUOUS_WITHIN_GENUS: Read hits multiple species within the same genus
            with similar bitscores (90% of top) but low ANI between them,
            used for Novel Genus candidates with genus-level competing hits
        CONSERVED_REGION: Very high placement ambiguity (U >= 5%) AND hits
            span multiple genera, indicating highly conserved gene (e.g., 16S)
        UNCLASSIFIED: Does not fit clear biological categories - may include
            strain variants or very high divergence (N > 25%)
    """

    KNOWN_SPECIES = "Known Species"
    NOVEL_SPECIES = "Novel Species"
    NOVEL_GENUS = "Novel Genus"
    SPECIES_BOUNDARY = "Species Boundary"
    AMBIGUOUS = "Ambiguous"
    AMBIGUOUS_WITHIN_GENUS = "Ambiguous Within Genus"
    CONSERVED_REGION = "Conserved Region"
    UNCLASSIFIED = "Unclassified"
    OFF_TARGET = "Off-target"

    @property
    def diversity_status(self) -> DiversityStatus:
        """Get the high-level diversity status for this taxonomic call."""
        return DiversityStatus(TAXONOMIC_TO_DIVERSITY[self.value])


class ReadClassification(BaseModel):
    """
    Classification result for a single metagenomic read.

    Contains the ANI-weighted placement metrics and final taxonomic call
    for detecting novel bacterial diversity in environmental samples.

    Attributes:
        read_id: Original read identifier from FASTQ/BLAST
        best_match_genome: Genome with highest BLAST bitscore hit
        top_hit_identity: Percent identity of best BLAST hit (0-100)
        novelty_index: 100 - top_hit_identity; measures divergence from reference
        placement_uncertainty: 100 - ANI(top_hit, secondary_hit); measures ambiguity
        num_ambiguous_hits: Number of hits within 95% of top bitscore
        taxonomic_call: Final classification based on thresholds
    """

    read_id: str = Field(description="Original read identifier")
    best_match_genome: str = Field(description="Genome with highest bitscore hit")
    top_hit_identity: float = Field(
        ge=0,
        le=100,
        description="Percent identity of best BLAST hit",
    )
    novelty_index: float = Field(
        ge=0,
        le=100,
        description="100 - top_hit_identity; measures divergence from reference",
    )
    placement_uncertainty: float = Field(
        ge=0,
        le=100,
        description="100 - ANI(top_hit, secondary_hit); measures placement ambiguity",
    )
    num_ambiguous_hits: int = Field(
        ge=0,
        description="Number of hits within 95% of top bitscore",
    )
    second_hit_identity: float | None = Field(
        default=None,
        ge=0,
        le=100,
        description="Percent identity of second-best hit to different genome (if exists)",
    )
    identity_gap: float | None = Field(
        default=None,
        ge=0,
        le=100,
        description="Gap between best and second-best hit identity (top - second)",
    )
    genus_uncertainty: float | None = Field(
        default=None,
        ge=0,
        le=100,
        description=(
            "Uncertainty based on genus-level hits (90% bitscore threshold). "
            "Higher values indicate competing hits from different species in "
            "the same genus with low ANI between them."
        ),
    )
    num_genus_hits: int | None = Field(
        default=None,
        ge=0,
        description="Number of hits within genus-level bitscore threshold (90%)",
    )
    confidence_score: float | None = Field(
        default=None,
        ge=0,
        le=100,
        description=(
            "Overall confidence score (0-100) integrating multiple quality factors: "
            "margin from classification threshold boundaries, alignment quality, "
            "bitscore gap to secondary hits, and number of competing placements. "
            "Scores below 50 indicate borderline classifications."
        ),
    )
    # Enhanced scoring fields (optional, enabled via --enhanced-scoring)
    inferred_uncertainty: float | None = Field(
        default=None,
        ge=0,
        le=100,
        description=(
            "Inferred placement uncertainty for single-hit reads based on "
            "novelty level. Only set when num_ambiguous_hits <= 1."
        ),
    )
    uncertainty_type: str | None = Field(
        default=None,
        description="Source of uncertainty: 'measured' (from ANI) or 'inferred' (from novelty)",
    )
    alignment_quality: float | None = Field(
        default=None,
        ge=0,
        le=100,
        description="Alignment quality score based on mismatch, gap, coverage, and evalue",
    )
    identity_confidence: float | None = Field(
        default=None,
        ge=0,
        le=100,
        description="Confidence in the identity measurement itself",
    )
    placement_confidence: float | None = Field(
        default=None,
        ge=0,
        le=100,
        description="Confidence in genome assignment, accounting for uncertainty source",
    )
    discovery_score: float | None = Field(
        default=None,
        ge=0,
        le=100,
        description="Priority score for novel discoveries (null for non-novel)",
    )
    taxonomic_call: TaxonomicCall = Field(description="Final classification")

    model_config = {"frozen": True}

    @computed_field
    @property
    def diversity_status(self) -> str:
        """Get the high-level diversity status (Known/Novel/Uncertain)."""
        return self.taxonomic_call.diversity_status.value

    @computed_field
    @property
    def is_novel(self) -> bool:
        """Check if this read represents novel diversity (dark matter)."""
        return self.taxonomic_call in (
            TaxonomicCall.NOVEL_SPECIES,
            TaxonomicCall.NOVEL_GENUS,
        )

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for DataFrame creation."""
        result = {
            "read_id": self.read_id,
            "best_match_genome": self.best_match_genome,
            "top_hit_identity": self.top_hit_identity,
            "novelty_index": self.novelty_index,
            "placement_uncertainty": self.placement_uncertainty,
            "num_ambiguous_hits": self.num_ambiguous_hits,
            "second_hit_identity": self.second_hit_identity,
            "identity_gap": self.identity_gap,
            "genus_uncertainty": self.genus_uncertainty,
            "num_genus_hits": self.num_genus_hits,
            "confidence_score": self.confidence_score,
            "taxonomic_call": self.taxonomic_call.value,
            "diversity_status": self.diversity_status,
            "is_novel": self.is_novel,
        }
        # Add enhanced scoring fields if they are set
        if self.inferred_uncertainty is not None:
            result["inferred_uncertainty"] = self.inferred_uncertainty
        if self.uncertainty_type is not None:
            result["uncertainty_type"] = self.uncertainty_type
        if self.alignment_quality is not None:
            result["alignment_quality"] = self.alignment_quality
        if self.identity_confidence is not None:
            result["identity_confidence"] = self.identity_confidence
        if self.placement_confidence is not None:
            result["placement_confidence"] = self.placement_confidence
        if self.discovery_score is not None:
            result["discovery_score"] = self.discovery_score
        return result


class TaxonomicSummary(BaseModel):
    """
    Summary statistics for a classified environmental sample.

    Aggregates read-level classifications to provide sample-level insights
    into microbial diversity and novel taxa detection.

    Attributes:
        total_reads: Total reads classified in this sample
        known_species: Reads classified as known species
        novel_species: Reads classified as novel species (dark matter)
        novel_genus: Reads classified as novel genus (dark matter)
        species_boundary: Reads at species boundary (U 2-5%, closely related species)
        ambiguous: Reads with high placement ambiguity within genus (U >= 5%)
        ambiguous_within_genus: Novel Genus candidates with genus-level competing hits
        conserved_regions: Reads matching conserved genes across genera
        unclassified: Reads not fitting clear biological categories
        diversity_known: Reads with confident known diversity
        diversity_novel: Reads with confident novel diversity
        diversity_uncertain: Reads with uncertain diversity status
        mean_novelty_index: Mean novelty index across all reads
        mean_placement_uncertainty: Mean placement uncertainty
        genome_hit_counts: Read counts per reference genome (top 50)
    """

    total_reads: int = Field(description="Total reads classified")
    known_species: int = Field(description="Reads classified as known species")
    novel_species: int = Field(description="Reads classified as novel species")
    novel_genus: int = Field(description="Reads classified as novel genus")
    species_boundary: int = Field(
        default=0,
        description="Reads at species boundary (U 2-5%, closely related species)",
    )
    ambiguous: int = Field(
        default=0,
        description="Reads with high placement ambiguity within genus (U >= 5%)",
    )
    ambiguous_within_genus: int = Field(
        default=0,
        description="Novel Genus candidates hitting multiple species within same genus",
    )
    conserved_regions: int = Field(description="Reads matching conserved genes across genera")
    unclassified: int = Field(default=0, description="Reads not fitting clear categories")
    off_target: int = Field(
        default=0,
        description="Reads classified as off-target (better match outside target family)",
    )
    # High-level diversity grouping
    diversity_known: int = Field(
        default=0,
        description="Reads with confident known diversity (Known Species)",
    )
    diversity_novel: int = Field(
        default=0,
        description="Reads with confident novel diversity (Novel Species + Novel Genus)",
    )
    diversity_uncertain: int = Field(
        default=0,
        description="Reads with uncertain diversity status",
    )
    mean_novelty_index: float = Field(description="Mean novelty index across all reads")
    mean_placement_uncertainty: float = Field(description="Mean placement uncertainty")
    genome_hit_counts: dict[str, int] = Field(
        default_factory=dict,
        description="Read counts per reference genome (top 50)",
    )
    species_hit_counts: dict[str, int] = Field(
        default_factory=dict,
        description="Read counts per species (requires metadata)",
    )

    @computed_field
    @property
    def diversity_known_pct(self) -> float:
        """Percentage of reads with known diversity."""
        if self.total_reads == 0:
            return 0.0
        return self.diversity_known / self.total_reads * 100

    @computed_field
    @property
    def diversity_novel_pct(self) -> float:
        """Percentage of reads with novel diversity."""
        if self.total_reads == 0:
            return 0.0
        return self.diversity_novel / self.total_reads * 100

    @computed_field
    @property
    def diversity_uncertain_pct(self) -> float:
        """Percentage of reads with uncertain diversity."""
        if self.total_reads == 0:
            return 0.0
        return self.diversity_uncertain / self.total_reads * 100

    @computed_field
    @property
    def known_species_pct(self) -> float:
        """Percentage of reads classified as known species."""
        if self.total_reads == 0:
            return 0.0
        return 100.0 * self.known_species / self.total_reads

    @computed_field
    @property
    def novel_species_pct(self) -> float:
        """Percentage of reads classified as novel species."""
        if self.total_reads == 0:
            return 0.0
        return 100.0 * self.novel_species / self.total_reads

    @computed_field
    @property
    def novel_genus_pct(self) -> float:
        """Percentage of reads classified as novel genus."""
        if self.total_reads == 0:
            return 0.0
        return 100.0 * self.novel_genus / self.total_reads

    @computed_field
    @property
    def novel_diversity_pct(self) -> float:
        """Total percentage of novel diversity detected in sample."""
        return self.novel_species_pct + self.novel_genus_pct

    @computed_field
    @property
    def species_count(self) -> int:
        """Number of unique species detected (requires metadata)."""
        return len(self.species_hit_counts)

    def to_json(self, path: Path) -> None:
        """Write summary to JSON file."""
        path.write_text(self.model_dump_json(indent=2))

    @classmethod
    def from_json(cls, path: Path) -> TaxonomicSummary:
        """Load summary from JSON file."""
        return cls.model_validate_json(path.read_text())
