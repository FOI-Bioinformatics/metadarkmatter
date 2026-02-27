"""
Data models for novel diversity clustering.

Defines the NovelCluster model representing a putative novel taxon cluster
inferred from grouped novel reads sharing similar characteristics.
"""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field, computed_field


class GenusDistance(BaseModel):
    """Distance from a novel cluster to a reference genus."""

    genus: str = Field(description="Genus name")
    representative_genome: str = Field(description="Closest genome in this genus")
    estimated_ani: float = Field(description="Estimated ANI from cluster to this genus")
    num_genomes_in_genus: int = Field(description="Reference genomes in this genus")

    model_config = {"frozen": True}


class PhylogeneticNeighborhood(BaseModel):
    """Phylogenetic neighborhood profile for a novel cluster."""

    cluster_id: str = Field(description="Cluster this neighborhood belongs to")
    nearest_genera: list[GenusDistance] = Field(
        description="Nearest genera sorted by ANI (descending), top 5"
    )
    placement_support: float = Field(
        ge=0, le=100, description="Placement support score (0-100)"
    )
    isolation_score: float = Field(
        ge=0, description="ANI gap between nearest and second-nearest genus"
    )
    neighborhood_density: int = Field(
        ge=0, description="Number of genera within genus_boundary + 5% ANI"
    )
    phylogenetic_context: str = Field(
        description="One-line human-readable placement text"
    )
    genus_boundary_ani: float | None = Field(
        default=None, description="Detected genus boundary ANI from GMM"
    )

    model_config = {"frozen": True}


class NovelCluster(BaseModel):
    """
    Represents a putative novel taxon cluster.

    Clusters are formed by grouping novel reads that share the same:
    - Best match genome (closest reference)
    - Novelty band (5% intervals of divergence)
    - Taxonomic call (Novel Species or Novel Genus)

    Attributes:
        cluster_id: Unique identifier (e.g., NSP_001, NGN_001)
        taxonomic_call: "Novel Species" or "Novel Genus"
        nearest_genome: Accession of the closest reference genome
        nearest_species: Species name of closest reference
        nearest_genus: Genus name of closest reference
        nearest_family: Family name of closest reference
        read_count: Number of reads in this cluster
        mean_novelty_index: Average novelty index of reads
        novelty_min: Minimum novelty index in cluster
        novelty_max: Maximum novelty index in cluster
        mean_placement_uncertainty: Average placement uncertainty
        mean_discovery_score: Average discovery score (if enhanced scoring)
        suggested_name: GTDB-style provisional name
        confidence: Cluster confidence rating
        phylogenetic_context: Inferred taxonomic placement description
        read_ids: List of read identifiers in this cluster
    """

    cluster_id: str = Field(
        description="Unique cluster identifier (NSP_001 for Novel Species, NGN_001 for Novel Genus)"
    )
    taxonomic_call: Literal["Novel Species", "Novel Genus"] = Field(
        description="Classification type for reads in this cluster"
    )
    nearest_genome: str = Field(
        description="Accession of the closest reference genome"
    )
    nearest_species: str = Field(
        description="Species name of closest reference genome"
    )
    nearest_genus: str = Field(
        description="Genus name of closest reference genome"
    )
    nearest_family: str = Field(
        description="Family name of closest reference genome"
    )
    novelty_band: int = Field(
        description="Novelty band (5-10, 10-15, 15-20, 20-25)"
    )
    read_count: int = Field(
        ge=1,
        description="Number of reads in this cluster"
    )
    mean_novelty_index: float = Field(
        ge=0,
        le=100,
        description="Average novelty index of reads in cluster"
    )
    novelty_min: float = Field(
        ge=0,
        le=100,
        description="Minimum novelty index in cluster"
    )
    novelty_max: float = Field(
        ge=0,
        le=100,
        description="Maximum novelty index in cluster"
    )
    mean_placement_uncertainty: float = Field(
        ge=0,
        le=100,
        description="Average placement uncertainty of reads"
    )
    mean_discovery_score: float | None = Field(
        default=None,
        ge=0,
        le=100,
        description="Average discovery score (if enhanced scoring enabled)"
    )
    suggested_name: str = Field(
        description="GTDB-style provisional name for the novel taxon"
    )
    confidence: Literal["High", "Medium", "Low"] = Field(
        description="Cluster confidence rating based on read count, uncertainty, and discovery score"
    )
    phylogenetic_context: str = Field(
        description="Inferred phylogenetic placement description"
    )
    read_ids: list[str] = Field(
        default_factory=list,
        description="List of read identifiers in this cluster"
    )
    mean_bayesian_confidence: float | None = Field(
        default=None,
        ge=0,
        le=100,
        description="Mean per-read Bayesian confidence_score (0-100) for reads in this cluster"
    )
    contributing_genomes: list[str] = Field(
        default_factory=list,
        description="Reference genomes that contributed reads to this cluster (after neighborhood merging)"
    )
    neighborhood: PhylogeneticNeighborhood | None = Field(
        default=None,
        description="Phylogenetic neighborhood profile (computed post-clustering)"
    )

    model_config = {"frozen": True}

    @computed_field
    @property
    def novelty_range(self) -> str:
        """Format novelty range as string."""
        return f"{self.novelty_min:.1f}-{self.novelty_max:.1f}"

    @computed_field
    @property
    def novelty_band_label(self) -> str:
        """Human-readable novelty band label."""
        band_labels = {
            5: "5-10% (Recently diverged)",
            10: "10-15% (Moderately divergent)",
            15: "15-20% (Highly divergent)",
            20: "20-25% (Novel genus candidate)",
        }
        return band_labels.get(self.novelty_band, f"{self.novelty_band}-{self.novelty_band + 5}%")

    @computed_field
    @property
    def estimated_ani(self) -> float:
        """Estimated ANI to nearest reference (100 - mean_novelty_index)."""
        return round(100.0 - self.mean_novelty_index, 1)

    @computed_field
    @property
    def closest_known_taxon(self) -> str:
        """
        Returns the closest known taxon at the appropriate level.

        For Novel Species: returns closest species name
        For Novel Genus: returns closest genus name
        """
        if self.taxonomic_call == "Novel Species":
            return self.nearest_species
        else:
            return self.nearest_genus

    @computed_field
    @property
    def closest_taxon_with_similarity(self) -> str:
        """
        Returns closest known taxon with estimated ANI similarity.

        For Novel Species: "Francisella tularensis (~92% ANI)"
        For Novel Genus: "Francisella (~78% ANI)"
        """
        taxon = self.closest_known_taxon
        return f"{taxon} (~{self.estimated_ani:.0f}% ANI)"

    def to_summary_dict(self) -> dict:
        """Convert to dictionary for summary output (excludes read_ids)."""
        result = {
            "cluster_id": self.cluster_id,
            "taxonomic_call": self.taxonomic_call,
            "read_count": self.read_count,
            "nearest_species": self.nearest_species,
            "nearest_genus": self.nearest_genus,
            "nearest_family": self.nearest_family,
            "closest_known_taxon": self.closest_known_taxon,
            "estimated_ani": self.estimated_ani,
            "mean_novelty": round(self.mean_novelty_index, 2),
            "novelty_range": self.novelty_range,
            "mean_uncertainty": round(self.mean_placement_uncertainty, 2),
            "mean_discovery_score": (
                round(self.mean_discovery_score, 1)
                if self.mean_discovery_score is not None
                else None
            ),
            "suggested_name": self.suggested_name,
            "confidence": self.confidence,
            "phylogenetic_context": self.phylogenetic_context,
        }
        if self.mean_bayesian_confidence is not None:
            result["mean_bayesian_confidence"] = round(self.mean_bayesian_confidence, 1)
        if self.contributing_genomes:
            result["contributing_genomes_count"] = len(self.contributing_genomes)
        return result


class NovelDiversitySummary(BaseModel):
    """
    Summary of novel diversity analysis for a sample.

    Provides aggregate statistics across all novel clusters,
    useful for report generation and quick assessment.

    Attributes:
        total_novel_reads: Total number of novel reads analyzed
        novel_species_reads: Number of Novel Species reads
        novel_genus_reads: Number of Novel Genus reads
        total_clusters: Total number of clusters formed
        novel_species_clusters: Number of Novel Species clusters
        novel_genus_clusters: Number of Novel Genus clusters
        high_confidence_clusters: Clusters with "High" confidence
        medium_confidence_clusters: Clusters with "Medium" confidence
        low_confidence_clusters: Clusters with "Low" confidence
        genera_with_novel_species: Unique genera with novel species
        families_with_novel_genera: Unique families with novel genera
        mean_cluster_size: Average reads per cluster
        largest_cluster_size: Size of largest cluster
    """

    total_novel_reads: int = Field(default=0)
    novel_species_reads: int = Field(default=0)
    novel_genus_reads: int = Field(default=0)
    total_clusters: int = Field(default=0)
    novel_species_clusters: int = Field(default=0)
    novel_genus_clusters: int = Field(default=0)
    high_confidence_clusters: int = Field(default=0)
    medium_confidence_clusters: int = Field(default=0)
    low_confidence_clusters: int = Field(default=0)
    genera_with_novel_species: int = Field(default=0)
    families_with_novel_genera: int = Field(default=0)
    mean_cluster_size: float = Field(default=0.0)
    largest_cluster_size: int = Field(default=0)

    @computed_field
    @property
    def novel_species_pct(self) -> float:
        """Percentage of novel reads that are Novel Species."""
        if self.total_novel_reads == 0:
            return 0.0
        return self.novel_species_reads / self.total_novel_reads * 100

    @computed_field
    @property
    def novel_genus_pct(self) -> float:
        """Percentage of novel reads that are Novel Genus."""
        if self.total_novel_reads == 0:
            return 0.0
        return self.novel_genus_reads / self.total_novel_reads * 100
