"""
Pydantic configuration models for metadarkmatter.

These models define configuration for the ANI-weighted placement uncertainty
algorithm, including scoring thresholds, file paths, and external tool settings.
Configuration can be loaded from YAML files, environment variables, or CLI arguments.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Literal, Self

from pydantic import BaseModel, Field, model_validator

logger = logging.getLogger(__name__)


# Default uniform priors over 6 Bayesian categories
_DEFAULT_PRIORS: dict[str, float] = {
    "known_species": 1.0 / 6,
    "novel_species": 1.0 / 6,
    "novel_genus": 1.0 / 6,
    "species_boundary": 1.0 / 6,
    "ambiguous": 1.0 / 6,
    "unclassified": 1.0 / 6,
}


class BayesianConfig(BaseModel):
    """Configuration for the Bayesian classification engine."""

    priors: dict[str, float] = Field(
        default_factory=lambda: dict(_DEFAULT_PRIORS),
        description="Prior probabilities for 6 Bayesian categories (sum to 1.0).",
    )
    identity_gap_boost: float = Field(
        default=2.0,
        ge=1.0,
        description=(
            "Prior modulation factor for Ambiguous when identity_gap < threshold "
            "and num_hits > 1."
        ),
    )
    single_hit_boost: float = Field(
        default=1.5,
        ge=1.0,
        description=(
            "Prior modulation factor for Ambiguous when single-hit reads have "
            "novelty in the novel-species range or higher."
        ),
    )
    entropy_threshold: float = Field(
        default=1.5,
        ge=0.0,
        description=(
            "Posterior entropy threshold for low-confidence flagging. "
            "Reads with entropy above this value are flagged as low confidence."
        ),
    )

    @model_validator(mode="after")
    def validate_priors_sum(self) -> Self:
        """Priors must sum to approximately 1.0."""
        total = sum(self.priors.values())
        if abs(total - 1.0) > 0.01:
            msg = f"Bayesian priors must sum to 1.0, got {total:.4f}"
            raise ValueError(msg)
        return self

    model_config = {"frozen": True}


class ScoringConfig(BaseModel):
    """
    Configuration for ANI-weighted placement scoring thresholds.

    These thresholds define the boundaries between different taxonomic
    classifications based on novelty index and placement uncertainty metrics.

    Literature Support:
        - 95-96% ANI is the widely accepted prokaryotic species boundary
          (Jain et al. 2018, Nature Communications; Goris et al. 2007, PNAS)
        - >96% identity = same species (stricter boundary per Jain et al.)
        - 70-80% ANI is the approximate genus boundary range

    Novelty Thresholds (N = 100 - pident, read-to-reference identity):
        - N < 4%: Known Species (pident >96%, at or above species boundary)
        - 4% <= N <= 20%: Novel Species (pident 80-96%)
        - 20% <= N <= 25%: Novel Genus (pident 75-80%)

    Note: Read-level BLAST identity can be 10-20% lower than genome-level ANI
    due to reads hitting divergent regions. The 20% novelty threshold for genus
    accounts for this discrepancy.

    Uncertainty Thresholds (U = 100 - ANI between competing genomes):
        - U < 1.5%: Confident placement (competitors >98.5% ANI, same species)
        - 1.5% <= U < 5%: Ambiguous (competitors 95-98.5% ANI, species boundary)
        - U >= 5%: Conserved region (competitors <95% ANI, different species)

    Enhanced scoring (always enabled):
        - Inferred uncertainty for single-hit reads (based on novelty level)
        - Alignment quality, identity confidence, placement confidence
        - Discovery score for prioritizing novel findings
        - Linear coverage weighting to reduce conserved-domain bias

    Protein-Level Classification:
        When alignment_mode is "protein", wider thresholds are used because
        amino acid sequences are more conserved than nucleotide sequences:
        - Protein N < 10%: Known Species (protein identity > 90%)
        - Protein 10% <= N <= 25%: Novel Species (protein identity 75-90%)
        - Protein 25% <= N <= 40%: Novel Genus (protein identity 60-75%)
    """

    # Alignment mode - determines which threshold set to use
    alignment_mode: Literal["nucleotide", "protein"] = Field(
        default="nucleotide",
        description=(
            "Alignment mode: 'nucleotide' for BLASTN results, 'protein' for BLASTX results. "
            "Protein mode uses wider novelty thresholds (10/25/40%) because amino acid "
            "sequences are more conserved than nucleotide sequences."
        ),
    )

    bitscore_threshold_pct: float = Field(
        default=95.0,
        ge=0,
        le=100,
        description="Percentage of top bitscore for ambiguous hit detection",
    )

    # Genus-aware bitscore threshold for detecting within-genus ambiguity
    # A lower threshold (90%) captures hits from other species in the same genus
    # that would otherwise be missed by the stricter 95% threshold
    genus_bitscore_threshold_pct: float = Field(
        default=90.0,
        ge=0,
        le=100,
        description=(
            "Percentage of top bitscore for genus-level ambiguity detection. "
            "When hits from multiple species in the same genus are within this "
            "threshold, the read is flagged as ambiguous within genus rather than "
            "novel. Default 90% captures same-genus hits missed by the 95% threshold."
        ),
    )

    # Novelty index thresholds
    # Based on: 95% pident = species boundary (Jain et al. 2018)
    # N = 100 - pident, so N < 5% means pident > 95% (same species)
    novelty_known_max: float = Field(
        default=4.0,
        ge=0,
        description="Maximum novelty index for known species (pident >96%)",
    )
    novelty_novel_species_min: float = Field(
        default=4.0,
        ge=0,
        description="Minimum novelty index for novel species (pident <=96%)",
    )
    novelty_novel_species_max: float = Field(
        default=20.0,
        ge=0,
        description="Maximum novelty index for novel species (pident >= 80%)",
    )
    novelty_novel_genus_min: float = Field(
        default=20.0,
        ge=0,
        description="Minimum novelty index for novel genus (pident < 80%)",
    )
    novelty_novel_genus_max: float = Field(
        default=25.0,
        ge=0,
        description="Maximum novelty index for novel genus",
    )

    # Placement uncertainty thresholds
    # Based on: 95-96% ANI = species boundary (Jain et al. 2018)
    # U < 2% means competing genomes share >98% ANI (same species, confident)
    # 2% <= U < 5% means 95-98% ANI (species boundary zone, ambiguous)
    # U >= 5% means <95% ANI (different species, conserved gene)
    uncertainty_known_max: float = Field(
        default=1.5,
        ge=0,
        description="Maximum placement uncertainty for known species (ANI >98.5%)",
    )
    uncertainty_novel_species_max: float = Field(
        default=1.5,
        ge=0,
        description="Maximum placement uncertainty for novel species (ANI >98.5%)",
    )
    uncertainty_novel_genus_max: float = Field(
        default=1.5,
        ge=0,
        description="Maximum placement uncertainty for novel genus (ANI >98.5%)",
    )
    uncertainty_conserved_min: float = Field(
        default=5.0,
        ge=0,
        description="Minimum uncertainty for conserved regions (ANI <95%)",
    )

    # Genus-level uncertainty threshold
    # When genus_uncertainty exceeds this value, the read hits multiple species
    # within the genus with low ANI, indicating ambiguous placement within genus
    genus_uncertainty_ambiguous_min: float = Field(
        default=10.0,
        ge=0,
        description=(
            "Minimum genus uncertainty to flag as ambiguous within genus. "
            "Genus uncertainty is calculated from ANI between competing species "
            "hits within the genus-level bitscore threshold. Values >= 10% indicate "
            "hits to species with <90% ANI (different species in same genus)."
        ),
    )

    # Identity gap threshold for ambiguity detection
    # When the gap between best hit and second-best hit identity is small,
    # the read matches multiple genomes with similar quality, indicating
    # ambiguous placement regardless of ANI-based uncertainty.
    identity_gap_ambiguous_max: float = Field(
        default=2.0,
        ge=0,
        le=100,
        description=(
            "Maximum identity gap (best - second best) before flagging as ambiguous. "
            "Gap < 2% means competing hits have similar identity, indicating ambiguity."
        ),
    )

    # Confidence threshold for quality flagging
    # Reads with confidence_score below this threshold are flagged as low_confidence
    # without affecting the classification category
    confidence_threshold: float = Field(
        default=50.0,
        ge=0,
        le=100,
        description=(
            "Minimum confidence score for high-confidence classification. "
            "Reads below this threshold are flagged but classification unchanged."
        ),
    )

    # Alignment quality filters (GTDB-compatible)
    # GTDB requires alignment fraction (AF) >= 50% for species assignment
    min_alignment_length: int = Field(
        default=100,
        ge=0,
        description="Minimum alignment length in bp (0 = no filter)",
    )
    min_alignment_fraction: float = Field(
        default=0.3,
        ge=0.0,
        le=1.0,
        description=(
            "Minimum fraction of read aligned (like GTDB's AF). "
            "Set to 0.5 for GTDB-compatible filtering. Default 0.3 filters "
            "short conserved domain hits."
        ),
    )

    # Coverage weighting for hit selection
    # Prioritizes longer alignments over short conserved domains
    coverage_weight_mode: Literal["none", "linear", "log", "sigmoid"] = Field(
        default="linear",
        description=(
            "Coverage weighting mode for hit selection. "
            "'linear' increases weight linearly with coverage (default). "
            "'none' uses raw bitscore only. "
            "'log' provides diminishing returns for increasing coverage. "
            "'sigmoid' applies sharp threshold around 60% coverage."
        ),
    )

    coverage_weight_strength: float = Field(
        default=0.5,
        ge=0.0,
        le=1.0,
        description=(
            "Coverage weight strength (0.0-1.0). "
            "Controls magnitude of coverage effect on scoring. "
            "Weight range = [1-strength, 1+strength]. "
            "At strength=0.5: weights range from 0.5x to 1.5x. "
            "Has no effect when coverage_weight_mode='none'."
        ),
    )

    # Uncertainty calculation mode
    uncertainty_mode: Literal["max", "second"] = Field(
        default="second",
        description=(
            "Mode for calculating placement uncertainty. "
            "'max' uses maximum ANI to any competing genome. "
            "'second' uses ANI to the second-best genome only (default)."
        ),
    )

    # Single-hit classification options
    # Single-hit reads (~70% of environmental data) get inferred uncertainty
    # and are gated by threshold to reduce novel species overestimation.
    single_hit_uncertainty_threshold: float = Field(
        default=10.0,
        ge=0,
        le=100,
        description=(
            "Inferred uncertainty threshold above which single-hit reads are flagged "
            "as Ambiguous. Inferred uncertainty is always computed for single-hit reads. "
            "Default 10% means reads with novelty ~7% or higher become Ambiguous. "
            "Lower values are stricter, higher values are more permissive."
        ),
    )

    # AAI (Average Amino Acid Identity) thresholds for genus-level classification
    # Based on Riesco & Trujillo 2024: 58-65% AAI genus boundary
    # AAI is more reliable than ANI for genus-level decisions because
    # ANI becomes unreliable below approximately 80% identity.
    aai_genus_boundary_low: float = Field(
        default=58.0,
        ge=0.0,
        le=100.0,
        description=(
            "Lower AAI boundary for genus classification. "
            "AAI below this value indicates different genus."
        ),
    )
    aai_genus_boundary_high: float = Field(
        default=65.0,
        ge=0.0,
        le=100.0,
        description=(
            "Upper AAI boundary for genus classification. "
            "AAI above this value indicates same genus."
        ),
    )
    use_aai_for_genus: bool = Field(
        default=True,
        description=(
            "Use AAI instead of ANI for genus-level classification when available. "
            "AAI provides more reliable genus boundaries than ANI at high divergence."
        ),
    )

    # Family validation for broad-database alignments
    # When target_family is set, reads with better external hits are flagged Off-target
    target_family: str | None = Field(
        default=None,
        description=(
            "Target taxonomic family for family validation (e.g., 'f__Francisellaceae'). "
            "When set, reads with substantially better hits outside the ANI matrix "
            "are classified as Off-target. None disables family validation."
        ),
    )
    family_ratio_threshold: float = Field(
        default=0.8,
        ge=0.0,
        le=1.0,
        description=(
            "Family bitscore ratio threshold for off-target detection. "
            "Reads with best_in_family_bitscore / best_overall_bitscore below "
            "this value are classified as Off-target. Default 0.8."
        ),
    )

    # Bayesian classification configuration
    bayesian: BayesianConfig = Field(
        default_factory=BayesianConfig,
        description="Configuration for the Bayesian classification engine.",
    )

    # Output control
    include_legacy_scores: bool = Field(
        default=False,
        description=(
            "Include legacy sub-scores (alignment_quality, identity_confidence, "
            "placement_confidence, discovery_score) in output. These are replaced "
            "by posterior entropy and confidence_score in the Bayesian-primary design."
        ),
    )

    @model_validator(mode="after")
    def validate_threshold_boundaries(self) -> Self:
        """
        Validate threshold relationships after all fields are set.

        This ensures proper ordering of classification boundaries regardless
        of field definition order in the model.

        Boundary logic:
            - Known Species: N < novelty_known_max (exclusive upper bound)
            - Novel Species: novelty_novel_species_min <= N <= novelty_novel_species_max
            - Novel Genus: novelty_novel_genus_min <= N <= novelty_novel_genus_max

        For continuous boundaries, novelty_known_max == novelty_novel_species_min.
        """
        # Ensure continuous boundary between known and novel species
        # Known: N < known_max, Novel Species: N >= novel_species_min
        # They should be equal for gapless classification
        if self.novelty_novel_species_min < self.novelty_known_max:
            msg = (
                f"novelty_novel_species_min ({self.novelty_novel_species_min}) "
                f"must be >= novelty_known_max ({self.novelty_known_max})"
            )
            raise ValueError(msg)

        # Ensure novel genus min equals novel species max for continuous boundaries
        if self.novelty_novel_genus_min != self.novelty_novel_species_max:
            msg = (
                f"novelty_novel_genus_min ({self.novelty_novel_genus_min}) "
                f"must equal novelty_novel_species_max ({self.novelty_novel_species_max})"
            )
            raise ValueError(msg)

        # Ensure uncertainty thresholds form a valid hierarchy
        if self.uncertainty_novel_species_max > self.uncertainty_novel_genus_max:
            msg = (
                f"uncertainty_novel_species_max ({self.uncertainty_novel_species_max}) "
                f"must be <= uncertainty_novel_genus_max ({self.uncertainty_novel_genus_max})"
            )
            raise ValueError(msg)

        # Coverage weighting validation
        if self.coverage_weight_mode != "none" and self.coverage_weight_strength == 0.0:
            msg = (
                "coverage_weight_strength cannot be 0.0 when coverage_weight_mode is enabled. "
                "Set coverage_weight_mode='none' or increase coverage_weight_strength > 0.0."
            )
            raise ValueError(msg)

        return self

    def get_effective_thresholds(self) -> dict[str, float | tuple[float, ...]]:
        """
        Get effective novelty, uncertainty, and confidence score thresholds.

        For protein mode, returns wider thresholds that account for the greater
        conservation of amino acid sequences compared to nucleotide sequences.
        Also includes confidence score scaling parameters for mode-specific
        confidence calculation.

        Returns:
            Dictionary with effective threshold values for classification,
            including confidence score scaling parameters.
        """
        if self.alignment_mode == "protein":
            # Protein-level thresholds (wider ranges)
            from metadarkmatter.core.protein_constants import (
                PROTEIN_CONFIDENCE_IDENTITY_GAP_THRESHOLDS,
                PROTEIN_CONFIDENCE_IDENTITY_SCORE_BASE,
                PROTEIN_CONFIDENCE_IDENTITY_SCORE_RANGE,
                PROTEIN_CONFIDENCE_MARGIN_DIVISOR_KNOWN,
                PROTEIN_CONFIDENCE_MARGIN_DIVISOR_NOVEL_GENUS,
                PROTEIN_CONFIDENCE_MARGIN_DIVISOR_NOVEL_SPECIES,
                PROTEIN_NOVELTY_KNOWN_MAX,
                PROTEIN_NOVELTY_NOVEL_GENUS_MAX,
                PROTEIN_NOVELTY_NOVEL_GENUS_MIN,
                PROTEIN_NOVELTY_NOVEL_SPECIES_MAX,
                PROTEIN_NOVELTY_NOVEL_SPECIES_MIN,
                PROTEIN_UNCERTAINTY_CONFIDENT_MAX,
                PROTEIN_UNCERTAINTY_CONSERVED_MIN,
            )
            return {
                # Classification thresholds
                "novelty_known_max": PROTEIN_NOVELTY_KNOWN_MAX,
                "novelty_novel_species_min": PROTEIN_NOVELTY_NOVEL_SPECIES_MIN,
                "novelty_novel_species_max": PROTEIN_NOVELTY_NOVEL_SPECIES_MAX,
                "novelty_novel_genus_min": PROTEIN_NOVELTY_NOVEL_GENUS_MIN,
                "novelty_novel_genus_max": PROTEIN_NOVELTY_NOVEL_GENUS_MAX,
                "uncertainty_known_max": PROTEIN_UNCERTAINTY_CONFIDENT_MAX,
                "uncertainty_novel_species_max": PROTEIN_UNCERTAINTY_CONFIDENT_MAX,
                "uncertainty_novel_genus_max": PROTEIN_UNCERTAINTY_CONFIDENT_MAX,
                "uncertainty_conserved_min": PROTEIN_UNCERTAINTY_CONSERVED_MIN,
                # Confidence score scaling parameters
                "margin_divisor_known": PROTEIN_CONFIDENCE_MARGIN_DIVISOR_KNOWN,
                "margin_divisor_novel_species": PROTEIN_CONFIDENCE_MARGIN_DIVISOR_NOVEL_SPECIES,
                "margin_divisor_novel_genus": PROTEIN_CONFIDENCE_MARGIN_DIVISOR_NOVEL_GENUS,
                "identity_gap_thresholds": PROTEIN_CONFIDENCE_IDENTITY_GAP_THRESHOLDS,
                "identity_score_base": PROTEIN_CONFIDENCE_IDENTITY_SCORE_BASE,
                "identity_score_range": PROTEIN_CONFIDENCE_IDENTITY_SCORE_RANGE,
            }
        else:
            # Nucleotide-level thresholds (default)
            from metadarkmatter.core.constants import (
                CONFIDENCE_IDENTITY_GAP_THRESHOLDS,
                CONFIDENCE_IDENTITY_SCORE_BASE,
                CONFIDENCE_IDENTITY_SCORE_RANGE,
                CONFIDENCE_MARGIN_DIVISOR_KNOWN,
                CONFIDENCE_MARGIN_DIVISOR_NOVEL_GENUS,
                CONFIDENCE_MARGIN_DIVISOR_NOVEL_SPECIES,
            )
            return {
                # Classification thresholds
                "novelty_known_max": self.novelty_known_max,
                "novelty_novel_species_min": self.novelty_novel_species_min,
                "novelty_novel_species_max": self.novelty_novel_species_max,
                "novelty_novel_genus_min": self.novelty_novel_genus_min,
                "novelty_novel_genus_max": self.novelty_novel_genus_max,
                "uncertainty_known_max": self.uncertainty_known_max,
                "uncertainty_novel_species_max": self.uncertainty_novel_species_max,
                "uncertainty_novel_genus_max": self.uncertainty_novel_genus_max,
                "uncertainty_conserved_min": self.uncertainty_conserved_min,
                # Confidence score scaling parameters
                "margin_divisor_known": CONFIDENCE_MARGIN_DIVISOR_KNOWN,
                "margin_divisor_novel_species": CONFIDENCE_MARGIN_DIVISOR_NOVEL_SPECIES,
                "margin_divisor_novel_genus": CONFIDENCE_MARGIN_DIVISOR_NOVEL_GENUS,
                "identity_gap_thresholds": CONFIDENCE_IDENTITY_GAP_THRESHOLDS,
                "identity_score_base": CONFIDENCE_IDENTITY_SCORE_BASE,
                "identity_score_range": CONFIDENCE_IDENTITY_SCORE_RANGE,
            }

    @classmethod
    def for_protein_mode(cls) -> ScoringConfig:
        """
        Create a ScoringConfig instance configured for protein-level classification.

        This factory method sets up the config with alignment_mode="protein",
        which causes get_effective_thresholds() to return protein-specific
        thresholds during classification.

        Returns:
            ScoringConfig configured for protein-level classification.
        """
        return cls(alignment_mode="protein")

    @classmethod
    def from_yaml(cls, path: Path) -> ScoringConfig:
        """
        Load scoring configuration from a YAML file.

        The YAML file uses a nested structure that maps to ScoringConfig fields.
        Unknown keys are ignored (forward compatibility). Nested sections
        (thresholds, bayesian, filters, etc.) are flattened to match model fields.

        Args:
            path: Path to YAML configuration file.

        Returns:
            ScoringConfig populated from YAML values merged with defaults.

        Raises:
            FileNotFoundError: If YAML file does not exist.
            ValueError: If YAML contains invalid threshold values.
        """
        import yaml

        raw = yaml.safe_load(path.read_text())
        if not isinstance(raw, dict):
            msg = f"YAML config must be a mapping, got {type(raw).__name__}"
            raise ValueError(msg)

        flat = _flatten_yaml_config(raw)
        return cls(**flat)

    def to_yaml(self, path: Path) -> None:
        """
        Write scoring configuration to a YAML file.

        Produces a structured YAML file matching the documented config format.

        Args:
            path: Output file path.
        """
        path.write_text(self.to_yaml_str())

    def to_yaml_str(self) -> str:
        """
        Serialize scoring configuration to a YAML string.

        Returns:
            YAML-formatted string with nested structure.
        """
        import yaml

        data = _build_yaml_structure(self)
        return yaml.dump(data, default_flow_style=False, sort_keys=False)

    model_config = {"frozen": True}


class BlastConfig(BaseModel):
    """
    Configuration for BLAST execution optimized for novel diversity detection.

    Parameters are tuned for detecting divergent sequences (75-95% identity)
    while maintaining specificity. For placement uncertainty calculation,
    max_target_seqs should be high enough to capture all competing hits.

    Attributes:
        num_threads: Number of CPU threads for BLAST (default: 4)
        word_size: Word size for seed matches (default: 7 for divergent sequences)
        max_target_seqs: Maximum aligned sequences (default: 500 for
            comprehensive competitive alignment)
        evalue: E-value threshold (default: 1e-5 balances sensitivity/specificity)
        perc_identity: Minimum percent identity (default: 0.0, no filter)
        outfmt: BLAST output format string
    """

    num_threads: int = Field(default=4, ge=1, description="Number of CPU threads")
    word_size: int = Field(
        default=7,
        ge=4,
        le=28,
        description="Word size for seed matches (7 optimal for divergent sequences)",
    )
    max_target_seqs: int = Field(
        default=500,
        ge=1,
        description="Maximum aligned sequences (high value for complete competitive alignment)",
    )
    evalue: float = Field(
        default=1e-5,
        gt=0,
        description="E-value threshold (1e-5 balances sensitivity for 75%+ identity)",
    )
    perc_identity: float = Field(
        default=0.0,
        ge=0,
        le=100,
        description="Minimum percent identity filter",
    )
    outfmt: str = Field(
        default=(
            "6 qseqid sseqid pident length mismatch gapopen "
            "qstart qend sstart send evalue bitscore"
        ),
        description="BLAST output format specification",
    )

    model_config = {"frozen": True}


def _flatten_yaml_config(raw: dict[str, Any]) -> dict[str, Any]:
    """
    Flatten nested YAML config structure into ScoringConfig keyword arguments.

    Maps the documented nested YAML structure:
        thresholds.novelty.known_max -> novelty_known_max
        bayesian.priors -> bayesian.priors (nested BayesianConfig)
        filters.min_alignment_length -> min_alignment_length
    """
    flat: dict[str, Any] = {}

    # Top-level scalars
    if "alignment_mode" in raw:
        flat["alignment_mode"] = raw["alignment_mode"]

    # Threshold section
    thresholds = raw.get("thresholds", {})
    novelty = thresholds.get("novelty", {})
    uncertainty = thresholds.get("uncertainty", {})
    identity_gap_sec = thresholds.get("identity_gap", {})
    genus_unc = thresholds.get("genus_uncertainty", {})

    _map_if_present(novelty, "known_max", flat, "novelty_known_max")
    _map_if_present(novelty, "novel_species_min", flat, "novelty_novel_species_min")
    _map_if_present(novelty, "novel_species_max", flat, "novelty_novel_species_max")
    _map_if_present(novelty, "novel_genus_min", flat, "novelty_novel_genus_min")
    _map_if_present(novelty, "novel_genus_max", flat, "novelty_novel_genus_max")

    _map_if_present(uncertainty, "known_max", flat, "uncertainty_known_max")
    _map_if_present(uncertainty, "novel_species_max", flat, "uncertainty_novel_species_max")
    _map_if_present(uncertainty, "novel_genus_max", flat, "uncertainty_novel_genus_max")
    _map_if_present(uncertainty, "conserved_min", flat, "uncertainty_conserved_min")

    _map_if_present(identity_gap_sec, "ambiguous_max", flat, "identity_gap_ambiguous_max")
    _map_if_present(genus_unc, "ambiguous_min", flat, "genus_uncertainty_ambiguous_min")

    # Bayesian section -> BayesianConfig
    bayesian_raw = raw.get("bayesian", {})
    if bayesian_raw:
        bay_kwargs: dict[str, Any] = {}
        if "priors" in bayesian_raw:
            bay_kwargs["priors"] = bayesian_raw["priors"]
        modulation = bayesian_raw.get("prior_modulation", {})
        _map_if_present(modulation, "identity_gap_boost", bay_kwargs, "identity_gap_boost")
        _map_if_present(modulation, "single_hit_boost", bay_kwargs, "single_hit_boost")
        _map_if_present(bayesian_raw, "entropy_threshold", bay_kwargs, "entropy_threshold")
        flat["bayesian"] = BayesianConfig(**bay_kwargs)

    # Filters section
    filters = raw.get("filters", {})
    _map_if_present(filters, "min_alignment_length", flat, "min_alignment_length")
    _map_if_present(filters, "min_alignment_fraction", flat, "min_alignment_fraction")
    _map_if_present(filters, "bitscore_threshold_pct", flat, "bitscore_threshold_pct")

    # Coverage weighting section
    cw = raw.get("coverage_weighting", {})
    _map_if_present(cw, "mode", flat, "coverage_weight_mode")
    _map_if_present(cw, "strength", flat, "coverage_weight_strength")

    # Uncertainty section
    unc = raw.get("uncertainty", {})
    _map_if_present(unc, "mode", flat, "uncertainty_mode")
    _map_if_present(unc, "single_hit_threshold", flat, "single_hit_uncertainty_threshold")

    # Family validation section
    fv = raw.get("family_validation", {})
    _map_if_present(fv, "target_family", flat, "target_family")
    _map_if_present(fv, "ratio_threshold", flat, "family_ratio_threshold")

    # Output section
    output_sec = raw.get("output", {})
    _map_if_present(output_sec, "include_legacy_scores", flat, "include_legacy_scores")

    return flat


def _map_if_present(
    source: dict[str, Any],
    source_key: str,
    target: dict[str, Any],
    target_key: str,
) -> None:
    """Copy value from source dict to target dict if key exists."""
    if source_key in source and source[source_key] is not None:
        target[target_key] = source[source_key]


def _build_yaml_structure(config: "ScoringConfig") -> dict[str, Any]:
    """Build nested YAML dict from a ScoringConfig instance."""
    return {
        "alignment_mode": config.alignment_mode,
        "thresholds": {
            "novelty": {
                "known_max": config.novelty_known_max,
                "novel_species_min": config.novelty_novel_species_min,
                "novel_species_max": config.novelty_novel_species_max,
                "novel_genus_min": config.novelty_novel_genus_min,
                "novel_genus_max": config.novelty_novel_genus_max,
            },
            "uncertainty": {
                "known_max": config.uncertainty_known_max,
                "novel_species_max": config.uncertainty_novel_species_max,
                "novel_genus_max": config.uncertainty_novel_genus_max,
                "conserved_min": config.uncertainty_conserved_min,
            },
            "identity_gap": {
                "ambiguous_max": config.identity_gap_ambiguous_max,
            },
            "genus_uncertainty": {
                "ambiguous_min": config.genus_uncertainty_ambiguous_min,
            },
        },
        "bayesian": {
            "priors": dict(config.bayesian.priors),
            "prior_modulation": {
                "identity_gap_boost": config.bayesian.identity_gap_boost,
                "single_hit_boost": config.bayesian.single_hit_boost,
            },
            "entropy_threshold": config.bayesian.entropy_threshold,
        },
        "filters": {
            "min_alignment_length": config.min_alignment_length,
            "min_alignment_fraction": config.min_alignment_fraction,
            "bitscore_threshold_pct": config.bitscore_threshold_pct,
        },
        "coverage_weighting": {
            "mode": config.coverage_weight_mode,
            "strength": config.coverage_weight_strength,
        },
        "uncertainty": {
            "mode": config.uncertainty_mode,
            "single_hit_threshold": config.single_hit_uncertainty_threshold,
        },
        "family_validation": {
            "target_family": config.target_family,
            "ratio_threshold": config.family_ratio_threshold,
        },
        "output": {
            "include_legacy_scores": config.include_legacy_scores,
        },
    }


