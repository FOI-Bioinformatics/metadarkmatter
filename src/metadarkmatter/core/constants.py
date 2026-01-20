"""
Constants used throughout the metadarkmatter package.

Centralizes magic strings, default values, and threshold constants
to improve maintainability and consistency.
"""

from __future__ import annotations

# =============================================================================
# Genome and Classification Constants
# =============================================================================

# Default placeholder for unknown or unresolved genome names
UNKNOWN_GENOME = "unknown"

# BLAST output format string for standard 12-column tabular output
BLAST_OUTFMT_12COL = (
    "6 qseqid sseqid pident length mismatch gapopen "
    "qstart qend sstart send evalue bitscore"
)

# =============================================================================
# ANI Thresholds - Literature-backed Values
#
# References:
# - Jain et al. 2018, Nature Communications: 95-96% ANI species boundary
# - Parks et al. 2020, Nature Biotechnology: GTDB taxonomy standards
# =============================================================================

# ANI values below which genomes are considered different genera
ANI_GENUS_BOUNDARY = 75.0

# ANI values at which genomes are considered different species
ANI_SPECIES_BOUNDARY_LOW = 95.0
ANI_SPECIES_BOUNDARY_HIGH = 96.0

# Default ANI for unrelated genomes (fallback when no ANI available)
ANI_DEFAULT_UNRELATED = 70.0

# Minimum ANI to consider genomes related enough for comparison
ANI_MIN_RELATED = 75.0

# =============================================================================
# AAI Thresholds - Literature-backed Values for Genus-Level Classification
#
# References:
# - Riesco & Trujillo 2024: 58-65% AAI genus boundary
# - Rodriguez-R & Konstantinidis 2014: AAI for taxonomic classification
# - Luo et al. 2014: AAI-based genus delineation
#
# AAI is more reliable than ANI for genus-level classification because
# ANI becomes unreliable below approximately 80% identity.
# =============================================================================

# AAI boundary range for same genus (values above 65% = same genus)
AAI_GENUS_BOUNDARY_HIGH = 65.0

# AAI boundary range for different genus (values below 58% = different genus)
AAI_GENUS_BOUNDARY_LOW = 58.0

# Default AAI for unrelated genomes (fallback when no AAI available)
AAI_DEFAULT_UNRELATED = 50.0

# Minimum AAI to consider genomes related enough for comparison
AAI_MIN_RELATED = 45.0

# =============================================================================
# Classification Categories
#
# These match the TaxonomicCall enum values in models/classification.py
# Use the enum directly when possible; these are for display/string matching.
# =============================================================================

CATEGORY_KNOWN_SPECIES = "Known Species"
CATEGORY_NOVEL_SPECIES = "Novel Species"
CATEGORY_NOVEL_GENUS = "Novel Genus"
CATEGORY_CONSERVED_REGION = "Conserved Region"
CATEGORY_AMBIGUOUS = "Ambiguous"
CATEGORY_UNCLASSIFIED = "Unclassified"

# =============================================================================
# Novelty Index Thresholds (percentage divergence from reference)
#
# Novelty Index = 100 - TopHitIdentity
# =============================================================================

# Maximum novelty for "Known Species" classification
NOVELTY_KNOWN_MAX = 5.0

# Range for "Novel Species" classification
# Note: Read-level BLAST identity can be 10-20% lower than genome-level ANI.
# The 20% threshold accounts for this gap (80-95% identity = novel species).
NOVELTY_NOVEL_SPECIES_MIN = 5.0
NOVELTY_NOVEL_SPECIES_MAX = 20.0

# Range for "Novel Genus" classification
# 75-80% identity corresponds to genus-level divergence.
NOVELTY_NOVEL_GENUS_MIN = 20.0
NOVELTY_NOVEL_GENUS_MAX = 25.0

# =============================================================================
# Placement Uncertainty Thresholds
#
# Placement Uncertainty = 100 - max(ANI between competing genomes)
# =============================================================================

# Maximum uncertainty for confident placement (Known/Novel Species/Genus)
UNCERTAINTY_CONFIDENT_MAX = 2.0

# Threshold above which hits are considered "Conserved Region"
UNCERTAINTY_CONSERVED_MIN = 5.0

# Range for "Ambiguous" classification
UNCERTAINTY_AMBIGUOUS_MIN = 2.0
UNCERTAINTY_AMBIGUOUS_MAX = 5.0


# =============================================================================
# Confidence Score Scaling Parameters (Nucleotide Mode)
#
# These parameters control how margin-from-boundary is converted to points.
# Protein mode uses wider scaling factors due to greater sequence conservation.
# =============================================================================

# Margin divisors for scaling novelty margin to points (0-40)
# Lower divisor = more sensitive to small margins
CONFIDENCE_MARGIN_DIVISOR_KNOWN = 5.0      # Known species: 5% margin -> 40 pts
CONFIDENCE_MARGIN_DIVISOR_NOVEL_SPECIES = 7.5  # Novel species: 7.5% margin -> 40 pts
CONFIDENCE_MARGIN_DIVISOR_NOVEL_GENUS = 5.0    # Novel genus: 5% margin -> 40 pts

# Identity gap thresholds for placement certainty scoring
# These are the minimum gaps (%) needed for each point tier
CONFIDENCE_IDENTITY_GAP_THRESHOLDS = (5.0, 2.0, 1.0, 0.5)  # -> 20, 15, 10, 5 pts

# Identity score scaling: (identity - base) / range * 20
# Nucleotide: 70% = 0 pts, 95% = 20 pts
CONFIDENCE_IDENTITY_SCORE_BASE = 70.0
CONFIDENCE_IDENTITY_SCORE_RANGE = 25.0


# =============================================================================
# Confidence Score Calculation
#
# Multi-factor confidence score (0-100) that quantifies classification
# reliability without changing the classification category.
# =============================================================================

def calculate_confidence_score(
    novelty_index: float,
    placement_uncertainty: float,
    num_ambiguous_hits: int,
    identity_gap: float | None,
    top_hit_identity: float,
    taxonomic_call: str,
    novelty_known_max: float = NOVELTY_KNOWN_MAX,
    novelty_novel_species_max: float = NOVELTY_NOVEL_SPECIES_MAX,
    novelty_novel_genus_max: float = NOVELTY_NOVEL_GENUS_MAX,
    uncertainty_confident_max: float = UNCERTAINTY_CONFIDENT_MAX,
    # Scaling parameters for mode-specific confidence calculation
    margin_divisor_known: float = CONFIDENCE_MARGIN_DIVISOR_KNOWN,
    margin_divisor_novel_species: float = CONFIDENCE_MARGIN_DIVISOR_NOVEL_SPECIES,
    margin_divisor_novel_genus: float = CONFIDENCE_MARGIN_DIVISOR_NOVEL_GENUS,
    identity_gap_thresholds: tuple[float, ...] = CONFIDENCE_IDENTITY_GAP_THRESHOLDS,
    identity_score_base: float = CONFIDENCE_IDENTITY_SCORE_BASE,
    identity_score_range: float = CONFIDENCE_IDENTITY_SCORE_RANGE,
) -> float:
    """
    Calculate a multi-factor confidence score for a read classification.

    The confidence score integrates:
    1. Margin from threshold boundaries (40 points max)
    2. Placement certainty from hit count and identity gap (40 points max)
    3. Alignment quality proxy from identity (20 points max)

    This provides uncertainty quantification without changing the
    classification category itself.

    Args:
        novelty_index: Read novelty (100 - top_hit_identity)
        placement_uncertainty: 100 - max(ANI between competing genomes)
        num_ambiguous_hits: Number of hits within bitscore threshold
        identity_gap: Gap between top hit and second-best different-genome hit
        top_hit_identity: Percent identity of best BLAST hit
        taxonomic_call: The classification category string
        novelty_known_max: Threshold for known species (default 5.0)
        novelty_novel_species_max: Upper bound for novel species (default 20.0)
        novelty_novel_genus_max: Upper bound for novel genus (default 25.0)
        uncertainty_confident_max: Max uncertainty for confident placement (default 2.0)
        margin_divisor_known: Divisor for known species margin scaling
        margin_divisor_novel_species: Divisor for novel species margin scaling
        margin_divisor_novel_genus: Divisor for novel genus margin scaling
        identity_gap_thresholds: Tuple of (high, mid-high, mid, low) gap thresholds
        identity_score_base: Base identity for alignment quality scoring
        identity_score_range: Range for identity score calculation

    Returns:
        Confidence score from 0-100 (higher = more confident)
    """
    score = 0.0

    # Component 1: Margin from threshold boundaries (0-40 points)
    # The further from classification boundaries, the more confident
    margin_score = 0.0

    if taxonomic_call == CATEGORY_KNOWN_SPECIES:
        # Known species: margin from upper novelty threshold
        margin_from_boundary = novelty_known_max - novelty_index
        margin_score = min(40.0, (margin_from_boundary / margin_divisor_known) * 40.0)

    elif taxonomic_call == CATEGORY_NOVEL_SPECIES:
        # Novel species: margin from both boundaries
        margin_lower = novelty_index - novelty_known_max
        margin_upper = novelty_novel_species_max - novelty_index
        min_margin = min(margin_lower, margin_upper)
        margin_score = min(40.0, (min_margin / margin_divisor_novel_species) * 40.0)

    elif taxonomic_call == CATEGORY_NOVEL_GENUS:
        # Novel genus: margin from both boundaries
        margin_lower = novelty_index - novelty_novel_species_max
        margin_upper = novelty_novel_genus_max - novelty_index
        min_margin = min(margin_lower, margin_upper)
        margin_score = min(40.0, (min_margin / margin_divisor_novel_genus) * 40.0)

    else:
        # Ambiguous/Unclassified: low base confidence
        margin_score = 10.0

    # Add uncertainty margin bonus (confident placement adds up to 10 pts)
    if placement_uncertainty < uncertainty_confident_max:
        uncertainty_margin = uncertainty_confident_max - placement_uncertainty
        margin_score = min(
            40.0, margin_score + (uncertainty_margin / uncertainty_confident_max) * 10.0
        )

    score += max(0.0, margin_score)

    # Component 2: Placement certainty (0-40 points)
    # Based on number of competing hits and identity gap
    placement_score = 0.0

    # Hit count factor: fewer hits = more confident (0-20 pts)
    if num_ambiguous_hits <= 1:
        placement_score += 20.0
    elif num_ambiguous_hits <= 3:
        placement_score += 15.0
    elif num_ambiguous_hits <= 5:
        placement_score += 10.0
    elif num_ambiguous_hits <= 10:
        placement_score += 5.0
    # >10 hits = 0 pts

    # Identity gap factor: larger gap = more confident (0-20 pts)
    # Use parameterized thresholds for mode-specific scoring
    if identity_gap is not None:
        gap_high, gap_mid_high, gap_mid, gap_low = identity_gap_thresholds
        if identity_gap >= gap_high:
            placement_score += 20.0  # Clear winner
        elif identity_gap >= gap_mid_high:
            placement_score += 15.0  # Good separation
        elif identity_gap >= gap_mid:
            placement_score += 10.0  # Moderate separation
        elif identity_gap >= gap_low:
            placement_score += 5.0   # Marginal separation
        # Below gap_low = 0 pts (very close competition)
    else:
        # No secondary hit = likely confident (single genome hit)
        if num_ambiguous_hits <= 1:
            placement_score += 20.0

    score += placement_score

    # Component 3: Alignment quality proxy (0-20 points)
    # Higher identity generally indicates better alignment quality
    identity_score = min(
        20.0,
        max(0.0, (top_hit_identity - identity_score_base) / identity_score_range * 20.0)
    )
    score += identity_score

    return round(min(100.0, max(0.0, score)), 1)
