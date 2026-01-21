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


# =============================================================================
# Enhanced Scoring Functions for Novel Diversity Detection
#
# These functions address limitations in the base scoring system:
# 1. Single-hit reads now get inferred uncertainty instead of 0%
# 2. Alignment quality uses mismatch/gap/coverage/evalue
# 3. Orthogonal confidence dimensions (identity vs placement)
# 4. Discovery score prioritizes novel findings for validation
# =============================================================================

# Inferred uncertainty thresholds
INFERRED_UNCERTAINTY_BASE = 5.0
INFERRED_UNCERTAINTY_NOVEL_SPECIES_SLOPE = 1.0
INFERRED_UNCERTAINTY_NOVEL_GENUS_SLOPE = 1.5
INFERRED_UNCERTAINTY_MAX = 35.0


def calculate_inferred_uncertainty(novelty_index: float) -> float:
    """
    Calculate inferred uncertainty for single-hit reads.

    For reads with only one hit (num_ambiguous_hits <= 1), we cannot measure
    placement uncertainty from ANI between competing genomes. Instead, we
    infer uncertainty from the novelty level: higher novelty suggests either
    novel diversity OR database incompleteness, both of which increase
    placement uncertainty.

    Args:
        novelty_index: Novelty index value (100 - top_hit_identity)

    Returns:
        Inferred uncertainty value (5-35%)

    Interpretation:
        - Low novelty (< 5%): Database likely complete for this species
        - Novel species range (5-15%): Uncertain if truly novel or database gap
        - Novel genus range (15-25%): High uncertainty about placement
        - Very high divergence (> 25%): Maximum uncertainty
    """
    if novelty_index < NOVELTY_KNOWN_MAX:
        # High identity: database likely complete for this species
        return INFERRED_UNCERTAINTY_BASE + novelty_index * 0.5  # 5-7.5%
    elif novelty_index < NOVELTY_NOVEL_SPECIES_MAX:
        # Novel species range: uncertain if truly novel or database gap
        return 7.5 + (novelty_index - NOVELTY_KNOWN_MAX) * INFERRED_UNCERTAINTY_NOVEL_SPECIES_SLOPE  # 7.5-17.5%
    elif novelty_index < NOVELTY_NOVEL_GENUS_MAX:
        # Novel genus range: high uncertainty
        return 17.5 + (novelty_index - NOVELTY_NOVEL_SPECIES_MAX) * INFERRED_UNCERTAINTY_NOVEL_GENUS_SLOPE  # 17.5-25%
    else:
        # Very high divergence: maximum uncertainty
        return INFERRED_UNCERTAINTY_MAX


def calculate_alignment_quality(
    mismatch: int,
    gapopen: int,
    length: int,
    coverage: float,
    evalue: float,
) -> float:
    """
    Calculate alignment quality score from BLAST statistics.

    Incorporates alignment statistics that are captured but not currently
    used in the base scoring system: mismatch count, gap opens, alignment
    coverage, and e-value.

    Args:
        mismatch: Number of mismatches in alignment
        gapopen: Number of gap openings in alignment
        length: Alignment length in bp
        coverage: Fraction of read aligned (0.0-1.0)
        evalue: BLAST e-value

    Returns:
        Alignment quality score (0-100)

    Components (each 0-25 points):
        - Mismatch density: Fewer mismatches = higher score
        - Gap complexity: Fewer gaps = higher score
        - Coverage: Higher coverage = higher score
        - E-value significance: Lower e-value = higher score
    """
    # Mismatch density penalty (0-25 pts)
    # Higher mismatch rate = lower score
    mismatch_rate = mismatch / max(1, length)
    mismatch_score = max(0.0, 25.0 - mismatch_rate * 50.0)

    # Gap complexity penalty (0-25 pts)
    # More gap openings per 100bp = lower score
    gap_rate = gapopen / max(1, length) * 100.0
    gap_score = max(0.0, 25.0 - gap_rate * 25.0)

    # Coverage bonus (0-25 pts)
    # Higher coverage = higher score
    coverage_score = min(25.0, coverage * 25.0)

    # E-value significance (0-25 pts)
    # Lower e-value = higher significance = higher score
    if evalue <= 1e-50:
        evalue_score = 25.0
    elif evalue <= 1e-20:
        evalue_score = 20.0
    elif evalue <= 1e-10:
        evalue_score = 15.0
    elif evalue <= 1e-5:
        evalue_score = 10.0
    else:
        evalue_score = 5.0

    return round(mismatch_score + gap_score + coverage_score + evalue_score, 1)


def calculate_identity_confidence(
    novelty_index: float,
    alignment_quality: float,
) -> float:
    """
    Calculate confidence in the identity measurement.

    This measures how reliable the percent identity value is, independent
    of placement confidence. High identity with good alignment quality
    means we can trust what this sequence is.

    Args:
        novelty_index: Novelty index value (100 - top_hit_identity)
        alignment_quality: Alignment quality score (0-100)

    Returns:
        Identity confidence score (0-100)

    Interpretation:
        - High score: Confident about the identity measurement
        - Low score: Identity measurement may be unreliable
    """
    # Base score from novelty (inverted - higher identity = more confident)
    if novelty_index < 5:
        base = 80.0  # Very high identity
    elif novelty_index < 15:
        base = 60.0 - (novelty_index - 5)  # 60-50
    elif novelty_index < 25:
        base = 50.0 - (novelty_index - 15) * 1.5  # 50-35
    else:
        base = 30.0

    # Alignment quality contribution (up to +20)
    quality_bonus = alignment_quality * 0.2

    return round(min(100.0, base + quality_bonus), 1)


def calculate_placement_confidence(
    uncertainty: float,
    uncertainty_type: str,
    identity_gap: float | None,
    num_ambiguous_hits: int,
) -> float:
    """
    Calculate confidence in genome assignment.

    This measures how confident we are about which genome this read
    belongs to, accounting for whether uncertainty was measured or
    inferred.

    Args:
        uncertainty: Placement uncertainty (measured OR inferred)
        uncertainty_type: "measured" (from ANI) or "inferred" (from novelty)
        identity_gap: Gap between best and second-best hit identity
        num_ambiguous_hits: Number of competing hits

    Returns:
        Placement confidence score (0-100)

    Key insight:
        Inferred uncertainty is penalized because we cannot directly
        measure competing placements - absence of data is not evidence
        of confident placement.
    """
    # Base score from uncertainty level
    if uncertainty < 2:
        base = 80.0
    elif uncertainty < 5:
        base = 60.0
    elif uncertainty < 10:
        base = 40.0
    elif uncertainty < 20:
        base = 25.0
    else:
        base = 10.0

    # Penalty for inferred (unmeasured) uncertainty
    # We can't be as confident without actual competing hit data
    if uncertainty_type == "inferred":
        base -= 15.0

    # Identity gap bonus (only meaningful for multi-hit reads)
    gap_bonus = 0.0
    if identity_gap is not None and num_ambiguous_hits > 1:
        if identity_gap >= 5:
            gap_bonus = 20.0
        elif identity_gap >= 2:
            gap_bonus = 10.0
        elif identity_gap >= 1:
            gap_bonus = 5.0

    return round(max(0.0, min(100.0, base + gap_bonus)), 1)


def calculate_discovery_score(
    taxonomic_call: str,
    novelty_index: float,
    identity_confidence: float,
    placement_confidence: float,
    alignment_quality: float,
) -> float | None:
    """
    Calculate discovery priority score for novel reads.

    This score helps prioritize which novel findings to investigate first.
    Only calculated for Novel Species and Novel Genus classifications.

    Args:
        taxonomic_call: Classification category string
        novelty_index: Novelty index value (100 - top_hit_identity)
        identity_confidence: Identity confidence score (0-100)
        placement_confidence: Placement confidence score (0-100)
        alignment_quality: Alignment quality score (0-100)

    Returns:
        Discovery score (0-100) for novel reads, None for non-novel

    Score interpretation:
        75-100: High-confidence discovery - prioritize for validation
        50-74: Moderate discovery signal - include in candidate list
        25-49: Low-confidence discovery - needs more evidence
        < 25: Unreliable signal - likely artifact
    """
    if taxonomic_call not in (CATEGORY_NOVEL_SPECIES, CATEGORY_NOVEL_GENUS):
        return None

    # Novelty component (0-40): Higher divergence = more interesting
    if taxonomic_call == CATEGORY_NOVEL_SPECIES:
        # Novel species: 5-20% novelty -> 15-30 points
        novelty_pts = 15.0 + (novelty_index - NOVELTY_NOVEL_SPECIES_MIN) * 1.5
    else:  # Novel Genus
        # Novel genus: 20-25% novelty -> 30-40 points
        novelty_pts = 30.0 + (novelty_index - NOVELTY_NOVEL_GENUS_MIN) * 2.0
    novelty_pts = min(40.0, novelty_pts)

    # Quality component (0-30): Better alignment = more trustworthy
    quality_pts = alignment_quality * 0.3

    # Confidence component (0-30): Average of identity and placement
    conf_pts = (identity_confidence + placement_confidence) / 2.0 * 0.3

    return round(novelty_pts + quality_pts + conf_pts, 1)
