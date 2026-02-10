"""
Protein-level constants for BLASTX-based read classification.

Protein identity operates in a wider range than nucleotide identity because
amino acid sequences are more conserved than nucleotide sequences. This module
defines thresholds for protein-level read classification that account for
this difference.

Threshold Rationale:
    Nucleotide-to-protein identity mapping:
    - Nucleotide 95% identity ~ Protein 90% identity (species boundary)
    - Nucleotide 85% identity ~ Protein 75% identity (novel species)
    - Nucleotide 75% identity ~ Protein 60% identity (novel genus)

    The wider ranges reflect the degeneracy of the genetic code - multiple
    codons can encode the same amino acid, so protein sequences change more
    slowly than nucleotide sequences during evolution.

References:
    - Konstantinidis & Tiedje 2005: ANI/AAI relationships
    - Qin et al. 2014: AAI thresholds for prokaryotic taxonomy
    - Rodriguez-R & Konstantinidis 2014: AAI for taxonomic classification
"""

from __future__ import annotations


# =============================================================================
# Protein Novelty Index Thresholds
#
# Novelty Index = 100 - TopHitIdentity
# These thresholds are wider than nucleotide thresholds because protein
# sequences are more conserved.
# =============================================================================

# Maximum novelty for "Known Species" classification at protein level
# Protein identity > 90% corresponds to nucleotide identity > 95%
PROTEIN_NOVELTY_KNOWN_MAX = 10.0

# Range for "Novel Species" classification at protein level
# Protein identity 75-90% corresponds to nucleotide identity ~85-95%
PROTEIN_NOVELTY_NOVEL_SPECIES_MIN = 10.0
PROTEIN_NOVELTY_NOVEL_SPECIES_MAX = 25.0

# Range for "Novel Genus" classification at protein level
# Protein identity 60-75% corresponds to nucleotide identity ~75-85%
PROTEIN_NOVELTY_NOVEL_GENUS_MIN = 25.0
PROTEIN_NOVELTY_NOVEL_GENUS_MAX = 40.0


# =============================================================================
# Protein Placement Uncertainty Thresholds
#
# These use AAI-based thresholds since BLASTX output reflects protein identity.
# Placement Uncertainty = 100 - max(AAI between competing genomes)
# =============================================================================

# Maximum uncertainty for confident placement at protein level
# Using 5% to account for wider protein identity ranges
PROTEIN_UNCERTAINTY_CONFIDENT_MAX = 5.0
PROTEIN_UNCERTAINTY_KNOWN_MAX = PROTEIN_UNCERTAINTY_CONFIDENT_MAX
PROTEIN_UNCERTAINTY_NOVEL_SPECIES_MAX = PROTEIN_UNCERTAINTY_CONFIDENT_MAX
PROTEIN_UNCERTAINTY_NOVEL_GENUS_MAX = PROTEIN_UNCERTAINTY_CONFIDENT_MAX

# Threshold above which hits are considered "Conserved Region"
PROTEIN_UNCERTAINTY_CONSERVED_MIN = 10.0


# =============================================================================
# Protein Confidence Score Scaling Parameters
#
# These parameters control how margin-from-boundary is converted to points.
# Protein mode uses wider scaling factors due to greater sequence conservation
# at the amino acid level compared to nucleotide sequences.
# =============================================================================

# Margin divisors for scaling novelty margin to points (0-40)
# Wider divisors reflect the larger threshold ranges in protein mode
PROTEIN_CONFIDENCE_MARGIN_DIVISOR_KNOWN = 10.0         # 10% margin -> 40 pts
PROTEIN_CONFIDENCE_MARGIN_DIVISOR_NOVEL_SPECIES = 15.0 # 15% margin -> 40 pts
PROTEIN_CONFIDENCE_MARGIN_DIVISOR_NOVEL_GENUS = 15.0   # 15% margin -> 40 pts

# Identity gap thresholds for placement certainty scoring
# Wider gaps expected at protein level due to greater conservation
PROTEIN_CONFIDENCE_IDENTITY_GAP_THRESHOLDS = (10.0, 5.0, 2.0, 1.0)  # -> 20, 15, 10, 5 pts

# Identity score scaling: (identity - base) / range * 20
# Protein: 50% = 0 pts, 90% = 20 pts (wider range due to divergent homology)
PROTEIN_CONFIDENCE_IDENTITY_SCORE_BASE = 50.0
PROTEIN_CONFIDENCE_IDENTITY_SCORE_RANGE = 40.0


# =============================================================================
# Confidence Score Calculation for Protein Mode
# =============================================================================

def calculate_protein_confidence_score(
    novelty_index: float,
    placement_uncertainty: float,
    num_ambiguous_hits: int,
    identity_gap: float | None,
    top_hit_identity: float,
    taxonomic_call: str,
    novelty_known_max: float = PROTEIN_NOVELTY_KNOWN_MAX,
    novelty_novel_species_max: float = PROTEIN_NOVELTY_NOVEL_SPECIES_MAX,
    novelty_novel_genus_max: float = PROTEIN_NOVELTY_NOVEL_GENUS_MAX,
    uncertainty_confident_max: float = PROTEIN_UNCERTAINTY_CONFIDENT_MAX,
) -> float:
    """
    Calculate confidence score for protein-level classification.

    Uses protein-specific thresholds that account for the wider identity
    ranges observed in amino acid sequences compared to nucleotides.

    This function delegates to the parameterized calculate_confidence_score
    with protein-specific scaling factors.

    Args:
        novelty_index: Read novelty (100 - top_hit_identity)
        placement_uncertainty: 100 - max(AAI between competing genomes)
        num_ambiguous_hits: Number of hits within bitscore threshold
        identity_gap: Gap between top hit and second-best different-genome hit
        top_hit_identity: Percent identity of best BLASTX hit
        taxonomic_call: The classification category string
        novelty_known_max: Threshold for known species (default 10.0)
        novelty_novel_species_max: Upper bound for novel species (default 25.0)
        novelty_novel_genus_max: Upper bound for novel genus (default 40.0)
        uncertainty_confident_max: Max uncertainty for confident placement (default 5.0)

    Returns:
        Confidence score from 0-100 (higher = more confident)
    """
    from metadarkmatter.core.constants import calculate_confidence_score

    return calculate_confidence_score(
        novelty_index=novelty_index,
        placement_uncertainty=placement_uncertainty,
        num_ambiguous_hits=num_ambiguous_hits,
        identity_gap=identity_gap,
        top_hit_identity=top_hit_identity,
        taxonomic_call=taxonomic_call,
        novelty_known_max=novelty_known_max,
        novelty_novel_species_max=novelty_novel_species_max,
        novelty_novel_genus_max=novelty_novel_genus_max,
        uncertainty_confident_max=uncertainty_confident_max,
        # Protein-specific scaling parameters
        margin_divisor_known=PROTEIN_CONFIDENCE_MARGIN_DIVISOR_KNOWN,
        margin_divisor_novel_species=PROTEIN_CONFIDENCE_MARGIN_DIVISOR_NOVEL_SPECIES,
        margin_divisor_novel_genus=PROTEIN_CONFIDENCE_MARGIN_DIVISOR_NOVEL_GENUS,
        identity_gap_thresholds=PROTEIN_CONFIDENCE_IDENTITY_GAP_THRESHOLDS,
        identity_score_base=PROTEIN_CONFIDENCE_IDENTITY_SCORE_BASE,
        identity_score_range=PROTEIN_CONFIDENCE_IDENTITY_SCORE_RANGE,
    )
