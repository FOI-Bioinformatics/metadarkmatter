"""
Reusable threshold classification logic.

Extracts classification decision rules from the vectorized classifier
so that sensitivity analysis can re-classify without re-computing metrics.
"""

from __future__ import annotations

import polars as pl

from metadarkmatter.models.classification import TaxonomicCall, TAXONOMIC_TO_DIVERSITY
from metadarkmatter.models.config import ScoringConfig


def apply_classification_thresholds(
    df: pl.DataFrame,
    config: ScoringConfig,
) -> pl.DataFrame:
    """
    Apply threshold-based classification to a DataFrame with pre-computed metrics.

    This function applies the classification decision rules without recomputing
    the underlying metrics (novelty_index, placement_uncertainty, etc.).
    This allows efficient re-classification with different threshold values
    for sensitivity analysis.

    Required columns in df:
        - novelty_index: float (100 - top_hit_identity)
        - placement_uncertainty: float (100 - max_ANI between competing genomes)
        - num_ambiguous_hits: int
        - identity_gap: float or null

    Optional columns (used when present):
        - genus_uncertainty: float
        - num_secondary_genomes: int
        - ambiguity_scope: str

    Args:
        df: DataFrame with pre-computed classification metrics.
        config: ScoringConfig with threshold values to apply.

    Returns:
        DataFrame with updated taxonomic_call and diversity_status columns.
    """
    # Build effective thresholds (handles protein vs nucleotide mode)
    eff = _effective_thresholds(config)

    # Compute inferred uncertainty breakpoints (continuous piecewise-linear)
    inferred_known_break = 5.0 + eff["novelty_known_max"] * 0.5
    inferred_species_break = (
        inferred_known_break
        + (eff["novelty_novel_species_max"] - eff["novelty_known_max"]) * 1.0
    )

    # Preserve Off-target reads (already classified by family validation)
    has_off_target = "taxonomic_call" in df.columns
    off_target_rows = None
    if has_off_target:
        off_target_mask = pl.col("taxonomic_call") == "Off-target"
        off_target_rows = df.filter(off_target_mask)
        df = df.filter(~off_target_mask)

    # Check which optional columns exist
    has_genus = "genus_uncertainty" in df.columns
    has_scope = "ambiguity_scope" in df.columns

    # Build classification expression
    classification = (
        # Rule 0: Identity gap check
        pl.when(
            (pl.col("identity_gap").is_not_null())
            & (pl.col("num_ambiguous_hits") > 1)
            & (pl.col("identity_gap") < config.identity_gap_ambiguous_max)
        )
        .then(pl.lit("Ambiguous"))
    )

    # Conserved Region and high uncertainty rules
    if has_scope:
        classification = classification.when(
            (pl.col("placement_uncertainty") >= eff["uncertainty_conserved_min"])
            & (pl.col("ambiguity_scope") == "across_genera")
        ).then(pl.lit("Conserved Region"))

    classification = (
        classification
        .when(pl.col("placement_uncertainty") >= eff["uncertainty_conserved_min"])
        .then(pl.lit("Ambiguous"))
        # Species boundary zone
        .when(pl.col("placement_uncertainty") >= eff["uncertainty_novel_genus_max"])
        .then(pl.lit("Species Boundary"))
        # Known Species
        .when(
            (pl.col("novelty_index") < eff["novelty_known_max"])
            & (pl.col("placement_uncertainty") < eff["uncertainty_known_max"])
        )
        .then(pl.lit("Known Species"))
        # Rule 0b: Single-hit inferred uncertainty check
        .when(
            (pl.col("num_ambiguous_hits") <= 1)
            & (pl.col("novelty_index") >= eff["novelty_novel_species_min"])
            & (
                pl.when(pl.col("novelty_index") < eff["novelty_known_max"])
                .then(5.0 + pl.col("novelty_index") * 0.5)
                .when(pl.col("novelty_index") < eff["novelty_novel_species_max"])
                .then(inferred_known_break + (pl.col("novelty_index") - eff["novelty_known_max"]) * 1.0)
                .when(pl.col("novelty_index") < eff["novelty_novel_genus_max"])
                .then(inferred_species_break + (pl.col("novelty_index") - eff["novelty_novel_species_max"]) * 1.5)
                .otherwise(35.0)
                >= pl.lit(config.single_hit_uncertainty_threshold)
            )
        )
        .then(pl.lit("Ambiguous"))
        # Novel Species
        .when(
            (pl.col("novelty_index") >= eff["novelty_novel_species_min"])
            & (pl.col("novelty_index") < eff["novelty_novel_species_max"])
            & (pl.col("placement_uncertainty") < eff["uncertainty_novel_species_max"])
        )
        .then(pl.lit("Novel Species"))
    )

    # Novel Genus with optional genus_uncertainty check
    if has_genus:
        classification = classification.when(
            (pl.col("novelty_index") >= eff["novelty_novel_genus_min"])
            & (pl.col("novelty_index") <= eff["novelty_novel_genus_max"])
            & (pl.col("placement_uncertainty") < eff["uncertainty_novel_genus_max"])
            & (
                (pl.col("genus_uncertainty") < config.genus_uncertainty_ambiguous_min)
                | (pl.col("num_secondary_genomes") <= 0)
            )
        ).then(pl.lit("Novel Genus"))
    else:
        classification = classification.when(
            (pl.col("novelty_index") >= eff["novelty_novel_genus_min"])
            & (pl.col("novelty_index") <= eff["novelty_novel_genus_max"])
            & (pl.col("placement_uncertainty") < eff["uncertainty_novel_genus_max"])
        ).then(pl.lit("Novel Genus"))

    classification = classification.otherwise(pl.lit("Unclassified")).alias("taxonomic_call")

    # Apply classification and derived columns
    result = df.with_columns([classification])
    result = result.with_columns([
        pl.col("taxonomic_call").replace(TAXONOMIC_TO_DIVERSITY).alias("diversity_status"),
        pl.col("taxonomic_call").is_in(["Novel Species", "Novel Genus"]).alias("is_novel"),
    ])

    # Re-add preserved Off-target rows
    if off_target_rows is not None and not off_target_rows.is_empty():
        off_target_rows = off_target_rows.with_columns([
            pl.lit("Off-target").alias("taxonomic_call"),
            pl.lit("Uncertain").alias("diversity_status"),
            pl.lit(False).alias("is_novel"),
        ])
        result = pl.concat([result, off_target_rows], how="diagonal_relaxed")

    return result


def _effective_thresholds(config: ScoringConfig) -> dict[str, float]:
    """
    Build effective threshold dictionary from config.

    Handles protein vs nucleotide mode by checking alignment_mode.
    """
    if config.alignment_mode == "protein":
        from metadarkmatter.core.protein_constants import (
            PROTEIN_NOVELTY_KNOWN_MAX,
            PROTEIN_NOVELTY_NOVEL_SPECIES_MIN,
            PROTEIN_NOVELTY_NOVEL_SPECIES_MAX,
            PROTEIN_NOVELTY_NOVEL_GENUS_MIN,
            PROTEIN_NOVELTY_NOVEL_GENUS_MAX,
            PROTEIN_UNCERTAINTY_KNOWN_MAX,
            PROTEIN_UNCERTAINTY_NOVEL_SPECIES_MAX,
            PROTEIN_UNCERTAINTY_NOVEL_GENUS_MAX,
            PROTEIN_UNCERTAINTY_CONSERVED_MIN,
        )
        return {
            "novelty_known_max": config.novelty_known_max if config.novelty_known_max != 4.0 else PROTEIN_NOVELTY_KNOWN_MAX,
            "novelty_novel_species_min": config.novelty_novel_species_min if config.novelty_novel_species_min != 4.0 else PROTEIN_NOVELTY_NOVEL_SPECIES_MIN,
            "novelty_novel_species_max": config.novelty_novel_species_max if config.novelty_novel_species_max != 20.0 else PROTEIN_NOVELTY_NOVEL_SPECIES_MAX,
            "novelty_novel_genus_min": config.novelty_novel_genus_min if config.novelty_novel_genus_min != 20.0 else PROTEIN_NOVELTY_NOVEL_GENUS_MIN,
            "novelty_novel_genus_max": config.novelty_novel_genus_max if config.novelty_novel_genus_max != 25.0 else PROTEIN_NOVELTY_NOVEL_GENUS_MAX,
            "uncertainty_known_max": PROTEIN_UNCERTAINTY_KNOWN_MAX,
            "uncertainty_novel_species_max": PROTEIN_UNCERTAINTY_NOVEL_SPECIES_MAX,
            "uncertainty_novel_genus_max": PROTEIN_UNCERTAINTY_NOVEL_GENUS_MAX,
            "uncertainty_conserved_min": PROTEIN_UNCERTAINTY_CONSERVED_MIN,
        }
    else:
        return {
            "novelty_known_max": config.novelty_known_max,
            "novelty_novel_species_min": config.novelty_novel_species_min,
            "novelty_novel_species_max": config.novelty_novel_species_max,
            "novelty_novel_genus_min": config.novelty_novel_genus_min,
            "novelty_novel_genus_max": config.novelty_novel_genus_max,
            "uncertainty_known_max": config.uncertainty_known_max,
            "uncertainty_novel_species_max": config.uncertainty_novel_species_max,
            "uncertainty_novel_genus_max": config.uncertainty_novel_genus_max,
            "uncertainty_conserved_min": config.uncertainty_conserved_min,
        }
