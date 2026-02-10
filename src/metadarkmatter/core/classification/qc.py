"""
Quality control metrics for classification pipeline.

Computes pre-classification and post-classification QC metrics
to help users identify problematic inputs and interpret results.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import polars as pl

if TYPE_CHECKING:
    from metadarkmatter.core.classification.ani_matrix import ANIMatrix
    from metadarkmatter.models.config import ScoringConfig


@dataclass
class QCMetrics:
    """Quality control metrics for a classification run."""

    # Pre-classification metrics
    total_alignments: int = 0
    filtered_alignments: int = 0
    filter_rate: float = 0.0
    genomes_in_alignment: int = 0
    genomes_in_ani: int = 0
    genome_coverage: float = 0.0  # fraction of ANI genomes seen in alignments
    missing_genomes: list[str] = field(default_factory=list)
    single_hit_fraction: float = 0.0
    mean_alignment_length: float = 0.0
    mean_identity: float = 0.0
    total_reads: int = 0

    # Post-classification metrics
    classification_counts: dict[str, int] = field(default_factory=dict)
    low_confidence_fraction: float = 0.0
    ambiguous_fraction: float = 0.0
    novel_fraction: float = 0.0

    # Warnings
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Convert to JSON-serializable dictionary."""
        return {
            "total_alignments": self.total_alignments,
            "filtered_alignments": self.filtered_alignments,
            "filter_rate": round(self.filter_rate, 4),
            "genomes_in_alignment": self.genomes_in_alignment,
            "genomes_in_ani": self.genomes_in_ani,
            "genome_coverage": round(self.genome_coverage, 4),
            "missing_genomes": self.missing_genomes,
            "single_hit_fraction": round(self.single_hit_fraction, 4),
            "mean_alignment_length": round(self.mean_alignment_length, 1),
            "mean_identity": round(self.mean_identity, 2),
            "total_reads": self.total_reads,
            "classification_counts": self.classification_counts,
            "low_confidence_fraction": round(self.low_confidence_fraction, 4),
            "ambiguous_fraction": round(self.ambiguous_fraction, 4),
            "novel_fraction": round(self.novel_fraction, 4),
            "warnings": self.warnings,
        }


def compute_pre_qc(
    raw_df: pl.DataFrame,
    filtered_df: pl.DataFrame,
    ani_matrix: ANIMatrix,
    config: ScoringConfig,
) -> QCMetrics:
    """
    Compute pre-classification QC metrics.

    Compares raw alignments to filtered alignments and checks
    genome coverage against the ANI matrix.

    Args:
        raw_df: Raw BLAST/alignment DataFrame before filtering.
            Must have columns: qseqid, sseqid, pident, length, genome_name.
        filtered_df: DataFrame after alignment quality filters applied.
        ani_matrix: ANI matrix used for classification.
        config: Scoring configuration.

    Returns:
        QCMetrics with pre-classification fields populated.
    """
    qc = QCMetrics()

    qc.total_alignments = raw_df.height
    qc.filtered_alignments = raw_df.height - filtered_df.height
    qc.filter_rate = qc.filtered_alignments / max(1, qc.total_alignments)

    # Genome coverage: how many ANI matrix genomes appear in alignments
    ani_genomes = set(ani_matrix.genomes)
    qc.genomes_in_ani = len(ani_genomes)

    if "genome_name" in filtered_df.columns and not filtered_df.is_empty():
        alignment_genomes = set(filtered_df["genome_name"].unique().to_list())
        qc.genomes_in_alignment = len(alignment_genomes)
        qc.genome_coverage = len(alignment_genomes & ani_genomes) / max(1, len(ani_genomes))
        qc.missing_genomes = sorted(ani_genomes - alignment_genomes)
    else:
        qc.genomes_in_alignment = 0
        qc.genome_coverage = 0.0
        qc.missing_genomes = sorted(ani_genomes)

    # Single-hit fraction
    if not filtered_df.is_empty() and "qseqid" in filtered_df.columns:
        reads_per_hit = filtered_df.group_by("qseqid").len()
        qc.total_reads = reads_per_hit.height
        single_hit_count = reads_per_hit.filter(pl.col("len") == 1).height
        qc.single_hit_fraction = single_hit_count / max(1, qc.total_reads)

    # Alignment statistics
    if not filtered_df.is_empty():
        if "length" in filtered_df.columns:
            qc.mean_alignment_length = filtered_df["length"].mean()
        if "pident" in filtered_df.columns:
            qc.mean_identity = filtered_df["pident"].mean()

    # Generate warnings
    if qc.filter_rate > 0.5:
        qc.warnings.append(
            f"High filter rate ({qc.filter_rate:.0%}): "
            f"over half of alignments were removed by quality filters. "
            f"Consider checking alignment parameters."
        )

    if qc.genome_coverage < 0.5 and qc.genomes_in_ani > 3:
        qc.warnings.append(
            f"Low genome coverage ({qc.genome_coverage:.0%}): "
            f"only {qc.genomes_in_alignment} of {qc.genomes_in_ani} ANI genomes "
            f"have alignments. Missing genomes may affect classification."
        )

    if qc.single_hit_fraction > 0.8:
        qc.warnings.append(
            f"High single-hit fraction ({qc.single_hit_fraction:.0%}): "
            f"most reads align to only one genome. "
            f"Placement uncertainty will be inferred rather than measured for most reads."
        )

    if qc.mean_identity < 80.0 and not filtered_df.is_empty():
        qc.warnings.append(
            f"Low mean identity ({qc.mean_identity:.1f}%): "
            f"reads may be too divergent from reference genomes for reliable classification."
        )

    return qc


def compute_post_qc(
    qc: QCMetrics,
    result_df: pl.DataFrame,
) -> QCMetrics:
    """
    Add post-classification QC metrics to existing QCMetrics.

    Args:
        qc: Pre-populated QCMetrics from compute_pre_qc.
        result_df: Classification result DataFrame.
            Expected columns: taxonomic_call, confidence_score.

    Returns:
        Updated QCMetrics with post-classification fields populated.
    """
    if result_df.is_empty():
        return qc

    # Classification distribution
    if "taxonomic_call" in result_df.columns:
        counts = result_df.group_by("taxonomic_call").len()
        qc.classification_counts = {
            row[0]: row[1] for row in counts.iter_rows()
        }

        total = result_df.height
        ambiguous_count = qc.classification_counts.get("Ambiguous", 0)
        qc.ambiguous_fraction = ambiguous_count / max(1, total)

        novel_species = qc.classification_counts.get("Novel Species", 0)
        novel_genus = qc.classification_counts.get("Novel Genus", 0)
        qc.novel_fraction = (novel_species + novel_genus) / max(1, total)

    # Confidence score analysis
    if "confidence_score" in result_df.columns:
        low_confidence = result_df.filter(pl.col("confidence_score") < 0.5)
        qc.low_confidence_fraction = low_confidence.height / max(1, result_df.height)

    # Post-classification warnings
    if qc.ambiguous_fraction > 0.5:
        qc.warnings.append(
            f"High ambiguous fraction ({qc.ambiguous_fraction:.0%}): "
            f"over half of reads could not be confidently classified. "
            f"This may indicate divergent taxa or insufficient reference coverage."
        )

    if qc.low_confidence_fraction > 0.3:
        qc.warnings.append(
            f"Low confidence reads ({qc.low_confidence_fraction:.0%}): "
            f"many reads have confidence scores below 0.5. "
            f"Results should be interpreted with caution."
        )

    return qc
