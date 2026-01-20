"""
Pydantic data models for metadarkmatter.

Provides type-safe models for BLAST results, classifications,
configuration, and sample metadata.
"""

from metadarkmatter.models.blast import BlastHit, BlastResult
from metadarkmatter.models.classification import (
    ReadClassification,
    TaxonomicCall,
    TaxonomicSummary,
)
from metadarkmatter.models.config import (
    BlastConfig,
    Bowtie2Config,
    GlobalConfig,
    KrakenConfig,
    ScoringConfig,
)
from metadarkmatter.models.genomes import AccessionList, GenomeAccession

__all__ = [
    "AccessionList",
    "BlastConfig",
    "BlastHit",
    "BlastResult",
    "Bowtie2Config",
    "GenomeAccession",
    "GlobalConfig",
    "KrakenConfig",
    "ReadClassification",
    "ScoringConfig",
    "TaxonomicCall",
    "TaxonomicSummary",
]
