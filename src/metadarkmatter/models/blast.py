"""
Pydantic models for BLAST results parsing.

These models represent BLAST tabular output (-outfmt 6) for competitive
read recruitment against reference genome databases.
"""

from __future__ import annotations

import re
from collections.abc import Iterator
from typing import ClassVar, Self

import numpy as np
from pydantic import BaseModel, Field, field_validator, model_validator


class BlastHit(BaseModel):
    """
    Single BLAST alignment hit from tabular output.

    Represents one alignment between a metagenomic read and a reference
    genome sequence. Multiple hits per read are expected in competitive
    recruitment scenarios.

    Attributes:
        qseqid: Query sequence identifier (read ID)
        sseqid: Subject sequence identifier (genome contig ID)
        pident: Percent identity (0-100)
        length: Alignment length in base pairs
        mismatch: Number of mismatches
        gapopen: Number of gap openings
        qstart: Start position in query
        qend: End position in query
        sstart: Start position in subject
        send: End position in subject
        evalue: Expectation value
        bitscore: Bit score
        qlen: Query sequence length (optional, None for backward compatibility)
    """

    qseqid: str = Field(description="Query sequence ID (read ID)")
    sseqid: str = Field(description="Subject sequence ID (genome contig)")
    pident: float = Field(ge=0, description="Percent identity (0-100, clamped)")
    length: int = Field(ge=0, description="Alignment length")
    mismatch: int = Field(ge=0, description="Number of mismatches")
    gapopen: int = Field(ge=0, description="Number of gap openings")
    qstart: int = Field(ge=1, description="Query start position")
    qend: int = Field(ge=1, description="Query end position")
    sstart: int = Field(ge=1, description="Subject start position")
    send: int = Field(ge=1, description="Subject end position")
    evalue: float = Field(ge=0, description="Expectation value")
    bitscore: float = Field(ge=0, description="Bit score")
    qlen: int | None = Field(
        default=None,
        ge=1,
        description="Query sequence length (optional, for accurate coverage calculation)",
    )

    @field_validator("pident", mode="before")
    @classmethod
    def clamp_pident(cls, v: float) -> float:
        """
        Clamp percent identity to valid range [0, 100].

        BLAST can occasionally report pident > 100 due to floating-point
        rounding artifacts in alignment scoring. This validator ensures
        the value is clamped to the valid range.
        """
        if isinstance(v, (int, float)):
            return max(0.0, min(100.0, float(v)))
        return v

    # Genome name extraction pattern for common formats
    # Examples:
    #   GCF_000123456.1_ASM12345v1_genomic -> GCF_000123456.1
    #   NZ_CP012345.1 -> NZ_CP012345.1
    #   contig_1234|genome_name -> genome_name
    GENOME_PATTERN: ClassVar[re.Pattern] = re.compile(
        r"^(?P<genome>GCF_\d+\.\d+|GCA_\d+\.\d+|NZ_[A-Z]+\d+\.\d+)"
        r"|(?P<alt>[^|]+\|(?P<name>[^|]+))"
    )

    @property
    def genome_name(self) -> str:
        """
        Extract genome identifier from subject sequence ID.

        Handles common RefSeq/GenBank formats and custom pipe-delimited formats.
        Falls back to full sseqid if no pattern matches.

        Returns:
            Genome identifier string, or "unknown" if sseqid is empty/invalid
        """
        # Handle empty or whitespace-only sseqid
        if not self.sseqid or not self.sseqid.strip():
            return "unknown"

        match = self.GENOME_PATTERN.match(self.sseqid)
        if match:
            if match.group("genome"):
                return match.group("genome")
            if match.group("name"):
                return match.group("name")

        # Fallback: use everything before first space or full sseqid
        parts = self.sseqid.split()
        if parts:
            return parts[0]

        # Final fallback for edge cases (e.g., only whitespace after strip)
        return "unknown"

    @model_validator(mode="after")
    def validate_query_positions(self) -> Self:
        """Ensure query end >= query start after all fields are set."""
        if self.qend < self.qstart:
            msg = f"qend ({self.qend}) must be >= qstart ({self.qstart})"
            raise ValueError(msg)
        return self

    def calculate_coverage(self, read_length: int | None = None) -> float:
        """
        Calculate alignment coverage as fraction of read length.

        Args:
            read_length: Explicit read length. If None, uses qlen if available,
                         otherwise falls back to qend proxy.

        Returns:
            Coverage fraction between 0.0 and 1.0

        Raises:
            ValueError: If effective read length is <= 0
        """
        # Determine effective read length
        if read_length is not None:
            effective_length = read_length
        elif self.qlen is not None:
            effective_length = self.qlen
        else:
            effective_length = self.qend

        if effective_length <= 0:
            msg = f"Effective read length must be positive, got {effective_length}"
            raise ValueError(msg)

        coverage = (self.qend - self.qstart + 1) / effective_length
        return min(1.0, max(0.0, coverage))

    def calculate_weighted_score(
        self,
        read_length: int | None = None,
        mode: str = "none",
        strength: float = 0.5,
    ) -> float:
        """
        Calculate coverage-weighted bitscore.

        Args:
            read_length: Explicit read length (if None, uses qlen or qend proxy)
            mode: Weighting mode ("none", "linear", "log", "sigmoid")
            strength: Weight strength parameter (0.0-1.0)

        Returns:
            Weighted bitscore (unmodified if mode="none")

        Raises:
            ValueError: If mode is unknown or read_length invalid
        """
        if mode == "none":
            return self.bitscore

        coverage = self.calculate_coverage(read_length)
        weight = _calculate_coverage_weight(coverage, mode, strength)
        return self.bitscore * weight

    model_config = {"frozen": True}

    @classmethod
    def from_blast_line(cls, line: str) -> BlastHit:
        """
        Parse a single line from BLAST tabular output.

        Expected format (outfmt 6):
        qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore [qlen]

        Args:
            line: Tab-delimited BLAST output line

        Returns:
            Parsed BlastHit instance

        Raises:
            ValueError: If line format is invalid (must be 12 or 13 fields)
        """
        fields = line.strip().split("\t")
        if len(fields) not in (12, 13):
            msg = f"Expected 12 or 13 fields in BLAST line, got {len(fields)}"
            raise ValueError(msg)

        return cls(
            qseqid=fields[0],
            sseqid=fields[1],
            pident=float(fields[2]),
            length=int(fields[3]),
            mismatch=int(fields[4]),
            gapopen=int(fields[5]),
            qstart=int(fields[6]),
            qend=int(fields[7]),
            sstart=int(fields[8]),
            send=int(fields[9]),
            evalue=float(fields[10]),
            bitscore=float(fields[11]),
            qlen=int(fields[12]) if len(fields) == 13 else None,
        )


def _calculate_coverage_weight(coverage: float, mode: str, strength: float) -> float:
    """
    Calculate coverage weight for given mode and strength.

    Args:
        coverage: Coverage fraction (0.0-1.0)
        mode: Weighting mode ("linear", "log", "sigmoid")
        strength: Weight strength parameter (0.0-1.0)

    Returns:
        Coverage weight multiplier

    Raises:
        ValueError: If mode is unknown
    """
    min_weight = 1.0 - strength
    max_weight = 1.0 + strength
    weight_range = max_weight - min_weight

    if mode == "linear":
        normalized = coverage
    elif mode == "log":
        normalized = np.log(1 + 9 * coverage) / np.log(10)
    elif mode == "sigmoid":
        normalized = 1.0 / (1.0 + np.exp(-10.0 * (coverage - 0.6)))
    else:
        msg = f"Unknown coverage weight mode: {mode}"
        raise ValueError(msg)

    return min_weight + weight_range * normalized


class BlastResult(BaseModel):
    """
    Collection of BLAST hits for a single read.

    Groups all alignment hits for one query sequence, sorted by bitscore
    in descending order. This model is immutable to ensure consistency
    with the frozen BlastHit model.

    Attributes:
        read_id: Query sequence identifier
        hits: List of BLAST hits sorted by bitscore (highest first)
    """

    read_id: str = Field(description="Query sequence ID")
    hits: tuple[BlastHit, ...] = Field(
        default_factory=tuple,
        description="BLAST hits sorted by bitscore descending",
    )

    @model_validator(mode="after")
    def ensure_sorted(self) -> Self:
        """
        Ensure hits are stored sorted by bitscore descending.

        Performance optimization: Only re-sorts if hits are not already sorted.
        When using StreamingBlastParser, hits are pre-sorted by Polars, so this
        check avoids redundant O(n log n) sorting operations.
        """
        if len(self.hits) <= 1:
            return self

        # Check if already sorted (O(n) check vs O(n log n) sort)
        is_sorted = all(
            self.hits[i].bitscore >= self.hits[i + 1].bitscore
            for i in range(len(self.hits) - 1)
        )

        if not is_sorted:
            sorted_hits = tuple(sorted(self.hits, key=lambda h: h.bitscore, reverse=True))
            # Use object.__setattr__ to bypass frozen restriction during validation
            object.__setattr__(self, "hits", sorted_hits)

        return self

    @property
    def best_hit(self) -> BlastHit | None:
        """Get the highest-scoring BLAST hit."""
        return self.hits[0] if self.hits else None

    @property
    def num_hits(self) -> int:
        """Total number of BLAST hits for this read."""
        return len(self.hits)

    @property
    def top_bitscore(self) -> float:
        """Bitscore of the best hit."""
        return self.best_hit.bitscore if self.best_hit else 0.0

    def iter_ambiguous_hits(self, threshold_pct: float = 95.0) -> Iterator[BlastHit]:
        """
        Iterate over hits within threshold percentage of top bitscore.

        Performance optimization: Uses early termination since hits are sorted
        by bitscore descending. Stops iteration as soon as bitscore drops
        below threshold, avoiding full list traversal.

        Used to identify competing genomes for placement uncertainty calculation.

        Args:
            threshold_pct: Percentage of top bitscore (default: 95%)

        Yields:
            BlastHit objects with bitscore >= (top_bitscore * threshold_pct / 100)
        """
        if not self.hits:
            return

        cutoff = self.top_bitscore * (threshold_pct / 100.0)

        # Hits are sorted by bitscore descending, so we can stop early
        for hit in self.hits:
            if hit.bitscore < cutoff:
                break  # Early termination - no more hits will qualify
            yield hit

    def get_ambiguous_hits(self, threshold_pct: float = 95.0) -> list[BlastHit]:
        """
        Get hits within threshold percentage of top bitscore.

        For memory efficiency with large hit counts, consider using
        iter_ambiguous_hits() instead.

        Args:
            threshold_pct: Percentage of top bitscore (default: 95%)

        Returns:
            List of hits with bitscore >= (top_bitscore * threshold_pct / 100)
        """
        return list(self.iter_ambiguous_hits(threshold_pct))

    def get_best_hit_weighted(
        self,
        read_length: int | None = None,
        mode: str = "none",
        strength: float = 0.5,
    ) -> BlastHit | None:
        """
        Get best hit using coverage-weighted scoring.

        Args:
            read_length: Explicit read length (if None, uses qlen or qend proxy)
            mode: Weighting mode ("none", "linear", "log", "sigmoid")
            strength: Weight strength parameter (0.0-1.0)

        Returns:
            BlastHit with highest weighted score, or None if no hits
        """
        if not self.hits:
            return None

        if mode == "none":
            return self.best_hit

        # Calculate weighted scores for all hits
        weighted_hits = [
            (hit, hit.calculate_weighted_score(read_length, mode, strength))
            for hit in self.hits
        ]

        # Sort by weighted score descending, then by pident descending as tiebreaker
        weighted_hits.sort(key=lambda x: (x[1], x[0].pident), reverse=True)

        return weighted_hits[0][0]

    def iter_ambiguous_hits_weighted(
        self,
        read_length: int | None = None,
        mode: str = "none",
        strength: float = 0.5,
        threshold_pct: float = 95.0,
    ) -> Iterator[BlastHit]:
        """
        Iterate over hits within threshold of top weighted bitscore.

        Args:
            read_length: Explicit read length (if None, uses qlen or qend proxy)
            mode: Weighting mode ("none", "linear", "log", "sigmoid")
            strength: Weight strength parameter (0.0-1.0)
            threshold_pct: Percentage of top weighted bitscore (default: 95%)

        Yields:
            BlastHit objects with weighted score >= (top_weighted * threshold_pct / 100)
        """
        if not self.hits:
            return

        if mode == "none":
            yield from self.iter_ambiguous_hits(threshold_pct)
            return

        best_hit = self.get_best_hit_weighted(read_length, mode, strength)
        if not best_hit:
            return

        top_weighted = best_hit.calculate_weighted_score(read_length, mode, strength)
        cutoff = top_weighted * (threshold_pct / 100.0)

        for hit in self.hits:
            weighted_score = hit.calculate_weighted_score(read_length, mode, strength)
            if weighted_score >= cutoff:
                yield hit

    model_config = {"frozen": True}
