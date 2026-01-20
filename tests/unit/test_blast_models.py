"""
Unit tests for BLAST data models.

Tests BlastHit and BlastResult Pydantic models including:
- Field validation and constraints
- Genome name extraction from various sseqid formats
- Hit sorting and ambiguous hit detection
- Line parsing from BLAST tabular output
"""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from metadarkmatter.models.blast import BlastHit, BlastResult


class TestBlastHit:
    """Tests for BlastHit Pydantic model."""

    def test_create_valid_hit(self, valid_blast_hit_dict):
        """BlastHit should be created from valid dictionary."""
        hit = BlastHit(**valid_blast_hit_dict)

        assert hit.qseqid == "read_001"
        assert hit.sseqid == "GCF_000123456.1_ASM123v1_genomic"
        assert hit.pident == 98.5
        assert hit.bitscore == 250.0

    def test_hit_is_frozen(self, valid_blast_hit_dict):
        """BlastHit should be immutable (frozen)."""
        hit = BlastHit(**valid_blast_hit_dict)

        with pytest.raises(ValidationError):
            hit.pident = 50.0

    def test_pident_clamped_above_100(self, valid_blast_hit_dict):
        """Percent identity > 100 should be clamped to 100."""
        valid_blast_hit_dict["pident"] = 105.0
        hit = BlastHit(**valid_blast_hit_dict)

        assert hit.pident == 100.0

    def test_pident_clamped_below_zero(self, valid_blast_hit_dict):
        """Percent identity < 0 should be clamped to 0."""
        valid_blast_hit_dict["pident"] = -5.0
        hit = BlastHit(**valid_blast_hit_dict)

        assert hit.pident == 0.0

    def test_pident_at_boundary_100(self, valid_blast_hit_dict):
        """Percent identity at exactly 100 should remain 100."""
        valid_blast_hit_dict["pident"] = 100.0
        hit = BlastHit(**valid_blast_hit_dict)

        assert hit.pident == 100.0

    def test_negative_bitscore_rejected(self, valid_blast_hit_dict):
        """Negative bitscore should raise validation error."""
        valid_blast_hit_dict["bitscore"] = -10.0

        with pytest.raises(ValidationError) as exc_info:
            BlastHit(**valid_blast_hit_dict)

        assert "bitscore" in str(exc_info.value)

    def test_negative_length_rejected(self, valid_blast_hit_dict):
        """Negative alignment length should raise validation error."""
        valid_blast_hit_dict["length"] = -1

        with pytest.raises(ValidationError) as exc_info:
            BlastHit(**valid_blast_hit_dict)

        assert "length" in str(exc_info.value)

    def test_qend_less_than_qstart_rejected(self, valid_blast_hit_dict):
        """Query end < query start should raise validation error."""
        valid_blast_hit_dict["qstart"] = 150
        valid_blast_hit_dict["qend"] = 1

        with pytest.raises(ValidationError) as exc_info:
            BlastHit(**valid_blast_hit_dict)

        assert "qend" in str(exc_info.value).lower() or "qstart" in str(exc_info.value).lower()

    def test_qstart_zero_rejected(self, valid_blast_hit_dict):
        """Query start position of 0 should raise validation error (1-based)."""
        valid_blast_hit_dict["qstart"] = 0

        with pytest.raises(ValidationError) as exc_info:
            BlastHit(**valid_blast_hit_dict)

        assert "qstart" in str(exc_info.value)


class TestBlastHitGenomeName:
    """Tests for genome name extraction from sseqid."""

    def test_refseq_format_gcf(self, valid_blast_hit_dict):
        """Should extract GCF accession from RefSeq format."""
        valid_blast_hit_dict["sseqid"] = "GCF_000123456.1_ASM123v1_genomic"
        hit = BlastHit(**valid_blast_hit_dict)

        assert hit.genome_name == "GCF_000123456.1"

    def test_refseq_format_gcf_version_2(self, valid_blast_hit_dict):
        """Should extract GCF accession with version 2."""
        valid_blast_hit_dict["sseqid"] = "GCF_000123456.2_ASM123v2_genomic"
        hit = BlastHit(**valid_blast_hit_dict)

        assert hit.genome_name == "GCF_000123456.2"

    def test_genbank_format_gca(self, valid_blast_hit_dict):
        """Should extract GCA accession from GenBank format."""
        valid_blast_hit_dict["sseqid"] = "GCA_000789012.1_ASM789v1_genomic"
        hit = BlastHit(**valid_blast_hit_dict)

        assert hit.genome_name == "GCA_000789012.1"

    def test_refseq_nz_prefix(self, valid_blast_hit_dict):
        """Should extract NZ_ prefixed accession."""
        valid_blast_hit_dict["sseqid"] = "NZ_CP012345.1"
        hit = BlastHit(**valid_blast_hit_dict)

        assert hit.genome_name == "NZ_CP012345.1"

    def test_refseq_nz_wgs(self, valid_blast_hit_dict):
        """Should extract NZ_ WGS accession format."""
        valid_blast_hit_dict["sseqid"] = "NZ_ABCD01000001.1"
        hit = BlastHit(**valid_blast_hit_dict)

        assert hit.genome_name == "NZ_ABCD01000001.1"

    def test_simple_accession_only(self, valid_blast_hit_dict):
        """Should return full accession when no suffix."""
        valid_blast_hit_dict["sseqid"] = "GCF_000123456.1"
        hit = BlastHit(**valid_blast_hit_dict)

        assert hit.genome_name == "GCF_000123456.1"

    def test_custom_format_fallback(self, valid_blast_hit_dict):
        """Should return full sseqid for non-standard formats."""
        valid_blast_hit_dict["sseqid"] = "custom_genome_contig_1"
        hit = BlastHit(**valid_blast_hit_dict)

        assert hit.genome_name == "custom_genome_contig_1"

    def test_sseqid_with_spaces(self, valid_blast_hit_dict):
        """Should return first part before space."""
        valid_blast_hit_dict["sseqid"] = "genome_name additional_info"
        hit = BlastHit(**valid_blast_hit_dict)

        assert hit.genome_name == "genome_name"

    def test_empty_sseqid(self, valid_blast_hit_dict):
        """Empty sseqid should return 'unknown'."""
        valid_blast_hit_dict["sseqid"] = ""
        hit = BlastHit(**valid_blast_hit_dict)

        assert hit.genome_name == "unknown"

    def test_whitespace_only_sseqid(self, valid_blast_hit_dict):
        """Whitespace-only sseqid should return 'unknown'."""
        valid_blast_hit_dict["sseqid"] = "   "
        hit = BlastHit(**valid_blast_hit_dict)

        assert hit.genome_name == "unknown"

    def test_pipe_delimited_format(self, valid_blast_hit_dict):
        """Should extract genome name from pipe-delimited format."""
        valid_blast_hit_dict["sseqid"] = "contig_1234|my_genome_name"
        hit = BlastHit(**valid_blast_hit_dict)

        # Based on regex pattern, should extract name after pipe
        assert hit.genome_name == "my_genome_name"


class TestBlastHitFromLine:
    """Tests for parsing BLAST tabular output lines."""

    def test_parse_valid_line(self, valid_blast_line):
        """Should parse valid BLAST line correctly."""
        hit = BlastHit.from_blast_line(valid_blast_line)

        assert hit.qseqid == "read_001"
        assert hit.sseqid == "GCF_000123456.1_ASM123v1_genomic"
        assert hit.pident == 98.5
        assert hit.length == 150
        assert hit.mismatch == 2
        assert hit.gapopen == 0
        assert hit.qstart == 1
        assert hit.qend == 150
        assert hit.sstart == 1000
        assert hit.send == 1150
        assert hit.evalue == 1e-50
        assert hit.bitscore == 250.0

    def test_parse_line_with_trailing_whitespace(self):
        """Should handle trailing whitespace."""
        line = "read_001\tGCF_000123456.1\t98.5\t150\t2\t0\t1\t150\t1000\t1150\t1e-50\t250.0\n"
        hit = BlastHit.from_blast_line(line)

        assert hit.qseqid == "read_001"
        assert hit.bitscore == 250.0

    def test_parse_line_too_few_fields(self):
        """Should raise ValueError for lines with too few fields."""
        line = "read_001\tGCF_000123456.1\t98.5\t150"

        with pytest.raises(ValueError) as exc_info:
            BlastHit.from_blast_line(line)

        assert "12 fields" in str(exc_info.value)

    def test_parse_line_too_many_fields(self):
        """Should raise ValueError for lines with too many fields."""
        line = "read_001\tGCF_000123456.1\t98.5\t150\t2\t0\t1\t150\t1000\t1150\t1e-50\t250.0\textra"

        with pytest.raises(ValueError) as exc_info:
            BlastHit.from_blast_line(line)

        assert "12 fields" in str(exc_info.value)

    def test_parse_line_invalid_number(self):
        """Should raise error for non-numeric fields."""
        line = "read_001\tGCF_000123456.1\tnot_a_number\t150\t2\t0\t1\t150\t1000\t1150\t1e-50\t250.0"

        with pytest.raises(ValueError):
            BlastHit.from_blast_line(line)

    def test_parse_scientific_notation_evalue(self):
        """Should correctly parse scientific notation evalue."""
        line = "read_001\tGCF_000123456.1\t98.5\t150\t2\t0\t1\t150\t1000\t1150\t1.5e-100\t250.0"
        hit = BlastHit.from_blast_line(line)

        assert hit.evalue == 1.5e-100


class TestBlastResult:
    """Tests for BlastResult collection model."""

    def test_create_empty_result(self):
        """Should create BlastResult with no hits."""
        result = BlastResult(read_id="read_001", hits=())

        assert result.read_id == "read_001"
        assert result.num_hits == 0
        assert result.best_hit is None

    def test_create_result_with_hits(self, multi_hit_blast_data):
        """Should create BlastResult with multiple hits."""
        hits = tuple(BlastHit(**d) for d in multi_hit_blast_data)
        result = BlastResult(read_id="read_001", hits=hits)

        assert result.num_hits == 3

    def test_hits_sorted_by_bitscore(self, multi_hit_blast_data):
        """Hits should be sorted by bitscore descending."""
        # Provide hits in wrong order
        hits = tuple(BlastHit(**d) for d in reversed(multi_hit_blast_data))
        result = BlastResult(read_id="read_001", hits=hits)

        # Best hit should have highest bitscore
        assert result.best_hit.bitscore == 250.0
        assert result.hits[0].bitscore >= result.hits[1].bitscore >= result.hits[2].bitscore

    def test_best_hit_property(self, multi_hit_blast_data):
        """best_hit should return highest bitscore hit."""
        hits = tuple(BlastHit(**d) for d in multi_hit_blast_data)
        result = BlastResult(read_id="read_001", hits=hits)

        assert result.best_hit.bitscore == 250.0
        assert result.best_hit.pident == 98.5

    def test_top_bitscore_property(self, multi_hit_blast_data):
        """top_bitscore should return bitscore of best hit."""
        hits = tuple(BlastHit(**d) for d in multi_hit_blast_data)
        result = BlastResult(read_id="read_001", hits=hits)

        assert result.top_bitscore == 250.0

    def test_top_bitscore_empty_result(self):
        """top_bitscore should return 0.0 for empty result."""
        result = BlastResult(read_id="read_001", hits=())

        assert result.top_bitscore == 0.0

    def test_result_is_frozen(self, multi_hit_blast_data):
        """BlastResult should be immutable."""
        hits = tuple(BlastHit(**d) for d in multi_hit_blast_data)
        result = BlastResult(read_id="read_001", hits=hits)

        with pytest.raises(ValidationError):
            result.read_id = "different_read"


class TestBlastResultAmbiguousHits:
    """Tests for ambiguous hit detection."""

    def test_get_ambiguous_hits_default_threshold(self, multi_hit_blast_data):
        """Should return hits within 95% of top bitscore."""
        hits = tuple(BlastHit(**d) for d in multi_hit_blast_data)
        result = BlastResult(read_id="read_001", hits=hits)

        # Top bitscore: 250.0, threshold at 95%: 237.5
        # Hit with 240.0 should be included, 180.0 should not
        ambiguous = result.get_ambiguous_hits(threshold_pct=95.0)

        assert len(ambiguous) == 2
        bitscores = [h.bitscore for h in ambiguous]
        assert 250.0 in bitscores
        assert 240.0 in bitscores
        assert 180.0 not in bitscores

    def test_get_ambiguous_hits_custom_threshold(self, multi_hit_blast_data):
        """Should respect custom threshold percentage."""
        hits = tuple(BlastHit(**d) for d in multi_hit_blast_data)
        result = BlastResult(read_id="read_001", hits=hits)

        # With 70% threshold: 250.0 * 0.7 = 175.0
        # All three hits should be included
        ambiguous = result.get_ambiguous_hits(threshold_pct=70.0)

        assert len(ambiguous) == 3

    def test_get_ambiguous_hits_strict_threshold(self, multi_hit_blast_data):
        """Very strict threshold should return only best hit."""
        hits = tuple(BlastHit(**d) for d in multi_hit_blast_data)
        result = BlastResult(read_id="read_001", hits=hits)

        # With 99.9% threshold, only exact matches
        ambiguous = result.get_ambiguous_hits(threshold_pct=99.9)

        assert len(ambiguous) == 1
        assert ambiguous[0].bitscore == 250.0

    def test_iter_ambiguous_hits_empty_result(self):
        """iter_ambiguous_hits should yield nothing for empty result."""
        result = BlastResult(read_id="read_001", hits=())

        ambiguous = list(result.iter_ambiguous_hits())

        assert ambiguous == []

    def test_iter_ambiguous_hits_single_hit(self, valid_blast_hit_dict):
        """Single hit should be returned as ambiguous."""
        hit = BlastHit(**valid_blast_hit_dict)
        result = BlastResult(read_id="read_001", hits=(hit,))

        ambiguous = list(result.iter_ambiguous_hits())

        assert len(ambiguous) == 1
        assert ambiguous[0] == hit

    def test_iter_ambiguous_hits_early_termination(self):
        """Iterator should stop when bitscore drops below threshold."""
        # Create many hits with decreasing bitscores
        hits = []
        for i in range(100):
            hits.append(
                BlastHit(
                    qseqid="read_001",
                    sseqid=f"GCF_{i:09d}.1",
                    pident=95.0,
                    length=150,
                    mismatch=7,
                    gapopen=0,
                    qstart=1,
                    qend=150,
                    sstart=1000,
                    send=1150,
                    evalue=1e-50,
                    bitscore=250.0 - i,  # Decreasing
                )
            )

        result = BlastResult(read_id="read_001", hits=tuple(hits))

        # With 95% threshold, only first few should qualify
        ambiguous = list(result.iter_ambiguous_hits(threshold_pct=95.0))

        # 250 * 0.95 = 237.5, so hits with bitscore >= 238 qualify
        # That's hits 0-12 (bitscores 250-238)
        assert len(ambiguous) == 13

    def test_sorted_by_bitscore_returns_self(self, multi_hit_blast_data):
        """sorted_by_bitscore should return self (already sorted)."""
        hits = tuple(BlastHit(**d) for d in multi_hit_blast_data)
        result = BlastResult(read_id="read_001", hits=hits)

        sorted_result = result.sorted_by_bitscore()

        assert sorted_result is result
