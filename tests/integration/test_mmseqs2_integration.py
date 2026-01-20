"""
Integration tests for MMseqs2 wrapper.

These tests require MMseqs2 to be installed and will be skipped if not available.
They validate real MMseqs2 execution and output format compatibility.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from metadarkmatter.external.mmseqs2 import MMseqs2


@pytest.fixture
def test_fasta(tmp_path: Path) -> Path:
    """Create a small test FASTA file."""
    fasta_path = tmp_path / "test_genomes.fasta"
    fasta_content = """>GCF_000001|contig1
ATCGATCGATCGATCGATCGATCG
>GCF_000002|contig1
GCTAGCTAGCTAGCTAGCTAGCTA
>GCF_000003|contig1
TTAATTAATTAATTAATTAATTAA
"""
    fasta_path.write_text(fasta_content)
    return fasta_path


@pytest.fixture
def test_query(tmp_path: Path) -> Path:
    """Create a small test query file."""
    query_path = tmp_path / "test_query.fasta"
    query_content = """>read1
ATCGATCGATCGATCG
>read2
GCTAGCTAGCTAGCTA
>read3
TTAATTAATTAATTAA
"""
    query_path.write_text(query_content)
    return query_path


@pytest.mark.requires_mmseqs2
class TestMMseqs2Integration:
    """Integration tests requiring actual MMseqs2 installation."""

    def test_mmseqs2_available(self):
        """Should detect MMseqs2 installation."""
        assert MMseqs2.check_available(), "MMseqs2 not found in PATH"

    def test_create_database(self, tmp_path: Path, test_fasta: Path):
        """Should create MMseqs2 database from FASTA."""
        mmseqs = MMseqs2()
        db_path = tmp_path / "test_db"

        # Create database
        mmseqs.create_database(
            input_fasta=test_fasta,
            database=db_path,
            dbtype=1,  # nucleotide
            timeout=60.0,
        )

        # Verify database files exist
        assert db_path.exists(), "Main database file should exist"
        assert db_path.with_suffix(".dbtype").exists(), "dbtype file should exist"

    def test_search_output_format(self, tmp_path: Path, test_fasta: Path, test_query: Path):
        """Should produce BLAST-compatible tabular output."""
        mmseqs = MMseqs2()

        # Create database
        db_path = tmp_path / "test_db"
        mmseqs.create_database(
            input_fasta=test_fasta,
            database=db_path,
            dbtype=1,
            timeout=60.0,
        )

        # Run search
        output_path = tmp_path / "results.tsv"
        mmseqs.search(
            query=test_query,
            database=db_path,
            output=output_path,
            threads=1,
            sensitivity=5.7,
            evalue=1e-3,
            timeout=60.0,
        )

        # Verify output exists
        assert output_path.exists(), "Output file should exist"

        # Verify output format
        lines = output_path.read_text().strip().split("\n")
        assert len(lines) > 0, "Should have at least one alignment"

        # Check first line has 12 columns (BLAST format)
        first_line = lines[0]
        columns = first_line.split("\t")
        assert len(columns) == 12, f"Expected 12 columns, got {len(columns)}"

        # Validate column types
        qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = columns

        # Query and subject IDs should be strings
        assert qseqid.startswith("read"), f"Query ID should start with 'read', got {qseqid}"
        assert sseqid.startswith("GCF_"), f"Subject ID should start with 'GCF_', got {sseqid}"

        # Numeric fields should be parseable
        float(pident)  # percent identity (0-100)
        int(length)    # alignment length
        int(mismatch)  # number of mismatches
        int(gapopen)   # number of gap openings
        int(qstart)    # query start position
        int(qend)      # query end position
        int(sstart)    # subject start position
        int(send)      # subject end position
        float(evalue)  # E-value
        float(bitscore)  # bit score

        # Validate ranges
        assert 0 <= float(pident) <= 100, "pident should be 0-100"
        assert float(bitscore) > 0, "bitscore should be positive"

    def test_parser_compatibility(self, tmp_path: Path, test_fasta: Path, test_query: Path):
        """Should be parseable by StreamingBlastParser."""
        from metadarkmatter.core.parsers import StreamingBlastParser

        mmseqs = MMseqs2()

        # Create database and run search
        db_path = tmp_path / "test_db"
        mmseqs.create_database(
            input_fasta=test_fasta,
            database=db_path,
            dbtype=1,
            timeout=60.0,
        )

        output_path = tmp_path / "results.tsv"
        mmseqs.search(
            query=test_query,
            database=db_path,
            output=output_path,
            threads=1,
            sensitivity=5.7,
            timeout=60.0,
        )

        # Parse with StreamingBlastParser
        parser = StreamingBlastParser(output_path)

        # Should parse without errors
        results = list(parser.iter_reads())

        # Should have results
        assert len(results) > 0, "Parser should extract results"

        # Validate first result structure
        result = results[0]
        assert result.read_id.startswith("read"), "Should have read ID"
        assert result.num_hits > 0, "Should have hits"

        # Validate hit structure
        hit = result.best_hit
        assert hit.qseqid.startswith("read"), "Hit should have query ID"
        assert hit.sseqid.startswith("GCF_"), "Hit should have subject ID"
        assert 0 <= hit.pident <= 100, "Percent identity in valid range"
        assert hit.bitscore > 0, "Bitscore should be positive"
        assert hit.genome_name.startswith("GCF_"), "Genome name should be extracted"

    def test_compressed_output(self, tmp_path: Path, test_fasta: Path, test_query: Path):
        """Should handle gzipped output correctly."""
        mmseqs = MMseqs2()

        # Create database
        db_path = tmp_path / "test_db"
        mmseqs.create_database(
            input_fasta=test_fasta,
            database=db_path,
            dbtype=1,
            timeout=60.0,
        )

        # Run search with gzipped output
        output_path = tmp_path / "results.tsv.gz"
        mmseqs.search(
            query=test_query,
            database=db_path,
            output=output_path,
            threads=1,
            timeout=60.0,
        )

        # Verify compressed output exists
        assert output_path.exists(), "Compressed output should exist"

        # Verify can decompress and parse
        with gzip.open(output_path, "rt") as f:
            lines = f.read().strip().split("\n")
            assert len(lines) > 0, "Should have results in compressed file"
            columns = lines[0].split("\t")
            assert len(columns) == 12, "Should have 12 columns"

    def test_sensitivity_parameter(self, tmp_path: Path, test_fasta: Path, test_query: Path):
        """Should respect sensitivity parameter."""
        mmseqs = MMseqs2()

        # Create database
        db_path = tmp_path / "test_db"
        mmseqs.create_database(
            input_fasta=test_fasta,
            database=db_path,
            dbtype=1,
            timeout=60.0,
        )

        # Run with low sensitivity (fast, may miss some hits)
        output_low = tmp_path / "results_low_sens.tsv"
        mmseqs.search(
            query=test_query,
            database=db_path,
            output=output_low,
            sensitivity=1.0,  # Lowest
            threads=1,
            timeout=60.0,
        )

        # Run with high sensitivity (slow, finds more hits)
        output_high = tmp_path / "results_high_sens.tsv"
        mmseqs.search(
            query=test_query,
            database=db_path,
            output=output_high,
            sensitivity=7.0,  # High
            threads=1,
            timeout=60.0,
        )

        # Both should produce valid output
        assert output_low.exists(), "Low sensitivity output should exist"
        assert output_high.exists(), "High sensitivity output should exist"

        # High sensitivity should find at least as many hits
        low_count = len(output_low.read_text().strip().split("\n"))
        high_count = len(output_high.read_text().strip().split("\n"))

        # Note: This may not always be true for very small test data,
        # but documents expected behavior
        assert high_count >= low_count, "Higher sensitivity should find more or equal hits"

    def test_min_identity_filter(self, tmp_path: Path, test_fasta: Path, test_query: Path):
        """Should filter results by minimum identity."""
        mmseqs = MMseqs2()

        # Create database
        db_path = tmp_path / "test_db"
        mmseqs.create_database(
            input_fasta=test_fasta,
            database=db_path,
            dbtype=1,
            timeout=60.0,
        )

        # Run with no identity filter
        output_all = tmp_path / "results_all.tsv"
        mmseqs.search(
            query=test_query,
            database=db_path,
            output=output_all,
            min_identity=None,
            threads=1,
            timeout=60.0,
        )

        # Run with high identity filter
        output_filtered = tmp_path / "results_filtered.tsv"
        mmseqs.search(
            query=test_query,
            database=db_path,
            output=output_filtered,
            min_identity=95.0,  # Very strict
            threads=1,
            timeout=60.0,
        )

        # Both should exist
        assert output_all.exists(), "Unfiltered output should exist"
        assert output_filtered.exists(), "Filtered output should exist"

        # Filtered should have fewer or equal results
        all_count = len(output_all.read_text().strip().split("\n"))
        filtered_count = len(output_filtered.read_text().strip().split("\n"))

        assert filtered_count <= all_count, "Filtering should reduce or maintain result count"

        # All filtered results should meet identity threshold
        if filtered_count > 0:
            for line in output_filtered.read_text().strip().split("\n"):
                columns = line.split("\t")
                pident = float(columns[2])
                assert pident >= 95.0, f"Identity {pident} should be >= 95.0"


@pytest.mark.requires_mmseqs2
class TestMMseqs2vsBlastComparison:
    """Comparative tests between MMseqs2 and BLAST."""

    def test_similar_results_to_blast(self, tmp_path: Path, test_fasta: Path, test_query: Path):
        """Should produce similar results to BLAST (if both installed)."""
        # Skip if BLAST not available
        from metadarkmatter.external.blast import BlastN, MakeBlastDb

        if not MakeBlastDb.check_available() or not BlastN.check_available():
            pytest.skip("BLAST not available for comparison")

        # Create BLAST database
        blast_db = tmp_path / "blast_db"
        blast_builder = MakeBlastDb()
        blast_builder.run_or_raise(
            input_fasta=test_fasta,
            output_db=blast_db,
            dbtype="nucl",
            parse_seqids=True,
            timeout=60.0,
        )

        # Run BLAST
        blast_output = tmp_path / "blast_results.tsv"
        blastn = BlastN()
        blastn.run_or_raise(
            query=test_query,
            database=blast_db,
            output=blast_output,
            outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
            word_size=7,
            evalue=1e-3,
            max_target_seqs=500,
            threads=1,
            timeout=60.0,
        )

        # Create MMseqs2 database
        mmseqs_db = tmp_path / "mmseqs_db"
        mmseqs = MMseqs2()
        mmseqs.create_database(
            input_fasta=test_fasta,
            database=mmseqs_db,
            dbtype=1,
            timeout=60.0,
        )

        # Run MMseqs2
        mmseqs_output = tmp_path / "mmseqs_results.tsv"
        mmseqs.search(
            query=test_query,
            database=mmseqs_db,
            output=mmseqs_output,
            sensitivity=5.7,
            evalue=1e-3,
            max_seqs=500,
            threads=1,
            timeout=60.0,
        )

        # Compare result counts
        blast_count = len(blast_output.read_text().strip().split("\n"))
        mmseqs_count = len(mmseqs_output.read_text().strip().split("\n"))

        # Should have similar number of hits (within 20%)
        ratio = min(blast_count, mmseqs_count) / max(blast_count, mmseqs_count)
        assert ratio >= 0.8, f"Result counts should be similar: BLAST={blast_count}, MMseqs2={mmseqs_count}"


# Pytest configuration for conditional test execution
def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line(
        "markers",
        "requires_mmseqs2: marks tests as requiring MMseqs2 installation",
    )


def pytest_collection_modifyitems(config, items):
    """Skip tests marked requires_mmseqs2 if MMseqs2 not available."""
    if not MMseqs2.check_available():
        skip_mmseqs2 = pytest.mark.skip(reason="MMseqs2 not installed")
        for item in items:
            if "requires_mmseqs2" in item.keywords:
                item.add_marker(skip_mmseqs2)
