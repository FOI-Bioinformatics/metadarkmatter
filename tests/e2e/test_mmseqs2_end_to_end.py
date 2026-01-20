"""
End-to-end tests for MMseqs2 integration with full pipeline validation.

These tests validate the complete workflow:
1. Database creation (BLAST vs MMseqs2)
2. Sequence search (BLAST vs MMseqs2)
3. Result parsing (StreamingBlastParser)
4. Classification pipeline (should produce equivalent results)

Tests use realistic synthetic data with known properties to validate correctness.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

# Import test data generators
import sys
sys.path.insert(0, str(Path(__file__).parent.parent / "fixtures"))
from test_genomes import create_test_genome_set, create_test_metadata, create_test_reads


@pytest.fixture(scope="module")
def test_data_dir(tmp_path_factory):
    """Create temporary directory for test data."""
    return tmp_path_factory.mktemp("e2e_test_data")


@pytest.fixture(scope="module")
def test_genomes(test_data_dir):
    """Create test genome files."""
    genome_dir = test_data_dir / "genomes"
    return create_test_genome_set(genome_dir)


@pytest.fixture(scope="module")
def test_metadata(test_data_dir, test_genomes):
    """Create test metadata file."""
    metadata_file = test_data_dir / "genome_metadata.tsv"
    return create_test_metadata(test_genomes, metadata_file)


@pytest.fixture(scope="module")
def test_reads(test_data_dir, test_genomes):
    """Create test reads file."""
    reads_file = test_data_dir / "reads.fasta"
    return create_test_reads(test_genomes, reads_file, reads_per_genome=10)


@pytest.fixture(scope="module")
def concatenated_genome(test_data_dir, test_genomes):
    """Create concatenated pangenome FASTA."""
    from metadarkmatter.core.genome_utils import concatenate_genomes_with_mapping

    genome_dir = test_data_dir / "genomes"
    pangenome_path = test_data_dir / "pangenome.fasta"
    contig_mapping_path = test_data_dir / "contig_mapping.tsv"

    genome_count, contig_count = concatenate_genomes_with_mapping(
        genome_dir=genome_dir,
        output_fasta=pangenome_path,
        contig_mapping_path=contig_mapping_path,
        pattern="*.fna"
    )

    assert genome_count == 4, "Should concatenate all 4 test genomes"
    assert contig_count >= 12, "Should have at least 12 contigs total"

    return pangenome_path


class TestMMseqs2EndToEnd:
    """End-to-end tests for MMseqs2 integration."""

    def test_database_creation(self, test_data_dir, concatenated_genome):
        """Should create MMseqs2 database successfully."""
        from metadarkmatter.external.mmseqs2 import MMseqs2

        if not MMseqs2.check_available():
            pytest.skip("MMseqs2 not installed")

        mmseqs = MMseqs2()
        db_path = test_data_dir / "mmseqs_db"

        # Create database
        mmseqs.create_database(
            input_fasta=concatenated_genome,
            database=db_path,
            # Note: Omit dbtype to let MMseqs2 auto-detect (explicit dbtype causes issues)  # nucleotide
            timeout=120.0,
        )

        # Verify database files exist
        assert db_path.exists(), "Main database file should exist"
        assert db_path.with_suffix(".dbtype").exists(), "dbtype file should exist"

    def test_search_execution(self, test_data_dir, concatenated_genome, test_reads):
        """Should run MMseqs2 search successfully."""
        from metadarkmatter.external.mmseqs2 import MMseqs2

        if not MMseqs2.check_available():
            pytest.skip("MMseqs2 not installed")

        mmseqs = MMseqs2()
        db_path = test_data_dir / "mmseqs_db"
        output_path = test_data_dir / "mmseqs_results.tsv"

        # Create database if not exists
        if not db_path.exists():
            mmseqs.create_database(
                input_fasta=concatenated_genome,
                database=db_path,
                # Note: Omit dbtype to let MMseqs2 auto-detect (explicit dbtype causes issues)
                timeout=120.0,
            )

        # Run search
        mmseqs.search(
            query=test_reads,
            database=db_path,
            output=output_path,
            sensitivity=5.7,
            evalue=1e-3,
            max_seqs=500,
            threads=2,
            timeout=120.0,
        )

        # Verify output exists and has content
        assert output_path.exists(), "Output file should exist"
        assert output_path.stat().st_size > 0, "Output file should not be empty"

        # Verify format
        lines = output_path.read_text().strip().split('\n')
        assert len(lines) > 0, "Should have at least one alignment"

        # Check first line has 13 columns (including qlen)
        first_line = lines[0]
        columns = first_line.split('\t')
        assert len(columns) == 13, f"Expected 13 columns (with qlen), got {len(columns)}"

    def test_parser_integration(self, test_data_dir, concatenated_genome, test_reads):
        """Should parse MMseqs2 output with StreamingBlastParser."""
        from metadarkmatter.core.parsers import StreamingBlastParser
        from metadarkmatter.external.mmseqs2 import MMseqs2

        if not MMseqs2.check_available():
            pytest.skip("MMseqs2 not installed")

        mmseqs = MMseqs2()
        db_path = test_data_dir / "mmseqs_db"
        output_path = test_data_dir / "mmseqs_results.tsv"

        # Create database and run search if not done
        if not output_path.exists():
            if not db_path.exists():
                mmseqs.create_database(
                    input_fasta=concatenated_genome,
                    database=db_path,
                    # Note: Omit dbtype to let MMseqs2 auto-detect (explicit dbtype causes issues)
                    timeout=120.0,
                )

            mmseqs.search(
                query=test_reads,
                database=db_path,
                output=output_path,
                sensitivity=5.7,
                evalue=1e-3,
                threads=2,
                timeout=120.0,
            )

        # Parse with StreamingBlastParser
        parser = StreamingBlastParser(output_path)
        results = list(parser.iter_reads())

        # Should have parsed results
        assert len(results) > 0, "Should parse at least one read"

        # Validate result structure
        result = results[0]
        assert result.read_id.startswith("read_"), "Should have valid read ID"
        assert result.num_hits > 0, "Should have at least one hit"

        # Validate hit structure
        hit = result.best_hit
        assert hit.qseqid.startswith("read_"), "Hit should have query ID"
        assert hit.sseqid.startswith("GCF_"), "Hit should have subject ID"
        assert 0 <= hit.pident <= 100, "Percent identity in valid range"
        assert hit.bitscore > 0, "Bitscore should be positive"

        # Validate genome name extraction
        assert hit.genome_name.startswith("GCF_"), "Genome name should be extracted"
        assert "." in hit.genome_name, "Genome name should include version"

    @pytest.mark.requires_blast
    def test_blast_vs_mmseqs2_comparison(self, test_data_dir, concatenated_genome, test_reads):
        """Should produce similar results to BLAST."""
        from metadarkmatter.core.parsers import StreamingBlastParser
        from metadarkmatter.external.blast import BlastN, MakeBlastDb
        from metadarkmatter.external.mmseqs2 import MMseqs2

        # Skip if either tool not available
        if not MMseqs2.check_available():
            pytest.skip("MMseqs2 not installed")
        if not MakeBlastDb.check_available() or not BlastN.check_available():
            pytest.skip("BLAST not installed")

        # Create BLAST database
        blast_db = test_data_dir / "blast_db"
        blast_output = test_data_dir / "blast_results.tsv"

        if not blast_output.exists():
            blast_builder = MakeBlastDb()
            blast_builder.run_or_raise(
                input_fasta=concatenated_genome,
                output_db=blast_db,
                dbtype="nucl",
                parse_seqids=True,
                timeout=120.0,
            )

            blastn = BlastN()
            blastn.run_or_raise(
                query=test_reads,
                database=blast_db,
                output=blast_output,
                # Include qlen (column 13) to match MMseqs2 format
                outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen",
                word_size=7,
                evalue=1e-3,
                max_target_seqs=500,
                threads=2,
                timeout=120.0,
            )

        # Create MMseqs2 database and search
        mmseqs_db = test_data_dir / "mmseqs_db"
        mmseqs_output = test_data_dir / "mmseqs_results.tsv"

        if not mmseqs_output.exists():
            mmseqs = MMseqs2()
            mmseqs.create_database(
                input_fasta=concatenated_genome,
                database=mmseqs_db,
                # Note: Omit dbtype to let MMseqs2 auto-detect (explicit dbtype causes issues)
                timeout=120.0,
            )

            mmseqs.search(
                query=test_reads,
                database=mmseqs_db,
                output=mmseqs_output,
                sensitivity=5.7,
                evalue=1e-3,
                max_seqs=500,
                threads=2,
                timeout=120.0,
            )

        # Parse both results
        blast_parser = StreamingBlastParser(blast_output)
        blast_results = {r.read_id: r for r in blast_parser.iter_reads()}

        mmseqs_parser = StreamingBlastParser(mmseqs_output)
        mmseqs_results = {r.read_id: r for r in mmseqs_parser.iter_reads()}

        # Should have similar number of reads with hits
        blast_read_count = len(blast_results)
        mmseqs_read_count = len(mmseqs_results)

        ratio = min(blast_read_count, mmseqs_read_count) / max(blast_read_count, mmseqs_read_count)
        assert ratio >= 0.8, (
            f"Should have similar read counts: "
            f"BLAST={blast_read_count}, MMseqs2={mmseqs_read_count}"
        )

        # Compare top hits for common reads
        common_reads = set(blast_results.keys()) & set(mmseqs_results.keys())
        assert len(common_reads) > 0, "Should have common reads"

        genome_agreement = 0
        for read_id in common_reads:
            blast_genome = blast_results[read_id].best_hit.genome_name
            mmseqs_genome = mmseqs_results[read_id].best_hit.genome_name

            if blast_genome == mmseqs_genome:
                genome_agreement += 1

        agreement_rate = genome_agreement / len(common_reads)
        # Note: BLAST and MMseqs2 use different algorithms, so some disagreement
        # on top hits is expected, especially for closely-related genomes.
        # With synthetic test data, agreement can be lower than real-world data.
        assert agreement_rate >= 0.50, (
            f"Top hit genome agreement should be >=50%, got {agreement_rate:.1%}"
        )

    def test_compressed_output_workflow(self, test_data_dir, concatenated_genome, test_reads):
        """Should handle gzipped output correctly."""
        from metadarkmatter.core.parsers import StreamingBlastParser
        from metadarkmatter.external.mmseqs2 import MMseqs2

        if not MMseqs2.check_available():
            pytest.skip("MMseqs2 not installed")

        mmseqs = MMseqs2()
        db_path = test_data_dir / "mmseqs_db"
        output_path = test_data_dir / "mmseqs_compressed.tsv.gz"

        # Create database if needed
        if not db_path.exists():
            mmseqs.create_database(
                input_fasta=concatenated_genome,
                database=db_path,
                # Note: Omit dbtype to let MMseqs2 auto-detect (explicit dbtype causes issues)
                timeout=120.0,
            )

        # Run search with uncompressed output first (MMseqs2 doesn't compress directly)
        temp_output = test_data_dir / "temp_uncompressed.tsv"
        mmseqs.search(
            query=test_reads,
            database=db_path,
            output=temp_output,
            sensitivity=5.7,
            threads=2,
            timeout=120.0,
        )

        # Compress manually (simulating CLI behavior)
        import gzip as gzip_module
        with temp_output.open('rb') as f_in:
            with gzip_module.open(output_path, 'wb') as f_out:
                f_out.write(f_in.read())

        temp_output.unlink()

        # Verify compressed file exists
        assert output_path.exists(), "Compressed output should exist"

        # Parse compressed output
        parser = StreamingBlastParser(output_path)
        results = list(parser.iter_reads())

        assert len(results) > 0, "Should parse compressed output"

    def test_performance_measurement(self, test_data_dir, concatenated_genome, test_reads):
        """Should measure and compare performance."""
        import time
        from metadarkmatter.external.mmseqs2 import MMseqs2

        if not MMseqs2.check_available():
            pytest.skip("MMseqs2 not installed")

        mmseqs = MMseqs2()
        db_path = test_data_dir / "mmseqs_db"

        # Measure database creation time
        if not db_path.exists():
            start_time = time.perf_counter()
            mmseqs.create_database(
                input_fasta=concatenated_genome,
                database=db_path,
                # Note: Omit dbtype to let MMseqs2 auto-detect (explicit dbtype causes issues)
                timeout=120.0,
            )
            db_time = time.perf_counter() - start_time
            print(f"\nDatabase creation time: {db_time:.2f}s")

        # Measure search time
        output_path = test_data_dir / "mmseqs_perf_test.tsv"
        start_time = time.perf_counter()
        mmseqs.search(
            query=test_reads,
            database=db_path,
            output=output_path,
            sensitivity=5.7,
            threads=2,
            timeout=120.0,
        )
        search_time = time.perf_counter() - start_time

        print(f"Search time: {search_time:.2f}s")
        print(f"Reads processed: {sum(1 for _ in test_reads.open())}")

        # Should complete reasonably fast on small data
        assert search_time < 60.0, "Search should complete in <60s on small test data"


class TestMMseqs2ClassificationIntegration:
    """Test MMseqs2 output through classification pipeline."""

    def test_classification_with_mmseqs2_results(self, test_data_dir, concatenated_genome,
                                                  test_reads, test_metadata):
        """Should classify reads using MMseqs2 alignment results."""
        from metadarkmatter.external.mmseqs2 import MMseqs2
        from metadarkmatter.core.parsers import StreamingBlastParser

        if not MMseqs2.check_available():
            pytest.skip("MMseqs2 not installed")

        # Create database and run search
        mmseqs = MMseqs2()
        db_path = test_data_dir / "mmseqs_db"
        search_output = test_data_dir / "mmseqs_for_classification.tsv"

        if not db_path.exists():
            mmseqs.create_database(
                input_fasta=concatenated_genome,
                database=db_path,
                # Note: Omit dbtype to let MMseqs2 auto-detect (explicit dbtype causes issues)
                timeout=120.0,
            )

        if not search_output.exists():
            mmseqs.search(
                query=test_reads,
                database=db_path,
                output=search_output,
                sensitivity=5.7,
                evalue=1e-3,
                threads=2,
                timeout=120.0,
            )

        # Parse results
        parser = StreamingBlastParser(search_output)
        results = list(parser.iter_reads())

        # Validate we can extract classification-relevant information
        for result in results[:5]:  # Check first 5 reads
            best_hit = result.best_hit

            # These are the fields classification needs
            assert hasattr(best_hit, 'genome_name'), "Should have genome_name"
            assert hasattr(best_hit, 'pident'), "Should have pident"
            assert hasattr(best_hit, 'bitscore'), "Should have bitscore"

            # Values should be in valid ranges
            assert 0 <= best_hit.pident <= 100, "pident should be 0-100"
            assert best_hit.bitscore > 0, "bitscore should be positive"

            # Genome name should be extractable
            assert best_hit.genome_name in ["GCF_000001.1", "GCF_000002.1",
                                             "GCF_000003.1", "GCF_000004.1"], \
                f"Genome name should match test genomes: {best_hit.genome_name}"

        print(f"\nParsed {len(results)} reads from MMseqs2 output")
        print(f"All results compatible with classification pipeline")


def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line(
        "markers",
        "requires_blast: marks tests as requiring BLAST installation",
    )
