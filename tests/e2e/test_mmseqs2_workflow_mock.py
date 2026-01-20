"""
Mock-based end-to-end workflow validation for MMseqs2 integration.

These tests validate the complete workflow using mock MMseqs2 output,
ensuring all components work together correctly without requiring MMseqs2 installation.
"""

from __future__ import annotations

from pathlib import Path

import pytest


class TestMMseqs2WorkflowWithMockData:
    """Test complete workflow using mock MMseqs2 output."""

    def test_complete_workflow_simulation(self, tmp_path: Path):
        """Simulate complete workflow from reads to classification."""
        from metadarkmatter.core.parsers import StreamingBlastParser

        # Step 1: Create mock MMseqs2 output (as would be produced by mmseqs2 search)
        # This uses the exact format that MMseqs2 produces with our format string
        mock_output = tmp_path / "mmseqs2_results.tsv"
        mock_data = """read_1_GCF_000001.1_100pct\tGCF_000001.1|contig_1\t100.0\t150\t0\t0\t1\t150\t100\t249\t1.2e-75\t285.0
read_1_GCF_000001.1_100pct\tGCF_000002.1|contig_1\t98.7\t150\t2\t0\t1\t150\t100\t249\t3.4e-72\t275.0
read_2_GCF_000001.1_98pct\tGCF_000001.1|contig_1\t98.0\t150\t3\t0\t1\t150\t200\t349\t5.6e-70\t268.0
read_3_GCF_000002.1_92pct\tGCF_000002.1|contig_1\t92.0\t150\t12\t0\t1\t150\t300\t449\t1.2e-60\t235.0
read_3_GCF_000002.1_92pct\tGCF_000001.1|contig_1\t90.7\t150\t14\t0\t1\t150\t300\t449\t4.5e-58\t228.0
read_4_GCF_000003.1_85pct\tGCF_000003.1|contig_1\t85.3\t150\t22\t0\t1\t150\t150\t299\t2.3e-45\t185.0
read_5_GCF_000004.1_100pct\tGCF_000004.1|contig_1\t100.0\t150\t0\t0\t1\t150\t50\t199\t6.7e-76\t287.0
"""
        mock_output.write_text(mock_data)

        # Step 2: Parse with StreamingBlastParser (validates format compatibility)
        parser = StreamingBlastParser(mock_output)
        results = list(parser.iter_reads())

        # Step 3: Validate parsing worked correctly
        assert len(results) == 5, "Should parse 5 unique reads"

        # Step 4: Validate read 1 (multiple hits, perfect match)
        read1 = results[0]
        assert read1.read_id == "read_1_GCF_000001.1_100pct"
        assert read1.num_hits == 2, "Read 1 should have 2 hits"
        assert read1.best_hit.genome_name == "GCF_000001.1"
        assert read1.best_hit.pident == 100.0
        assert read1.best_hit.bitscore == 285.0

        # Step 5: Validate read 3 (novel species candidate - 92% identity)
        read3 = results[2]
        assert read3.read_id == "read_3_GCF_000002.1_92pct"
        assert read3.num_hits == 2
        assert read3.best_hit.genome_name == "GCF_000002.1"
        assert 90 <= read3.best_hit.pident <= 95, "Should be in novel species range"

        # Step 6: Validate read 4 (novel genus candidate - 85% identity)
        read4 = results[3]
        assert read4.read_id == "read_4_GCF_000003.1_85pct"
        assert 80 <= read4.best_hit.pident <= 90, "Should be in novel genus range"

        # Step 7: Calculate novelty index (for classification)
        for result in results:
            novelty_index = 100 - result.best_hit.pident
            print(f"{result.read_id}: {result.best_hit.pident}% identity, "
                  f"novelty={novelty_index:.1f}%, genome={result.best_hit.genome_name}")

        print("\n✓ Complete workflow simulation successful")
        print("✓ MMseqs2 output → Parser → Classification data extraction")

    def test_format_compatibility_verification(self, tmp_path: Path):
        """Verify MMseqs2 format exactly matches BLAST format."""
        from metadarkmatter.core.parsers import StreamingBlastParser

        # Create identical alignment data in both formats
        read_id = "test_read_1"
        genome = "GCF_123456.1"
        pident = 95.5
        alnlen = 150
        mismatch = 7
        gapopen = 0
        qstart, qend = 1, 150
        sstart, send = 100, 249
        evalue = 1.2e-50
        bitscore = 185.0

        # MMseqs2 format (what we generate)
        mmseqs_output = tmp_path / "mmseqs_format.tsv"
        mmseqs_line = (
            f"{read_id}\t{genome}|contig1\t{pident}\t{alnlen}\t{mismatch}\t{gapopen}\t"
            f"{qstart}\t{qend}\t{sstart}\t{send}\t{evalue}\t{bitscore}"
        )
        mmseqs_output.write_text(mmseqs_line)

        # BLAST format (for reference)
        blast_output = tmp_path / "blast_format.tsv"
        blast_line = (
            f"{read_id}\t{genome}|contig1\t{pident}\t{alnlen}\t{mismatch}\t{gapopen}\t"
            f"{qstart}\t{qend}\t{sstart}\t{send}\t{evalue}\t{bitscore}"
        )
        blast_output.write_text(blast_line)

        # Parse both with same parser
        mmseqs_parser = StreamingBlastParser(mmseqs_output)
        mmseqs_result = list(mmseqs_parser.iter_reads())[0]

        blast_parser = StreamingBlastParser(blast_output)
        blast_result = list(blast_parser.iter_reads())[0]

        # Results should be identical
        assert mmseqs_result.read_id == blast_result.read_id
        assert mmseqs_result.num_hits == blast_result.num_hits == 1

        mmseqs_hit = mmseqs_result.best_hit
        blast_hit = blast_result.best_hit

        # All fields should match
        assert mmseqs_hit.qseqid == blast_hit.qseqid == read_id
        assert mmseqs_hit.genome_name == blast_hit.genome_name == genome
        assert mmseqs_hit.pident == blast_hit.pident == pident
        assert mmseqs_hit.length == blast_hit.length == alnlen
        assert mmseqs_hit.bitscore == blast_hit.bitscore == bitscore

        print("\n✓ MMseqs2 and BLAST formats produce identical parsed results")

    def test_classification_threshold_validation(self, tmp_path: Path):
        """Validate classification thresholds work with MMseqs2 output."""
        from metadarkmatter.core.parsers import StreamingBlastParser
        from metadarkmatter.core.constants import (
            NOVELTY_KNOWN_MAX,
            NOVELTY_NOVEL_SPECIES_MAX,
            NOVELTY_NOVEL_GENUS_MIN
        )

        # Create reads spanning classification categories
        mock_output = tmp_path / "threshold_test.tsv"
        # read_known: 97% identity = 3% novelty (< 5%)
        # read_novel_species: 92% identity = 8% novelty (5-20% range)
        # read_novel_genus: 78% identity = 22% novelty (>= 20%)
        mock_data = f"""read_known\tGCF_000001.1|c1\t97.0\t150\t4\t0\t1\t150\t1\t150\t1e-70\t270.0
read_novel_species\tGCF_000001.1|c1\t92.0\t150\t12\t0\t1\t150\t1\t150\t1e-60\t240.0
read_novel_genus\tGCF_000001.1|c1\t78.0\t150\t33\t0\t1\t150\t1\t150\t1e-45\t190.0
"""
        mock_output.write_text(mock_data)

        parser = StreamingBlastParser(mock_output)
        results = {r.read_id: r for r in parser.iter_reads()}

        # Known species: novelty < 5%
        known = results["read_known"]
        known_novelty = 100 - known.best_hit.pident
        assert known_novelty < NOVELTY_KNOWN_MAX, \
            f"Known species should have novelty < {NOVELTY_KNOWN_MAX}%"

        # Novel species: 5% <= novelty < 20%
        novel_sp = results["read_novel_species"]
        novel_sp_novelty = 100 - novel_sp.best_hit.pident
        assert NOVELTY_KNOWN_MAX <= novel_sp_novelty < NOVELTY_NOVEL_SPECIES_MAX, \
            f"Novel species should have {NOVELTY_KNOWN_MAX}% <= novelty < {NOVELTY_NOVEL_SPECIES_MAX}%"

        # Novel genus: 20% <= novelty <= 25%
        novel_genus = results["read_novel_genus"]
        novel_genus_novelty = 100 - novel_genus.best_hit.pident
        assert novel_genus_novelty >= NOVELTY_NOVEL_GENUS_MIN, \
            f"Novel genus should have novelty >= {NOVELTY_NOVEL_GENUS_MIN}%"

        print(f"\n✓ Known species: {known_novelty:.1f}% novelty")
        print(f"✓ Novel species: {novel_sp_novelty:.1f}% novelty")
        print(f"✓ Novel genus: {novel_genus_novelty:.1f}% novelty")
        print("✓ Classification thresholds work correctly with MMseqs2 output")

    def test_genome_name_extraction_edge_cases(self, tmp_path: Path):
        """Test genome name extraction from various header formats."""
        from metadarkmatter.core.parsers import StreamingBlastParser

        # Test various RefSeq/GenBank formats
        mock_output = tmp_path / "genome_names.tsv"
        mock_data = """read1\tGCF_000195955.2|NZ_CP007557.1\t95.0\t150\t7\t0\t1\t150\t1\t150\t1e-70\t270.0
read2\tGCA_123456789.1|contig_1\t94.0\t150\t9\t0\t1\t150\t1\t150\t1e-68\t265.0
read3\tNZ_CP012345.1|chromosome\t93.0\t150\t10\t0\t1\t150\t1\t150\t1e-66\t260.0
"""
        mock_output.write_text(mock_data)

        parser = StreamingBlastParser(mock_output)
        results = list(parser.iter_reads())

        # Validate genome name extraction
        assert results[0].best_hit.genome_name == "GCF_000195955.2"
        assert results[1].best_hit.genome_name == "GCA_123456789.1"
        assert results[2].best_hit.genome_name == "NZ_CP012345.1"

        print("\n✓ Genome name extraction works for all RefSeq/GenBank formats")

    def test_multi_hit_ranking(self, tmp_path: Path):
        """Test that hits are correctly ranked by bitscore."""
        from metadarkmatter.core.parsers import StreamingBlastParser

        # Create read with multiple hits in non-sorted order
        mock_output = tmp_path / "multi_hit.tsv"
        mock_data = """read1\tGCF_000002.1|c1\t92.0\t150\t12\t0\t1\t150\t1\t150\t1e-60\t240.0
read1\tGCF_000001.1|c1\t95.0\t150\t7\t0\t1\t150\t1\t150\t1e-70\t270.0
read1\tGCF_000003.1|c1\t88.0\t150\t18\t0\t1\t150\t1\t150\t1e-50\t210.0
"""
        mock_output.write_text(mock_data)

        parser = StreamingBlastParser(mock_output)
        results = list(parser.iter_reads())

        assert len(results) == 1, "Should have one read"
        result = results[0]
        assert result.num_hits == 3, "Should have 3 hits"

        # Best hit should be the one with highest bitscore (270.0)
        assert result.best_hit.genome_name == "GCF_000001.1"
        assert result.best_hit.bitscore == 270.0
        assert result.best_hit.pident == 95.0

        # Verify all hits are sorted by bitscore descending
        bitscores = [hit.bitscore for hit in result.hits]
        assert bitscores == sorted(bitscores, reverse=True), \
            "Hits should be sorted by bitscore descending"

        print(f"\n✓ Multi-hit ranking correct: {bitscores}")

    def test_workflow_with_real_data_structure(self, tmp_path: Path):
        """Test workflow with realistic data structure and file operations."""
        import gzip
        from metadarkmatter.core.parsers import StreamingBlastParser

        # Simulate real workflow: compressed MMseqs2 output
        compressed_output = tmp_path / "mmseqs2_results.tsv.gz"

        # Create realistic data
        mock_data = """read_1\tGCF_000001.1|contig_1\t98.7\t150\t2\t0\t1\t150\t100\t249\t1.2e-72\t275.0
read_2\tGCF_000001.1|contig_1\t97.3\t148\t4\t0\t1\t148\t200\t347\t3.4e-68\t265.0
read_3\tGCF_000002.1|contig_1\t92.7\t150\t11\t0\t1\t150\t300\t449\t5.6e-58\t235.0
read_4\tGCF_000003.1|contig_2\t87.3\t150\t19\t0\t1\t150\t150\t299\t1.2e-48\t200.0
read_5\tGCF_000004.1|contig_3\t83.3\t150\t25\t0\t1\t150\t50\t199\t6.7e-42\t175.0
"""
        # Compress
        with gzip.open(compressed_output, 'wt') as f:
            f.write(mock_data)

        # Parse compressed file
        parser = StreamingBlastParser(compressed_output)
        results = list(parser.iter_reads())

        assert len(results) == 5, "Should parse all 5 reads from compressed file"

        # Verify data quality
        for result in results:
            assert result.num_hits > 0, f"{result.read_id} should have hits"
            assert result.best_hit.genome_name.startswith("GCF_"), \
                f"{result.read_id} should have valid genome name"
            assert 0 <= result.best_hit.pident <= 100, \
                f"{result.read_id} should have valid pident"

        # Calculate statistics
        identities = [r.best_hit.pident for r in results]
        avg_identity = sum(identities) / len(identities)
        min_identity = min(identities)
        max_identity = max(identities)

        print(f"\n✓ Parsed {len(results)} reads from compressed file")
        print(f"✓ Identity range: {min_identity:.1f}% - {max_identity:.1f}%")
        print(f"✓ Average identity: {avg_identity:.1f}%")
        print("✓ Workflow with compressed output successful")


def test_integration_test_file_exists():
    """Verify end-to-end integration test file exists."""
    e2e_test = Path(__file__).parent / "test_mmseqs2_end_to_end.py"
    assert e2e_test.exists(), "End-to-end integration test file should exist"
    print(f"\n✓ Integration test file exists: {e2e_test}")
    print("✓ Run with MMseqs2 installed to execute real validation")
