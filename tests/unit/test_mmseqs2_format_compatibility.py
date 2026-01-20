"""
Unit tests for MMseqs2 output format compatibility with BLAST parser.

These tests validate format compatibility without requiring MMseqs2 installation.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metadarkmatter.external.mmseqs2 import MMseqs2


class TestMMseqs2FormatCompatibility:
    """Tests for BLAST parser format compatibility."""

    def test_column_count_matches_blast(self):
        """Should produce 12 columns to match BLAST format."""
        mmseqs = MMseqs2()

        # BLAST format has 12 columns
        expected_columns = 12

        # Count columns in MMseqs2 format string
        format_columns = mmseqs.BLAST_FORMAT_COLUMNS
        assert len(format_columns) == expected_columns, (
            f"MMseqs2 format should have {expected_columns} columns, "
            f"got {len(format_columns)}"
        )

    def test_format_string_structure(self):
        """Should have correct format string for MMseqs2."""
        mmseqs = MMseqs2()

        format_string = ",".join(mmseqs.BLAST_FORMAT_COLUMNS)

        # Should be comma-separated
        assert "," in format_string, "Format string should be comma-separated"

        # Should contain key fields
        required_fields = ["query", "target", "pident", "evalue", "bits"]
        for field in required_fields:
            assert field in format_string, f"Format should include {field}"

    def test_column_order_matches_blast_parser(self):
        """Should match column order expected by StreamingBlastParser."""
        from metadarkmatter.core.parsers import StreamingBlastParser

        mmseqs = MMseqs2()

        # BLAST parser expects these columns in this order
        blast_columns = StreamingBlastParser.COLUMN_NAMES

        # MMseqs2 format columns (positional mapping)
        mmseqs_columns = mmseqs.BLAST_FORMAT_COLUMNS

        # Create mapping of what MMseqs2 columns map to
        column_mapping = {
            "query": "qseqid",
            "target": "sseqid",
            "pident": "pident",
            "alnlen": "length",
            "mismatch": "mismatch",
            "gapopen": "gapopen",
            "qstart": "qstart",
            "qend": "qend",
            "tstart": "sstart",
            "tend": "send",
            "evalue": "evalue",
            "bits": "bitscore",
        }

        # Verify each position maps correctly
        for idx, (mmseqs_col, blast_col) in enumerate(zip(mmseqs_columns, blast_columns)):
            expected_blast_col = column_mapping[mmseqs_col]
            assert blast_col == expected_blast_col, (
                f"Position {idx}: MMseqs2 column '{mmseqs_col}' maps to "
                f"BLAST column '{expected_blast_col}', but parser expects '{blast_col}'"
            )

    def test_mock_output_parses_correctly(self, tmp_path: Path):
        """Should parse mock MMseqs2 output with BLAST parser."""
        from metadarkmatter.core.parsers import StreamingBlastParser

        # Create mock MMseqs2 output (12 columns, no header)
        # Using proper RefSeq format with version numbers
        mock_output = tmp_path / "mock_mmseqs2.tsv"
        mock_data = """read1\tGCF_000001.1|contig1\t95.5\t100\t4\t0\t1\t100\t1\t100\t1.2e-50\t185.0
read2\tGCF_000002.1|contig1\t88.3\t120\t14\t0\t1\t120\t1\t120\t3.4e-40\t152.0
read2\tGCF_000001.1|contig1\t85.0\t115\t17\t1\t1\t115\t1\t115\t5.6e-35\t140.0
read3\tGCF_000003.1|contig1\t100.0\t50\t0\t0\t1\t50\t1\t50\t1.0e-25\t95.0
"""
        mock_output.write_text(mock_data)

        # Parse with StreamingBlastParser
        parser = StreamingBlastParser(mock_output)
        results = list(parser.iter_reads())

        # Should parse 3 reads
        assert len(results) == 3, "Should parse all reads"

        # Validate first read
        read1 = results[0]
        assert read1.read_id == "read1"
        assert read1.num_hits == 1
        assert read1.best_hit.genome_name == "GCF_000001.1"
        assert read1.best_hit.pident == 95.5
        assert read1.best_hit.bitscore == 185.0

        # Validate second read (multiple hits)
        read2 = results[1]
        assert read2.read_id == "read2"
        assert read2.num_hits == 2
        assert read2.best_hit.genome_name == "GCF_000002.1"  # Higher bitscore
        assert read2.best_hit.bitscore == 152.0

    def test_genome_name_extraction(self, tmp_path: Path):
        """Should extract genome names from standardized headers."""
        from metadarkmatter.core.parsers import StreamingBlastParser

        # Create mock output with various header formats
        mock_output = tmp_path / "mock_headers.tsv"
        mock_data = """read1\tGCF_000195955.2|NZ_CP007557.1\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t180.0
read2\tGCA_123456789.1|contig_1\t90.0\t100\t10\t0\t1\t100\t1\t100\t1e-40\t160.0
"""
        mock_output.write_text(mock_data)

        parser = StreamingBlastParser(mock_output)
        results = list(parser.iter_reads())

        # Should extract genome accessions correctly
        assert results[0].best_hit.genome_name == "GCF_000195955.2"
        assert results[1].best_hit.genome_name == "GCA_123456789.1"

    def test_numeric_field_ranges(self, tmp_path: Path):
        """Should validate numeric field ranges."""
        from metadarkmatter.core.parsers import StreamingBlastParser

        # Create output with edge case values
        mock_output = tmp_path / "mock_ranges.tsv"
        mock_data = """read1\tGCF_000001.1|c1\t0.0\t10\t10\t0\t1\t10\t1\t10\t1.0\t20.0
read2\tGCF_000002.1|c1\t100.0\t100\t0\t0\t1\t100\t1\t100\t1e-100\t500.0
read3\tGCF_000003.1|c1\t50.0\t50\t25\t0\t1\t50\t1\t50\t0.001\t100.0
"""
        mock_output.write_text(mock_data)

        parser = StreamingBlastParser(mock_output)
        results = list(parser.iter_reads())

        # Validate ranges
        for result in results:
            hit = result.best_hit
            assert 0 <= hit.pident <= 100, "pident should be 0-100"
            assert hit.bitscore > 0, "bitscore should be positive"
            assert hit.evalue >= 0, "evalue should be non-negative"

    def test_format_mode_parameter(self, tmp_path: Path):
        """Should use format-mode 0 for simple tabular output."""
        mmseqs = MMseqs2()
        MMseqs2._executable_cache["mmseqs"] = Path("/usr/bin/mmseqs")

        query = tmp_path / "query.fasta"
        query.write_text(">read1\nACGT\n")

        cmd = mmseqs.build_command(
            mode="search",
            query=query,
            database=tmp_path / "db",
            output=tmp_path / "out.tsv",
        )

        # Should include format-mode 0
        assert "--format-mode" in cmd
        format_mode_idx = cmd.index("--format-mode")
        assert cmd[format_mode_idx + 1] == "0", "Should use format-mode 0"

    def test_format_output_parameter(self, tmp_path: Path):
        """Should specify exact column format."""
        mmseqs = MMseqs2()
        MMseqs2._executable_cache["mmseqs"] = Path("/usr/bin/mmseqs")

        query = tmp_path / "query.fasta"
        query.write_text(">read1\nACGT\n")

        cmd = mmseqs.build_command(
            mode="search",
            query=query,
            database=tmp_path / "db",
            output=tmp_path / "out.tsv",
        )

        # Should include format-output
        assert "--format-output" in cmd
        format_output_idx = cmd.index("--format-output")

        # Next argument should be the format string
        format_str = cmd[format_output_idx + 1]

        # Should match BLAST_FORMAT_COLUMNS
        expected_format = ",".join(mmseqs.BLAST_FORMAT_COLUMNS)
        assert format_str == expected_format, (
            f"Format string mismatch: expected '{expected_format}', got '{format_str}'"
        )

    def test_documentation_examples(self):
        """Should validate format strings match documentation."""
        mmseqs = MMseqs2()

        # Format string from docstring/documentation
        documented_format = (
            "query,target,pident,alnlen,mismatch,gapopen,"
            "qstart,qend,tstart,tend,evalue,bits"
        )

        actual_format = ",".join(mmseqs.BLAST_FORMAT_COLUMNS)

        assert actual_format == documented_format, (
            "Format string should match documentation"
        )


class TestMMseqs2BlastEquivalence:
    """Tests for parameter equivalence between MMseqs2 and BLAST."""

    def test_evalue_parameter_equivalence(self, tmp_path: Path):
        """Should use same E-value scale as BLAST."""
        from metadarkmatter.external.blast import BlastN

        mmseqs = MMseqs2()
        blastn = BlastN()

        # Cache executables
        MMseqs2._executable_cache["mmseqs"] = Path("/usr/bin/mmseqs")
        BlastN._executable_cache["blastn"] = Path("/usr/bin/blastn")

        query = tmp_path / "query.fasta"
        query.write_text(">read1\nACGT\n")

        # Build commands with same E-value
        mmseqs_cmd = mmseqs.build_command(
            mode="search",
            query=query,
            database=tmp_path / "db",
            output=tmp_path / "out.tsv",
            evalue=1e-3,
        )

        blast_cmd = blastn.build_command(
            query=query,
            database=tmp_path / "db",
            output=tmp_path / "out.tsv",
            evalue=1e-3,
        )

        # Both should have E-value parameter
        assert "-e" in mmseqs_cmd
        assert "-evalue" in blast_cmd

        # MMseqs2 E-value should be in scientific notation or decimal
        mmseqs_evalue_idx = mmseqs_cmd.index("-e")
        mmseqs_evalue = mmseqs_cmd[mmseqs_evalue_idx + 1]
        assert "0.001" in mmseqs_evalue or "1e-03" in mmseqs_evalue

    def test_max_targets_parameter_equivalence(self, tmp_path: Path):
        """Should use equivalent max targets parameter."""
        from metadarkmatter.external.blast import BlastN

        mmseqs = MMseqs2()
        blastn = BlastN()

        MMseqs2._executable_cache["mmseqs"] = Path("/usr/bin/mmseqs")
        BlastN._executable_cache["blastn"] = Path("/usr/bin/blastn")

        query = tmp_path / "query.fasta"
        query.write_text(">read1\nACGT\n")

        # Build commands with same max targets
        mmseqs_cmd = mmseqs.build_command(
            mode="search",
            query=query,
            database=tmp_path / "db",
            output=tmp_path / "out.tsv",
            max_seqs=500,
        )

        blast_cmd = blastn.build_command(
            query=query,
            database=tmp_path / "db",
            output=tmp_path / "out.tsv",
            max_target_seqs=500,
        )

        # MMseqs2 uses --max-seqs
        assert "--max-seqs" in mmseqs_cmd
        max_seqs_idx = mmseqs_cmd.index("--max-seqs")
        assert mmseqs_cmd[max_seqs_idx + 1] == "500"

        # BLAST uses -max_target_seqs
        assert "-max_target_seqs" in blast_cmd
