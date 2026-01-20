"""
Kraken2 and KrakenTools wrapper classes.

Provides Python interfaces for:
- Kraken2: Taxonomic classification of metagenomic sequences
- ExtractKrakenReads: Extraction of reads assigned to specific taxa
"""

from __future__ import annotations

from pathlib import Path

from metadarkmatter.external.base import ExternalTool


class Kraken2(ExternalTool):
    """Wrapper for Kraken2 taxonomic classifier.

    Kraken2 assigns taxonomic labels to short DNA sequences using
    exact k-mer matches to a database of reference genomes.

    Example:
        >>> kraken = Kraken2()
        >>> result = kraken.run(
        ...     reads_1=Path("sample_R1.fastq.gz"),
        ...     reads_2=Path("sample_R2.fastq.gz"),
        ...     database=Path("/path/to/kraken_db"),
        ...     output=Path("sample.kraken"),
        ...     report=Path("sample.kreport"),
        ... )
    """

    TOOL_NAME = "kraken2"
    INSTALL_HINT = "conda install -c bioconda kraken2"

    def build_command(
        self,
        *,
        reads_1: Path,
        database: Path,
        output: Path,
        report: Path,
        reads_2: Path | None = None,
        threads: int = 4,
        confidence: float = 0.0,
        use_names: bool = True,
        memory_mapping: bool = False,
        report_zero_counts: bool = True,
        gzip_compressed: bool | None = None,
        bzip2_compressed: bool = False,
    ) -> list[str]:
        """Build Kraken2 command.

        Args:
            reads_1: Forward reads file (FASTQ/FASTA, optionally gzipped).
            database: Path to Kraken2 database directory.
            output: Output file for per-read classification.
            report: Output file for classification report.
            reads_2: Reverse reads file for paired-end data.
            threads: Number of threads to use.
            confidence: Confidence threshold (0-1) for classification.
            use_names: Include taxon names in output.
            memory_mapping: Use memory mapping (slower but lower RAM).
            report_zero_counts: Include taxa with zero counts in report.
            gzip_compressed: Force gzip decompression (auto-detected if None).
            bzip2_compressed: Input is bzip2 compressed.

        Returns:
            Command as list of strings.
        """
        exe = str(self.get_executable())
        cmd = [exe]

        # Database
        cmd.extend(["--db", str(database)])

        # Output files
        cmd.extend(["--output", str(output)])
        cmd.extend(["--report", str(report)])

        # Processing options
        cmd.extend(["--threads", str(threads)])

        if confidence > 0:
            cmd.extend(["--confidence", str(confidence)])

        if use_names:
            cmd.append("--use-names")

        if memory_mapping:
            cmd.append("--memory-mapping")

        if report_zero_counts:
            cmd.append("--report-zero-counts")

        # Compression handling
        if gzip_compressed is True:
            cmd.append("--gzip-compressed")
        if bzip2_compressed:
            cmd.append("--bzip2-compressed")

        # Input files
        if reads_2 is not None:
            cmd.append("--paired")
            cmd.extend([str(reads_1), str(reads_2)])
        else:
            cmd.append(str(reads_1))

        return cmd


class ExtractKrakenReads(ExternalTool):
    """Wrapper for KrakenTools extract_kraken_reads.py.

    Extracts reads assigned to specific taxonomic IDs from Kraken2 output.
    Can include all descendants of a taxon (--include-children).

    Example:
        >>> extractor = ExtractKrakenReads()
        >>> result = extractor.run(
        ...     kraken_output=Path("sample.kraken"),
        ...     kraken_report=Path("sample.kreport"),
        ...     reads_1=Path("sample_R1.fastq.gz"),
        ...     reads_2=Path("sample_R2.fastq.gz"),
        ...     taxid=1224,  # Proteobacteria
        ...     output_1=Path("proteobacteria_R1.fastq"),
        ...     output_2=Path("proteobacteria_R2.fastq"),
        ...     include_children=True,
        ... )
    """

    TOOL_NAME = "extract_kraken_reads.py"
    TOOL_ALIASES = ("extractKrakenReads.py", "extract-kraken-reads")
    INSTALL_HINT = "pip install krakentools  # or: conda install -c bioconda krakentools"

    def build_command(
        self,
        *,
        kraken_output: Path,
        kraken_report: Path,
        reads_1: Path,
        output_1: Path,
        taxid: int | list[int],
        reads_2: Path | None = None,
        output_2: Path | None = None,
        include_children: bool = True,
        include_parents: bool = False,
        exclude: bool = False,
        fastq_output: bool = True,
    ) -> list[str]:
        """Build extract_kraken_reads.py command.

        Args:
            kraken_output: Kraken2 output file (.kraken).
            kraken_report: Kraken2 report file (.kreport).
            reads_1: Forward reads file used in Kraken2 run.
            output_1: Output file for extracted forward reads.
            taxid: Taxonomic ID(s) to extract.
            reads_2: Reverse reads file (for paired-end).
            output_2: Output file for extracted reverse reads.
            include_children: Include reads from child taxa.
            include_parents: Include reads from parent taxa.
            exclude: Exclude (instead of include) the specified taxa.
            fastq_output: Output in FASTQ format (vs FASTA).

        Returns:
            Command as list of strings.
        """
        exe = str(self.get_executable())
        cmd = [exe]

        # Input files
        cmd.extend(["-k", str(kraken_output)])
        cmd.extend(["-r", str(kraken_report)])
        cmd.extend(["-s1", str(reads_1)])

        if reads_2 is not None:
            cmd.extend(["-s2", str(reads_2)])

        # Taxonomic ID(s)
        if isinstance(taxid, int):
            cmd.extend(["-t", str(taxid)])
        else:
            cmd.extend(["-t", " ".join(str(t) for t in taxid)])

        # Options
        if include_children:
            cmd.append("--include-children")

        if include_parents:
            cmd.append("--include-parents")

        if exclude:
            cmd.append("--exclude")

        if fastq_output:
            cmd.append("--fastq-output")

        # Output files
        cmd.extend(["-o", str(output_1)])

        if output_2 is not None:
            cmd.extend(["-o2", str(output_2)])

        return cmd


class KrakenReport:
    """Parser for Kraken2 report files.

    Parses .kreport files to extract taxonomic counts and lineage information.
    """

    def __init__(self, report_path: Path):
        """Initialize with path to kreport file.

        Args:
            report_path: Path to Kraken2 report file.
        """
        self.report_path = report_path
        self._data: list[dict[str, object]] | None = None

    def parse(self) -> list[dict[str, object]]:
        """Parse the report file.

        Returns:
            List of dicts with keys: pct, clade_reads, taxon_reads,
            rank, taxid, name.
        """
        if self._data is not None:
            return self._data

        self._data = []

        with self.report_path.open() as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 6:
                    self._data.append({
                        "pct": float(parts[0]),
                        "clade_reads": int(parts[1]),
                        "taxon_reads": int(parts[2]),
                        "rank": parts[3],
                        "taxid": int(parts[4]),
                        "name": parts[5].strip(),
                    })

        return self._data

    def get_taxon_reads(self, taxid: int) -> int:
        """Get number of reads assigned directly to a taxon.

        Args:
            taxid: Taxonomic ID to look up.

        Returns:
            Number of reads assigned to this taxon (not including children).
        """
        for entry in self.parse():
            if entry["taxid"] == taxid:
                return int(entry["taxon_reads"])
        return 0

    def get_clade_reads(self, taxid: int) -> int:
        """Get number of reads assigned to a taxon and its descendants.

        Args:
            taxid: Taxonomic ID to look up.

        Returns:
            Number of reads in this clade (including all children).
        """
        for entry in self.parse():
            if entry["taxid"] == taxid:
                return int(entry["clade_reads"])
        return 0

    def get_family_taxids(self) -> list[int]:
        """Get all family-level taxonomic IDs.

        Returns:
            List of taxids with rank 'F' (family).
        """
        return [
            int(entry["taxid"])
            for entry in self.parse()
            if entry["rank"] == "F"
        ]
