"""
Test data factories for E2E testing.

Provides deterministic, seeded mock data generation for BLAST files
and ANI matrices to enable reproducible end-to-end testing.
"""

from __future__ import annotations

import random
from pathlib import Path
from typing import NamedTuple

import polars as pl


class BlastRecord(NamedTuple):
    """A single BLAST alignment record."""

    qseqid: str
    sseqid: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float


class BlastDataFactory:
    """
    Factory for generating synthetic BLAST data with controlled characteristics.

    Creates reads with specific classification outcomes for testing:
    - Known species (high identity, unambiguous)
    - Novel species (moderate identity divergence)
    - Novel genus (high identity divergence)
    - Conserved regions (ambiguous multi-genome hits)

    All data generation is seeded for reproducibility.
    """

    # Standard RefSeq-style genome prefixes for realistic test data
    GENOME_PREFIXES = [
        "GCF_000123456.1",
        "GCF_000789012.1",
        "GCA_000111222.1",
        "GCF_000333444.1",
        "GCA_000555666.1",
        "GCF_000777888.1",
        "GCA_000999000.1",
        "GCF_001111222.1",
    ]

    def __init__(self, seed: int = 42):
        """Initialize with reproducible seed."""
        self._rng = random.Random(seed)
        self._read_counter = 0

    def _next_read_id(self, prefix: str = "read") -> str:
        """Generate unique read ID."""
        self._read_counter += 1
        return f"{prefix}_{self._read_counter:06d}"

    def _genome_sseqid(self, genome_id: str) -> str:
        """Convert genome ID to full sseqid format.

        Uses the standardized pipe-delimited format: {accession}|{contig_id}
        This matches the format produced by blast makedb header rewriting.
        """
        return f"{genome_id}|ASM1v1_genomic"

    def create_known_species_read(
        self,
        read_id: str | None = None,
        genome: str | None = None,
        pident: float = 99.0,
        bitscore: float = 280.0,
    ) -> list[BlastRecord]:
        """
        Create a read with a single high-identity hit (known species).

        Args:
            read_id: Read identifier (auto-generated if None)
            genome: Target genome ID (random selection if None)
            pident: Percent identity (default 99.0 for known species)
            bitscore: Alignment bitscore

        Returns:
            List with single BlastRecord for this read
        """
        if read_id is None:
            read_id = self._next_read_id("known")
        if genome is None:
            genome = self._rng.choice(self.GENOME_PREFIXES)

        start_pos = self._rng.randint(1000, 100000)
        length = self._rng.randint(100, 300)

        return [
            BlastRecord(
                qseqid=read_id,
                sseqid=self._genome_sseqid(genome),
                pident=pident,
                length=length,
                mismatch=int(length * (100 - pident) / 100),
                gapopen=0,
                qstart=1,
                qend=length,
                sstart=start_pos,
                send=start_pos + length - 1,
                evalue=1e-80,
                bitscore=bitscore,
            )
        ]

    def create_novel_species_read(
        self,
        read_id: str | None = None,
        genome: str | None = None,
        pident: float = 90.0,
        bitscore: float = 220.0,
    ) -> list[BlastRecord]:
        """
        Create a read with moderate identity (novel species, 5-20% divergence).

        Args:
            read_id: Read identifier (auto-generated if None)
            genome: Target genome ID (random selection if None)
            pident: Percent identity (default 90.0 for novel species)
            bitscore: Alignment bitscore

        Returns:
            List with single BlastRecord for this read
        """
        if read_id is None:
            read_id = self._next_read_id("novel_sp")
        if genome is None:
            genome = self._rng.choice(self.GENOME_PREFIXES)

        start_pos = self._rng.randint(1000, 100000)
        length = self._rng.randint(100, 300)

        return [
            BlastRecord(
                qseqid=read_id,
                sseqid=self._genome_sseqid(genome),
                pident=pident,
                length=length,
                mismatch=int(length * (100 - pident) / 100),
                gapopen=self._rng.randint(0, 2),
                qstart=1,
                qend=length,
                sstart=start_pos,
                send=start_pos + length - 1,
                evalue=1e-50,
                bitscore=bitscore,
            )
        ]

    def create_novel_genus_read(
        self,
        read_id: str | None = None,
        genome: str | None = None,
        pident: float = 82.0,
        bitscore: float = 180.0,
    ) -> list[BlastRecord]:
        """
        Create a read with low identity (novel genus, 20-25% divergence).

        Args:
            read_id: Read identifier (auto-generated if None)
            genome: Target genome ID (random selection if None)
            pident: Percent identity (default 82.0 for novel genus)
            bitscore: Alignment bitscore

        Returns:
            List with single BlastRecord for this read
        """
        if read_id is None:
            read_id = self._next_read_id("novel_gen")
        if genome is None:
            genome = self._rng.choice(self.GENOME_PREFIXES)

        start_pos = self._rng.randint(1000, 100000)
        length = self._rng.randint(100, 300)

        return [
            BlastRecord(
                qseqid=read_id,
                sseqid=self._genome_sseqid(genome),
                pident=pident,
                length=length,
                mismatch=int(length * (100 - pident) / 100),
                gapopen=self._rng.randint(1, 4),
                qstart=1,
                qend=length,
                sstart=start_pos,
                send=start_pos + length - 1,
                evalue=1e-30,
                bitscore=bitscore,
            )
        ]

    def create_conserved_region_read(
        self,
        read_id: str | None = None,
        genome1: str | None = None,
        genome2: str | None = None,
        pident1: float = 99.0,
        pident2: float = 98.5,
        bitscore1: float = 280.0,
        bitscore2: float = 275.0,
    ) -> list[BlastRecord]:
        """
        Create a read with multiple high-scoring hits (conserved region).

        These represent reads mapping to conserved regions across multiple
        genomes with similar scores, creating placement ambiguity.

        Args:
            read_id: Read identifier (auto-generated if None)
            genome1: First target genome ID
            genome2: Second target genome ID
            pident1: Percent identity for first hit
            pident2: Percent identity for second hit
            bitscore1: Bitscore for first hit
            bitscore2: Bitscore for second hit

        Returns:
            List with two BlastRecords for this read
        """
        if read_id is None:
            read_id = self._next_read_id("conserved")

        available = self.GENOME_PREFIXES.copy()
        if genome1 is None:
            genome1 = self._rng.choice(available)
        available.remove(genome1)
        if genome2 is None:
            genome2 = self._rng.choice(available)

        length = self._rng.randint(100, 300)
        start1 = self._rng.randint(1000, 100000)
        start2 = self._rng.randint(1000, 100000)

        return [
            BlastRecord(
                qseqid=read_id,
                sseqid=self._genome_sseqid(genome1),
                pident=pident1,
                length=length,
                mismatch=int(length * (100 - pident1) / 100),
                gapopen=0,
                qstart=1,
                qend=length,
                sstart=start1,
                send=start1 + length - 1,
                evalue=1e-80,
                bitscore=bitscore1,
            ),
            BlastRecord(
                qseqid=read_id,
                sseqid=self._genome_sseqid(genome2),
                pident=pident2,
                length=length,
                mismatch=int(length * (100 - pident2) / 100),
                gapopen=0,
                qstart=1,
                qend=length,
                sstart=start2,
                send=start2 + length - 1,
                evalue=1e-78,
                bitscore=bitscore2,
            ),
        ]

    def create_mixed_sample(
        self,
        n_known: int = 10,
        n_novel_species: int = 5,
        n_novel_genus: int = 3,
        n_conserved: int = 5,
    ) -> list[BlastRecord]:
        """
        Create a mixed sample with various read types.

        Args:
            n_known: Number of known species reads
            n_novel_species: Number of novel species reads
            n_novel_genus: Number of novel genus reads
            n_conserved: Number of conserved region reads

        Returns:
            List of all BlastRecords, sorted by read ID then bitscore
        """
        records: list[BlastRecord] = []

        for _ in range(n_known):
            records.extend(self.create_known_species_read())

        for _ in range(n_novel_species):
            records.extend(self.create_novel_species_read())

        for _ in range(n_novel_genus):
            records.extend(self.create_novel_genus_read())

        for _ in range(n_conserved):
            records.extend(self.create_conserved_region_read())

        # Sort by qseqid, then by bitscore descending
        records.sort(key=lambda r: (r.qseqid, -r.bitscore))
        return records

    def write_blast_file(self, records: list[BlastRecord], path: Path) -> Path:
        """
        Write BLAST records to a TSV file.

        Args:
            records: List of BlastRecords to write
            path: Output file path

        Returns:
            Path to the written file
        """
        lines = []
        for r in records:
            line = (
                f"{r.qseqid}\t{r.sseqid}\t{r.pident}\t{r.length}\t"
                f"{r.mismatch}\t{r.gapopen}\t{r.qstart}\t{r.qend}\t"
                f"{r.sstart}\t{r.send}\t{r.evalue}\t{r.bitscore}"
            )
            lines.append(line)

        path.write_text("\n".join(lines) + "\n")
        return path


class ANIMatrixFactory:
    """
    Factory for generating ANI matrices with controlled characteristics.

    Creates matrices representing different evolutionary scenarios:
    - High ANI (same species, >95% ANI)
    - Moderate ANI (same genus, 80-95% ANI)
    - Low ANI (different genera, <80% ANI)
    """

    def __init__(self, seed: int = 42):
        """Initialize with reproducible seed."""
        self._rng = random.Random(seed)

    def create_matrix(
        self,
        genomes: list[str],
        species_groups: dict[str, list[str]] | None = None,
        default_ani: float = 75.0,
        intra_species_ani: float = 98.0,
        intra_genus_ani: float = 88.0,
    ) -> dict[str, dict[str, float]]:
        """
        Create an ANI matrix with specified genome relationships.

        Args:
            genomes: List of genome identifiers
            species_groups: Dict mapping species name to genome IDs in that species
            default_ani: ANI value for unrelated genomes
            intra_species_ani: ANI value for same-species pairs
            intra_genus_ani: ANI value for same-genus pairs (unused for now)

        Returns:
            Nested dict {genome1: {genome2: ani_value}}
        """
        matrix: dict[str, dict[str, float]] = {}

        for g1 in genomes:
            matrix[g1] = {}
            for g2 in genomes:
                if g1 == g2:
                    matrix[g1][g2] = 100.0
                else:
                    matrix[g1][g2] = default_ani

        # Apply species groupings
        if species_groups:
            for _species, members in species_groups.items():
                for g1 in members:
                    for g2 in members:
                        if g1 != g2 and g1 in matrix and g2 in matrix[g1]:
                            # Add small variation for realism
                            variation = self._rng.uniform(-0.5, 0.5)
                            matrix[g1][g2] = intra_species_ani + variation

        return matrix

    def create_high_ani_matrix(
        self,
        genomes: list[str],
        min_ani: float = 92.0,
        max_ani: float = 98.0,
    ) -> dict[str, dict[str, float]]:
        """
        Create a matrix where all genomes are closely related (high ANI).

        Args:
            genomes: List of genome identifiers
            min_ani: Minimum ANI between any two genomes
            max_ani: Maximum ANI between any two genomes

        Returns:
            Nested dict {genome1: {genome2: ani_value}}
        """
        matrix: dict[str, dict[str, float]] = {}

        for g1 in genomes:
            matrix[g1] = {}
            for g2 in genomes:
                if g1 == g2:
                    matrix[g1][g2] = 100.0
                else:
                    matrix[g1][g2] = self._rng.uniform(min_ani, max_ani)

        return matrix

    def create_low_ani_matrix(
        self,
        genomes: list[str],
        min_ani: float = 70.0,
        max_ani: float = 78.0,
    ) -> dict[str, dict[str, float]]:
        """
        Create a matrix where genomes are distantly related (low ANI).

        Args:
            genomes: List of genome identifiers
            min_ani: Minimum ANI between any two genomes
            max_ani: Maximum ANI between any two genomes

        Returns:
            Nested dict {genome1: {genome2: ani_value}}
        """
        matrix: dict[str, dict[str, float]] = {}

        for g1 in genomes:
            matrix[g1] = {}
            for g2 in genomes:
                if g1 == g2:
                    matrix[g1][g2] = 100.0
                else:
                    matrix[g1][g2] = self._rng.uniform(min_ani, max_ani)

        return matrix

    def write_ani_file(
        self,
        matrix: dict[str, dict[str, float]],
        path: Path,
    ) -> Path:
        """
        Write ANI matrix to a CSV file.

        Args:
            matrix: ANI matrix as nested dict
            path: Output file path

        Returns:
            Path to the written file
        """
        genomes = sorted(matrix.keys())

        data: dict[str, list[float | str]] = {"genome": genomes}
        for g in genomes:
            data[g] = [matrix[genome][g] for genome in genomes]

        df = pl.DataFrame(data)
        df.write_csv(path)
        return path


class E2ETestDataset:
    """
    Generates complete test datasets for E2E testing.

    Combines BlastDataFactory and ANIMatrixFactory to create
    matched BLAST files and ANI matrices for testing the full
    classification pipeline.
    """

    def __init__(self, temp_dir: Path, seed: int = 42):
        """
        Initialize dataset generator.

        Args:
            temp_dir: Directory for storing generated files
            seed: Random seed for reproducibility
        """
        self.temp_dir = temp_dir
        self.seed = seed
        self.blast_factory = BlastDataFactory(seed=seed)
        self.ani_factory = ANIMatrixFactory(seed=seed)

    def create_standard_dataset(
        self,
        n_known: int = 20,
        n_novel_species: int = 10,
        n_novel_genus: int = 5,
        n_conserved: int = 8,
    ) -> tuple[Path, Path]:
        """
        Create a standard mixed dataset for testing.

        Args:
            n_known: Number of known species reads
            n_novel_species: Number of novel species reads
            n_novel_genus: Number of novel genus reads
            n_conserved: Number of conserved region reads

        Returns:
            Tuple of (blast_path, ani_path)
        """
        # Create BLAST data
        records = self.blast_factory.create_mixed_sample(
            n_known=n_known,
            n_novel_species=n_novel_species,
            n_novel_genus=n_novel_genus,
            n_conserved=n_conserved,
        )

        blast_path = self.temp_dir / "standard_sample.blast.tsv"
        self.blast_factory.write_blast_file(records, blast_path)

        # Create matching ANI matrix
        genomes = BlastDataFactory.GENOME_PREFIXES[:8]
        species_groups = {
            "species_A": [genomes[0], genomes[1]],
            "species_B": [genomes[2], genomes[3]],
        }
        matrix = self.ani_factory.create_matrix(
            genomes=genomes,
            species_groups=species_groups,
        )

        ani_path = self.temp_dir / "standard_sample.ani.csv"
        self.ani_factory.write_ani_file(matrix, ani_path)

        return blast_path, ani_path

    def create_batch_dataset(
        self,
        n_samples: int = 3,
        reads_per_sample: int = 20,
    ) -> tuple[Path, Path]:
        """
        Create multiple BLAST files for batch processing tests.

        Args:
            n_samples: Number of sample files to create
            reads_per_sample: Approximate reads per sample

        Returns:
            Tuple of (blast_dir, ani_path)
        """
        blast_dir = self.temp_dir / "batch_blast"
        blast_dir.mkdir(exist_ok=True)

        for i in range(n_samples):
            # Vary composition slightly per sample
            n_known = reads_per_sample // 2
            n_novel = reads_per_sample // 4
            n_genus = reads_per_sample // 8
            n_conserved = reads_per_sample - n_known - n_novel - n_genus

            records = self.blast_factory.create_mixed_sample(
                n_known=n_known,
                n_novel_species=n_novel,
                n_novel_genus=n_genus,
                n_conserved=n_conserved,
            )

            blast_path = blast_dir / f"sample_{i + 1:02d}.blast.tsv"
            self.blast_factory.write_blast_file(records, blast_path)

        # Create shared ANI matrix
        genomes = BlastDataFactory.GENOME_PREFIXES[:8]
        matrix = self.ani_factory.create_matrix(genomes=genomes)

        ani_path = self.temp_dir / "batch.ani.csv"
        self.ani_factory.write_ani_file(matrix, ani_path)

        return blast_dir, ani_path

    def create_empty_blast_file(self) -> Path:
        """Create an empty BLAST file for edge case testing."""
        path = self.temp_dir / "empty.blast.tsv"
        path.write_text("")
        return path

    def create_single_read_dataset(self) -> tuple[Path, Path]:
        """Create a dataset with a single read for edge case testing."""
        records = self.blast_factory.create_known_species_read()

        blast_path = self.temp_dir / "single_read.blast.tsv"
        self.blast_factory.write_blast_file(records, blast_path)

        # Minimal ANI matrix
        genomes = BlastDataFactory.GENOME_PREFIXES[:2]
        matrix = self.ani_factory.create_matrix(genomes=genomes)

        ani_path = self.temp_dir / "single_read.ani.csv"
        self.ani_factory.write_ani_file(matrix, ani_path)

        return blast_path, ani_path

    def create_boundary_novelty_dataset(self) -> tuple[Path, Path]:
        """
        Create dataset with reads at classification boundary values.

        Includes reads at exactly the novelty thresholds (2.0, 5.0, 15.0).
        """
        records = []

        # Boundary at 2.0 (98% identity) - Known vs Novel Species boundary
        records.extend(
            self.blast_factory.create_known_species_read(
                read_id="boundary_2.0",
                pident=98.0,  # novelty = 2.0 exactly
            )
        )

        # Boundary at 5.0 (95% identity) - Novel Species lower bound
        records.extend(
            self.blast_factory.create_novel_species_read(
                read_id="boundary_5.0",
                pident=95.0,  # novelty = 5.0 exactly
            )
        )

        # Boundary at 15.0 (85% identity) - Novel Genus lower bound
        records.extend(
            self.blast_factory.create_novel_genus_read(
                read_id="boundary_15.0",
                pident=85.0,  # novelty = 15.0 exactly
            )
        )

        blast_path = self.temp_dir / "boundary.blast.tsv"
        self.blast_factory.write_blast_file(records, blast_path)

        genomes = BlastDataFactory.GENOME_PREFIXES[:4]
        matrix = self.ani_factory.create_matrix(genomes=genomes)

        ani_path = self.temp_dir / "boundary.ani.csv"
        self.ani_factory.write_ani_file(matrix, ani_path)

        return blast_path, ani_path
