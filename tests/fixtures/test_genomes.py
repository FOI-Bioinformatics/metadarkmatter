"""
Generate realistic test genomes and reads for end-to-end testing.

Creates synthetic bacterial genomes with realistic properties:
- RefSeq-style accessions
- Multiple contigs per genome
- Realistic sequence lengths
- Known sequence relationships for testing ANI/classification
"""

from __future__ import annotations

from pathlib import Path


def generate_test_genome(accession: str, num_contigs: int = 3, contig_length: int = 1000) -> str:
    """Generate a synthetic genome with multiple contigs.

    Args:
        accession: RefSeq accession (e.g., "GCF_000001.1")
        num_contigs: Number of contigs to generate
        contig_length: Length of each contig in bp

    Returns:
        FASTA string with all contigs
    """
    # Use deterministic sequence generation based on accession
    # This ensures consistent sequences across test runs
    seed = sum(ord(c) for c in accession)

    contigs = []
    for i in range(num_contigs):
        contig_id = f"contig_{i+1}"
        # Generate pseudo-random but deterministic sequence
        seq = generate_sequence(seed + i, contig_length)
        contigs.append(f">{accession}|{contig_id}\n{seq}")

    return "\n".join(contigs)


def generate_sequence(seed: int, length: int) -> str:
    """Generate a deterministic DNA sequence.

    Args:
        seed: Random seed for reproducibility
        length: Sequence length in bp

    Returns:
        DNA sequence string (multiline, 80 bp per line)
    """
    import random
    rng = random.Random(seed)
    bases = ['A', 'T', 'G', 'C']

    seq = ''.join(rng.choice(bases) for _ in range(length))

    # Format as 80 bp per line (standard FASTA)
    lines = []
    for i in range(0, len(seq), 80):
        lines.append(seq[i:i+80])

    return '\n'.join(lines)


def generate_read_from_genome(genome_seq: str, read_id: str, start_pos: int, read_length: int = 150,
                               identity: float = 100.0) -> str:
    """Generate a synthetic read from a genome sequence.

    Args:
        genome_seq: Source genome sequence (without newlines)
        read_id: Read identifier
        start_pos: Starting position in genome
        read_length: Read length in bp
        identity: Percent identity to genome (100 = exact match, 95 = 5% mutations)

    Returns:
        FASTA format read
    """
    import random

    # Extract subsequence
    clean_seq = genome_seq.replace('\n', '')
    if start_pos + read_length > len(clean_seq):
        start_pos = len(clean_seq) - read_length

    read_seq = clean_seq[start_pos:start_pos + read_length]

    # Introduce mutations to achieve target identity
    if identity < 100.0:
        mutation_rate = (100.0 - identity) / 100.0
        read_list = list(read_seq)
        bases = ['A', 'T', 'G', 'C']

        for i in range(len(read_list)):
            if random.random() < mutation_rate:
                # Mutate to different base
                current = read_list[i]
                read_list[i] = random.choice([b for b in bases if b != current])

        read_seq = ''.join(read_list)

    return f">{read_id}\n{read_seq}"


def create_test_genome_set(output_dir: Path) -> dict[str, Path]:
    """Create a set of test genomes representing a bacterial family.

    Creates genomes with known relationships:
    - Same species: High ANI (>95%)
    - Same genus: Moderate ANI (85-95%)
    - Different genus: Low ANI (<85%)

    Args:
        output_dir: Directory to write genome files

    Returns:
        Dict mapping genome accession to file path
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    genomes = {
        # Species A - strain 1 (reference)
        "GCF_000001.1": generate_test_genome("GCF_000001.1", num_contigs=3, contig_length=2000),

        # Species A - strain 2 (>96% ANI to strain 1)
        "GCF_000002.1": generate_test_genome("GCF_000002.1", num_contigs=3, contig_length=2000),

        # Species B - same genus (85-95% ANI)
        "GCF_000003.1": generate_test_genome("GCF_000003.1", num_contigs=3, contig_length=2000),

        # Different genus (< 85% ANI)
        "GCF_000004.1": generate_test_genome("GCF_000004.1", num_contigs=4, contig_length=1500),
    }

    paths = {}
    for accession, fasta_content in genomes.items():
        path = output_dir / f"{accession}.fna"
        path.write_text(fasta_content + "\n")
        paths[accession] = path

    return paths


def create_test_reads(genome_paths: dict[str, Path], output_file: Path,
                      reads_per_genome: int = 10) -> Path:
    """Create synthetic reads from test genomes.

    Generates reads with varying identity to test classification:
    - 100% identity: Known species
    - 95-99% identity: Known species (with variation)
    - 90-94% identity: Novel species candidates
    - 80-89% identity: Novel genus candidates

    Args:
        genome_paths: Dict mapping accession to genome file path
        output_file: Output FASTA file for reads
        reads_per_genome: Number of reads to generate per genome

    Returns:
        Path to created reads file
    """
    import random

    reads = []
    read_counter = 1

    for accession, genome_path in genome_paths.items():
        # Load genome sequence
        genome_fasta = genome_path.read_text()
        lines = genome_fasta.split('\n')
        genome_seq = ''.join(line for line in lines if not line.startswith('>'))

        # Generate reads with different identities
        identities = [100.0] * 3 + [98.0] * 3 + [92.0] * 2 + [85.0] * 2

        for i in range(reads_per_genome):
            identity = identities[i % len(identities)]
            start_pos = random.randint(0, len(genome_seq) - 200)

            read = generate_read_from_genome(
                genome_seq,
                read_id=f"read_{read_counter}_{accession}_{int(identity)}pct",
                start_pos=start_pos,
                read_length=150,
                identity=identity
            )
            reads.append(read)
            read_counter += 1

    # Write all reads to file
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text('\n'.join(reads) + '\n')

    return output_file


def create_test_metadata(genome_paths: dict[str, Path], output_file: Path) -> Path:
    """Create genome metadata file for test genomes.

    Args:
        genome_paths: Dict mapping accession to genome file path
        output_file: Output TSV file path

    Returns:
        Path to created metadata file
    """
    # Create realistic GTDB-style metadata
    metadata_lines = [
        "accession\tspecies\tgenus\tfamily\tgtdb_taxonomy"
    ]

    genome_metadata = {
        "GCF_000001.1": {
            "species": "Testbacterium testii",
            "genus": "Testbacterium",
            "family": "Testbacteriaceae",
            "taxonomy": "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Testales;f__Testbacteriaceae;g__Testbacterium;s__Testbacterium testii"
        },
        "GCF_000002.1": {
            "species": "Testbacterium testii",  # Same species
            "genus": "Testbacterium",
            "family": "Testbacteriaceae",
            "taxonomy": "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Testales;f__Testbacteriaceae;g__Testbacterium;s__Testbacterium testii"
        },
        "GCF_000003.1": {
            "species": "Testbacterium alterum",  # Different species, same genus
            "genus": "Testbacterium",
            "family": "Testbacteriaceae",
            "taxonomy": "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Testales;f__Testbacteriaceae;g__Testbacterium;s__Testbacterium alterum"
        },
        "GCF_000004.1": {
            "species": "Differentium novum",  # Different genus
            "genus": "Differentium",
            "family": "Testbacteriaceae",
            "taxonomy": "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Testales;f__Testbacteriaceae;g__Differentium;s__Differentium novum"
        },
    }

    for accession in genome_paths.keys():
        meta = genome_metadata[accession]
        metadata_lines.append(
            f"{accession}\t{meta['species']}\t{meta['genus']}\t{meta['family']}\t{meta['taxonomy']}"
        )

    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text('\n'.join(metadata_lines) + '\n')

    return output_file
