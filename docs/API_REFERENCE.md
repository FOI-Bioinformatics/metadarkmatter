# Metadarkmatter API Reference

Complete API documentation for the metadarkmatter Python package. This reference covers all public interfaces for the ANI-weighted placement uncertainty algorithm and supporting data structures.

**Package Version**: 0.1.0
**Last Updated**: 2026-02-20

---

## Table of Contents

1. [Core Placement Algorithm](#core-placement-algorithm)
2. [Parsers](#parsers)
3. [Data Models](#data-models)
4. [Configuration](#configuration)
5. [Performance Characteristics](#performance-characteristics)
6. [Error Handling](#error-handling)

---

## Core Placement Algorithm

The core module implements the ANI-weighted placement uncertainty algorithm for detecting novel microbial diversity.

### ANIMatrix

High-performance ANI matrix using NumPy arrays with integer indexing.

#### Constructor

```python
ANIMatrix(ani_dict: dict[str, dict[str, float]]) -> ANIMatrix
```

**Description**:
Initializes an ANI matrix from a nested dictionary of genome-to-genome similarity values. The matrix is optimized using NumPy arrays with integer indexing for fast lookups, providing 18x memory reduction and 16x faster access compared to nested dictionaries.

**Parameters**:
- `ani_dict` (dict[str, dict[str, float]]): Nested dictionary where outer keys are genome names and inner dictionaries map genome names to ANI values. Example: `{"GCF_001": {"GCF_002": 85.5, ...}, ...}`

**Returns**:
ANIMatrix instance ready for classification tasks

**Example**:
```python
from metadarkmatter.core.classification.ani_matrix import ANIMatrix

ani_dict = {
    "GCF_000001": {"GCF_000001": 100.0, "GCF_000002": 85.3},
    "GCF_000002": {"GCF_000001": 85.3, "GCF_000002": 100.0},
}
ani_matrix = ANIMatrix(ani_dict)
print(f"Contains {len(ani_matrix)} genomes")
print(f"Memory usage: {ani_matrix.memory_usage_bytes() / 1e6:.1f} MB")
```

**Performance Notes**:
- Initialization time: O(n^2) where n is number of genomes (only occurs once)
- Memory usage: ~4 MB for 1000 genomes (dense matrix)
- Lookup time: O(1) after initial dict lookup

---

#### ANIMatrix.from_file()

```python
@classmethod
from_file(cls, path: Path) -> ANIMatrix
```

**Description**:
Load ANI matrix from a CSV/TSV file using the ANIMatrixParser.

**Parameters**:
- `path` (Path): File path to ANI matrix (CSV or TSV format)

**Returns**:
ANIMatrix instance loaded from file

**Raises**:
- `FileNotFoundError`: If file does not exist
- `ValueError`: If matrix is not valid (non-square, missing genomes, invalid values)

**Example**:
```python
from pathlib import Path
from metadarkmatter.core.classification.ani_matrix import ANIMatrix

ani_matrix = ANIMatrix.from_file(Path("ani_matrix.csv"))
```

---

#### ANIMatrix.genomes

```python
@property
genomes(self) -> set[str]
```

**Description**:
Get the set of all genome names in the matrix.

**Returns**:
Set of genome identifiers

**Example**:
```python
genomes = ani_matrix.genomes
print(f"Total genomes: {len(genomes)}")
print(f"First 5: {sorted(genomes)[:5]}")
```

---

#### ANIMatrix.__len__()

```python
def __len__(self) -> int
```

**Description**:
Get the number of genomes in the matrix.

**Returns**:
Integer count of genomes

**Example**:
```python
num_genomes = len(ani_matrix)
```

---

#### ANIMatrix.get_ani()

```python
def get_ani(self, genome1: str, genome2: str) -> float
```

**Description**:
Get ANI value between two genomes. Uses integer-indexed NumPy array for O(1) access after initial genome name lookup.

**Parameters**:
- `genome1` (str): First genome identifier
- `genome2` (str): Second genome identifier

**Returns**:
ANI value (0-100), or 0.0 if genomes not found. Returns 100.0 if genome1 == genome2 (diagonal).

**Example**:
```python
ani_value = ani_matrix.get_ani("GCF_000001", "GCF_000002")
print(f"ANI between genomes: {ani_value:.1f}%")
```

---

#### ANIMatrix.get_ani_by_idx()

```python
def get_ani_by_idx(self, idx1: int, idx2: int) -> float
```

**Description**:
Get ANI value using pre-computed integer indices. This is the fastest lookup method when genome indices are already known, avoiding string-based dictionary lookups.

**Parameters**:
- `idx1` (int): First genome index
- `idx2` (int): Second genome index

**Returns**:
ANI value (0-100)

**Performance Notes**:
- Fastest method for repeated lookups: pure NumPy array access O(1)
- Use with `get_genome_idx()` to pre-compute indices

**Example**:
```python
# Pre-compute indices for hot path
idx1 = ani_matrix.get_genome_idx("GCF_000001")
idx2 = ani_matrix.get_genome_idx("GCF_000002")

# Fast lookups in loop
for _ in range(1_000_000):
    ani = ani_matrix.get_ani_by_idx(idx1, idx2)  # 10x faster
```

---

#### ANIMatrix.get_genome_idx()

```python
def get_genome_idx(self, genome: str) -> int | None
```

**Description**:
Get integer index for a genome name. Useful for pre-computing indices before bulk operations.

**Parameters**:
- `genome` (str): Genome identifier

**Returns**:
Integer index (0 to num_genomes-1), or None if genome not found

**Example**:
```python
idx = ani_matrix.get_genome_idx("GCF_000001")
if idx is not None:
    print(f"Genome index: {idx}")
```

---

#### ANIMatrix.has_genome()

```python
def has_genome(self, genome: str) -> bool
```

**Description**:
Check if a genome is present in the ANI matrix.

**Parameters**:
- `genome` (str): Genome identifier

**Returns**:
True if genome is in matrix, False otherwise

**Example**:
```python
if ani_matrix.has_genome("GCF_000001"):
    ani = ani_matrix.get_ani("GCF_000001", "GCF_000002")
```

---

#### ANIMatrix.memory_usage_bytes()

```python
def memory_usage_bytes(self) -> int
```

**Description**:
Estimate total memory usage in bytes, including both the NumPy array and index dictionaries.

**Returns**:
Memory usage in bytes

**Example**:
```python
mem_mb = ani_matrix.memory_usage_bytes() / 1e6
print(f"ANI matrix memory: {mem_mb:.1f} MB")
```

---

### ANIWeightedClassifier

Main classifier for detecting novel microbial diversity using ANI-weighted placement. This class provides the core single-read classification algorithm used by higher-level classifiers.

#### Constructor

```python
ANIWeightedClassifier(
    ani_matrix: ANIMatrix,
    config: ScoringConfig | None = None,
) -> ANIWeightedClassifier
```

**Description**:
Initialize the ANI-weighted classifier with an ANI matrix and optional scoring configuration. Uses default scoring thresholds if no config provided.

**Parameters**:
- `ani_matrix` (ANIMatrix): Precomputed ANI matrix for genome comparisons
- `config` (ScoringConfig, optional): Scoring configuration with classification thresholds. Defaults to ScoringConfig() with standard parameters.

**Returns**:
ANIWeightedClassifier instance ready for classification

**Example**:
```python
from metadarkmatter.core.classification.classifiers.base import ANIWeightedClassifier
from metadarkmatter.core.classification.ani_matrix import ANIMatrix
from metadarkmatter.models.config import ScoringConfig

ani_matrix = ANIMatrix.from_file(Path("ani_matrix.csv"))

# Use default configuration
classifier = ANIWeightedClassifier(ani_matrix)

# Or customize scoring thresholds
config = ScoringConfig(
    bitscore_threshold_pct=95.0,
    novelty_novel_species_min=4.0,
    novelty_novel_species_max=20.0,
)
classifier = ANIWeightedClassifier(ani_matrix, config)
```

---

#### ANIWeightedClassifier.classify_read()

```python
def classify_read(self, blast_result: BlastResult) -> ReadClassification | None
```

**Description**:
Classify a single read based on its BLAST hits. Calculates novelty index and placement uncertainty metrics to determine taxonomic classification.

**Parameters**:
- `blast_result` (BlastResult): BLAST results for one read containing all hits sorted by bitscore

**Returns**:
ReadClassification if classification successful, None if no hits

**Algorithm**:
1. Extract top BLAST hit (highest bitscore)
2. Calculate Novelty Index: N = 100 - top_hit_identity
3. Find secondary hits within 95% of top bitscore
4. Calculate Placement Uncertainty: U = 100 - ANI(top_hit, secondary_hits)
5. Classify based on N and U thresholds

**Example**:
```python
from metadarkmatter.models.blast import BlastResult

# Parse BLAST results
blast_result = BlastResult(...)

# Classify single read
classification = classifier.classify_read(blast_result)
if classification:
    print(f"{classification.read_id}: {classification.taxonomic_call}")
    print(f"  Novelty: {classification.novelty_index:.1f}")
    print(f"  Uncertainty: {classification.placement_uncertainty:.1f}")
```

---

### VectorizedClassifier

Fully vectorized classifier using Polars for maximum performance. This is the sole classifier used by the CLI.

#### Constructor

```python
VectorizedClassifier(
    ani_matrix: ANIMatrix,
    config: ScoringConfig | None = None,
) -> VectorizedClassifier
```

**Description**:
Initialize vectorized classifier using Polars' native Rust backend. All operations run in Polars with automatic CPU parallelization, avoiding Python loops entirely.

**Parameters**:
- `ani_matrix` (ANIMatrix): Precomputed ANI matrix
- `config` (ScoringConfig, optional): Scoring configuration (uses defaults if None)

**Returns**:
VectorizedClassifier instance

**Example**:
```python
from metadarkmatter.core.classification.classifiers.vectorized import VectorizedClassifier
from metadarkmatter.core.classification.ani_matrix import ANIMatrix

ani_matrix = ANIMatrix.from_file(Path("ani_matrix.csv"))
vectorized = VectorizedClassifier(ani_matrix)
```

---

#### VectorizedClassifier.classify_file()

```python
def classify_file(self, blast_path: Path) -> pl.DataFrame
```

**Description**:
Classify BLAST file using fully vectorized Polars operations. All operations run in Polars' Rust backend with automatic parallelization.

**Parameters**:
- `blast_path` (Path): Path to BLAST tabular output

**Returns**:
Polars DataFrame with classification results

**Example**:
```python
df = vectorized.classify_file(Path("blast.tsv"))
print(f"Classified {len(df)} reads")
```

---

#### VectorizedClassifier.stream_to_file()

```python
def stream_to_file(
    self,
    blast_path: Path,
    output_path: Path,
    output_format: str = "parquet",
    partition_size: int = 5_000_000,
    progress_callback: Callable[[int, int, float], None] | None = None,
) -> int
```

**Description**:
Stream classification results to file for very large BLAST files. Processes input in partitions to maintain bounded memory usage.

**Parameters**:
- `blast_path` (Path): Path to BLAST tabular output
- `output_path` (Path): Output file path
- `output_format` (str): 'parquet' or 'csv'
- `partition_size` (int): Number of alignments per partition (default: 5,000,000)
- `progress_callback` (Callable, optional): Function called as `progress_callback(rows_processed, total_reads, elapsed_secs)`

**Returns**:
Total number of reads classified

**Memory Usage**:
- Bounded by partition_size
- ~2 GB per 5M alignments (default)

**Example**:
```python
def progress(rows, reads, elapsed):
    print(f"Processed {rows:,} rows -> {reads:,} reads in {elapsed:.1f}s")

total = vectorized.stream_to_file(
    Path("blast.tsv"),
    Path("classifications.parquet"),
    progress_callback=progress
)
```

---

---

## Parsers

Streaming parsers for efficient handling of large bioinformatics files.

### StreamingBlastParser

Memory-efficient streaming parser for BLAST tabular output.

#### Constructor

```python
StreamingBlastParser(blast_path: Path, chunk_size: int = 1_000_000) -> StreamingBlastParser
```

**Description**:
Initialize streaming BLAST parser. Uses Polars with chunked batch reading to process BLAST results without loading entire files into memory.

**Parameters**:
- `blast_path` (Path): Path to BLAST tabular output file (may be gzipped)
- `chunk_size` (int): Number of rows to read per batch (default: 1,000,000). For 100M records, this processes in ~100 batches.

**Returns**:
StreamingBlastParser instance

**Expected Format**:
BLAST tabular output (-outfmt 6):
```
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
```

**Example**:
```python
from pathlib import Path
from metadarkmatter.core.parsers import StreamingBlastParser

parser = StreamingBlastParser(Path("blast.tsv.gz"))
```

---

#### StreamingBlastParser.parse_lazy()

```python
def parse_lazy(self) -> pl.LazyFrame
```

**Description**:
Parse BLAST file into Polars LazyFrame for deferred execution. No data is loaded into memory until explicitly collected.

**Returns**:
Polars LazyFrame with BLAST results

**Example**:
```python
lazy_df = parser.parse_lazy()

# Filter and aggregate without loading entire file
top_hits = (
    lazy_df
    .group_by("qseqid")
    .agg(pl.col("bitscore").max().alias("max_bitscore"))
    .collect(streaming=True)
)
```

---

#### StreamingBlastParser.parse_eager()

```python
def parse_eager(self) -> pl.DataFrame
```

**Description**:
Parse BLAST file into materialized Polars DataFrame. Loads entire file into memory.

**Returns**:
Materialized DataFrame with BLAST results

**Warning**:
Use only for small-to-medium files or when immediate computation is required. For large files, prefer parse_lazy().

**Example**:
```python
df = parser.parse_eager()
```

---

#### StreamingBlastParser.iter_reads()

```python
def iter_reads(self) -> Iterator[BlastResult]
```

**Description**:
Iterate over BLAST results grouped by read ID. Uses true chunked streaming to maintain bounded memory usage.

**Yields**:
BlastResult objects with hits sorted by bitscore descending

**Example**:
```python
for result in parser.iter_reads():
    if result.num_hits >= 10:
        print(f"{result.read_id}: {result.best_hit.genome_name}")
```

---

#### StreamingBlastParser.iter_reads_fast()

```python
def iter_reads_fast(self) -> Iterator[BlastResultFast]
```

**Description**:
High-performance iterator using lightweight NamedTuple data structures. Approximately 50x faster object creation than iter_reads().

**Performance Optimizations**:
- Uses NamedTuples instead of Pydantic models (~50x faster)
- Vectorized genome extraction in Polars (~100x faster)
- Minimal fields (only what classification needs)

**Yields**:
BlastResultFast objects with pre-extracted genome names

**Example**:
```python
for result in parser.iter_reads_fast():
    best = result.hits[0]
    print(f"{result.read_id}: {best.genome_name} @ {best.pident}%")
```

---

#### StreamingBlastParser.get_top_hits_per_read()

```python
def get_top_hits_per_read(self) -> pl.DataFrame
```

**Description**:
Extract only the top BLAST hit for each read using lazy evaluation.

**Returns**:
DataFrame with one row per read (highest bitscore hit only)

**Example**:
```python
top_hits = parser.get_top_hits_per_read()
```

---

#### StreamingBlastParser.get_ambiguous_hits()

```python
def get_ambiguous_hits(self, bitscore_threshold_pct: float = 95.0) -> pl.DataFrame
```

**Description**:
Get all hits within threshold percentage of top bitscore per read.

**Parameters**:
- `bitscore_threshold_pct` (float): Percentage of max bitscore (default: 95%)

**Returns**:
DataFrame with ambiguous hits, including 'max_bitscore' column

**Example**:
```python
ambiguous = parser.get_ambiguous_hits(bitscore_threshold_pct=95.0)
print(f"Total competitive placements: {len(ambiguous)}")
```

---

#### StreamingBlastParser.count_hits_per_read()

```python
def count_hits_per_read(self) -> pl.DataFrame
```

**Description**:
Count number of BLAST hits per read.

**Returns**:
DataFrame with columns: qseqid, hit_count

**Example**:
```python
hit_counts = parser.count_hits_per_read()
mean_hits = hit_counts["hit_count"].mean()
print(f"Mean hits per read: {mean_hits:.1f}")
```

---

#### StreamingBlastParser.get_summary_stats()

```python
def get_summary_stats(self) -> dict[str, float]
```

**Description**:
Calculate summary statistics for BLAST results.

**Returns**:
Dictionary with keys:
- `total_hits` (int): Total number of alignments
- `unique_reads` (int): Number of unique query sequences
- `mean_pident` (float): Mean percent identity
- `mean_bitscore` (float): Mean bit score

**Example**:
```python
stats = parser.get_summary_stats()
print(f"Total hits: {stats['total_hits']:,}")
print(f"Unique reads: {stats['unique_reads']:,}")
print(f"Mean identity: {stats['mean_pident']:.1f}%")
```

---

### ANIMatrixParser

Parser for precomputed ANI (Average Nucleotide Identity) matrices.

#### Constructor

```python
ANIMatrixParser(ani_path: Path) -> ANIMatrixParser
```

**Description**:
Initialize ANI matrix parser.

**Parameters**:
- `ani_path` (Path): Path to ANI matrix file (CSV/TSV)

**Returns**:
ANIMatrixParser instance

**Expected Format**:
- CSV or TSV with genome names as first column and header row
- Symmetric matrix with ANI values (0-100)

**Example**:
```python
from metadarkmatter.core.parsers import ANIMatrixParser

parser = ANIMatrixParser(Path("ani_matrix.csv"))
```

---

#### ANIMatrixParser.parse()

```python
def parse(self) -> pl.DataFrame
```

**Description**:
Parse ANI matrix into Polars DataFrame with validation.

**Returns**:
Polars DataFrame with genome names as first column and ANI values

**Raises**:
- `ValueError`: If matrix is not symmetric, contains invalid values, or row/column names don't match

**Example**:
```python
df = parser.parse()
print(f"Parsed {len(df)} genomes")
```

---

#### ANIMatrixParser.to_dict()

```python
def to_dict(self) -> dict[str, dict[str, float]]
```

**Description**:
Convert ANI matrix to nested dictionary for fast lookups.

**Returns**:
Nested dict: `{genome1: {genome2: ani_value, ...}, ...}`

**Example**:
```python
ani_dict = parser.to_dict()
print(f"ANI between GCF_001 and GCF_002: {ani_dict['GCF_001']['GCF_002']:.1f}")
```

---

### Data Structures

#### BlastHitFast

Lightweight BLAST hit for high-performance internal processing using NamedTuple.

```python
class BlastHitFast(NamedTuple):
    qseqid: str           # Query sequence identifier
    sseqid: str           # Subject sequence identifier
    pident: float         # Percent identity (0-100)
    bitscore: float       # Bit score
    genome_name: str      # Pre-extracted genome name
```

**Description**:
No validation overhead - use only with pre-validated data from Polars. Approximately 50x faster to create than BlastHit.

---

#### BlastResultFast

Lightweight BLAST result for high-performance internal processing.

```python
class BlastResultFast(NamedTuple):
    read_id: str                          # Query sequence identifier
    hits: tuple[BlastHitFast, ...]        # Pre-sorted hits by bitscore
```

---

---

## Data Models

Pydantic models for representing classification results and configuration.

### ReadClassification

Classification result for a single metagenomic read.

```python
class ReadClassification(BaseModel):
    read_id: str
    best_match_genome: str
    top_hit_identity: float              # [0-100]
    novelty_index: float                 # [0-100]
    placement_uncertainty: float         # [0-100]
    num_ambiguous_hits: int
    second_hit_identity: float | None
    identity_gap: float | None
    genus_uncertainty: float | None
    num_genus_hits: int | None
    confidence_score: float | None
    inferred_uncertainty: float | None
    uncertainty_type: str | None
    alignment_quality: float | None
    identity_confidence: float | None
    placement_confidence: float | None
    discovery_score: float | None
    taxonomic_call: TaxonomicCall
```

**Description**:
Contains ANI-weighted placement metrics and final taxonomic call for detecting novel bacterial diversity.

**Core Attributes**:
- `read_id` (str): Original read identifier from FASTQ/BLAST
- `best_match_genome` (str): Genome with highest BLAST bitscore hit
- `top_hit_identity` (float): Percent identity of best BLAST hit (0-100)
- `novelty_index` (float): 100 - top_hit_identity; measures divergence from reference
- `placement_uncertainty` (float): 100 - ANI(top_hit, secondary_hit); measures ambiguity
- `num_ambiguous_hits` (int): Number of hits within 95% of top bitscore
- `taxonomic_call` (TaxonomicCall): Final classification result

**Enhanced Scoring Attributes** (always computed):
- `second_hit_identity` (float | None): Percent identity of second-best hit to different genome
- `identity_gap` (float | None): Gap between best and second-best hit identity
- `genus_uncertainty` (float | None): Uncertainty from genus-level hits (90% bitscore threshold)
- `num_genus_hits` (int | None): Number of hits within genus-level bitscore threshold
- `confidence_score` (float | None): Overall confidence score (0-100) integrating multiple quality factors
- `inferred_uncertainty` (float | None): Inferred placement uncertainty for single-hit reads
- `uncertainty_type` (str | None): Source of uncertainty: 'measured' or 'inferred'
- `alignment_quality` (float | None): Alignment quality score
- `identity_confidence` (float | None): Confidence in the identity measurement
- `placement_confidence` (float | None): Confidence in genome assignment
- `discovery_score` (float | None): Priority score for novel discoveries (None for non-novel)

**Computed Properties**:
- `is_novel` (bool): True if taxonomic_call is NOVEL_SPECIES or NOVEL_GENUS
- `diversity_status` (str): High-level status: 'Known', 'Novel', or 'Uncertain'

**Methods**:
```python
def to_dict(self) -> dict[str, Any]
```
Convert to dictionary for DataFrame creation.

**Example**:
```python
from metadarkmatter.models.classification import ReadClassification, TaxonomicCall

classification = ReadClassification(
    read_id="read_001",
    best_match_genome="GCF_000001",
    top_hit_identity=92.5,
    novelty_index=7.5,
    placement_uncertainty=0.3,
    num_ambiguous_hits=2,
    taxonomic_call=TaxonomicCall.NOVEL_SPECIES,
)

print(f"Novel: {classification.is_novel}")
print(f"Dict: {classification.to_dict()}")
```

---

### TaxonomicCall

Enum representing taxonomic classification categories.

```python
class TaxonomicCall(str, Enum):
    KNOWN_SPECIES = "Known Species"
    NOVEL_SPECIES = "Novel Species"
    NOVEL_GENUS = "Novel Genus"
    SPECIES_BOUNDARY = "Species Boundary"
    AMBIGUOUS = "Ambiguous"
    AMBIGUOUS_WITHIN_GENUS = "Ambiguous Within Genus"
    CONSERVED_REGION = "Conserved Region"
    UNCLASSIFIED = "Unclassified"
    OFF_TARGET = "Off-target"
```

**Description**:
Classification categories based on novelty index and placement uncertainty thresholds from competitive read recruitment methodology.

**Values**:
- `KNOWN_SPECIES`: High identity, low ambiguity (N < 4%, U < 1.5%)
- `NOVEL_SPECIES`: Moderate divergence, low ambiguity (4% <= N < 20%, U < 1.5%)
- `NOVEL_GENUS`: High divergence, low ambiguity (20% <= N <= 25%, U < 1.5%)
- `SPECIES_BOUNDARY`: Moderate placement ambiguity (1.5% <= U < 5%), indicating the read matches multiple closely-related species at the species boundary zone
- `AMBIGUOUS`: High placement ambiguity (U >= 5%) within genus, or identity gap ambiguity
- `AMBIGUOUS_WITHIN_GENUS`: Read hits multiple species within the same genus with similar bitscores but low ANI between them
- `CONSERVED_REGION`: Very high placement ambiguity (U >= 5%) and hits span multiple genera, indicating a highly conserved gene (e.g., 16S rRNA)
- `UNCLASSIFIED`: Does not fit clear biological categories (strain variants or very high divergence N > 25%)
- `OFF_TARGET`: Read has substantially better hits outside the ANI matrix (family validation)

**Example**:
```python
from metadarkmatter.models.classification import TaxonomicCall

call = TaxonomicCall.NOVEL_SPECIES
print(f"Classification: {call.value}")
print(f"Is novel: {call in (TaxonomicCall.NOVEL_SPECIES, TaxonomicCall.NOVEL_GENUS)}")
```

---

### TaxonomicSummary

Summary statistics for a classified environmental sample.

```python
class TaxonomicSummary(BaseModel):
    total_reads: int
    known_species: int
    novel_species: int
    novel_genus: int
    conserved_regions: int
    mean_novelty_index: float
    mean_placement_uncertainty: float
    genome_hit_counts: dict[str, int]
```

**Description**:
Aggregates read-level classifications to provide sample-level insights into microbial diversity and novel taxa detection.

**Attributes**:
- `total_reads` (int): Total reads classified in this sample
- `known_species` (int): Reads classified as known species
- `novel_species` (int): Reads classified as novel species (dark matter)
- `novel_genus` (int): Reads classified as novel genus (dark matter)
- `conserved_regions` (int): Reads in conserved/ambiguous regions
- `mean_novelty_index` (float): Mean novelty index across all reads
- `mean_placement_uncertainty` (float): Mean placement uncertainty
- `genome_hit_counts` (dict[str, int]): Read counts per reference genome (top 50)

**Computed Properties**:
- `known_species_pct` (float): Percentage of reads classified as known species
- `novel_species_pct` (float): Percentage of reads classified as novel species
- `novel_genus_pct` (float): Percentage of reads classified as novel genus
- `novel_diversity_pct` (float): Total percentage of novel diversity detected

**Methods**:
```python
def to_json(self, path: Path) -> None
```
Write summary to JSON file.

```python
@classmethod
def from_json(cls, path: Path) -> TaxonomicSummary
```
Load summary from JSON file.

**Example**:
```python
from metadarkmatter.models.classification import TaxonomicSummary

summary = TaxonomicSummary(
    total_reads=10000,
    known_species=7500,
    novel_species=1200,
    novel_genus=800,
    conserved_regions=500,
    mean_novelty_index=4.5,
    mean_placement_uncertainty=0.8,
    genome_hit_counts={"GCF_001": 5000, "GCF_002": 2500},
)

print(f"Novel diversity: {summary.novel_diversity_pct:.1f}%")
summary.to_json(Path("summary.json"))
```

---

### BlastHit

Single BLAST alignment hit from tabular output.

```python
class BlastHit(BaseModel):
    qseqid: str              # Query sequence identifier
    sseqid: str              # Subject sequence identifier
    pident: float            # Percent identity [0-100]
    length: int              # Alignment length (bp)
    mismatch: int            # Number of mismatches
    gapopen: int             # Number of gap openings
    qstart: int              # Start position in query
    qend: int                # End position in query
    sstart: int              # Start position in subject
    send: int                # End position in subject
    evalue: float            # Expectation value
    bitscore: float          # Bit score
```

**Description**:
Represents one alignment between a metagenomic read and a reference genome sequence. Multiple hits per read are expected in competitive recruitment scenarios.

**Properties**:
```python
@property
def genome_name(self) -> str
```
Extract genome identifier from subject sequence ID. Handles common RefSeq/GenBank formats and custom pipe-delimited formats.

**Methods**:
```python
@classmethod
def from_blast_line(cls, line: str) -> BlastHit
```
Parse a single line from BLAST tabular output (outfmt 6).

**Example**:
```python
from metadarkmatter.models.blast import BlastHit

# Parse from BLAST line
blast_line = "read_001\tGCF_000001.1_contig_1\t95.5\t100\t4\t1\t1\t100\t1\t100\t1e-30\t150.5"
hit = BlastHit.from_blast_line(blast_line)
print(f"Genome: {hit.genome_name}")
print(f"Identity: {hit.pident:.1f}%")
```

---

### BlastResult

Collection of BLAST hits for a single read.

```python
class BlastResult(BaseModel):
    read_id: str
    hits: tuple[BlastHit, ...]
```

**Description**:
Groups all alignment hits for one query sequence, sorted by bitscore in descending order.

**Properties**:
```python
@property
def best_hit(self) -> BlastHit | None
```
Get the highest-scoring BLAST hit.

```python
@property
def num_hits(self) -> int
```
Total number of BLAST hits for this read.

```python
@property
def top_bitscore(self) -> float
```
Bitscore of the best hit.

**Methods**:
```python
def iter_ambiguous_hits(self, threshold_pct: float = 95.0) -> Iterator[BlastHit]
```
Iterate over hits within threshold percentage of top bitscore with early termination.

```python
def get_ambiguous_hits(self, threshold_pct: float = 95.0) -> list[BlastHit]
```
Get hits within threshold percentage (returns list).

```python
def sorted_by_bitscore(self) -> BlastResult
```
Return new BlastResult sorted by bitscore (already sorted on creation).

**Example**:
```python
from metadarkmatter.models.blast import BlastResult, BlastHit

hits = (
    BlastHit.from_blast_line(line1),
    BlastHit.from_blast_line(line2),
)

result = BlastResult(read_id="read_001", hits=hits)
print(f"Best hit: {result.best_hit.genome_name}")
print(f"Total hits: {result.num_hits}")

# Get ambiguous hits
for hit in result.iter_ambiguous_hits(threshold_pct=95.0):
    print(f"  {hit.genome_name}: {hit.bitscore:.1f}")
```

---

---

## Configuration

Pydantic configuration models for pipeline parameters.

### ScoringConfig

Configuration for ANI-weighted placement scoring thresholds.

```python
class ScoringConfig(BaseModel):
    alignment_mode: Literal["nucleotide", "protein"] = "nucleotide"
    bitscore_threshold_pct: float = 95.0
    genus_bitscore_threshold_pct: float = 90.0
    novelty_known_max: float = 4.0
    novelty_novel_species_min: float = 4.0
    novelty_novel_species_max: float = 20.0
    novelty_novel_genus_min: float = 20.0
    novelty_novel_genus_max: float = 25.0
    uncertainty_known_max: float = 1.5
    uncertainty_novel_species_max: float = 1.5
    uncertainty_novel_genus_max: float = 1.5
    uncertainty_conserved_min: float = 5.0
    genus_uncertainty_ambiguous_min: float = 10.0
    identity_gap_ambiguous_max: float = 2.0
    confidence_threshold: float = 50.0
    min_alignment_length: int = 100
    min_alignment_fraction: float = 0.3
    coverage_weight_mode: Literal["none", "linear", "log", "sigmoid"] = "linear"
    coverage_weight_strength: float = 0.5
    uncertainty_mode: Literal["max", "second"] = "second"
    single_hit_uncertainty_threshold: float = 10.0
    aai_genus_boundary_low: float = 58.0
    aai_genus_boundary_high: float = 65.0
    use_aai_for_genus: bool = True
    target_family: str | None = None
    family_ratio_threshold: float = 0.8
```

**Description**:
Thresholds define boundaries between different taxonomic classifications based on novelty index and placement uncertainty metrics.

**Nucleotide Mode Thresholds (default)**:
- Known Species: N < 4%, U < 1.5%
- Novel Species: 4% <= N < 20%, U < 1.5%
- Novel Genus: 20% <= N <= 25%, U < 1.5%

**Protein Mode Thresholds** (`alignment_mode="protein"`):
- Known Species: N < 10%, U < 5%
- Novel Species: 10% <= N < 25%, U < 5%
- Novel Genus: 25% <= N <= 40%, U < 5%

**Key Parameters**:
- `alignment_mode` (Literal): 'nucleotide' for BLASTN, 'protein' for BLASTX results
- `bitscore_threshold_pct` (float): Percentage of top bitscore for ambiguous hit detection (default: 95%)
- `genus_bitscore_threshold_pct` (float): Percentage threshold for genus-level ambiguity (default: 90%)
- `novelty_known_max` (float): Maximum novelty index for known species (default: 4.0)
- `novelty_novel_species_min` (float): Minimum novelty index for novel species (default: 4.0)
- `novelty_novel_species_max` (float): Maximum novelty index for novel species (default: 20.0)
- `novelty_novel_genus_min` (float): Minimum novelty index for novel genus (default: 20.0)
- `novelty_novel_genus_max` (float): Maximum novelty index for novel genus (default: 25.0)
- `uncertainty_known_max` (float): Maximum placement uncertainty for known species (default: 1.5)
- `uncertainty_novel_species_max` (float): Maximum uncertainty for novel species (default: 1.5)
- `uncertainty_novel_genus_max` (float): Maximum uncertainty for novel genus (default: 1.5)
- `uncertainty_conserved_min` (float): Minimum uncertainty for conserved regions (default: 5.0)
- `genus_uncertainty_ambiguous_min` (float): Minimum genus uncertainty for ambiguous-within-genus flag (default: 10.0)
- `identity_gap_ambiguous_max` (float): Maximum identity gap before ambiguity flag (default: 2.0)
- `confidence_threshold` (float): Minimum confidence score for high-confidence classification (default: 50.0)
- `min_alignment_length` (int): Minimum alignment length in bp (default: 100)
- `min_alignment_fraction` (float): Minimum fraction of read aligned (default: 0.3)
- `coverage_weight_mode` (Literal): Coverage weighting mode for hit selection (default: 'linear')
- `coverage_weight_strength` (float): Magnitude of coverage effect on scoring (default: 0.5)
- `uncertainty_mode` (Literal): How to calculate placement uncertainty - 'max' or 'second' (default: 'second')
- `single_hit_uncertainty_threshold` (float): Inferred uncertainty threshold for single-hit reads (default: 10.0)
- `aai_genus_boundary_low` (float): Lower AAI boundary for genus classification (default: 58.0)
- `aai_genus_boundary_high` (float): Upper AAI boundary for genus classification (default: 65.0)
- `use_aai_for_genus` (bool): Use AAI instead of ANI for genus-level decisions when available (default: True)
- `target_family` (str | None): Target family for off-target detection (default: None)
- `family_ratio_threshold` (float): Bitscore ratio threshold for off-target reads (default: 0.8)

**Validation**:
Configuration is validated to ensure threshold relationships:
- novelty_novel_species_min >= novelty_known_max
- novelty_novel_genus_min == novelty_novel_species_max (continuous boundaries)
- uncertainty_novel_species_max <= uncertainty_novel_genus_max

**Methods**:
```python
def get_effective_thresholds(self) -> dict[str, float | tuple[float, ...]]
```
Get effective thresholds for the current alignment mode. In protein mode, returns wider thresholds from `protein_constants.py`.

```python
@classmethod
def for_protein_mode(cls) -> ScoringConfig
```
Factory method for protein-level classification configuration.

**Example**:
```python
from metadarkmatter.models.config import ScoringConfig

# Default configuration (nucleotide mode)
config = ScoringConfig()

# Custom thresholds
config = ScoringConfig(
    bitscore_threshold_pct=95.0,
    novelty_novel_species_min=4.0,
    novelty_novel_species_max=20.0,
    uncertainty_novel_species_max=1.5,
)

# Protein mode
protein_config = ScoringConfig.for_protein_mode()
```

---

### BlastConfig

Configuration for BLAST execution parameters.

```python
class BlastConfig(BaseModel):
    num_threads: int = 4
    max_target_seqs: int = 5000
    evalue: float = 1e-5
    perc_identity: float = 0.0
    outfmt: str = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
```

**Parameters**:
- `num_threads` (int): Number of CPU threads for BLAST (default: 4)
- `max_target_seqs` (int): Maximum number of aligned sequences to keep (default: 5000)
- `evalue` (float): Expectation value threshold (default: 1e-5)
- `perc_identity` (float): Minimum percent identity (default: 0.0)
- `outfmt` (str): BLAST output format specification

---

---

## Performance Characteristics

### Memory Usage

| Component | Dataset Size | Memory | Notes |
|-----------|--------------|--------|-------|
| ANIMatrix | 1000 genomes | 4 MB | Dense NumPy matrix |
| ANIMatrix | 10000 genomes | 400 MB | Dense matrix for large sets |
| BLAST Parsing | 100M records | ~16 GB | Polars streaming keeps memory constant |
| Full Classification | 100M records | ~32 GB | Streaming writes to disk in chunks |

### Classification Speed

| Method | Dataset | Speed | Notes |
|--------|---------|-------|-------|
| VectorizedClassifier.classify_file() | 10M reads | ~3 sec | Fully vectorized, Polars Rust backend |
| VectorizedClassifier.stream_to_file() | 100M reads | ~5 min | Partitioned processing, bounded memory |

### Bitscore Threshold Effect

The `bitscore_threshold_pct` parameter affects performance and results:

- **95% (default)**: Most reads have 1-5 ambiguous hits, moderate uncertainty
- **90%**: More competitive placements, higher uncertainty values
- **85%**: Very competitive, many reads marked conserved region
- **50%**: Maximum sensitivity, many false ambiguities

Recommended: Keep at 95% for standard competitive recruitment.

---

## Error Handling

### Common Errors and Solutions

#### FileNotFoundError

```python
# BLAST file not found
FileNotFoundError: BLAST file not found: /path/to/blast.tsv
```

**Solution**: Ensure BLAST file path is correct and file exists.

```python
from pathlib import Path
blast_path = Path("blast.tsv")
assert blast_path.exists(), f"BLAST file not found: {blast_path}"
```

#### ValueError: ANI matrix is not square

```python
ValueError: ANI matrix is not square: 1000 rows, 999 value columns
```

**Solution**: Check that ANI matrix has same number of rows and columns (genome names must match both axes).

#### ValueError: ANI matrix contains null/missing values

```python
ValueError: ANI matrix contains null/missing values
```

**Solution**: Replace NaN values in ANI matrix with valid ANI values before loading.

#### ValueError: threshold_pct validation failed

```python
ValueError: bitscore_threshold_pct must be between 0 and 100
```

**Solution**: Ensure threshold percentage is in valid range.

```python
config = ScoringConfig(bitscore_threshold_pct=95.0)  # Valid
```

### Validation Examples

```python
from metadarkmatter.core.classification.ani_matrix import ANIMatrix
from metadarkmatter.core.classification.classifiers.base import ANIWeightedClassifier
from metadarkmatter.models.config import ScoringConfig
from pathlib import Path

# Validate ANI matrix
try:
    ani_matrix = ANIMatrix.from_file(Path("ani.csv"))
except ValueError as e:
    print(f"Invalid ANI matrix: {e}")
    exit(1)

# Validate configuration
try:
    config = ScoringConfig(novelty_novel_species_min=3.0)  # Invalid: must be >= 4.0
except ValueError as e:
    print(f"Invalid config: {e}")
    exit(1)

# Classify with VectorizedClassifier (recommended)
from metadarkmatter.core.classification.classifiers.vectorized import VectorizedClassifier

vectorized = VectorizedClassifier(ani_matrix)
df = vectorized.classify_file(Path("blast.tsv"))
print(f"Successfully classified {len(df)} reads")
```

---

## Quick Reference

### Most Common Workflows

#### 1. Classify BLAST file with VectorizedClassifier

```python
from pathlib import Path
from metadarkmatter.core.classification.ani_matrix import ANIMatrix
from metadarkmatter.core.classification.classifiers.vectorized import VectorizedClassifier

# Load ANI matrix
ani_matrix = ANIMatrix.from_file(Path("ani_matrix.csv"))

# Create classifier and classify
vectorized = VectorizedClassifier(ani_matrix)
df = vectorized.classify_file(Path("blast.tsv"))
print(f"Classified {len(df):,} reads")

# Save results
df.write_csv("classifications.csv")
```

#### 2. Stream large files to disk

```python
# Use stream_to_file for 100M+ reads
total = vectorized.stream_to_file(
    Path("blast.tsv"),
    Path("classifications.parquet"),
    partition_size=5_000_000
)
print(f"Classified {total:,} reads to Parquet")
```

#### 3. Analyze classification results

```python
import polars as pl

# Load classifications
df = pl.read_csv("classifications.csv")

# Summary statistics
summary = {
    "total_reads": len(df),
    "novel_reads": df.filter(pl.col("is_novel")).height,
    "novel_pct": 100 * df["is_novel"].sum() / len(df),
    "mean_novelty": df["novelty_index"].mean(),
}

print(f"Novel reads: {summary['novel_reads']:,} ({summary['novel_pct']:.1f}%)")
print(f"Mean novelty: {summary['mean_novelty']:.2f}")
```

#### 4. Process with custom scoring thresholds

```python
from metadarkmatter.models.config import ScoringConfig

# More stringent novel species threshold
config = ScoringConfig(
    novelty_novel_species_min=7.0,    # Require higher divergence
    novelty_known_max=7.0,            # Must match for continuous boundaries
    uncertainty_novel_species_max=1.0,  # Require lower ambiguity
)

vectorized = VectorizedClassifier(ani_matrix, config)
df = vectorized.classify_file(Path("blast.tsv"))
```

#### 5. Programmatic single-read classification

```python
from metadarkmatter.core.classification.classifiers.base import ANIWeightedClassifier
from metadarkmatter.models.blast import BlastResult

# ANIWeightedClassifier exposes single-read classification
classifier = ANIWeightedClassifier(ani_matrix)
classification = classifier.classify_read(blast_result)
if classification:
    print(f"{classification.read_id}: {classification.taxonomic_call.value}")
```

---

## Further Reading

- See `CLAUDE.md` for project background and methodology
- See `docs/METHODS.md` for the scientific methods documentation
- See `docs/REFERENCE.md` for CLI reference and algorithm details
- See workflow examples in `src/metadarkmatter/cli/` for CLI usage patterns
