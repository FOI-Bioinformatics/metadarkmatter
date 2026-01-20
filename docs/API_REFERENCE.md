# Metadarkmatter API Reference

Complete API documentation for the metadarkmatter Python package. This reference covers all public interfaces for the ANI-weighted placement uncertainty algorithm and supporting data structures.

**Package Version**: 0.1.0
**Last Updated**: 2025-12-31

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
from metadarkmatter.core.ani_placement import ANIMatrix

ani_dict = {
    "GCF_000001": {"GCF_000001": 100.0, "GCF_000002": 85.3},
    "GCF_000002": {"GCF_000001": 85.3, "GCF_000002": 100.0},
}
ani_matrix = ANIMatrix(ani_dict)
print(f"Contains {len(ani_matrix)} genomes")
print(f"Memory usage: {ani_matrix.memory_usage_bytes() / 1e6:.1f} MB")
```

**Performance Notes**:
- Initialization time: O(n²) where n is number of genomes (only occurs once)
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
from metadarkmatter.core.ani_placement import ANIMatrix

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

Main classifier for detecting novel microbial diversity using ANI-weighted placement.

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
from metadarkmatter.core.ani_placement import ANIMatrix, ANIWeightedClassifier
from metadarkmatter.models.config import ScoringConfig

ani_matrix = ANIMatrix.from_file(Path("ani_matrix.csv"))

# Use default configuration
classifier = ANIWeightedClassifier(ani_matrix)

# Or customize scoring thresholds
config = ScoringConfig(
    bitscore_threshold_pct=95.0,
    novelty_novel_species_min=5.0,
    novelty_novel_species_max=15.0,
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

#### ANIWeightedClassifier.classify_blast_file()

```python
def classify_blast_file(
    self,
    blast_path: Path,
) -> Iterator[ReadClassification]
```

**Description**:
Classify all reads from a BLAST file using memory-efficient streaming. Yields results one at a time without loading entire file into memory.

**Parameters**:
- `blast_path` (Path): Path to BLAST tabular output

**Yields**:
ReadClassification objects for each successfully classified read

**Memory Usage**:
- Constant memory regardless of file size
- Processes in configurable chunks (default: 1M records per chunk)

**Example**:
```python
# Stream results without loading entire file
for classification in classifier.classify_blast_file(Path("blast.tsv")):
    if classification.is_novel:
        print(f"Novel: {classification.read_id} in {classification.best_match_genome}")

# Or collect with early termination
novel_only = [c for c in classifier.classify_blast_file(Path("blast.tsv"))
              if c.is_novel]
```

---

#### ANIWeightedClassifier.classify_to_dataframe()

```python
def classify_to_dataframe(
    self,
    blast_path: Path,
) -> pl.DataFrame
```

**Description**:
Classify all reads and return results as a Polars DataFrame. Use for medium-sized files where full results fit in memory.

**Parameters**:
- `blast_path` (Path): Path to BLAST tabular output

**Returns**:
Polars DataFrame with columns:
- `read_id` (Utf8): Read identifier
- `best_match_genome` (Utf8): Top hit genome name
- `top_hit_identity` (Float64): Percent identity of best hit (0-100)
- `novelty_index` (Float64): Divergence metric (0-100)
- `placement_uncertainty` (Float64): Ambiguity metric (0-100)
- `num_ambiguous_hits` (Int64): Number of competitive hits
- `taxonomic_call` (Utf8): Classification result
- `is_novel` (Boolean): True if Novel Species or Novel Genus

**Example**:
```python
df = classifier.classify_to_dataframe(Path("blast.tsv"))

# Analyze results
novel_count = (df["is_novel"]).sum()
mean_novelty = df["novelty_index"].mean()
print(f"Novel reads: {novel_count} / {len(df)}")
print(f"Mean novelty: {mean_novelty:.2f}")

# Filter and export
df.filter(pl.col("is_novel")).write_csv("novel_only.csv")
```

---

#### ANIWeightedClassifier.write_classifications()

```python
def write_classifications(
    self,
    blast_path: Path,
    output_path: Path,
    output_format: str = "csv",
) -> int
```

**Description**:
Classify reads and write results directly to file. Useful for large files or when further analysis is not needed in memory.

**Parameters**:
- `blast_path` (Path): Path to BLAST tabular output
- `output_path` (Path): Path for output file
- `output_format` (str): Output format - "csv" or "parquet" (default: "csv"). Parquet is 10x smaller and 10x faster for large datasets.

**Returns**:
Number of reads classified

**Example**:
```python
# Write to CSV (human-readable)
num_reads = classifier.write_classifications(
    Path("blast.tsv"),
    Path("classifications.csv"),
    output_format="csv"
)
print(f"Classified {num_reads} reads")

# Write to Parquet (efficient for 100M+ reads)
num_reads = classifier.write_classifications(
    Path("blast.tsv"),
    Path("classifications.parquet"),
    output_format="parquet"
)
```

---

#### ANIWeightedClassifier.classify_read_fast()

```python
def classify_read_fast(
    self,
    result: BlastResultFast,
) -> dict | None
```

**Description**:
Classify a read using lightweight data structures. Approximately 10x faster than classify_read() for large-scale processing.

**Performance Optimizations**:
- Uses BlastResultFast (NamedTuple) instead of Pydantic model
- Pre-extracted genome names avoid regex in hot path
- Returns dict directly, bypassing Pydantic on output

**Parameters**:
- `result` (BlastResultFast): BLAST result with pre-extracted genome names

**Returns**:
Dictionary with classification results, or None if no hits:
```python
{
    "read_id": str,
    "best_match_genome": str,
    "top_hit_identity": float,
    "novelty_index": float,
    "placement_uncertainty": float,
    "num_ambiguous_hits": int,
    "taxonomic_call": str,
    "is_novel": bool,
}
```

**Example**:
```python
parser = StreamingBlastParser(Path("blast.tsv"))
for result in parser.iter_reads_fast():
    classification = classifier.classify_read_fast(result)
    if classification and classification["is_novel"]:
        print(f"Novel: {classification['read_id']}")
```

---

#### ANIWeightedClassifier.classify_blast_file_fast()

```python
def classify_blast_file_fast(
    self,
    blast_path: Path,
) -> Iterator[dict]
```

**Description**:
High-performance streaming classification using lightweight data structures and vectorized genome extraction. Approximately 10x faster than classify_blast_file().

**Performance Optimizations**:
- Uses iter_reads_fast() with NamedTuples (~50x faster object creation)
- Vectorized genome extraction in Polars (~100x faster)
- Returns dicts instead of Pydantic models

**Parameters**:
- `blast_path` (Path): Path to BLAST tabular output

**Yields**:
Dict with classification results for each read

**Example**:
```python
# High-performance streaming for large files
for classification in classifier.classify_blast_file_fast(Path("blast.tsv")):
    if classification["is_novel"]:
        print(f"Novel: {classification['read_id']}")
```

---

#### ANIWeightedClassifier.classify_to_dataframe_fast()

```python
def classify_to_dataframe_fast(
    self,
    blast_path: Path,
) -> pl.DataFrame
```

**Description**:
Classify reads and return Polars DataFrame directly using fast methods. This is the fastest method for processing large BLAST files that fit in memory.

**Performance vs classify_to_dataframe()**:
- ~10x faster for large files (100M+ reads)
- ~5x less memory usage
- Streaming-compatible

**Parameters**:
- `blast_path` (Path): Path to BLAST tabular output

**Returns**:
Polars DataFrame with classification results (same schema as classify_to_dataframe)

**Example**:
```python
# Fast processing for 10M+ read files
df = classifier.classify_to_dataframe_fast(Path("blast.tsv"))
print(f"Classified {len(df)} reads in ~30 seconds")
```

---

#### ANIWeightedClassifier.stream_to_file_fast()

```python
def stream_to_file_fast(
    self,
    blast_path: Path,
    output_path: Path,
    output_format: str = "parquet",
    chunk_size: int = 100_000,
    progress_callback: Callable[[int, int], None] | None = None,
) -> int
```

**Description**:
Stream classification results directly to disk in chunks. Suitable for very large files (100M+ reads) where memory is constrained.

**Parameters**:
- `blast_path` (Path): Path to BLAST tabular output
- `output_path` (Path): Output file path
- `output_format` (str): 'parquet' or 'csv' (default: "parquet")
- `chunk_size` (int): Number of classifications per batch (default: 100,000)
- `progress_callback` (Callable, optional): Function called as `progress_callback(processed, total)` for progress reporting

**Returns**:
Total number of reads classified

**Memory Usage**:
- Bounded by chunk_size: ~200 bytes per record
- Default chunk_size (100K) uses ~20 MB per batch

**Example**:
```python
def show_progress(done, total):
    print(f"Processed {done:,} / {total:,} reads...")

total = classifier.stream_to_file_fast(
    Path("blast.tsv"),
    Path("classifications.parquet"),
    chunk_size=100_000,
    progress_callback=show_progress
)
print(f"Total classified: {total:,} reads")
```

---

### ParallelClassifier

Parallel classifier for very large BLAST files using multiprocessing.

#### Constructor

```python
ParallelClassifier(
    ani_matrix: ANIMatrix,
    config: ScoringConfig | None = None,
    num_workers: int | None = None,
    chunk_size: int = 50_000,
) -> ParallelClassifier
```

**Description**:
Initialize parallel classifier for multicore processing. The ANI matrix is shared across processes using shared memory, minimizing memory overhead.

**Parameters**:
- `ani_matrix` (ANIMatrix): Precomputed ANI matrix
- `config` (ScoringConfig, optional): Scoring configuration (uses defaults if None)
- `num_workers` (int, optional): Number of worker processes (default: CPU count - 1)
- `chunk_size` (int): Reads per worker chunk (default: 50,000)

**Returns**:
ParallelClassifier instance

**Example**:
```python
ani_matrix = ANIMatrix.from_file(Path("ani_matrix.csv"))
parallel = ParallelClassifier(ani_matrix, num_workers=8, chunk_size=50_000)
```

---

#### ParallelClassifier.classify_file()

```python
def classify_file(
    self,
    blast_path: Path,
    progress_callback: Callable[[int], None] | None = None,
) -> pl.DataFrame
```

**Description**:
Classify BLAST file using parallel processing across multiple CPU cores.

**Parameters**:
- `blast_path` (Path): Path to BLAST tabular output
- `progress_callback` (Callable, optional): Function called as `progress_callback(chunks_completed)` for progress

**Returns**:
Polars DataFrame with classification results

**Performance**:
- Ideal for files with 10M+ reads where single-threaded processing becomes a bottleneck
- Memory efficient through chunked processing
- Automatic CPU core utilization

**Example**:
```python
def progress(chunks):
    print(f"Completed {chunks} chunks...")

df = parallel.classify_file(
    Path("blast.tsv"),
    progress_callback=progress
)
print(f"Classified {len(df)} reads using {parallel.num_workers} workers")
```

---

#### ParallelClassifier.stream_to_file()

```python
def stream_to_file(
    self,
    blast_path: Path,
    output_path: Path,
    output_format: str = "parquet",
    progress_callback: Callable[[int], None] | None = None,
) -> int
```

**Description**:
Classify and stream results to file using parallel processing. Combines parallel classification with chunked output for maximum scalability.

**Parameters**:
- `blast_path` (Path): Path to BLAST tabular output
- `output_path` (Path): Output file path
- `output_format` (str): 'parquet' or 'csv'
- `progress_callback` (Callable, optional): Function called as `progress_callback(chunks_completed)`

**Returns**:
Total reads classified

**Example**:
```python
total = parallel.stream_to_file(
    Path("blast.tsv"),
    Path("classifications.parquet")
)
```

---

### VectorizedClassifier

Fully vectorized classifier using Polars for maximum performance.

#### Constructor

```python
VectorizedClassifier(
    ani_matrix: ANIMatrix,
    config: ScoringConfig | None = None,
) -> VectorizedClassifier
```

**Description**:
Initialize vectorized classifier using Polars' native Rust backend. All operations run in Polars with automatic CPU parallelization, avoiding Python loops entirely.

**Performance vs Other Methods**:
- 5-10x faster than classify_to_dataframe_fast() for large files
- Uses Polars' internal parallelism (no multiprocessing overhead)
- Memory-efficient through lazy evaluation

**Parameters**:
- `ani_matrix` (ANIMatrix): Precomputed ANI matrix
- `config` (ScoringConfig, optional): Scoring configuration (uses defaults if None)

**Returns**:
VectorizedClassifier instance

**Example**:
```python
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

### SparseANIMatrix

Sparse ANI matrix for very large genome sets (10K+ genomes).

#### Constructor

```python
SparseANIMatrix(
    ani_dict: dict[str, dict[str, float]],
    default_ani: float = 70.0,
    min_ani: float = 75.0,
) -> SparseANIMatrix
```

**Description**:
Initialize sparse ANI matrix for large genome databases. Stores only non-zero ANI values, reducing memory from 400 MB (dense) to ~20 MB (5% density) for 10K genomes.

**Parameters**:
- `ani_dict` (dict[str, dict[str, float]]): Nested dictionary of ANI values
- `default_ani` (float): Default ANI for missing pairs (default: 70.0)
- `min_ani` (float): Minimum ANI value to store; below this uses default (default: 75.0)

**Returns**:
SparseANIMatrix instance

**Example**:
```python
from metadarkmatter.core.ani_placement import SparseANIMatrix

# Use sparse representation for 10K+ genomes
sparse_matrix = SparseANIMatrix(ani_dict, default_ani=70.0)
classifier = ANIWeightedClassifier(sparse_matrix)
```

---

#### SparseANIMatrix.from_file()

```python
@classmethod
from_file(cls, path: Path, **kwargs) -> SparseANIMatrix
```

**Description**:
Load sparse ANI matrix from file.

**Parameters**:
- `path` (Path): File path to ANI matrix
- `**kwargs`: Additional arguments passed to constructor (default_ani, min_ani)

**Returns**:
SparseANIMatrix instance

**Example**:
```python
sparse_matrix = SparseANIMatrix.from_file(Path("ani_matrix.csv"), default_ani=70.0)
```

---

#### SparseANIMatrix.get_ani()

```python
def get_ani(self, genome1: str, genome2: str) -> float
```

**Description**:
Get ANI between two genomes, using default ANI for missing pairs.

**Parameters**:
- `genome1` (str): First genome identifier
- `genome2` (str): Second genome identifier

**Returns**:
ANI value (0-100), or default_ani if genomes not found

**Example**:
```python
ani = sparse_matrix.get_ani("GCF_001", "GCF_002")  # Returns default if not stored
```

---

#### SparseANIMatrix.density()

```python
def density(self) -> float
```

**Description**:
Calculate matrix density (fraction of non-default values stored).

**Returns**:
Density as float between 0 and 1

**Example**:
```python
density = sparse_matrix.density()
print(f"Matrix density: {density:.1%}")  # e.g., "Matrix density: 5.2%"
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
    taxonomic_call: TaxonomicCall
```

**Description**:
Contains ANI-weighted placement metrics and final taxonomic call for detecting novel bacterial diversity.

**Attributes**:
- `read_id` (str): Original read identifier from FASTQ/BLAST
- `best_match_genome` (str): Genome with highest BLAST bitscore hit
- `top_hit_identity` (float): Percent identity of best BLAST hit (0-100)
- `novelty_index` (float): 100 - top_hit_identity; measures divergence from reference
- `placement_uncertainty` (float): 100 - ANI(top_hit, secondary_hit); measures ambiguity
- `num_ambiguous_hits` (int): Number of hits within 95% of top bitscore
- `taxonomic_call` (TaxonomicCall): Final classification result

**Properties**:
- `is_novel` (bool, computed): True if taxonomic_call is NOVEL_SPECIES or NOVEL_GENUS

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
    CONSERVED_REGION = "Conserved Region"
```

**Description**:
Classification categories based on novelty index and placement uncertainty thresholds from competitive read recruitment methodology.

**Values**:
- `KNOWN_SPECIES`: Low novelty, low uncertainty (N < 2, U < 0.5)
- `NOVEL_SPECIES`: Moderate novelty, low uncertainty (5 <= N <= 15, U < 0.5)
- `NOVEL_GENUS`: High novelty, moderate uncertainty (15 <= N <= 25, U < 2)
- `CONSERVED_REGION`: High uncertainty or no clear classification (U >= 5)

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
    bitscore_threshold_pct: float = 95.0
    novelty_known_max: float = 2.0
    novelty_novel_species_min: float = 5.0
    novelty_novel_species_max: float = 15.0
    novelty_novel_genus_min: float = 15.0
    novelty_novel_genus_max: float = 25.0
    uncertainty_known_max: float = 0.5
    uncertainty_novel_species_max: float = 0.5
    uncertainty_novel_genus_max: float = 2.0
    uncertainty_conserved_min: float = 5.0
```

**Description**:
Thresholds define boundaries between different taxonomic classifications based on novelty index and placement uncertainty metrics.

**Parameters**:
- `bitscore_threshold_pct` (float): Percentage of top bitscore for ambiguous hits (default: 95%)
- `novelty_known_max` (float): Maximum novelty index for known species (default: 2.0)
- `novelty_novel_species_min` (float): Minimum novelty index for novel species (default: 5.0)
- `novelty_novel_species_max` (float): Maximum novelty index for novel species (default: 15.0)
- `novelty_novel_genus_min` (float): Minimum novelty index for novel genus (default: 15.0)
- `novelty_novel_genus_max` (float): Maximum novelty index for novel genus (default: 25.0)
- `uncertainty_known_max` (float): Maximum placement uncertainty for known species (default: 0.5)
- `uncertainty_novel_species_max` (float): Maximum uncertainty for novel species (default: 0.5)
- `uncertainty_novel_genus_max` (float): Maximum uncertainty for novel genus (default: 2.0)
- `uncertainty_conserved_min` (float): Minimum uncertainty for conserved regions (default: 5.0)

**Validation**:
Configuration is validated to ensure threshold relationships:
- novelty_novel_species_min > novelty_known_max
- novelty_novel_genus_min == novelty_novel_species_max (continuous boundaries)
- uncertainty_novel_species_max <= uncertainty_novel_genus_max

**Example**:
```python
from metadarkmatter.models.config import ScoringConfig

# Default configuration
config = ScoringConfig()

# Custom thresholds
config = ScoringConfig(
    bitscore_threshold_pct=95.0,
    novelty_novel_species_min=5.0,
    novelty_novel_species_max=15.0,
    uncertainty_novel_species_max=0.5,
)
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

### Bowtie2Config

Configuration for Bowtie2 competitive mapping.

```python
class Bowtie2Config(BaseModel):
    num_threads: int = 4
    mode: Literal["local", "end-to-end"] = "local"
    max_alignments: int = 100
```

**Parameters**:
- `num_threads` (int): Number of CPU threads (default: 4)
- `mode` (Literal["local", "end-to-end"]): Alignment mode (default: "local")
- `max_alignments` (int): Maximum alignments per read (-k parameter, default: 100)

---

### KrakenConfig

Configuration for Kraken2 taxonomic classification.

```python
class KrakenConfig(BaseModel):
    database: Path
    num_threads: int = 4
    confidence: float = 0.0
    minimum_base_quality: int = 0
```

**Parameters**:
- `database` (Path): Path to Kraken2 database (validated to exist)
- `num_threads` (int): Number of CPU threads (default: 4)
- `confidence` (float): Confidence score threshold (default: 0.0)
- `minimum_base_quality` (int): Minimum base quality for reads (default: 0)

---

### GlobalConfig

Global configuration for metadarkmatter pipeline.

```python
class GlobalConfig(BaseSettings):
    project_name: str = "metadarkmatter_analysis"
    output_dir: Path = Path("./results")
    scoring: ScoringConfig = Field(default_factory=ScoringConfig)
    blast: BlastConfig = Field(default_factory=BlastConfig)
    bowtie2: Bowtie2Config = Field(default_factory=Bowtie2Config)
    kraken: KrakenConfig | None = None
    num_threads: int = 4
    verbose: bool = False
```

**Description**:
Global configuration using pydantic-settings to load from YAML, environment variables, or CLI arguments.

**Environment Variables**:
- Prefix: `MDM_`
- Delimiter: `__` (for nested fields)

**Methods**:
```python
def to_yaml(self, path: Path) -> None
```
Write configuration to YAML file.

```python
@classmethod
def from_yaml(cls, path: Path) -> GlobalConfig
```
Load configuration from YAML file.

**Example**:
```python
from metadarkmatter.models.config import GlobalConfig

# Load from YAML
config = GlobalConfig.from_yaml(Path("config.yaml"))

# Save to YAML
config.to_yaml(Path("config_out.yaml"))

# Access nested config
print(f"BLAST evalue: {config.blast.evalue}")
print(f"Novel species min novelty: {config.scoring.novelty_novel_species_min}")
```

---

---

## Performance Characteristics

### Memory Usage Comparison

| Component | Dataset Size | Dense ANI | Sparse ANI | Notes |
|-----------|--------------|-----------|------------|-------|
| ANIMatrix | 1000 genomes | 4 MB | 4 MB | Minimal difference for small sets |
| ANIMatrix | 10000 genomes | 400 MB | ~20 MB (5% density) | Sparse recommended for large sets |
| BLAST Parsing | 100M records | 40+ GB | ~16 GB | Polars streaming keeps memory constant |
| Full Classification | 100M records | RAM overflow | ~32 GB | Streaming writes to disk in chunks |

### Classification Speed Comparison

| Method | Dataset | Speed | Notes |
|--------|---------|-------|-------|
| ANIWeightedClassifier.classify_to_dataframe() | 10M reads | ~5 min | Single-threaded, medium files |
| ANIWeightedClassifier.classify_to_dataframe_fast() | 10M reads | ~30 sec | Lightweight structures, ~10x faster |
| VectorizedClassifier.classify_file() | 10M reads | ~3 sec | Fully vectorized, ~100x faster |
| ParallelClassifier.classify_file() | 10M reads | ~1 min | 8 workers, with multiprocessing overhead |
| ParallelClassifier (16 workers) | 100M reads | ~2 hours | True horizontal scaling for massive files |

### Bitscore Threshold Effect

The `bitscore_threshold_pct` parameter significantly affects performance and results:

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
from metadarkmatter.core.ani_placement import ANIMatrix, ANIWeightedClassifier
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
    config = ScoringConfig(novelty_novel_species_min=3.0)  # Invalid: must be > 2.0
except ValueError as e:
    print(f"Invalid config: {e}")
    exit(1)

# Validate BLAST file
classifier = ANIWeightedClassifier(ani_matrix)
try:
    count = 0
    for result in classifier.classify_blast_file(Path("blast.tsv")):
        count += 1
except Exception as e:
    print(f"Error classifying reads: {e}")
    exit(1)

print(f"Successfully classified {count} reads")
```

---

## Quick Reference

### Most Common Workflows

#### 1. Classify BLAST file with default settings

```python
from pathlib import Path
from metadarkmatter.core.ani_placement import ANIMatrix, ANIWeightedClassifier

# Load ANI matrix
ani_matrix = ANIMatrix.from_file(Path("ani_matrix.csv"))

# Create classifier
classifier = ANIWeightedClassifier(ani_matrix)

# Classify and save results
num_reads = classifier.write_classifications(
    Path("blast.tsv"),
    Path("classifications.csv")
)
print(f"Classified {num_reads:,} reads")
```

#### 2. High-performance classification for 100M+ reads

```python
# Use vectorized classifier for maximum speed
vectorized = VectorizedClassifier(ani_matrix)

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
    uncertainty_novel_species_max=0.2,  # Require lower ambiguity
)

classifier = ANIWeightedClassifier(ani_matrix, config)
df = classifier.classify_to_dataframe_fast(Path("blast.tsv"))
```

---

## Further Reading

- See `IMPLEMENTATION_STRATEGY.md` for architectural details
- See `CLAUDE.md` for project background and methodology
- See workflow examples in `src/metadarkmatter/cli/` for CLI usage patterns
