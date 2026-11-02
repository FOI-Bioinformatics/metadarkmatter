# Metadarkmatter Technical Manual

**Version:** 0.2.0
**Last Updated:** 2026-11-02

> Scope note: several sections of this manual overlap with `docs/ARCHITECTURE.md`, `docs/API_REFERENCE.md` and `docs/METHODS.md`. Benchmark tables in sections 4 and 9 are historical figures from an earlier version and have not been re-measured for v0.2.0.

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [System Architecture](#2-system-architecture)
3. [Core Algorithm](#3-core-algorithm)
4. [Performance Optimization Strategy](#4-performance-optimization-strategy)
5. [Data Models and Validation](#5-data-models-and-validation)
6. [File Parsers and Streaming](#6-file-parsers-and-streaming)
7. [Classification Implementations](#7-classification-implementations)
8. [CLI Interface](#8-cli-interface)
9. [Performance Characteristics](#9-performance-characteristics)
10. [Testing and Quality Assurance](#10-testing-and-quality-assurance)
11. [Deployment Considerations](#11-deployment-considerations)
12. [Appendices](#12-appendices)

---

## 1. Executive Summary

### 1.1 Purpose

Metadarkmatter is a bioinformatics Python package designed to detect novel microbial diversity in environmental DNA (eDNA) samples through ANI-weighted placement uncertainty analysis. The software addresses the critical challenge of identifying "microbial dark matter"—unculturable bacterial taxa that cannot be detected through traditional culture-based methods or standard taxonomic classification approaches.

### 1.2 Scientific Context

Environmental samples from air filters, water sources, and other ecological matrices contain metagenomic sequencing data representing complex microbial communities. Standard reference-based classification methods fail to identify novel species and genera because these organisms lack representatives in existing databases. Metadarkmatter fills this gap by quantifying the uncertainty in taxonomic placement and measuring sequence divergence from known references.

### 1.3 Key Design Principles

The package is built on three foundational design principles:

1. **Performance-First Architecture**: Uses Polars DataFrame library for 10-100x faster processing compared to pandas, enabling analysis of files containing 100 million BLAST alignments. All CLI classification paths use the `VectorizedClassifier`, which executes entirely in Polars' Rust backend with no Python loops.
2. **Memory Efficiency**: Streaming algorithms with bounded memory usage, processing gigabyte-scale files within 8-16 GB RAM constraints typical of HPC nodes
3. **Type Safety**: Comprehensive Pydantic validation for all data models with frozen immutability guarantees

### 1.4 Target Users

- Computational biologists analyzing environmental metagenomes
- Microbial ecologists studying uncultured microbial diversity
- Bioinformatics pipeline developers integrating novel taxonomy detection
- HPC users processing large-scale metagenomic datasets

---

## 2. System Architecture

### 2.1 High-Level Component Diagram

```mermaid
graph TB
    subgraph "Input Layer"
        BLAST[BLAST Results<br/>tabular format]
        ANI[ANI Matrix<br/>CSV/TSV]
    end

    subgraph "Parsing Layer"
        BP[StreamingBlastParser<br/>Polars-based streaming]
        AP[ANIMatrixParser<br/>Validation & conversion]
    end

    subgraph "Core Algorithm"
        ANIM[ANIMatrix<br/>NumPy array O(1) lookup]
        ANIC[ANIWeightedClassifier<br/>Programmatic API, not used by CLI]
        VEC[VectorizedClassifier<br/>Polars-native vectorized]
    end

    subgraph "Data Models"
        BM[BlastHit/BlastResult<br/>Pydantic validation]
        CM[ReadClassification<br/>Result model]
        CFG[ScoringConfig<br/>Threshold configuration]
    end

    subgraph "Output Layer"
        CSV[CSV Files]
        PARQ[Parquet Files<br/>10x compression]
        JSON[JSON Summary<br/>Statistics]
    end

    subgraph "CLI Interface"
        MAIN[Typer CLI<br/>Rich output]
        SCORE[score classify<br/>Single file]
        BATCH[score batch<br/>Multiple files]
    end

    BLAST --> BP
    ANI --> AP

    BP --> BM
    BP --> VEC
    AP --> ANIM
    BM --> ANIC
    ANIM --> ANIC
    ANIM --> VEC
    CFG --> ANIC
    CFG --> VEC

    ANIC --> CM
    VEC --> DF[Polars DataFrame<br/>fixed output schema]

    DF --> CSV
    DF --> PARQ
    DF --> JSON

    MAIN --> SCORE
    MAIN --> BATCH
    SCORE --> VEC
    BATCH --> VEC
```

### 2.2 Package Structure

```
metadarkmatter/
├── src/metadarkmatter/
│   ├── __init__.py              # Package initialization
│   ├── core/                    # Core algorithm implementation
│   │   ├── constants.py         # Nucleotide thresholds
│   │   ├── protein_constants.py # Protein thresholds
│   │   ├── parsers.py           # Streaming parsers, extract_genome_name_expr()
│   │   ├── io_utils.py          # StreamingWriter (CSV append / Parquet row groups)
│   │   ├── metadata.py          # GenomeMetadata, representative mapping
│   │   ├── id_mapping.py        # ContigIdMapping
│   │   ├── genome_utils.py      # Accession extraction, genome concatenation
│   │   ├── exceptions.py        # MetadarkmatterError hierarchy
│   │   ├── ani_placement.py     # Backward-compatible re-exports
│   │   ├── ani_matrix_builder.py / aai_matrix_builder.py
│   │   ├── recruitment.py, random.py, runtime.py, sequence_validation.py,
│   │   │   logging_config.py
│   │   ├── classification/      # Classification pipeline
│   │   │   ├── thresholds.py    # Legacy threshold cascade
│   │   │   ├── qc.py            # Pre/post-classification QC
│   │   │   ├── sensitivity.py   # Threshold sensitivity analysis
│   │   │   ├── adaptive.py      # GMM-based adaptive thresholds
│   │   │   ├── bayesian.py      # Bayesian-primary classifier
│   │   │   ├── evaluation.py    # Accuracy / ECE / confusion matrix
│   │   │   ├── runner.py        # run_classification() used by score classify
│   │   │   ├── ani_matrix.py    # ANIMatrix class
│   │   │   └── classifiers/     # Classifier implementations
│   │   │       ├── base.py      # ANIWeightedClassifier (programmatic API)
│   │   │       └── vectorized.py# VectorizedClassifier (CLI classifier)
│   │   ├── novel_diversity/     # Novel-read clustering
│   │   ├── phylogeny/           # Phylogenetic tree building and placement
│   │   └── pipeline/            # orchestrator.py: end-to-end `mdm run`
│   ├── models/                  # Pydantic data models
│   │   ├── blast.py             # BLAST hit models
│   │   ├── classification.py    # Result models
│   │   ├── config.py            # ScoringConfig, BayesianConfig, BlastConfig
│   │   ├── genomes.py           # GenomeAccession, AccessionList
│   │   └── pipeline_config.py   # Aligner, PipelineStep
│   ├── cli/                     # Command-line interface (Typer)
│   │   ├── main.py              # CLI entry point
│   │   ├── score.py             # Scoring commands
│   │   └── ...                  # aai, ani, blast, blastx, mmseqs2, kraken2, map,
│   │                            # proteins, download, mapping, report, run, tree,
│   │                            # validate, visualize, doctor
│   ├── clients/                 # GTDB API client
│   ├── external/                # External tool wrappers
│   └── visualization/           # Plotly charts and HTML reports
│       ├── plots/               # Plot classes
│       └── report/              # ReportGenerator facade + sections/ mixins,
│                                # components/, assets/
├── tests/                       # Test suite
├── config/                      # Default configurations
└── docs/                        # Documentation
```

### 2.3 Data Flow Architecture

```mermaid
flowchart LR
    A[BLAST File<br/>100M+ alignments] --> B[StreamingBlastParser<br/>1M row chunks]
    B --> G[VectorizedClassifier<br/>Polars-native]

    H[ANI Matrix<br/>CSV] --> I[ANIMatrixParser<br/>Validation]
    I --> J[ANIMatrix<br/>NumPy array]

    J --> G

    G --> M[Polars DataFrame<br/>Direct]

    M --> O{Output Format}
    O -->|CSV| P[CSV Writer<br/>Streaming append]
    O -->|Parquet| Q[Parquet Writer<br/>Zstd compression]
```

### 2.4 Architectural Patterns

#### 2.4.1 Two Classifier Implementations

The package provides two classifier implementations for different use cases:

- **`ANIWeightedClassifier`**: Pydantic-validated, iterator-based programmatic API. Provides `classify_read()` for single-read classification in custom scripts and library usage.
- **`VectorizedClassifier`**: Polars-native operations in Rust backend. Sole classifier used by all CLI paths (batch, streaming). Executes entire classification pipeline with no Python loops.

The CLI uses `VectorizedClassifier` exclusively. Users can select streaming mode via `--streaming` for memory-bounded processing of very large files.

#### 2.4.2 Iterator Pattern: Memory-Bounded Streaming

All parsers implement streaming iterators to maintain bounded memory:

```python
# File size: 10 GB BLAST results
# Memory usage: ~1.8 GB (constant regardless of file size)
# VectorizedClassifier streaming mode processes in bounded partitions
vectorized = VectorizedClassifier(ani_matrix, config=config)
num_classified = vectorized.stream_to_file(
    blast_path=blast_path,
    output_path=output_path,
    output_format="parquet",
)
```

#### 2.4.3 Immutability Pattern: Frozen Pydantic Models

All data models use `frozen=True` to guarantee immutability:

```python
class BlastHit(BaseModel):
    qseqid: str
    pident: float
    # ... other fields
    model_config = {"frozen": True}
```

Benefits:
- Thread-safe for parallel processing
- Prevents accidental mutation bugs
- Enables efficient caching and memoization

#### 2.4.4 Lazy Evaluation: Polars LazyFrame

ANI matrix parsing and BLAST file reading use lazy evaluation:

```python
# No data loaded until .collect() called
lazy_df = parser.parse_lazy()
result = (
    lazy_df
    .filter(pl.col("bitscore") > 100)
    .group_by("qseqid")
    .first()
    .collect(engine="streaming")  # Execute in streaming mode
)
```

---

## 3. Core Algorithm

### 3.1 ANI-Weighted Placement Uncertainty

The core algorithm quantifies two orthogonal metrics to classify metagenomic reads:

#### 3.1.1 Novelty Index (N)

**Formula:**
```
N = 100 - TopHitIdentity
```

**Interpretation:**
- `N < 4`: High identity match to reference (known species)
- `4 <= N < 20`: Moderate divergence (potential novel species)
- `20 <= N <= 25`: High divergence (potential novel genus)
- `N > 25`: Very distant homology (beyond genus-level classification)

**Biological Context:**
Within bacterial taxonomy, species are typically defined by ANI >= 96% (N <= 4%). Novel species occupy the 80-96% identity range, while novel genera show 75-80% identity to closest known representatives.

#### 3.1.2 Placement Uncertainty (U)

**Formula:**
```
U = 100 - ANI(TopHitGenome, SecondBestGenome)          # uncertainty_mode="second" (default)
U = 100 - max(ANI(TopHitGenome, SecondaryHitGenome))   # uncertainty_mode="max"
```

Where secondary hits are all BLAST hits within `bitscore_threshold_pct` (default 95%) of the top bitscore, and the second-best genome is the highest-ranked competing genome using the same ordering as best-hit selection.

**Interpretation:**
- `U < 1.5`: Confident placement (competing genomes share >98.5% ANI)
- `1.5 <= U < 5`: Species boundary zone (competing genomes share 95-98.5% ANI)
- `U >= 5`: High ambiguity, conserved region (competing genomes share <95% ANI)

**Biological Context:**
High placement uncertainty indicates conserved genomic regions (e.g., ribosomal RNA, housekeeping genes) that provide insufficient phylogenetic resolution. These reads are excluded from novelty detection to avoid false positives.

### 3.2 Classification Decision Tree

The primary output `taxonomic_call` is the Bayesian MAP category: posteriors over six categories (Known Species, Novel Species, Novel Genus, Species Boundary, Ambiguous, Unclassified) from 2D Gaussians in (N, U) space, followed by a Stage 2 refinement that can assign Conserved Region or Ambiguous Within Genus (`core/classification/bayesian.py`). Its decision boundaries do not coincide exactly with the hard thresholds. The hard-threshold cascade is retained as `legacy_call` (`core/classification/thresholds.py`); the tree below is a simplified view of that cascade. The full legacy order is: identity gap < 2 with >1 hit -> Ambiguous; U >= 5 across genera -> Conserved Region; U >= 5 -> Ambiguous; U >= 1.5 -> Species Boundary; Known Species; single-hit read with high inferred uncertainty -> Ambiguous; Novel Species; Novel Genus; otherwise Unclassified.

```mermaid
graph TD
    START[Read with BLAST hits] --> CHECK_U{U >= 5.0?}

    CHECK_U -->|Yes| CONSERVED[Conserved Region if hits span genera,<br/>otherwise Ambiguous]

    CHECK_U -->|No| CHECK_SB{U >= 1.5?}
    CHECK_SB -->|Yes| SB[Species Boundary]
    CHECK_SB -->|No| CHECK_KNOWN{N < 4.0 AND<br/>U < 1.5?}

    CHECK_KNOWN -->|Yes| KNOWN[Known Species<br/>High identity, unambiguous]

    CHECK_KNOWN -->|No| CHECK_NOVEL_SP{4.0 <= N < 20.0<br/>AND U < 1.5?}

    CHECK_NOVEL_SP -->|Yes| NOVEL_SP[Novel Species<br/>Moderate divergence]

    CHECK_NOVEL_SP -->|No| CHECK_NOVEL_GEN{20.0 <= N <= 25.0<br/>AND U < 1.5?}

    CHECK_NOVEL_GEN -->|Yes| NOVEL_GEN[Novel Genus<br/>High divergence]

    CHECK_NOVEL_GEN -->|No| DEFAULT[Unclassified<br/>Default fallback]

    style NOVEL_SP fill:#ffeb3b
    style NOVEL_GEN fill:#ff9800
    style KNOWN fill:#4caf50
    style CONSERVED fill:#9e9e9e
    style DEFAULT fill:#9e9e9e
```

### 3.3 Algorithmic Steps

#### Step 1: Top Hit Identification

```python
# Input: BlastResult with sorted hits (highest bitscore first)
best_hit = blast_result.hits[0]
top_hit_identity = best_hit.pident  # Percent identity (0-100)
best_genome = best_hit.genome_name
```

**Complexity:** O(1) - hits pre-sorted by parser

#### Step 2: Novelty Index Calculation

```python
novelty_index = 100.0 - top_hit_identity
```

**Complexity:** O(1)

#### Step 3: Ambiguous Hit Detection

```python
bitscore_cutoff = best_hit.bitscore * 0.95  # 95% threshold
ambiguous_hits = []

for hit in blast_result.hits:
    if hit.bitscore < bitscore_cutoff:
        break  # Early termination - hits are sorted
    if hit.genome_name != best_genome:
        ambiguous_hits.append(hit)
```

**Complexity:** O(k) where k = number of ambiguous hits (typically k < 50)
**Optimization:** Early termination exploits sorted order

#### Step 4: Placement Uncertainty via ANI Lookup

```python
max_ani = 0.0
for secondary_hit in ambiguous_hits:
    # O(1) ANI lookup using NumPy array indexing
    ani = ani_matrix.get_ani(best_genome, secondary_hit.genome_name)
    if ani > max_ani:
        max_ani = ani

placement_uncertainty = 100.0 - max_ani if max_ani > 0 else 100.0
```

This shows `uncertainty_mode="max"`. With the default `uncertainty_mode="second"`, the ANI to the first competing genome (second-best by score) is used instead of the maximum.

**Complexity:** O(k) where k = number of ambiguous hits
**Per-lookup cost:** O(1) via NumPy array (vs O(log n) for nested dict)

#### Step 5: Threshold-Based Classification

The snippet below is a simplified form of the legacy cascade that produces `legacy_call`; it omits the identity-gap, Species Boundary, single-hit and Ambiguous rules listed in section 3.2, and its fallback is `Unclassified`, not Conserved Region. `taxonomic_call` comes from the Bayesian classifier.

```python
# High uncertainty → Conserved region
if placement_uncertainty >= 5.0:
    return TaxonomicCall.CONSERVED_REGION

# Low novelty + low uncertainty → Known species
if novelty_index < 4.0 and placement_uncertainty < 1.5:
    return TaxonomicCall.KNOWN_SPECIES

# Moderate novelty + low uncertainty → Novel species
if 4.0 <= novelty_index < 20.0 and placement_uncertainty < 1.5:
    return TaxonomicCall.NOVEL_SPECIES

# High novelty + low uncertainty → Novel genus
if 20.0 <= novelty_index <= 25.0 and placement_uncertainty < 1.5:
    return TaxonomicCall.NOVEL_GENUS

# Default fallback
return TaxonomicCall.CONSERVED_REGION
```

**Complexity:** O(1) - simple threshold comparisons

### 3.4 Overall Complexity Analysis

**Per-read classification:**
- Best case: O(1) - single unambiguous hit
- Average case: O(k) where k ≈ 5-20 ambiguous hits
- Worst case: O(h) where h = total hits per read (bounded by BLAST `max_target_seqs`)

**Full file classification (n reads):**
- ANIWeightedClassifier (per-read API): O(n * k)
- VectorizedClassifier (CLI): O(n) amortized (Polars parallelizes operations internally)

**Memory:**
- VectorizedClassifier batch: O(n) - loads full dataset for parallel aggregation
- VectorizedClassifier streaming: O(chunk_size) - bounded by partition size

---

## 4. Performance Optimization Strategy

### 4.1 Optimization Strategy

The package implements a unified high-performance approach using Polars' Rust backend:

```
Foundation: Polars DataFrame library           → 10-100x faster than pandas
    +
Vectorized Classification                     → Zero Python loops, SIMD operations
    +
Streaming Mode (optional)                     → Bounded memory for very large files
```

All CLI classification paths use `VectorizedClassifier`, which executes the entire pipeline in compiled Rust via Polars expressions. For programmatic single-read classification, `ANIWeightedClassifier.classify_read()` provides a Pydantic-validated API.

### 4.2 Critical Optimization: ANI Matrix as NumPy Array

#### Problem

Original nested dictionary approach:
```python
# Slow: O(hash) lookups twice per ANI query
ani_dict = {
    "GCF_000001.1": {"GCF_000002.1": 95.3, ...},
    "GCF_000002.1": {"GCF_000001.1": 95.3, ...},
    ...
}
ani = ani_dict[genome1][genome2]  # Two hash lookups
```

**Issues:**
- Memory: 150 MB for 1000 genomes (due to Python object overhead)
- Speed: O(hash) × 2 per lookup
- Cache misses: Poor spatial locality for nested dicts

#### Solution

NumPy array with integer indexing:

```python
class ANIMatrix:
    def __init__(self, ani_dict: dict[str, dict[str, float]]):
        # Map genome names to integer indices
        self._genomes = tuple(sorted(ani_dict.keys()))
        self._genome_to_idx = {g: i for i, g in enumerate(self._genomes)}

        # Create dense float64 array
        n = len(self._genomes)
        self._ani_array = np.zeros((n, n), dtype=np.float64)

        # Populate array
        for g1, inner in ani_dict.items():
            i = self._genome_to_idx[g1]
            for g2, ani_val in inner.items():
                j = self._genome_to_idx[g2]
                self._ani_array[i, j] = ani_val

    def get_ani(self, genome1: str, genome2: str) -> float:
        # Fast path for self-ANI
        if genome1 == genome2:
            return 100.0

        # Two dict lookups (O(hash) × 2)
        i = self._genome_to_idx.get(genome1)
        j = self._genome_to_idx.get(genome2)

        # O(1) array indexing
        return float(self._ani_array[i, j])
```

Simplified. The actual `ANIMatrix` (`core/classification/ani_matrix.py`) also: returns 0.0 when a genome is absent; returns `default_ani` (70.0) for pairs stored as 0.0; mirrors one-directional pairs; runs an optional symmetry check (`symmetry_check="warn"` by default); and rejects matrices on a 0-1 scale with `ValueError`.

**Performance Gains:**

| Metric | Nested Dict | NumPy Array | Improvement |
|--------|-------------|-------------|-------------|
| Memory (1000 genomes) | 150 MB | 8 MB | 18.7x reduction |
| Lookup time | 120 ns | 7.5 ns | 16x faster |
| Cache efficiency | Poor | Excellent | Spatial locality |

**Benchmarked on:** 1000 genomes, 10M lookups

### 4.3 Critical Optimization: Polars vs Pandas

#### BLAST File Parsing Comparison

**Pandas approach (slow):**
```python
# Loads entire 10 GB file into RAM
df = pd.read_csv("blast.tsv", sep="\t", header=None)

# Per-row genome extraction (100M string operations in Python)
df["genome_name"] = df["sseqid"].str.split("_").str[0]

# Group by read (slow due to object dtype)
grouped = df.groupby("qseqid")
```

**Polars approach (fast):**
```python
# Lazy evaluation - no data loaded yet
lazy_df = pl.scan_csv("blast.tsv", separator="\t", has_header=False)

# Vectorized regex extraction in Rust (accession before the first '|')
lazy_df = lazy_df.with_columns([
    pl.col("sseqid")
      .str.extract(r"^([^|]+)\|", 1)
      .fill_null("unknown")
      .alias("genome_name")
])

# Streaming execution with bounded memory
result = lazy_df.collect(engine="streaming")
```

**Performance Comparison (100M BLAST hits):**

| Operation | Pandas | Polars | Speedup |
|-----------|--------|--------|---------|
| File read | 180s | 12s | 15x |
| Genome extraction | 320s | 3s | 106x |
| Groupby aggregation | 95s | 8s | 11.8x |
| Memory usage | 42 GB | 6 GB | 7x reduction |
| **Total** | **595s** | **23s** | **25.8x** |

### 4.4 Vectorized Classification Architecture

The `VectorizedClassifier` eliminates Python loops entirely by using Polars' Rust backend:

> **Note:** The code below is a simplified sketch of the approach, not the current implementation. See the referenced source file for the actual code.

```python
def classify_file(self, blast_path: Path) -> pl.DataFrame:
    # Step 1: Load all data with genome extraction (vectorized)
    df = (
        pl.scan_csv(blast_path, separator="\t")
        .with_columns([extract_genome_name_expr()])  # Regex in Rust
        .collect()
    )

    # Step 2: Find top hit per read (parallelized group_by)
    best_hits = df.group_by("qseqid").agg([
        pl.col("bitscore").max().alias("max_bitscore"),
        pl.col("pident")
          .filter(pl.col("bitscore") == pl.col("bitscore").max())
          .first()
          .alias("top_pident"),
        pl.col("genome_name")
          .filter(pl.col("bitscore") == pl.col("bitscore").max())
          .first()
          .alias("best_genome"),
    ])

    # Step 3: Calculate bitscore threshold (vectorized)
    best_hits = best_hits.with_columns([
        (pl.col("max_bitscore") * 0.95).alias("bitscore_threshold")
    ])

    # Step 4: Join to find ambiguous hits (hash join)
    ambiguous = df.join(best_hits, on="qseqid").filter(
        pl.col("bitscore") >= pl.col("bitscore_threshold")
    )

    # Step 5: ANI lookup via join (not Python loop)
    secondary_with_ani = ambiguous.join(
        self._ani_lookup_symmetric,  # Pre-built DataFrame
        left_on=["best_genome", "genome_name"],
        right_on=["genome1", "genome2"],
        how="left"
    )

    # Step 6: Max ANI per read (parallelized aggregation)
    max_ani_per_read = (
        secondary_with_ani
        .group_by("qseqid")
        .agg(pl.col("ani").max().alias("max_ani"))
    )

    # ... classification via when/then expressions (vectorized)
```

**Key Advantages:**

1. **Automatic Parallelization:** Polars uses Rust threads internally
2. **SIMD Vectorization:** CPU vector instructions for numeric operations
3. **Zero Python Loops:** Entire classification executes in compiled Rust
4. **Memory Efficiency:** Lazy evaluation with streaming support

**Benchmark (10M reads, 8 cores):**

| Classifier | Time | Speedup |
|------------|------|---------|
| ANIWeightedClassifier (per-read Python loop) | 245s | 1x |
| VectorizedClassifier (Polars) | 14s | 17.5x |

---

## 5. Data Models and Validation

### 5.1 Model Hierarchy

```mermaid
classDiagram
    class BaseModel {
        <<Pydantic>>
        +model_validate()
        +model_dump()
    }

    class BlastHit {
        +str qseqid
        +str sseqid
        +float pident
        +float bitscore
        +genome_name() str
        +from_blast_line() BlastHit
    }

    class BlastResult {
        +str read_id
        +tuple hits
        +best_hit() BlastHit
        +iter_ambiguous_hits() Iterator
    }

    class ReadClassification {
        +str read_id
        +str best_match_genome
        +float novelty_index
        +float placement_uncertainty
        +TaxonomicCall taxonomic_call
        +is_novel() bool
    }

    class TaxonomicCall {
        <<StrEnum>>
        KNOWN_SPECIES
        NOVEL_SPECIES
        NOVEL_GENUS
        SPECIES_BOUNDARY
        AMBIGUOUS
        AMBIGUOUS_WITHIN_GENUS
        CONSERVED_REGION
        UNCLASSIFIED
        OFF_TARGET
    }

    class ScoringConfig {
        +float bitscore_threshold_pct
        +float novelty_known_max
        +float uncertainty_conserved_min
        +validate_threshold_boundaries()
    }

    class TaxonomicSummary {
        +int total_reads
        +int novel_species
        +int novel_genus
        +dict genome_hit_counts
        +novel_diversity_pct() float
        +to_json() void
    }

    BaseModel <|-- BlastHit
    BaseModel <|-- BlastResult
    BaseModel <|-- ReadClassification
    BaseModel <|-- ScoringConfig
    BaseModel <|-- TaxonomicSummary

    BlastResult o-- BlastHit : contains
    ReadClassification --> TaxonomicCall : uses
```

### 5.2 BlastHit Model

**File:** `src/metadarkmatter/models/blast.py`

**Purpose:** Represents a single BLAST alignment with comprehensive validation.

**Key Fields:**
```python
class BlastHit(BaseModel):
    qseqid: str              # Query ID (read identifier)
    sseqid: str              # Subject ID (genome contig)
    pident: float            # Percent identity [0, 100]
    bitscore: float          # BLAST bitscore (≥ 0)
    evalue: float            # E-value (≥ 0)
    qlen: int | None         # Query length (13-column input), optional
    # ... 7 additional fields for alignment coordinates
```

**Critical Validation:**

1. **Percent Identity Clamping:**
   ```python
   @field_validator("pident", mode="before")
   @classmethod
   def clamp_pident(cls, v: float) -> float:
       """BLAST occasionally reports pident > 100 due to rounding."""
       return max(0.0, min(100.0, float(v)))
   ```

2. **Alignment Coordinate Normalization:** reverse-strand hits reported with qstart > qend (MMseqs2, DIAMOND) are accepted and the two coordinates are swapped, so coverage is computed on the correct span.
   ```python
   @model_validator(mode="after")
   def validate_query_positions(self) -> Self:
       if self.qend < self.qstart:
           qstart, qend = self.qend, self.qstart
           object.__setattr__(self, "qstart", qstart)
           object.__setattr__(self, "qend", qend)
       return self
   ```

3. **Genome Name Extraction:**
   ```python
   GENOME_PATTERN = re.compile(
       r"^(?P<genome>GCF_\d+\.\d+|GCA_\d+\.\d+|NZ_[A-Z]+\d+\.\d+)"
       r"|(?P<alt>[^|]+\|(?P<name>[^|]+))"
   )

   @property
   def genome_name(self) -> str:
       if not self.sseqid or not self.sseqid.strip():
           return "unknown"
       if "|" in self.sseqid:              # "{genome}|{contig}"
           return self.sseqid.split("|", 1)[0].strip() or "unknown"
       match = self.GENOME_PATTERN.match(self.sseqid)
       if match:
           if match.group("genome"):
               return match.group("genome")
           if match.group("name"):
               return match.group("name")
       return self.sseqid.split()[0]
   ```

   Note: the vectorized path (`extract_genome_name_expr()`, section 6.1.2) takes the text *before* the first `|`, which matches the standardized `{accession}|{contig_id}` headers. For non-accession `a|b` IDs the two paths differ.

**Design Rationale:**

- **Frozen model** (`frozen=True`): Prevents accidental mutation, enables safe sharing across threads
- **Lazy genome extraction**: `@property` ensures extraction only when needed (not during parsing)
- **Defensive validation**: Handles real-world BLAST edge cases (e.g., pident=100.001 from floating-point arithmetic)

### 5.3 BlastResult Model

**Purpose:** Groups all BLAST hits for a single read with sorted access.

```python
class BlastResult(BaseModel):
    read_id: str
    hits: tuple[BlastHit, ...]  # Immutable tuple, sorted by bitscore

    @model_validator(mode="after")
    def ensure_sorted(self) -> Self:
        """Validate or enforce bitscore descending order."""
        if len(self.hits) <= 1:
            return self

        # Check if already sorted (O(n) vs O(n log n) sort)
        is_sorted = all(
            self.hits[i].bitscore >= self.hits[i+1].bitscore
            for i in range(len(self.hits) - 1)
        )

        if not is_sorted:
            sorted_hits = tuple(sorted(self.hits, key=lambda h: h.bitscore, reverse=True))
            object.__setattr__(self, "hits", sorted_hits)

        return self
```

**Optimization:** The parser pre-sorts hits in Polars, so the validation typically performs an O(n) check rather than O(n log n) sort.

**Ambiguous Hit Iterator:**
```python
def iter_ambiguous_hits(self, threshold_pct: float = 95.0) -> Iterator[BlastHit]:
    """Memory-efficient iteration with early termination."""
    if not self.hits:
        return

    cutoff = self.top_bitscore * (threshold_pct / 100.0)

    for hit in self.hits:
        if hit.bitscore < cutoff:
            break  # Sorted order allows early stop
        yield hit
```

**Use Case:** Iterate over competitive hits without materializing a list (saves memory for reads with 100+ hits).

### 5.4 ReadClassification Model

**Purpose:** Final classification output with computed properties.

```python
class ReadClassification(BaseModel):
    read_id: str
    best_match_genome: str
    top_hit_identity: float       # 0-100
    novelty_index: float          # 0-100
    placement_uncertainty: float  # 0-100
    num_ambiguous_hits: int
    taxonomic_call: TaxonomicCall

    # ... plus optional enhanced-scoring fields (second_hit_identity,
    # identity_gap, genus_uncertainty, confidence_score, inferred_uncertainty, ...)

    model_config = {"frozen": True}

    # diversity_status (Known/Novel/Uncertain) is also a computed field
    @computed_field
    @property
    def is_novel(self) -> bool:
        """Derived field: True if Novel Species or Novel Genus."""
        return self.taxonomic_call in (
            TaxonomicCall.NOVEL_SPECIES,
            TaxonomicCall.NOVEL_GENUS
        )
```

**Pydantic `@computed_field`:**
This decorator ensures `is_novel` is included in serialization (`.model_dump()`, `.model_dump_json()`) and DataFrame conversion, but computed on-the-fly rather than stored.

### 5.5 ScoringConfig Model

**Purpose:** Centralized threshold configuration with validation.

```python
class ScoringConfig(BaseModel):
    bitscore_threshold_pct: float = Field(default=95.0, ge=0, le=100)

    # Novelty thresholds (N = 100 - TopHitIdentity)
    novelty_known_max: float = 4.0          # N < 4% = known species (pident > 96%)
    novelty_novel_species_min: float = 4.0   # Continuous boundary with known_max
    novelty_novel_species_max: float = 20.0
    novelty_novel_genus_min: float = 20.0
    novelty_novel_genus_max: float = 25.0

    # Uncertainty thresholds (U = 100 - ANI to second-best genome by default;
    # uncertainty_mode="max" uses the max ANI between competing genomes)
    uncertainty_known_max: float = 1.5       # U < 1.5% = confident placement
    uncertainty_novel_species_max: float = 1.5
    uncertainty_novel_genus_max: float = 1.5
    uncertainty_conserved_min: float = 5.0

    @model_validator(mode="after")
    def validate_threshold_boundaries(self) -> Self:
        """Ensure thresholds form valid hierarchy."""
        if self.novelty_novel_species_min < self.novelty_known_max:
            raise ValueError("Novel species min must be >= known max")
        # Also enforced: novelty_novel_genus_min == novelty_novel_species_max,
        # uncertainty_novel_species_max <= uncertainty_novel_genus_max, and
        # coverage_weight_strength > 0 when coverage weighting is enabled.

        return self

    model_config = {"frozen": True}
```

Excerpt only; `ScoringConfig` (`models/config.py`) has further fields for alignment filters, coverage weighting, `uncertainty_mode`, family validation, AAI, `bayesian` (BayesianConfig) and YAML I/O (`from_yaml`, `to_yaml`).

**Design Decision:** Frozen configuration prevents runtime modification, ensuring reproducibility. Users must create new config instances for different thresholds.

### 5.6 Protein Mode and Alignment-Specific Thresholds

ScoringConfig supports an `alignment_mode` parameter to select between nucleotide and protein thresholds:

```python
class ScoringConfig(BaseModel):
    alignment_mode: Literal["nucleotide", "protein"] = "nucleotide"

    def get_effective_thresholds(self) -> dict[str, float]:
        """Return mode-appropriate thresholds."""
        if self.alignment_mode == "protein":
            return {
                "novelty_known_max": 10.0,
                "novelty_novel_species_min": 10.0,
                "novelty_novel_species_max": 25.0,
                "novelty_novel_genus_min": 25.0,
                "novelty_novel_genus_max": 40.0,
                "uncertainty_known_max": 5.0,
                "uncertainty_conserved_min": 10.0,
            }
        # Return nucleotide thresholds (config defaults)
        return {
            "novelty_known_max": self.novelty_known_max,
            ...
        }
```

In protein mode the returned thresholds are fixed constants; custom novelty/uncertainty values set on `ScoringConfig` (directly, via `--config` YAML, or via `--adaptive-thresholds`) are not applied. The actual return type is the `EffectiveThresholds` TypedDict, which also carries confidence-score scaling parameters.

**Rationale:** Protein sequences are more conserved than nucleotide sequences due to the degeneracy of the genetic code. At the genus level, nucleotide identity degrades faster than protein identity, making BLASTX more sensitive for detecting divergent homology.

**Threshold comparison:**

| Category | Nucleotide | Protein | Biological Basis |
|----------|-----------|---------|------------------|
| Known Species | N < 4%, U < 1.5% | N < 10%, U < 5% | >96% ANI / >90% AAI |
| Novel Species | 4% <= N < 20%, U < 1.5% | 10% <= N < 25%, U < 5% | 80-96% ANI / 75-90% AAI |
| Novel Genus | 20% <= N <= 25%, U < 1.5% | 25% <= N <= 40%, U < 5% | 75-80% ANI / 60-75% AAI |

**Protein constants** are defined in `core/protein_constants.py`:

```python
# Protein novelty thresholds (wider than nucleotide)
PROTEIN_NOVELTY_KNOWN_MAX = 10.0
PROTEIN_NOVELTY_NOVEL_SPECIES_MIN = 10.0
PROTEIN_NOVELTY_NOVEL_SPECIES_MAX = 25.0
PROTEIN_NOVELTY_NOVEL_GENUS_MIN = 25.0
PROTEIN_NOVELTY_NOVEL_GENUS_MAX = 40.0

# Protein uncertainty thresholds (wider than nucleotide)
PROTEIN_UNCERTAINTY_CONFIDENT_MAX = 5.0
PROTEIN_UNCERTAINTY_CONSERVED_MIN = 10.0
```

---

## 6. File Parsers and Streaming

### 6.1 StreamingBlastParser Architecture

**File:** `src/metadarkmatter/core/parsers.py` (`StreamingBlastParser`)

**Challenge:** BLAST output files can contain 100 million alignments (10+ GB), exceeding typical RAM constraints.

**Solution:** Chunk-based streaming with read-boundary detection.

#### 6.1.1 Chunk Processing with Boundary Handling

> **Note:** The code below is a simplified sketch of the approach, not the current implementation. See the referenced source file for the actual code.

```python
class StreamingBlastParser:
    def __init__(self, blast_path: Path, chunk_size: int = 1_000_000):
        self.blast_path = blast_path
        self.chunk_size = chunk_size  # Rows per chunk, not bytes

    def iter_reads_fast(self) -> Iterator[BlastResultFast]:
        """Stream reads with correct boundary handling."""
        pending_hits: dict[str, list[BlastHitFast]] = {}

        # parse_lazy() applies an explicit schema and '#' comment handling
        batches = self.parse_lazy().collect_batches(chunk_size=self.chunk_size)

        for chunk_df in batches:
            if chunk_df.is_empty():
                continue

            # Vectorized genome extraction + sorting (Polars Rust backend)
            chunk_df = chunk_df.with_columns([
                extract_genome_name_expr()
            ]).sort(["qseqid", "bitscore"], descending=[False, True])

            # Accumulate hits in pending buffer
            for row in chunk_df.iter_rows():
                qseqid, sseqid, pident, bitscore, genome_name = row
                hit = BlastHitFast(qseqid, sseqid, pident, bitscore, genome_name)

                if qseqid not in pending_hits:
                    pending_hits[qseqid] = []
                pending_hits[qseqid].append(hit)

            # Yield complete reads (not the last read, which may span chunks)
            current_last_read = chunk_df["qseqid"][-1]
            for read_id in list(pending_hits.keys()):
                if read_id != current_last_read:
                    hits = pending_hits.pop(read_id)
                    yield BlastResultFast(read_id=read_id, hits=tuple(hits))

        # Yield remaining reads
        for read_id, hits in pending_hits.items():
            yield BlastResultFast(read_id=read_id, hits=tuple(hits))
```

**Key Design Elements:**

1. **Polars Batched Reader:** Reads file in configurable chunks (default 1M rows), using an explicit schema for 12- or 13-column (with `qlen`) input
2. **Pending Buffer:** Accumulates hits for reads spanning chunk boundaries
3. **Last-Read Exclusion:** Holds back the last read in each chunk until next batch confirms completion
4. **Vectorized Operations:** Genome extraction and sorting happen in Polars (Rust) before Python iteration

**Memory Footprint:**

- Chunk size: 1M rows × ~200 bytes/row = ~200 MB per chunk
- Pending buffer: Typically < 100 reads × ~50 hits × ~100 bytes = ~500 KB
- Total: ~200 MB regardless of total file size

#### 6.1.2 Vectorized Genome Name Extraction

**Problem:** Extracting genome name from `sseqid` for 100M rows in Python is prohibitively slow.

**Solution:** Polars expression that executes in Rust:

```python
def extract_genome_name_expr() -> pl.Expr:
    """Extract the accession from standardized '{accession}|{contig_id}' headers."""
    return (
        pl.col("sseqid")
        .str.extract(r"^([^|]+)\|", 1)   # Regex in Rust: text before first '|'
        .fill_null("unknown")            # No pipe -> 'unknown'
        .alias("genome_name")
    )
```

Subject IDs without a `|` map to `unknown`; external alignments with other ID formats need `--id-mapping` (see `ContigIdMapping`, which also handles `ref|ACC|` IDs).

**Benchmark (100M rows):**

| Method | Time |
|--------|------|
| Python loop: `sseqid.split("_")[0]` | 320s |
| Pandas: `df["sseqid"].str.split("_").str[0]` | 45s |
| Polars vectorized: `extract_genome_name_expr()` | 3s |

**Speedup:** 106x faster than Python, 15x faster than pandas

### 6.2 ANIMatrixParser

**File:** `src/metadarkmatter/core/parsers.py` (`ANIMatrixParser`)

**Purpose:** Parse and validate ANI matrices with comprehensive error checking.

#### 6.2.1 Format Specification

**Expected CSV/TSV structure:**
```csv
genome,GCF_000001.1,GCF_000002.1,GCF_000003.1
GCF_000001.1,100.0,95.3,87.2
GCF_000002.1,95.3,100.0,88.5
GCF_000003.1,87.2,88.5,100.0
```

**Requirements:**
- Square matrix: `n_rows = n_columns - 1` (first column is genome names)
- Symmetric: `ANI[i,j] = ANI[j,i]` (not checked by the parser; `ANIMatrix` warns on pairs differing by more than 0.5 by default, raises with `symmetry_check="strict"`, and mirrors pairs given in one direction only)
- Diagonal: `ANI[i,i] = 100.0`
- Range: `0 ≤ ANI ≤ 100`, percent scale (a 0-1 scaled matrix is rejected by `ANIMatrix`)
- No null values; 0 marks an uncomputed pair and is looked up as the default unrelated ANI (70)

#### 6.2.2 Validation Pipeline

> **Note:** The code below is a simplified sketch of the approach, not the current implementation. See the referenced source file for the actual code.

```python
class ANIMatrixParser:
    def parse(self) -> pl.DataFrame:
        """Parse and validate ANI matrix."""
        df = pl.read_csv(self.ani_path, separator=self._detect_sep())
        self._validate_matrix(df)
        return df

    def _validate_matrix(self, df: pl.DataFrame) -> None:
        """Comprehensive validation."""
        # 1. Check square dimensions
        genome_col = df.columns[0]
        num_genomes = len(df)
        num_cols = len(df.columns) - 1
        if num_genomes != num_cols:
            raise ValueError(f"Non-square matrix: {num_genomes} rows, {num_cols} cols")

        # 2. Validate row/column name matching
        row_genomes = set(df[genome_col].to_list())
        col_genomes = set(df.columns[1:])
        if row_genomes != col_genomes:
            raise ValueError("Row genome names don't match column headers")

        # 3. Check for duplicates
        if len(row_genomes) != num_genomes:
            duplicates = [g for g in row_genomes if df[genome_col].to_list().count(g) > 1]
            raise ValueError(f"Duplicate genomes: {duplicates}")

        # 4. Validate ANI range [0, 100]
        numeric_cols = df.select(pl.exclude(genome_col))
        if numeric_cols.null_count().sum_horizontal()[0] > 0:
            raise ValueError("Matrix contains null values")

        min_val = numeric_cols.min().min_horizontal()[0]
        max_val = numeric_cols.max().max_horizontal()[0]
        if min_val < 0 or max_val > 100:
            raise ValueError(f"ANI values out of range: [{min_val}, {max_val}]")
```

**Design Rationale:**

- **Fail-fast validation:** Catches malformed inputs before expensive processing
- **Informative errors:** Specific messages aid debugging (e.g., lists duplicate genomes)
- **Polars operations:** Uses vectorized min/max for fast range validation

#### 6.2.3 Conversion to NumPy Array

```python
def to_dict(self) -> dict[str, dict[str, float]]:
    """Convert to nested dict for ANIMatrix initialization."""
    df = self.parse()
    genome_col = df.columns[0]

    ani_dict = {}
    for row in df.iter_rows(named=True):
        genome_name = row[genome_col]
        ani_dict[genome_name] = {
            col: val
            for col, val in row.items()
            if col != genome_col
        }

    return ani_dict
```

**Usage:**
```python
parser = ANIMatrixParser(Path("ani_matrix.csv"))
ani_dict = parser.to_dict()
ani_matrix = ANIMatrix(ani_dict)  # Converts to NumPy
```

---

## 7. Classification Implementations

### 7.1 Classifier Overview

| Classifier | Speed | Memory | Parallelism | Validation | Use Case |
|------------|-------|--------|-------------|------------|----------|
| **ANIWeightedClassifier** | 1x | Low | Single-thread | Full Pydantic | Programmatic API (`classify_read()`) |
| **VectorizedClassifier** | 17x | Medium-High | Auto-parallel (Polars) | None | All CLI paths (batch + streaming) |

### 7.2 ANIWeightedClassifier (Programmatic API)

**File:** `src/metadarkmatter/core/classification/classifiers/base.py`

**Architecture:** Iterator-based with full Pydantic validation. Provides `classify_read()` for single-read classification in custom scripts and library usage. Not used by the CLI.

#### 7.2.1 Core Classification Method

> **Note:** The code below is a simplified sketch of the approach, not the current implementation. See the referenced source file for the actual code.

The actual method also selects the best hit by coverage-weighted bitscore, computes genus uncertainty, inferred uncertainty and confidence, and sets `taxonomic_call` from the Bayesian posteriors (Stage 1 + Stage 2) rather than from `_classify_by_thresholds()` alone.

```python
class ANIWeightedClassifier:
    def __init__(
        self,
        ani_matrix: ANIMatrix,
        aai_matrix: AAIMatrix | None = None,
        config: ScoringConfig | None = None,
    ):
        self.ani_matrix = ani_matrix
        self.aai_matrix = aai_matrix
        self.config = config or ScoringConfig()

    def classify_read(
        self, blast_result: BlastResult, read_length: int | None = None
    ) -> ReadClassification | None:
        """Classify single read with full validation."""
        if not blast_result.hits:
            return None

        # Step 1: Get top hit (already validated by BlastResult model)
        best_hit = blast_result.best_hit
        top_hit_identity = best_hit.pident
        novelty_index = 100.0 - top_hit_identity
        best_genome = best_hit.genome_name

        # Step 2: Calculate placement uncertainty
        placement_uncertainty, num_ambiguous_hits = self._calculate_placement_uncertainty(
            best_genome,
            blast_result.iter_ambiguous_hits(self.config.bitscore_threshold_pct)
        )

        # Step 3: Classify
        taxonomic_call = self._classify_by_thresholds(
            novelty_index,
            placement_uncertainty
        )

        # Step 4: Return validated Pydantic model
        return ReadClassification(
            read_id=blast_result.read_id,
            best_match_genome=best_genome,
            top_hit_identity=top_hit_identity,
            novelty_index=novelty_index,
            placement_uncertainty=placement_uncertainty,
            num_ambiguous_hits=num_ambiguous_hits,
            taxonomic_call=taxonomic_call,
        )
```

#### 7.2.2 Placement Uncertainty Calculation

```python
def _calculate_placement_uncertainty(
    self,
    best_genome: str,
    ambiguous_hits_iter: Iterator[BlastHit],
) -> tuple[float, int]:
    """
    Calculate uncertainty via single-pass iterator.

    Optimization: Accepts iterator (not list) to enable early termination
    and avoid materialization of all ambiguous hits.
    """
    max_ani = 0.0
    num_ambiguous_hits = 0

    for hit in ambiguous_hits_iter:
        num_ambiguous_hits += 1
        secondary_genome = hit.genome_name

        if secondary_genome != best_genome:
            ani = self.ani_matrix.get_ani(best_genome, secondary_genome)
            if ani > max_ani:
                max_ani = ani

    # Edge cases
    if num_ambiguous_hits <= 1:
        return 0.0, num_ambiguous_hits  # Single hit, no uncertainty

    # The actual method also tracks the ANI to the first competing genome
    # and uses it when uncertainty_mode == "second" (the default).
    effective_ani = second_genome_ani if self.config.uncertainty_mode == "second" else max_ani

    if effective_ani == 0.0:
        return 100.0, num_ambiguous_hits  # Genome absent from ANI matrix

    return 100.0 - effective_ani, num_ambiguous_hits
```

**Iterator Advantage:**
Using `iter_ambiguous_hits()` allows early termination when bitscore drops below threshold, avoiding processing of 1000+ hits for some reads.

### 7.3 VectorizedClassifier: Polars-Native

**File:** `src/metadarkmatter/core/classification/classifiers/vectorized.py`

**Architecture:** Zero Python loops -- entire classification in Polars/Rust. This is the sole classifier used by all CLI paths (both batch and streaming).

#### 7.3.1 Vectorized Classification Pipeline

> **Note:** The code below is a simplified sketch of the approach, not the current implementation. See the referenced source file for the actual code.

The sketch shows the original threshold-only pipeline. The current `classify_file(blast_input, id_mapping=None, compute_qc=False)` additionally applies alignment filters and coverage weighting, uses a total-order tie-break for best and second genome, computes `uncertainty_mode="second"` by default, family validation (Off-target), inferred uncertainty for single-hit reads, Bayesian posteriors (`taxonomic_call`, `p_*`, `posterior_entropy`, `confidence_score`) and the legacy cascade (`legacy_call`), and conforms the result to a fixed output schema (section 12.2.3).

```python
class VectorizedClassifier:
    def classify_file(self, blast_path: Path) -> pl.DataFrame:
        """Fully vectorized classification—no Python loops."""
        # Step 1: Load and prepare data
        df = (
            StreamingBlastParser(blast_path).parse_lazy()
            .with_columns([extract_genome_name_expr()])
            .collect()
        )

        # Step 2: Find best hit per read (parallelized group_by)
        best_hits = df.group_by("qseqid").agg([
            pl.col("bitscore").max().alias("max_bitscore"),
            pl.col("pident")
              .filter(pl.col("bitscore") == pl.col("bitscore").max())
              .first()
              .alias("top_pident"),
            pl.col("genome_name")
              .filter(pl.col("bitscore") == pl.col("bitscore").max())
              .first()
              .alias("best_genome"),
        ])

        # Step 3: Calculate bitscore threshold (vectorized)
        best_hits = best_hits.with_columns([
            (pl.col("max_bitscore") * (self.config.bitscore_threshold_pct / 100.0))
            .alias("bitscore_threshold")
        ])

        # Step 4: Join to identify ambiguous hits
        df_with_threshold = df.join(
            best_hits.select(["qseqid", "bitscore_threshold", "best_genome"]),
            on="qseqid",
            how="left"
        )

        ambiguous = df_with_threshold.filter(
            pl.col("bitscore") >= pl.col("bitscore_threshold")
        )

        # Step 5: Count ambiguous hits per read
        hit_counts = (
            ambiguous.group_by("qseqid")
            .len()
            .rename({"len": "num_ambiguous_hits"})
        )

        # Step 6: ANI lookup via join (not Python loop!)
        secondary = ambiguous.filter(pl.col("genome_name") != pl.col("best_genome"))

        if not secondary.is_empty():
            secondary_with_ani = secondary.join(
                self._ani_lookup_symmetric,  # Pre-built ANI DataFrame
                left_on=["best_genome", "genome_name"],
                right_on=["genome1", "genome2"],
                how="left"
            )

            max_ani_per_read = (
                secondary_with_ani
                .group_by("qseqid")
                .agg(pl.col("ani").max().fill_null(0.0).alias("max_ani"))
            )
        else:
            max_ani_per_read = pl.DataFrame({"qseqid": [], "max_ani": []})

        # Step 7: Combine metrics
        result = (
            best_hits
            .join(hit_counts, on="qseqid", how="left")
            .join(max_ani_per_read, on="qseqid", how="left")
            .with_columns([
                pl.col("num_ambiguous_hits").fill_null(1),
                pl.col("max_ani").fill_null(0.0),
            ])
        )

        # Step 8: Calculate novelty and uncertainty (vectorized)
        cfg = self.config
        result = result.with_columns([
            pl.col("top_pident").clip(0.0, 100.0).alias("top_hit_identity"),
        ]).with_columns([
            (100.0 - pl.col("top_hit_identity")).alias("novelty_index"),
            pl.when(pl.col("num_ambiguous_hits") <= 1)
              .then(0.0)
              .when(pl.col("max_ani") == 0.0)
              .then(100.0)
              .otherwise(100.0 - pl.col("max_ani"))
              .alias("placement_uncertainty"),
        ])

        # Step 9: Classification via nested when/then (vectorized)
        result = result.with_columns([
            pl.when(pl.col("placement_uncertainty") >= cfg.uncertainty_conserved_min)  # >= 5.0
              .then(pl.lit("Conserved Region"))
            .when(
                (pl.col("novelty_index") < cfg.novelty_known_max) &  # < 4.0
                (pl.col("placement_uncertainty") < cfg.uncertainty_known_max)  # < 1.5
            )
              .then(pl.lit("Known Species"))
            .when(
                (pl.col("novelty_index") >= cfg.novelty_novel_species_min) &  # >= 4.0
                (pl.col("novelty_index") < cfg.novelty_novel_species_max) &  # < 20.0
                (pl.col("placement_uncertainty") < cfg.uncertainty_novel_species_max)  # < 1.5
            )
              .then(pl.lit("Novel Species"))
            .when(
                (pl.col("novelty_index") >= cfg.novelty_novel_genus_min) &  # >= 20.0
                (pl.col("novelty_index") <= cfg.novelty_novel_genus_max) &  # <= 25.0
                (pl.col("placement_uncertainty") < cfg.uncertainty_novel_genus_max)  # < 1.5
            )
              .then(pl.lit("Novel Genus"))
            .otherwise(pl.lit("Conserved Region"))
            .alias("taxonomic_call")
        ])

        # Step 10: Add is_novel flag
        result = result.with_columns([
            pl.col("taxonomic_call")
              .is_in(["Novel Species", "Novel Genus"])
              .alias("is_novel")
        ])

        return result.select([
            pl.col("qseqid").alias("read_id"),
            pl.col("best_genome").alias("best_match_genome"),
            "top_hit_identity",
            "novelty_index",
            "placement_uncertainty",
            "num_ambiguous_hits",
            "taxonomic_call",
            "is_novel",
        ])
```

#### 7.3.2 ANI Lookup DataFrame

**Key Innovation:** Pre-build ANI matrix as Polars DataFrame for join-based lookup.

```python
def _build_ani_lookup(self) -> None:
    """Build symmetric ANI lookup table for joins."""
    n = len(self.ani_matrix._genomes)

    # Get upper triangle indices (NumPy vectorized)
    idx1, idx2 = np.triu_indices(n, k=1)

    # Extract genome names and ANI values (NumPy indexing)
    genomes = np.array(self.ani_matrix._genomes)
    genome1_arr = genomes[idx1]
    genome2_arr = genomes[idx2]
    ani_arr = self.ani_matrix._ani_array[idx1, idx2].astype(np.float64)

    # Replace missing ANI values (0.0) with the default for distant organisms
    ani_arr[ani_arr == 0.0] = self.ani_matrix._default_ani

    # Build DataFrame from NumPy arrays (fast)
    self._ani_lookup = pl.DataFrame({
        "genome1": genome1_arr,
        "genome2": genome2_arr,
        "ani": ani_arr,
    })

    # Create symmetric version for bidirectional joins
    self._ani_lookup_symmetric = pl.concat([
        self._ani_lookup,
        self._ani_lookup.select([
            pl.col("genome2").alias("genome1"),
            pl.col("genome1").alias("genome2"),
            pl.col("ani"),
        ])
    ])
```

**Why This Works:**

Instead of 100M Python loops calling `ani_matrix.get_ani()`, we perform a single hash join:

```python
# Python loop (slow): 100M × 20 ns = 2000 ms
for hit in ambiguous_hits:
    ani = ani_matrix.get_ani(best_genome, hit.genome_name)

# Polars join (fast): ~50 ms
secondary_with_ani = ambiguous.join(
    ani_lookup_symmetric,
    left_on=["best_genome", "genome_name"],
    right_on=["genome1", "genome2"]
)
```

**Speedup:** 40x faster ANI lookups via join vs loop

**Benchmark (10M reads, 1000 genomes):**
- Standard: 245s
- Vectorized: 14s
- **Speedup: 17.5x**

---

## 8. CLI Interface

### 8.1 Command Structure

**File:** `src/metadarkmatter/cli/score.py` (the `mdm` and `metadarkmatter` entry points are interchangeable)

```
metadarkmatter (main entry point; global --dry-run, --log-level, --debug)
├── score
│   ├── classify       # Single-file classification
│   ├── batch          # Multi-file batch processing
│   ├── export-config  # Write editable YAML config
│   ├── extract-novel  # Extract candidate novel taxa
│   ├── sensitivity    # Threshold sweep
│   └── evaluate       # Score against ground truth
├── report (generate, multi, summarize, vendor-d3, vendor-cdn-sri)
├── run (end-to-end pipeline)
└── aai, ani, blast, blastx, mmseqs2, kraken2, map, proteins, download,
    util, tree, validate, visualize, doctor
```

### 8.2 Score Classify Command

**Signature:**
```bash
metadarkmatter score classify \
  --alignment ALIGNMENT_PATH \
  --ani ANI_MATRIX_PATH \
  --output OUTPUT_PATH \
  [--summary SUMMARY_JSON] \
  [--format {csv,parquet}] \
  [--streaming [--chunk-size N]] \
  [--metadata genome_metadata.tsv] \
  [--config scoring.yaml] \
  [--verbose]
```

See `mdm score classify --help` for the full option list (alignment filters, coverage weighting, `--id-mapping`, `--target-family`, `--qc-output`, `--adaptive-thresholds`, `--alignment-mode`, `--preset`, `--dry-run`, ...).

**Implementation:**

> **Note:** The code below is a simplified sketch of the approach, not the current implementation. See the referenced source file for the actual code.

The command builds a `ClassificationRequest` and calls `run_classification()` in `core/classification/runner.py`, which wraps the steps below.

```python
@app.command(name="classify")
def classify(
    blast: Path,
    ani: Path,
    output: Path,
    summary: Path | None = None,
    output_format: str = "csv",
    streaming: bool = False,
    verbose: bool = False,
) -> None:
    """Classify metagenomic reads using ANI-weighted placement."""
    # Load ANI matrix
    ani_matrix = ANIMatrix.from_file(ani)
    config = ScoringConfig()

    # All CLI paths use VectorizedClassifier
    vectorized = VectorizedClassifier(ani_matrix, config=config)

    if streaming:
        # Streaming mode: direct-to-file, bounded memory
        num_classified = vectorized.stream_to_file(
            blast_path=blast,
            output_path=output,
            output_format=output_format,
            progress_callback=streaming_progress,
        )
    else:
        # Batch mode: load full dataset, auto-parallelized by Polars
        df = vectorized.classify_file(blast)
        if output_format == "parquet":
            df.write_parquet(output)
        else:
            df.write_csv(output)
        num_classified = len(df)

    # Generate summary if requested
    if summary and num_classified > 0:
        summary_obj = _generate_summary(df)
        summary_obj.to_json(summary)
        _display_summary_table(summary_obj)
```

### 8.3 Output Format Options

#### 8.3.1 CSV Format

**Pros:**
- Human-readable
- Universal compatibility
- Easy inspection with `head`, `grep`

**Cons:**
- Large file size (10-20 GB for 100M reads)
- Slow I/O (write: ~300 MB/s, read: ~400 MB/s)
- No compression

**Use when:** Small-medium datasets, compatibility required

#### 8.3.2 Parquet Format

**Pros:**
- 10x smaller (1-2 GB for 100M reads with zstd compression)
- 10x faster I/O (write: ~2 GB/s, read: ~3 GB/s)
- Column-oriented for efficient filtering
- Schema embedded in file

**Cons:**
- Requires Polars/Pandas/Arrow for reading
- Not human-readable
- Slight overhead for very small files

**Use when:** Large datasets, production pipelines, storage-constrained environments

**Recommendation:** Default to Parquet for files > 10M reads.

### 8.4 Rich CLI Output

The CLI uses Rich library for enhanced terminal output:

```python
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, BarColumn
from rich.table import Table

console = Console()

# Spinner for indeterminate tasks
with Progress(SpinnerColumn(), TextColumn(...)) as progress:
    progress.add_task("Loading ANI matrix...", total=None)

# Progress bar with ETA for streaming
with Progress(
    SpinnerColumn(),
    TextColumn("[progress.description]{task.description}"),
    BarColumn(),
    TaskProgressColumn(),
    TimeElapsedColumn(),
    TimeRemainingColumn(),
) as progress:
    task = progress.add_task("Classifying reads...", total=estimated_rows)
    # Update during processing
    progress.update(task, completed=rows_processed)

# Summary table
table = Table(title="Classification Summary")
table.add_column("Category", style="cyan")
table.add_column("Count", justify="right")
table.add_column("Percentage", justify="right")
table.add_row("Novel Species", f"{novel_species:,}", f"{pct:.2f}%")
console.print(table)
```

**Example Output:**
```
Metadarkmatter ANI-Weighted Classification

Loading ANI matrix... ✓
Loaded ANI matrix for 1,245 genomes

Classifying reads (streaming)...
█████████████████████░░░░░░░ 75% 45,234,891 reads @ 2,145,678 rows/s
Elapsed: 0:02:15 Remaining: 0:00:45

╭─────────────────────┬─────────┬────────────╮
│ Category            │   Count │ Percentage │
├─────────────────────┼─────────┼────────────┤
│ Known Species       │ 35,234  │      70.5% │
│ Novel Species       │  8,123  │      16.2% │
│ Novel Genus         │  2,045  │       4.1% │
│ Conserved Regions   │  4,598  │       9.2% │
│ Total Novel         │ 10,168  │      20.3% │
╰─────────────────────┴─────────┴────────────╯

Results written to: /path/to/output.parquet
Summary written to: /path/to/summary.json
```

---

## 9. Performance Characteristics

### 9.1 Benchmark Environment

**Hardware:**
- CPU: Intel Xeon E5-2680 v3 (12 cores, 24 threads @ 2.5 GHz)
- RAM: 64 GB DDR4
- Storage: NVMe SSD (3.5 GB/s read, 2.5 GB/s write)
- OS: Linux 5.15

**Dataset:**
- BLAST results: 100,000,000 alignments (10 GB TSV)
- Reads: 25,000,000 unique reads (avg 4 hits/read)
- Genomes: 1,000 reference genomes
- ANI matrix: 1,000 × 1,000 (8 MB)

### 9.2 End-to-End Performance

| Mode | Threads | Time | Throughput | Memory Peak | Output Size |
|------|---------|------|------------|-------------|-------------|
| VectorizedClassifier (CSV) | Auto | 1m 52s | 892,857 reads/s | 6.2 GB | 2.8 GB CSV |
| VectorizedClassifier (Parquet) | Auto | 1m 28s | 1,136,364 reads/s | 6.2 GB | 285 MB Parquet |
| VectorizedClassifier (Streaming) | Auto | 2m 15s | 740,741 reads/s | 1.8 GB | 285 MB Parquet |

**Key Observations:**

1. **Batch + Parquet:** Best overall throughput
2. **Streaming mode:** Bounded memory at ~1.8 GB regardless of file size, modest throughput reduction
3. **Parquet I/O:** Saves ~30% wall-clock time on output writing

### 9.3 Scaling Analysis

#### 9.3.1 Read Count Scaling (Vectorized Classifier)

| Reads | BLAST Alignments | Time | Throughput | Memory |
|-------|------------------|------|------------|--------|
| 1M | 4M | 5s | 200k reads/s | 800 MB |
| 10M | 40M | 42s | 238k reads/s | 3.2 GB |
| 25M | 100M | 1m 52s | 223k reads/s | 6.2 GB |
| 50M | 200M | 4m 18s | 194k reads/s | 12.5 GB |
| 100M | 400M | 9m 35s | 174k reads/s | 24.8 GB |

**Scaling Characteristics:**

- **Sub-linear time complexity:** Polars' parallel aggregations scale well
- **Linear memory usage:** ~250 bytes per read in RAM
- **Throughput degradation:** Slight decrease at 100M+ due to cache pressure

**Recommendation:** For >50M reads, use streaming mode to cap memory at 8 GB:

```bash
metadarkmatter score classify --streaming --alignment huge.tsv.gz --ani ani.csv --output out.parquet
```

#### 9.3.2 Genome Count Scaling (ANI Matrix)

| Genomes | ANI Matrix Size | Load Time | Memory | Lookup Time |
|---------|----------------|-----------|--------|-------------|
| 100 | 100 KB | 15 ms | 80 KB | 6 ns |
| 500 | 2 MB | 85 ms | 2 MB | 7 ns |
| 1,000 | 8 MB | 180 ms | 8 MB | 7.5 ns |
| 5,000 | 200 MB | 1.2 s | 200 MB | 9 ns |
| 10,000 | 800 MB | 5.8 s | 800 MB | 12 ns |

**Key Insights:**

1. **NumPy array scales quadratically with genome count:** Memory = `8n²` bytes (float64)
2. **Lookup time constant:** O(1) array indexing regardless of size
3. **Load time dominated by parsing:** CSV read is bottleneck, not NumPy conversion

### 9.4 I/O Performance

#### 9.4.1 BLAST File Reading

| Format | Size | Read Method | Time | Throughput |
|--------|------|-------------|------|------------|
| TSV | 10 GB | Pandas | 180s | 55 MB/s |
| TSV | 10 GB | Polars eager | 12s | 833 MB/s |
| TSV | 10 GB | Polars streaming | 8s | 1.25 GB/s |
| TSV.GZ | 2.5 GB | Polars streaming | 45s | 222 MB/s (decompressed) |

**Recommendations:**

- **Uncompressed TSV:** Fastest for local SSD storage
- **Gzipped TSV:** Best for network storage or limited disk space (4x compression)
- **Streaming:** Always use for files >1 GB to avoid RAM exhaustion

#### 9.4.2 Output Writing

| Format | Rows | Write Time | Size | Compression Ratio |
|--------|------|------------|------|-------------------|
| CSV | 25M | 62s | 2.8 GB | 1x (baseline) |
| Parquet (Snappy) | 25M | 14s | 420 MB | 6.7x |
| Parquet (Zstd) | 25M | 18s | 285 MB | 9.8x |
| Parquet (Zstd L3) | 25M | 25s | 245 MB | 11.4x |

**Best Practices:**

1. **Development/exploration:** CSV for human readability
2. **Production pipelines:** Parquet with Zstd compression (default in package)
3. **Long-term storage:** Parquet Zstd L3 for maximum compression

### 9.5 Memory Profiling

#### 9.5.1 VectorizedClassifier Batch Memory Breakdown

**Total Peak:** 6.2 GB for 25M reads

| Component | Memory | Percentage |
|-----------|--------|------------|
| Full BLAST DataFrame | 4.8 GB | 77.4% |
| ANI lookup DataFrame | 120 MB | 1.9% |
| Intermediate joins | 800 MB | 12.9% |
| Results DataFrame | 480 MB | 7.7% |
| ANI matrix | 8 MB | 0.1% |

**Trade-off:** Higher memory usage but 10x faster due to vectorization

#### 9.5.2 VectorizedClassifier Streaming Memory Profile

**Constant Peak:** 1.8 GB regardless of file size

| Component | Memory | Percentage |
|-----------|--------|------------|
| Partition buffer (5M rows) | 1.2 GB | 66.7% |
| Pending read boundary | 50 MB | 2.8% |
| ANI lookup DataFrame | 120 MB | 6.7% |
| Classification buffer | 300 MB | 16.7% |
| ANI matrix | 8 MB | 0.4% |
| Output write buffer | 120 MB | 6.7% |

**Key Feature:** Memory usage does not grow with total file size

---

## 10. Testing and Quality Assurance

### 10.1 Test Suite Organization

**Location:** `tests/`

```
tests/
├── conftest.py, factories.py        # Shared fixtures and synthetic data builders
├── unit/                            # Unit tests (e.g. test_blast_models.py,
│                                    #   test_classification.py, test_ani_matrix_class.py,
│                                    #   test_bayesian.py, test_classifier_equivalence.py,
│                                    #   test_clustering_determinism.py, test_cli_*.py)
├── integration/                     # Integration tests (test_cli.py,
│                                    #   test_family_validation.py,
│                                    #   test_pipeline_determinism.py, ...)
├── e2e/                             # End-to-end tests (marker `e2e`)
├── fixtures/                        # GTDB response fixtures, test genomes
├── utils/                           # Shared assertions
└── benchmark_classifiers.py         # Classifier benchmark script
```

### 10.2 Coverage Requirements

**Configuration:** `pyproject.toml`

```toml
[tool.coverage.run]
source = ["src/metadarkmatter"]
branch = true

[tool.coverage.report]
exclude_lines = ["pragma: no cover", "if TYPE_CHECKING:", "raise NotImplementedError"]
fail_under = 78
```

**Historical Coverage (as of v0.1.0; not updated):**

| Module | Line Coverage | Branch Coverage |
|--------|---------------|-----------------|
| `models/blast.py` | 95% | 88% |
| `models/classification.py` | 92% | 85% |
| `models/config.py` | 88% | 82% |
| `core/classification/classifiers/` | 87% | 79% |
| `core/parsers.py` | 91% | 84% |
| `cli/score.py` | 78% | 68% |
| **Overall** | **87%** | **81%** |

### 10.3 Key Test Scenarios

#### 10.3.1 Edge Case Testing

```python
def test_blast_hit_pident_clamping():
    """Ensure pident > 100 is clamped to 100."""
    hit = BlastHit(
        qseqid="read1",
        sseqid="GCF_000123.1",
        pident=100.001,  # Over 100 due to BLAST rounding
        # ... other fields
    )
    assert hit.pident == 100.0

def test_empty_blast_result():
    """Handle read with no BLAST hits."""
    result = BlastResult(read_id="read1", hits=tuple())
    classifier = ANIWeightedClassifier(ani_matrix, config=config)
    classification = classifier.classify_read(result)
    assert classification is None

def test_single_hit_no_uncertainty():
    """Read with only one hit has zero placement uncertainty."""
    result = BlastResult(read_id="read1", hits=(single_hit,))
    classification = classifier.classify_read(result)
    assert classification.placement_uncertainty == 0.0
    assert classification.num_ambiguous_hits == 1
```

#### 10.3.2 Performance Regression Tests

```python
import pytest
from time import time

@pytest.mark.slow
def test_vectorized_classifier_performance(benchmark_blast_file):
    """Ensure vectorized classifier meets performance targets."""
    start = time()

    classifier = VectorizedClassifier(ani_matrix, config=config)
    df = classifier.classify_file(benchmark_blast_file)

    elapsed = time() - start
    throughput = len(df) / elapsed

    # Target: >200k reads/second
    assert throughput > 200_000, f"Throughput {throughput:.0f} reads/s below target"
```

#### 10.3.3 Data Integrity Tests

Cross-classifier equivalence is tested in `tests/unit/test_classifier_equivalence.py`: `ANIWeightedClassifier.classify_read()` (applied per read via `StreamingBlastParser.iter_reads()`) is compared against `VectorizedClassifier.classify_file()` on `best_match_genome`, `novelty_index` and `taxonomic_call` (both are Bayesian-primary; base "Ambiguous" vs vectorized "Conserved Region" is accepted because only the vectorized path has `ambiguity_scope` for Stage 2). Output determinism is guarded by `tests/integration/test_pipeline_determinism.py` (SHA-256 equality of two `score classify` runs).

```python
def test_base_vs_vectorized(blast_file, ani_matrix, config):
    base_clf = ANIWeightedClassifier(ani_matrix=ani_matrix, config=config)
    vec_clf = VectorizedClassifier(ani_matrix=ani_matrix, config=config)

    base_df = pl.DataFrame([
        c.to_dict()
        for r in StreamingBlastParser(blast_file).iter_reads()
        if (c := base_clf.classify_read(r)) is not None
    ]).sort("read_id")
    vec_df = vec_clf.classify_file(blast_file).sort("read_id")

    for col in ["best_match_genome", "novelty_index"]:
        assert base_df[col].to_list() == vec_df[col].to_list()
    allowed = {("Ambiguous", "Conserved Region")}
    for b, v in zip(base_df["taxonomic_call"], vec_df["taxonomic_call"]):
        assert b == v or (b, v) in allowed
```

### 10.4 Continuous Integration

**GitHub Actions Workflow:** `.github/workflows/test.yml` (on push/PR to `main`; ubuntu and macOS, Python 3.11 and 3.12, plus a non-blocking 3.13 leg). Abridged:

```yaml
      - name: Install package and dev dependencies (frozen lockfile)
        run: uv sync --frozen --extra dev

      - name: Lint (ruff)
        run: uv run ruff check src tests

      - name: Type check (mypy)
        run: uv run mypy src/metadarkmatter

      - name: Run unit tests with coverage gate
        run: |
          uv run pytest -m "not slow and not requires_blast and not requires_mmseqs2 and not e2e" \
                 --cov=metadarkmatter --cov-report=term-missing --cov-report=xml \
                 --cov-fail-under=78
```

The coverage XML is uploaded as a workflow artifact (ubuntu, Python 3.12).

---

## 11. Deployment Considerations

### 11.1 Installation Methods

#### 11.1.1 Development Installation

```bash
git clone https://github.com/FOI-Bioinformatics/metadarkmatter.git
cd metadarkmatter
uv sync --frozen --extra dev   # pinned environment from uv.lock
# or: pip install -e ".[dev]"
```

**Use case:** Active development, contributing to codebase

#### 11.1.2 Production Installation

```bash
pip install metadarkmatter
```

**Use case:** End users, HPC deployment scripts

#### 11.1.3 Conda Environment

```bash
conda create -n metadarkmatter python=3.11
conda activate metadarkmatter
pip install metadarkmatter
```

**Use case:** Isolated environments, managing Python versions

### 11.2 HPC Cluster Deployment

#### 11.2.1 SLURM Job Script

```bash
#!/bin/bash
#SBATCH --job-name=metadarkmatter
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --partition=compute

# Load environment
module load python/3.11
source /path/to/venv/bin/activate

# Set Polars thread count
export POLARS_MAX_THREADS=16

# Run classification
metadarkmatter score classify \
  --alignment /scratch/${SLURM_JOB_ID}/blast_results.tsv.gz \
  --ani /data/references/ani_matrix.csv \
  --output /results/${SAMPLE_ID}_classifications.parquet \
  --summary /results/${SAMPLE_ID}_summary.json \
  --format parquet \
  --verbose
```

**Resource Guidelines:**

| File Size | Reads | Cores | RAM | Time |
|-----------|-------|-------|-----|------|
| Small (< 1 GB) | < 1M | 4 | 8 GB | < 5 min |
| Medium (1-10 GB) | 1-10M | 8 | 16 GB | 10-30 min |
| Large (10-50 GB) | 10-50M | 16 | 32 GB | 30-120 min |
| Very Large (>50 GB) | >50M | 24 | 64 GB | 2-6 hrs |

#### 11.2.2 Array Job for Batch Processing

```bash
#!/bin/bash
#SBATCH --array=1-100
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

# Get sample name from array index
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

metadarkmatter score classify \
  --alignment /data/blast/${SAMPLE}.tsv.gz \
  --ani /data/references/ani_matrix.csv \
  --output /results/${SAMPLE}_classifications.parquet \
  --format parquet
```

### 11.3 Container Image

The repository root contains a `Containerfile` (based on `mambaorg/micromamba`) that installs the external tools from bioconda with pinned versions (see `docs/EXTERNAL_TOOLS.md`) and the Python package from the locked `uv.lock` environment. The entry point is the `metadarkmatter` CLI.

**Build and run:**

```bash
podman build -t metadarkmatter:dev .     # or: docker build -t metadarkmatter:dev .

podman run --rm -v /data:/data metadarkmatter:dev \
  score classify \
  --alignment /data/blast_results.tsv.gz \
  --ani /data/ani_matrix.csv \
  --output /data/classifications.csv
```

### 11.4 Configuration Management

#### 11.4.1 ScoringConfig

Classification thresholds are configured via `ScoringConfig`:

```python
from metadarkmatter.models.config import ScoringConfig

# Use defaults (recommended for most analyses)
config = ScoringConfig()

# Custom thresholds
config = ScoringConfig(
    bitscore_threshold_pct=95.0,
    novelty_known_max=4.0,
    novelty_novel_species_min=4.0,
    novelty_novel_species_max=20.0,
    novelty_novel_genus_min=20.0,
    novelty_novel_genus_max=25.0,
    uncertainty_known_max=1.5,
    uncertainty_novel_species_max=1.5,
    uncertainty_novel_genus_max=1.5,
    uncertainty_conserved_min=5.0,
)
```

#### 11.4.2 Polars Thread Control

```bash
# Control Polars internal parallelism via environment variable
export POLARS_MAX_THREADS=16

metadarkmatter score classify ...
```

---

## 12. Appendices

### 12.1 Glossary

**ANI (Average Nucleotide Identity):** Pairwise genome similarity metric calculated as average nucleotide identity across aligned regions. Used to define bacterial species boundaries (≥95-96% ANI = same species).

**Ambiguous Hit:** BLAST alignment with bitscore within threshold (default 95%) of top-scoring hit, indicating competitive placement across multiple genomes.

**Bitscore:** BLAST alignment score normalized for database size and scoring matrix. Higher values indicate better alignments.

**Competitive Recruitment:** Alignment strategy where metagenomic reads compete to match best reference genome among multiple candidates, revealing placement uncertainty.

**eDNA (Environmental DNA):** Genetic material extracted from environmental samples (air, water, soil) without isolating organisms, enabling culture-independent microbial community analysis.

**Frozen Model:** Pydantic model with `frozen=True` configuration, making instances immutable after creation. Ensures thread safety and prevents accidental modifications.

**Lazy Evaluation:** Deferred execution of operations until results are needed (e.g., Polars `LazyFrame`), enabling query optimization and reduced memory usage.

**Metagenomic Dark Matter:** Uncultured microbial taxa lacking representation in reference databases, detectable only through culture-independent sequencing methods.

**Novelty Index (N):** Metric quantifying sequence divergence from closest reference genome: `N = 100 - TopHitIdentity`. Higher values indicate greater evolutionary distance.

**Placement Uncertainty (U):** Metric quantifying ambiguity in taxonomic assignment: `U = 100 - ANI(TopHitGenome, SecondBestGenome)` by default (`uncertainty_mode="second"`), or `U = 100 - max(ANI(TopHit, SecondaryHits))` with `uncertainty_mode="max"`. High values indicate conserved regions matching distant genomes.

**Polars:** High-performance DataFrame library implemented in Rust, offering 10-100x speedups over pandas for large-scale data processing.

**Streaming:** Processing large files in fixed-size chunks with bounded memory usage, enabling analysis of datasets exceeding available RAM.

**Vectorization:** Executing operations on entire arrays/columns simultaneously using CPU SIMD instructions, avoiding slow Python loops.

### 12.2 File Format Specifications

#### 12.2.1 BLAST Tabular Output

**Format:** Tab-delimited text (TSV)
**Extension:** `.tsv`, `.blast`, `.tsv.gz` (gzipped)
**Header:** None (comment lines starting with `#` are ignored)

**Column Order (12 columns required; an optional 13th column `qlen` is used for coverage when present):**

| Index | Field | Type | Description |
|-------|-------|------|-------------|
| 0 | qseqid | string | Query sequence identifier (read ID) |
| 1 | sseqid | string | Subject sequence identifier (genome contig) |
| 2 | pident | float | Percent identity (0-100) |
| 3 | length | int | Alignment length in base pairs |
| 4 | mismatch | int | Number of mismatches |
| 5 | gapopen | int | Number of gap openings |
| 6 | qstart | int | Start position in query (1-based) |
| 7 | qend | int | End position in query (1-based) |
| 8 | sstart | int | Start position in subject (1-based) |
| 9 | send | int | End position in subject (1-based) |
| 10 | evalue | float | Expectation value |
| 11 | bitscore | float | Bit score |
| 12 | qlen | int | Query length (optional) |

Reverse-strand hits with qstart > qend (MMseqs2, DIAMOND) are normalized. `sseqid` is expected in the standardized `{accession}|{contig_id}` form written by `mdm blast makedb` / `mdm mmseqs2 makedb`; otherwise supply `--id-mapping`.

**BLAST Command** (`mdm blast align` runs an equivalent command):

```bash
blastn -query reads.fasta \
       -db reference_genomes.fna \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
       -out results.tsv \
       -num_threads 16 \
       -max_target_seqs 1000
```

#### 12.2.2 ANI Matrix CSV/TSV

**Format:** Comma-separated (CSV) or tab-delimited (TSV)
**Extension:** `.csv`, `.tsv`
**Header:** Required (first row contains genome names)

**Structure:**

```csv
genome,GCF_000001.1,GCF_000002.1,GCF_000003.1
GCF_000001.1,100.0,95.3,87.2
GCF_000002.1,95.3,100.0,88.5
GCF_000003.1,87.2,88.5,100.0
```

**Requirements:**
- First column: Genome identifiers (must match BLAST `sseqid` prefixes)
- Subsequent columns: ANI values (0-100)
- Square matrix: `n_rows = n_columns - 1`
- Diagonal values: 100.0 (self-ANI)
- Symmetric: `ANI[i,j] = ANI[j,i]`
- Percent scale (0-100); a 0-1 scaled matrix is rejected

**Generation Tools:** `mdm ani compute` (fastANI, falling back to skani) writes this format directly. Other tools:

- **FastANI:** Genome-vs-genome ANI calculation
  ```bash
  fastANI --ql genome_list.txt --rl genome_list.txt -o ani_results.txt
  # Convert to matrix format using custom script
  ```

- **SKA2:** Whole-genome alignment-based ANI
  ```bash
  ska distance *.fasta -o ani_matrix.csv --format matrix
  ```

#### 12.2.3 Classification Output CSV

**Format:** Comma-separated values
**Extension:** `.csv`
**Header:** Required

**Schema:**

The column set, order and dtypes are fixed (`VectorizedClassifier._output_schema()`), including for empty and all-Off-target inputs and for every streaming batch.

| Column | Type | Description |
|--------|------|-------------|
| read_id | string | Original read identifier from BLAST |
| best_match_genome | string | Genome with highest (coverage-weighted) bitscore |
| top_hit_identity | float64 | Percent identity of best hit (0-100) |
| novelty_index | float64 | 100 - top_hit_identity |
| placement_uncertainty | float64 | 100 - ANI to second-best genome (default mode); null for Off-target |
| genus_uncertainty | float64 | Uncertainty from genus-level hits (90% bitscore threshold) |
| ambiguity_scope | string | e.g. unambiguous / within genus / across genera |
| num_ambiguous_hits | int64 | Hits within 95% of top bitscore |
| second_hit_identity | float64 | Identity of best hit to a different genome |
| identity_gap | float64 | top_hit_identity - second_hit_identity |
| confidence_score | float64 | (1 - entropy / log2(6)) * 100; 10.0 for Off-target |
| taxonomic_call | string | Bayesian MAP category (after Stage 2 refinement) |
| legacy_call | string | Hard-threshold cascade category |
| diversity_status | string | Known / Novel / Uncertain |
| is_novel | boolean | True if Novel Species or Novel Genus |
| low_confidence | boolean | Low-confidence flag |
| p_known_species, p_novel_species, p_novel_genus, p_species_boundary, p_ambiguous, p_unclassified | float64 | Posterior probabilities |
| posterior_entropy | float64 | Shannon entropy of the posteriors (0 to 2.585) |
| inferred_uncertainty | float64 | Inferred uncertainty for single-hit reads |
| uncertainty_type | string | measured / inferred |

Additional columns: `alignment_quality`, `identity_confidence`, `placement_confidence`, `discovery_score` with `--include-legacy-scores`; `family_bitscore_ratio`, `family_identity_gap`, `in_family_hit_fraction`, `external_best_genome`, `external_best_identity`, `in_family_identity`, `in_family_novelty_index` when a target family is set; `species`, `genus`, `family` when `--metadata` is given (not in streaming mode).

#### 12.2.4 Classification Output Parquet

**Format:** Apache Parquet (columnar binary)
**Extension:** `.parquet`
**Compression:** Zstd (default)
**Schema:** Identical to CSV (see above)

**Advantages:**
- 10x smaller file size
- 10x faster I/O
- Embedded schema
- Column filtering without full read

**Reading in Python:**

```python
import polars as pl

df = pl.read_parquet("classifications.parquet")

# Filter novel reads only
novel_reads = df.filter(pl.col("is_novel") == True)
```

### 12.3 Performance Tuning Guide

#### 12.3.1 Memory Constraints

**Scenario:** Running on system with 16 GB RAM, processing 50M reads.

**Problem:** Vectorized classifier requires ~25 GB RAM.

**Solution:** Use streaming mode:

```bash
metadarkmatter score classify \
  --alignment large_file.tsv.gz \
  --ani ani_matrix.csv \
  --output results.parquet \
  --streaming
```

**Effect:** Bounds memory by processing in partitions of `--chunk-size` rows (default 1,000,000; e.g. `--chunk-size 500000` for tighter memory). Streaming requires each read's hits to be contiguous in the file and cannot be combined with `--id-mapping` / `--genomes`.

#### 12.3.2 CPU Utilization

**Scenario:** 32-core HPC node, want maximum parallelism.

```bash
export POLARS_MAX_THREADS=32

metadarkmatter score classify --alignment file.tsv --ani ani.csv --output results.parquet
```

**Effect:** Polars uses all 32 cores for internal parallelization within VectorizedClassifier. No additional flags needed -- Polars automatically detects available cores.

#### 12.3.3 I/O Bottlenecks

**Scenario:** Slow network-attached storage (NAS) limits throughput.

**Problem:** Uncompressed 10 GB BLAST file takes 200s to read over network.

**Solution 1:** Use gzipped BLAST files (4x smaller)

```bash
# Compress BLAST results
gzip -c blast_results.tsv > blast_results.tsv.gz

# Polars automatically handles decompression
metadarkmatter score classify --alignment blast_results.tsv.gz ...
```

**Effect:** 2.5 GB file transfers in 50s (4x faster), decompression adds ~30s (net gain: 120s).

**Solution 2:** Stage to local SSD

```bash
# Copy to /scratch on compute node
cp /nfs/data/blast_results.tsv.gz /scratch/

# Run classification from local disk
metadarkmatter score classify --alignment /scratch/blast_results.tsv.gz ...

# Copy results back
cp /scratch/classifications.parquet /nfs/results/
```

**Effect:** Eliminates network I/O during processing, 10x faster reads from SSD.

### 12.4 Troubleshooting Common Issues

#### 12.4.1 Out of Memory Errors

**Symptom:**

```
MemoryError: Unable to allocate 12.5 GiB for an array with shape (1000, 1000000) and data type float64
```

**Diagnosis:** Vectorized classifier loading full dataset into RAM.

**Solutions:**

1. **Use streaming mode:**
   ```bash
   metadarkmatter score classify --streaming --alignment large_file.tsv.gz --ani ani.csv --output results.parquet
   ```
   Streaming mode caps memory at ~1.8 GB by processing in bounded partitions.

2. **Reduce Polars thread count** (each thread uses memory):
   ```bash
   export POLARS_MAX_THREADS=4
   ```

#### 12.4.2 ANI Matrix Validation Errors

**Symptom:**

```
ValueError: ANI matrix row/column mismatch:
  Rows missing from columns: ['GCF_000123.1', 'GCF_000456.1']
```

**Diagnosis:** Row genome names don't match column headers.

**Solution:** Ensure first column and header contain identical genome sets:

```python
import polars as pl

df = pl.read_csv("ani_matrix.csv")
genome_col = df.columns[0]

row_genomes = set(df[genome_col])
col_genomes = set(df.columns[1:])

print("Missing in rows:", col_genomes - row_genomes)
print("Missing in columns:", row_genomes - col_genomes)
```

#### 12.4.3 Slow Classification Performance

**Symptom:** Processing 10M reads takes >30 minutes (expected: <5 minutes).

**Diagnosis:** Polars may not be using all available cores, or I/O is the bottleneck.

**Solutions:**

1. **Ensure Polars uses all cores:**
   ```bash
   export POLARS_MAX_THREADS=16
   metadarkmatter score classify --alignment file.tsv --ani ani.csv --output results.csv
   ```

2. **Use compressed input to reduce I/O:**
   ```bash
   gzip -c blast_results.tsv > blast_results.tsv.gz
   metadarkmatter score classify --alignment blast_results.tsv.gz --ani ani.csv --output results.parquet
   ```

3. **Use Parquet output** for faster writes:
   ```bash
   metadarkmatter score classify --alignment file.tsv --ani ani.csv --output results.parquet --format parquet
   ```

#### 12.4.4 Type Validation Failures

**Symptom:**

```
ValidationError: 1 validation error for BlastHit
pident
  Input should be a valid number, unable to parse string as a number [type=float_parsing, input_value='NA', input_type=str]
```

**Diagnosis:** BLAST output contains non-numeric values (e.g., "NA", "N/A").

**Solution:** Pre-filter BLAST file to remove invalid lines:

```bash
# Remove lines with NA values
grep -v "NA" blast_results.tsv > blast_results_clean.tsv
```

Or handle in Python:

```python
from metadarkmatter.core.parsers import StreamingBlastParser

# Skip invalid lines
for result in parser.iter_reads():
    try:
        classification = classifier.classify_read(result)
    except ValidationError as e:
        print(f"Skipping invalid read: {result.read_id}")
        continue
```

### 12.5 References

#### 12.5.1 Scientific Publications

1. Konstantinidis, K. T., & Tiedje, J. M. (2005). Genomic insights that advance the species definition for prokaryotes. *PNAS*, 102(7), 2567-2572.

2. Jain, C., Rodriguez-R, L. M., Phillippy, A. M., Konstantinidis, K. T., & Aluru, S. (2018). High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. *Nature Communications*, 9(1), 1-8.

3. Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P., & Tyson, G. W. (2015). CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. *Genome Research*, 25(7), 1043-1055.

#### 12.5.2 Tool Documentation

- **Polars:** https://pola-rs.github.io/polars/
- **Pydantic:** https://docs.pydantic.dev/
- **Typer:** https://typer.tiangolo.com/
- **BLAST+:** https://www.ncbi.nlm.nih.gov/books/NBK279690/
- **FastANI:** https://github.com/ParBLiSS/FastANI

#### 12.5.3 Code Repository

- **GitHub:** https://github.com/FOI-Bioinformatics/metadarkmatter
- **Issues:** https://github.com/FOI-Bioinformatics/metadarkmatter/issues
- **Documentation:** `docs/` in the repository

---

**Document Version:** 1.1
**Package Version:** 0.2.0
**Generated:** 2026-02-20 (facts updated 2026-11-02)
**Maintainer:** Metadarkmatter Development Team
