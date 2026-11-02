# Performance Guide

Optimizing metadarkmatter for different dataset sizes and hardware configurations.

## Classification Modes

All classification uses the Polars-based vectorized engine (`VectorizedClassifier`) with automatic parallelization across CPU cores.

### Default Mode

```bash
metadarkmatter score classify --alignment sample.blast.tsv.gz --ani ani.csv --output out.csv
```

| Metric | Value |
|--------|-------|
| Best for | Up to ~5M alignments, or more when RAM allows |
| RAM usage | Scales with data; pre-run estimate about 1.5x the uncompressed alignment size |
| Speed | Automatic multi-core parallelization |

Loads all data into memory, performs vectorized classification in Polars' Rust backend. Before loading, the classifier estimates memory as about 1.5x the uncompressed file size (about 15x a gzipped file) and logs a warning if this may exceed available memory.

### Streaming Mode

```bash
metadarkmatter score classify --alignment sample.blast.tsv.gz --ani ani.csv --output out.csv \
    --streaming --chunk-size 500000
```

| Metric | Value |
|--------|-------|
| Best for | >5M alignments or memory-constrained systems |
| RAM usage | Bounded by `--chunk-size` |
| Speed | Comparable to default for large files |

Processes data in batches of `--chunk-size` rows (default 1,000,000), writing results incrementally, so memory use does not grow with file size. Restrictions: each read's hits must be contiguous in the input (BLAST and MMseqs2 output satisfies this; otherwise the run stops with an error), `--id-mapping`/`--genomes` are rejected, and with `--metadata` the species/genus/family columns are omitted (representative mapping and family inference still apply). Classifications are expected to match the in-memory path (checked by `tests/e2e/test_mode_equivalence.py`).

### When to Use Streaming

```
Dataset size?
    |
    +--< ~5M alignments --> Default (loads all into memory)
    |
    +-->= ~5M alignments or limited RAM --> --streaming
```

## Output Format Selection

### CSV (Default)

```bash
metadarkmatter score classify --output results.csv
```

- Human-readable
- Compatible with Excel, R, any tool
- Larger file size

### Parquet

```bash
metadarkmatter score classify --output results.parquet --format parquet
```

- 10x smaller file size
- 10x faster write speed
- 15x faster read speed
- Requires Parquet-compatible tools

**Comparison (10M reads):**

| Format | File Size | Write Time | Read Time |
|--------|-----------|------------|-----------|
| CSV | 1.2 GB | 45s | 30s |
| Parquet | 120 MB | 4s | 2s |

## Hardware Recommendations

| Dataset | CPU Cores | RAM | Storage | Runtime |
|---------|-----------|-----|---------|---------|
| Small (< 1M) | 2-4 | 8 GB | 10 GB | 2-5 min |
| Medium (1-10M) | 4-8 | 16 GB | 50 GB | 5-15 min |
| Large (10-100M) | 8-16 | 32 GB | 200 GB | 15-45 min |
| Very large (100M+) | 16-32 | 64 GB | 1 TB | 1-2 hours |

## Memory Management

### If Running Out of Memory

1. **Switch to streaming mode:**
   ```bash
   metadarkmatter score classify --streaming ...
   ```

2. **Use compressed input:**
   ```bash
   gzip blast_results.tsv
   metadarkmatter score classify --alignment blast_results.tsv.gz ...
   ```

3. **Reduce BLAST hits per read:**
   ```bash
   metadarkmatter blast align --max-targets 100 ...  # Instead of the default 500
   ```

### Memory Usage by Component

| Component | Typical Usage |
|-----------|---------------|
| ANI matrix (500 genomes) | 2 GB |
| BLAST data (vectorized) | 4-8 GB |
| Classification buffer | 2-4 GB |

## Batch Processing Optimization

### Inefficient (repeated ANI loading)

```bash
for f in samples/*.blast.tsv.gz; do
  metadarkmatter score classify --alignment "$f" --ani ani.csv --output "${f%.blast.tsv.gz}.csv"
done
```

### Efficient (single ANI load)

```bash
metadarkmatter score batch \
  --alignment-dir samples/ \
  --ani ani.csv \
  --output-dir results/
```

The batch command loads the ANI matrix once and reuses it for all samples.

## Alignment Performance

### BLAST vs MMseqs2

| Reads | BLAST Time | MMseqs2 Time | Recommendation |
|-------|-----------|--------------|----------------|
| <10K | 10s-10min | Slower | Use BLAST |
| 100K | 30-60 min | 5-10 min | Either tool works |
| 1M+ | 3-6 hours | 15-30 min | MMseqs2 recommended |

### BLAST Alignment Parameters

| Parameter | Sensitive (default) | Balanced | Fast |
|-----------|-----------|----------|------|
| task | blastn | dc-megablast | megablast |
| word_size (task default) | 7 | 11 | 28 |
| evalue | 1e-5 | 1e-5 | 1e-10 |
| max_targets | 500 | 100 | 50 |

`--word-size` defaults to the value for the chosen `--task`; pass it explicitly to override. For closely related references, `--task dc-megablast` has been about 10x faster than `blastn` with similar sensitivity in the 80-100% identity range.

### Parallel BLAST

```bash
# Use multiple threads
metadarkmatter blast align --threads 16 ...

# Or split query and run in parallel
split -n l/4 reads.fasta reads_part_
parallel metadarkmatter blast align --query {} --database db --output {}.blast.tsv ::: reads_part_*
cat reads_part_*.blast.tsv > combined.blast.tsv
```

Concatenating per-part files keeps each read's hits contiguous, so the result can still be classified with `--streaming`.

External tools run without a default timeout, so long BLAST/MMseqs2/ANI runs are not cut off.

## Tips for Large Datasets

1. **Use compressed files:** BLAST output compresses 5-10x
2. **Use Parquet output:** 10x smaller, 10x faster I/O
3. **Use MMseqs2:** 5-100x faster than BLAST for >100K reads
4. **Use streaming mode:** Bounded memory for very large files
5. **Batch process:** Use `score batch` for multiple samples

## Monitoring Progress

All commands show progress indicators:

```
Metadarkmatter ANI-Weighted Classification

Loading ANI matrix...
Classifying reads...  [################............]  60%  ETA: 2:30

Classification Summary
Total reads: 1,234,567
Known species: 876,543 (71.0%)
Novel species: 234,567 (19.0%)
```

Use `--quiet` to suppress output for scripted workflows.
