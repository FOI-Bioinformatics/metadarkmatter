# Performance Guide

Optimizing metadarkmatter for different dataset sizes and hardware configurations.

## Processing Modes

The `score classify` command offers four processing modes optimized for different scenarios:

### Standard Mode (Default)

```bash
metadarkmatter score classify --alignment sample.blast.tsv.gz --ani ani.csv --output out.csv
```

| Metric | Value |
|--------|-------|
| Best for | < 1M alignments |
| RAM usage | 2-4 GB |
| Speed | Baseline |

### Fast Mode

```bash
metadarkmatter score classify --alignment sample.blast.tsv.gz --ani ani.csv --output out.csv --fast
```

| Metric | Value |
|--------|-------|
| Best for | 1-10M alignments |
| RAM usage | 4-8 GB |
| Speed | ~3x faster |

Optimized single-threaded processing path with reduced object creation overhead.

### Parallel Mode (Recommended)

```bash
metadarkmatter score classify --alignment sample.blast.tsv.gz --ani ani.csv --output out.csv --parallel
```

| Metric | Value |
|--------|-------|
| Best for | 10-100M alignments |
| RAM usage | 8-16 GB |
| Speed | ~16x faster |

Uses Polars vectorized operations with automatic parallelization across all CPU cores.

### Streaming Mode

```bash
metadarkmatter score classify --alignment sample.blast.tsv.gz --ani ani.csv --output out.csv --streaming
```

| Metric | Value |
|--------|-------|
| Best for | 100M+ alignments |
| RAM usage | 8-16 GB (constant) |
| Speed | ~10x faster than standard |

Processes data in 5M alignment chunks, writing results incrementally. Memory usage remains constant regardless of file size.

## Mode Selection Guide

```
                    Dataset Size
                         |
         +---------------+---------------+
         |               |               |
      < 1M           1-10M          > 10M
         |               |               |
    Standard          Fast         Memory OK?
                                        |
                                   +----+----+
                                   |         |
                                  Yes        No
                                   |         |
                               Parallel  Streaming
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

## Benchmark Results

Tested on 64-core HPC node with 128 GB RAM:

| Mode | 1M reads | 10M reads | 100M reads |
|------|----------|-----------|------------|
| Standard | 2.1 min | 21.5 min | 215 min |
| Fast | 0.7 min | 7.2 min | 72 min |
| Parallel | 0.4 min | 3.8 min | 38 min |
| Streaming | 0.8 min | 8.1 min | 45 min |

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
   metadarkmatter blast align --max-targets 50 ...  # Instead of 100
   ```

### Memory Usage by Component

| Component | Typical Usage |
|-----------|---------------|
| ANI matrix (500 genomes) | 2 GB |
| BLAST parsing (parallel) | 4-8 GB |
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
  --output-dir results/ \
  --parallel
```

The batch command loads the ANI matrix once and reuses it for all samples.

## BLAST Optimization

### Alignment Parameters

| Parameter | Sensitive | Balanced | Fast |
|-----------|-----------|----------|------|
| word_size | 7 | 11 | 16 |
| evalue | 1e-3 | 1e-5 | 1e-10 |
| max_targets | 100 | 50 | 20 |

### Parallel BLAST

```bash
# Use multiple threads
metadarkmatter blast align --threads 16 ...

# Or split query and run in parallel
split -n l/4 reads.fasta reads_part_
parallel metadarkmatter blast align --query {} --database db --output {}.blast.tsv ::: reads_part_*
cat reads_part_*.blast.tsv > combined.blast.tsv
```

## Tips for Large Datasets

1. **Use compressed files:** BLAST output compresses 5-10x
2. **Use Parquet output:** 10x smaller, 10x faster I/O
3. **Use parallel mode:** Best balance of speed and memory
4. **Process on HPC:** Use high-memory nodes for 100M+ reads
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
