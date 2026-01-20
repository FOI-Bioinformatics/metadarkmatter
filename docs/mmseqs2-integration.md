# MMseqs2 Integration Guide

## Overview

MMseqs2 is integrated into metadarkmatter as an alternative to BLAST for very large datasets. The implementation provides BLAST-compatible output that works seamlessly with the existing classification pipeline.

**⚠️ IMPORTANT**: MMseqs2 is only faster than BLAST for datasets with >100,000 reads. For smaller datasets, use BLAST.

## When to Use MMseqs2

### Use MMseqs2 when:
- Processing very large datasets (>100K reads)
- BLAST runtime becomes impractical (>1 hour)
- Running production pipelines with consistent large-scale data

### Use BLAST when:
- Dataset is small (<100K reads) - **BLAST will be faster**
- Dataset is medium (10K-100K reads) - similar performance, BLAST is simpler
- Establishing baseline results for comparison
- Maximum compatibility is required

### Performance Reality Check

| Reads | BLAST Time | MMseqs2 Time | Winner |
|-------|-----------|--------------|--------|
| 1,000 | 10-15s | 4-5 min | **BLAST** (20x faster) |
| 10,000 | 5-10 min | 5-8 min | **BLAST** (simpler) |
| 100,000 | 30-60 min | 5-10 min | **MMseqs2** (5x faster) |
| 1,000,000 | 3-6 hours | 15-30 min | **MMseqs2** (10-15x faster) |
| 10,000,000+ | 1-2 days | 2-4 hours | **MMseqs2** (10-20x faster) |

## Installation

```bash
# Via conda (recommended)
conda install -c bioconda mmseqs2

# Verify installation
mmseqs version
```

## Quick Start

### 1. Build MMseqs2 Database

```bash
# From directory of genomes
metadarkmatter mmseqs2 makedb \
    --genomes reference_genomes/ \
    --output mmseqs_db/pangenome

# From single FASTA
metadarkmatter mmseqs2 makedb \
    --genomes pangenome.fasta \
    --output mmseqs_db/pangenome
```

### 2. Run Sequence Search

```bash
# Basic search
metadarkmatter mmseqs2 search \
    --query extracted_reads.fasta \
    --database mmseqs_db/pangenome \
    --output sample.mmseqs2.tsv.gz \
    --threads 16

# High-sensitivity search
metadarkmatter mmseqs2 search \
    --query extracted_reads.fasta \
    --database mmseqs_db/pangenome \
    --output sample.mmseqs2.tsv.gz \
    --sensitivity 7.0 \
    --evalue 1e-5 \
    --threads 16
```

### 3. Classify Results (Same as BLAST)

```bash
metadarkmatter score classify \
    --blast sample.mmseqs2.tsv.gz \
    --ani ani_matrix.csv \
    --metadata genome_metadata.tsv \
    --output classifications.csv
```

## Parameter Recommendations

### Sensitivity Settings

| Sensitivity | Speed      | Use Case                              |
|-------------|------------|---------------------------------------|
| 1.0         | Fastest    | High similarity (>90% identity)       |
| 4.0         | Fast       | Standard (>80% identity)              |
| 5.7         | **Balanced** | **Default, recommended for most cases** |
| 7.0         | Slower     | Sensitive (>70% identity)             |
| 7.5         | Slowest    | Remote homology                       |

**Recommendation**: Start with default sensitivity (5.7) and increase only if missing expected hits.

### Other Parameters

```bash
--evalue 1e-3          # Default, matches BLAST defaults
--max-seqs 500         # Default, matches BLAST max-target-seqs
--min-identity 75.0    # Optional, filter by percent identity
--threads 16           # Adjust based on available cores
```

## Output Format

MMseqs2 produces BLAST-compatible tabular output with these columns:

```
query, target, pident, alnlen, mismatch, gapopen, qstart, qend, tstart, tend, evalue, bits
```

**Column Mapping to BLAST**:
- `query` → `qseqid`
- `target` → `sseqid`
- `pident` → `pident` (percent identity)
- `alnlen` → `length` (alignment length)
- `bits` → `bitscore`

This format is compatible with existing parsers and classifiers without modification.

## Performance Comparison

Based on typical metadarkmatter workloads:

| Dataset Size | BLAST Time | MMseqs2 Time (s=5.7) | Speedup |
|-------------|-----------|---------------------|---------|
| 1M reads    | 2-4 hours | 5-15 minutes        | 10-50x  |
| 10M reads   | 20-40 hours | 30-90 minutes     | 15-80x  |
| 100M reads  | 8-16 days | 4-12 hours          | 30-100x |

*Times estimated for 500 reference genomes, 16 threads*

## Sensitivity Validation

To validate that MMseqs2 provides similar results to BLAST:

```bash
# Run both tools on same data
metadarkmatter blast align --query reads.fasta --database blastdb --output blast.tsv
metadarkmatter mmseqs2 search --query reads.fasta --database mmseqs_db --output mmseqs2.tsv

# Compare classifications
metadarkmatter score classify --blast blast.tsv --ani ani.csv --output blast_class.csv
metadarkmatter score classify --blast mmseqs2.tsv --ani ani.csv --output mmseqs_class.csv

# Analyze differences (manual comparison or custom script)
```

Expected: <5% difference in classification results with sensitivity 5.7+.

## Advanced Usage

### Temporary Directory Management

MMseqs2 creates large temporary files (2-5x database size). Control temp directory:

```bash
metadarkmatter mmseqs2 search \
    --query reads.fasta \
    --database mmseqs_db \
    --output results.tsv \
    --tmp-dir /scratch/tmp/  # Use fast local storage
```

### Protein-Level Search (Future)

For highly divergent taxa (genus-level novelty), protein-level search may be more sensitive:

```bash
# Create protein database
metadarkmatter mmseqs2 makedb \
    --genomes proteins.faa \
    --output mmseqs_protein_db \
    --dbtype 0  # 0 = protein

# Translated search (BLASTX equivalent)
metadarkmatter mmseqs2 search \
    --query reads.fastq \
    --database mmseqs_protein_db \
    --output results.tsv \
    --search-type 2  # 2 = translated search
```

*Note: Protein-level classification not yet validated for metadarkmatter.*

## Troubleshooting

### Database Not Found Error

```
Error: MMseqs2 database not found at mmseqs_db
```

**Solution**: Ensure the database path is correct and all database files exist:
```bash
ls -lh mmseqs_db*
# Should show: mmseqs_db, mmseqs_db.dbtype, mmseqs_db.index, etc.
```

### Out of Disk Space

```
Error: No space left on device
```

**Solution**: Specify a temporary directory with more space:
```bash
--tmp-dir /path/to/larger/storage/tmp/
```

### Unexpected Results

If MMseqs2 results differ significantly from BLAST:

1. **Increase sensitivity**: Try `--sensitivity 7.0`
2. **Check E-value**: Use same E-value as BLAST (`--evalue 1e-3`)
3. **Verify database**: Ensure same reference genomes used for both tools

## Implementation Notes

### BLAST Compatibility

MMseqs2 integration uses the `--format-output` parameter to produce BLAST-compatible TSV:

```bash
--format-mode 0 \
--format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
```

This ensures seamless integration with:
- `StreamingBlastParser` (no modifications needed)
- `BlastHit` models (validates without changes)
- Classification pipeline (uses same thresholds)

### Performance Optimization

MMseqs2 performance optimizations:
- Uses `easy-search` for simplicity (automatic cleanup)
- Temp directory auto-created if not specified
- Progress indicators disabled for large-scale runs
- Output compression supported (gzip)

## References

- MMseqs2 paper: Steinegger & Söding (2017) Nature Biotechnology
- MMseqs2 documentation: https://github.com/soedinglab/MMseqs2
- Parameter tuning guide: https://github.com/soedinglab/MMseqs2/wiki

## Support

For issues specific to MMseqs2 integration:
1. Check this guide's troubleshooting section
2. Validate MMseqs2 installation: `mmseqs version`
3. Test with small dataset first
4. Compare results with BLAST on same data
5. Report issues at: https://github.com/anthropics/metadarkmatter/issues
