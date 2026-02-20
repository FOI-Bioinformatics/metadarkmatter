# Troubleshooting Guide

Solutions for common issues when using metadarkmatter.

## Installation Issues

### Python version too old

**Error:**
```
ERROR: metadarkmatter requires Python '>=3.11'
```

**Solution:**
```bash
# Using conda
conda create -n metadarkmatter python=3.11
conda activate metadarkmatter
pip install -e .

# Using pyenv
pyenv install 3.11.0
pyenv local 3.11.0
pip install -e .
```

### Dependency conflicts

**Error:**
```
ERROR: Cannot install metadarkmatter because these package versions conflict
```

**Solution:**
```bash
# Create fresh virtual environment
python3.11 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -e .
```

### Missing external tools

**Error:**
```
Error: Required tool not found: blastn (BLAST+)
```

**Solution:**
```bash
conda install -c bioconda blast kraken2 krakentools bowtie2 samtools fastani
```

---

## Input File Issues

### BLAST file not found

**Error:**
```
Error: BLAST file not found: results.blast.tsv.gz
```

**Solution:**
```bash
# Check file exists
ls -lh results.blast.tsv.gz

# Use absolute path if needed
metadarkmatter score classify --alignment /absolute/path/to/results.blast.tsv.gz ...
```

### Invalid BLAST format

**Error:**
```
ValueError: Expected 12 fields in BLAST line, got 6
```

**Solution:**
```bash
# Verify BLAST output format
zcat results.blast.tsv.gz | head -n 1
# Should have 12 tab-separated columns

# Re-run BLAST with correct format
blastn -query reads.fasta -db genomes.fasta -outfmt 6 -out results.blast.tsv
```

### ANI matrix format incorrect

**Error:**
```
Error loading ANI matrix: CSV parse error
```

**Solution:**
```bash
# Check format
head -n 5 ani_matrix.csv

# Required format:
# genome,GCF_001,GCF_002,GCF_003
# GCF_001,100.0,95.2,87.3
# GCF_002,95.2,100.0,88.1
```

### Genome name mismatch

**Error:**
```
KeyError: 'GCF_000195955.1' not found in ANI matrix
```

**Solution:**
```bash
# Check genome IDs in BLAST output
zcat results.blast.tsv.gz | cut -f2 | sed 's/_.*//g' | sort -u > blast_genomes.txt

# Check genome IDs in ANI matrix
head -n 1 ani_matrix.csv | tr ',' '\n' | tail -n +2 > ani_genomes.txt

# Compare
diff blast_genomes.txt ani_genomes.txt
```

Ensure genome identifiers match exactly between BLAST subject IDs and ANI matrix.

---

## Runtime Issues

### Out of memory

**Error:**
```
MemoryError: Unable to allocate array
```

**Solutions:**

1. Use streaming mode:
   ```bash
   metadarkmatter score classify --streaming ...
   ```

2. Reduce BLAST hits per read:
   ```bash
   metadarkmatter blast align --max-targets 50 ...  # Instead of 100
   ```

3. Split BLAST file:
   ```bash
   zcat huge_file.blast.tsv.gz | split -l 10000000 - chunk_
   for chunk in chunk_*; do
     metadarkmatter score classify --alignment "$chunk" --ani ani.csv --output "${chunk}.csv"
   done
   cat chunk_*.csv > combined.csv
   ```

### Processing too slow

**Solutions:**

1. Use Parquet output:
   ```bash
   metadarkmatter score classify --format parquet ...
   ```

### BLAST database not found

**Error:**
```
Error: BLAST database not found at blastdb/pangenome
```

**Solution:**
```bash
# Verify database files exist
ls blastdb/pangenome.{nhr,nin,nsq}

# Rebuild if missing
metadarkmatter blast makedb --genomes genomes/ --output blastdb/pangenome
```

---

## Output Issues

### Empty output file

**Symptoms:** Classification file has only header, no data rows.

**Causes:**
1. No reads matched criteria
2. All reads filtered out by bitscore threshold

**Solution:**
```bash
# Count BLAST alignments
zcat results.blast.tsv.gz | wc -l

# Check bitscores
zcat results.blast.tsv.gz | cut -f12 | sort -n | tail -n 10

# Try lower bitscore threshold
metadarkmatter score classify --bitscore-threshold 90.0 ...
```

### All reads classified as "Conserved Region"

**Symptoms:** > 80% reads as conserved regions, very low novel diversity.

**Causes:**
1. BLAST max_target_seqs too high
2. Reference genomes are highly similar
3. Reads from conserved regions (rRNA)

**Solutions:**

1. Reduce BLAST targets:
   ```bash
   metadarkmatter blast align --max-targets 50 ...
   ```

2. Filter rRNA reads before BLAST:
   ```bash
   # Use SortMeRNA or similar
   sortmerna --ref rRNA_databases/silva-bac-16s-id90.fasta \
             --reads family_reads.fasta \
             --aligned rRNA_reads \
             --other non_rRNA_reads
   ```

### Unexpected classification

**Issue:** Read with 98% identity classified as "Novel Species" instead of "Known Species".

**Cause:** High placement uncertainty due to multiple similar hits.

**Solution:**
```bash
# Examine read's metrics
grep "read_00123" classifications.csv

# Check:
# - placement_uncertainty: Should be < 0.5 for Known Species
# - num_ambiguous_hits: High values indicate conserved region
```

---

## Data Quality Issues

### High novel diversity in pure culture

**Symptoms:** > 20% novel diversity in pure culture sample.

**Causes:**
1. Reference genome not in database
2. Sample contamination
3. Sequencing errors

**Solutions:**

1. Verify reference genome:
   ```bash
   grep "GCF_XXXXXX" ani_matrix.csv
   ```

2. Quality filter reads:
   ```bash
   fastp -i reads_R1.fastq -I reads_R2.fastq \
         -o filtered_R1.fastq -O filtered_R2.fastq
   ```

### Low novel diversity in environmental sample

**Symptoms:** < 5% novel diversity in soil/water sample.

**Causes:**
1. Reference database has excellent coverage (rare)
2. BLAST parameters too stringent

**Solution:**
```bash
# Verify BLAST parameters
# Correct (sensitive):
metadarkmatter blast align --word-size 7 --evalue 1e-3 ...

# Too stringent (may miss divergent hits):
metadarkmatter blast align --word-size 11 --evalue 1e-10 ...
```

---

## Report Issues

### Report not generating

**Error:**
```
Error: No classification data found
```

**Solution:**
```bash
# Check classifications file has data
wc -l classifications.csv

# Check file format
head classifications.csv
```

### Plots not rendering

**Cause:** Missing visualization dependencies.

**Solution:**
```bash
pip install -e ".[viz]"
# Includes: kaleido, matplotlib
```

---

## Getting Help

If issues persist:

1. **Enable verbose output:**
   ```bash
   metadarkmatter score classify --verbose ...
   ```

2. **Check with minimal example:**
   ```bash
   zcat full_file.blast.tsv.gz | head -n 1000 > test.blast.tsv
   metadarkmatter score classify --alignment test.blast.tsv ...
   ```

3. **Report issue on GitHub:**
   https://github.com/metadarkmatter/metadarkmatter/issues

   Include:
   - Command executed
   - Error message
   - File sizes and read counts
   - System info (OS, RAM, Python version)
