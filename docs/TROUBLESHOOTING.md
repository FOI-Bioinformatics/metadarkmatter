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
Error: Required tool not found: blastn
```

**Solution:**
```bash
# Show which external tools are found and their versions
mdm doctor

# Core tools (skani or fastani for ANI; mmseqs2/diamond/prodigal optional)
conda install -c bioconda blast kraken2 krakentools skani fastani
# bowtie2 + samtools are only needed for 'mdm map' / recruitment plots
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
Error: Invalid input: Unexpected BLAST format: expected 12, 13, or 14 columns, got 6 in results.blast.tsv.gz
```

**Solution:**
```bash
# Verify BLAST output format
zcat results.blast.tsv.gz | head -n 1 | awk -F'\t' '{print NF}'
# Should print 12, 13 (with qlen) or 14 (with qlen and qcov)

# Or validate the file
mdm validate blast --input results.blast.tsv.gz

# Re-run BLAST with the standard 12 columns plus qlen
blastn -query reads.fasta -db blastdb/pangenome \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" \
  -out results.blast.tsv
```

### ANI matrix format incorrect

**Error (examples):**
```
Error: Data processing error: found more fields than defined in 'Schema'
ValueError: ANI matrix is not square: 4 rows, 2 value columns
Error: Invalid input: ANI matrix values appear to be fractions (max off-diagonal value 0.975 <= 1.0); expected percent identity on a 0-100 scale.
```

**Solution:**
```bash
# Check format
head -n 5 ani_matrix.csv

# Or validate the file
mdm validate ani --input ani_matrix.csv

# Required format (square matrix, values 0-100):
# genome,GCF_001,GCF_002,GCF_003
# GCF_001,100.0,95.2,87.3
# GCF_002,95.2,100.0,88.1
```

Long-format (`genome1,genome2,ani`) tables are not accepted, and 0-1
scaled matrices must be multiplied by 100.

### Genome name mismatch

**Symptom:** classification runs, but prints a warning such as:
```
Warning: Low genome coverage (0.0%)
  Only 0 of 3 genomes in BLAST file are in ANI matrix.
  This may indicate mismatched input files.
```
Hits to genomes missing from the matrix fall back to the default 70% ANI,
which inflates placement uncertainty.

**Solution:**
```bash
# Check genome IDs in BLAST output (sseqid format: {accession}|{contig_id})
zcat results.blast.tsv.gz | cut -f2 | cut -d'|' -f1 | sort -u > blast_genomes.txt

# Check genome IDs in ANI matrix
head -n 1 ani_matrix.csv | tr ',' '\n' | tail -n +2 | sort > ani_genomes.txt

# Compare
diff blast_genomes.txt ani_genomes.txt
```

Ensure genome identifiers match exactly between BLAST subject IDs and ANI matrix.
`mdm ani validate --ani ani_matrix.csv --blast results.blast.tsv.gz` reports
the coverage directly. For alignments made outside metadarkmatter, use
`--id-mapping` (see `util generate-mapping`). When `--metadata` has a
`representative` column, only the representatives need to be in the matrix.

---

## Runtime Issues

### Out of memory

**Error:**
```
MemoryError: Unable to allocate array
```

**Solutions:**

1. Use streaming mode (lower `--chunk-size`, default 1,000,000 rows, to cap memory further):
   ```bash
   metadarkmatter score classify --streaming --chunk-size 500000 ...
   ```
   Streaming requires each read's hits to be contiguous (as BLAST/MMseqs2
   write them) and does not support `--id-mapping`/`--genomes`.

2. Reduce BLAST hits per read:
   ```bash
   metadarkmatter blast align --max-targets 50 ...  # Default is 500
   ```

Splitting the alignment file by line count is not recommended: a split
can separate one read's hits across files, which changes its
classification.

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

**Issue:** A read's `taxonomic_call` differs from what the novelty/uncertainty thresholds suggest.

**Cause:** `taxonomic_call` is the Bayesian MAP category; its boundaries
differ from the hard thresholds used by `legacy_call` (for example, with
U = 0 the Known -> Novel Species transition is near N ~ 6.2% rather than
4%). For single-hit reads, the inferred uncertainty only affects
`legacy_call`.

**Solution:**
```bash
# Examine read's metrics
grep "read_00123" classifications.csv

# Check:
# - legacy_call vs taxonomic_call
# - p_* posterior columns and posterior_entropy / confidence_score
# - placement_uncertainty: legacy Known Species requires < 1.5
# - num_ambiguous_hits / ambiguity_scope: many hits across genera indicate a conserved region
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

Note that the `blast align` default `--evalue` is 1e-5 and `--word-size`
defaults per task (7 for blastn, 11 for dc-megablast, 28 for megablast).

---

## Report Issues

### Report not generating

**Symptom:** `report generate` fails or the report is nearly empty.

Empty and all-Off-target classification files produce a report (with
empty sections), so an empty report usually means no reads were
classified upstream.

**Solution:**
```bash
# Check classifications file has data
wc -l classifications.csv

# Check file format
mdm validate classifications --input classifications.csv
```

### Plots not rendering

**Cause:** Reports built with `--report-mode cdn` load Plotly/D3 from a CDN
and need network access to render. The default `--report-mode offline`
inlines Plotly; the phylogenetic tree uses D3, which is inlined only if a
vendored copy exists (`mdm report vendor-d3`), otherwise loaded from a CDN.

**Solution:**
```bash
# Regenerate in offline mode (default) and vendor D3 for the tree panel
mdm report vendor-d3
mdm report generate --classifications classifications.csv --output report.html

# Static image export (kaleido) and matplotlib plots need the viz extra
pip install -e ".[viz]"
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
   https://github.com/FOI-Bioinformatics/metadarkmatter/issues

   Include:
   - Command executed
   - Error message
   - File sizes and read counts
   - System info (OS, RAM, Python version)
