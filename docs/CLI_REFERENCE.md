# Metadarkmatter CLI Reference

Complete command-line interface reference for the metadarkmatter tool.

## Installation

```bash
pip install metadarkmatter
```

## Quick Start

```bash
# Classify alignment results using ANI-weighted placement
metadarkmatter score classify \
    --alignment results.blast.tsv.gz \
    --ani ani_matrix.csv \
    --output classifications.csv

# Batch process multiple files
metadarkmatter score batch \
    --alignment-dir ./blast_results/ \
    --ani ani_matrix.csv \
    --output-dir ./classifications/
```

---

## Global Commands

### `metadarkmatter --version`

Display the installed version.

```bash
metadarkmatter --version
# Output: metadarkmatter 0.1.0
```

### `metadarkmatter --help`

Show available commands and options.

---

## Blast Subcommand

The `blast` subcommand wraps BLAST+ tools for database creation and competitive alignment.

### `metadarkmatter blast makedb`

Create a BLAST nucleotide database from reference genomes.

#### Required Options

| Option | Short | Type | Description |
|--------|-------|------|-------------|
| `--genomes` | `-g` | PATH | Path to genome FASTA file or directory of genomes |
| `--output` | `-o` | PATH | Output database prefix |

#### Optional Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--title` | `-t` | TEXT | output name | Database title |
| `--genome-pattern` | `-p` | TEXT | *.fna | Glob pattern for genome files (when --genomes is directory) |
| `--skip-if-exists` | | FLAG | False | Skip if database already exists |
| `--verbose` | `-v` | FLAG | False | Enable verbose output |
| `--quiet` | `-q` | FLAG | False | Suppress progress output |
| `--dry-run` | | FLAG | False | Show commands without executing |

#### Examples

**From directory of genomes**
```bash
metadarkmatter blast makedb \
    --genomes reference_genomes/ \
    --output blastdb/pangenome
```

**From single FASTA**
```bash
metadarkmatter blast makedb \
    --genomes pangenome.fasta \
    --output blastdb/pangenome
```

---

### `metadarkmatter blast align`

Run BLASTN competitive alignment for ANI-weighted placement.

#### Required Options

| Option | Short | Type | Description |
|--------|-------|------|-------------|
| `--query` | `-q` | PATH | Query reads file (FASTA/FASTQ) |
| `--database` | `-d` | PATH | BLAST database prefix |
| `--output` | `-o` | PATH | Output file for BLAST results |

#### Optional Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--word-size` | `-w` | INT | 7 | Word size for seeds (smaller = more sensitive) |
| `--evalue` | `-e` | FLOAT | 1e-3 | E-value threshold |
| `--max-targets` | `-m` | INT | 100 | Maximum target sequences per query |
| `--min-identity` | | FLOAT | None | Minimum percent identity filter |
| `--task` | | TEXT | blastn | BLAST task (blastn, megablast, dc-megablast) |
| `--threads` | `-p` | INT | 4 | Number of threads |
| `--compress/--no-compress` | | FLAG | True | Compress output with gzip |
| `--skip-if-exists` | | FLAG | False | Skip if output file exists |
| `--verbose` | `-v` | FLAG | False | Enable verbose output |
| `--quiet` | `-q` | FLAG | False | Suppress progress output |
| `--dry-run` | | FLAG | False | Show commands without executing |

#### Examples

**Basic alignment**
```bash
metadarkmatter blast align \
    --query extracted_reads.fasta \
    --database blastdb/pangenome \
    --output sample.blast.tsv.gz
```

**High-sensitivity alignment**
```bash
metadarkmatter blast align \
    --query extracted_reads.fasta \
    --database blastdb/pangenome \
    --output sample.blast.tsv.gz \
    --word-size 7 \
    --evalue 1e-3 \
    --max-targets 100 \
    --threads 16
```

**Output format**

The output is tabular format (outfmt 6) with 12 columns:
```
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
```

---

## BLASTX Subcommand

The `blastx` subcommand wraps Diamond for protein-level read classification using translated nucleotide queries against protein databases. This is useful for detecting divergent homology that nucleotide-based alignment may miss.

### `metadarkmatter blastx makedb`

Create a Diamond protein database from reference genome protein files (.faa).

#### Required Options

| Option | Short | Type | Description |
|--------|-------|------|-------------|
| `--proteins` | `-p` | PATH | Path to protein FASTA file or directory of .faa files |
| `--output` | `-o` | PATH | Output database prefix |

#### Optional Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--protein-pattern` | | TEXT | *.faa | Glob pattern for protein files (when --proteins is directory) |
| `--threads` | `-t` | INT | 4 | Number of threads for database building |
| `--skip-if-exists` | | FLAG | False | Skip if database already exists |
| `--verbose` | `-v` | FLAG | False | Enable verbose output |
| `--quiet` | `-q` | FLAG | False | Suppress progress output |
| `--dry-run` | | FLAG | False | Show commands without executing |

#### Examples

**From directory of protein files**
```bash
metadarkmatter blastx makedb \
    --proteins reference_proteins/ \
    --output blastdb/panproteome \
    --threads 16
```

**From single protein FASTA**
```bash
metadarkmatter blastx makedb \
    --proteins all_proteins.faa \
    --output blastdb/panproteome
```

---

### `metadarkmatter blastx align`

Run Diamond BLASTX alignment (DNA query vs protein database) for protein-level read classification.

#### Required Options

| Option | Short | Type | Description |
|--------|-------|------|-------------|
| `--query` | `-q` | PATH | Query reads file (FASTA/FASTQ, optionally gzipped) |
| `--database` | `-d` | PATH | Diamond protein database prefix |
| `--output` | `-o` | PATH | Output file for BLASTX results |

#### Optional Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--evalue` | `-e` | FLOAT | 1e-5 | E-value threshold |
| `--min-identity` | | FLOAT | 30.0 | Minimum percent identity filter |
| `--max-targets` | `-m` | INT | 500 | Maximum target sequences per query |
| `--threads` | `-p` | INT | 4 | Number of threads |
| `--sensitive/--fast` | | FLAG | True | Use sensitive mode (recommended for divergent sequences) |
| `--compress/--no-compress` | | FLAG | True | Compress output with gzip |
| `--skip-if-exists` | | FLAG | False | Skip if output file exists |
| `--verbose` | `-v` | FLAG | False | Enable verbose output |
| `--quiet` | `-q` | FLAG | False | Suppress progress output |
| `--dry-run` | | FLAG | False | Show commands without executing |

#### Examples

**Basic BLASTX alignment**
```bash
metadarkmatter blastx align \
    --query extracted_reads.fastq.gz \
    --database blastdb/panproteome \
    --output sample.blastx.tsv.gz
```

**High-sensitivity alignment for divergent sequences**
```bash
metadarkmatter blastx align \
    --query extracted_reads.fastq.gz \
    --database blastdb/panproteome \
    --output sample.blastx.tsv.gz \
    --evalue 1e-3 \
    --min-identity 30.0 \
    --max-targets 500 \
    --sensitive \
    --threads 16
```

**Output format**

The output is 12-column tabular format matching BLASTN:
```
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
```

**Note:** BLASTX translates DNA queries in all six reading frames before alignment against the protein database. Use `--alignment-mode protein` when running `score classify` on BLASTX output to use protein-calibrated thresholds.

---

## Kraken2 Subcommand

The `kraken2` subcommand provides Kraken2 taxonomic classification and read extraction for target taxa.

### `metadarkmatter kraken2 classify`

Run Kraken2 taxonomic classification on metagenomic reads.

#### Required Options

| Option | Short | Type | Description |
|--------|-------|------|-------------|
| `--reads-1` | `-1` | PATH | Forward reads file (FASTQ, optionally gzipped) |
| `--kraken-db` | `-k` | PATH | Path to Kraken2 database directory |
| `--output` | `-o` | PATH | Output directory for Kraken2 files |

#### Optional Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--reads-2` | `-2` | PATH | None | Reverse reads file for paired-end data |
| `--sample-name` | `-n` | TEXT | auto | Sample name (default: derived from reads filename) |
| `--confidence` | | FLOAT | 0.0 | Kraken2 confidence threshold (0-1) |
| `--minimum-hit-groups` | | INT | 2 | Minimum hit groups for classification |
| `--skip-if-exists` | | FLAG | False | Skip if output files already exist |
| `--threads` | `-p` | INT | 4 | Number of threads for Kraken2 |
| `--verbose` | `-v` | FLAG | False | Enable verbose output |
| `--quiet` | `-q` | FLAG | False | Suppress progress output |
| `--dry-run` | | FLAG | False | Show commands without executing |

#### Examples

**Basic Classification**
```bash
metadarkmatter kraken2 classify \
    --reads-1 sample_R1.fastq.gz \
    --reads-2 sample_R2.fastq.gz \
    --kraken-db /path/to/kraken_db \
    --output kraken_output/
```

**With Custom Settings**
```bash
metadarkmatter kraken2 classify \
    --reads-1 sample_R1.fastq.gz \
    --kraken-db /path/to/kraken_db \
    --output kraken_output/ \
    --confidence 0.1 \
    --threads 16
```

**Output files:**
- `{sample}.kraken` - Per-read classification output
- `{sample}.kreport` - Hierarchical taxonomic report

---

### `metadarkmatter kraken2 extract`

Extract reads belonging to a specific taxon from Kraken2 output.

#### Required Options

| Option | Short | Type | Description |
|--------|-------|------|-------------|
| `--kraken-output` | `-k` | PATH | Kraken2 output file (.kraken) from 'kraken2 classify' |
| `--reads-1` | `-1` | PATH | Original forward reads file (FASTQ) |
| `--taxid` | `-t` | INT | Target taxonomic ID to extract |
| `--output` | `-o` | PATH | Output directory for extracted reads |

#### Optional Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--reads-2` | `-2` | PATH | None | Original reverse reads file for paired-end data |
| `--sample-name` | `-n` | TEXT | auto | Sample name for output files |
| `--include-children/--exact-match` | | FLAG | True | Include child taxa (default) or exact match only |
| `--skip-if-exists` | | FLAG | False | Skip if output files already exist |
| `--verbose` | `-v` | FLAG | False | Enable verbose output |
| `--quiet` | `-q` | FLAG | False | Suppress progress output |
| `--dry-run` | | FLAG | False | Show commands without executing |

#### Examples

**Extract Family Reads (include child taxa)**
```bash
metadarkmatter kraken2 extract \
    --kraken-output kraken_output/sample.kraken \
    --reads-1 sample_R1.fastq.gz \
    --reads-2 sample_R2.fastq.gz \
    --taxid 262 \
    --output extracted/
```

**Exact Match Only**
```bash
metadarkmatter kraken2 extract \
    --kraken-output kraken_output/sample.kraken \
    --reads-1 sample_R1.fastq.gz \
    --taxid 262 \
    --output extracted/ \
    --exact-match
```

**Output files:**
- `{sample}_taxid{taxid}_R1.fastq` - Extracted forward reads
- `{sample}_taxid{taxid}_R2.fastq` - Extracted reverse reads (if paired)

---

## Score Subcommand

The `score` subcommand performs ANI-weighted placement classification on BLAST results.

### `metadarkmatter score classify`

Classify metagenomic reads from an alignment file (BLAST or MMseqs2).

#### Required Options

| Option | Short | Type | Description |
|--------|-------|------|-------------|
| `--alignment` | `-b` | PATH | Path to alignment file (.tsv or .tsv.gz) - accepts BLAST or MMseqs2 tabular format |
| `--ani` | `-a` | PATH | Path to ANI matrix file (CSV or TSV) |
| `--output` | `-o` | PATH | Output path for classification results |

#### Optional Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--summary` | `-s` | PATH | None | Output path for summary statistics (JSON) |
| `--alignment-mode` | | TEXT | nucleotide | Alignment type: 'nucleotide' for BLASTN, 'protein' for BLASTX |
| `--bitscore-threshold` | | FLOAT | 95.0 | Percentage of top bitscore for ambiguous hits (0-100) |
| `--format` | `-f` | TEXT | csv | Output format: 'csv' or 'parquet' |
| `--verbose` | `-v` | FLAG | False | Enable verbose output |
| `--quiet` | `-q` | FLAG | False | Suppress progress output (for scripting) |
| `--dry-run` | | FLAG | False | Validate inputs without processing |

#### Coverage Weighting Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--coverage-weight-mode` | TEXT | none | Coverage weighting mode: 'none', 'linear', 'log', 'sigmoid' |
| `--coverage-weight-strength` | FLOAT | 0.5 | Coverage weight strength (0.0-1.0) |

**Coverage Weighting Modes:**
- `none` (default): Use raw bitscore for hit selection
- `linear`: Gradual penalty for low coverage
- `log`: Rewards any coverage, diminishing returns
- `sigmoid`: Sharp threshold around 60% coverage

#### Alignment Mode

The `--alignment-mode` option selects threshold calibration:

| Mode | Source | Thresholds |
|------|--------|------------|
| `nucleotide` | BLASTN output | Known: N<5%, Novel Species: 5-20%, Novel Genus: 20-25% |
| `protein` | BLASTX output | Known: N<10%, Novel Species: 10-25%, Novel Genus: 25-40% |

Protein thresholds are wider because amino acid sequences diverge more slowly than nucleotide sequences.

#### Processing Mode Options (Mutually Exclusive)

| Option | Description | Use Case |
|--------|-------------|----------|
| `--fast` | Single-threaded optimized (~3x faster) | Small to medium files |
| `--parallel` | Polars vectorized (~16x faster) | Large files, recommended |
| `--streaming` | Memory-efficient for 100M+ alignments | Very large files |

#### Examples

**Basic Classification**
```bash
metadarkmatter score classify \
    --alignment sample.blast.tsv.gz \
    --ani genomes.ani.csv \
    --output sample_classifications.csv
```

**With Summary and Parallel Processing**
```bash
metadarkmatter score classify \
    --alignment sample.blast.tsv.gz \
    --ani genomes.ani.csv \
    --output sample_classifications.parquet \
    --summary sample_summary.json \
    --format parquet \
    --parallel
```

**Coverage-Weighted Classification**
```bash
metadarkmatter score classify \
    --alignment sample.blast.tsv.gz \
    --ani genomes.ani.csv \
    --output sample_classifications.csv \
    --coverage-weight-mode linear \
    --coverage-weight-strength 0.5
```

**Dry Run to Validate Inputs**
```bash
metadarkmatter score classify \
    --alignment sample.blast.tsv.gz \
    --ani genomes.ani.csv \
    --output sample_classifications.csv \
    --dry-run
```

**Quiet Mode for Scripts**
```bash
metadarkmatter score classify \
    --alignment sample.blast.tsv.gz \
    --ani genomes.ani.csv \
    --output sample_classifications.csv \
    --quiet
```

**Protein Mode (for BLASTX output)**
```bash
metadarkmatter score classify \
    --alignment sample.blastx.tsv.gz \
    --ani genomes.ani.csv \
    --output sample_classifications.csv \
    --alignment-mode protein \
    --parallel
```

---

### `metadarkmatter score batch`

Batch classify multiple alignment files in a directory.

#### Required Options

| Option | Short | Type | Description |
|--------|-------|------|-------------|
| `--alignment-dir` | `-b` | PATH | Directory containing alignment files (BLAST or MMseqs2) |
| `--ani` | `-a` | PATH | Path to ANI matrix file (CSV or TSV) |
| `--output-dir` | `-o` | PATH | Output directory for classification results |

#### Optional Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--pattern` | `-p` | TEXT | *.blast.tsv.gz | Glob pattern for alignment files |
| `--bitscore-threshold` | | FLOAT | 95.0 | Percentage of top bitscore for ambiguous hits |
| `--format` | `-f` | TEXT | csv | Output format: 'csv' or 'parquet' |
| `--workers` | `-w` | INT | CPU-1 | Number of worker processes |
| `--verbose` | `-v` | FLAG | False | Enable verbose output |
| `--coverage-weight-mode` | | TEXT | none | Coverage weighting mode |
| `--coverage-weight-strength` | | FLOAT | 0.5 | Coverage weight strength (0.0-1.0) |

#### Processing Mode Options

| Option | Description |
|--------|-------------|
| `--fast` | Single-threaded optimized per file |
| `--parallel` | Multi-core vectorized per file |

#### Examples

**Basic Batch Processing**
```bash
metadarkmatter score batch \
    --alignment-dir ./blast_results/ \
    --ani genomes.ani.csv \
    --output-dir ./classifications/
```

**Custom Pattern with Parquet Output**
```bash
metadarkmatter score batch \
    --alignment-dir ./blast_results/ \
    --ani genomes.ani.csv \
    --output-dir ./classifications/ \
    --pattern "*.tsv" \
    --format parquet \
    --parallel
```

**With Coverage Weighting**
```bash
metadarkmatter score batch \
    --alignment-dir ./blast_results/ \
    --ani genomes.ani.csv \
    --output-dir ./classifications/ \
    --coverage-weight-mode linear \
    --parallel
```

---

## Input File Formats

### Alignment Tabular Format

Both BLAST and MMseqs2 produce compatible tabular output. The standard format includes 13 columns:

```
qseqid  sseqid  pident  length  mismatch  gapopen  qstart  qend  sstart  send  evalue  bitscore  qlen
```

| Column | Description | Example |
|--------|-------------|---------|
| qseqid | Query sequence ID | read_001 |
| sseqid | Subject sequence ID | GCF_000123456.1\|contig1 |
| pident | Percent identity | 98.5 |
| length | Alignment length | 150 |
| mismatch | Number of mismatches | 2 |
| gapopen | Number of gap openings | 0 |
| qstart | Query start position | 1 |
| qend | Query end position | 150 |
| sstart | Subject start position | 1000 |
| send | Subject end position | 1150 |
| evalue | E-value | 1e-80 |
| bitscore | Bit score | 280 |
| qlen | Query sequence length | 150 |

**Note:** The `qlen` column (column 13) is used for coverage-weighted hit selection. Older 12-column output (without qlen) is still supported - coverage is estimated from qend in this case.

### ANI Matrix Format

CSV or TSV with genome names as first column and header row:

```csv
genome,GCF_000123456.1,GCF_000789012.1,GCA_000111222.1
GCF_000123456.1,100.0,95.5,82.3
GCF_000789012.1,95.5,100.0,83.1
GCA_000111222.1,82.3,83.1,100.0
```

---

## Output File Formats

### Classification Output (CSV/Parquet)

| Column | Type | Description |
|--------|------|-------------|
| read_id | string | Query sequence identifier |
| best_match_genome | string | Top-scoring genome |
| novelty_index | float | 100 - percent_identity (0-100) |
| placement_uncertainty | float | ANI-based ambiguity (0-100) |
| taxonomic_call | string | Classification result |
| is_novel | boolean | True if Novel Species or Novel Genus |

**Taxonomic Call Values:**
- `Known Species`: High identity, unambiguous placement
- `Novel Species`: Moderate divergence (5-15% novelty)
- `Novel Genus`: High divergence (15-25% novelty)
- `Conserved Region`: Ambiguous placement across genomes

### Summary Output (JSON)

```json
{
  "total_reads": 10000,
  "known_species": 6500,
  "novel_species": 2000,
  "novel_genus": 500,
  "conserved_regions": 1000,
  "known_species_pct": 65.0,
  "novel_species_pct": 20.0,
  "novel_genus_pct": 5.0,
  "novel_diversity_pct": 25.0,
  "mean_novelty_index": 3.5,
  "mean_placement_uncertainty": 1.2,
  "genome_hit_counts": {
    "GCF_000123456.1": 3500,
    "GCF_000789012.1": 2800
  }
}
```

---

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Error (invalid input, processing failure) |
| 2 | Missing required argument |

---

## Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `METADARKMATTER_CONFIG` | Path to YAML config file | None |
| `POLARS_MAX_THREADS` | Maximum threads for Polars | CPU count |

---

## Performance Guide

### Classifier Selection

| Scenario | Recommended Mode | Command |
|----------|-----------------|---------|
| < 1M alignments | Standard | (default) |
| 1-10M alignments | Fast | `--fast` |
| 10-100M alignments | Parallel | `--parallel` |
| > 100M alignments | Streaming | `--streaming` |

### Memory Usage

| Mode | Memory per 10M Alignments |
|------|--------------------------|
| Standard | ~4 GB |
| Fast | ~3 GB |
| Parallel | ~2 GB |
| Streaming | ~500 MB (bounded) |

---

## Troubleshooting

### Common Errors

**"Processing modes are mutually exclusive"**
- Use only one of `--fast`, `--parallel`, or `--streaming`

**"No BLAST files found matching pattern"**
- Check the `--pattern` option matches your files
- Verify the `--blast-dir` path is correct

**"Low genome coverage warning"**
- Many genomes in BLAST file are not in ANI matrix
- Consider regenerating ANI matrix with all reference genomes

**"Invalid format"**
- Use only 'csv' or 'parquet' with `--format`

### Getting Help

```bash
# Command help
metadarkmatter score classify --help
metadarkmatter score batch --help

# Validate inputs without processing
metadarkmatter score classify --blast file.tsv --ani matrix.csv --output out.csv --dry-run
```

---

## Report Subcommand

The `report` subcommand generates unified HTML reports with interactive visualizations.

### `metadarkmatter report generate`

Generate a unified HTML report from classification results.

#### Required Options

| Option | Short | Type | Description |
|--------|-------|------|-------------|
| `--classifications` | `-c` | PATH | Classification results file (CSV or Parquet) |
| `--output` | `-o` | PATH | Output HTML report path |

#### Optional Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--sample-name` | `-n` | TEXT | filename | Sample name for report header |
| `--title` | `-t` | TEXT | "Metadarkmatter Classification Report" | Report title |
| `--ani` | `-a` | PATH | None | ANI matrix file for heatmap visualization |
| `--recruitment` | `-r` | PATH | None | Recruitment data file for recruitment plots |
| `--bam` | `-b` | PATH | None | BAM file for generating recruitment data on-the-fly |
| `--theme` | | TEXT | light | Color theme: 'light' or 'dark' |
| `--max-points` | | INT | 50000 | Maximum points in scatter plots (for performance) |
| `--max-table-rows` | | INT | 10000 | Maximum rows in data table (for file size) |
| `--verbose` | `-v` | FLAG | False | Enable verbose output |
| `--quiet` | `-q` | FLAG | False | Suppress progress output |

#### Report Sections

The generated HTML report includes 6 tabbed sections:

| Section | Contents |
|---------|----------|
| **Overview** | Metric cards, classification donut chart, bar chart |
| **Distributions** | Novelty histogram, uncertainty histogram, 2D scatter plot |
| **Recruitment** | Read recruitment plots (if BAM/data provided) |
| **Genomes** | Per-genome read count bar chart, identity box plots |
| **ANI Matrix** | Genome-genome ANI heatmap (if provided) |
| **Data** | Interactive sortable/filterable classification table with export |

#### Examples

**Basic Report**
```bash
metadarkmatter report generate \
    --classifications results.csv \
    --output report.html \
    --sample-name "Sample_001"
```

**Full Report with ANI and Recruitment**
```bash
metadarkmatter report generate \
    --classifications results.csv \
    --output report.html \
    --ani ani_matrix.csv \
    --bam mapped.bam \
    --sample-name "Environmental Sample"
```

**Dark Theme Report**
```bash
metadarkmatter report generate \
    --classifications results.csv \
    --output report.html \
    --theme dark
```

**Performance Optimization for Large Datasets**
```bash
metadarkmatter report generate \
    --classifications results.parquet \
    --output report.html \
    --max-points 30000 \
    --max-table-rows 5000
```

---

### `metadarkmatter report multi`

Generate a multi-sample comparison report.

#### Required Options

| Option | Short | Type | Description |
|--------|-------|------|-------------|
| `--input-dir` | `-i` | PATH | Directory containing classification result files |
| `--output` | `-o` | PATH | Output HTML report path |

#### Optional Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--pattern` | `-p` | TEXT | *classifications*.csv | Glob pattern to match classification files |
| `--title` | `-t` | TEXT | "Multi-Sample Comparison Report" | Report title |
| `--theme` | | TEXT | light | Color theme: 'light' or 'dark' |
| `--verbose` | `-v` | FLAG | False | Enable verbose output |
| `--quiet` | `-q` | FLAG | False | Suppress progress output |

#### Multi-Sample Visualizations

| Visualization | Description |
|--------------|-------------|
| **Summary Table** | Statistics for each sample (reads, classifications, novelty) |
| **Stacked Bar Chart** | Classification proportions by sample |
| **Grouped Bar Chart** | Classification comparison across samples |
| **Heatmap** | Sample-category classification proportions |
| **Diversity Scatter** | Mean novelty vs uncertainty per sample |
| **Box Plots** | Novelty index distribution by sample |
| **Trend Line** | Novel diversity percentage across samples |

#### Examples

**Basic Multi-Sample Report**
```bash
metadarkmatter report multi \
    --input-dir ./results/ \
    --output comparison.html
```

**Custom File Pattern**
```bash
metadarkmatter report multi \
    --input-dir ./results/ \
    --output comparison.html \
    --pattern "*_classifications.csv"
```

**Verbose Multi-Sample Processing**
```bash
metadarkmatter report multi \
    --input-dir ./results/ \
    --output comparison.html \
    --verbose
```

---

## Report Output Details

### HTML Report Structure

The generated HTML report is self-contained with:

- **Embedded CSS**: Professional styling with light/dark themes
- **Embedded Plotly.js**: Interactive charts (loaded from CDN)
- **Tab Navigation**: JavaScript-powered section switching
- **Data Table**: Sortable, filterable, searchable with pagination

### File Size Estimates

| Dataset | Reads | Report Size |
|---------|-------|-------------|
| Small | 10K | 0.5 MB |
| Medium | 100K | 1-2 MB |
| Large | 500K | 3-5 MB |
| Very Large | 1M+ | 5-10 MB |

### Browser Compatibility

Reports are compatible with modern browsers:
- Chrome 80+
- Firefox 75+
- Safari 13+
- Edge 80+
