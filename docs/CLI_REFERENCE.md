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
# Output: metadarkmatter version 0.2.0
```

`mdm` is an alias for `metadarkmatter`; both entry points are interchangeable.

### `metadarkmatter --help`

Show available commands and options.

### Global Options

These options go before the subcommand (e.g. `metadarkmatter --dry-run score classify ...`).

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--version` | `-v` | FLAG | off | Show version and exit |
| `--log-level` | | TEXT | WARNING | Logging level (DEBUG, INFO, WARNING, ERROR) |
| `--log-format` | | TEXT | text | Log format: text or json |
| `--log-file` | | TEXT | - | Log file path |
| `--dry-run` | | FLAG | off | Process-wide dry-run; subcommands that support dry-run honour it. Equivalent to `MDM_DRY_RUN=1` |
| `--debug` | | FLAG | off | Print full tracebacks on error instead of a short message. Equivalent to `MDM_DEBUG=1` |

### `metadarkmatter doctor`

Report Python version, dependency versions, external tool paths/versions, and the project environment variables. Run it before a long pipeline invocation.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--kraken-db` | | PATH | - | Optionally validate a Kraken2 database directory |
| `--strict` | | FLAG | off | Exit non-zero if a required tool or package is missing (blastn, makeblastdb, kraken2, skani or fastANI, core Python packages) or `--kraken-db` is invalid. Without it, doctor reports problems and exits 0. |

Optional tools are listed with the command that needs them: bowtie2/bowtie2-build and samtools (`map reads`, `visualize recruitment`), prodigal (`proteins predict`).

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

When `--genomes` is a directory, FASTA headers are rewritten to `{accession}|{contig_id}` and a contig mapping is written to `{output}_contig_mapping.tsv` (columns `contig_id`, `genome_accession`, `original_header`). If any standardized ID is longer than 50 characters, `makeblastdb` is run without `-parse_seqids` (subject IDs in the alignment output are unchanged).

---

### `metadarkmatter blast align`

Run BLASTN competitive alignment for ANI-weighted placement.

#### Required Options

| Option | Short | Type | Description |
|--------|-------|------|-------------|
| `--query` | | PATH | Single-end query reads (FASTA/FASTQ). Use this OR `--query-1`/`--query-2` |
| `--query-1` | `-1` | PATH | Forward (R1) reads for paired-end input (gzip OK) |
| `--query-2` | `-2` | PATH | Reverse reads for paired-end input (gzip OK) |
| `--database` | `-d` | PATH | BLAST database prefix |
| `--output` | `-o` | PATH | Output file for BLAST results |

#### Optional Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--word-size` | `-w` | INT | per task | Word size for seeds (smaller = more sensitive); default 7 for blastn, 11 for dc-megablast, 28 for megablast |
| `--evalue` | `-e` | FLOAT | 1e-5 | E-value threshold |
| `--max-targets` | `-m` | INT | 500 | Maximum target sequences per query (high for competitive alignment) |
| `--min-identity` | | FLOAT | None | Minimum percent identity filter |
| `--task` | | TEXT | blastn | BLAST task (blastn, megablast, dc-megablast) |
| `--threads` | `-p` | INT | 4 | Number of threads |
| `--compress/--no-compress` | | FLAG | auto | Compress output with gzip (default: only when `--output` ends in `.gz`; `--compress` appends `.gz`) |
| `--skip-if-exists` | | FLAG | False | Skip if output file exists |
| `--verbose` | `-v` | FLAG | False | Enable verbose output |
| `--quiet` | `-q` | FLAG | False | Suppress progress output |
| `--dry-run` | | FLAG | False | Show commands without executing (the global `--dry-run` is also honoured) |

Paired-end input (`--query-1`/`--query-2`) aligns both mates in one BLAST run and writes one combined output file. As with `mmseqs2 search`, each mate ID receives a `/1` or `/2` suffix unless it already ends in that suffix, so mates are classified as separate reads; R1 records precede R2 records, and every read's hits are contiguous, so the output can be classified with `--streaming`. Both mates must have the same format (FASTA or FASTQ).

Gzipped queries are decompressed to a temporary file before BLASTN runs; FASTQ input is converted to FASTA. The query format is detected from the first record, so gzipped FASTA is not mistaken for FASTQ. Multi-volume databases and database prefixes containing dots are accepted.

#### Examples

**Basic alignment**
```bash
metadarkmatter blast align \
    --query extracted_reads.fasta \
    --database blastdb/pangenome \
    --output sample.blast.tsv.gz
```

**Paired-end alignment (both mates, one output)**
```bash
metadarkmatter blast align \
    --query-1 extracted_R1.fastq.gz \
    --query-2 extracted_R2.fastq.gz \
    --database blastdb/pangenome \
    --output sample.blast.tsv.gz
```

**Faster alignment for closely related references**
```bash
metadarkmatter blast align \
    --query extracted_reads.fasta \
    --database blastdb/pangenome \
    --output sample.blast.tsv.gz \
    --task dc-megablast \
    --threads 16
```

With `--task dc-megablast` the word size defaults to 11; pass `--word-size` to override it.

**Output format**

The output is tabular format (outfmt 6) with 13 columns (`qlen` is used for coverage weighting):
```
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen
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
| `--skip-if-exists` | | FLAG | False | Skip if database already exists |
| `--verbose` | `-v` | FLAG | False | Enable verbose output |
| `--quiet` | `-q` | FLAG | False | Suppress progress output |
| `--dry-run` | | FLAG | False | Show commands without executing |

#### Examples

**From directory of protein files**
```bash
metadarkmatter blastx makedb \
    --proteins reference_proteins/ \
    --output blastdb/panproteome
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
| `--query` | | PATH | Query reads file (FASTA/FASTQ, optionally gzipped) |
| `--database` | `-d` | PATH | Diamond protein database prefix |
| `--output` | `-o` | PATH | Output file for BLASTX results |

#### Optional Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--evalue` | `-e` | FLOAT | 1e-5 | E-value threshold |
| `--min-identity` | | FLOAT | 30.0 | Minimum percent identity filter |
| `--max-targets` | `-k` | INT | 500 | Maximum target sequences per query |
| `--threads` | `-p` | INT | 4 | Number of threads |
| `--sensitive/--fast` | | FLAG | True | Use sensitive mode (recommended for divergent sequences) |
| `--compress/--no-compress` | | FLAG | auto | Compress output with gzip (default: only when `--output` ends in `.gz`; `--compress` appends `.gz`) |
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

**Note:** BLASTX translates DNA queries in all six reading frames before alignment against the protein database. Use `--alignment-mode protein` when running `score classify` on BLASTX output to use protein-calibrated thresholds. Genome accessions are derived from the protein file names with the same rules as nucleotide genomes (`extract_accession_from_filename`). Reverse-frame hits (qstart > qend) are normalised when the alignment is parsed.

---

## MMseqs2 Subcommand

The `mmseqs2` subcommand wraps MMseqs2 as a faster alternative to BLASTN. Output is BLAST-compatible 13-column tabular format (including `qlen`). See [mmseqs2-integration.md](mmseqs2-integration.md) for details.

### `metadarkmatter mmseqs2 makedb`

Create an MMseqs2 nucleotide database from reference genomes.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--genomes` | `-g` | PATH | - | Path to genome FASTA file or directory of genome files (required) |
| `--output` | `-o` | PATH | - | Output database path (e.g., 'mmseqs_db/pangenome') (required) |
| `--genome-pattern` | `-p` | TEXT | *.fna | Glob pattern for genome files when --genomes is a directory |
| `--dbtype` |  | INT | - | Database type (mmseqs createdb): 0=auto, 1=amino acid, 2=nucleotide (default: auto) |
| `--skip-if-exists` |  | FLAG | off | Skip database creation if output files already exist |
| `--verbose` | `-v` | FLAG | off | Enable verbose output |
| `--quiet` |  | FLAG | off | Suppress progress output |
| `--dry-run` |  | FLAG | off | Show commands without executing |

`--dbtype` is passed to `mmseqs createdb`: 0 = auto, 1 = amino acid, 2 = nucleotide.

### `metadarkmatter mmseqs2 search`

Run MMseqs2 fast sequence search with BLAST-compatible output.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--query` | `-q` | PATH | - | Query reads file (FASTA/FASTQ). Use this OR --query-1/--query-2. |
| `--query-1` | `-1` | PATH | - | Forward reads for paired-end data (FASTA/FASTQ, gzipped OK) |
| `--query-2` | `-2` | PATH | - | Reverse reads for paired-end data (FASTA/FASTQ, gzipped OK) |
| `--database` | `-d` | PATH | - | MMseqs2 database path (from 'mmseqs2 makedb') (required) |
| `--output` | `-o` | PATH | - | Output file for search results (TSV, optionally .gz) (required) |
| `--query-db` |  | PATH | - | Pre-computed query database for caching (optional, speeds up repeated searches) |
| `--sensitivity` | `-s` | FLOAT | 5.7 | Sensitivity: 1.0 (fast) to 7.5 (sensitive), default 5.7 (balanced) |
| `--evalue` | `-e` | FLOAT | 0.001 | E-value threshold (default 1e-3) |
| `--max-seqs` | `-m` | INT | 500 | Max target sequences per query (high for competitive alignment) |
| `--min-identity` |  | FLOAT | - | Minimum percent identity filter (0-100 scale, optional) |
| `--search-type` |  | INT | 3 | Search type: 2=translated, 3=nucleotide (default: 3) |
| `--threads` | `-p` | INT | 4 | Number of threads |
| `--tmp-dir` |  | PATH | - | Temporary directory for MMseqs2 (auto-created if not specified) |
| `--compress/--no-compress` |  | FLAG | - | Compress output with gzip (default: auto, compress only when --output ends in .gz; --compress appends .gz) |
| `--skip-if-exists` |  | FLAG | off | Skip search if output file already exists |
| `--verbose` | `-v` | FLAG | off | Enable verbose output |
| `--quiet` |  | FLAG | off | Suppress progress output |
| `--dry-run` |  | FLAG | off | Show commands without executing |

Notes:

- Use `--query` for single-end input, or `--query-1` with `--query-2` for paired-end input. Paired files are concatenated and each mate ID receives a `/1` or `/2` suffix (unless it already ends in that suffix), so mates are classified as separate reads. Gzip compression is detected per file.
- `--search-type 2` runs a translated search; the query database is then built as nucleotide.
- Output compression follows the `--output` extension (`.gz` => gzip); `--compress` forces it and appends `.gz`.

```bash
metadarkmatter mmseqs2 makedb --genomes genomes/ --output mmseqs_db/pangenome
metadarkmatter mmseqs2 search \
    --query-1 sample_taxid34064_R1.fastq.gz --query-2 sample_taxid34064_R2.fastq.gz \
    --database mmseqs_db/pangenome --output sample.mmseqs2.tsv.gz --threads 16
```

---

## ANI and AAI Subcommands

### `metadarkmatter ani compute`

Compute ANI matrix from reference genomes.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--genomes` | `-g` | PATH | - | Directory containing genome FASTA files (required) |
| `--output` | `-o` | PATH | - | Output CSV file for ANI matrix (required) |
| `--metadata` | `-m` | PATH | - | Path to genome_metadata.tsv (required for --representatives-only) |
| `--representatives-only` |  | FLAG | off | Compute ANI only for species representative genomes (from metadata). Reduces O(N^2) comparisons for large families. |
| `--backend` | `-b` | CHOICE | auto | ANI computation backend: fastani, skani, or auto Choices: fastani, skani, auto. |
| `--fastani-batches` |  | INT | 1 | Run fastANI as N parallel subprocess batches over query chunks. Default 1 uses a single fastANI process with internal threading. Higher values (e.g. 4 or 8) often speed up large all-vs-all runs because each batch owns its own memory. Ignored by the skani backend. |
| `--sensitivity` | `-s` | CHOICE | sensitive | Sensitivity preset: default (standard), sensitive (captures distant genomes) Choices: default, sensitive, custom. |
| `--min-af` |  | FLOAT | - | Minimum aligned fraction (0-100). Lower = more distant comparisons. Default: 15 (default mode), 5 (sensitive mode) |
| `--screen` |  | FLOAT | - | Pre-filter threshold for ANI (0-100). Lower = attempt more comparisons. Default: 80 (default mode), 50 (sensitive mode). Only applies to skani. |
| `--skani-mode` |  | CHOICE | default | skani sensitivity mode: default, medium (`-c 70`) or slow (`-c 30`). More sensitive modes detect more distant pairs but are slower. Only applies to skani. |
| `--skani-c` |  | INT | - | skani compression factor (`-c`); overrides `--skani-mode`. Only applies to skani. |
| `--threads` | `-t` | INT | 4 | Number of threads for parallel execution |
| `--genome-pattern` | `-p` | TEXT | *.fna | Glob pattern for genome files |
| `--verbose` | `-v` | FLAG | off | Enable verbose output |
| `--quiet` | `-q` | FLAG | off | Suppress progress output |
| `--dry-run` |  | FLAG | off | Show commands without executing |

With `--backend auto`, skani is used when it is on `PATH`, otherwise fastANI. `--representatives-only` requires `--metadata` and is honoured by both backends. Matrix values are written on a 0-100 scale. fastANI asymmetric pairs are averaged.

skani reports nothing for pairs it does not detect (written as 0 = not computed). When more than 20% of genome pairs are undetected, `ani compute` prints a warning with the count. On a 10-genome Francisellaceae set, skani left 22 of 45 pairs undetected (also with `--skani-mode slow`) while fastANI computed all of them (lowest 78%); undetected pairs are treated as 70% ANI during classification, so prefer `--backend fastani` when many pairs are missing.

### `metadarkmatter ani validate`

Validate ANI matrix coverage against BLAST results.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--ani` | `-a` | PATH | - | ANI matrix CSV file to validate (required) |
| `--blast` | `-b` | PATH | - | BLAST results file to check coverage against (required) |
| `--metadata` | `-m` | PATH | - | Path to genome_metadata.tsv with representative column. When provided, validates that representative genomes are in ANI matrix rather than requiring all BLAST genomes. |
| `--sample-rows` |  | INT | 100000 | Number of BLAST rows to sample for coverage check |
| `--verbose` | `-v` | FLAG | off | List missing genomes |
| `--quiet` | `-q` | FLAG | off | Suppress progress output |

### `metadarkmatter aai compute`

Compute AAI matrix from genome protein files.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--genomes` | `-g` | PATH | - | Directory containing genome protein FASTA files (.faa) (required) |
| `--output` | `-o` | PATH | - | Output CSV file for AAI matrix (required) |
| `--threads` | `-t` | INT | 4 | Number of threads for Diamond |
| `--min-identity` | `-i` | FLOAT | 30.0 | Minimum percent identity for Diamond hits |
| `--min-coverage` | `-c` | FLOAT | 0.5 | Minimum query coverage fraction for reciprocal best hits |
| `--evalue` | `-e` | FLOAT | 1e-05 | E-value threshold for Diamond |
| `--protein-pattern` | `-p` | TEXT | *.faa | Glob pattern for protein FASTA files |
| `--keep-intermediates` | `-k` | FLAG | off | Keep intermediate files (Diamond output, RBH pairs) |
| `--sensitive` | `-s` | FLAG | off | Use Diamond sensitive mode (slower but more accurate for divergent sequences) |
| `--verbose` | `-v` | FLAG | off | Enable verbose output |
| `--quiet` | `-q` | FLAG | off | Suppress progress output |
| `--dry-run` |  | FLAG | off | Show what would be computed without executing |

### `metadarkmatter aai validate`

Validate AAI matrix coverage against ANI matrix.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--aai` | `-a` | PATH | - | AAI matrix CSV file to validate (required) |
| `--ani` |  | PATH | - | ANI matrix CSV file to compare coverage against (required) |
| `--verbose` | `-v` | FLAG | off | List genomes with missing AAI/ANI pairs |
| `--quiet` | `-q` | FLAG | off | Suppress progress output |

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

**Output files** (`--output` is a directory):
- `{output}/{sample}.kraken` - Per-read classification output
- `{output}/{sample}.kreport` - Hierarchical taxonomic report

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
| `--compress/--no-compress` | | FLAG | True | Compress output files with gzip (default: compress) |
| `--verbose` | `-v` | FLAG | False | Enable verbose output |
| `--quiet` | `-q` | FLAG | False | Suppress progress output |
| `--dry-run` | | FLAG | False | Show commands without executing |

The Kraken2 report must sit next to the `.kraken` file: either `{name}.kreport` with the same stem, or `{sample}.kreport` in the same directory (as written by `kraken2 classify`). The sample name defaults to the `.kraken` file name with common suffixes removed.

#### Examples

**Extract Family Reads (include child taxa; 34064 = Francisellaceae)**
```bash
metadarkmatter kraken2 extract \
    --kraken-output kraken_output/sample.kraken \
    --reads-1 sample_R1.fastq.gz \
    --reads-2 sample_R2.fastq.gz \
    --taxid 34064 \
    --output extracted/
```

**Exact Match Only**
```bash
metadarkmatter kraken2 extract \
    --kraken-output kraken_output/sample.kraken \
    --reads-1 sample_R1.fastq.gz \
    --taxid 34064 \
    --output extracted/ \
    --exact-match
```

**Output files** (`.fastq` instead of `.fastq.gz` with `--no-compress`):
- `{sample}_taxid{taxid}_R1.fastq.gz` - Extracted forward reads
- `{sample}_taxid{taxid}_R2.fastq.gz` - Extracted reverse reads (if paired)

---

## Score Subcommand

The `score` subcommand performs ANI-weighted placement classification on BLAST results.

### `metadarkmatter score classify`

Classify metagenomic reads from an alignment file (BLAST or MMseqs2).

#### Required Options

| Option | Short | Type | Description |
|--------|-------|------|-------------|
| `--alignment` | | PATH | Path to alignment file (.tsv or .tsv.gz) - accepts BLAST or MMseqs2 tabular format |
| `--ani` | `-a` | PATH | Path to ANI matrix file (CSV or TSV) |
| `--output` | `-o` | PATH | Output path for classification results |

#### Optional Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--aai` | | PATH | None | AAI matrix (CSV) for genus-level classification |
| `--summary` | `-s` | PATH | None | Output path for summary statistics (JSON) |
| `--metadata` | `-m` | PATH | None | Genome metadata TSV: species/genus/family columns, representative mapping, and target-family inference |
| `--genomes` | `-g` | PATH | None | Genome FASTA directory for auto-generating ID mapping (external results) |
| `--id-mapping` | | PATH | None | Pre-computed ID mapping TSV (alternative to --genomes) |
| `--config` | | PATH | None | YAML configuration (generate a template with `score export-config`); explicitly passed flags override YAML values |
| `--preset` | | TEXT | None | Threshold preset: `default`, `gtdb-strict`, `gtdb-relaxed`, `conservative`, `literature-strict`, `coverage-strict`, `adaptive`. Explicitly passed flags override preset values |
| `--alignment-mode` | | TEXT | nucleotide | Alignment type: 'nucleotide' for BLASTN/MMseqs2, 'protein' for BLASTX |
| `--bitscore-threshold` | | FLOAT | 95.0 | Percentage of top bitscore for ambiguous hits (0-100) |
| `--min-alignment-length` | | INT | 100 | Minimum alignment length in bp (0 = no filter). In `--alignment-mode protein` it is compared with the query span in nucleotides (`qend - qstart + 1`), because BLASTX/DIAMOND `length` is in amino acids. |
| `--min-alignment-fraction` | | FLOAT | 0.3 | Minimum fraction of the read aligned (uses `qlen` when present); 0.5 is GTDB-compatible |
| `--max-evalue` | | FLOAT | 0.0 | Maximum E-value for hit filtering (0 disables) |
| `--min-percent-identity` | | FLOAT | 0.0 | Minimum percent identity for hit filtering (0 disables) |
| `--min-bitscore` | | FLOAT | 0.0 | Minimum bitscore for hit filtering (0 disables) |
| `--min-read-length` | | INT | 0 | Minimum read length in bp (0 disables) |
| `--min-query-coverage` | | FLOAT | 0.0 | Minimum query coverage percentage, 0-100 (0 disables) |
| `--single-hit-threshold` | | FLOAT | 10.0 | Inferred-uncertainty threshold above which single-hit reads become Ambiguous in `legacy_call` |
| `--format` | `-f` | TEXT | csv | Output format: 'csv' or 'parquet' |
| `--strict-ani` | | FLAG | False | Refuse ANI matrices whose forward/reverse entries differ by more than 0.5 (default: warn) |
| `--include-legacy-scores` | | FLAG | False | Add legacy sub-scores (alignment_quality, identity_confidence, placement_confidence, discovery_score) |
| `--clade-novel-species/--no-clade-novel-species` | | FLAG | on | Stage 2 rule (nucleotide mode, on by default; `--no-clade-novel-species` disables it): call Ambiguous/Species Boundary reads Novel Species when all competitors are within one genus and the read is at least 5% divergent from both the best and the second-best reference (and below 18.6%). Adds a `clade_novel_species` column. See METHODS.md section 10.5 and CLADE_NOVELTY_VALIDATION.md. |
| `--single-hit-novel-species/--no-single-hit-novel-species` | | FLAG | on | Stage 2 rule (nucleotide mode, on by default): call Known Species reads with no competing reference Novel Species when their novelty is >= the Known Species threshold (`novelty_known_max`, 4%). Adds a `single_hit_novel_species` column. See METHODS.md section 10.6. |
| `--verbose` | `-v` | FLAG | False | Enable verbose output |
| `--provenance/--no-provenance` | | FLAG | on | Write `<output>.provenance.json` (see [Provenance sidecar](#provenance-sidecar)) and embed it under `"provenance"` in `--summary` |
| `--quiet` | `-q` | FLAG | False | Suppress progress output (for scripting) |
| `--target-family` | | TEXT | None | Override the target family for off-target detection. Auto-inferred from `--metadata` when omitted (e.g., `f__Francisellaceae`). |
| `--family-ratio-threshold` | | FLOAT | 0.8 | Bitscore ratio threshold for off-target detection (0.0-1.0) |
| `--dry-run` | | FLAG | False | Validate inputs without processing (honours `MDM_DRY_RUN=1` env var and the global top-level `--dry-run`) |

#### Family Validation Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--target-family` | TEXT | None | Override the target family for off-target read detection. Optional - leave unset and family validation will auto-infer from `--metadata`. |
| `--family-ratio-threshold` | FLOAT | 0.8 | Bitscore ratio below which a read is classified as Off-target |

Family validation activates whenever metadarkmatter knows the target family: either explicitly via `--target-family` or auto-inferred from the most common entry in the `family` column of `--metadata`. The classifier partitions BLAST hits into in-family (genome present in the ANI matrix) and external hits, then flags reads whose best in-family bitscore is substantially lower than their best overall bitscore. Off-target reads carry `placement_uncertainty = null` and `confidence_score = 10.0` with `low_confidence = True`; filter on `taxonomic_call == "Off-target"` downstream.

#### Advanced Classification Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--bayesian` | FLAG | False | Deprecated, no effect beyond a warning: Bayesian posteriors are always computed |
| `--adaptive-thresholds` | FLAG | False | Detect species boundary from ANI matrix via GMM (requires scikit-learn) |
| `--qc-output` | PATH | None | Output path for QC metrics JSON (not available with `--streaming`) |
| `--uncertainty-mode` | TEXT | second | Uncertainty calculation: 'second' (second-best genome) or 'max' (maximum ANI to any competing genome) |

**Bayesian posteriors (always on):** `taxonomic_call` is the maximum a posteriori category from 2D Gaussian likelihoods over novelty and uncertainty, followed by a discrete refinement step (Conserved Region, Ambiguous Within Genus). The output carries six posterior columns (`p_known_species`, `p_novel_species`, `p_novel_genus`, `p_species_boundary`, `p_ambiguous`, `p_unclassified`), `posterior_entropy` (0 = confident, log2(6) = 2.585 = uniform), and `confidence_score = (1 - entropy / log2(6)) * 100`. The hard-threshold cascade is kept in `legacy_call`. See [METHODS.md](METHODS.md) for detailed formulas.

**Adaptive Thresholds (`--adaptive-thresholds`):** Fits a 2-component Gaussian Mixture Model to the ANI matrix distribution to detect the natural species boundary. Overrides the default `novelty_known_max` threshold with the detected value. Falls back to the default 96% threshold if the GMM does not converge, the components are not well separated, or the detected boundary is outside the plausible species range of 92-99% ANI (for example when the matrix holds no within-species pairs). Requires scikit-learn: `pip install metadarkmatter[adaptive]`. It has no effect in protein mode, where the fixed protein thresholds are always used.

**QC Metrics (`--qc-output`):** Writes a JSON file with pre-classification metrics (filter rate, genome coverage, single-hit fraction, mean identity) and post-classification metrics (ambiguous fraction, low-confidence fraction). Includes automated warnings for problematic inputs (high filter rate, low genome coverage).

**Uncertainty Mode (`--uncertainty-mode`):** Controls how placement uncertainty is calculated. `second` (default) uses ANI to the second-best genome only. `max` uses the maximum ANI across all competing genomes. The `second` mode is generally preferred as it better reflects the specificity of placement.

#### Coverage Weighting Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--coverage-weight-mode` | TEXT | linear | Coverage weighting mode: 'none', 'linear', 'log', 'sigmoid' |
| `--coverage-weight-strength` | FLOAT | 0.5 | Coverage weight strength (0.0-1.0) |

**Coverage Weighting Modes:**
- `none`: Use raw bitscore for hit selection
- `linear`: Gradual penalty for low coverage
- `log`: Rewards any coverage, diminishing returns
- `sigmoid`: Sharp threshold around 60% coverage

#### Alignment Mode

The `--alignment-mode` option selects threshold calibration:

| Mode | Source | Thresholds (`legacy_call`) |
|------|--------|------------|
| `nucleotide` | BLASTN / MMseqs2 output | Known: N<4%, Novel Species: 4-20%, Novel Genus: 20-25% (U < 1.5%) |
| `protein` | BLASTX output | Known: N<10%, Novel Species: 10-25%, Novel Genus: 25-40% (U < 5%) |

Protein thresholds are wider because amino acid sequences diverge more slowly than nucleotide sequences. In protein mode these fixed values are used regardless of `--preset`, `--config`, or `--adaptive-thresholds` novelty settings.

These hard thresholds produce the `legacy_call` column. The primary `taxonomic_call` comes from the Bayesian posteriors, whose decision boundaries do not coincide exactly with the table. With the default parameters and U = 0, the measured nucleotide-mode transitions are approximately N = 6.2 (Known -> Novel Species), 18.6 (Novel Species -> Novel Genus) and 27.6 (Novel Genus -> Unclassified); at N = 1, Species Boundary becomes the MAP call from about U = 3.8. In protein mode the corresponding novelty transitions are about 11.3, 25.1 and 40.4. Single-hit reads use inferred uncertainty only in `legacy_call`.

#### Processing Mode Options

| Option | Description | Use Case |
|--------|-------------|----------|
| `--streaming` | Memory-efficient, processes the file in batches | Very large files (>5M alignments) or constrained memory |
| `--chunk-size` | Rows per streaming batch (default 1,000,000; range 100 to 100,000,000) | Lower it (e.g. 500000) to cap memory |

By default, all classification uses the VectorizedClassifier (Polars-based, auto-parallelized). `--streaming` processes the alignment file in batches to keep memory bounded. Restrictions:

- `--genomes` and `--id-mapping` are rejected with an error in `--streaming` mode; pre-transform the subject IDs or use the in-memory path.
- With `--metadata`, the representative mapping and target-family inference still apply, but the `species`/`genus`/`family` output columns are omitted.
- All hits for a read must be contiguous in the file (as BLAST and MMseqs2 write them); otherwise the run stops with an error. Sort with e.g. `sort -k1,1 -s` if needed.

#### Examples

**Basic Classification**
```bash
metadarkmatter score classify \
    --alignment sample.blast.tsv.gz \
    --ani genomes.ani.csv \
    --output sample_classifications.csv
```

**With Summary and Parquet Output**
```bash
metadarkmatter score classify \
    --alignment sample.blast.tsv.gz \
    --ani genomes.ani.csv \
    --output sample_classifications.parquet \
    --summary sample_summary.json \
    --format parquet
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

**External BLAST Results with ID Mapping**
```bash
# Using pre-generated mapping file
metadarkmatter score classify \
    --alignment external_results.tsv.gz \
    --ani ani_matrix.csv \
    --id-mapping id_mapping.tsv \
    --output classifications.csv

# Auto-generate mapping from genome directory
metadarkmatter score classify \
    --alignment external_results.tsv.gz \
    --ani ani_matrix.csv \
    --genomes reference_genomes/ \
    --output classifications.csv
```

See the [Workflow Guide](WORKFLOW.md#importing-external-alignment-results) for the complete external import workflow.

**Protein Mode (for BLASTX output)**
```bash
metadarkmatter score classify \
    --alignment sample.blastx.tsv.gz \
    --ani genomes.ani.csv \
    --output sample_classifications.csv \
    --alignment-mode protein
```

**Family Validation (broad-database alignment)**
```bash
# Recommended: auto-infer the target family from --metadata.
metadarkmatter score classify \
    --alignment broad_results.tsv.gz \
    --ani family_ani_matrix.csv \
    --metadata genome_metadata.tsv \
    --family-ratio-threshold 0.8 \
    --output classifications.csv

# Override: pin the target family explicitly.
metadarkmatter score classify \
    --alignment broad_results.tsv.gz \
    --ani family_ani_matrix.csv \
    --target-family "f__Francisellaceae" \
    --family-ratio-threshold 0.8 \
    --output classifications.csv
```

**YAML Configuration**
```bash
metadarkmatter score export-config --output my_config.yaml
# edit my_config.yaml, then:
metadarkmatter score classify \
    --alignment sample.blast.tsv.gz \
    --ani genomes.ani.csv \
    --output sample_classifications.csv \
    --config my_config.yaml
```

**Adaptive Thresholds with QC Output**
```bash
metadarkmatter score classify \
    --alignment sample.blast.tsv.gz \
    --ani genomes.ani.csv \
    --output sample_classifications.csv \
    --adaptive-thresholds \
    --qc-output sample_qc.json
```

**Full Advanced Classification**
```bash
# Target family auto-inferred from --metadata; --target-family
# can be added to override.
metadarkmatter score classify \
    --alignment sample.blast.tsv.gz \
    --ani genomes.ani.csv \
    --metadata genome_metadata.tsv \
    --output sample_classifications.csv \
    --adaptive-thresholds \
    --qc-output sample_qc.json
```

#### Provenance sidecar

`score classify` writes `<output>.provenance.json` next to the classification
table (for example `classifications.csv.provenance.json`; disable with
`--no-provenance`). The `--summary` JSON also carries the same object under the
key `"provenance"`. The classification table itself is unchanged (no extra
columns). Schema `metadarkmatter.provenance/1`:

```json
{
  "schema": "metadarkmatter.provenance/1",
  "metadarkmatter_version": "0.2.0",
  "created_utc": "2026-11-05T17:59:10Z",
  "command": ["mdm", "score", "classify", "..."],
  "python": "3.13.1",
  "platform": "Linux-...",
  "seed": 42,
  "inputs": {
    "alignment":  {"path": "/abs/sample.blast.tsv.gz", "size": 651400, "sha256": "..."},
    "ani":        {"path": "...", "size": 1099, "sha256": "..."},
    "metadata":   {"path": "...", "size": 2284, "sha256": "..."},
    "aai": null,
    "id_mapping": null
  },
  "config": {"...": "ScoringConfig as applied (adaptive thresholds and inferred target family included)"},
  "outputs": {"classifications": "/abs/classifications.csv", "summary": "/abs/summary.json"},
  "results": {"total_reads": 3973, "taxonomic_call_counts": {"Ambiguous": 639, "Known Species": 2205}}
}
```

Input files are hashed by streaming, so large alignments are not loaded into
memory. Everything except `created_utc` is identical between two runs with the
same inputs, options and `METADARKMATTER_SEED`. A `--genomes` directory is not
hashed (use `--id-mapping` if the mapping must be pinned).

---

### `metadarkmatter score batch`

Batch classify multiple alignment files in a directory. Each file goes through
the same code path as `score classify` (`run_classification` and the same
config builder), so for the same options a sample's table has exactly the same
columns and values as a single `score classify` run.

#### Required Options

| Option | Short | Type | Description |
|--------|-------|------|-------------|
| `--alignment-dir` | `-b` | PATH | Directory containing alignment files (BLAST or MMseqs2) |
| `--ani` | `-a` | PATH | Path to ANI matrix file (CSV or TSV) |
| `--output-dir` | `-o` | PATH | Output directory for classification results |

#### Optional Options

Accepts the classify options below with the same meaning and defaults:
`--aai`, `--metadata`, `--genomes`, `--id-mapping`, `--config`, `--preset`,
`--alignment-mode`, `--target-family`, `--family-ratio-threshold`,
`--bitscore-threshold`, `--min-alignment-length`, `--min-alignment-fraction`,
`--max-evalue`, `--min-percent-identity`, `--min-bitscore`,
`--min-read-length`, `--min-query-coverage`, `--coverage-weight-mode`,
`--coverage-weight-strength`, `--uncertainty-mode`,
`--single-hit-threshold`, `--include-legacy-scores`,
`--clade-novel-species/--no-clade-novel-species`,
`--single-hit-novel-species/--no-single-hit-novel-species`, `--streaming`,
`--chunk-size`, `--strict-ani`, `--adaptive-thresholds`, `--format`,
`--verbose`, `--provenance/--no-provenance`. Config precedence is the same as
for classify (`--config` > flags > `--preset` > defaults). `--qc-output` and
`--dry-run` are not available in batch.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--pattern` | `-p` | TEXT | `*.blast.tsv.gz`, `*.mmseqs2.tsv.gz`, `*.blast.tsv`, `*.mmseqs2.tsv` | Glob pattern(s) for alignment files; repeat the option to give several |
| `--workers` | `-w` | INT | 0 | Deprecated and ignored (hidden); a warning is printed when set |

For each input file `{sample}_classifications.csv` (or `.parquet`),
`{sample}_summary.json` and (unless `--no-provenance`)
`{sample}_classifications.csv.provenance.json` are written to `--output-dir`.

**Sample names.** The name is the file name with `.blast.tsv.gz`,
`.mmseqs2.tsv.gz`, `.blast.tsv`, `.mmseqs2.tsv`, `.tsv.gz` or `.tsv` removed.
If two files map to the same name (for example `s1.blast.tsv.gz` and
`s1.mmseqs2.tsv.gz`) the command stops before processing anything and lists
the files; narrow `--pattern` or rename a file.

**Sequential processing.** Samples run one after another, in sorted file
order. `--workers` was never implemented and is now deprecated rather than
implemented: each sample already uses all cores through Polars, and N parallel
processes would multiply peak memory (about 18 GB at 10M hits in vectorised
mode) for little gain. Use `--streaming` to bound memory, or run several
`score batch` invocations on disjoint directories.

#### Examples

**Basic Batch Processing (BLAST and MMseqs2 outputs)**
```bash
metadarkmatter score batch \
    --alignment-dir ./alignments/ \
    --ani genomes.ani.csv \
    --metadata genome_metadata.tsv \
    --output-dir ./classifications/
```

**Only MMseqs2 outputs, Parquet**
```bash
metadarkmatter score batch \
    --alignment-dir ./alignments/ \
    --ani genomes.ani.csv \
    --output-dir ./classifications/ \
    --pattern "*.mmseqs2.tsv.gz" \
    --format parquet
```

**Streaming mode for large samples**
```bash
metadarkmatter score batch \
    --alignment-dir ./alignments/ \
    --ani genomes.ani.csv \
    --metadata genome_metadata.tsv \
    --output-dir ./classifications/ \
    --streaming --chunk-size 500000
```

---

### `metadarkmatter score export-config`

Export an editable YAML configuration file.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--output` | `-o` | PATH | - | Output path for the YAML config file (required) |
| `--preset` | `-p` | TEXT | - | Base preset to export: 'default', 'gtdb-strict', 'gtdb-relaxed', 'conservative', 'literature-strict', 'coverage-strict', 'adaptive'. Defaults to built-in defaults. |

### `metadarkmatter score extract-novel`

Extract candidate novel species or genera from classification results.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--classifications` | `-c` | PATH | - | Path to classification results file (CSV or Parquet) (required) |
| `--output` | `-o` | PATH | - | Output path for novel species candidates (CSV) (required) |
| `--read-ids` | `-r` | PATH | - | Output path for read IDs only (for seqtk extraction) |
| `--category` | `-t` | TEXT | all | Category to extract: 'species', 'genus', or 'all' (both) |
| `--min-novelty` |  | FLOAT | - | Minimum novelty index threshold (overrides category defaults) |
| `--max-uncertainty` |  | FLOAT | 2.0 | Maximum placement uncertainty (exclude ambiguous placements). Reads flagged `clade_novel_species` are exempt: their U reflects distant within-genus competitors, not uncertain novelty. |
| `--group-by-genome/--no-group` |  | FLAG | on | Group results by closest reference genome (default: group) |
| `--min-reads` |  | INT | 1 | Minimum reads per reference genome; applies to the grouped and per-read (`--no-group`) output, `--read-ids` and the read export |
| `--reads-1` |  | PATH | - | Original FASTQ/FASTA(.gz) (R1 or single-end) to extract the selected reads from; needs an output below |
| `--reads-2` |  | PATH | - | R2 FASTQ/FASTA(.gz); outputs become `{stem}_R1` / `{stem}_R2` |
| `--fastq-out` |  | PATH | - | Write selected reads as FASTQ (gzip if the path ends in `.gz`) |
| `--fasta-out` |  | PATH | - | Write selected reads as FASTA (gzip if `.gz`; FASTQ input is converted) |
| `--per-lineage-dir` |  | PATH | - | Also write one gzipped file per selected lineage (closest reference genome) |
| `--verbose` | `-v` | FLAG | off | Enable verbose output |
| `--quiet` | `-q` | FLAG | off | Suppress progress output |

#### Exporting the selected reads

With `--reads-1` the command streams the original reads once (gzip handled
transparently, memory bounded by the set of IDs, not by the FASTQ) and writes
the records whose ID was selected. The selection is the same as for
`--read-ids` (including the `--min-reads` lineage filter).

- **ID matching.** Only the first whitespace-delimited token of each header
  is used and a trailing `/1` or `/2` is removed on both sides. IDs classified
  from a paired MMseqs2 search (which appends `/1` and `/2`) therefore match
  FASTQ headers with or without the suffix.
- **Pairs.** A pair is written when either mate's ID is selected, so
  `{stem}_R1` and `{stem}_R2` stay in sync. With `--fastq-out novel.fastq.gz`
  and `--reads-2` the files are `novel_R1.fastq.gz` and `novel_R2.fastq.gz`.
- **Per lineage.** `--per-lineage-dir` writes `{genome}.fastq.gz` (or
  `{genome}_R1/_R2.fastq.gz`) per selected `best_match_genome`; a read whose
  mates hit different genomes goes to the alphabetically first one.
- **Report.** The command prints requested versus found IDs per input file and
  warns about missing IDs (usually a sign that the reads are not the files that
  were aligned). `--fastq-out` requires FASTQ input.

```bash
metadarkmatter score extract-novel \
    --classifications classifications.csv --output novel.csv \
    --reads-1 sample_R1.fastq.gz --reads-2 sample_R2.fastq.gz \
    --fastq-out novel/novel.fastq.gz --per-lineage-dir novel/by_lineage/
```

### `metadarkmatter score sensitivity`

Sweep classification thresholds and report category counts.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--classifications` | `-c` | PATH | - | Classification CSV/Parquet produced by 'score classify'. Must contain the columns novelty_index, placement_uncertainty, num_ambiguous_hits, and identity_gap. (required) |
| `--output` | `-o` | PATH | - | Output file path (TSV or JSON based on --format). (required) |
| `--config` |  | PATH | - | Base ScoringConfig YAML to vary. Defaults to built-in defaults. Threshold ranges sweep around the novelty_known_max and uncertainty_known_max fields. |
| `--novelty-min` |  | FLOAT | 2.0 | Lower bound of the novelty_known_max sweep (default 2.0). |
| `--novelty-max` |  | FLOAT | 8.0 | Upper bound of the novelty_known_max sweep (default 8.0). |
| `--uncertainty-min` |  | FLOAT | 0.5 | Lower bound of the uncertainty_known_max sweep (default 0.5). |
| `--uncertainty-max` |  | FLOAT | 3.0 | Upper bound of the uncertainty_known_max sweep (default 3.0). |
| `--steps` |  | INT | 9 | Number of threshold points to evaluate (default 9). |
| `--format` | `-f` | TEXT | tsv | Output format: 'tsv' (long form) or 'json'. |
| `--quiet` | `-q` | FLAG | off | Suppress progress output. |

### `metadarkmatter score evaluate`

Score classifications against a labelled ground-truth column.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--predictions` | `-p` | PATH | - | Classification CSV/Parquet produced by 'score classify'. Must contain a 'taxonomic_call' column; 'confidence_score' is required for calibration metrics. (required) |
| `--truth-column` |  | TEXT | true_category | Name of the ground-truth column in the predictions file. |
| `--output` | `-o` | PATH | - | Optional JSON path for the full evaluation result; without it only the console summary is printed. |
| `--bins` |  | INT | 10 | Number of equal-width confidence bins for ECE (default 10). |
| `--quiet` | `-q` | FLAG | off | Suppress progress and summary output. |

`score sensitivity` re-applies the hard-threshold cascade (the rules behind `legacy_call`) to a saved classification file at each threshold point and writes one row per threshold point and category.

---

## Util Subcommand

The `util` subcommand provides helper commands for data preparation, particularly for importing external alignment results.

### `metadarkmatter util generate-mapping`

Generate an ID mapping file from a genome FASTA directory. This maps original contig IDs to their parent genome accessions, enabling classification of externally-run BLAST or MMseqs2 results.

#### Required Options

| Option | Short | Type | Description |
|--------|-------|------|-------------|
| `--genomes` | `-g` | PATH | Directory containing genome FASTA files |
| `--output` | `-o` | PATH | Output TSV file for ID mapping |

#### Optional Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--pattern` | `-p` | TEXT | *.fna | Glob pattern for genome files |
| `--allow-duplicates` | | FLAG | False | Tolerate contig IDs that occur in more than one genome file (by default duplicates are an error) |
| `--quiet` | `-q` | FLAG | False | Suppress progress output |

#### Examples

**Basic usage**
```bash
metadarkmatter util generate-mapping \
    --genomes reference_genomes/ \
    --output id_mapping.tsv
```

**Non-standard file extensions**
```bash
metadarkmatter util generate-mapping \
    --genomes reference_genomes/ \
    --output id_mapping.tsv \
    --pattern "*.fasta.gz"
```

**Output format**

Tab-separated file with two columns:

```
original_contig_id	genome_accession
NZ_CP007557.1	GCF_000195955.2
NC_006570.2	GCF_000008985.1
```

---

### `metadarkmatter util validate-mapping`

Validate an ID mapping file and optionally check coverage against a BLAST results file.

#### Required Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `MAPPING_FILE` | PATH | Path to ID mapping TSV file |

#### Optional Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--blast` | `-b` | PATH | None | BLAST TSV file to check mapping coverage |
| `--min-coverage` | | FLOAT | None | With `--blast`: exit non-zero when the percentage of subject IDs covered by the mapping (or already in `accession|contig` form) is below this value |
| `--quiet` | `-q` | FLAG | False | Suppress progress output |

A contig ID mapped to two different accessions is an error (exit 1, conflicts listed). Subject IDs already in the standardized `accession|contig` form count as covered; NCBI-style `ref|ID|` IDs are looked up by their bare ID.

#### Examples

**Validate mapping file format**
```bash
metadarkmatter util validate-mapping id_mapping.tsv
```

**Check coverage against BLAST results**
```bash
metadarkmatter util validate-mapping id_mapping.tsv \
    --blast external_results.blast.tsv.gz
```

Reports the number of unique subject IDs in the BLAST file and what fraction can be mapped. Warns if coverage falls below 95%.

---

### `metadarkmatter util generate-metadata`

Generate genome_metadata.tsv from a directory of genome FASTA files.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--genomes` | `-g` | PATH | - | Directory containing genome FASTA files (required) |
| `--output` | `-o` | PATH | - | Output path for genome_metadata.tsv (required) |
| `--pattern` | `-p` | TEXT | *.fna | Glob pattern for genome files |
| `--gtdb-metadata` |  | PATH | - | GTDB metadata file for representative assignment |
| `--family` |  | TEXT | - | Override family name (instead of parsing from filename) |
| `--quiet` | `-q` | FLAG | off | Suppress progress output |

### `metadarkmatter util select-representatives`

Assign or update species representative genomes in metadata.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--metadata` | `-m` | PATH | - | Path to genome_metadata.tsv file (required) |
| `--output` | `-o` | PATH | - | Output path for updated metadata TSV (required) |
| `--quiet` | `-q` | FLAG | off | Suppress progress output |

---

## Other Subcommands

### `metadarkmatter download genomes list`

Query GTDB API for representative genome accessions.

Positional argument: `TAXON`

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--output` | `-o` | PATH | - | Output TSV file for accession list (required) |
| `--reps-only/--all-genomes` |  | FLAG | on | Only include representative genomes (default) or all genomes. With `--all-genomes`, each species gets exactly one representative in `genome_metadata.tsv` (the GTDB representative, or the first accession in sorted order when none was returned, with a warning). |
| `--include-metadata` |  | FLAG | off | Include additional metadata columns (genome_size, etc.) |
| `--verbose` | `-v` | FLAG | off | Show detailed breakdown by genus |
| `--quiet` | `-q` | FLAG | off | Suppress progress output |
| `--dry-run` |  | FLAG | off | Show what would be queried without making requests |

### `metadarkmatter download genomes fetch`

Download genome FASTAs from NCBI using accession list.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--accessions` | `-a` | PATH | - | TSV file with accession list (from 'list' command) (required) |
| `--output-dir` | `-o` | PATH | - | Output directory for genome FASTA files (required) |
| `--include-protein` |  | FLAG | off | Also download protein FASTA files (.faa) for AAI computation |
| `--skip-if-exists/--redownload` |  | FLAG | on | Skip already downloaded genomes (default) or redownload all |
| `--decompress/--keep-compressed` |  | FLAG | on | Decompress downloaded .gz files (default) or keep compressed |
| `--batch-size` |  | INT | 100 | Number of accessions per NCBI request (max 1000) |
| `--verbose` | `-v` | FLAG | off | Show detailed progress |
| `--quiet` | `-q` | FLAG | off | Suppress progress output |
| `--dry-run` |  | FLAG | off | Show what would be downloaded without making requests |
| `--require-all` |  | FLAG | off | Exit non-zero if any genome could not be downloaded (failures are listed in `download_failures.tsv`); by default partial downloads exit 0 with a summary |

### `metadarkmatter proteins predict`

Predict protein sequences from genome FASTA files using Prodigal.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--genomes` | `-g` | PATH | - | Directory containing genome FASTA files (.fna) (required) |
| `--output-dir` | `-o` | PATH | - | Output directory for protein files (default: same as genomes) |
| `--missing-only` | `-m` | FLAG | off | Only predict proteins for genomes without existing .faa files |
| `--threads` | `-t` | INT | 4 | Number of parallel Prodigal processes |
| `--procedure` | `-p` | TEXT | single | Prodigal procedure: 'single' for isolate genomes, 'meta' for metagenomes |
| `--genome-pattern` |  | TEXT | *.fna | Glob pattern for genome files |
| `--verbose` | `-v` | FLAG | off | Show detailed progress for each genome |
| `--quiet` | `-q` | FLAG | off | Suppress progress output |
| `--dry-run` |  | FLAG | off | Show what would be predicted without running |

### `metadarkmatter map reads`

Competitively map reads against reference genomes.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--reads-1` | `-1` | PATH | - | Forward reads file (FASTQ, optionally gzipped) (required) |
| `--reads-2` | `-2` | PATH | - | Reverse reads file for paired-end data |
| `--genomes` | `-g` | PATH | - | Reference genome directory (FASTA files) or pangenome FASTA (required) |
| `--output` | `-o` | PATH | - | Output BAM file path (required) |
| `--index` | `-i` | PATH | - | Existing Bowtie2 index prefix (skip index building) |
| `--genome-pattern` |  | TEXT | *.fna | Glob pattern for genome files in directory |
| `--mode` | `-m` | TEXT | local | Alignment mode: 'local' or 'end-to-end' |
| `--very-sensitive/--fast` |  | FLAG | on | Use very-sensitive preset for divergent sequences (default: very-sensitive) |
| `--max-alignments` | `-k` | INT | 500 | Report up to this many alignments per read (high value captures all competing hits) |
| `--no-unal/--include-unal` |  | FLAG | on | Exclude unaligned reads from output (default: exclude) |
| `--no-discordant` |  | FLAG | off | Do not report discordant paired-end alignments |
| `--skip-if-exists` |  | FLAG | off | Skip processing if output BAM already exists |
| `--threads` | `-p` | INT | 4 | Number of threads |
| `--verbose` | `-v` | FLAG | off | Enable verbose output |
| `--quiet` | `-q` | FLAG | off | Suppress progress output |
| `--dry-run` |  | FLAG | off | Show commands without executing |

### `metadarkmatter tree build`

Build a phylogenetic tree.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--method` | `-m` | CHOICE | - | Tree building method: nj, upgma, or mashtree Choices: nj, upgma, mashtree. (required) |
| `--ani` | `-a` | PATH | - | ANI matrix CSV file (required for nj/upgma methods) |
| `--genomes` | `-g` | PATH | - | Directory containing genome FASTA files (required for mashtree method) |
| `--output` | `-o` | PATH | - | Output Newick tree file (required) |
| `--threads` | `-t` | INT | 4 | Number of threads (mashtree only) |
| `--genome-pattern` | `-p` | TEXT | *.fna | Glob pattern for genome files (mashtree only) |
| `--verbose` | `-v` | FLAG | off | Enable verbose output |
| `--quiet` | `-q` | FLAG | off | Suppress progress output |

### `metadarkmatter validate ani`

Validate ANI matrix format and values.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--input` | `-i` | PATH | - | ANI matrix CSV (required) |

### `metadarkmatter validate blast`

Validate BLAST tabular output format (outfmt 6). Accepts 12, 13 (+qlen) or 14 (+qlen, qcov) columns, as `score classify` does; all sampled lines must have the same column count.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--input` | `-i` | PATH | - | BLAST tabular output file (.tsv or .tsv.gz) (required) |

### `metadarkmatter validate classifications`

Validate classification output file.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--input` | `-i` | PATH | - | Classification output CSV (required) |

### `metadarkmatter validate fasta`

Validate a FASTA file (streaming, dependency-free).

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--input` | `-i` | PATH | - | FASTA file (.fa/.fasta/.fna, optionally .gz) (required) |
| `--max-records` |  | INT | 0 | Stop after this many records (0 = read whole file). |
| `--sequence-type` | `-t` | TEXT | nucleotide | Allowed-character set: 'nucleotide' or 'protein'. |

### `metadarkmatter validate fastq`

Validate a FASTQ file (4-line record structure).

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--input` | `-i` | PATH | - | FASTQ file (.fq/.fastq, optionally .gz) (required) |
| `--max-records` |  | INT | 0 | Stop after this many records (0 = read whole file). |

### `metadarkmatter visualize recruitment`

Create recruitment plot from BAM alignment data.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--bam` | `-b` | PATH | - | Input BAM file from competitive mapping (required) |
| `--output` | `-o` | PATH | - | Output plot file path (required) |
| `--format` | `-f` | TEXT | html | Output format: 'html' (interactive), 'png' (static), or 'json' (data) |
| `--genome` | `-g` | TEXT | - | Specific genome(s) to plot (can be repeated) |
| `--min-identity` |  | FLOAT | 70.0 | Minimum percent identity to include |
| `--min-mapq` |  | INT | 0 | Minimum mapping quality to include |
| `--show-bands/--hide-bands` |  | FLAG | on | Show identity threshold bands |
| `--known-threshold` |  | FLOAT | 96.0 | Identity threshold for known species (100 - novelty_known_max) |
| `--novel-species-threshold` |  | FLOAT | 80.0 | Identity threshold for novel species (100 - novelty_novel_species_max) |
| `--novel-genus-threshold` |  | FLOAT | 75.0 | Identity threshold for novel genus |
| `--max-points` |  | INT | 100000 | Maximum points to plot (subsamples if exceeded) |
| `--width` |  | INT | 1200 | Plot width in pixels |
| `--height` |  | INT | 800 | Plot height in pixels |
| `--export-tsv` |  | PATH | - | Export data as TSV for anvi'o compatibility |
| `--multi-panel` |  | FLAG | off | Create multi-panel plot with one subplot per genome |
| `--verbose` | `-v` | FLAG | off | Enable verbose output |
| `--quiet` | `-q` | FLAG | off | Suppress progress output |

### `metadarkmatter visualize summary`

Create summary visualization from classification results.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--classifications` | `-c` | PATH | - | Classification results file (CSV or Parquet) (required) |
| `--output` | `-o` | PATH | - | Output plot file path (required) |
| `--format` | `-f` | TEXT | html | Output format: 'html' or 'png' |
| `--width` |  | INT | 1000 | Plot width in pixels |
| `--height` |  | INT | 600 | Plot height in pixels |
| `--verbose` | `-v` | FLAG | off | Enable verbose output |
| `--quiet` | `-q` | FLAG | off | Suppress progress output |

### `metadarkmatter run pipeline`

Run the complete pipeline from reads (or extracted reads) to an HTML report.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--config` | `-c` | PATH | - | Pipeline config YAML (CLI flags override its values). |
| `--family` |  | TEXT | - | GTDB family to download references for (e.g. f__Francisellaceae). |
| `--family-taxid` |  | INT | - | NCBI taxid for the family (needed for Kraken2 extraction). |
| `--reads-1` | `-1` | PATH | - | R1 reads (FASTQ[.gz]). |
| `--reads-2` | `-2` | PATH | - | R2 reads for paired-end input. |
| `--kraken-db` |  | PATH | - | Kraken2 database directory. |
| `--genomes` | `-g` | PATH | - | Pre-downloaded reference genomes dir (skips download). |
| `--kraken-output` |  | PATH | - | Pre-computed .kraken file (skips kraken classify). |
| `--extracted-reads-1` |  | PATH | - | Pre-extracted family R1 reads (skips both kraken steps). |
| `--extracted-reads-2` |  | PATH | - | Pre-extracted family R2 reads (paired-end; both aligners align both mates). |
| `--output-dir` | `-o` | PATH | - | Root output directory for the run (required unless set in --config). |
| `--aligner` |  | CHOICE | - | Aligner for the alignment step [default: blast]. Choices: blast, mmseqs2. |
| `--ani-backend` |  | TEXT | - | ANI backend: auto, fastani or skani [default: auto, or the --config value]. |
| `--no-report` |  | FLAG | off | Skip the HTML report step. |
| `--from` |  | CHOICE | - | Resume from this step (earlier outputs must already exist). Choices: download, kraken-classify, kraken-extract, makedb, align, ani, classify, report. |
| `--force` |  | FLAG | off | Re-run every step even if its outputs already exist. |
| `--threads` | `-t` | INT | - | Threads passed to each tool step [default: 8, or the --config value]. |
| `--quiet` | `-q` | FLAG | off | Suppress progress output (useful for scripting). |
| `--dry-run` |  | FLAG | off | Show what would run without executing. Honours the global --dry-run too. |

`download genomes list` writes both the accession TSV given by `--output` and `genome_metadata.tsv` (accession, species, genus, family, representative, gtdb_taxonomy) next to it.

For `run pipeline`, values in the `--config` YAML (e.g. aligner, threads, output_dir) apply unless the corresponding flag is passed explicitly. `--from STEP` re-runs that step and all later steps even when their outputs exist. The run manifest records start and finish times for each step. With either aligner, R2 reads (`--reads-2` / `--extracted-reads-2`) are aligned together with R1 (mates get `/1` and `/2` ID suffixes). Unknown keys in the `--config` YAML produce a warning.

**Checkpointing.** A step is skipped on re-run only if a completion marker (`.mdm_steps/<step>.done.json`, written after the step exits 0) exists and its fingerprint (the step's command line plus size and modification time of its input files) matches. A changed fingerprint, or a step that was interrupted, re-runs that step and every later step. The previous `run_manifest.json` is kept as `run_manifest.<timestamp>.json`.


### `metadarkmatter run cohort`

Run the pipeline over many samples that share one reference.

```bash
# samples.tsv (tab-separated; relative paths resolve against the sheet):
# sample    reads_1         reads_2
# S1        S1_R1.fq.gz     S1_R2.fq.gz
# S2        S2_R1.fq.gz     S2_R2.fq.gz
mdm run cohort --samples samples.tsv --family f__Francisellaceae --family-taxid 34064 \
    --kraken-db db/ --output-dir results/cohort01
```

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--samples` | `-s` | PATH | required | Sample sheet TSV: `sample`, `reads_1`, optional `reads_2`; or `extracted_reads_1`, optional `extracted_reads_2` (skips Kraken2). `#` lines are ignored. |
| `--output-dir` | `-o` | PATH | required | Cohort root directory. |
| `--family`, `--family-taxid`, `--kraken-db`, `--genomes` | | | - | As for `run pipeline`. |
| `--aligner` | | CHOICE | blast | blast or mmseqs2. |
| `--ani-backend` | | TEXT | auto | auto, fastani or skani. |
| `--no-report` | | FLAG | off | Skip per-sample reports and the cohort report. |
| `--no-sample-reports` | | FLAG | off | Skip only the per-sample HTML reports. |
| `--from` | | CHOICE | - | Resume from this step, for the reference and for every sample. |
| `--force` | | FLAG | off | Re-run every step. |
| `--threads` | `-t` | INT | 8 | Threads per tool step. |
| `--quiet` | `-q` | FLAG | off | Suppress progress output. |
| `--dry-run` | | FLAG | off | Plan all steps for the reference and every sample without executing (global `--dry-run` also honoured). |

Layout:

```
cohort01/
  reference/               genomes, genome_metadata.tsv, blastdb/ (or mmseqs_db/), ani_matrix.csv
  samples/<sample>/        <sample>.kraken, extraction/, <sample>.blast.tsv.gz,
                           classifications.csv, summary.json, reports/report.html,
                           run_manifest.json, .mdm_steps/
  cohort_summary.tsv       one row per sample
  cohort_report.html       multi-sample report (via `report multi`)
  cohort_manifest.json     per-sample status
  multi_inputs/            symlinks to each completed sample's classifications.csv
```

Sample names (sheet column `sample`, also used as the Kraken2 `--sample-name` and as file stems) must be unique, start with a letter or digit, and contain only letters, digits, `_`, `.` and `-`; names containing `_classifications`, `.kraken`, `.blast` or `.tsv` are rejected because `report multi` would alter them. All read files must exist. Samples are processed sequentially in sheet order, each with the same per-step fingerprint checkpointing as `run pipeline`, so re-running the command only redoes steps whose inputs changed. A failing reference build stops the run; a failing sample does not stop the other samples, but the command exits with status 1 and the sample is marked `failed` (with the error) in `cohort_summary.tsv`.

`cohort_summary.tsv` columns: `sample`, `status` (`completed` or `failed`), `total_reads`, one count column per `taxonomic_call` category present (sorted by category name), `novel_fraction` ((Novel Species + Novel Genus) / total reads), `error`.

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

**Note:** The `qlen` column (column 13) is used for coverage-weighted hit selection and `--min-alignment-fraction`. Older 12-column output (without qlen) is still supported - `qend` is used as a proxy for read length in this case. Hits with `qstart > qend` (minus-strand hits as reported by MMseqs2 and DIAMOND) are normalised at parse time.

### ANI Matrix Format

CSV or TSV with genome names as first column and header row:

```csv
genome,GCF_000123456.1,GCF_000789012.1,GCA_000111222.1
GCF_000123456.1,100.0,95.5,82.3
GCF_000789012.1,95.5,100.0,83.1
GCA_000111222.1,82.3,83.1,100.0
```

Values must be percent identity on a 0-100 scale; a matrix whose values are all <= 1.0 is rejected with an error. Pairs present in only one direction are mirrored. Missing pairs (0) are treated as the default unrelated ANI (70%) for lookups. Forward/reverse disagreement above 0.5 produces a warning (an error with `--strict-ani`).

---

## Output File Formats

### Classification Output (CSV/Parquet)

Columns written by `score classify` (in this order):

| Column | Type | Description |
|--------|------|-------------|
| read_id | string | Query sequence identifier |
| best_match_genome | string | Top-scoring genome |
| top_hit_identity | float | Percent identity of the best hit |
| novelty_index | float | 100 - top_hit_identity |
| placement_uncertainty | float | 100 - ANI to competing genome(s), per `--uncertainty-mode`; null for Off-target reads |
| genus_uncertainty | float | Based on minimum ANI to secondary genomes |
| ambiguity_scope | string | `unambiguous`, `within_species`, `within_genus`, or `across_genera` |
| num_ambiguous_hits | int | Hits within the bitscore threshold of the top hit |
| second_hit_identity | float | Identity of the best hit to a different genome |
| identity_gap | float | top_hit_identity - second_hit_identity |
| confidence_score | float | `(1 - posterior_entropy / log2(6)) * 100`; 10.0 for Off-target reads |
| taxonomic_call | string | Bayesian MAP category after Stage 2 refinement (primary call) |
| legacy_call | string | Hard-threshold cascade result |
| diversity_status | string | Known, Novel, or Uncertain |
| is_novel | boolean | True if taxonomic_call is Novel Species or Novel Genus |
| low_confidence | boolean | Low-confidence flag |
| p_known_species, p_novel_species, p_novel_genus, p_species_boundary, p_ambiguous, p_unclassified | float | Posterior probabilities (sum to 1) |
| posterior_entropy | float | Shannon entropy of the posteriors (0 to 2.585) |
| inferred_uncertainty | float | Inferred uncertainty for single-hit reads |
| uncertainty_type | string | `measured`, `inferred` (single-hit reads), or `off_target` |

When family validation is active, `family_bitscore_ratio`, `family_identity_gap`, `in_family_hit_fraction`, `external_best_genome`, `external_best_identity`, `in_family_identity` and `in_family_novelty_index` are added. With `--metadata` (in-memory mode), `species`, `genus` and `family` are appended. `--include-legacy-scores` adds `alignment_quality`, `identity_confidence`, `placement_confidence` and `discovery_score`. In nucleotide mode the clade-novelty and single-hit novelty rules (both on by default) add the boolean columns `clade_novel_species` and `single_hit_novel_species`; each is absent when its rule is disabled and in protein mode.

**Taxonomic Call Values:** `Known Species`, `Novel Species`, `Novel Genus`, `Species Boundary`, `Ambiguous`, `Ambiguous Within Genus`, `Conserved Region`, `Unclassified`, and `Off-target`. See the Alignment Mode section above for how the hard thresholds behind `legacy_call` relate to the Bayesian `taxonomic_call`.

### Summary Output (JSON)

```json
{
  "total_reads": 10000,
  "known_species": 6500,
  "novel_species": 2000,
  "novel_genus": 500,
  "species_boundary": 300,
  "ambiguous": 400,
  "ambiguous_within_genus": 50,
  "conserved_regions": 200,
  "unclassified": 50,
  "off_target": 0,
  "diversity_known": 6500,
  "diversity_novel": 2500,
  "diversity_uncertain": 1000,
  "mean_novelty_index": 3.5,
  "mean_placement_uncertainty": 1.2,
  "genome_hit_counts": {"GCF_000123456.1": 3500, "GCF_000789012.1": 2800},
  "species_hit_counts": {},
  "diversity_known_pct": 65.0,
  "diversity_novel_pct": 25.0,
  "diversity_uncertain_pct": 10.0,
  "known_species_pct": 65.0,
  "novel_species_pct": 20.0,
  "novel_genus_pct": 5.0,
  "novel_diversity_pct": 25.0,
  "species_count": 0
}
```

Counts are taken from `taxonomic_call`; the `diversity_*` fields group them by the per-read `diversity_status` (Known = Known Species; Novel = Novel Species + Novel Genus; Uncertain = the remaining in-family categories). `species_hit_counts` and `species_count` are filled only when `--metadata` is given and are empty with `--streaming`.

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
| `MDM_DRY_RUN` | Set to `1` for a process-wide dry run (same as global `--dry-run`) | unset |
| `MDM_DEBUG` | Set to `1` to print full tracebacks (same as global `--debug`) | unset |
| `METADARKMATTER_SEED` | Seed for every randomness-using step (plot subsampling, adaptive GMM fits) | 42 |
| `METADARKMATTER_GTDB_CACHE_DIR` | On-disk cache for GTDB API responses (7-day TTL); empty disables | unset |
| `POLARS_MAX_THREADS` | Maximum threads for Polars (read by Polars itself) | CPU count |

External tools (BLAST, MMseqs2, Kraken2, fastANI, etc.) run without a default timeout; a timeout, when set programmatically, is reported as a tool error.

---

## Performance Guide

### Classifier Selection

| Scenario | Recommended Mode | Command |
|----------|-----------------|---------|
| Up to ~5M alignments | Default (VectorizedClassifier) | (default) |
| > ~5M alignments or limited RAM | Streaming | `--streaming --chunk-size 500000` |

All classification uses the VectorizedClassifier by default, which is Polars-based and auto-parallelized. The `--streaming` mode is available for large files where memory is a constraint.

### Memory Usage

In-memory (vectorized) mode holds the whole alignment table. Before running it estimates its memory use as about 1.5x the uncompressed alignment size (about 15x for a gzipped file) and warns if that may exceed available memory; in the [SRR25038281 case study](CASE_STUDIES/SRR25038281_chesapeake_cb52.md) this estimate was about 18 GB for a 1.3 GB gzipped BLAST file. Streaming mode memory is bounded by `--chunk-size` (default 1,000,000 rows). These figures are rough guidance, not measurements.

---

## Troubleshooting

### Common Errors

**"No BLAST files found matching pattern"**
- Check the `--pattern` option matches your files
- Verify the `--alignment-dir` path is correct

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
metadarkmatter score classify --alignment file.tsv --ani matrix.csv --output out.csv --dry-run
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
| `--ani` | `-a` | PATH | None | ANI matrix file for heatmap visualization and the phylogeny tree |
| `--aai` | | PATH | None | AAI matrix file (CSV) for genus-level heatmap |
| `--metadata` | `-m` | PATH | None | Genome metadata TSV for species-level breakdown |
| `--recruitment` | `-r` | PATH | None | Recruitment data file for recruitment plots |
| `--bam` | `-b` | PATH | None | BAM file for generating recruitment data on-the-fly |
| `--theme` | | TEXT | light | Color theme: 'light' or 'dark' |
| `--alignment-mode` | | TEXT | nucleotide | 'nucleotide' (uses ANI) or 'protein' (uses AAI) |
| `--max-points` | | INT | 50000 | Maximum points in scatter plots (for performance) |
| `--max-table-rows` | | INT | no limit | Maximum rows embedded in the Data tab table. When set, a deterministic sample balanced across `taxonomic_call` categories is kept (rare categories in full; highest novelty first within a category; original order preserved) and the report states that the table is truncated. Summary counts and plots always use all reads |
| `--provenance` | | PATH | `<classifications>.provenance.json` if present | Provenance sidecar written by `score classify`; feeds the Summary tab's Run information panel |
| `--max-phylo-clusters` | | INT | 20 | Maximum novel clusters in the phylogenetic context heatmap |
| `--max-phylo-references` | | INT | 50 | Maximum reference genomes in the phylogenetic context heatmap |
| `--tree` | | PATH | None | Newick tree file; if omitted, a neighbor-joining tree is built from the ANI matrix |
| `--no-phylogeny` | | FLAG | False | Skip the phylogeny view (faster for large datasets) |
| `--cluster-grouping` | | TEXT | auto | How the Novel Diversity tab groups novel reads: `auto` (genus when metadata gives genus names, otherwise reference), `genus` (GTDB genus of the nearest reference; within a genus, references below 80% ANI are kept apart), `reference` (nearest reference) or `legacy` (ANI neighbourhood x 5% novelty band x call, the previous method). See METHODS.md section 4.1 and NOVEL_LINEAGE_GROUPING.md. |
| `--report-mode` | | TEXT | offline | 'offline' inlines Plotly (and vendored D3) for a self-contained file; 'cdn' uses remote script tags |
| `--verbose` | `-v` | FLAG | False | Enable verbose output |
| `--quiet` | `-q` | FLAG | False | Suppress progress output |

#### Report Sections

The generated HTML report has the following tabs, plus a collapsible Methods panel below them:

| Section | Contents |
|---------|----------|
| **Summary** | Metric cards, key findings, top species (with metadata) |
| **Classification** | Category charts, novelty/uncertainty distributions and scatter, Bayesian confidence panels |
| **Novel Diversity** | Novel clusters and phylogenetic context (shown only when novel reads are present) |
| **Reference** | ANI/AAI heatmaps, recruitment plots (if BAM/data provided), species bar chart (with metadata) |
| **Data** | Interactive sortable/filterable classification table with export |

Off-target rows are handled; null values are shown as "-" in the table. Apart from the generation timestamp, a report is byte-identical across runs with the same inputs.

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

### `metadarkmatter report summarize`

Summarize novel diversity into novel lineage groups (one row per group).

Positional argument: `CLASSIFICATIONS`

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--output` | `-o` | PATH | - | Output file path (TSV or JSON based on --format) (required) |
| `--metadata` | `-m` | PATH | - | Genome metadata file (TSV) for species/genus/family lookups |
| `--band-size` | `-b` | FLOAT | 5.0 | Size of novelty bands in percent, a whole number (`--grouping legacy` only) |
| `--min-cluster-size` | `-s` | INT | 3 | Minimum reads to form a group (default 3) |
| `--ani` | `-a` | PATH | - | ANI matrix (CSV). With genus grouping, keeps references of one genus below 80% ANI apart; with legacy grouping, enables neighbourhood merging |
| `--grouping` | `-g` | TEXT | auto | `auto` (genus when metadata gives genus names, otherwise reference), `genus` (GTDB genus of the nearest reference; within a genus, references below 80% ANI are kept apart), `reference` (nearest reference) or `legacy` (ANI neighbourhood x 5% novelty band x call, the previous method). See METHODS.md section 4.1 and NOVEL_LINEAGE_GROUPING.md. |
| `--format` | `-f` | TEXT | tsv | Output format: 'tsv' or 'json' |
| `--json` | `-j` | PATH | - | Additional JSON output path with full metadata |
| `--verbose` | `-v` | FLAG | off | Enable verbose output |
| `--quiet` | `-q` | FLAG | off | Suppress progress output |

### `metadarkmatter report vendor-d3`

Download d3.v7.min.js into the package's report assets directory.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--force` |  | FLAG | off | Overwrite the cached D3 asset if it already exists. |
| `--timeout` |  | FLOAT | 30.0 | HTTP timeout in seconds. |

### `metadarkmatter report vendor-cdn-sri`

Download and hash the CDN assets so the report can attach SRI integrity.

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--timeout` |  | FLOAT | 30.0 | HTTP timeout in seconds. |

---

## Report Output Details

### HTML Report Structure

The generated HTML report is self-contained with:

- **Embedded CSS**: Styling with light/dark themes
- **Plotly.js**: Inlined by default (`--report-mode offline`); loaded from a CDN with `--report-mode cdn` (with SRI attributes after `report vendor-cdn-sri`)
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
