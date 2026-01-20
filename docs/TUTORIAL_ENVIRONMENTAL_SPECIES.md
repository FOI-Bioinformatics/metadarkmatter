# Tutorial: Discovering Environmental Species in Your Favorite Bacterial Family

A practical guide to using metadarkmatter for identifying novel bacterial diversity in metagenomic samples.

---

## Quick Start

For experienced users, here is the complete workflow:

```bash
# Setup
mkdir my_analysis && cd my_analysis
FAMILY="Francisellaceae"
TAXID=34064
SAMPLE="water_sample"
THREADS=16

# 1. Download reference genomes (auto-creates genome_metadata.tsv)
metadarkmatter download genomes list "f__${FAMILY}" --output genomes.tsv
metadarkmatter download genomes fetch --accessions genomes.tsv --output-dir genomes/ --include-protein

# 1b. Generate missing protein files with Prodigal (if any)
metadarkmatter proteins predict --genomes genomes/ --missing-only --threads $THREADS

# 2. Build ANI matrix
metadarkmatter ani compute --genomes genomes/ --output ani_matrix.csv --threads $THREADS

# 2b. Build AAI matrix (optional, for genus-level disambiguation)
metadarkmatter aai compute --genomes genomes/ --output aai_matrix.csv --threads $THREADS

# 3. Build BLAST database (standardizes headers for multi-contig genomes)
metadarkmatter blast makedb --genomes genomes/ --output blastdb/pangenome

# 4. Classify metagenome and extract family reads
metadarkmatter kraken2 classify \
  --reads-1 ${SAMPLE}_R1.fastq.gz --reads-2 ${SAMPLE}_R2.fastq.gz \
  --kraken-db /path/to/kraken2_db --output kraken_out/ --threads $THREADS

metadarkmatter kraken2 extract \
  --kraken-output kraken_out/${SAMPLE}.kraken \
  --reads-1 ${SAMPLE}_R1.fastq.gz --reads-2 ${SAMPLE}_R2.fastq.gz \
  --taxid $TAXID --output extracted/

# 5. Sequence alignment with BLAST
metadarkmatter blast align \
  --query extracted/${SAMPLE}_taxid${TAXID}_R1.fastq.gz \
  --database blastdb/pangenome --output ${SAMPLE}.blast.tsv.gz --threads $THREADS

# 6. ANI-weighted classification (with species-level tracking and AAI)
metadarkmatter score classify \
  --blast ${SAMPLE}.blast.tsv.gz --ani ani_matrix.csv --aai aai_matrix.csv \
  --metadata genome_metadata.tsv \
  --output classifications.csv --summary summary.json --parallel

# 7. Extract novel candidates and generate report (with species breakdown)
metadarkmatter score extract-novel \
  --classifications classifications.csv \
  --output novel_candidates.csv --read-ids novel_reads.txt

metadarkmatter report generate \
  --classifications classifications.csv \
  --metadata genome_metadata.tsv \
  --ani ani_matrix.csv --aai aai_matrix.csv \
  --output report.html
```

---

## Introduction

### The Problem

Standard taxonomic classifiers like Kraken2 assign metagenomic reads to known taxa, but they cannot tell you *how similar* those reads are to their assigned references. A read classified as *Pseudomonas* might be:

- **Nearly identical** to a sequenced strain (known species)
- **Moderately divergent** representing a potentially new species
- **Highly divergent** representing a potentially new genus
- **Misclassified** actually belonging to a different family

This distinction matters for ecological studies. Environmental samples often contain organisms that have never been cultured or sequenced - the so-called "microbial dark matter."

### The Solution

Metadarkmatter addresses this by calculating two metrics for each read:

| Metric | Definition | What It Tells You |
|--------|------------|-------------------|
| **Novelty Index** | 100 - percent identity to best match | How divergent the read is from known references |
| **Placement Uncertainty** | 100 - max ANI between competing matches | Whether the read could belong to multiple taxa |

These metrics enable classification into biologically meaningful categories:

| Category | Novelty | Uncertainty | Interpretation |
|----------|---------|-------------|----------------|
| Known Species | < 5% | < 2% | Matches a sequenced species |
| Novel Species | 5-15% | < 2% | Likely a new species within a known genus |
| Novel Genus | 15-25% | < 2% | Likely a new genus within the family |
| Conserved Region | any | > 5% | Conserved gene, placement ambiguous |

### Threshold Basis

The thresholds are based on established prokaryotic taxonomy:
- **95-96% ANI** defines the species boundary (Jain et al. 2018)
- **85% identity** approximates genus-level divergence
- **75% identity** approximates family-level divergence

---

## Prerequisites

### Software Installation

```bash
# Python 3.11+ required
python --version

# Install metadarkmatter
pip install metadarkmatter
# Or from source: pip install -e /path/to/metadarkmatter
```

### External Software Dependencies

Metadarkmatter requires several external bioinformatics tools:

**Required tools** (install all):
```bash
conda install -c bioconda kraken2 krakentools blast skani
```

| Tool | Purpose |
|------|---------|
| Kraken2 | Taxonomic classification of reads |
| KrakenTools | Extract reads by taxid |
| BLAST+ | Nucleotide alignment (blastn) |
| skani | Fast ANI computation |

**Optional tools**:
```bash
# MMseqs2 - only for datasets >100K reads
conda install -c bioconda mmseqs2

# seqtk - only for assembly workflows
conda install -c bioconda seqtk

# fastANI - alternative to skani for ANI computation
conda install -c bioconda fastani
```

**Note:** BLAST accepts FASTQ files directly (automatic conversion), so seqtk is only needed for specialized workflows.

```bash

### Data Requirements

| Requirement | Description |
|-------------|-------------|
| Metagenomic reads | Paired-end FASTQ files (gzipped OK) |
| Kraken2 database | Standard, PlusPF, or custom database |
| Target family TaxID | NCBI Taxonomy ID (find at ncbi.nlm.nih.gov/taxonomy) |

### Hardware Recommendations

| Dataset Size | RAM | Disk | Runtime (BLAST) |
|--------------|-----|------|-----------------|
| < 1M reads | 16 GB | 20 GB | 1-2 hours |
| 1-10M reads | 32 GB | 100 GB | 4-8 hours |
| 10-100M reads | 64 GB | 500 GB | 1-3 days |

**Note:** For very large datasets (>100K reads), see the MMseqs2 section in Step 7 for faster alignment options.

---

## Step 1: Choose Your Target Family

Select a bacterial family based on:

**Scientific criteria:**
- Ecological relevance to your sample type
- Known to harbor uncultured diversity
- Contains organisms of medical, agricultural, or industrial interest

**Practical criteria:**
- Has sufficient reference genomes in GTDB (ideally > 30)
- Is represented in your Kraken2 database

### Example Families

| Family | TaxID | Typical Genomes | Notes |
|--------|-------|-----------------|-------|
| Francisellaceae | 34064 | ~80 | Pathogens + environmental |
| Pseudomonadaceae | 135621 | ~2000 | Very diverse, many environmental |
| Methylococcaceae | 135618 | ~150 | Methanotrophs |
| Nitrospiraceae | 1760 | ~50 | Nitrite oxidizers |

For this tutorial, we use **Francisellaceae** as an example.

---

## Step 2: Set Up Your Analysis

Create a project directory and define variables:

```bash
# Create project structure
mkdir -p francisella_analysis
cd francisella_analysis

# Define analysis parameters
FAMILY="Francisellaceae"
TAXID=34064
SAMPLE="env_sample_01"
THREADS=16
KRAKEN_DB="/path/to/kraken2_standard_db"
```

---

## Step 3: Download Reference Genomes

### 3a. Query GTDB for Available Genomes

```bash
metadarkmatter download genomes list "f__${FAMILY}" \
  --output genomes_list.tsv
```

This command does two things:
1. Creates `genomes_list.tsv` with accession numbers for downloading
2. **Automatically creates `genome_metadata.tsv`** with species-level taxonomy

Examine the results:

```bash
# Count genomes
wc -l genomes_list.tsv
# 82 genomes_list.tsv

# Preview accession list
head -5 genomes_list.tsv
```

Expected output:

```
accession        species                      genus           gtdb_taxonomy
GCF_000195955.2  s__Francisella tularensis    g__Francisella  d__Bacteria;p__Proteobacteria;...
GCF_000242755.1  s__Francisella philomiragia  g__Francisella  d__Bacteria;p__Proteobacteria;...
GCF_000297895.1  s__Francisella noatunensis   g__Francisella  d__Bacteria;p__Proteobacteria;...
GCF_001971545.1  s__Francisella opportunistica g__Francisella d__Bacteria;p__Proteobacteria;...
```

Check the metadata file:

```bash
# Preview genome metadata (used later for species-level tracking)
head -3 genome_metadata.tsv
```

```
accession	species	genus	family	gtdb_taxonomy
GCF_000195955.2	Francisella tularensis	Francisella	Francisellaceae	d__Bacteria;p__Proteobacteria;...
GCF_000242755.1	Francisella philomiragia	Francisella	Francisellaceae	d__Bacteria;p__Proteobacteria;...
```

This metadata file links genome accessions to species names and is used throughout the pipeline for species-level aggregation.

### 3b. Download Genomes

```bash
metadarkmatter download genomes fetch \
  --accessions genomes_list.tsv \
  --output-dir genomes/ \
  --include-protein
```

The `--include-protein` flag downloads both genome nucleotide files (.fna) and protein files (.faa) from NCBI. This enables AAI computation without needing a gene prediction step.

For 80 genomes, expect 5-15 minutes.

```bash
# Verify downloads
ls genomes/*.fna | wc -l
# 82

ls genomes/*.faa | wc -l
# 82

ls -lh genomes/ | head -5
# -rw-r--r--  1 user  staff  1.9M  GCF_000195955.2.fna
# -rw-r--r--  1 user  staff  589K  GCF_000195955.2.faa
# -rw-r--r--  1 user  staff  2.1M  GCF_000242755.1.fna
# -rw-r--r--  1 user  staff  645K  GCF_000242755.1.faa
```

**Without AAI:**
If you don't plan to use AAI, you can omit the `--include-protein` flag to save disk space and download time:

```bash
metadarkmatter download genomes fetch \
  --accessions genomes_list.tsv \
  --output-dir genomes/
```

**Missing Protein Files:**
NCBI doesn't provide protein annotations for all genomes (typically 5-15% are missing). This is common for:
- GenBank-only assemblies (no RefSeq equivalent)
- Draft genomes without automated annotation
- Older assemblies predating annotation pipelines

If you get a warning about missing protein files:

```bash
Warning: 5 genomes missing protein files
This is common for GenBank-only assemblies or draft genomes.

Missing protein files for:
  - GCA_027621155.1      Francisella halioticida
  - GCA_038141505.1      Francisella sp. nov.
  ...
```

Use the built-in Prodigal integration to fill in the gaps:

```bash
# Generate protein files for missing genomes only (recommended)
metadarkmatter proteins predict \
  --genomes genomes/ \
  --missing-only \
  --threads 8
```

Or use the shell script approach:

```bash
cd genomes/
for fna in *.fna; do
  base=$(basename $fna .fna)
  [ ! -f ${base}.faa ] && prodigal -i $fna -a ${base}.faa -q
done
```

Install Prodigal if needed:
```bash
conda install -c bioconda prodigal
```

**Alternative: Predict proteins for ALL genomes:**
If NCBI annotations are inconsistent or you prefer uniform annotations:

```bash
metadarkmatter proteins predict \
  --genomes genomes/ \
  --threads 8
```

### Tips for Large Families

For families with > 500 genomes (e.g., Pseudomonadaceae), download only representative species:

```bash
metadarkmatter download genomes list "f__Pseudomonadaceae" \
  --output genomes_list.tsv \
  --representative-only
```

---

## Step 4: Build the ANI Matrix

ANI (Average Nucleotide Identity) quantifies how similar two genomes are. The ANI matrix is essential for calculating placement uncertainty.

```bash
metadarkmatter ani compute \
  --genomes genomes/ \
  --output ani_matrix.csv \
  --threads $THREADS
```

The command:
- Auto-detects the best available tool (prefers skani over fastANI)
- Computes all pairwise ANI values
- Outputs a symmetric matrix in CSV format

Expected output:

```
Metadarkmatter ANI Matrix Builder

Step 1: Detecting ANI tool...
  Using: skani (faster, recommended)

Step 2: Computing pairwise ANI...
  Genomes: 82
  Comparisons: 3,321
  Completed in 2m 34s

Step 3: Building matrix...
  Matrix saved to ani_matrix.csv
```

### Backend Selection

You can explicitly choose the ANI backend:

```bash
# Use skani (faster, lower memory)
metadarkmatter ani compute --genomes genomes/ --output ani_matrix.csv --backend skani

# Use fastANI (original tool, widely validated)
metadarkmatter ani compute --genomes genomes/ --output ani_matrix.csv --backend fastani
```

| Backend | Speed | Best For |
|---------|-------|----------|
| skani | ~5x faster | Large genome sets (> 100) |
| fastANI | Moderate | Publication validation, smaller sets |

### Validate the Matrix

```bash
metadarkmatter ani validate \
  --ani ani_matrix.csv \
  --genomes genomes/
```

This ensures all genomes are represented in the matrix.

---

## Step 4b: Build the AAI Matrix (Optional but Recommended)

AAI (Average Amino Acid Identity) provides genus-level resolution for reads with 20-25% novelty. It's particularly useful for:
- Samples with high genus-level diversity
- Disambiguating borderline genus/species reads
- Cross-genus validation

### When to Use AAI

| Scenario | Use AAI? | Reason |
|----------|----------|--------|
| Well-characterized family | Optional | Most reads are species-level matches |
| Novel diversity expected | **Recommended** | Helps distinguish novel species from novel genus |
| Cross-genus analysis | **Essential** | AAI > 65% confirms same genus |
| High ambiguity (>40%) | **Recommended** | May help resolve genus-level placement |

### Compute AAI Matrix

```bash
metadarkmatter aai compute \
  --genomes genomes/ \
  --output aai_matrix.csv \
  --threads $THREADS
```

The command:
- Extracts protein sequences from all genomes
- Runs Diamond BLASTP for all-vs-all protein alignment
- Computes reciprocal best hits (RBH) for ortholog identification
- Calculates mean AAI from RBH identities
- Outputs a symmetric matrix in CSV format

Expected output:

```
Metadarkmatter AAI Matrix Builder

Input: 82 protein files from genomes
Threads: 16
Parameters: min_identity=30.0%, min_coverage=0.5

Processing 82 genomes with 8 parallel workers...
Completed 82/82 genomes

AAI matrix computed successfully!

Output: aai_matrix.csv
Genomes: 82
Size: 0.15 MB

Matrix shape: 82 x 82
AAI range: 38.5% - 97.2%
AAI mean: 62.8%
```

### Performance

**With optimization (parallel processing):**
- 42 genomes: ~1 minute
- 82 genomes: ~2-3 minutes
- 200 genomes: ~5-10 minutes

**Thread allocation:**
- 8 threads = 4 genomes in parallel, 2 threads each
- 16 threads = 8 genomes in parallel, 2 threads each

### AAI Interpretation (Riesco & Trujillo 2024)

| AAI Range | Interpretation |
|-----------|---------------|
| > 65% | Same genus |
| 58-65% | Genus boundary zone |
| 45-58% | Different genus |
| < 45% | Highly divergent (family level) |

### Impact on Classification

With full AAI coverage, expect:
- **3-5% of genus-boundary reads** reclassified more accurately
- Novel Genus reads (20-25% novelty) refined to Novel Species when AAI > 65%
- Improved confidence scores for borderline classifications

**Example from validation:**
- 16 reads reclassified (3.3% of total)
- 8 reads: Novel Genus → Novel Species (AAI confirmed same genus)
- Novel Species count increased by 11.4%

---

## Step 5: Build the BLAST Database

Create a BLAST database from your reference genomes:

```bash
metadarkmatter blast makedb \
  --genomes genomes/ \
  --output blastdb/pangenome
```

This command:
1. Concatenates all genomes into a single FASTA file
2. **Standardizes FASTA headers** to format: `{accession}|{contig_id}`
3. Creates a contig mapping file for tracking multi-contig draft genomes
4. Builds the BLAST nucleotide database

```bash
# Verify
ls blastdb/
# pangenome.nhr  pangenome.nin  pangenome.nsq
# pangenome_pangenome.fasta  pangenome_contig_mapping.tsv
```

### Standardized Headers

Many reference genomes are draft assemblies with multiple contigs. The standardized header format ensures reliable genome identification from BLAST hits:

```bash
# Original header (varies by genome)
>NZ_CP000439.1 Francisella tularensis subsp. holarctica complete genome

# Standardized header (consistent format)
>GCF_000195955.2|NZ_CP000439.1
```

The contig mapping file tracks the relationship between contigs and genomes:

```bash
head -3 blastdb/pangenome_contig_mapping.tsv
```

```
contig_id	genome_accession	original_header
GCF_000195955.2|NZ_CP000439.1	GCF_000195955.2	NZ_CP000439.1 Francisella tularensis...
GCF_000242755.1|NZ_CP000937.1	GCF_000242755.1	NZ_CP000937.1 Francisella philomiragia...
```

---

## Step 6: Extract Family Reads from Metagenome

### 6a. Run Kraken2 Classification

```bash
metadarkmatter kraken2 classify \
  --reads-1 ${SAMPLE}_R1.fastq.gz \
  --reads-2 ${SAMPLE}_R2.fastq.gz \
  --kraken-db $KRAKEN_DB \
  --output kraken_output/ \
  --threads $THREADS
```

### 6b. Check Family Abundance

```bash
# Find your target family in the report
grep -i "francisella" kraken_output/${SAMPLE}.kreport
```

Example output:

```
  0.24   48521   1205   F       34064         Francisellaceae
  0.21   42316   8432   G       262             Francisella
  0.08   16234   16234  S       263               Francisella tularensis
  0.05   10150   10150  S       119857            Francisella philomiragia
```

This shows ~48,500 reads classified to Francisellaceae.

### 6c. Extract Family Reads

```bash
metadarkmatter kraken2 extract \
  --kraken-output kraken_output/${SAMPLE}.kraken \
  --reads-1 ${SAMPLE}_R1.fastq.gz \
  --reads-2 ${SAMPLE}_R2.fastq.gz \
  --taxid $TAXID \
  --output extraction/
```

Output (compressed by default):

```bash
ls -lh extraction/
# -rw-r--r--  1 user  staff  12M  env_sample_01_taxid34064_R1.fastq.gz
# -rw-r--r--  1 user  staff  13M  env_sample_01_taxid34064_R2.fastq.gz
```

---

## Step 7: Run Sequence Alignment

Align extracted reads against the reference database using BLAST.

### BLAST Alignment

Run BLAST to map reads to reference genomes:

```bash
metadarkmatter blast align \
  --query extraction/${SAMPLE}_taxid${TAXID}_R1.fastq.gz \
  --database blastdb/pangenome \
  --output ${SAMPLE}.blast.tsv.gz \
  --word-size 7 \
  --evalue 1e-3 \
  --max-targets 100 \
  --threads $THREADS
```

**Input formats**: FASTA, FASTQ, or gzipped (`.gz`) versions of either format. FASTQ files are automatically converted.

**Parameter Guidance:**

| Parameter | Default | Sensitive | Purpose |
|-----------|---------|-----------|---------|
| --word-size | 7 | 6 | Smaller = more sensitive, slower |
| --evalue | 1e-3 | 1e-2 | Higher = more remote homologs |
| --max-targets | 100 | 200 | More targets = better uncertainty estimation |

**For highly divergent samples:**

```bash
metadarkmatter blast align \
  --query extraction/${SAMPLE}_taxid${TAXID}_R1.fastq.gz \
  --database blastdb/pangenome \
  --output ${SAMPLE}.blast.tsv.gz \
  --word-size 6 --evalue 1e-2 --max-targets 200 --threads $THREADS
```

**Runtime Expectations:**

| Extracted Reads | Approximate Time |
|-----------------|------------------|
| 1,000 | ~10 seconds |
| 10,000 | 5-10 minutes |
| 100,000 | 30-60 minutes |
| 1,000,000 | 3-6 hours |
| 10,000,000 | 1-2 days |

### Advanced: MMseqs2 for Very Large Datasets

**⚠️ Use MMseqs2 only for datasets with >100,000 reads. For smaller datasets, use BLAST (faster and simpler).**

For very large datasets, MMseqs2 provides significant speedup at the cost of additional complexity:

```bash
# Install MMseqs2
conda install -c bioconda mmseqs2

# Create database (one-time)
metadarkmatter mmseqs2 makedb --genomes genomes/ --output mmseqs_db/pangenome

# Run search
metadarkmatter mmseqs2 search \
  --query extraction/${SAMPLE}_taxid${TAXID}_R1.fastq.gz \
  --database mmseqs_db/pangenome \
  --output ${SAMPLE}.blast.tsv.gz \
  --query-db query_cache/reads \
  --threads $THREADS
```

**When MMseqs2 is faster:**

| Reads | BLAST | MMseqs2 | Speedup |
|-------|-------|---------|---------|
| 1K-10K | 10s-10min | **Slower** | ❌ Don't use |
| 100K | 30-60 min | 5-10 min | ✓ 5x faster |
| 1M+ | 3-6 hours | 15-30 min | ✓ 10-15x faster |

**Key points:**
- Output is BLAST-compatible (works seamlessly with classification)
- Use `--query-db` to cache query database for repeated searches
- Requires 2-5x database size in temporary disk space
- Best for production pipelines with >100K reads

---

## Step 8: ANI-Weighted Classification (with AAI)

Run the core classification algorithm with species-level tracking and genus-level disambiguation:

```bash
metadarkmatter score classify \
  --blast ${SAMPLE}.blast.tsv.gz \
  --ani ani_matrix.csv \
  --aai aai_matrix.csv \
  --metadata genome_metadata.tsv \
  --output classifications.csv \
  --summary summary.json \
  --parallel
```

**Key options:**
- `--aai`: Enables genus-level disambiguation for reads with 20-25% novelty
- `--metadata`: Adds species/genus columns and enables species-level aggregation
- `--parallel`: Recommended for datasets with > 10M alignments

**AAI Impact:**
- AAI is only consulted for reads in the genus boundary range (20-25% novelty)
- If AAI > 65%: Read stays in same genus, classified as Novel Species
- If AAI < 60%: Read represents different genus, classified as Novel Genus
- Improves accuracy by 3-5% for genus-boundary reads

**Without AAI:**
If you didn't compute an AAI matrix, simply omit the `--aai` flag:

```bash
metadarkmatter score classify \
  --blast ${SAMPLE}.blast.tsv.gz \
  --ani ani_matrix.csv \
  --metadata genome_metadata.tsv \
  --output classifications.csv \
  --summary summary.json \
  --parallel
```

### Processing Modes

| Mode | Flag | Best For |
|------|------|----------|
| Standard | (none) | < 1M alignments |
| Fast | --fast | 1-10M alignments |
| Parallel | --parallel | 10-100M alignments (recommended) |
| Streaming | --streaming | > 100M alignments, memory-constrained |

Note: `--streaming` mode does not support `--metadata` for species tracking.

### Output Files

**classifications.csv** - Per-read results with species information:

```csv
read_id,best_match_genome,top_hit_identity,novelty_index,placement_uncertainty,taxonomic_call,species,genus
read_00001,GCF_000195955.2,98.7,1.3,0.5,Known Species,Francisella tularensis,Francisella
read_00002,GCF_000242755.1,91.2,8.8,1.2,Novel Species,Francisella philomiragia,Francisella
read_00003,GCF_000297895.1,82.4,17.6,0.8,Novel Genus,Francisella noatunensis,Francisella
```

**summary.json** - Sample statistics with species breakdown:

```json
{
  "total_reads": 48521,
  "known_species": 31250,
  "novel_species": 12340,
  "novel_genus": 2890,
  "conserved_regions": 1541,
  "unclassified": 500,
  "novel_diversity_pct": 31.4,
  "mean_novelty_index": 5.8,
  "species_hit_counts": {
    "Francisella tularensis": 28450,
    "Francisella philomiragia": 12340,
    "Francisella noatunensis": 5210,
    "Unknown species": 2521
  },
  "species_count": 4
}
```

### Console Output

When using `--metadata`, the console displays species information:

```
Mean Novelty Index: 5.80
Mean Placement Uncertainty: 1.24
Unique Species: 4

Top Species:
  Francisella tularensis: 28,450 (58.6%)
  Francisella philomiragia: 12,340 (25.4%)
  Francisella noatunensis: 5,210 (10.7%)
  Unknown species: 2,521 (5.2%)
```

---

## Step 9: Extract Novel Species Candidates

Identify and extract reads representing novel diversity:

```bash
# Extract all novel candidates (species + genus)
metadarkmatter score extract-novel \
  --classifications classifications.csv \
  --output novel_candidates.csv \
  --read-ids novel_read_ids.txt

# Extract only novel species
metadarkmatter score extract-novel \
  --classifications classifications.csv \
  --output novel_species.csv \
  --category species \
  --read-ids novel_species_reads.txt

# Extract only novel genus candidates
metadarkmatter score extract-novel \
  --classifications classifications.csv \
  --output novel_genus.csv \
  --category genus

# Custom threshold: extract reads with > 20% novelty
metadarkmatter score extract-novel \
  --classifications classifications.csv \
  --output highly_novel.csv \
  --min-novelty 20 \
  --read-ids highly_novel_reads.txt
```

### Output: Candidate Lineages

The output groups novel reads by their closest reference genome:

```csv
best_match_genome,read_count,mean_novelty,min_novelty,max_novelty,mean_identity,dominant_category
GCF_000242755.1,3521,8.4,5.1,14.2,91.6,Novel Species
GCF_000297895.1,1842,18.2,15.3,24.1,81.8,Novel Genus
GCF_001971545.1,1205,7.2,5.0,11.8,92.8,Novel Species
```

Interpretation:
- **GCF_000242755.1**: 3,521 reads at ~92% identity = candidate novel species related to *F. philomiragia*
- **GCF_000297895.1**: 1,842 reads at ~82% identity = candidate novel genus related to *F. noatunensis*

### Extract Reads for Assembly

Use seqtk to extract specific reads for targeted genome assembly:

```bash
# Install seqtk if not already installed
conda install -c bioconda seqtk

# Extract novel reads from original FASTQ
seqtk subseq extraction/${SAMPLE}_taxid${TAXID}_R1.fastq.gz novel_read_ids.txt > novel_R1.fastq
seqtk subseq extraction/${SAMPLE}_taxid${TAXID}_R2.fastq.gz novel_read_ids.txt > novel_R2.fastq

# Assemble with SPAdes or MEGAHIT
megahit -1 novel_R1.fastq -2 novel_R2.fastq -o novel_assembly/
```

---

## Step 10: Generate Reports

### Single-Sample Report (with ANI and AAI Matrices)

```bash
metadarkmatter report generate \
  --classifications classifications.csv \
  --metadata genome_metadata.tsv \
  --output report.html \
  --ani ani_matrix.csv \
  --aai aai_matrix.csv \
  --sample-name "Environmental Sample 01"
```

**Matrix visualizations:**
- `--ani`: Adds ANI Matrix heatmap showing nucleotide-level genome relationships
- `--aai`: Adds AAI Matrix heatmap showing protein-level genus relationships
- Both matrices are interactive with hover tooltips showing pairwise values

**Without AAI matrix:**
If you didn't compute AAI, simply omit the `--aai` flag:

```bash
metadarkmatter report generate \
  --classifications classifications.csv \
  --metadata genome_metadata.tsv \
  --ani ani_matrix.csv \
  --output report.html
```

Open `report.html` in a browser. The interactive report includes:

1. **Overview** - Classification breakdown pie chart and statistics
2. **Distributions** - Novelty index and uncertainty histograms
3. **Scatter Plot** - Reads in novelty vs. uncertainty space
4. **Recruitment** - Read recruitment plots (if BAM provided)
5. **Species** - Species-level breakdown (when `--metadata` is provided)
6. **Genomes** - Top reference genome matches
7. **ANI Matrix** - Interactive heatmap of nucleotide-level genome similarities (if `--ani` provided)
8. **AAI Matrix** - Interactive heatmap of protein-level genus relationships (if `--aai` provided)
9. **Data Table** - Searchable classification results

### Species Breakdown Tab

When you provide `--metadata`, the report includes a dedicated "Species Breakdown" tab with:

- **Species Composition Pie Chart** - Shows proportion of reads assigned to each species (top 10)
- **Species Read Count Bar Chart** - Horizontal bar chart of top 20 species by read count
- **Summary Metrics** - Unique species detected vs. reference species in database

This is particularly useful for:
- Understanding which known species dominate your sample
- Identifying species with high read counts but novel diversity (potential subspecies)
- Comparing species composition across samples

### Multi-Sample Comparison

Compare multiple samples:

```bash
# Assuming you have classifications for multiple samples
metadarkmatter report multi \
  --input-dir results/ \
  --output comparison_report.html \
  --pattern "*_classifications.csv"
```

---

## Step 11: Interpret Results

### Classification Summary

Examine your `summary.json`:

| Metric | Example Value | Interpretation |
|--------|---------------|----------------|
| novel_diversity_pct | 31.4% | ~1/3 of reads represent novel diversity |
| mean_novelty_index | 5.8 | Average divergence from references |
| novel_species | 12,340 | Reads from potential new species |
| novel_genus | 2,890 | Reads from potential new genera |
| species_count | 4 | Number of reference species detected |

### Species-Level Insights

When using `--metadata`, examine the species breakdown:

```python
import polars as pl

# Load classifications with species information
df = pl.read_csv("classifications.csv")

# Species distribution
species_stats = (
    df.group_by("species")
    .agg([
        pl.len().alias("read_count"),
        pl.col("novelty_index").mean().alias("mean_novelty"),
        pl.col("taxonomic_call").filter(pl.col("taxonomic_call") == "Novel Species").len().alias("novel_reads"),
    ])
    .sort("read_count", descending=True)
)

print(species_stats)
```

Example output:

```
species                      read_count  mean_novelty  novel_reads
Francisella tularensis       28450       2.1           1205
Francisella philomiragia     12340       8.4           9876
Francisella noatunensis      5210        12.8          4102
Unknown species              2521        18.2          2521
```

**Interpreting species patterns:**

| Pattern | What It Means |
|---------|---------------|
| High reads + low novelty | Known species, well-characterized in your sample |
| High reads + high novelty | Potential new subspecies or strain of known species |
| "Unknown species" | Reads didn't match metadata (check accession mapping) |
| Many species with similar counts | Diverse community, no dominant species |

### What Does Novel Diversity % Mean?

| Range | Interpretation |
|-------|----------------|
| < 10% | Well-characterized; reference database covers most diversity |
| 10-30% | Moderate discovery potential; typical for many environments |
| 30-50% | High discovery potential; undersampled ecosystem |
| > 50% | Very high; may indicate database gaps or misclassification |

### Validating Novel Candidates

High novelty scores require validation:

1. **Check for misclassification** - BLAST novel reads against GTDB to confirm they belong to your target family
2. **Examine read quality** - Low-quality reads can show inflated divergence
3. **Look for consistency** - Novel hits present across multiple samples are more reliable
4. **Consider contamination** - Unexpected taxa may indicate sample issues
5. **Compare to species breakdown** - Novel reads from known species may represent strain variation

---

## Advanced: Hybrid Validation Strategy

For publication-quality results, validate novel hits against the full GTDB database:

```bash
# Stage 1: Extract highly novel reads
metadarkmatter score extract-novel \
  --classifications classifications.csv \
  --output highly_novel.csv \
  --min-novelty 15 \
  --read-ids validation_reads.txt

# Stage 2: BLAST against GTDB (if you have a GTDB BLAST database)
# Note: BLAST accepts FASTQ directly - no conversion needed
metadarkmatter blast align \
  --query extraction/${SAMPLE}_taxid${TAXID}_R1.fastq.gz \
  --database /path/to/gtdb_blast_db \
  --output validation.blast.tsv.gz \
  --threads $THREADS

# Stage 3: Check results
# Reads matching non-Francisellaceae at > 95% identity are misclassifications
zcat validation.blast.tsv.gz | awk '$3 > 95' | head
```

---

## Advanced: Python API for Species Analysis

For custom analysis, use the `GenomeMetadata` class directly:

```python
from pathlib import Path
import polars as pl
from metadarkmatter.core.metadata import GenomeMetadata

# Load genome metadata
metadata = GenomeMetadata.from_file(Path("genome_metadata.tsv"))

# Print summary
print(f"Loaded metadata for {metadata.genome_count} genomes")
print(f"Covering {metadata.species_count} species in {metadata.genus_count} genera")

# Load classifications (without species columns)
df = pl.read_csv("classifications_no_metadata.csv")

# Enrich with species information
enriched_df = metadata.join_classifications(df)
print(enriched_df.columns)
# ['read_id', 'best_match_genome', ..., 'species', 'genus']

# Aggregate by species
species_summary = metadata.aggregate_by_species(enriched_df)
print(species_summary)
# Columns: species, read_count, mean_novelty, mean_identity, mean_uncertainty, genome_count

# Aggregate by genus
genus_summary = metadata.aggregate_by_genus(enriched_df)
print(genus_summary)
# Columns: genus, read_count, mean_novelty, mean_identity, mean_uncertainty, genome_count, species_count

# Look up individual genomes
species = metadata.get_species("GCF_000195955.2")
genus = metadata.get_genus("GCF_000195955.2")
print(f"GCF_000195955.2 -> {species} ({genus})")
# GCF_000195955.2 -> Francisella tularensis (Francisella)

# Access the underlying DataFrame for custom queries
metadata_df = metadata.dataframe
family_counts = (
    metadata_df.group_by("family")
    .len()
    .sort("len", descending=True)
)
print(family_counts)
```

### Combining Multiple Samples with Species Data

```python
import polars as pl
from pathlib import Path

# Load multiple classification files
samples = {}
for f in Path("results").glob("*_classifications.csv"):
    sample_name = f.stem.replace("_classifications", "")
    samples[sample_name] = pl.read_csv(f)

# Combine and compare species across samples
combined = []
for sample_name, df in samples.items():
    species_counts = (
        df.group_by("species")
        .len()
        .with_columns(pl.lit(sample_name).alias("sample"))
    )
    combined.append(species_counts)

all_samples = pl.concat(combined)

# Pivot to create species x sample matrix
species_matrix = all_samples.pivot(
    values="len",
    index="species",
    on="sample"
).fill_null(0)

print(species_matrix)
```

---

## Complete Workflow Script

Save this as `analyze_family.sh`:

```bash
#!/bin/bash
set -euo pipefail

# Usage: ./analyze_family.sh <sample> <family> <taxid> <kraken_db> [threads]
# Example: ./analyze_family.sh water_sample Francisellaceae 34064 /db/kraken2 16

SAMPLE=$1
FAMILY=$2
TAXID=$3
KRAKEN_DB=$4
THREADS=${5:-16}

echo "=== Metadarkmatter Analysis ==="
echo "Sample: $SAMPLE"
echo "Family: $FAMILY (TaxID: $TAXID)"
echo "Threads: $THREADS"
echo ""

# Create project structure
mkdir -p ${SAMPLE}_analysis/{genomes,blastdb,kraken,extraction,results}
cd ${SAMPLE}_analysis

# Step 1: Download reference genomes (auto-creates genome_metadata.tsv)
echo "Step 1: Downloading reference genomes..."
metadarkmatter download genomes list "f__${FAMILY}" --output genomes.tsv
metadarkmatter download genomes fetch --accessions genomes.tsv --output-dir genomes/ --include-protein

# Step 1b: Generate missing protein files with Prodigal (if any)
echo "Step 1b: Generating missing protein files..."
metadarkmatter proteins predict --genomes genomes/ --missing-only --threads $THREADS

# Step 2: Build ANI matrix
echo "Step 2: Computing ANI matrix..."
metadarkmatter ani compute --genomes genomes/ --output ani_matrix.csv --threads $THREADS

# Step 2b: Build AAI matrix (for genus-level disambiguation)
echo "Step 2b: Computing AAI matrix..."
metadarkmatter aai compute --genomes genomes/ --output aai_matrix.csv --threads $THREADS

# Step 3: Build BLAST database (standardizes headers for multi-contig genomes)
echo "Step 3: Building BLAST database..."
metadarkmatter blast makedb --genomes genomes/ --output blastdb/pangenome

# Step 4: Kraken2 classification
echo "Step 4: Running Kraken2 classification..."
metadarkmatter kraken2 classify \
  --reads-1 ../${SAMPLE}_R1.fastq.gz \
  --reads-2 ../${SAMPLE}_R2.fastq.gz \
  --kraken-db $KRAKEN_DB \
  --output kraken/ \
  --threads $THREADS

# Step 5: Extract family reads
echo "Step 5: Extracting family reads..."
metadarkmatter kraken2 extract \
  --kraken-output kraken/${SAMPLE}.kraken \
  --reads-1 ../${SAMPLE}_R1.fastq.gz \
  --reads-2 ../${SAMPLE}_R2.fastq.gz \
  --taxid $TAXID \
  --output extraction/

# Step 6: BLAST alignment
echo "Step 6: Running BLAST alignment..."
metadarkmatter blast align \
  --query extraction/${SAMPLE}_taxid${TAXID}_R1.fastq.gz \
  --database blastdb/pangenome \
  --output results/${SAMPLE}.blast.tsv.gz \
  --threads $THREADS

# Step 7: ANI-weighted classification with species tracking
echo "Step 7: Classifying reads..."
metadarkmatter score classify \
  --blast results/${SAMPLE}.blast.tsv.gz \
  --ani ani_matrix.csv \
  --metadata genome_metadata.tsv \
  --output results/classifications.csv \
  --summary results/summary.json \
  --parallel

# Step 8: Extract novel candidates
echo "Step 8: Extracting novel candidates..."
metadarkmatter score extract-novel \
  --classifications results/classifications.csv \
  --output results/novel_candidates.csv \
  --read-ids results/novel_read_ids.txt

# Step 9: Generate report with species breakdown
echo "Step 9: Generating report..."
metadarkmatter report generate \
  --classifications results/classifications.csv \
  --metadata genome_metadata.tsv \
  --output results/report.html \
  --ani ani_matrix.csv \
  --sample-name "$SAMPLE"

echo ""
echo "=== Analysis Complete ==="
echo "Results: ${SAMPLE}_analysis/results/"
echo "Report:  ${SAMPLE}_analysis/results/report.html"
echo ""

# Display summary with species information
echo "Summary:"
cat results/summary.json | python -m json.tool

# Display species breakdown
echo ""
echo "Files created:"
echo "  - genome_metadata.tsv (species/genus taxonomy)"
echo "  - blastdb/pangenome_contig_mapping.tsv (contig-to-genome mapping)"
echo "  - results/classifications.csv (with species columns)"
echo "  - results/report.html (with Species Breakdown tab)"
```

Run the script:

```bash
chmod +x analyze_family.sh
./analyze_family.sh water_sample Francisellaceae 34064 /path/to/kraken2_db 16
```

---

## Troubleshooting

### Very Few Reads Extracted

**Symptoms:** < 100 reads extracted for your target family

**Solutions:**
1. Check family abundance in Kraken report: `grep "YourFamily" kraken_output/*.kreport`
2. Verify TaxID is correct (check NCBI Taxonomy)
3. Consider using a more comprehensive Kraken2 database

### High "Conserved Region" Percentage

**Symptoms:** > 50% of reads classified as Conserved Region

**Solutions:**
1. Reduce `--max-targets` to 50 in BLAST alignment
2. Check if reference genomes are too similar (redundant strains)
3. This is expected for highly conserved genes (16S, housekeeping genes)

### Unexpectedly High Novel Diversity

**Symptoms:** > 60% novel species/genus

**Solutions:**
1. Validate with hybrid database strategy (BLAST against GTDB)
2. Check for Kraken2 misclassification
3. May be genuine - some environments have high novel diversity

### BLAST Alignment Too Slow

**Symptoms:** Alignment taking > 6 hours for very large datasets (>100K reads)

**Solutions:**

1. **Use more threads**: `--threads 32` (or more if available)
2. **Reduce sensitivity** (if appropriate): `--word-size 11 --evalue 1e-5`
3. **Run on HPC cluster** with job parallelization
4. **For >100K reads only**: Consider MMseqs2 (see Step 7 Advanced section)
   - ⚠️ Only use MMseqs2 for datasets >100,000 reads
   - For smaller datasets, MMseqs2 is actually slower than BLAST
   - See "Decision Guide: BLAST vs MMseqs2" in Step 7 for details

---

## Summary

This tutorial demonstrated how to use metadarkmatter to discover novel bacterial diversity within a target family. The workflow:

1. **Downloads reference genomes** from GTDB (with automatic metadata extraction)
2. **Builds an ANI matrix** for uncertainty calculation
3. **Builds sequence databases** with standardized headers for multi-contig genome tracking
4. **Extracts family reads** using Kraken2
5. **Aligns reads** competitively against references using BLAST
6. **Classifies reads** by novelty, placement uncertainty, and species assignment
7. **Extracts novel candidates** for further analysis
8. **Generates reports** with species breakdown visualizations

### Key Features

- **Species-Level Tracking**: Genome metadata propagates through the pipeline, enabling species-level aggregation in results and reports
- **Multi-Contig Support**: Standardized FASTA headers (`{accession}|{contig_id}`) ensure reliable genome identification from draft assemblies
- **Scalable Performance**: BLAST handles most datasets efficiently; MMseqs2 available for very large datasets (>100K reads)
- **Integrated Reporting**: HTML reports include a dedicated "Species Breakdown" tab showing composition charts and read counts per species

### Output Files Summary

| File | Purpose |
|------|---------|
| `genome_metadata.tsv` | Species/genus taxonomy for all reference genomes |
| `contig_mapping.tsv` | Links contigs to parent genomes for draft assemblies |
| `classifications.csv` | Per-read results with species and genus columns |
| `summary.json` | Sample statistics including species_hit_counts |
| `report.html` | Interactive report with Species Breakdown tab |

The ANI-weighted approach provides a quantitative framework for distinguishing known species from potential novel species and genera, enabling targeted discovery of microbial dark matter.

---

## References

- Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S. (2018). High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. *Nature Communications*, 9:5114.
- Parks DH, et al. (2022). GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy. *Nucleic Acids Research*, 50:D1.
- Wood DE, Lu J, Langmead B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20:257.
- Shaw J, Yu YW. (2023). Fast and robust metagenomic sequence comparison through sparse chaining with skani. *Nature Methods*, 20:1661-1665.
- Steinegger M, Söding J. (2017). MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. *Nature Biotechnology*, 35:1026-1028.
