# Metadarkmatter Workflow Guide

A step-by-step guide to analyze metagenomic samples for novel bacterial diversity.

## Database Strategy: Family vs GTDB

Before starting, you must choose a BLAST database strategy. This decision affects sensitivity, specificity, and computational requirements.

### Option A: Family BLAST Database (Recommended for most analyses)

**Use genomes only from your target bacterial family.**

```
family_blast_db/
└── Contains: 100-500 genomes from target family only
```

**Advantages:**
- Fast alignment (30 min - 2 hr for 10M reads)
- All hits are within family context
- Novelty Index directly measures divergence from known family members
- Lower computational resources (~50 GB storage, 32 GB RAM)

**Disadvantages:**
- Reads misclassified by Kraken2 will be forced to match family genomes
- These will appear as "highly novel" - potential false positives
- Cannot distinguish true novel genera from Kraken2 errors

**Best for:**
- Well-characterized families with good Kraken2 database coverage
- Samples where target family dominance is expected
- Rapid screening of multiple samples

### Option B: GTDB Representative Database

**Use all ~85,000 representative bacterial genomes from GTDB.**

```
gtdb_blast_db/
└── Contains: All GTDB r220 representative genomes
```

**Advantages:**
- Comprehensive context - reads can match ANY bacterial family
- Validates that reads truly belong to target family
- Reduces false positives from Kraken2 misclassification
- Better for discovering deep-branching novel lineages

**Disadvantages:**
- Large database (~500 GB storage)
- Slow alignment (10-50x longer than family DB)
- Requires post-filtering to focus on family of interest
- Higher memory requirements (64-128 GB RAM)

**Best for:**
- Discovery-focused studies where false positives are costly
- Samples with uncertain taxonomic composition
- Validation of surprising results from family DB analysis

### Option C: Hybrid Two-Stage Approach (Recommended for publication)

**Combine both strategies for maximum confidence.**

```
Stage 1: Family DB (fast screening)
    ↓
Stage 2: GTDB DB (validation of novel hits only)
```

**Workflow:**

```bash
# Stage 1: Primary analysis with family database
metadarkmatter blast align \
  --query extracted_reads.fasta \
  --database family_blast_db/pangenome \
  --output stage1.blast.tsv.gz

metadarkmatter score classify \
  --blast stage1.blast.tsv.gz \
  --ani family_ani_matrix.csv \
  --output stage1_classifications.csv

# Stage 2: Extract reads classified as highly novel
python << 'EOF'
import polars as pl
df = pl.read_csv("stage1_classifications.csv")
novel = df.filter(
    (pl.col("taxonomic_call") == "Novel Genus") |
    (pl.col("novelty_index") > 20)
)
novel.select("read_id").write_csv("novel_read_ids.txt", include_header=False)
EOF

# Extract novel reads for validation
seqtk subseq extracted_reads.fasta novel_read_ids.txt > novel_reads.fasta

# Stage 2: Validate against GTDB
metadarkmatter blast align \
  --query novel_reads.fasta \
  --database gtdb_blast_db/gtdb_reps \
  --output stage2_validation.blast.tsv.gz

# Check if "novel" reads match known species outside family
python << 'EOF'
import polars as pl
validation = pl.read_csv("stage2_validation.blast.tsv.gz", separator="\t",
                         has_header=False,
                         new_columns=["qseqid", "sseqid", "pident", "length",
                                     "mismatch", "gapopen", "qstart", "qend",
                                     "sstart", "send", "evalue", "bitscore"])

# Get best hit per read
best_hits = validation.group_by("qseqid").agg(pl.all().sort_by("bitscore").last())

# Check: if best hit is >95% identity to non-family genome, it's likely misclassification
high_identity = best_hits.filter(pl.col("pident") > 95)
print(f"Reads with >95% identity to GTDB genomes: {len(high_identity)}")
print("These may be Kraken2 misclassifications, not true novel diversity")
EOF
```

**Interpretation:**
- Reads that match GTDB genomes at >95% identity are likely Kraken2 errors
- Reads with no good GTDB matches (<85% identity) are more likely truly novel
- Reads matching family genomes best in GTDB confirm family assignment

### Decision Matrix

| Scenario | Recommended Strategy |
|----------|---------------------|
| Rapid screening, many samples | Family DB only |
| Single sample, publication-quality | Hybrid two-stage |
| Unknown sample composition | GTDB DB with post-filtering |
| Well-characterized ecosystem | Family DB only |
| Discovery of deep-branching lineages | GTDB DB or hybrid |
| Limited computational resources | Family DB only |

### Building the Databases

**Family database:**
```bash
metadarkmatter download genomes list "f__YourFamily" --output family_genomes.tsv
metadarkmatter download genomes fetch --accessions family_genomes.tsv --output-dir family_genomes/
metadarkmatter blast makedb --genomes family_genomes/ --output family_blast_db/pangenome
```

**GTDB database (requires significant resources):**
```bash
# Download GTDB representative genomes (~50 GB compressed)
wget https://data.gtdb.ecogenomic.org/releases/release220/220.0/genomic_files_reps/gtdb_genomes_reps_r220.tar.gz
tar -xzf gtdb_genomes_reps_r220.tar.gz

# Build database (requires ~500 GB disk, several hours)
metadarkmatter blast makedb --genomes gtdb_genomes_reps/ --output gtdb_blast_db/gtdb_reps
```

### ANI Matrix Considerations

For family DB analysis, build ANI matrix only for family genomes:
```bash
fastANI --ql family_genomes.txt --rl family_genomes.txt -o family_ani.txt
```

For GTDB validation, ANI matrix is typically not needed - you're checking for presence/absence of good matches, not nuanced placement uncertainty.

---

### Option D: Protein-Level Read Classification (BLASTX)

**Use translated nucleotide alignment for divergent sequences.**

For highly divergent taxa where nucleotide alignment fails to detect homology, protein-level classification using Diamond BLASTX can capture divergent homology that BLASTN misses.

**Why BLASTX?**

| Scenario | BLASTN | BLASTX |
|----------|--------|--------|
| Known species (>95% ANI) | Excellent | Good |
| Novel species (80-95% ANI) | Good | Good |
| Novel genus (60-80% ANI) | Poor | **Better** |
| Deep divergence (<60% ANI) | Fails | **Detects** |

**Protein workflow:**

```bash
# Build protein database from reference proteomes
metadarkmatter blastx makedb \
  --proteins reference_proteins/ \
  --output blastdb/panproteome \
  --threads 16

# Run BLASTX alignment (DNA reads vs protein database)
metadarkmatter blastx align \
  --query extracted_reads.fastq.gz \
  --database blastdb/panproteome \
  --output sample.blastx.tsv.gz \
  --threads 16

# Classify with protein-calibrated thresholds
metadarkmatter score classify \
  --blast sample.blastx.tsv.gz \
  --ani ani_matrix.csv \
  --alignment-mode protein \
  --output classifications.csv \
  --parallel
```

**Key differences from nucleotide mode:**

| Threshold | Nucleotide | Protein |
|-----------|-----------|---------|
| Known Species | N < 5% | N < 10% |
| Novel Species | 5-20% | 10-25% |
| Novel Genus | 20-25% | 25-40% |

**Best for:**
- Genus-level novelty detection
- Samples with highly divergent sequences
- When BLASTN fails to find homology

---

## Prerequisites

Before starting, ensure you have:

- [ ] Python 3.11 or higher
- [ ] metadarkmatter installed (`pip install -e .`)
- [ ] Kraken2 database for your target ecosystem
- [ ] Reference genomes for your bacterial family of interest
- [ ] Precomputed ANI matrix for reference genomes (see Step 2)

External tools (install via conda):
```bash
conda install -c bioconda kraken2 krakentools blast bowtie2 samtools fastani diamond
```

## Workflow Overview

```
Metagenome reads
       |
       v
[1. Download genomes] --> Reference genomes
       |                        |
       v                        v
[2. Build ANI matrix]    [3. Extract reads] --> Family reads
       |                        |
       v                        v
   ANI matrix            [4. Build BLAST DB]
       |                        |
       |                        v
       |                 [5. BLAST alignment]
       |                        |
       +------------------------+
                   |
                   v
          [6. Classification]
                   |
                   v
          [7. Report generation]
```

---

## Step 1: Download Reference Genomes

Query GTDB for genomes in your target bacterial family:

```bash
# List available genomes (example: Francisellaceae)
metadarkmatter download genomes list "f__Francisellaceae" \
  --output genomes.tsv \
  --include-metadata

# Download genomes from NCBI
metadarkmatter download genomes fetch \
  --accessions genomes.tsv \
  --output-dir reference_genomes/
```

**Expected output:**
- `genomes.tsv` - List of genome accessions
- `reference_genomes/` - Directory with genome FASTA files

---

## Step 2: Build ANI Matrix

Compute pairwise ANI values between all reference genomes:

```bash
# Create genome list
ls reference_genomes/*.fna > genome_list.txt

# Run fastANI (all vs all)
fastANI \
  --ql genome_list.txt \
  --rl genome_list.txt \
  -o ani_output.txt \
  -t 16
```

Convert to matrix format:

```python
import pandas as pd

# Read fastANI output
ani = pd.read_csv('ani_output.txt', sep='\t',
                  names=['genome1', 'genome2', 'ANI', 'frag1', 'frag2'])

# Extract genome IDs from paths
ani['genome1'] = ani['genome1'].str.extract(r'(GCF_\d+\.\d+)')
ani['genome2'] = ani['genome2'].str.extract(r'(GCF_\d+\.\d+)')

# Pivot to matrix
genomes = sorted(set(ani['genome1']) | set(ani['genome2']))
matrix = pd.DataFrame(100.0, index=genomes, columns=genomes)

for _, row in ani.iterrows():
    matrix.loc[row['genome1'], row['genome2']] = row['ANI']
    matrix.loc[row['genome2'], row['genome1']] = row['ANI']

matrix.to_csv('ani_matrix.csv')
```

**Expected output:**
- `ani_matrix.csv` - Symmetric ANI matrix

---

## Step 3: Extract Family Reads

This step has two parts: taxonomic classification with Kraken2, then extraction of reads for your target family.

### 3a: Run Kraken2 Classification

```bash
metadarkmatter kraken2 classify \
  --reads-1 sample_R1.fastq.gz \
  --reads-2 sample_R2.fastq.gz \
  --kraken-db /path/to/kraken2_db \
  --output kraken_output/ \
  --threads 16
```

**Parameters:**
- `--confidence`: Kraken2 confidence threshold (0-1, default: 0)
- `--minimum-hit-groups`: Minimum hit groups (default: 2)

**Expected output:**
- `kraken_output/sample.kraken` - Per-read classification
- `kraken_output/sample.kreport` - Hierarchical taxonomic report

### 3b: Extract Target Family Reads

```bash
metadarkmatter kraken2 extract \
  --kraken-output kraken_output/sample.kraken \
  --reads-1 sample_R1.fastq.gz \
  --reads-2 sample_R2.fastq.gz \
  --taxid 262 \
  --output extraction/
```

**Parameters:**
- `--taxid`: NCBI TaxID for your family (find at taxonomy.ncbi.nlm.nih.gov)
- `--include-children/--exact-match`: Include all child taxa (default) or exact match only

**Expected output:**
- `extraction/sample_taxid262_R1.fastq` - Extracted forward reads
- `extraction/sample_taxid262_R2.fastq` - Extracted reverse reads

---

## Step 4: Build BLAST Database

Create BLAST database from reference genomes:

```bash
metadarkmatter blast makedb \
  --genomes reference_genomes/ \
  --output blastdb/pangenome
```

**Expected output:**
- `blastdb/pangenome.nhr`
- `blastdb/pangenome.nin`
- `blastdb/pangenome.nsq`

---

## Step 5: Run BLAST Alignment

Competitive alignment with high-sensitivity parameters:

```bash
metadarkmatter blast align \
  --query extraction/sample_family_R1.fastq \
  --database blastdb/pangenome \
  --output sample.blast.tsv.gz \
  --word-size 7 \
  --evalue 1e-3 \
  --max-targets 100 \
  --threads 16
```

**Parameters for divergent sequence detection:**
- `--word-size 7`: Smaller seeds for sensitive detection (default: 7)
- `--evalue 1e-3`: Relaxed threshold for remote homology
- `--max-targets 100`: Capture multiple hits per read

**Expected output:**
- `sample.blast.tsv.gz` - Compressed BLAST results

**Runtime estimate:** 3-6 hours for 10M reads

---

## Step 6: ANI-Weighted Classification

Run the core classification algorithm:

```bash
# For most datasets (10-100M alignments)
metadarkmatter score classify \
  --blast sample.blast.tsv.gz \
  --ani ani_matrix.csv \
  --output classifications.csv \
  --summary summary.json \
  --parallel
```

**Processing modes:**

| Dataset Size | Flag | RAM Usage |
|--------------|------|-----------|
| < 1M alignments | (default) | 2-4 GB |
| 1-10M | `--fast` | 4-8 GB |
| 10-100M | `--parallel` | 8-16 GB |
| 100M+ | `--streaming` | 8-16 GB |

**Expected output:**
- `classifications.csv` - Per-read classifications
- `summary.json` - Sample statistics

---

## Step 7: Generate Reports

Create HTML report for visualization:

```bash
# Single sample report
metadarkmatter report generate \
  --classifications classifications.csv \
  --output report.html \
  --ani ani_matrix.csv \
  --sample-name "Sample01"
```

For multiple samples:

```bash
# Process all samples first, then compare
metadarkmatter report multi \
  --input-dir results/ \
  --output comparison.html \
  --pattern "*_classifications.csv"
```

**Expected output:**
- `report.html` - Interactive HTML report with:
  - Classification distribution
  - Novelty vs uncertainty scatter plot
  - Top genome matches
  - ANI heatmap

---

## Interpreting Results

### Classification Categories

| Category | Novelty Index | Placement Uncertainty | Interpretation |
|----------|---------------|----------------------|----------------|
| Known Species | < 2% | < 0.5% | Matches reference genome |
| Novel Species | 5-15% | < 0.5% | New species in known genus |
| Novel Genus | 15-25% | < 2% | New genus in family |
| Conserved Region | any | > 5% | Ambiguous (rRNA, etc.) |

### Key Metrics

- **Novel diversity %**: `(novel_species + novel_genus) / total_reads * 100`
  - < 10%: Well-characterized community
  - 10-30%: Typical environmental sample
  - \> 30%: High discovery potential

### Report Sections

1. **Overview**: Summary statistics and classification pie chart
2. **Distributions**: Novelty and uncertainty histograms
3. **Scatter Plot**: Novelty vs uncertainty with classification regions
4. **Genomes**: Top matched reference genomes
5. **ANI Matrix**: Genome-genome similarity heatmap

---

## Batch Processing

For multiple samples with the same reference database:

```bash
# Process all samples
for sample in sample_A sample_B sample_C; do
  # Run Kraken2 classification
  metadarkmatter kraken2 classify \
    --reads-1 data/${sample}_R1.fastq.gz \
    --kraken-db /path/to/kraken_db \
    --output kraken_output/${sample}/

  # Extract target family reads
  metadarkmatter kraken2 extract \
    --kraken-output kraken_output/${sample}/${sample}.kraken \
    --reads-1 data/${sample}_R1.fastq.gz \
    --taxid 262 \
    --output extraction/${sample}/

  # BLAST alignment
  metadarkmatter blast align \
    --query extraction/${sample}/${sample}_taxid262_R1.fastq \
    --database blastdb/pangenome \
    --output results/${sample}.blast.tsv.gz
done

# Batch classification
metadarkmatter score batch \
  --blast-dir results/ \
  --ani ani_matrix.csv \
  --output-dir classifications/ \
  --parallel

# Multi-sample comparison
metadarkmatter report multi \
  --input-dir classifications/ \
  --output comparison_report.html
```

---

## Troubleshooting

### Common Issues

**"BLAST database not found"**
```bash
# Verify database files exist
ls blastdb/pangenome.{nhr,nin,nsq}
```

**"Genome not found in ANI matrix"**
```bash
# Check genome IDs match between BLAST results and ANI matrix
zcat sample.blast.tsv.gz | cut -f2 | sort -u | head
head -1 ani_matrix.csv
```

**Out of memory during classification**
```bash
# Use streaming mode for large files
metadarkmatter score classify --streaming ...
```

**High percentage of "Conserved Region"**
- Reduce `--max-targets` in BLAST step
- May indicate rRNA or mobile element reads

See `docs/TROUBLESHOOTING.md` for more solutions.
