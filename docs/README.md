# Metadarkmatter Documentation

Welcome to the metadarkmatter documentation.

## External Software Dependencies

Metadarkmatter requires external bioinformatics tools (non-Python software):

**Required:**
```bash
conda install -c bioconda kraken2 krakentools blast skani
```

| Tool | Purpose |
|------|---------|
| Kraken2 | Taxonomic classification |
| KrakenTools | Read extraction by taxid |
| BLAST+ | Nucleotide alignment |
| skani | Fast ANI computation |

**Optional:**
- **MMseqs2**: For datasets >100K reads only (`conda install -c bioconda mmseqs2`)
- **seqtk**: For assembly workflows (`conda install -c bioconda seqtk`)

---

## Getting Started

**New users start here:**

1. **[../README.md](../README.md)** - Installation and quick start
2. **[TUTORIAL_ENVIRONMENTAL_SPECIES.md](TUTORIAL_ENVIRONMENTAL_SPECIES.md)** - Complete hands-on tutorial
3. **[USER_GUIDE.md](USER_GUIDE.md)** - Practical usage guide with output interpretation

## Core Documentation

**Workflows and Usage:**
- **[WORKFLOW.md](WORKFLOW.md)** - Workflow decisions and database strategies
- **[USER_GUIDE.md](USER_GUIDE.md)** - Input formats, output interpretation, common workflows
- **[CLI_REFERENCE.md](CLI_REFERENCE.md)** - Complete command-line reference
- **[REFERENCE.md](REFERENCE.md)** - Algorithm overview and technical reference
- **[ALGORITHM_DETAILED.md](ALGORITHM_DETAILED.md)** - Complete algorithm explanation with ANI/AAI decision tree
- **[ALIGNMENT_OUTPUT_STATISTICS.md](ALIGNMENT_OUTPUT_STATISTICS.md)** - BLAST/MMseqs2 output statistics explained

**Optimization and Troubleshooting:**
- **[PERFORMANCE.md](PERFORMANCE.md)** - Performance optimization guide
- **[TROUBLESHOOTING.md](TROUBLESHOOTING.md)** - Common issues and solutions

## Advanced Topics

- **[mmseqs2-integration.md](mmseqs2-integration.md)** - MMseqs2 for very large datasets (>100K reads)
- **[CLASSIFICATION_STATISTICS.md](CLASSIFICATION_STATISTICS.md)** - Statistical framework and literature

## Technical Documentation

**For developers and contributors:**
- **[ARCHITECTURE.md](ARCHITECTURE.md)** - System architecture and design
- **[TECHNICAL_MANUAL.md](TECHNICAL_MANUAL.md)** - Algorithm implementation details
- **[API_REFERENCE.md](API_REFERENCE.md)** - Python API documentation
- **[CHANGELOG.md](CHANGELOG.md)** - Project changes and updates

---

## Quick Links by Task

### Find novel diversity in environmental samples
**→ [TUTORIAL_ENVIRONMENTAL_SPECIES.md](TUTORIAL_ENVIRONMENTAL_SPECIES.md)**

### Classify metagenomic reads
```bash
metadarkmatter score classify \
    --alignment sample.blast.tsv.gz \
    --ani ani_matrix.csv \
    --metadata genome_metadata.tsv \
    --output classifications.csv
```
**See:** [CLI_REFERENCE.md](CLI_REFERENCE.md) or [USER_GUIDE.md](USER_GUIDE.md)

### Interpret classification output
**→ [USER_GUIDE.md](USER_GUIDE.md#output-interpretation)**

### Run complete workflow
**→ [TUTORIAL_ENVIRONMENTAL_SPECIES.md](TUTORIAL_ENVIRONMENTAL_SPECIES.md)** or **[WORKFLOW.md](WORKFLOW.md)**

### Optimize for large datasets
**→ [PERFORMANCE.md](PERFORMANCE.md)** and **[mmseqs2-integration.md](mmseqs2-integration.md)**

### Fix common errors
**→ [TROUBLESHOOTING.md](TROUBLESHOOTING.md)**

### Use Python API
**→ [API_REFERENCE.md](API_REFERENCE.md)**

---

## Key Concepts

**Novelty Index (N)**: `100 - percent_identity`
- N < 5%: Known species
- 5-20%: Novel species
- 20-25%: Novel genus

**Placement Uncertainty (U)**: `100 - max(ANI between competing genomes)`
- U < 2%: Clear placement
- U > 2%: Conserved region (ambiguous)

**Species-Level Tracking**: Use `--metadata genome_metadata.tsv` to track diversity by species/genus in classifications and reports.

---

## Documentation Structure

```
docs/
├── README.md                          # This file - documentation index
│
├── Getting Started
│   ├── TUTORIAL_ENVIRONMENTAL_SPECIES.md    # Hands-on tutorial (start here)
│   ├── USER_GUIDE.md                        # Practical usage guide
│   └── WORKFLOW.md                          # Workflow strategies
│
├── Reference
│   ├── CLI_REFERENCE.md                     # Complete command reference
│   ├── REFERENCE.md                         # Algorithm overview
│   ├── ALGORITHM_DETAILED.md                # Detailed algorithm with ANI/AAI decision tree
│   ├── ALIGNMENT_OUTPUT_STATISTICS.md       # BLAST/MMseqs2 statistics explained
│   ├── PERFORMANCE.md                       # Performance optimization
│   └── TROUBLESHOOTING.md                   # Common issues
│
├── Advanced
│   ├── mmseqs2-integration.md               # MMseqs2 for large datasets
│   └── CLASSIFICATION_STATISTICS.md         # Statistical framework
│
└── Technical
    ├── ARCHITECTURE.md                      # System design
    ├── TECHNICAL_MANUAL.md                  # Implementation details
    ├── API_REFERENCE.md                     # Python API
    └── CHANGELOG.md                         # Project changes
```

---

## Important Notes

**BLAST Input Formats:**
- BLAST accepts FASTQ, FASTA, or gzipped versions (.gz)
- FASTQ files are automatically converted to FASTA (transparent)
- No manual conversion with seqtk required

**MMseqs2 Usage:**
- Only use MMseqs2 for datasets >100,000 reads
- For smaller datasets, BLAST is faster and simpler
- See [mmseqs2-integration.md](mmseqs2-integration.md) for details

**Classification Modes:**
- Default (nucleotide): For BLASTN results
- Protein mode (`--alignment-mode protein`): For BLASTX results with different thresholds

---

## Support

GitHub Issues: https://github.com/metadarkmatter/metadarkmatter/issues
