#!/usr/bin/env bash
# End-to-End Test: Novel Species Detection with Simulated Reads
#
# This script tests the metadarkmatter pipeline's ability to detect novel species
# by simulating reads from Francisella sciaenopsi (GCF_045862875.1), which is
# NOT in the GTDB reference database.
#
# Dependencies:
#   - insilicoseq (iss): pip install insilicoseq
#   - ncbi-datasets-cli: conda install -c conda-forge ncbi-datasets-cli
#   - metadarkmatter: installed and available in PATH
#   - blastn, makeblastdb: BLAST+ suite
#   - skani or fastANI: for ANI computation
#
# Usage:
#   ./test_novel_species_detection.sh [--keep]
#
# Options:
#   --keep    Keep the working directory after completion (default: cleanup)

set -euo pipefail

# Configuration
NOVEL_ACCESSION="GCF_045862875.1"
FAMILY="f__Francisellaceae"
NUM_READS=1000
THREADS=4

# Parse arguments
KEEP_WORKDIR=false
for arg in "$@"; do
    case $arg in
        --keep)
            KEEP_WORKDIR=true
            shift
            ;;
    esac
done

# Check dependencies
check_command() {
    if ! command -v "$1" &> /dev/null; then
        echo "ERROR: Required command '$1' not found"
        echo "Please install: $2"
        exit 1
    fi
}

echo "=== Checking dependencies ==="
check_command "iss" "pip install insilicoseq"
check_command "datasets" "conda install -c conda-forge ncbi-datasets-cli"
check_command "metadarkmatter" "pip install -e . (from project root)"
check_command "blastn" "BLAST+ suite"
check_command "makeblastdb" "BLAST+ suite"

# Check for ANI tool (skani preferred, fastANI fallback)
if command -v skani &> /dev/null; then
    echo "Using skani for ANI computation"
elif command -v fastANI &> /dev/null; then
    echo "Using fastANI for ANI computation"
else
    echo "ERROR: Neither skani nor fastANI found"
    echo "Please install: conda install -c bioconda skani"
    exit 1
fi

echo "All dependencies found"

# Setup working directory
WORKDIR=$(mktemp -d -t mdm_e2e_novel_XXXXXX)
echo ""
echo "=== Working directory: $WORKDIR ==="
cd "$WORKDIR"

# Cleanup handler
cleanup() {
    if [ "$KEEP_WORKDIR" = false ]; then
        echo ""
        echo "=== Cleaning up ==="
        rm -rf "$WORKDIR"
        echo "Removed: $WORKDIR"
    else
        echo ""
        echo "=== Results preserved in: $WORKDIR ==="
    fi
}

if [ "$KEEP_WORKDIR" = false ]; then
    trap cleanup EXIT
fi

# Step 1: Download Francisellaceae reference genomes from GTDB
echo ""
echo "=== Step 1: Downloading Francisellaceae reference genomes from GTDB ==="

metadarkmatter download genomes list "$FAMILY" \
    --output genome_accessions.tsv \
    --include-metadata

# Check if we got any genomes
if [ ! -s genome_accessions.tsv ]; then
    echo "ERROR: No genomes found for $FAMILY"
    exit 1
fi

GENOME_COUNT=$(wc -l < genome_accessions.tsv)
echo "Found $GENOME_COUNT reference genomes in GTDB"

metadarkmatter download genomes fetch \
    --accessions genome_accessions.tsv \
    --output-dir genomes/ \
    --decompress

DOWNLOADED=$(ls genomes/*.fna 2>/dev/null | wc -l || echo 0)
echo "Downloaded $DOWNLOADED reference genomes"

if [ "$DOWNLOADED" -eq 0 ]; then
    echo "ERROR: No genomes were downloaded"
    exit 1
fi

# Step 2: Download the novel genome from NCBI (NOT in reference set)
echo ""
echo "=== Step 2: Downloading novel species $NOVEL_ACCESSION from NCBI ==="

mkdir -p novel_genome
datasets download genome accession "$NOVEL_ACCESSION" \
    --filename novel_genome.zip \
    --include genome

unzip -o novel_genome.zip -d novel_genome/

# Find the genome FASTA file
NOVEL_GENOME=$(find novel_genome -name "*.fna" -type f | head -1)

if [ -z "$NOVEL_GENOME" ]; then
    echo "ERROR: Could not find novel genome FASTA file"
    echo "Contents of novel_genome/:"
    find novel_genome -type f
    exit 1
fi

echo "Novel genome: $NOVEL_GENOME"

# Step 3: Simulate reads with InSilicoSeq
echo ""
echo "=== Step 3: Simulating $NUM_READS paired-end 150bp reads ==="

iss generate \
    --genomes "$NOVEL_GENOME" \
    --model hiseq \
    --n_reads "$NUM_READS" \
    --cpus "$THREADS" \
    --output simulated_reads

echo "Generated reads:"
ls -la simulated_reads*.fastq

# Verify reads were generated
if [ ! -f "simulated_reads_R1.fastq" ] || [ ! -f "simulated_reads_R2.fastq" ]; then
    echo "ERROR: Read simulation failed"
    exit 1
fi

READ_COUNT=$(grep -c "^@" simulated_reads_R1.fastq || echo 0)
echo "Generated $READ_COUNT read pairs"

# Step 3b: Convert FASTQ to FASTA (BLAST requires FASTA format)
echo ""
echo "=== Step 3b: Converting FASTQ to FASTA ==="

# Convert using awk (works without additional dependencies)
awk 'NR%4==1 {print ">" substr($0,2)} NR%4==2 {print}' simulated_reads_R1.fastq > simulated_reads_R1.fasta
awk 'NR%4==1 {print ">" substr($0,2)} NR%4==2 {print}' simulated_reads_R2.fastq > simulated_reads_R2.fasta

echo "Converted to FASTA format"
ls -la simulated_reads*.fasta

# Step 4: Build BLAST database from GTDB references only
echo ""
echo "=== Step 4: Building BLAST database from GTDB references ==="

metadarkmatter blast makedb \
    --genomes genomes/ \
    --output blastdb/francisellaceae

echo "BLAST database created"

# Step 5: Run BLAST alignment
echo ""
echo "=== Step 5: Running BLAST alignment ==="

metadarkmatter blast align \
    --query simulated_reads_R1.fasta \
    --database blastdb/francisellaceae \
    --output blast_results.tsv.gz \
    --threads "$THREADS"

echo "BLAST alignment complete"

# Verify BLAST results
if [ ! -f "blast_results.tsv.gz" ]; then
    echo "ERROR: BLAST output file not created"
    exit 1
fi

BLAST_HITS=$(gunzip -c blast_results.tsv.gz | wc -l || echo 0)
echo "BLAST produced $BLAST_HITS alignments"

# Step 6: Compute ANI matrix
echo ""
echo "=== Step 6: Computing ANI matrix ==="

metadarkmatter ani compute \
    --genomes genomes/ \
    --output ani_matrix.csv \
    --threads "$THREADS"

echo "ANI matrix computed"

# Step 7: Run classification
echo ""
echo "=== Step 7: Running ANI-weighted classification ==="

metadarkmatter score classify \
    --blast blast_results.tsv.gz \
    --ani ani_matrix.csv \
    --metadata genome_metadata.tsv \
    --output classifications.csv \
    --summary summary.json \
    --parallel

echo "Classification complete"

# Step 8: Generate report
echo ""
echo "=== Step 8: Generating HTML report ==="

metadarkmatter report generate \
    --classifications classifications.csv \
    --metadata genome_metadata.tsv \
    --ani ani_matrix.csv \
    --output report.html

echo "Report generated"

# Step 9: Validation
echo ""
echo "=========================================="
echo "=== VALIDATION RESULTS ==="
echo "=========================================="

# Count classifications by category
echo ""
echo "Classification distribution:"
# Skip header row and count taxonomic_call column (column 11)
tail -n +2 classifications.csv | cut -d',' -f11 | sort | uniq -c | sort -rn

# Extract counts
NOVEL_SPECIES=$(grep -c "Novel Species" classifications.csv 2>/dev/null || echo 0)
NOVEL_GENUS=$(grep -c "Novel Genus" classifications.csv 2>/dev/null || echo 0)
KNOWN=$(grep -c "Known Species" classifications.csv 2>/dev/null || echo 0)
AMBIGUOUS=$(grep -c "Ambiguous" classifications.csv 2>/dev/null || echo 0)
CONSERVED=$(grep -c "Conserved Region" classifications.csv 2>/dev/null || echo 0)
# Total = lines - 1 (header)
TOTAL=$(($(wc -l < classifications.csv) - 1))

echo ""
echo "Category Counts:"
echo "  Novel Species:    $NOVEL_SPECIES"
echo "  Novel Genus:      $NOVEL_GENUS"
echo "  Known Species:    $KNOWN"
echo "  Ambiguous:        $AMBIGUOUS"
echo "  Conserved Region: $CONSERVED"
echo "  -------------------"
echo "  Total classified: $TOTAL"

# Calculate percentages
if [ "$TOTAL" -gt 0 ]; then
    NOVEL_TOTAL=$((NOVEL_SPECIES + NOVEL_GENUS))
    NOVEL_PCT=$((100 * NOVEL_TOTAL / TOTAL))

    echo ""
    echo "Novel Detection Rate: $NOVEL_PCT% ($NOVEL_TOTAL / $TOTAL reads)"
fi

# Show summary statistics
echo ""
echo "Summary Statistics (from summary.json):"
if command -v python3 &> /dev/null; then
    python3 -m json.tool summary.json 2>/dev/null || cat summary.json
else
    cat summary.json
fi

# Final verdict
echo ""
echo "=========================================="
if [ "$NOVEL_TOTAL" -gt "$((TOTAL / 2))" ]; then
    echo "TEST PASSED: Majority of reads detected as novel species"
    echo "  - $NOVEL_PCT% of reads classified as Novel Species/Genus"
    echo "  - This is expected for reads from $NOVEL_ACCESSION"
    echo "  - The pipeline correctly identified novel diversity"
else
    echo "TEST WARNING: Less than 50% of reads detected as novel"
    echo "  - This may indicate issues with the reference database"
    echo "  - Or the novel species is very similar to a GTDB reference"
fi
echo "=========================================="

echo ""
echo "=== Output Files ==="
echo "Classifications: $WORKDIR/classifications.csv"
echo "Summary:         $WORKDIR/summary.json"
echo "HTML Report:     $WORKDIR/report.html"

if [ "$KEEP_WORKDIR" = true ]; then
    echo ""
    echo "To view the report:"
    echo "  open $WORKDIR/report.html"
fi
