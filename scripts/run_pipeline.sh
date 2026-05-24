#!/usr/bin/env bash
# Canonical end-to-end metadarkmatter pipeline for an internal lab run.
#
# Invocation:
#   scripts/run_pipeline.sh <family> <reads_R1.fastq.gz> [<reads_R2.fastq.gz>] <kraken_db> <out_dir>
#
# Example:
#   scripts/run_pipeline.sh f__Francisellaceae \
#       data/sample_R1.fastq.gz data/sample_R2.fastq.gz \
#       /opt/db/kraken2 results/run01
#
# Steps:
#   1. Resolve the target taxid for the family from Kraken2 / GTDB.
#   2. Download reference genomes for the family from NCBI Datasets via
#      the metadarkmatter GTDB client (cached on disk).
#   3. Extract reads classified to the family by Kraken2.
#   4. Build a pangenome BLAST database from the genomes.
#   5. Align the extracted reads against the pangenome.
#   6. Compute the ANI matrix for the representative genomes.
#   7. Classify the alignment with the ANI-weighted Bayesian classifier.
#   8. Generate the HTML report.
#   9. Print a final summary.
#
# Notes:
# - Idempotent steps reuse cached outputs in $OUT_DIR; rm -rf to redo.
# - Set METADARKMATTER_SEED if you care about bit-reproducible
#   subsamples in the report; default 42.
# - The script aborts on the first error (set -e).

set -euo pipefail

if (( $# < 4 )); then
    echo "usage: $0 <family> <reads_R1.fastq.gz> [<reads_R2.fastq.gz>] <kraken_db> <out_dir>" >&2
    exit 2
fi

FAMILY="$1"
READS_R1="$2"
if (( $# >= 5 )); then
    READS_R2="$3"
    KRAKEN_DB="$4"
    OUT_DIR="$5"
else
    READS_R2=""
    KRAKEN_DB="$3"
    OUT_DIR="$4"
fi

mkdir -p "$OUT_DIR"/{genomes,blastdb,extraction,reports}

THREADS="${THREADS:-8}"
MDM="${MDM:-metadarkmatter}"

echo "[1/9] Family resolution: $FAMILY"
# The taxid for the family is consumed by Kraken2 in step 3. For
# Francisellaceae it is 34064; users running other families should
# look it up once and export FAMILY_TAXID.
: "${FAMILY_TAXID:?Set FAMILY_TAXID to the NCBI taxid for $FAMILY before running}"

echo "[2/9] Downloading reference genomes"
"$MDM" download genomes list "$FAMILY" --output "$OUT_DIR/genomes.tsv"
"$MDM" download genomes fetch \
    --accessions "$OUT_DIR/genomes.tsv" \
    --output-dir "$OUT_DIR/genomes/"

echo "[3/9] Extracting family reads with Kraken2"
KRAKEN_OUT="$OUT_DIR/sample.kraken"
if [[ ! -s "$KRAKEN_OUT" ]]; then
    if [[ -n "$READS_R2" ]]; then
        kraken2 --db "$KRAKEN_DB" --threads "$THREADS" --paired \
            --output "$KRAKEN_OUT" "$READS_R1" "$READS_R2"
    else
        kraken2 --db "$KRAKEN_DB" --threads "$THREADS" \
            --output "$KRAKEN_OUT" "$READS_R1"
    fi
fi
"$MDM" kraken2 extract \
    --kraken-output "$KRAKEN_OUT" \
    --reads-1 "$READS_R1" \
    ${READS_R2:+--reads-2 "$READS_R2"} \
    --taxid "$FAMILY_TAXID" \
    --output "$OUT_DIR/extraction/"

echo "[4/9] Building BLAST database"
"$MDM" blast makedb \
    --genomes "$OUT_DIR/genomes/" \
    --output "$OUT_DIR/blastdb/pangenome"

echo "[5/9] Aligning reads"
"$MDM" blast align \
    --query "$OUT_DIR/extraction/reads_R1.fastq.gz" \
    --database "$OUT_DIR/blastdb/pangenome" \
    --output "$OUT_DIR/sample.blast.tsv.gz" \
    --threads "$THREADS"

echo "[6/9] Computing ANI matrix"
"$MDM" ani compute \
    --genomes "$OUT_DIR/genomes/" \
    --output "$OUT_DIR/ani_matrix.csv" \
    --threads "$THREADS"

echo "[7/9] Classifying"
"$MDM" score classify \
    --alignment "$OUT_DIR/sample.blast.tsv.gz" \
    --ani "$OUT_DIR/ani_matrix.csv" \
    --output "$OUT_DIR/classifications.csv"

echo "[8/9] Generating report"
"$MDM" report generate \
    --classifications "$OUT_DIR/classifications.csv" \
    --ani "$OUT_DIR/ani_matrix.csv" \
    --output "$OUT_DIR/reports/report.html" \
    --report-mode offline

echo "[9/9] Done"
echo "    Classifications: $OUT_DIR/classifications.csv"
echo "    Report:          $OUT_DIR/reports/report.html"
