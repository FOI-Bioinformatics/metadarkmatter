#!/usr/bin/env bash
# End-to-end per-read-labelled Bayesian calibration experiment.
#
# Wires the existing per-step scripts into a single invocation so the
# operator can reproduce the experiment described in
# docs/STATISTICS_AUDIT.md (Step 3): build a labelled corpus from one
# family, simulate reads from the target genomes, run the standard
# classifier against the held-out reference set, attach per-read truth
# from the genome-vs-source ANI, fit a calibrated YAML, then evaluate.
#
# Usage:
#   scripts/run_per_read_calibration.sh <family> <out_dir> [threads]
#
# Example:
#   scripts/run_per_read_calibration.sh f__Francisellaceae \
#       corpora/francisella 16
#
# Inputs:
#   - InSilicoSeq, BLAST+, fastANI or skani must be on PATH.
#   - METADARKMATTER_SEED controls every random draw (default 42).
#
# Outputs (under <out_dir>):
#   - reads_pooled_R1.fastq.gz / R2 - pooled simulated reads
#   - classifications_default.csv - baseline run (hand-tuned defaults)
#   - metrics_per_read.tsv - per-read truth + classifier output
#   - bayesian_per_read.yaml - fitted config (not in configs/, by design)
#   - classifications_calibrated.csv - re-run with the fitted config
#   - eval_default.json / eval_calibrated.json - score evaluate output
#   - summary.txt - one-line accuracy delta
#
# The shipped hand-tuned defaults remain the reference; this script
# produces the artefacts needed to decide whether a calibrated YAML
# clears the cross-family gate in docs/CALIBRATION_RESULTS.md before
# anyone proposes shipping one.

set -euo pipefail

if (( $# < 2 )); then
    echo "usage: $0 <family> <out_dir> [threads]" >&2
    exit 2
fi

FAMILY="$1"
OUT_DIR="$2"
THREADS="${3:-8}"
SEED="${METADARKMATTER_SEED:-42}"

mkdir -p "$OUT_DIR"

echo "[1/7] Build corpus (download + ANI + partition + emit simulate.sh)"
if [ ! -f "$OUT_DIR/target_to_label.tsv" ]; then
    uv run python scripts/build_corpus.py \
        --family "$FAMILY" \
        --output-dir "$OUT_DIR" \
        --threads "$THREADS"
fi

echo "[2/7] Simulate reads (InSilicoSeq via emitted simulate.sh)"
if [ ! -f "$OUT_DIR/reads_pooled_R1.fastq.gz" ]; then
    bash "$OUT_DIR/simulate.sh"
fi

REF_DB="$OUT_DIR/reference_blastdb/pangenome"
echo "[3/7] Build reference BLAST DB from held-out reference set"
if [ ! -f "${REF_DB}.nhr" ]; then
    REF_DIR="$OUT_DIR/reference_genomes_dir"
    mkdir -p "$REF_DIR"
    while read -r acc; do
        src="$OUT_DIR/genomes/${acc}.fna.gz"
        [ -f "$src" ] && cp "$src" "$REF_DIR/"
    done < "$OUT_DIR/reference_genomes.txt"
    mkdir -p "$(dirname "$REF_DB")"
    uv run metadarkmatter blast makedb \
        --genomes "$REF_DIR" --output "$REF_DB"
fi

echo "[4/7] Align simulated reads against reference DB (single-end R1)"
ALIGN="$OUT_DIR/alignment.tsv.gz"
if [ ! -f "$ALIGN" ]; then
    uv run metadarkmatter blast align \
        --query "$OUT_DIR/reads_pooled_R1.fastq.gz" \
        --database "$REF_DB" \
        --output "$ALIGN" --threads "$THREADS"
fi

REF_ANI="$OUT_DIR/ani_reference.csv"
echo "[5/7] Compute reference-only ANI matrix for classification"
if [ ! -f "$REF_ANI" ]; then
    uv run metadarkmatter ani compute \
        --genomes "$OUT_DIR/reference_genomes_dir" \
        --output "$REF_ANI" --threads "$THREADS"
fi

echo "[6/7] Baseline classify (hand-tuned defaults)"
BASE_CSV="$OUT_DIR/classifications_default.csv"
uv run metadarkmatter score classify \
    --alignment "$ALIGN" --ani "$REF_ANI" \
    --output "$BASE_CSV" --strict-ani

echo "       Attach per-read truth (per_read mode uses ANI to source)"
METRICS="$OUT_DIR/metrics_per_read.tsv"
uv run python scripts/build_metrics_tsv.py \
    --classifications "$BASE_CSV" \
    --labels "$OUT_DIR/target_to_label.tsv" \
    --label-mode per_read \
    --ani-matrix "$OUT_DIR/ani_full.csv" \
    --output "$METRICS"

echo "[7/7] Fit per-read Bayesian config (outside configs/ - no guardrail trip)"
FIT_YAML="$OUT_DIR/bayesian_per_read.yaml"
uv run python scripts/calibrate_bayesian.py \
    --benchmark "$METRICS" \
    --output "$FIT_YAML" \
    --label-mode per_read

echo "       Re-classify with fitted config and evaluate both"
CAL_CSV="$OUT_DIR/classifications_calibrated.csv"
uv run metadarkmatter score classify \
    --alignment "$ALIGN" --ani "$REF_ANI" \
    --config "$FIT_YAML" \
    --output "$CAL_CSV" --strict-ani

# build_metrics_tsv to attach truth to the calibrated run too
CAL_METRICS="$OUT_DIR/metrics_per_read_calibrated.tsv"
uv run python scripts/build_metrics_tsv.py \
    --classifications "$CAL_CSV" \
    --labels "$OUT_DIR/target_to_label.tsv" \
    --label-mode per_read \
    --ani-matrix "$OUT_DIR/ani_full.csv" \
    --output "$CAL_METRICS"

uv run metadarkmatter score evaluate \
    --predictions "$METRICS" --truth-column true_category \
    --output "$OUT_DIR/eval_default.json"
uv run metadarkmatter score evaluate \
    --predictions "$CAL_METRICS" --truth-column true_category \
    --output "$OUT_DIR/eval_calibrated.json"

python - <<'PYEOF' "$OUT_DIR/eval_default.json" "$OUT_DIR/eval_calibrated.json" \
    > "$OUT_DIR/summary.txt"
import json, sys
d = json.loads(open(sys.argv[1]).read())
c = json.loads(open(sys.argv[2]).read())
ad = d.get("accuracy", float("nan"))
ac = c.get("accuracy", float("nan"))
print(f"default accuracy:    {ad:.4f}")
print(f"calibrated accuracy: {ac:.4f}")
print(f"delta:               {ac - ad:+.4f}")
print()
print("Ship the calibrated YAML only if this delta is >= 0 here AND")
print("on at least one other family per docs/CALIBRATION_RESULTS.md,")
print("with no family worse than -0.02.")
PYEOF

cat "$OUT_DIR/summary.txt"
echo
echo "All artefacts in $OUT_DIR. Per-read fitted YAML is at $FIT_YAML."
echo "Do NOT copy it into configs/ without running the cross-family gate."
