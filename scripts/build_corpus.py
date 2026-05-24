#!/usr/bin/env python
"""
Build a labelled benchmark corpus for Bayesian calibration.

Phase 1 of the calibration validation plan. For one family:

  1. Download GTDB representative genomes.
  2. Compute the full pairwise ANI matrix.
  3. Partition genomes 75/25 into a reference set and a target set,
     stratifying so every ANI quintile is represented in both.
  4. Assign a true_category to every target genome based on its
     nearest-reference ANI (the labels are derived from genome-level
     ANI, never from running the classifier under test).
  5. Emit a shell script (simulate.sh) that, when run, calls
     InSilicoSeq to synthesise paired-end reads from each target
     genome and pools them with read headers tagged by source
     accession.

The genome download, ANI computation, and BLAST DB build are invoked
via the metadarkmatter CLI. Read simulation is emitted as a separate
shell script because it depends on InSilicoSeq being installed and is
the slowest step; the lab can substitute mason/wgsim by editing the
emitted script.

Usage:

    python scripts/build_corpus.py \\
        --family f__Francisellaceae \\
        --output-dir corpora/francisella \\
        --threads 16

The output directory layout after the Python phase:

    corpora/francisella/
        genomes.tsv                 # full accession list
        genomes/                    # downloaded FASTAs
        ani_full.csv                # full-family ANI matrix
        reference_genomes.txt       # 75% subset
        target_genomes.txt          # 25% subset
        target_to_label.tsv         # per-target ground truth label
        simulate.sh                 # next-step read simulation script
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import polars as pl

logger = logging.getLogger(__name__)

# Label thresholds. These mirror the package defaults and are written here
# (not imported) so the corpus labels stay fixed even if the package
# thresholds are tuned later.
ANI_KNOWN_SPECIES = 96.0
ANI_NOVEL_SPECIES = 80.0
ANI_NOVEL_GENUS = 75.0
# A target is labelled "Species Boundary" when its best and second-best
# reference ANIs are within this band - signals ambiguous placement.
SPECIES_BOUNDARY_GAP = 1.5
AMBIGUOUS_GAP = 0.5  # Within this gap and at low ANI: Ambiguous

DEFAULT_REFERENCE_FRACTION = 0.75
DEFAULT_QUINTILES = 5


def _run(cmd: list[str]) -> None:
    """Run a subprocess command, streaming output to the parent terminal."""
    logger.info("$ %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def download_and_fetch_genomes(family: str, out_dir: Path, threads: int) -> Path:
    """Step 1.1-1.2: list and fetch the family's GTDB representatives."""
    out_dir.mkdir(parents=True, exist_ok=True)
    accession_tsv = out_dir / "genomes.tsv"
    genome_dir = out_dir / "genomes"

    if accession_tsv.exists() and genome_dir.exists():
        logger.info("Skipping download (%s already exists).", accession_tsv)
        return genome_dir

    _run([
        "metadarkmatter", "download", "genomes", "list", family,
        "--output", str(accession_tsv),
    ])
    _run([
        "metadarkmatter", "download", "genomes", "fetch",
        "--accessions", str(accession_tsv),
        "--output-dir", str(genome_dir),
    ])
    return genome_dir


def compute_full_ani(genome_dir: Path, out_dir: Path, threads: int) -> Path:
    """Step 1.2 continued: full-family ANI matrix."""
    ani_path = out_dir / "ani_full.csv"
    if ani_path.exists():
        logger.info("Skipping ANI compute (%s already exists).", ani_path)
        return ani_path
    _run([
        "metadarkmatter", "ani", "compute",
        "--genomes", str(genome_dir),
        "--output", str(ani_path),
        "--threads", str(threads),
        "--fastani-batches", "4",
    ])
    return ani_path


def _load_ani_matrix(ani_path: Path) -> tuple[list[str], np.ndarray]:
    """Return (genomes, ANI matrix) from a metadarkmatter ANI CSV."""
    df = pl.read_csv(ani_path)
    genomes = df.columns[1:]
    matrix = df.select(genomes).to_numpy()
    return list(genomes), matrix


def stratified_split(
    genomes: list[str],
    matrix: np.ndarray,
    reference_fraction: float,
    quintiles: int,
    rng: np.random.Generator,
) -> tuple[list[str], list[str]]:
    """Partition genomes into reference/target stratifying by ANI rank.

    For each genome g, compute its mean ANI to all other genomes; bin
    those means into ``quintiles`` strata; sample ``reference_fraction``
    of each stratum into the reference set, the remainder into target.
    Stratification prevents the target set from being dominated by
    only-close or only-far genomes.
    """
    n = len(genomes)
    mean_ani = np.zeros(n)
    for i in range(n):
        others = np.delete(matrix[i], i)
        # Treat 0 (no signal from fastANI) as missing for the stratification.
        mask = others > 0
        mean_ani[i] = float(others[mask].mean()) if mask.any() else 0.0

    # Quantile-bin into strata.
    edges = np.quantile(mean_ani, np.linspace(0, 1, quintiles + 1))
    edges[0] -= 1e-9  # ensure the lowest value lands in bin 0

    references: list[str] = []
    targets: list[str] = []
    for q in range(quintiles):
        in_bin = [
            g
            for g, a in zip(genomes, mean_ani, strict=True)
            if edges[q] < a <= edges[q + 1]
        ]
        if not in_bin:
            continue
        rng.shuffle(in_bin)
        cut = max(1, int(round(len(in_bin) * reference_fraction)))
        references.extend(in_bin[:cut])
        targets.extend(in_bin[cut:])

    # Guarantee at least one target if rounding ate them all.
    if not targets and references:
        targets.append(references.pop())
    return sorted(references), sorted(targets)


def categorise(best_ani: float, second_ani: float | None) -> str:
    """Map nearest-reference ANI to one of the six Bayesian categories.

    The label is derived from genome-level ANI alone, independent of
    the classifier under test. Species Boundary fires when the second-
    best reference is very close to the best (ambiguous placement) and
    the best is itself near the species line; Ambiguous fires at lower
    ANI with a similar gap pattern.
    """
    gap = (best_ani - second_ani) if second_ani is not None else float("inf")
    near_boundary = abs(best_ani - ANI_KNOWN_SPECIES) <= 2.0

    if best_ani >= ANI_KNOWN_SPECIES:
        if near_boundary and gap < SPECIES_BOUNDARY_GAP:
            return "Species Boundary"
        return "Known Species"
    if best_ani >= ANI_NOVEL_SPECIES:
        if gap < AMBIGUOUS_GAP and best_ani < ANI_KNOWN_SPECIES - 2.0:
            return "Ambiguous"
        return "Novel Species"
    if best_ani >= ANI_NOVEL_GENUS:
        return "Novel Genus"
    return "Unclassified"


def assign_labels(
    references: list[str],
    targets: list[str],
    genomes: list[str],
    matrix: np.ndarray,
) -> pl.DataFrame:
    """For each target, find nearest reference ANI and label the row."""
    idx = {g: i for i, g in enumerate(genomes)}
    rows = []
    for target in targets:
        t_i = idx[target]
        ref_anis = sorted(
            ((r, float(matrix[t_i, idx[r]])) for r in references if matrix[t_i, idx[r]] > 0),
            key=lambda kv: kv[1],
            reverse=True,
        )
        if not ref_anis:
            best_ani = 0.0
            second_ani = None
            nearest = "none"
        else:
            nearest, best_ani = ref_anis[0]
            second_ani = ref_anis[1][1] if len(ref_anis) > 1 else None
        label = categorise(best_ani, second_ani)
        rows.append(
            {
                "target_accession": target,
                "nearest_ref_accession": nearest,
                "nearest_ref_ani": best_ani,
                "second_ref_ani": second_ani if second_ani is not None else float("nan"),
                "true_category": label,
            }
        )
    return pl.DataFrame(rows)


def emit_simulate_script(
    out_dir: Path,
    genome_dir: Path,
    targets: list[str],
    depth: int,
    seed: int,
) -> Path:
    """Write a shell script that synthesises tagged reads with InSilicoSeq."""
    script = out_dir / "simulate.sh"
    target_dir = out_dir / "simulated"
    pooled_r1 = out_dir / "reads_pooled_R1.fastq.gz"
    pooled_r2 = out_dir / "reads_pooled_R2.fastq.gz"

    lines = [
        "#!/usr/bin/env bash",
        "# Phase 1.4: simulate paired-end reads per target genome and pool.",
        "# Requires InSilicoSeq on PATH (conda install -c bioconda insilicoseq).",
        "# Replace iss with mason or wgsim by editing the per-target block.",
        "set -euo pipefail",
        f"mkdir -p {target_dir}",
        "",
    ]
    for t in targets:
        # The genome file may have any of several extensions; the
        # download CLI writes .fna.gz by default.
        lines.extend(
            [
                f"echo '[simulate] {t}'",
                f"target_fa={genome_dir}/{t}.fna",
                f"if [ ! -f \"$target_fa\" ] && [ -f \"$target_fa.gz\" ]; then",
                f"    gunzip -k \"$target_fa.gz\"",
                "fi",
                f"iss generate \\",
                f"    --genomes \"$target_fa\" \\",
                f"    --model hiseq \\",
                f"    --n_reads {depth * 1000} \\",
                f"    --seed {seed} \\",
                f"    --output {target_dir}/{t}",
                "",
                "# Tag each read with the source accession (read_id prefix)",
                f"for mate in R1 R2; do",
                f"    awk -v acc='{t}' 'NR%4==1{{sub(/^@/, \"@\"acc\"__\")}}1' \\",
                f"        {target_dir}/{t}_${{mate}}.fastq \\",
                f"        | gzip >> {out_dir}/reads_pooled_${{mate}}.fastq.gz",
                "    rm -f " + f"{target_dir}/{t}_${{mate}}.fastq",
                "done",
                "",
            ]
        )
    lines.extend(
        [
            f"echo 'Pooled reads written to {pooled_r1} and {pooled_r2}'",
            "",
        ]
    )
    script.write_text("\n".join(lines))
    script.chmod(0o755)
    return script


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family", required=True, help="GTDB family, e.g. f__Francisellaceae")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument(
        "--reference-fraction",
        type=float,
        default=DEFAULT_REFERENCE_FRACTION,
        help="Fraction of genomes used as reference set (default 0.75).",
    )
    parser.add_argument(
        "--quintiles",
        type=int,
        default=DEFAULT_QUINTILES,
        help="Strata to balance the split across (default 5).",
    )
    parser.add_argument(
        "--depth",
        type=int,
        default=10,
        help="Approximate read depth per target genome (default 10).",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Reuse existing genomes/ANI under --output-dir.",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    if not args.skip_download:
        if shutil.which("metadarkmatter") is None:
            logger.error("metadarkmatter CLI not on PATH; install the package first.")
            return 2
        genome_dir = download_and_fetch_genomes(args.family, out_dir, args.threads)
        ani_path = compute_full_ani(genome_dir, out_dir, args.threads)
    else:
        genome_dir = out_dir / "genomes"
        ani_path = out_dir / "ani_full.csv"
        if not ani_path.exists():
            logger.error("--skip-download passed but %s is missing.", ani_path)
            return 2

    genomes, matrix = _load_ani_matrix(ani_path)
    logger.info("Loaded %d genomes from %s.", len(genomes), ani_path)

    rng = np.random.default_rng(args.seed)
    references, targets = stratified_split(
        genomes, matrix, args.reference_fraction, args.quintiles, rng
    )
    logger.info("Split: %d reference, %d target.", len(references), len(targets))

    (out_dir / "reference_genomes.txt").write_text("\n".join(references) + "\n")
    (out_dir / "target_genomes.txt").write_text("\n".join(targets) + "\n")

    labels = assign_labels(references, targets, genomes, matrix)
    labels.write_csv(out_dir / "target_to_label.tsv", separator="\t")
    logger.info(
        "Wrote %d target labels:\n%s",
        len(labels),
        labels.group_by("true_category").len().sort("true_category"),
    )

    script = emit_simulate_script(
        out_dir, genome_dir, targets, depth=args.depth, seed=args.seed
    )
    logger.info(
        "Next: run %s to synthesise reads, then continue with the "
        "alignment/classification pipeline. See docs/CALIBRATION.md.",
        script,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
