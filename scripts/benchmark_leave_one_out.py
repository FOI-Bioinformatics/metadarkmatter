#!/usr/bin/env python
"""
Leave-one-out novelty benchmark for the ANI-weighted classifier.

Idea: reads simulated from genome *g* are classified against a reference
set from which *g* has been removed (its hits are dropped and its row and
column are removed from the ANI matrix). The nearest remaining relative of
*g* then plays the role of the closest known genome of a novel organism,
and the true category is known from the full ANI matrix:

    nearest ANI >= 96  -> Known Species (another genome of the same species)
    nearest ANI >= 80  -> Novel Species
    otherwise          -> Novel Genus

A control run (nothing removed, truth "Known Species" for every read)
measures the false-novel rate.

Pipeline (existing outputs in --out are reused):

    1. simulate single-end reads from every genome (read IDs ``r{N}_{accession}``)
    2. ``mdm ani compute``            -> ani_matrix.csv
    3. ``mdm mmseqs2 makedb/search``  (or ``mdm blast makedb/align``)
    4. per genome: leave-one-out classification for each ScoringConfig
    5. metrics.json, per_genome.tsv, results.parquet and a printed summary

Usage:

    python scripts/benchmark_leave_one_out.py \
        --genomes genomes/ --out loo_bench/ --threads 8

    # Compare the default config with the clade-level Novel Species rule
    python scripts/benchmark_leave_one_out.py \
        --genomes genomes/ --out loo_bench/ --compare-clade-rule

    # Compare arbitrary ScoringConfig YAMLs
    python scripts/benchmark_leave_one_out.py \
        --genomes genomes/ --out loo_bench/ \
        --config-yaml base=configs/a.yaml --config-yaml tuned=configs/b.yaml

Genome files must be named ``{accession}.fna`` or ``{accession}.fna.gz``.
The simulated reads are error-prone substitutions only (no indels) and are
intended for relative comparison of configurations, not as a realistic
sequencing model.
"""

from __future__ import annotations

import argparse
import gzip
import json
import random
import subprocess
import sys
from pathlib import Path
from typing import Any

import polars as pl

KNOWN_ANI = 96.0
NOVEL_SPECIES_ANI = 80.0
BOUNDARY_LOW = 94.0
NOVEL_CALLS = ("Novel Species", "Novel Genus")
ACCESSION_RE = r"^r\d+_(.+)$"


# ---------------------------------------------------------------------------
# Read simulation
# ---------------------------------------------------------------------------


def read_genome_sequence(path: Path) -> str:
    """Return the concatenated upper-case sequence of a (gzipped) FASTA."""
    opener = gzip.open if path.suffix == ".gz" else open
    parts: list[str] = []
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                parts.append(line.strip().upper())
    return "".join(parts)


def accession_from_path(path: Path) -> str:
    """``GCF_1.1.fna.gz`` -> ``GCF_1.1``."""
    name = path.name
    for suffix in (".fna.gz", ".fna"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def mutate(seq: str, error_rate: float, rng: random.Random) -> str:
    """Substitute each base with a random base with probability error_rate."""
    if error_rate <= 0:
        return seq
    return "".join(
        rng.choice("ACGT") if rng.random() < error_rate else c for c in seq
    )


def simulate_reads(
    sequence: str,
    accession: str,
    n_reads: int,
    read_length: int,
    error_rate: float,
    rng: random.Random,
    start_index: int = 1,
) -> list[tuple[str, str]]:
    """Simulate single-end reads; returns ``(read_id, sequence)`` pairs.

    Read IDs are ``r{N}_{accession}``. Reads containing ambiguous bases
    are skipped (so fewer than ``n_reads`` may be returned); the counter
    only advances for emitted reads.
    """
    out: list[tuple[str, str]] = []
    if len(sequence) < read_length:
        return out
    idx = start_index
    for _ in range(n_reads):
        pos = rng.randrange(0, len(sequence) - read_length + 1)
        frag = sequence[pos : pos + read_length]
        read = mutate(frag, error_rate, rng)
        if set(read) - set("ACGT"):
            continue
        out.append((f"r{idx}_{accession}", read))
        idx += 1
    return out


def write_fastq(reads: list[tuple[str, str]], path: Path) -> None:
    """Write reads to a gzipped FASTQ with constant quality."""
    with gzip.open(path, "wt") as fh:
        for rid, seq in reads:
            fh.write(f"@{rid}\n{seq}\n+\n{'I' * len(seq)}\n")


def simulate_all(
    genome_files: dict[str, Path],
    reads_per_genome: int,
    read_length: int,
    error_rate: float,
    seed: int,
) -> list[tuple[str, str]]:
    """Simulate reads for every genome with one seeded RNG (deterministic)."""
    rng = random.Random(seed)
    reads: list[tuple[str, str]] = []
    for acc in sorted(genome_files):
        seq = read_genome_sequence(genome_files[acc])
        reads.extend(
            simulate_reads(
                seq, acc, reads_per_genome, read_length, error_rate, rng,
                start_index=len(reads) + 1,
            )
        )
    return reads


# ---------------------------------------------------------------------------
# Truth and leave-one-out helpers
# ---------------------------------------------------------------------------


def assign_truth(nearest_ani: float) -> str:
    """True category from the nearest remaining relative's ANI."""
    if nearest_ani >= KNOWN_ANI:
        return "Known Species"
    if nearest_ani >= NOVEL_SPECIES_ANI:
        return "Novel Species"
    return "Novel Genus"


def load_ani_dict(path: Path) -> dict[str, dict[str, float]]:
    """Load a square ANI CSV (first column = genome label) into nested dicts."""
    df = pl.read_csv(path)
    label = df.columns[0]
    cols = [c for c in df.columns if c != label]
    return {
        r[label]: {c: float(r[c]) for c in cols} for r in df.iter_rows(named=True)
    }


def ani_without(
    full: dict[str, dict[str, float]], genome: str
) -> dict[str, dict[str, float]]:
    """ANI matrix with ``genome`` removed (row and column)."""
    return {
        a: {b: v for b, v in row.items() if b != genome}
        for a, row in full.items()
        if a != genome
    }


def nearest_ani(full: dict[str, dict[str, float]], genome: str) -> float:
    """Highest ANI between ``genome`` and any other genome."""
    return max(v for b, v in full[genome].items() if b != genome)


def prepare_alignment(aln: pl.DataFrame) -> pl.DataFrame:
    """Add ``src`` (source genome of the read) and ``genome_name`` columns."""
    return aln.with_columns(
        pl.col("qseqid").str.extract(ACCESSION_RE, 1).alias("src"),
        pl.col("sseqid").str.extract(r"^([^|]+)\|", 1).alias("genome_name"),
    )


def leave_one_out_hits(aln: pl.DataFrame, genome: str) -> pl.DataFrame:
    """Hits of ``genome``'s reads with hits to ``genome`` itself dropped."""
    return aln.filter(
        (pl.col("src") == genome) & (pl.col("genome_name") != genome)
    ).drop("src")


# ---------------------------------------------------------------------------
# Configs
# ---------------------------------------------------------------------------


def build_configs(
    yaml_specs: list[str], compare_clade_rule: bool
) -> dict[str, Any]:
    """Return ordered ``{name: ScoringConfig}`` from CLI options."""
    from metadarkmatter.models.config import ScoringConfig

    configs: dict[str, Any] = {}
    if compare_clade_rule:
        if "clade_novel_species" not in ScoringConfig.model_fields:
            raise SystemExit(
                "ScoringConfig has no 'clade_novel_species' field in this checkout"
            )
        configs["default"] = ScoringConfig()
        configs["clade_rule"] = ScoringConfig(clade_novel_species=True)
    for spec in yaml_specs:
        name, sep, path = spec.partition("=")
        if not sep or not name or not path:
            raise SystemExit(f"--config-yaml expects NAME=PATH, got {spec!r}")
        if name in configs:
            raise SystemExit(f"duplicate config name {name!r}")
        configs[name] = ScoringConfig.from_yaml(Path(path))
    if not configs:
        configs["default"] = ScoringConfig()
    return configs


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------


def classify_loo(
    aln: pl.DataFrame,
    full_ani: dict[str, dict[str, float]],
    config: Any,
    config_name: str,
) -> pl.DataFrame:
    """Run leave-one-out for every genome plus the control for one config."""
    from metadarkmatter.core.classification.ani_matrix import ANIMatrix
    from metadarkmatter.core.classification.classifiers.vectorized import (
        VectorizedClassifier,
    )

    frames: list[pl.DataFrame] = []
    for g in sorted(full_ani):
        reads = leave_one_out_hits(aln, g)
        if reads.height == 0:
            continue
        near = nearest_ani(full_ani, g)
        res = VectorizedClassifier(
            ANIMatrix(ani_without(full_ani, g)), config=config
        ).classify_file(reads)
        frames.append(
            res.with_columns(
                pl.lit(g).alias("src"),
                pl.lit(near).alias("true_nearest_ani"),
                pl.lit(assign_truth(near)).alias("truth"),
                pl.lit("loo").alias("setting"),
            )
        )
    control_in = aln.drop("src")
    res = VectorizedClassifier(ANIMatrix(full_ani), config=config).classify_file(
        control_in
    )
    frames.append(
        res.with_columns(
            pl.col("read_id").str.extract(ACCESSION_RE, 1).alias("src"),
            pl.lit(100.0).alias("true_nearest_ani"),
            pl.lit("Known Species").alias("truth"),
            pl.lit("control").alias("setting"),
        )
    )
    return pl.concat(frames, how="diagonal_relaxed").with_columns(
        pl.lit(config_name).alias("config")
    )


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------


def _ratio(num: int, den: int) -> float | None:
    return num / den if den else None


def _frac_called(df: pl.DataFrame, calls: tuple[str, ...]) -> float | None:
    return _ratio(df.filter(pl.col("taxonomic_call").is_in(calls)).height, df.height)


def compute_metrics(results: pl.DataFrame) -> dict[str, dict[str, Any]]:
    """Per-config metrics from the per-read results frame.

    Required columns: config, setting ("loo"/"control"), truth,
    true_nearest_ani, taxonomic_call. Fractions with an empty denominator
    are ``None``.
    """
    metrics: dict[str, dict[str, Any]] = {}
    for cfg in results["config"].unique(maintain_order=True).to_list():
        d = results.filter(pl.col("config") == cfg)
        loo = d.filter(pl.col("setting") == "loo")
        control = d.filter(pl.col("setting") == "control")
        ani = pl.col("true_nearest_ani")
        ns_group = loo.filter(ani.is_between(NOVEL_SPECIES_ANI, BOUNDARY_LOW, closed="left"))
        boundary = loo.filter(ani.is_between(BOUNDARY_LOW, KNOWN_ANI, closed="left"))
        ng_group = loo.filter(pl.col("truth") == "Novel Genus")

        called_ns = d.filter(pl.col("taxonomic_call") == "Novel Species")
        called_novel = d.filter(pl.col("taxonomic_call").is_in(NOVEL_CALLS))
        truth_novel = d.filter(pl.col("truth") != "Known Species")

        metrics[cfg] = {
            "n_reads_loo": loo.height,
            "n_reads_control": control.height,
            "novel_species_recall": _frac_called(ns_group, ("Novel Species",)),
            "n_novel_species_reads": ns_group.height,
            "boundary_called_novel_species": _frac_called(boundary, ("Novel Species",)),
            "n_boundary_reads": boundary.height,
            "novel_genus_recall": _frac_called(ng_group, ("Novel Genus",)),
            "novel_genus_called_novel_species": _frac_called(ng_group, ("Novel Species",)),
            "n_novel_genus_reads": ng_group.height,
            "novel_species_precision": _ratio(
                called_ns.filter(pl.col("truth") == "Novel Species").height,
                called_ns.height,
            ),
            "any_novel_recall": _ratio(
                truth_novel.filter(pl.col("taxonomic_call").is_in(NOVEL_CALLS)).height,
                truth_novel.height,
            ),
            "any_novel_precision": _ratio(
                called_novel.filter(pl.col("truth") != "Known Species").height,
                called_novel.height,
            ),
            "control_false_novel": control.filter(
                pl.col("taxonomic_call").is_in(NOVEL_CALLS)
            ).height,
        }
    return metrics


def per_genome_table(results: pl.DataFrame) -> pl.DataFrame:
    """One row per genome: nearest ANI, truth, n reads, recall per config.

    Recall is the fraction of the genome's leave-one-out reads whose
    ``taxonomic_call`` equals the truth.
    """
    loo = results.filter(pl.col("setting") == "loo")
    if loo.height == 0:
        return pl.DataFrame()
    agg = (
        loo.group_by("config", "src")
        .agg(
            pl.col("true_nearest_ani").first(),
            pl.col("truth").first(),
            pl.len().alias("n_reads"),
            (pl.col("taxonomic_call") == pl.col("truth")).mean().alias("recall"),
        )
        .sort("config", "src")
    )
    wide = agg.pivot(on="config", index="src", values="recall")
    cfgs = agg["config"].unique(maintain_order=True).to_list()
    wide = wide.rename({c: f"recall_{c}" for c in cfgs})
    info = (
        agg.group_by("src")
        .agg(
            pl.col("true_nearest_ani").first(),
            pl.col("truth").first(),
            pl.col("n_reads").first(),
        )
    )
    return (
        info.join(wide, on="src")
        .rename({"src": "accession", "true_nearest_ani": "nearest_ani"})
        .sort("nearest_ani", "accession")
    )


def format_summary(metrics: dict[str, dict[str, Any]]) -> str:
    """Compact fixed-width table (metrics as rows, configs as columns)."""
    cfgs = list(metrics)
    rows = [
        ("novel_species_recall", "NS recall [80,94)"),
        ("boundary_called_novel_species", "boundary [94,96) -> NS"),
        ("novel_genus_recall", "NG recall"),
        ("novel_genus_called_novel_species", "NG called NS"),
        ("novel_species_precision", "NS precision"),
        ("any_novel_recall", "any-novel recall"),
        ("any_novel_precision", "any-novel precision"),
        ("control_false_novel", "control false-novel"),
    ]
    width = max(14, *(len(c) + 2 for c in cfgs))
    lines = ["metric".ljust(26) + "".join(c.rjust(width) for c in cfgs)]
    for key, label in rows:
        cells = []
        for c in cfgs:
            v = metrics[c][key]
            if v is None:
                s = "n/a"
            elif key == "control_false_novel":
                s = f"{v}/{metrics[c]['n_reads_control']}"
            else:
                s = f"{v:.3f}"
            cells.append(s.rjust(width))
        lines.append(label.ljust(26) + "".join(cells))
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# External steps
# ---------------------------------------------------------------------------


def run_mdm(args: list[str]) -> None:
    """Run a metadarkmatter CLI subcommand in a subprocess."""
    cmd = [sys.executable, "-m", "metadarkmatter.cli.main", *args]
    print("+ " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def find_genomes(genome_dir: Path) -> dict[str, Path]:
    files = sorted(list(genome_dir.glob("*.fna")) + list(genome_dir.glob("*.fna.gz")))
    return {accession_from_path(f): f for f in files}


def ensure_inputs(args: argparse.Namespace, genomes: dict[str, Path]) -> tuple[Path, Path]:
    """Create (or reuse) reads, ANI matrix, database and alignment."""
    out: Path = args.out
    out.mkdir(parents=True, exist_ok=True)
    reads_path = out / "reads.fastq.gz"
    ani_path = out / "ani_matrix.csv"
    aln_path = out / "alignment.tsv.gz"
    db_dir = out / "db"
    db_prefix = db_dir / "pangenome"

    if reads_path.exists():
        print(f"reusing {reads_path}")
    else:
        reads = simulate_all(
            genomes, args.reads_per_genome, args.read_length, args.error_rate, args.seed
        )
        write_fastq(reads, reads_path)
        print(f"simulated {len(reads)} reads -> {reads_path}")

    if ani_path.exists():
        print(f"reusing {ani_path}")
    else:
        run_mdm([
            "ani", "compute", "--genomes", str(args.genomes), "--output", str(ani_path),
            "--backend", args.ani_backend, "--threads", str(args.threads),
        ])

    if aln_path.exists():
        print(f"reusing {aln_path}")
        return ani_path, aln_path

    db_dir.mkdir(parents=True, exist_ok=True)
    tool = "mmseqs2" if args.aligner == "mmseqs2" else "blast"
    run_mdm([
        tool, "makedb", "--genomes", str(args.genomes), "--output", str(db_prefix),
        "--skip-if-exists",
    ])
    sub = "search" if tool == "mmseqs2" else "align"
    run_mdm([
        tool, sub, "--query", str(reads_path), "--database", str(db_prefix),
        "--output", str(aln_path), "--threads", str(args.threads),
    ])
    return ani_path, aln_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genomes", required=True, type=Path,
                        help="Directory of {accession}.fna[.gz] genome FASTAs.")
    parser.add_argument("--out", required=True, type=Path,
                        help="Output directory (existing intermediates are reused).")
    parser.add_argument("--reads-per-genome", type=int, default=400,
                        help="Reads simulated per genome (default 400).")
    parser.add_argument("--read-length", type=int, default=150,
                        help="Simulated read length (default 150).")
    parser.add_argument("--error-rate", type=float, default=0.005,
                        help="Per-base substitution rate (default 0.005).")
    parser.add_argument("--seed", type=int, default=42,
                        help="RNG seed for read simulation (default 42).")
    parser.add_argument("--threads", type=int, default=4,
                        help="Threads for ANI and alignment (default 4).")
    parser.add_argument("--ani-backend", default="fastani",
                        choices=["fastani", "skani", "auto"],
                        help="ANI backend (default fastani).")
    parser.add_argument("--aligner", default="mmseqs2", choices=["mmseqs2", "blast"],
                        help="Aligner (default mmseqs2).")
    parser.add_argument("--config-yaml", action="append", default=[], metavar="NAME=PATH",
                        help="ScoringConfig YAML to evaluate; repeatable.")
    parser.add_argument("--compare-clade-rule", action="store_true",
                        help="Compare ScoringConfig() with "
                             "ScoringConfig(clade_novel_species=True).")
    args = parser.parse_args()

    genomes = find_genomes(args.genomes)
    if len(genomes) < 3:
        raise SystemExit(f"need at least 3 genomes in {args.genomes}, found {len(genomes)}")
    configs = build_configs(args.config_yaml, args.compare_clade_rule)

    ani_path, aln_path = ensure_inputs(args, genomes)

    from metadarkmatter.core.parsers import StreamingBlastParser

    full_ani = load_ani_dict(ani_path)
    aln = prepare_alignment(StreamingBlastParser(aln_path).parse_eager())

    frames = []
    for name, cfg in configs.items():
        print(f"classifying with config {name!r} ...", flush=True)
        frames.append(classify_loo(aln, full_ani, cfg, name))
    results = pl.concat(frames, how="diagonal_relaxed")

    metrics = compute_metrics(results)
    table = per_genome_table(results)

    results.write_parquet(args.out / "results.parquet")
    (args.out / "metrics.json").write_text(json.dumps(metrics, indent=2))
    table.write_csv(args.out / "per_genome.tsv", separator="\t")

    print()
    print(format_summary(metrics))
    print(f"\nwrote results.parquet, metrics.json, per_genome.tsv to {args.out}")


if __name__ == "__main__":
    main()
