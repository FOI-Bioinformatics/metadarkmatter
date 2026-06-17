"""Regression test: streaming and in-memory classification filter identically.

Previously VectorizedClassifier.classify_file only applied alignment-quality
filters on the file-path branch, so the streaming path (which passes DataFrame
partitions) silently bypassed min_alignment_length etc. and classified more
reads than the in-memory path. These tests lock the two modes to the same
filtered read set.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from metadarkmatter.core.classification.ani_matrix import ANIMatrix
from metadarkmatter.core.classification.classifiers.vectorized import VectorizedClassifier
from metadarkmatter.models.config import ScoringConfig


@pytest.fixture
def blast_and_ani(tmp_path: Path) -> tuple[Path, ANIMatrix]:
    """A BLAST file where one read's only hit is below the length filter."""
    # columns: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    rows = [
        # read_long: 200 bp alignment -> kept under min_alignment_length=100
        ("read_long", "GCF_000001.1", 99.0, 200, 2, 0, 1, 200, 1, 200, 1e-90, 380.0),
        # read_short: 50 bp alignment -> dropped under min_alignment_length=100
        ("read_short", "GCF_000001.1", 98.0, 50, 1, 0, 1, 50, 1, 50, 1e-20, 95.0),
        # read_mid: 150 bp -> kept
        ("read_mid", "GCF_000002.1", 92.0, 150, 12, 0, 1, 150, 1, 150, 1e-55, 240.0),
    ]
    blast_path = tmp_path / "sample.blast.tsv"
    lines = ["\t".join(str(c) for c in r) for r in rows]
    blast_path.write_text("\n".join(lines) + "\n")

    ani_df = pl.DataFrame(
        {
            "genome": ["GCF_000001.1", "GCF_000002.1"],
            "GCF_000001.1": [100.0, 85.0],
            "GCF_000002.1": [85.0, 100.0],
        }
    )
    ani_path = tmp_path / "ani.csv"
    ani_df.write_csv(ani_path)
    return blast_path, ANIMatrix.from_file(ani_path)


def test_streaming_and_vectorized_filter_identically(
    blast_and_ani: tuple[Path, ANIMatrix], tmp_path: Path
):
    blast_path, ani = blast_and_ani
    config = ScoringConfig(min_alignment_length=100)
    clf = VectorizedClassifier(ani_matrix=ani, config=config)

    vec = clf.classify_file(blast_path)
    out = tmp_path / "streamed.parquet"
    # Tiny partitions to force the DataFrame-partition path.
    clf.stream_to_file(blast_path, out, output_format="parquet", partition_size=1)
    streamed = pl.read_parquet(out)

    vec_reads = set(vec["read_id"].to_list())
    streamed_reads = set(streamed["read_id"].to_list())

    # read_short must be filtered out in BOTH modes; the sets must match.
    assert "read_short" not in vec_reads
    assert "read_short" not in streamed_reads
    assert vec_reads == streamed_reads == {"read_long", "read_mid"}


def test_no_filter_keeps_all_reads_in_both_modes(
    blast_and_ani: tuple[Path, ANIMatrix], tmp_path: Path
):
    blast_path, ani = blast_and_ani
    # Default config has min_alignment_length=100; disable it to keep all reads.
    config = ScoringConfig(min_alignment_length=0, min_alignment_fraction=0.0)
    clf = VectorizedClassifier(ani_matrix=ani, config=config)

    vec = clf.classify_file(blast_path)
    out = tmp_path / "streamed.parquet"
    clf.stream_to_file(blast_path, out, output_format="parquet", partition_size=1)
    streamed = pl.read_parquet(out)

    assert set(vec["read_id"].to_list()) == set(streamed["read_id"].to_list())
    assert "read_short" in set(vec["read_id"].to_list())
