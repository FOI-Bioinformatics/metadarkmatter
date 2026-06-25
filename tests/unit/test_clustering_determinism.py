"""Regression tests: novel-cluster assignment must be order-independent.

Polars ``group_by``/``unique`` do not guarantee row order, so sequential
cluster IDs and neighborhood representatives must not depend on input order
or process hash seed. These tests shuffle the inputs and assert the
cluster_id -> genome mapping is stable.
"""

from __future__ import annotations

import polars as pl

from metadarkmatter.core.classification.ani_matrix import ANIMatrix
from metadarkmatter.core.novel_diversity import NovelDiversityAnalyzer
from metadarkmatter.core.phylogeny.placement import extract_novel_clusters


def _novel_df() -> pl.DataFrame:
    """Novel reads where two genomes tie on read count (the bug trigger)."""
    rows = []
    # NSP on three genomes: GCF_b and GCF_c tie at 8 reads, GCF_a has 9.
    spec = [("GCF_000a.1", 9), ("GCF_000b.1", 8), ("GCF_000c.1", 8)]
    rid = 0
    for genome, n in spec:
        for _ in range(n):
            rid += 1
            rows.append({
                "read_id": f"r{rid}",
                "best_match_genome": genome,
                "top_hit_identity": 92.0,
                "novelty_index": 8.0,
                "placement_uncertainty": 0.5,
                "taxonomic_call": "Novel Species",
            })
    return pl.DataFrame(rows)


def test_extract_novel_clusters_id_assignment_is_order_independent() -> None:
    """Shuffling input rows must not change which genome gets which cluster ID."""
    df = _novel_df()

    def id_map(frame: pl.DataFrame) -> dict[str, str]:
        clusters = extract_novel_clusters(frame, min_reads=3)
        return {c.cluster_id: c.best_match_genome for c in clusters}

    baseline = id_map(df)
    assert baseline  # non-empty
    # Several deterministic reorderings of the same rows.
    for seed in range(5):
        shuffled = df.sample(fraction=1.0, shuffle=True, seed=seed)
        assert id_map(shuffled) == baseline


def test_neighborhood_representative_is_order_independent() -> None:
    """Tied read counts must resolve to the same representative every time."""
    # Two genomes tie on novel-read count; ANI links them into one neighborhood.
    ani = pl.DataFrame({
        "genome": ["GCF_000b.1", "GCF_000c.1"],
        "GCF_000b.1": [100.0, 99.0],
        "GCF_000c.1": [99.0, 100.0],
    })
    ani_matrix = ANIMatrix.from_dataframe(ani)
    rows = []
    rid = 0
    for genome in ("GCF_000b.1", "GCF_000c.1"):
        for _ in range(8):
            rid += 1
            rows.append({
                "read_id": f"r{rid}",
                "best_match_genome": genome,
                "top_hit_identity": 92.0,
                "novelty_index": 8.0,
                "placement_uncertainty": 0.5,
                "taxonomic_call": "Novel Species",
            })
    df = pl.DataFrame(rows)

    reps = set()
    for seed in range(5):
        shuffled = df.sample(fraction=1.0, shuffle=True, seed=seed)
        analyzer = NovelDiversityAnalyzer(
            classifications=shuffled, metadata=None, ani_matrix=ani_matrix,
            novelty_band_size=5.0, min_cluster_size=3,
        )
        clusters = analyzer.cluster_novel_reads()
        reps.add(tuple(sorted(c.nearest_genome for c in clusters)))
    assert len(reps) == 1, f"representative selection not deterministic: {reps}"
