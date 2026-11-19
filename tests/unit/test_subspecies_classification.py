"""Classifier wiring for per-read NCBI species labels (closest_subspecies)."""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from metadarkmatter.core.classification.ani_matrix import ANIMatrix
from metadarkmatter.core.classification.classifiers.vectorized import VectorizedClassifier
from metadarkmatter.core.subspecies import UNRESOLVED

GENOME_A = "GCF_000001.1"  # labelled Brucella abortus
GENOME_S = "GCF_000002.1"  # labelled Brucella suis
GENOME_X = "GCF_000003.1"  # no label

SUBSPECIES_MAP = {GENOME_A: "Brucella abortus", GENOME_S: "Brucella suis"}


@pytest.fixture
def alignment_and_ani(tmp_path: Path) -> tuple[Path, ANIMatrix]:
    # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    rows = [
        # read_tie: equal best hits to A and S -> two labels -> UNRESOLVED
        ("read_tie", GENOME_A, 99.8, 200, 0, 0, 1, 200, 1, 200, 1e-100, 380.0),
        ("read_tie", GENOME_S, 99.8, 200, 0, 0, 1, 200, 1, 200, 1e-100, 380.0),
        # read_a: best hit A, weaker hit X -> Brucella abortus
        ("read_a", GENOME_A, 99.9, 200, 0, 0, 1, 200, 1, 200, 1e-100, 390.0),
        ("read_a", GENOME_X, 99.0, 200, 2, 0, 1, 200, 1, 200, 1e-90, 360.0),
        # read_x: only the unlabelled genome -> null
        ("read_x", GENOME_X, 99.9, 200, 0, 0, 1, 200, 1, 200, 1e-100, 390.0),
    ]
    # Subject ids use the standardized "{accession}|{contig}" form.
    rows = [(q, f"{g}|contig_1", *rest) for q, g, *rest in rows]
    aln = tmp_path / "sample.blast.tsv"
    aln.write_text("\n".join("\t".join(str(c) for c in r) for r in rows) + "\n")

    genomes = [GENOME_A, GENOME_S, GENOME_X]
    ani_df = pl.DataFrame(
        {"genome": genomes}
        | {g: [100.0 if g == h else 99.5 for h in genomes] for g in genomes}
    )
    ani_path = tmp_path / "ani.csv"
    ani_df.write_csv(ani_path)
    return aln, ANIMatrix.from_file(ani_path)


def test_classifier_adds_subspecies_columns(alignment_and_ani):
    aln, ani = alignment_and_ani
    clf = VectorizedClassifier(ani_matrix=ani, subspecies_map=SUBSPECIES_MAP)
    df = clf.classify_file(aln)
    assert isinstance(df, pl.DataFrame)

    assert df.columns[-2:] == ["closest_subspecies", "subspecies_candidates"]
    assert "top_hit_genomes" not in df.columns
    assert set(df["taxonomic_call"].to_list()) == {"Known Species"}

    rows = {r["read_id"]: r for r in df.to_dicts()}
    assert rows["read_tie"]["closest_subspecies"] == UNRESOLVED
    assert rows["read_tie"]["subspecies_candidates"] == "Brucella abortus;Brucella suis"
    assert rows["read_a"]["closest_subspecies"] == "Brucella abortus"
    assert rows["read_a"]["subspecies_candidates"] == "Brucella abortus"
    assert rows["read_x"]["closest_subspecies"] is None
    assert rows["read_x"]["subspecies_candidates"] is None


def test_classifier_without_subspecies_map_has_no_columns(alignment_and_ani):
    aln, ani = alignment_and_ani
    baseline = VectorizedClassifier(ani_matrix=ani).classify_file(aln)
    empty_map = VectorizedClassifier(ani_matrix=ani, subspecies_map={}).classify_file(aln)
    assert isinstance(baseline, pl.DataFrame)
    assert isinstance(empty_map, pl.DataFrame)

    for df in (baseline, empty_map):
        assert "closest_subspecies" not in df.columns
        assert "subspecies_candidates" not in df.columns
        assert "top_hit_genomes" not in df.columns
    assert baseline.equals(empty_map)


def test_streaming_matches_in_memory_subspecies(alignment_and_ani, tmp_path: Path):
    aln, ani = alignment_and_ani
    clf = VectorizedClassifier(ani_matrix=ani, subspecies_map=SUBSPECIES_MAP)
    in_memory = clf.classify_file(aln)
    assert isinstance(in_memory, pl.DataFrame)

    out = tmp_path / "streamed.csv"
    clf.stream_to_file(aln, out, output_format="csv", partition_size=2)
    streamed = pl.read_csv(out)

    assert streamed.columns[-2:] == ["closest_subspecies", "subspecies_candidates"]
    cols = ["read_id", "closest_subspecies", "subspecies_candidates"]
    assert streamed.select(cols).sort("read_id").equals(
        in_memory.select(cols).sort("read_id")
    )


def test_off_target_reads_have_null_subspecies(alignment_and_ani, tmp_path: Path):
    from metadarkmatter.models.config import ScoringConfig

    aln, ani = alignment_and_ani
    # read_off: best hit to a genome outside the ANI matrix -> Off-target
    with aln.open("a") as fh:
        fh.write(
            "read_off\tGCF_000009.1|contig_1\t99.9\t200\t0\t0\t1\t200\t1\t200\t1e-100\t390.0\n"
            "read_off\tGCF_000001.1|contig_1\t90.0\t200\t20\t0\t1\t200\t1\t200\t1e-50\t200.0\n"
        )
    config = ScoringConfig(target_family="f__Brucellaceae")
    clf = VectorizedClassifier(ani_matrix=ani, config=config, subspecies_map=SUBSPECIES_MAP)
    df = clf.classify_file(aln)
    assert isinstance(df, pl.DataFrame)

    assert df.columns[-2:] == ["closest_subspecies", "subspecies_candidates"]
    rows = {r["read_id"]: r for r in df.to_dicts()}
    assert rows["read_off"]["taxonomic_call"] == "Off-target"
    assert rows["read_off"]["closest_subspecies"] is None
    assert rows["read_off"]["subspecies_candidates"] is None
    assert rows["read_a"]["closest_subspecies"] == "Brucella abortus"

    # Off-target only: the early-return path also carries the columns.
    off_only = tmp_path / "off_only.tsv"
    off_only.write_text(
        "read_off\tGCF_000009.1|contig_1\t99.9\t200\t0\t0\t1\t200\t1\t200\t1e-100\t390.0\n"
    )
    df_off = clf.classify_file(off_only)
    assert isinstance(df_off, pl.DataFrame)
    assert df_off.columns[-2:] == ["closest_subspecies", "subspecies_candidates"]
    assert df_off["closest_subspecies"].to_list() == [None]


@pytest.mark.parametrize("streaming", [False, True])
def test_runner_passes_subspecies_map(alignment_and_ani, tmp_path: Path, streaming: bool):
    from metadarkmatter.core.classification.runner import (
        ClassificationRequest,
        run_classification,
    )
    from metadarkmatter.models.config import ScoringConfig

    aln, _ = alignment_and_ani
    metadata = tmp_path / "genome_metadata.tsv"
    metadata.write_text(
        "accession\tspecies\tgenus\tsubspecies\n"
        f"{GENOME_A}\tBrucella melitensis\tBrucella\tBrucella abortus\n"
        f"{GENOME_S}\tBrucella melitensis\tBrucella\tBrucella suis\n"
        f"{GENOME_X}\tBrucella melitensis\tBrucella\t\n"
    )
    output = tmp_path / "out.csv"
    messages: list[str] = []
    run_classification(
        ClassificationRequest(
            alignment=aln,
            ani=tmp_path / "ani.csv",
            output=output,
            config=ScoringConfig(),
            metadata=metadata,
            streaming=streaming,
        ),
        log=messages.append,
    )
    out = pl.read_csv(output)
    rows = {r["read_id"]: r for r in out.to_dicts()}
    assert rows["read_a"]["closest_subspecies"] == "Brucella abortus"
    assert rows["read_tie"]["closest_subspecies"] == UNRESOLVED
    if streaming:
        assert out.columns[-2:] == ["closest_subspecies", "subspecies_candidates"]
        assert "Subspecies columns" in "\n".join(messages)


def test_assign_subspecies_empty_map_returns_null_utf8():
    from metadarkmatter.core.subspecies import assign_subspecies

    df = pl.DataFrame(
        {"taxonomic_call": ["Known Species"], "top_hit_genomes": [[GENOME_A]]}
    )
    out = assign_subspecies(df, {})
    assert out.columns == ["taxonomic_call", "closest_subspecies", "subspecies_candidates"]
    assert out.schema["closest_subspecies"] == pl.Utf8
    assert out.schema["subspecies_candidates"] == pl.Utf8
    assert out["closest_subspecies"].to_list() == [None]


@pytest.fixture
def rep_only_ani(tmp_path: Path) -> tuple[ANIMatrix, dict[str, str]]:
    """ANI matrix with the representative only; S maps to A."""
    ani_path = tmp_path / "ani_rep.csv"
    pl.DataFrame({"genome": [GENOME_A, GENOME_X],
                  GENOME_A: [100.0, 99.0], GENOME_X: [99.0, 100.0]}).write_csv(ani_path)
    return ANIMatrix.from_file(ani_path), {GENOME_A: GENOME_A, GENOME_S: GENOME_A}


def _write_aln(path: Path, rows: list[tuple]) -> Path:
    rows = [(q, f"{g}|contig_1", *rest) for q, g, *rest in rows]
    path.write_text("\n".join("\t".join(str(c) for c in r) for r in rows) + "\n")
    return path


def test_non_representative_best_hit_labelled_with_representative_mapping(
    tmp_path: Path, rep_only_ani
):
    ani, rep_map = rep_only_ani
    aln = _write_aln(
        tmp_path / "s.blast.tsv",
        [
            # single best hit is the non-representative S
            ("read_s", GENOME_S, 99.9, 200, 0, 0, 1, 200, 1, 200, 1e-100, 390.0),
            ("read_s", GENOME_X, 98.0, 200, 4, 0, 1, 200, 1, 200, 1e-80, 340.0),
            # tie between A and S with different labels
            ("read_tie", GENOME_A, 99.8, 200, 0, 0, 1, 200, 1, 200, 1e-100, 380.0),
            ("read_tie", GENOME_S, 99.8, 200, 0, 0, 1, 200, 1, 200, 1e-100, 380.0),
        ],
    )
    clf = VectorizedClassifier(
        ani_matrix=ani, representative_mapping=rep_map, subspecies_map=SUBSPECIES_MAP
    )
    df = clf.classify_file(aln)
    assert isinstance(df, pl.DataFrame)
    rows = {r["read_id"]: r for r in df.to_dicts()}
    assert rows["read_s"]["closest_subspecies"] == "Brucella suis"
    assert rows["read_tie"]["closest_subspecies"] == UNRESOLVED

    baseline = VectorizedClassifier(ani_matrix=ani, representative_mapping=rep_map)
    base_df = baseline.classify_file(aln)
    assert isinstance(base_df, pl.DataFrame)
    core = [c for c in df.columns if c not in ("closest_subspecies", "subspecies_candidates")]
    assert df.select(core).equals(base_df.select(core))
