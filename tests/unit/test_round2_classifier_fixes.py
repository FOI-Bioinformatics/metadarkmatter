"""Regression tests for the round-2 audit fixes (classification core)."""

from __future__ import annotations

import logging
from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.main import app
from metadarkmatter.core.classification import ANIMatrix, VectorizedClassifier
from metadarkmatter.core.classification.evaluation import evaluate_classifications
from metadarkmatter.core.classification.runner import (
    ClassificationRequest,
    run_classification,
)
from metadarkmatter.core.exceptions import ConfigurationError
from metadarkmatter.core.metadata import GenomeMetadata
from metadarkmatter.models.config import ScoringConfig

COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
    "qend", "sstart", "send", "evalue", "bitscore", "qlen",
]
GENOMES = ["A", "B", "C"]


def mkdf(rows, length=150):
    df = pl.DataFrame(
        [(q, s, p, length, 0, 0, 1, 150, 1, 150, 1e-30, b, 150) for q, s, p, b in rows],
        schema=COLS,
        orient="row",
    )
    return df.with_columns(pl.col("sseqid").str.extract(r"^([^|]+)\|", 1).alias("genome_name"))


def write_blast(path: Path, rows) -> Path:
    lines = [
        "\t".join(str(x) for x in (q, s, p, 150, 0, 0, 1, 150, 1, 150, "1e-30", b, 150))
        for q, s, p, b in rows
    ]
    path.write_text("\n".join(lines) + "\n")
    return path


def write_ani(path: Path, genomes=GENOMES, vals=None) -> Path:
    vals = vals or {("A", "B"): 90.0, ("A", "C"): 85.0, ("B", "C"): 88.0}
    rows = []
    for g in genomes:
        row = {"genome": g}
        for h in genomes:
            row[h] = 100.0 if g == h else vals.get((g, h), vals.get((h, g), 0.0))
        rows.append(row)
    pl.DataFrame(rows).write_csv(path)
    return path


@pytest.fixture
def ani():
    return ANIMatrix({"A": {"B": 90.0, "C": 85.0}, "B": {"A": 90.0, "C": 88.0},
                      "C": {"A": 85.0, "B": 88.0}})


@pytest.fixture
def cfg():
    return ScoringConfig(min_alignment_length=0, min_alignment_fraction=0.0)


# 1. protein-mode alignment length ------------------------------------------
def test_protein_mode_length_filter_uses_nucleotide_span(ani):
    # BLASTX: length in amino acids (50), query span 150 nt.
    df = mkdf([("r1", "A|c", 80.0, 100.0)], length=50)
    prot = ScoringConfig(alignment_mode="protein", min_alignment_fraction=0.0)
    assert len(VectorizedClassifier(ani, config=prot).classify_file(df)) == 1
    nucl = ScoringConfig(min_alignment_fraction=0.0)
    assert len(VectorizedClassifier(ani, config=nucl).classify_file(df)) == 0


def test_all_rows_filtered_warns(ani, caplog):
    df = mkdf([("r1", "A|c", 80.0, 100.0)], length=50)
    with caplog.at_level(logging.WARNING):
        out = VectorizedClassifier(ani, config=ScoringConfig(min_alignment_fraction=0.0)).classify_file(df)
    assert len(out) == 0
    assert "min_alignment_length=100" in caplog.text


# 2. zero classified ---------------------------------------------------------
def test_zero_classified_writes_header_only_and_warns(tmp_path):
    blast = write_blast(tmp_path / "a.tsv", [("r1", "A|c", 99.0, 250.0)])
    ani_p = write_ani(tmp_path / "ani.csv")
    out = tmp_path / "out.csv"
    out.write_text("stale\n")
    cfg = ScoringConfig(min_alignment_length=1000)
    msgs: list[str] = []
    res = run_classification(
        ClassificationRequest(alignment=blast, ani=ani_p, output=out, config=cfg),
        log=msgs.append,
    )
    assert res.num_classified == 0
    header = out.read_text().strip()
    assert header.startswith("read_id,") and "stale" not in header
    assert any("0 reads were classified" in m for m in msgs)


# 3. subject IDs without '|' -------------------------------------------------
def test_ids_without_pipe_raise_configuration_error(tmp_path):
    blast = write_blast(tmp_path / "a.tsv", [("r1", "A_contig", 99.0, 250.0)])
    ani_p = write_ani(tmp_path / "ani.csv")
    req = ClassificationRequest(
        alignment=blast, ani=ani_p, output=tmp_path / "o.csv",
        config=ScoringConfig(min_alignment_length=0, min_alignment_fraction=0.0),
    )
    with pytest.raises(ConfigurationError, match="unknown"):
        run_classification(req)


# D. 0% overlap --------------------------------------------------------------
def test_zero_overlap_with_ani_matrix_is_error(tmp_path):
    blast = write_blast(tmp_path / "a.tsv", [("r1", "Z|c", 99.0, 250.0)])
    ani_p = write_ani(tmp_path / "ani.csv")
    req = ClassificationRequest(
        alignment=blast, ani=ani_p, output=tmp_path / "o.csv",
        config=ScoringConfig(min_alignment_length=0, min_alignment_fraction=0.0),
    )
    with pytest.raises(ConfigurationError, match="0% overlap"):
        run_classification(req)


# 4. missing ANI genomes warning ---------------------------------------------
def test_missing_ani_genomes_reported(tmp_path, cfg, ani):
    blast = write_blast(
        tmp_path / "a.tsv",
        [("r1", "A|c", 99.0, 250.0), ("r1", "Z|c", 98.5, 248.0),
         ("r2", "B|c", 99.0, 250.0)],
    )
    ani_p = write_ani(tmp_path / "ani.csv")
    msgs: list[str] = []
    run_classification(
        ClassificationRequest(alignment=blast, ani=ani_p, output=tmp_path / "o.csv", config=cfg),
        log=msgs.append,
    )
    joined = "\n".join(msgs)
    assert "absent from the ANI matrix" in joined
    assert "1 hit genome(s)" in joined


def test_no_missing_ani_warning_when_all_present(tmp_path, cfg):
    blast = write_blast(tmp_path / "a.tsv", [("r1", "A|c", 99.0, 250.0), ("r1", "B|c", 98.5, 248.0)])
    ani_p = write_ani(tmp_path / "ani.csv")
    msgs: list[str] = []
    run_classification(
        ClassificationRequest(alignment=blast, ani=ani_p, output=tmp_path / "o.csv", config=cfg),
        log=msgs.append,
    )
    assert not any("absent from the ANI matrix" in m for m in msgs)


# 5. blank representative ----------------------------------------------------
def test_blank_representative_falls_back_to_accession(tmp_path):
    p = tmp_path / "m.tsv"
    p.write_text(
        "accession\tspecies\tgenus\tfamily\trepresentative\n"
        "A\tsp1\tg\tF\tA\nB\tsp1\tg\tF\t\nC\tsp2\tg\tF\tC\n"
    )
    md = GenomeMetadata.from_file(p)
    assert md.build_representative_mapping()["B"] == "B"


def test_none_representative_keeps_competitor(ani, cfg):
    df = mkdf([("r1", "A|c", 99.0, 250.0), ("r1", "B|c", 98.5, 248.0)])
    out = VectorizedClassifier(
        ani, config=cfg, representative_mapping={"A": "A", "B": None, "C": "C"}
    ).classify_file(df)
    assert out["placement_uncertainty"][0] == pytest.approx(10.0)


# 6. NaN / NA ANI --------------------------------------------------------------
@pytest.mark.parametrize("token", ["nan", "NA", ""])
def test_nan_ani_treated_as_missing(tmp_path, token, caplog):
    p = tmp_path / "ani.csv"
    p.write_text(f"genome,A,B\nA,100,{token}\nB,{token},100\n")
    with caplog.at_level(logging.WARNING):
        m = ANIMatrix.from_file(p)
    assert m.get_ani("A", "B") == pytest.approx(70.0)  # default for missing
    assert "missing" in caplog.text


def test_nonfinite_novelty_is_unclassified(cfg):
    from metadarkmatter.core.classification.bayesian import BayesianClassifier

    df = pl.DataFrame({
        "novelty_index": [float("nan"), 1.0],
        "placement_uncertainty": [0.0, 0.0],
        "num_ambiguous_hits": [1, 1],
    })
    out = BayesianClassifier(cfg).classify_primary(df)
    assert out["taxonomic_call"][0] == "Unclassified"
    assert out["taxonomic_call"][1] == "Known Species"


# 7. single-genome matrix ------------------------------------------------------
def test_single_genome_matrix_competitor_outside_matrix(cfg):
    a1 = ANIMatrix({"A": {}})
    out = VectorizedClassifier(a1, config=cfg).classify_file(
        mkdf([("r1", "A|c", 98.0, 250.0), ("r1", "B|d", 97.5, 249.0)])
    )
    assert out["placement_uncertainty"][0] == pytest.approx(100.0)


# 8. metadata ------------------------------------------------------------------
def test_duplicate_metadata_rows_do_not_duplicate_reads(caplog):
    md = pl.DataFrame({"accession": ["A", "A"], "species": ["s1", "s2"], "genus": ["g", "g"]})
    with caplog.at_level(logging.WARNING):
        gm = GenomeMetadata(md)
    reads = pl.DataFrame({"read_id": ["r1"], "best_match_genome": ["A"]})
    joined = gm.join_classifications(reads)
    assert joined.height == 1 and joined["species"][0] == "s1"
    assert "duplicate" in caplog.text


def test_infer_target_family_tiebreak_deterministic():
    md = pl.DataFrame({
        "accession": ["1", "2", "3", "4"], "species": ["a"] * 4, "genus": ["g"] * 4,
        "family": ["Zeta", "Alpha", "Zeta", "Alpha"],
    })
    assert GenomeMetadata(md).infer_target_family() == "Alpha"
    assert GenomeMetadata(md.reverse()).infer_target_family() == "Alpha"


# 9. base vs vectorized agree on the "second" genome ---------------------------
def test_base_classifier_uses_weighted_order_for_second_genome():
    from metadarkmatter.core.classification import ANIWeightedClassifier
    from metadarkmatter.models.blast import BlastHit, BlastResult

    ani = ANIMatrix({"A": {"B": 99.0, "C": 80.0}, "B": {"A": 99.0, "C": 85.0},
                     "C": {"A": 80.0, "B": 85.0}})
    cfg = ScoringConfig(min_alignment_length=0, min_alignment_fraction=0.0)
    # A is best (full coverage). B has a slightly higher raw bitscore than C
    # but poor coverage, so after weighting C is the "second" genome.
    rows = [
        ("r1", "A|c", 99.0, 200, 1, 150, 200.0),
        ("r1", "B|c", 98.0, 60, 1, 60, 199.0),
        ("r1", "C|c", 98.0, 150, 1, 150, 198.0),
    ]
    hits = tuple(
        BlastHit(qseqid=q, sseqid=s, pident=p, length=ln, mismatch=0, gapopen=0,
                 qstart=qs, qend=qe, sstart=1, send=qe, evalue=1e-30, bitscore=bs, qlen=150)
        for q, s, p, ln, qs, qe, bs in rows
    )
    base = ANIWeightedClassifier(ani, config=cfg).classify_read(BlastResult(read_id="r1", hits=hits))
    df = pl.DataFrame(
        [(q, s, p, ln, 0, 0, qs, qe, 1, qe, 1e-30, bs, 150) for q, s, p, ln, qs, qe, bs in rows],
        schema=COLS, orient="row",
    ).with_columns(pl.col("sseqid").str.extract(r"^([^|]+)\|", 1).alias("genome_name"))
    vec = VectorizedClassifier(ani, config=cfg).classify_file(df)
    assert base.placement_uncertainty == pytest.approx(vec["placement_uncertainty"][0])


# 10. evaluation with null truth -------------------------------------------------
def test_evaluation_null_truth_rows_dropped():
    df = pl.DataFrame({
        "taxonomic_call": ["Known Species", "Novel Genus", "Known Species"],
        "true_category": ["Known Species", None, "Novel Species"],
    })
    res = evaluate_classifications(df, truth_column="true_category",
                                   prediction_column="taxonomic_call", confidence_column=None)
    assert res.n_rows == 2


# E. YAML unknown keys -------------------------------------------------------------
def test_yaml_unknown_keys_warn(tmp_path, caplog):
    p = tmp_path / "c.yaml"
    p.write_text("thresholds:\n  novelty:\n    known_maxx: 3\n    known_max: 4\nbogus: 1\n")
    with caplog.at_level(logging.WARNING):
        cfg = ScoringConfig.from_yaml(p)
    assert cfg.novelty_known_max == 4
    assert "thresholds.novelty.known_maxx" in caplog.text and "bogus" in caplog.text


def test_yaml_known_keys_do_not_warn(tmp_path, caplog):
    p = tmp_path / "c.yaml"
    ScoringConfig().to_yaml(p)
    with caplog.at_level(logging.WARNING):
        ScoringConfig.from_yaml(p)
    assert "unknown" not in caplog.text.lower()


# CLI: batch / sensitivity / extract-novel -------------------------------------------
runner = CliRunner()


def test_batch_exits_nonzero_when_a_file_fails(tmp_path):
    d = tmp_path / "in"
    d.mkdir()
    write_blast(d / "good.tsv", [("r1", "A|c", 99.0, 250.0)])
    (d / "bad.tsv").write_text("not\ta\tblast\n")
    ani_p = write_ani(tmp_path / "ani.csv")
    res = runner.invoke(app, ["score", "batch", "--alignment-dir", str(d), "--ani", str(ani_p),
                              "--output-dir", str(tmp_path / "o"), "--pattern", "*.tsv",
                              "--min-alignment-length", "0", "--min-alignment-fraction", "0"])
    assert res.exit_code != 0
    assert (tmp_path / "o" / "good_classifications.csv").exists()


def _class_csv(path: Path) -> Path:
    pl.DataFrame({
        "read_id": ["r1", "r2", "r3", "r4"],
        "best_match_genome": ["A", "A", "B", "B"],
        "taxonomic_call": ["Novel Species", "Novel Species", "Novel Genus", "Novel Species"],
        "novelty_index": [8.0, 12.0, 22.0, 9.0],
        "placement_uncertainty": [0.5, 3.0, 0.5, 0.5],
        "top_hit_identity": [92.0, 88.0, 78.0, 91.0],
        "num_ambiguous_hits": [1, 2, 1, 1],
        "identity_gap": [None, 1.0, None, None],
    }).write_csv(path)
    return path


def test_extract_novel_filters_intersect_and_empty_outputs(tmp_path):
    c = _class_csv(tmp_path / "c.csv")
    o, ids = tmp_path / "o.csv", tmp_path / "ids.txt"
    res = runner.invoke(app, ["score", "extract-novel", "-c", str(c), "-o", str(o), "-r", str(ids),
                              "-t", "species", "--min-novelty", "5", "--min-reads", "2"])
    assert res.exit_code == 0
    # r2 excluded by max-uncertainty (default 2.0); r1 alone in A (<2 reads), B has r4 only.
    assert ids.read_text().strip() == ""
    assert o.exists() and o.read_text().startswith("best_match_genome")

    # category + min-novelty intersect: genus reads are never returned for species
    res = runner.invoke(app, ["score", "extract-novel", "-c", str(c), "-o", str(o), "-r", str(ids),
                              "-t", "species", "--min-novelty", "5", "--no-group"])
    assert sorted(ids.read_text().split()) == ["r1", "r4"]

    # nothing matches -> header-only output still written
    o2, ids2 = tmp_path / "o2.csv", tmp_path / "ids2.txt"
    res = runner.invoke(app, ["score", "extract-novel", "-c", str(c), "-o", str(o2),
                              "-r", str(ids2), "--min-novelty", "50"])
    assert res.exit_code == 0 and o2.exists() and ids2.exists()


def test_sensitivity_labels_legacy_call_and_sums(tmp_path):
    c = tmp_path / "c.csv"
    pl.DataFrame({
        "read_id": ["r1", "r2", "r3"],
        "taxonomic_call": ["Known Species", "Off-target", "Novel Species"],
        "novelty_index": [1.0, 5.0, 10.0],
        "placement_uncertainty": [0.0, None, 0.5],
        "num_ambiguous_hits": [1, 1, 1],
        "identity_gap": [None, None, 5.0],
    }).write_csv(c)
    o = tmp_path / "s.tsv"
    res = runner.invoke(app, ["score", "sensitivity", "-c", str(c), "-o", str(o),
                              "--steps", "3", "-f", "TSV"])
    assert res.exit_code == 0, res.output
    df = pl.read_csv(o, separator="\t")
    assert "legacy_call" in df.columns and "taxonomic_call" not in df.columns
    assert (df.group_by("novelty_known_max").agg(pl.col("count").sum())["count"] == 3).all()
    assert "Off-target" in df["legacy_call"].to_list()
