"""Regression tests for classifier / summary bugs found in the 2026-11 docs audit."""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from metadarkmatter.cli.score import _generate_summary
from metadarkmatter.core.aai_matrix_builder import AAIMatrix
from metadarkmatter.core.classification.ani_matrix import ANIMatrix
from metadarkmatter.core.classification.bayesian import BayesianClassifier
from metadarkmatter.core.classification.classifiers.base import ANIWeightedClassifier
from metadarkmatter.core.classification.classifiers.vectorized import (
    VectorizedClassifier,
)
from metadarkmatter.core.classification.runner import (
    ClassificationRequest,
    run_classification,
)
from metadarkmatter.core.parsers import StreamingBlastParser
from metadarkmatter.models.config import ScoringConfig


def _write(rows: list[tuple[str, str, float, float]], path: Path) -> Path:
    """rows: (qseqid, sseqid, pident, bitscore); 150 bp full-length hits."""
    path.write_text(
        "".join(
            f"{q}\t{s}\t{p}\t150\t0\t0\t1\t150\t1\t150\t1e-50\t{b}\t150\n"
            for q, s, p, b in rows
        )
    )
    return path


def _ani(pairs: dict[tuple[str, str], float]) -> ANIMatrix:
    genomes = sorted({g for pair in pairs for g in pair})
    d = {g: {h: (100.0 if g == h else 0.0) for h in genomes} for g in genomes}
    for (a, b), v in pairs.items():
        d[a][b] = d[b][a] = v
    return ANIMatrix(d)


class TestAmbiguityScopeCountsGenomes:
    def test_multiple_contigs_of_one_competitor(self, tmp_path):
        """Two contig hits to one same-species competitor are still within_species.

        same_species_count used to count hit rows (2) while
        num_secondary_genomes counted genomes (1), so the read fell through to
        across_genera.
        """
        ani = _ani({("A", "B"): 98.0})
        rows = [
            ("r1", "A|c1", 99.0, 280.0),
            ("r1", "B|c1", 96.0, 275.0),
            ("r1", "B|c2", 96.0, 274.0),
        ]
        out = VectorizedClassifier(ani, config=ScoringConfig()).classify_file(
            _write(rows, tmp_path / "a.tsv")
        )
        assert out["ambiguity_scope"][0] == "within_species"


class TestAAIGenusOverride:
    @pytest.mark.parametrize(
        ("aai", "expected"),
        [
            (70.0, "Novel Species"),
            (60.0, "Ambiguous Within Genus"),
            (50.0, "Novel Genus"),
        ],
    )
    def test_vectorized_applies_aai_and_matches_base(self, tmp_path, aai, expected):
        """--aai used to be loaded but never applied by the CLI classifier."""
        ani = _ani({("A", "B"): 99.0})
        aai_m = AAIMatrix({"A": {"A": 100.0, "B": aai}, "B": {"A": aai, "B": 100.0}})
        path = _write(
            [("r1", "A|c1", 78.0, 200.0), ("r1", "B|c1", 77.5, 197.0)],
            tmp_path / "a.tsv",
        )
        v = VectorizedClassifier(ani, aai_matrix=aai_m, config=ScoringConfig()).classify_file(
            path
        )
        b = ANIWeightedClassifier(ani, aai_matrix=aai_m, config=ScoringConfig()).classify_read(
            next(StreamingBlastParser(path).iter_reads())
        )
        assert v["taxonomic_call"][0] == expected
        assert v["legacy_call"][0] == expected
        assert b.taxonomic_call.value == expected

    def test_no_override_outside_genus_range(self, tmp_path):
        """A Known Species read is not touched by the AAI rule."""
        ani = _ani({("A", "B"): 99.0})
        aai_m = AAIMatrix({"A": {"A": 100.0, "B": 50.0}, "B": {"A": 50.0, "B": 100.0}})
        path = _write(
            [("r1", "A|c1", 99.5, 280.0), ("r1", "B|c1", 99.0, 278.0)],
            tmp_path / "a.tsv",
        )
        with_aai = VectorizedClassifier(
            ani, aai_matrix=aai_m, config=ScoringConfig()
        ).classify_file(path)
        without = VectorizedClassifier(ani, config=ScoringConfig()).classify_file(path)
        assert with_aai["taxonomic_call"][0] == without["taxonomic_call"][0]


class TestBaseClassifierParity:
    def test_hits_to_same_genome_are_confident(self, tmp_path):
        """All ambiguous hits to the best genome: U=0, not U=100."""
        ani = _ani({("GA", "GB"): 85.0})
        path = _write(
            [("r1", "GA|c1", 99.0, 250.0), ("r1", "GA|c2", 98.5, 248.0)],
            tmp_path / "a.tsv",
        )
        b = ANIWeightedClassifier(ani, config=ScoringConfig()).classify_read(
            next(StreamingBlastParser(path).iter_reads())
        )
        v = VectorizedClassifier(ani, config=ScoringConfig()).classify_file(path)
        assert b.placement_uncertainty == 0.0
        assert v["placement_uncertainty"][0] == 0.0
        assert b.taxonomic_call.value == v["taxonomic_call"][0] == "Known Species"

    def test_scalar_single_hit_boost_uses_protein_threshold(self):
        """Protein single-hit boost starts at N >= 10, not the nucleotide 4."""
        nuc = BayesianClassifier(ScoringConfig())
        prot = BayesianClassifier(ScoringConfig(alignment_mode="protein"))
        prot_cfg = ScoringConfig(alignment_mode="protein")
        n = 7.0  # between the nucleotide (4) and protein (10) thresholds
        p_one = prot.compute_posteriors(n, 0.0, num_hits=1).p_ambiguous
        p_two = prot.compute_posteriors(n, 0.0, num_hits=2, identity_gap=10.0).p_ambiguous
        assert p_one == pytest.approx(p_two)
        assert prot_cfg.get_effective_thresholds()["novelty_novel_species_min"] > n
        # nucleotide mode still boosts at the same N
        assert nuc.compute_posteriors(n, 0.0, num_hits=1).p_ambiguous > nuc.compute_posteriors(
            n, 0.0, num_hits=2, identity_gap=10.0
        ).p_ambiguous


class TestSummary:
    def test_diversity_and_species_counts(self, tmp_path):
        """--summary used to report 0 diversity reads and no species counts."""
        rows =[(f"k{i}", "GCF_000000001.1|c", 99.5, 280.0) for i in range(4)]
        rows += [(f"n{i}", "GCF_000000002.1|c", 90.0, 200.0) for i in range(3)]
        aln = _write(rows, tmp_path / "a.tsv")
        ani_csv = tmp_path / "ani.csv"
        ani_csv.write_text(
            "genome,GCF_000000001.1,GCF_000000002.1\n"
            "GCF_000000001.1,100,85\nGCF_000000002.1,85,100\n"
        )
        meta = tmp_path / "meta.tsv"
        meta.write_text(
            "accession\tspecies\tgenus\tfamily\n"
            "GCF_000000001.1\tSp one\tGen\tFam\n"
            "GCF_000000002.1\tSp two\tGen\tFam\n"
        )
        result = run_classification(
            ClassificationRequest(
                alignment=aln, ani=ani_csv, output=tmp_path / "o.csv",
                metadata=meta, config=ScoringConfig(),
            ),
            log=lambda *a, **k: None,
        )
        assert result.classification_df is not None
        assert "species" in result.classification_df.columns
        s = _generate_summary(result.classification_df)
        assert s.diversity_known + s.diversity_novel + s.diversity_uncertain + s.off_target == 7
        assert s.diversity_known == s.known_species
        assert s.diversity_novel == s.novel_species + s.novel_genus
        assert sum(s.species_hit_counts.values()) == 7

    def test_summary_without_diversity_column(self):
        df = pl.DataFrame({
            "taxonomic_call": ["Known Species"],
            "best_match_genome": ["G"],
            "novelty_index": [1.0],
            "placement_uncertainty": [0.0],
        })
        assert _generate_summary(df).diversity_known == 0
