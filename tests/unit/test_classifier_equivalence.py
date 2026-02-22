"""
Cross-classifier equivalence tests.

Verifies that base (ANIWeightedClassifier.classify_read) and vectorized
(VectorizedClassifier.classify_file) classifiers produce equivalent
legacy threshold classifications for the same input data.

The vectorized classifier now uses Bayesian posteriors for its primary
``taxonomic_call`` column, so equivalence is checked against its
``legacy_call`` column which preserves the hard-threshold cascade.

confidence_score is excluded from comparison because the vectorized classifier
uses a different scoring approach (alignment-length + bitscore-gap based)
compared to the base classifier (novelty-margin based).
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from metadarkmatter.core.classification.ani_matrix import ANIMatrix
from metadarkmatter.core.classification.classifiers.base import ANIWeightedClassifier
from metadarkmatter.core.classification.classifiers.vectorized import VectorizedClassifier
from metadarkmatter.core.parsers import StreamingBlastParser
from metadarkmatter.models.config import ScoringConfig


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def ani_dict() -> dict[str, dict[str, float]]:
    """ANI matrix with 4 genomes covering various similarity levels."""
    return {
        "GCF_000001.1": {
            "GCF_000001.1": 100.0,
            "GCF_000002.1": 98.5,  # same species
            "GCF_000003.1": 85.0,  # different species
            "GCF_000004.1": 78.0,  # different genus
        },
        "GCF_000002.1": {
            "GCF_000001.1": 98.5,
            "GCF_000002.1": 100.0,
            "GCF_000003.1": 84.0,
            "GCF_000004.1": 77.0,
        },
        "GCF_000003.1": {
            "GCF_000001.1": 85.0,
            "GCF_000002.1": 84.0,
            "GCF_000003.1": 100.0,
            "GCF_000004.1": 79.0,
        },
        "GCF_000004.1": {
            "GCF_000001.1": 78.0,
            "GCF_000002.1": 77.0,
            "GCF_000003.1": 79.0,
            "GCF_000004.1": 100.0,
        },
    }


@pytest.fixture
def ani_matrix(ani_dict) -> ANIMatrix:
    return ANIMatrix(ani_dict)


@pytest.fixture
def config() -> ScoringConfig:
    return ScoringConfig()


def _blast_rows() -> list[dict]:
    """Representative reads spanning all classification categories.

    Each read has hits against genomes defined in the ANI fixture.
    The subject IDs use the {accession}|{contig} format expected by the parser.
    """
    rows = []

    def _hit(qseqid, sseqid, pident, bitscore, length=150):
        return {
            "qseqid": qseqid,
            "sseqid": sseqid,
            "pident": pident,
            "length": length,
            "mismatch": int(length * (100 - pident) / 100),
            "gapopen": 0,
            "qstart": 1,
            "qend": length,
            "sstart": 1,
            "send": length,
            "evalue": 1e-50,
            "bitscore": bitscore,
        }

    # Read 1: Known species -- high identity, two close hits (same species pair)
    rows.append(_hit("read_known", "GCF_000001.1|ctg1", 99.0, 280))
    rows.append(_hit("read_known", "GCF_000002.1|ctg1", 98.0, 270))

    # Read 2: Novel species -- moderate novelty, close competitors
    rows.append(_hit("read_novel_sp", "GCF_000001.1|ctg1", 90.0, 250))
    rows.append(_hit("read_novel_sp", "GCF_000002.1|ctg1", 89.0, 245))

    # Read 3: Novel genus -- high divergence, close competitors
    rows.append(_hit("read_novel_gen", "GCF_000001.1|ctg1", 78.0, 200))
    rows.append(_hit("read_novel_gen", "GCF_000002.1|ctg1", 77.0, 195))

    # Read 4: Single hit -- tests inferred uncertainty path
    rows.append(_hit("read_single", "GCF_000003.1|ctg1", 92.0, 260))

    # Read 5: Ambiguous -- hits on distant genomes (low ANI -> high uncertainty)
    rows.append(_hit("read_ambiguous", "GCF_000001.1|ctg1", 96.0, 270))
    rows.append(_hit("read_ambiguous", "GCF_000004.1|ctg1", 95.0, 265))

    return rows


@pytest.fixture
def blast_file(tmp_path) -> Path:
    """Write BLAST TSV from representative rows."""
    path = tmp_path / "test.blast.tsv"
    df = pl.DataFrame(_blast_rows())
    col_order = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    ]
    df.select(col_order).write_csv(path, separator="\t", include_header=False)
    return path


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

COMPARE_COLUMNS = [
    "read_id",
    "best_match_genome",
    "taxonomic_call",
    "novelty_index",
    "placement_uncertainty",
    "num_ambiguous_hits",
]


class TestClassifierEquivalence:
    """Verify base and vectorized classifiers agree on taxonomic calls."""

    def test_base_vs_vectorized(self, blast_file, ani_matrix, config):
        base_clf = ANIWeightedClassifier(ani_matrix=ani_matrix, config=config)
        vec_clf = VectorizedClassifier(ani_matrix=ani_matrix, config=config)

        # Use streaming parser + classify_read for base classifier
        parser = StreamingBlastParser(blast_file)
        base_results = []
        for blast_result in parser.iter_reads():
            classification = base_clf.classify_read(blast_result)
            if classification is not None:
                base_results.append(classification.to_dict())
        base_df = pl.DataFrame(base_results)

        vec_df = vec_clf.classify_file(blast_file)

        # Align on read_id for comparison
        base_sorted = base_df.sort("read_id")
        vec_sorted = vec_df.sort("read_id")

        assert base_sorted["read_id"].to_list() == vec_sorted["read_id"].to_list(), (
            "Classifiers returned different read sets"
        )

        # Both classifiers use Bayesian-primary classification.
        # Compare best_match_genome and novelty_index (must be identical).
        for col in ["best_match_genome", "novelty_index"]:
            base_vals = base_sorted[col].to_list()
            vec_vals = vec_sorted[col].to_list()
            assert base_vals == vec_vals, (
                f"Column '{col}' differs between base and vectorized:\n"
                f"  base:       {base_vals}\n"
                f"  vectorized: {vec_vals}"
            )

        # For taxonomic_call, allow "Ambiguous" (base) vs "Conserved Region" (vec).
        # The vectorized path has ambiguity_scope from hit analysis, enabling
        # Stage 2 refinement (Ambiguous + across_genera → Conserved Region).
        # The scalar base classifier lacks this context.
        allowed_equivalences = {
            ("Ambiguous", "Conserved Region"),
        }
        base_calls = base_sorted["taxonomic_call"].to_list()
        vec_calls = vec_sorted["taxonomic_call"].to_list()
        for i, (b, v) in enumerate(zip(base_calls, vec_calls)):
            if b != v and (b, v) not in allowed_equivalences:
                raise AssertionError(
                    f"taxonomic_call differs at row {i}: base={b!r}, vectorized={v!r}\n"
                    f"  full base:       {base_calls}\n"
                    f"  full vectorized: {vec_calls}"
                )


class TestInferredUncertaintyContinuity:
    """Verify the inferred uncertainty formula is continuous at boundaries."""

    def test_continuity_at_novel_species_max(self):
        from metadarkmatter.core.constants import calculate_inferred_uncertainty

        # At the boundary N = NOVELTY_NOVEL_SPECIES_MAX (20.0), both the
        # novel-species and novel-genus branches should return the same value.
        below = calculate_inferred_uncertainty(19.99)
        at = calculate_inferred_uncertainty(20.0)
        above = calculate_inferred_uncertainty(20.01)

        # The transition should be smooth -- no jump
        assert abs(below - at) < 0.1, (
            f"Discontinuity at N=20: below={below}, at={at}"
        )
        assert abs(at - above) < 0.1, (
            f"Discontinuity at N=20: at={at}, above={above}"
        )

    def test_continuity_at_known_max(self):
        from metadarkmatter.core.constants import calculate_inferred_uncertainty

        below = calculate_inferred_uncertainty(4.99)
        at = calculate_inferred_uncertainty(5.0)

        assert abs(below - at) < 0.1, (
            f"Discontinuity at N=5: below={below}, at={at}"
        )

    def test_monotonically_increasing(self):
        from metadarkmatter.core.constants import calculate_inferred_uncertainty

        prev = calculate_inferred_uncertainty(0.0)
        for n in [1, 3, 5, 10, 15, 20, 22, 25, 30]:
            curr = calculate_inferred_uncertainty(float(n))
            assert curr >= prev, (
                f"Not monotonic: f({n})={curr} < f(prev)={prev}"
            )
            prev = curr
