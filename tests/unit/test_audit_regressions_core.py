"""Regression tests for classifier / parser defects found in the 2026-11 audit.

Each test pins one concrete failure that was reproduced against real or
synthetic data before the fix.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import polars as pl
import pytest

from metadarkmatter.clients.gtdb import GTDBClient
from metadarkmatter.core.classification.adaptive import (
    build_adaptive_config,
    detect_species_boundary,
)
from metadarkmatter.core.classification.ani_matrix import ANIMatrix
from metadarkmatter.core.classification.classifiers.vectorized import (
    VectorizedClassifier,
)
from metadarkmatter.core.classification.qc import QCMetrics, compute_post_qc
from metadarkmatter.core.classification.sensitivity import run_sensitivity_analysis
from metadarkmatter.core.genome_utils import concatenate_genomes_with_mapping
from metadarkmatter.core.id_mapping import ContigIdMapping
from metadarkmatter.core.parsers import StreamingBlastParser
from metadarkmatter.core.sequence_validation import validate_fasta
from metadarkmatter.models.config import ScoringConfig

COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen",
]


def _ani(pairs: dict[tuple[str, str], float], genomes: list[str]) -> ANIMatrix:
    d = {g: {h: (100.0 if g == h else 0.0) for h in genomes} for g in genomes}
    for (a, b), v in pairs.items():
        d[a][b] = v
        d[b][a] = v
    return ANIMatrix(d)


def _hit(q, g, pid, bits, qs=1, qe=150, qlen=150, contig="c1"):
    return {
        "qseqid": q, "sseqid": f"{g}|{contig}", "pident": float(pid),
        "length": abs(qe - qs) + 1, "mismatch": 0, "gapopen": 0,
        "qstart": qs, "qend": qe, "sstart": 1, "send": 150,
        "evalue": 1e-50, "bitscore": float(bits), "qlen": qlen,
    }


def _write(rows, path: Path, header: str | None = None) -> Path:
    lines = []
    if header:
        lines.append(header)
    for r in rows:
        lines.append("\t".join(str(r[c]) for c in COLS))
    path.write_text("\n".join(lines) + "\n")
    return path


GENOMES = ["GA", "GB", "GC", "GD"]
ANI = _ani(
    {("GA", "GB"): 95.21, ("GA", "GC"): 92.76, ("GB", "GC"): 94.26,
     ("GA", "GD"): 85.0, ("GB", "GD"): 85.0, ("GC", "GD"): 85.0},
    GENOMES,
)


class TestReverseStrandCoordinates:
    def test_mmseqs_minus_strand_best_hit_is_kept(self, tmp_path):
        """MMseqs2 reports minus-strand hits as qstart=150, qend=1.

        Before the fix the coverage term was negative, the perfect hit to
        the true genome was filtered/down-weighted and a 90% hit to another
        genome won (reads from a reference genome became 'Novel Species').
        """
        rows = [
            _hit("r1", "GA", 99.3, 266, qs=150, qe=1),
            _hit("r1", "GD", 90.9, 169, qs=1, qe=121),
        ]
        out = VectorizedClassifier(ani_matrix=ANI, config=ScoringConfig()).classify_file(
            _write(rows, tmp_path / "a.tsv")
        )
        row = out.row(0, named=True)
        assert row["best_match_genome"] == "GA"
        assert row["top_hit_identity"] == pytest.approx(99.3)
        assert row["taxonomic_call"] == "Known Species"


class TestDeterministicSecondGenome:
    def test_tied_secondaries_independent_of_batch(self, tmp_path):
        """With tied bitscores the 'second' genome must not depend on batch.

        Default uncertainty_mode='second' picked the competitor via an
        unstable sort on raw bitscore, so the same read got different
        placement_uncertainty in streaming vs in-memory runs.
        """
        target = [
            _hit("t", "GA", 82.6, 154),
            _hit("t", "GB", 82.6, 154),
            _hit("t", "GC", 82.6, 154),
        ]
        filler = [
            _hit(f"f{i:04d}", g, 99.0, 250) for i in range(300) for g in ("GA", "GD")
        ]
        clf = VectorizedClassifier(ani_matrix=ANI, config=ScoringConfig())
        alone = clf.classify_file(_write(target, tmp_path / "a.tsv"))
        mixed = clf.classify_file(_write(filler + target + filler, tmp_path / "b.tsv"))
        u_alone = alone.filter(pl.col("read_id") == "t")["placement_uncertainty"][0]
        u_mixed = mixed.filter(pl.col("read_id") == "t")["placement_uncertainty"][0]
        assert u_alone == u_mixed
        # best = GA (genome-name tiebreak); second = GB (same tiebreak)
        assert u_alone == pytest.approx(100.0 - 95.21)


class TestStreaming:
    def _rows(self):
        rows = []
        for i in range(60):
            # integer-looking bitscores first, decimals later
            rows.append(_hit(f"r{i:03d}", "GA", 99.0, 250))
            rows.append(_hit(f"r{i:03d}", "GB", 97.5, 45.4 if i > 40 else 240))
        return rows

    def test_streaming_matches_in_memory(self, tmp_path):
        path = _write(self._rows(), tmp_path / "a.tsv", header="# comment line")
        clf = VectorizedClassifier(ani_matrix=ANI, config=ScoringConfig())
        mem = clf.classify_file(path).sort("read_id")
        out = tmp_path / "s.csv"
        clf.stream_to_file(path, out, output_format="csv", partition_size=7)
        streamed = pl.read_csv(out).sort("read_id")
        assert streamed.columns == mem.columns
        assert streamed["taxonomic_call"].to_list() == mem["taxonomic_call"].to_list()
        assert streamed["placement_uncertainty"].to_list() == pytest.approx(
            mem["placement_uncertainty"].to_list()
        )

    def test_streaming_parquet_schema_stable(self, tmp_path):
        """Batches with and without competing hits used to differ in dtype."""
        rows = [_hit(f"s{i:03d}", "GA", 99.0, 250) for i in range(20)]
        rows += self._rows()
        path = _write(rows, tmp_path / "a.tsv")
        out = tmp_path / "s.parquet"
        VectorizedClassifier(ani_matrix=ANI, config=ScoringConfig()).stream_to_file(
            path, out, output_format="parquet", partition_size=10
        )
        df = pl.read_parquet(out)
        assert df.height == 80
        assert df.schema["placement_uncertainty"] == pl.Float64

    def test_streaming_family_validation_columns_aligned(self, tmp_path):
        """An all-Off-target batch must use the same columns as the others."""
        rows = [_hit(f"k{i}", "GA", 99.0, 250) for i in range(5)]
        rows += [_hit(f"x{i}", "EXT", 99.0, 250) for i in range(5)]
        path = _write(rows, tmp_path / "a.tsv")
        cfg = ScoringConfig(target_family="f__Test")
        out = tmp_path / "s.csv"
        VectorizedClassifier(ani_matrix=ANI, config=cfg).stream_to_file(
            path, out, output_format="csv", partition_size=5
        )
        df = pl.read_csv(out)
        off = df.filter(pl.col("read_id").str.starts_with("x"))
        assert off.height == 5
        assert set(off["taxonomic_call"]) == {"Off-target"}
        assert set(off["legacy_call"]) == {"Off-target"}
        assert off["p_known_species"].is_null().all()

    def test_noncontiguous_reads_rejected(self, tmp_path):
        rows = [_hit("a", "GA", 99, 250)] + [
            _hit(f"m{i}", "GA", 99, 250) for i in range(10)
        ] + [_hit("a", "GB", 98, 240)]
        path = _write(rows, tmp_path / "a.tsv")
        with pytest.raises(ValueError, match="not contiguous"):
            VectorizedClassifier(ani_matrix=ANI, config=ScoringConfig()).stream_to_file(
                path, tmp_path / "o.csv", output_format="csv", partition_size=4
            )


class TestANIMatrixInput:
    def test_fraction_scale_rejected(self):
        with pytest.raises(ValueError, match="fractions"):
            ANIMatrix({"A": {"A": 1.0, "B": 0.99}, "B": {"A": 0.99, "B": 1.0}})

    def test_one_directional_pair_mirrored(self):
        m = ANIMatrix({"A": {"A": 100.0}, "B": {"A": 97.0, "B": 100.0}},
                      symmetry_check="off")
        assert m.get_ani("A", "B") == pytest.approx(97.0)
        assert m.get_ani("B", "A") == pytest.approx(97.0)

    def test_uncertainty_exact_for_two_decimal_ani(self, tmp_path):
        ani = _ani({("GA", "GB"): 98.6}, ["GA", "GB"])
        rows = [_hit("r", "GA", 99.0, 250), _hit("r", "GB", 99.0, 249)]
        out = VectorizedClassifier(ani_matrix=ani, config=ScoringConfig()).classify_file(
            _write(rows, tmp_path / "a.tsv")
        )
        assert out["placement_uncertainty"][0] == pytest.approx(1.4, abs=1e-9)


class TestAdaptiveAndSensitivity:
    def test_missing_pairs_do_not_clamp_boundary(self):
        genomes = [f"S{s}_{k}" for s in range(6) for k in range(4)]
        d = {g: {h: (100.0 if g == h else 0.0) for h in genomes} for g in genomes}
        for i, a in enumerate(genomes):
            for b in genomes[i + 1:]:
                same = a.split("_")[0] == b.split("_")[0]
                if same:
                    v = 99.0
                elif (hash(a + b) % 2) == 0:
                    continue  # uncomputed distant pair
                else:
                    v = 86.0
                d[a][b] = d[b][a] = v
        res = detect_species_boundary(ANIMatrix(d))
        assert res.species_boundary > 90.0

    def test_build_adaptive_config_preserves_settings(self):
        base = ScoringConfig(
            alignment_mode="protein", target_family="f__X",
            min_alignment_length=50, include_legacy_scores=True,
        )
        res = detect_species_boundary(ANI, min_genomes=100)  # fallback
        cfg = build_adaptive_config(base, res)
        assert cfg.alignment_mode == "protein"
        assert cfg.target_family == "f__X"
        assert cfg.min_alignment_length == 50
        assert cfg.include_legacy_scores is True

    def test_sensitivity_on_classify_output(self, tmp_path):
        rows = [_hit(f"r{i}", "GA", 99 - i * 0.5, 250) for i in range(20)]
        out = VectorizedClassifier(ani_matrix=ANI, config=ScoringConfig()).classify_file(
            _write(rows, tmp_path / "a.tsv")
        )
        assert "num_secondary_genomes" not in out.columns
        res = run_sensitivity_analysis(out, ScoringConfig(), steps=3)
        assert sum(v[0] for v in res.counts.values()) == 20

    def test_qc_low_confidence_uses_percent_scale(self):
        df = pl.DataFrame({
            "taxonomic_call": ["Known Species"] * 2,
            "confidence_score": [30.0, 90.0],
        })
        assert compute_post_qc(QCMetrics(), df).low_confidence_fraction == 0.5


class TestParsersAndIO:
    def test_chunk_boundary_read_not_duplicated(self, tmp_path):
        rows = [_hit("zz", "GA", 99, 250)] * 1500 + [_hit("aa", "GA", 99, 250)] * 100
        path = _write(rows, tmp_path / "a.tsv")
        p = StreamingBlastParser(path, chunk_size=1000)
        assert sorted((r.read_id, r.num_hits) for r in p.iter_reads_chunked()) == [
            ("aa", 100), ("zz", 1500),
        ]
        assert sorted((r.read_id, len(r.hits)) for r in p.iter_reads_fast()) == [
            ("aa", 100), ("zz", 1500),
        ]

    def test_parse_lazy_normalizes_reverse_coordinates(self, tmp_path):
        path = _write([_hit("r", "GA", 99, 250, qs=150, qe=1)], tmp_path / "a.tsv")
        df = StreamingBlastParser(path).parse_eager()
        assert (df["qstart"][0], df["qend"][0]) == (1, 150)

    def test_id_mapping_ncbi_style_and_numeric_ids(self, tmp_path):
        tsv = tmp_path / "m.tsv"
        tsv.write_text(
            "original_contig_id\tgenome_accession\n"
            "NC_006570.2\tGCF_000008985.1\n007\tGCF_1.1\n"
        )
        m = ContigIdMapping.from_tsv(tsv)
        assert m.transform_sseqid("ref|NC_006570.2|") == "GCF_000008985.1|NC_006570.2"
        df = pl.DataFrame({"sseqid": ["ref|NC_006570.2|", "007", "GCF_9.1|x", "zz"]})
        assert m.transform_column(df)["sseqid"].to_list() == [
            "GCF_000008985.1|NC_006570.2", "GCF_1.1|007", "GCF_9.1|x", "zz",
        ]

    def test_concatenate_without_trailing_newline(self, tmp_path):
        gdir = tmp_path / "g"
        gdir.mkdir()
        (gdir / "GCF_000000001.1.fna").write_text(">c1\nACGT")  # no newline
        with gzip.open(gdir / "GCF_000000002.1.fna.gz", "wt") as fh:
            fh.write(">c2\nTTTT\n")
        out = tmp_path / "pan.fna"
        concatenate_genomes_with_mapping(gdir, out, tmp_path / "map.tsv", pattern="*.fna*")
        lines = out.read_text().splitlines()
        assert ">GCF_000000002.1|c2" in lines
        assert "ACGT" in lines

    def test_gtdb_placeholder_taxa_accepted(self):
        assert GTDBClient.validate_taxon_format("f__CAG-138")
        assert GTDBClient.validate_taxon_format("g__UBA1234")
        assert not GTDBClient.validate_taxon_format("Francisellaceae")

    def test_validate_fasta_max_records(self, tmp_path):
        fa = tmp_path / "x.fa"
        fa.write_text(">a\nACGT\n>b\nACGT\n>c\nACGT\n")
        res = validate_fasta(fa, max_records=2)
        assert res.record_count == 2
        assert not res.issues
