"""Regression tests for the round-2 novel-diversity / phylogeny / wrapper fixes."""

from __future__ import annotations

import gzip
import json
import os
import re
import time
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import polars as pl
import pytest
from typer.testing import CliRunner

from metadarkmatter.core.classification.ani_matrix import ANIMatrix
from metadarkmatter.core.metadata import GenomeMetadata
from metadarkmatter.core.novel_diversity import NovelDiversityAnalyzer
from metadarkmatter.core.novel_diversity.models import NovelCluster
from metadarkmatter.core.novel_diversity.neighborhood import (
    PhylogeneticNeighborhoodAnalyzer,
)

REP = "GCF_000000001.1"
NONREP = "GCF_000000009.1"
OTHERS = ["GCF_000000002.1", "GCF_000000003.1", "GCF_000000004.1"]


def _metadata() -> GenomeMetadata:
    rows = [
        (REP, "Genus a", "Genus", "Fam", REP),
        (NONREP, "Genus a", "Genus", "Fam", REP),
        (OTHERS[0], "Genus b", "Genus", "Fam", OTHERS[0]),
        (OTHERS[1], "Genus c", "Genus2", "Fam", OTHERS[1]),
        (OTHERS[2], "Genus d", "Genus2", "Fam", OTHERS[2]),
    ]
    df = pl.DataFrame(
        rows,
        schema=["accession", "species", "genus", "family", "representative"],
        orient="row",
    )
    return GenomeMetadata(df)


def _ani_df() -> pl.DataFrame:
    genomes = [REP, *OTHERS]
    vals = {
        (REP, OTHERS[0]): 85.0,
        (REP, OTHERS[1]): 75.0,
        (REP, OTHERS[2]): 74.0,
        (OTHERS[0], OTHERS[1]): 76.0,
        (OTHERS[0], OTHERS[2]): 75.0,
        (OTHERS[1], OTHERS[2]): 88.0,
    }
    data: dict[str, list] = {"genome": genomes}
    for g2 in genomes:
        col = []
        for g1 in genomes:
            if g1 == g2:
                col.append(100.0)
            else:
                col.append(vals.get((g1, g2), vals.get((g2, g1))))
        data[g2] = col
    return pl.DataFrame(data)


def _novel_df(genome: str = NONREP, n: int = 6, unc: float = 0.0) -> pl.DataFrame:
    return pl.DataFrame({
        "read_id": [f"r{i}" for i in range(n)],
        "best_match_genome": [genome] * n,
        "top_hit_identity": [90.0] * n,
        "novelty_index": [10.0 + 0.1 * i for i in range(n)],
        "placement_uncertainty": [unc] * n,
        "num_ambiguous_hits": [1] * n,
        "taxonomic_call": ["Novel Species"] * n,
        "confidence_score": [80.0] * n,
    })


# ---------------------------------------------------------------------------
# 1. Representative mapping
# ---------------------------------------------------------------------------


class TestRepresentativeMapping:
    def test_cluster_uses_representative_and_keeps_raw_genome(self):
        analyzer = NovelDiversityAnalyzer(
            _novel_df(), metadata=_metadata(),
            ani_matrix=ANIMatrix.from_dataframe(_ani_df()),
        )
        clusters = analyzer.cluster_novel_reads()
        assert len(clusters) == 1
        assert clusters[0].nearest_genome == REP
        assert clusters[0].contributing_genomes == [NONREP]

    def test_neighborhood_computed_for_nonrepresentative_hits(self):
        ani = ANIMatrix.from_dataframe(_ani_df())
        clusters = NovelDiversityAnalyzer(
            _novel_df(), metadata=_metadata(), ani_matrix=ani,
        ).cluster_novel_reads()
        nbr = PhylogeneticNeighborhoodAnalyzer(ani_matrix=ani)
        out = nbr.analyze(clusters)
        assert out[0].neighborhood is not None
        assert nbr.unplaced_cluster_ids == []

    def test_neighborhood_skips_genome_missing_from_matrix(self):
        ani = ANIMatrix.from_dataframe(_ani_df())
        clusters = NovelDiversityAnalyzer(
            _novel_df(), metadata=None, ani_matrix=ani,
        ).cluster_novel_reads()  # no metadata -> raw non-representative
        nbr = PhylogeneticNeighborhoodAnalyzer(ani_matrix=ani)
        out = nbr.analyze(clusters)
        assert out[0].neighborhood is None
        assert nbr.unplaced_cluster_ids == [clusters[0].cluster_id]


def _generate_report(tmp_path: Path, classifications: pl.DataFrame, metadata: bool) -> str:
    from metadarkmatter.cli.report import app

    cls_path = tmp_path / "cls.csv"
    ani_path = tmp_path / "ani.csv"
    out_path = tmp_path / "report.html"
    classifications.write_csv(cls_path)
    _ani_df().write_csv(ani_path)
    args = [
        "generate", "--classifications", str(cls_path), "--ani", str(ani_path),
        "--output", str(out_path),
    ]
    if metadata:
        meta_path = tmp_path / "meta.tsv"
        pl.DataFrame(
            {
                "accession": [REP, NONREP, *OTHERS],
                "species": ["Genus a", "Genus a", "Genus b", "Genus c", "Genus d"],
                "genus": ["Genus", "Genus", "Genus", "Genus2", "Genus2"],
                "family": ["Fam"] * 5,
                "representative": [REP, REP, *OTHERS],
            }
        ).write_csv(meta_path, separator="\t")
        args += ["--metadata", str(meta_path)]
    result = CliRunner().invoke(app, args)
    assert result.exit_code == 0, result.output
    return out_path.read_text()


def _tree_data(html: str) -> dict:
    m = re.search(r"TREE_DATA\s*=\s*(\{.*?\});\s*\n", html, re.S)
    assert m, "TREE_DATA not found"
    return json.loads(m.group(1))


class TestReportPaths:
    def test_report_places_nonrepresentative_and_ids_match(self, tmp_path):
        html = _generate_report(tmp_path, _novel_df(), metadata=True)
        tree = _tree_data(html)
        assert tree["novel_count"] == 1
        ann = tree["annotations"]["NSP_001"]
        assert ann["nearest_ref"] == REP
        # Novel Diversity tab uses the same ID and read count
        assert "NSP_001" in html
        assert ann["read_count"] == 6
        # Confidence in the tree is the analyzer's (confidence_score=80 -> Medium
        # with 6 reads), not the legacy tree-only rule.
        assert ann["confidence"] == "Medium"
        assert "could not be placed" not in html

    def test_unplaced_clusters_are_reported_visibly(self, tmp_path):
        # Without metadata the non-representative hit cannot be placed.
        html = _generate_report(tmp_path, _novel_df(), metadata=False)
        tree = _tree_data(html)
        assert tree["novel_count"] == 0
        assert "1 of 1 novel cluster(s) could not be placed" in tree["source_note"]
        assert "1 of 1 cluster(s) have a nearest genome that is not in the ANI matrix" in html


# ---------------------------------------------------------------------------
# 2. Tree uses the analyzer clusters; branch length rule
# ---------------------------------------------------------------------------


class TestTreeClusters:
    def test_branch_length_same_rule_for_species_and_genus(self):
        from metadarkmatter.core.phylogeny.placement import _estimate_branch_length

        assert _estimate_branch_length(20.0, "novel_genus") == pytest.approx(20.0)
        assert _estimate_branch_length(20.0, "novel_species") == pytest.approx(20.0)

    def test_from_diversity_cluster_copies_fields(self):
        from metadarkmatter.core.phylogeny import from_diversity_cluster

        c = NovelDiversityAnalyzer(
            _novel_df(), metadata=_metadata(),
        ).cluster_novel_reads()[0]
        p = from_diversity_cluster(c)
        assert (p.cluster_id, p.read_count, p.confidence_rating) == (
            c.cluster_id, c.read_count, c.confidence,
        )
        assert p.best_match_genome == c.nearest_genome
        assert p.mean_novelty == pytest.approx(c.mean_novelty_index, abs=0.01)

    def test_place_reports_unplaced(self):
        from metadarkmatter.core.phylogeny import ani_to_newick, place_novel_clusters
        from metadarkmatter.core.phylogeny.placement import NovelCluster as PCluster

        ani = _ani_df().to_pandas().set_index("genome")
        nw = ani_to_newick(ani)
        assert nw is not None
        bad = PCluster("NSP_001", "novel_species", "GCF_missing", 90, 10, 0, 5, "Low")
        unplaced: list[str] = []
        place_novel_clusters(nw, [bad], ani, unplaced=unplaced)
        assert unplaced == ["NSP_001"]


# ---------------------------------------------------------------------------
# 3. Cluster uncertainty and band label
# ---------------------------------------------------------------------------


def _cluster(**kw):
    base = {
        "cluster_id": "NSP_001", "taxonomic_call": "Novel Species",
        "nearest_genome": "G", "nearest_species": "s", "nearest_genus": "g",
        "nearest_family": "f", "novelty_band": 5, "read_count": 5,
        "mean_novelty_index": 12.0, "novelty_min": 5.0, "novelty_max": 18.0,
        "mean_placement_uncertainty": 0.0, "suggested_name": "x",
        "confidence": "Low", "phylogenetic_context": "c",
    }
    base.update(kw)
    return NovelCluster(**base)


class TestUncertaintyAndBand:
    def test_single_hit_uncertainty_not_fabricated(self):
        c = NovelDiversityAnalyzer(_novel_df(unc=0.0), metadata=None).cluster_novel_reads()[0]
        assert c.mean_placement_uncertainty == 0.0
        assert c.mean_inferred_uncertainty is None

    def test_merged_band_label_reflects_range(self):
        c = _cluster()  # lowest band 5 but range spans 5-20, mean 12
        assert c.novelty_band_label == "5-20% (Moderately divergent)"

    def test_single_band_label_unchanged_and_ge_25(self):
        assert _cluster(novelty_min=5.2, novelty_max=9.8, mean_novelty_index=7.5
                        ).novelty_band_label == "5-10% (Recently diverged)"
        assert _cluster(novelty_band=25, novelty_min=25.0, novelty_max=29.0,
                        mean_novelty_index=27.0).novelty_band_label == "25-30%"

    def test_real_merge_produces_range_label(self):
        df = pl.DataFrame({
            "read_id": [f"r{i}" for i in range(6)],
            "best_match_genome": ["G"] * 6,
            "novelty_index": [9.0, 9.5, 9.8, 10.2, 10.5, 10.9],
            "placement_uncertainty": [1.0] * 6,
            "taxonomic_call": ["Novel Species"] * 6,
        })
        c = NovelDiversityAnalyzer(df, min_cluster_size=3).cluster_novel_reads()[0]
        assert c.novelty_band == 5  # lowest merged band
        assert c.novelty_band_label.startswith("5-15%")


# ---------------------------------------------------------------------------
# 4. Novel tab labels / plot
# ---------------------------------------------------------------------------


class TestNovelTabLabels:
    def test_scatter_y_axis_informative_and_effective_thresholds(self):
        from metadarkmatter.visualization.report.novel_section import (
            build_cluster_scatter_figure,
        )

        cl = [
            _cluster(cluster_id="NSP_001", read_count=10, mean_bayesian_confidence=72.0,
                     confidence="High"),
            _cluster(cluster_id="NSP_002", read_count=4, mean_bayesian_confidence=55.0,
                     confidence="Low"),
        ]
        fig = build_cluster_scatter_figure(cl, novel_species_min=10.0, novel_genus_min=25.0)
        ys = [y for t in fig.data for y in t.y]
        assert len(set(ys)) > 1
        assert fig.layout.yaxis.title.text == "Mean Confidence Score"
        labels = [a.text for a in fig.layout.annotations]
        assert "Novel Species (10%)" in labels
        assert "Novel Genus (25%)" in labels

    def test_scatter_falls_back_to_read_count(self):
        from metadarkmatter.visualization.report.novel_section import (
            build_cluster_scatter_figure,
        )

        cl = [_cluster(cluster_id="NSP_001", read_count=10),
              _cluster(cluster_id="NSP_002", read_count=4)]
        fig = build_cluster_scatter_figure(cl)
        assert fig.layout.yaxis.title.text == "Reads in Cluster"

    def test_genus_boundary_trust(self):
        from metadarkmatter.core.classification.adaptive import AdaptiveGenusThreshold
        from metadarkmatter.visualization.report.sections.novel import NovelMixin

        def res(method, rng):
            return AdaptiveGenusThreshold(82.0, 18.0, 0.5, method, rng, 20)

        trusted = NovelMixin._genus_boundary_trusted
        assert trusted(res("gmm_3component", (78.0, 86.0)))
        assert not trusted(res("fallback", None))
        assert not trusted(res("gmm_3component", (60.0, 70.0)))  # diverges > 5

    def test_report_labels_units(self, tmp_path):
        html = _generate_report(tmp_path, _novel_df(), metadata=True)
        assert "Novel Species (reads)" in html
        assert "Novel Species clusters" in html


# ---------------------------------------------------------------------------
# 5. Tree builder missing ANI
# ---------------------------------------------------------------------------


class TestTreeMissingAni:
    def test_zero_is_missing_not_distance_100(self):
        from io import StringIO

        from Bio import Phylo

        from metadarkmatter.core.phylogeny.tree_builder import (
            ani_to_newick,
            count_missing_pairs,
        )

        g = ["A", "B", "C", "D"]
        m = pd.DataFrame(
            [[100, 95, 90, 0], [95, 100, 91, 80], [90, 91, 100, 82], [0, 80, 82, 100]],
            index=g, columns=g, dtype=float,
        )
        assert count_missing_pairs(m) == 1
        filled = m.copy()
        filled.loc["A", "D"] = filled.loc["D", "A"] = 70.0
        t0 = Phylo.read(StringIO(ani_to_newick(m)), "newick")
        t70 = Phylo.read(StringIO(ani_to_newick(filled)), "newick")
        for a in g:
            for b in g:
                if a < b:
                    assert t0.distance(a, b) == pytest.approx(t70.distance(a, b))

    def test_mirrored_pair_not_missing(self):
        from metadarkmatter.core.phylogeny.tree_builder import count_missing_pairs

        g = ["A", "B", "C"]
        m = pd.DataFrame(
            [[100, 90, 0], [90, 100, 85], [80, 85, 100]], index=g, columns=g, dtype=float
        )
        assert count_missing_pairs(m) == 0

    def test_tree_note_mentions_missing_pairs(self, tmp_path):
        ani = _ani_df().to_pandas().set_index("genome")
        ani.loc[REP, OTHERS[2]] = ani.loc[OTHERS[2], REP] = 0.0
        pl.from_pandas(ani.reset_index()).write_csv(tmp_path / "ani.csv")
        from metadarkmatter.cli.report import app

        cls = _novel_df(genome=REP)
        cls.write_csv(tmp_path / "cls.csv")
        r = CliRunner().invoke(app, [
            "generate", "--classifications", str(tmp_path / "cls.csv"),
            "--ani", str(tmp_path / "ani.csv"), "--output", str(tmp_path / "r.html"),
        ])
        assert r.exit_code == 0, r.output
        tree = _tree_data((tmp_path / "r.html").read_text())
        assert "1 genome pair(s) had no computed ANI" in tree["source_note"]


# ---------------------------------------------------------------------------
# 6. ANI CLI: skani mode, undetected-pair warning, directory expansion
# ---------------------------------------------------------------------------


class TestAniCli:
    def test_count_undetected_pairs(self):
        from metadarkmatter.cli.ani import count_undetected_pairs

        d = {
            "A": {"A": 100.0, "B": 90.0, "C": 0.0},
            "B": {"A": 90.0, "B": 100.0, "C": 0.0},
            "C": {"A": 0.0, "B": 0.0, "C": 100.0},
        }
        assert count_undetected_pairs(d) == (2, 3)

    def test_skani_build_command_modes(self, tmp_path):
        from metadarkmatter.external.skani import Skani

        files = [tmp_path / "b.fna", tmp_path / "a.fna"]
        for f in files:
            f.write_text(">x\nACGT\n")
        (tmp_path / "notes.txt").write_text("x")
        with patch.object(Skani, "get_executable", return_value=Path("skani")):
            cmd = Skani().build_command(genomes=tmp_path, output=tmp_path / "o", mode="slow")
            assert "*" not in " ".join(cmd)
            assert cmd[2:4] == [str(tmp_path / "a.fna"), str(tmp_path / "b.fna")]
            assert "--slow" in cmd
            cmd_c = Skani().build_command(
                genomes=files, output=tmp_path / "o", mode="slow", compression_factor=50
            )
            assert "--slow" not in cmd_c and cmd_c[cmd_c.index("-c") + 1] == "50"
            assert "--slow" not in Skani().build_command(genomes=files, output=tmp_path / "o")

    def test_skani_large_directory_uses_list_file(self, tmp_path):
        from metadarkmatter.external.skani import Skani

        gdir = tmp_path / "g"
        gdir.mkdir()
        for i in range(5):
            (gdir / f"g{i}.fna").write_text(">x\nA\n")
        with patch.object(Skani, "get_executable", return_value=Path("skani")), patch.object(
            Skani, "MAX_ARGV_GENOMES", 2
        ):
            cmd = Skani().build_command(genomes=gdir, output=tmp_path / "o.txt")
        assert "-l" in cmd
        lines = Path(cmd[cmd.index("-l") + 1]).read_text().splitlines()
        assert lines == sorted(lines) and len(lines) == 5


# ---------------------------------------------------------------------------
# 7. fastANI averaging, AAI mirroring / parsing / splitting / threads
# ---------------------------------------------------------------------------


class TestMatrixBuilders:
    def test_fastani_asymmetric_pairs_averaged(self, tmp_path):
        from metadarkmatter.core.ani_matrix_builder import parse_fastani_output

        f = tmp_path / "fa.out"
        f.write_text(
            "/x/GCF_000000001.1.fna\t/x/GCF_000000002.1.fna\t90.0\t10\t12\n"
            "/x/GCF_000000002.1.fna\t/x/GCF_000000001.1.fna\t90.5\t10\t12\n"
            "/x/GCF_000000001.1.fna\t/x/GCF_000000003.1.fna\t85.0\t10\t12\n"
        )
        d = parse_fastani_output(f)
        assert d["GCF_000000001.1"]["GCF_000000002.1"] == pytest.approx(90.25)
        assert d["GCF_000000002.1"]["GCF_000000001.1"] == pytest.approx(90.25)
        assert d["GCF_000000003.1"]["GCF_000000001.1"] == pytest.approx(85.0)
        assert d["GCF_000000002.1"]["GCF_000000003.1"] == 0.0

    def test_aai_matrix_mirrors_one_directional_pairs(self):
        from metadarkmatter.core.aai_matrix_builder import AAIMatrix

        m = AAIMatrix({"A": {"A": 100.0, "B": 60.0}, "B": {"A": 0.0, "B": 100.0}})
        assert m.get_aai("B", "A") == pytest.approx(60.0)
        assert m.get_aai("A", "B") == pytest.approx(60.0)

    def test_parse_aai_csv_matches_and_is_fast(self, tmp_path):
        import numpy as np

        from metadarkmatter.core.aai_matrix_builder import parse_aai_csv

        n = 300
        names = [f"G{i:04d}" for i in range(n)]
        rng = np.random.default_rng(0)
        arr = rng.uniform(40, 100, (n, n))
        df = pl.DataFrame({"genome": names, **{g: arr[:, j] for j, g in enumerate(names)}})
        p = tmp_path / "aai.csv"
        df.write_csv(p)
        t0 = time.perf_counter()
        d = parse_aai_csv(p)
        assert time.perf_counter() - t0 < 2.0
        assert d["G0003"]["G0007"] == pytest.approx(arr[3, 7])
        assert len(d) == n and len(d["G0000"]) == n

    def test_split_proteins_single_pass(self, tmp_path):
        from metadarkmatter.core.aai_matrix_builder import split_proteins_by_genome

        fa = tmp_path / "all.faa"
        fa.write_text(">A|p1\nMK\n>A|p2\nMM\n>B|p1\nMA\n>A|p3\nMT\n>C|p1\nMM\n")
        out = tmp_path / "out"
        out.mkdir()
        counts = split_proteins_by_genome(fa, ["A", "B", "C", "D"], out)
        assert counts == {"A": 3, "B": 1, "C": 1, "D": 0}
        assert (out / "A.faa").read_text().count(">") == 3
        assert not (out / "D.faa").exists()

    @pytest.mark.parametrize("threads", [1, 2, 8, 16, 64])
    def test_compute_aai_threads_not_oversubscribed(self, tmp_path, threads):
        import metadarkmatter.core.aai_matrix_builder as mod

        pdir = tmp_path / "prot"
        pdir.mkdir()
        for g in ("GCF_000000001.1", "GCF_000000002.1"):
            (pdir / f"{g}.faa").write_text(f">{g}|p1\nMKV\n")
        seen: dict[str, int] = {}

        class FakeDiamond:
            def make_database(self, *a, **k):
                pass

            def blastp(self, *, output, threads, **k):
                seen["tpj"] = threads
                Path(output).write_text("")

        real_pool = mod.ThreadPoolExecutor

        def pool(max_workers=None, **k):
            seen["workers"] = max_workers
            return real_pool(max_workers=max_workers, **k)

        with patch("metadarkmatter.external.diamond.Diamond", FakeDiamond), patch.object(
            mod, "ThreadPoolExecutor", pool
        ):
            try:
                mod.compute_aai_matrix(pdir, tmp_path / "o.csv", threads=threads)
            except Exception:  # empty hits -> downstream RBH may reject; thread plan is set already
                pass
        assert seen["tpj"] * seen["workers"] <= threads


# ---------------------------------------------------------------------------
# 8. Recruitment
# ---------------------------------------------------------------------------


def _sam(flag=0, cigar="100M", tags=()):
    return "\t".join(
        ["r1", str(flag), "GCF_000000001.1_scaffold_1", "10", "60", cigar, "*", "0", "0",
         "A" * 100, "I" * 100, *tags]
    ) + "\n"


class TestRecruitment:
    def test_identity_uses_nm_tag(self):
        from metadarkmatter.core.recruitment import calculate_percent_identity, parse_sam_line

        assert calculate_percent_identity("100M", nm_tag=10) == pytest.approx(90.0)
        rec = parse_sam_line(_sam(tags=("NM:i:20", "MD:Z:40A19T39")))
        assert rec is not None and rec.percent_identity == pytest.approx(80.0)

    def test_identity_falls_back_to_md_then_cigar(self):
        from metadarkmatter.core.recruitment import calculate_percent_identity

        # 2 mismatches + 1 deleted base (^A not counted as mismatch) + 3I
        assert calculate_percent_identity(
            "50M1D47M3I", md_tag="10A39^A20C26"
        ) == pytest.approx(100 * (1 - (2 + 3 + 1) / 101))
        assert calculate_percent_identity("100M") == 100.0

    def test_secondary_and_supplementary_skipped(self):
        from metadarkmatter.core.recruitment import parse_sam_line

        assert parse_sam_line(_sam(flag=256)) is None
        assert parse_sam_line(_sam(flag=2048)) is None
        assert parse_sam_line(_sam(flag=16)) is not None

    def test_stream_terminates_process_when_stopped_early(self, tmp_path):
        from metadarkmatter.core.recruitment import stream_bam_alignments

        bam = tmp_path / "x.bam"
        bam.touch()
        proc = MagicMock()
        proc.stdout = iter([_sam(), _sam(), _sam()])
        proc.poll.return_value = None
        proc.wait.return_value = -15
        with patch("metadarkmatter.core.recruitment.subprocess.Popen", return_value=proc):
            gen = stream_bam_alignments(bam)
            next(gen)
            gen.close()
        proc.terminate.assert_called_once()
        proc.wait.assert_called()

    def test_stream_raises_on_samtools_failure(self, tmp_path):
        from metadarkmatter.core.recruitment import stream_bam_alignments

        bam = tmp_path / "x.bam"
        bam.touch()
        proc = MagicMock()
        proc.stdout = iter([])
        proc.wait.return_value = 1
        with patch("metadarkmatter.core.recruitment.subprocess.Popen", return_value=proc):
            with pytest.raises(RuntimeError, match="samtools view failed"):
                list(stream_bam_alignments(bam))


# ---------------------------------------------------------------------------
# 9. Wrappers
# ---------------------------------------------------------------------------


class TestWrapperArgv:
    def test_bowtie2_interleaved(self):
        from metadarkmatter.external.bowtie2 import Bowtie2

        with patch.object(Bowtie2, "get_executable", return_value=Path("bowtie2")):
            cmd = Bowtie2().build_command(
                index_prefix=Path("idx"), reads_1=Path("r.fq"),
                output_sam=Path("o.sam"), interleaved=True,
            )
        i = cmd.index("--interleaved")
        assert cmd[i + 1] == "r.fq"
        assert "-U" not in cmd

    def test_kraken_taxid_list_is_separate_args(self):
        from metadarkmatter.external.kraken import ExtractKrakenReads

        with patch.object(
            ExtractKrakenReads, "get_executable", return_value=Path("extract_kraken_reads.py")
        ):
            cmd = ExtractKrakenReads().build_command(
                kraken_output=Path("a.kraken"), kraken_report=Path("a.kreport"),
                reads_1=Path("r.fq"), output_1=Path("o.fq"), taxid=[10, 20, 30],
            )
        i = cmd.index("-t")
        assert cmd[i + 1:i + 4] == ["10", "20", "30"]


# ---------------------------------------------------------------------------
# map reads: index reuse and empty input
# ---------------------------------------------------------------------------


class TestMapReads:
    def _make_index(self, prefix: Path):
        for sfx in (".1", ".2", ".3", ".4", ".rev.1", ".rev.2"):
            Path(f"{prefix}{sfx}.bt2").write_bytes(b"x")

    def test_index_is_current_requires_all_files_and_freshness(self, tmp_path):
        from metadarkmatter.cli.map import _index_is_current

        genome = tmp_path / "g.fna"
        genome.write_text(">a\nACGT\n")
        prefix = tmp_path / "idx"
        assert not _index_is_current(prefix, [genome])
        self._make_index(prefix)
        os.utime(genome, (time.time() - 100, time.time() - 100))
        assert _index_is_current(prefix, [genome])
        Path(f"{prefix}.rev.2.bt2").unlink()
        assert not _index_is_current(prefix, [genome])
        self._make_index(prefix)
        os.utime(genome, None)  # input newer than a (now older) index
        for f in tmp_path.glob("idx.*.bt2"):
            os.utime(f, (time.time() - 100, time.time() - 100))
        assert not _index_is_current(prefix, [genome])

    def test_empty_reads_detection(self, tmp_path):
        from metadarkmatter.cli.map import _reads_file_is_empty

        e = tmp_path / "e.fastq.gz"
        with gzip.open(e, "wb"):
            pass
        z = tmp_path / "z.fastq"
        z.write_bytes(b"")
        ok = tmp_path / "ok.fastq"
        ok.write_text("@r\nA\n+\nI\n")
        assert _reads_file_is_empty(e) and _reads_file_is_empty(z)
        assert not _reads_file_is_empty(ok)

    def test_cli_errors_on_empty_reads(self, tmp_path):
        from metadarkmatter.cli.map import app

        e = tmp_path / "e.fastq.gz"
        with gzip.open(e, "wb"):
            pass
        g = tmp_path / "g.fna"
        g.write_text(">a\nACGT\n")
        with patch("metadarkmatter.cli.map._check_tools_available", return_value=(True, [])):
            r = CliRunner().invoke(app, [
                "-1", str(e), "-g", str(g), "-o", str(tmp_path / "o.bam"),
            ])
        assert r.exit_code == 1
        assert "empty" in r.output.lower()
        assert "Mapping complete" not in r.output

    def test_cli_reuses_existing_index(self, tmp_path):
        from metadarkmatter.cli import map as map_mod

        reads = tmp_path / "s_R1.fastq"
        reads.write_text("@r\nACGT\n+\nIIII\n")
        gdir = tmp_path / "genomes"
        gdir.mkdir()
        gfile = gdir / "a.fna"
        gfile.write_text(">c\nACGTACGT\n")
        os.utime(gfile, (time.time() - 100, time.time() - 100))
        out = tmp_path / "out" / "s.bam"
        out.parent.mkdir()
        pangenome = out.parent / f"{map_mod._derive_output_name(reads)}_pangenome.fasta"
        pangenome.write_text(">a|c\nACGTACGT\n")
        self._make_index(pangenome.with_suffix(""))

        builder = MagicMock()
        aligner = MagicMock()
        aligner.run_or_raise.return_value = MagicMock(elapsed_seconds=0.1, stderr="",
                                                       command_string="b")
        samtools = MagicMock()
        r_ = MagicMock(elapsed_seconds=0.1)
        samtools.sam_to_sorted_bam.return_value = (r_, r_, r_)
        concat = MagicMock()
        with patch.object(map_mod, "_check_tools_available", return_value=(True, [])), patch.object(
            map_mod, "Bowtie2Build", return_value=builder
        ), patch.object(map_mod, "Bowtie2", return_value=aligner), patch.object(
            map_mod, "Samtools", return_value=samtools
        ), patch.object(map_mod, "concatenate_genomes", concat):
            r = CliRunner().invoke(map_mod.app, [
                "-1", str(reads), "-g", str(gdir), "-o", str(out),
            ])
        assert r.exit_code == 0, r.output
        builder.run_or_raise.assert_not_called()
        concat.assert_not_called()
