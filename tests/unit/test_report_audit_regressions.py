"""Regression tests for report / novel-diversity / phylogeny audit findings.

Each test class corresponds to one audited bug:

1. Off-target rows (placement_uncertainty=null) crashed plots and the table.
2. Report output differed across processes (hash-seed dependent ordering).
3. Metadata species/genus names were injected unescaped into HTML/JS.
4. min_cluster_size was applied before adjacent-band merging; chained
   band merges stopped after the first merge.
5. "Ambiguous Within Genus" missing from summary grid; scatter plotted
   only four categories.
6. Tree branch-length units and label-based matrix indexing.
7. D3-missing fallback targeted a non-existent element id.
8. Empty / all-Off-target classification CSVs crashed the summary.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
import textwrap
from pathlib import Path

import pandas as pd
import polars as pl
import pytest

pytest.importorskip("plotly")

GENOMES = ["GCF_000001.1", "GCF_000002.1", "GCF_000003.1", "GCF_000004.1"]


def _ani_matrix() -> pl.DataFrame:
    return pl.DataFrame({
        "genome": GENOMES,
        "GCF_000001.1": [100.0, 95.5, 85.2, 78.0],
        "GCF_000002.1": [95.5, 100.0, 84.1, 77.5],
        "GCF_000003.1": [85.2, 84.1, 100.0, 92.0],
        "GCF_000004.1": [78.0, 77.5, 92.0, 100.0],
    })


def _classifications(with_off_target: bool = True) -> pl.DataFrame:
    n_known, n_sp, n_gen, n_awg = 20, 12, 10, 5
    n_ot = 4 if with_off_target else 0
    calls = (
        ["Known Species"] * n_known
        + ["Novel Species"] * n_sp
        + ["Novel Genus"] * n_gen
        + ["Ambiguous Within Genus"] * n_awg
        + ["Off-target"] * n_ot
    )
    n = len(calls)
    genomes = (
        ["GCF_000001.1"] * n_known
        + ["GCF_000002.1"] * n_sp
        + ["GCF_000003.1"] * n_gen
        + ["GCF_000004.1"] * n_awg
        + ["GCA_999999.1"] * n_ot
    )
    identity = [99.0] * n_known + [91.0] * n_sp + [78.0] * n_gen + [97.0] * n_awg + [99.0] * n_ot
    uncertainty: list[float | None] = (
        [0.5] * n_known + [0.8] * n_sp + [1.2] * n_gen + [6.0] * n_awg + [None] * n_ot
    )
    return pl.DataFrame({
        "read_id": [f"read_{i:03d}" for i in range(n)],
        "best_match_genome": genomes,
        "top_hit_identity": identity,
        "novelty_index": [100.0 - x for x in identity],
        "placement_uncertainty": uncertainty,
        "num_ambiguous_hits": [1, 2] * (n // 2) + [1] * (n % 2),
        "taxonomic_call": calls,
        "is_novel": [c in ("Novel Species", "Novel Genus") for c in calls],
    }, schema_overrides={"placement_uncertainty": pl.Float64})


def _metadata(species: dict[str, str] | None = None, genus: dict[str, str] | None = None):
    from metadarkmatter.core.metadata import GenomeMetadata

    species = species or {g: f"Francisella sp{i}" for i, g in enumerate(GENOMES)}
    genus = genus or dict.fromkeys(GENOMES, "Francisella")
    return GenomeMetadata(pl.DataFrame({
        "accession": GENOMES,
        "species": [species[g] for g in GENOMES],
        "genus": [genus[g] for g in GENOMES],
        "family": ["Francisellaceae"] * len(GENOMES),
        "gtdb_taxonomy": ["d__Bacteria"] * len(GENOMES),
    }))


def _render(df: pl.DataFrame, tmp_path: Path, **kwargs) -> str:
    from metadarkmatter.visualization.report.generator import ReportGenerator

    out = tmp_path / "report.html"
    ReportGenerator(df, **kwargs).generate(out)
    return out.read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# 1. Null placement_uncertainty (Off-target rows)
# ---------------------------------------------------------------------------


class TestNullUncertainty:
    def test_distribution_plots_ignore_nulls(self):
        from metadarkmatter.visualization.plots.distributions import (
            CombinedDistributionPlot,
            IdentityHistogram,
            NoveltyHistogram,
            UncertaintyHistogram,
        )

        df = _classifications()
        for cls in (NoveltyHistogram, UncertaintyHistogram, IdentityHistogram,
                    CombinedDistributionPlot):
            assert cls(df).create_figure() is not None

    def test_distribution_plots_all_null(self):
        from metadarkmatter.visualization.plots.distributions import UncertaintyHistogram

        df = pl.DataFrame({"placement_uncertainty": [None, None]},
                          schema={"placement_uncertainty": pl.Float64})
        assert UncertaintyHistogram(df).create_figure() is not None

    def test_report_with_off_target_embeds_null_uncertainty(self, tmp_path):
        from tests.unit.report_helpers import embedded_json

        html_out = _render(_classifications(), tmp_path, ani_matrix=_ani_matrix())
        payload = embedded_json(html_out, "mdm-table-data")
        calls = payload["dicts"]["classification"]
        off = calls.index("Off-target")
        rows = [
            i for i, c in enumerate(payload["cols"]["classification"]) if c == off
        ]
        assert rows
        # Off-target rows carry null uncertainty; the client renders "-"
        assert all(payload["cols"]["uncertainty"][i] is None for i in rows)


# ---------------------------------------------------------------------------
# 2. Cross-process determinism
# ---------------------------------------------------------------------------


_DETERMINISM_SCRIPT = textwrap.dedent('''
    import sys
    from pathlib import Path
    import polars as pl
    from metadarkmatter.core.metadata import GenomeMetadata
    from metadarkmatter.visualization.report.generator import ReportGenerator

    genomes = [f"GCF_{i:06d}.1" for i in range(1, 9)]
    # Many tied ANI values so any unordered iteration shows up in output.
    rows = {"genome": genomes}
    for j, b in enumerate(genomes):
        rows[b] = [100.0 if i == j else (90.0 if i // 2 == j // 2 else 76.0)
                   for i in range(len(genomes))]
    ani = pl.DataFrame(rows)
    read_genomes, ident, calls = [], [], []
    for i, g in enumerate(genomes):
        for _ in range(6):  # equal counts -> ties in every count sort
            read_genomes.append(g)
            ident.append(91.0 if i % 2 == 0 else 78.0)
            calls.append("Novel Species" if i % 2 == 0 else "Novel Genus")
    n = len(read_genomes)
    df = pl.DataFrame({
        "read_id": [f"r{i:04d}" for i in range(n)],
        "best_match_genome": read_genomes,
        "top_hit_identity": ident,
        "novelty_index": [100.0 - x for x in ident],
        "placement_uncertainty": [0.5] * n,
        "num_ambiguous_hits": [2] * n,
        "taxonomic_call": calls,
        "is_novel": [True] * n,
    })
    md = GenomeMetadata(pl.DataFrame({
        "accession": genomes,
        "species": [f"Sp {chr(65 + i // 2)}" for i in range(len(genomes))],
        "genus": [f"Genus{i // 4}" for i in range(len(genomes))],
        "family": ["Fam"] * len(genomes),
        "gtdb_taxonomy": ["d__Bacteria"] * len(genomes),
    }))
    ReportGenerator(df, ani_matrix=ani, genome_metadata=md).generate(Path(sys.argv[1]))
''')


class TestCrossProcessDeterminism:
    @pytest.mark.slow
    def test_report_identical_across_hash_seeds(self, tmp_path):
        script = tmp_path / "gen_report.py"
        script.write_text(_DETERMINISM_SCRIPT)
        outputs = []
        for seed in ("0", "12345"):
            out = tmp_path / f"report_{seed}.html"
            env = {**os.environ, "PYTHONHASHSEED": seed}
            subprocess.run(
                [sys.executable, str(script), str(out)],
                check=True, env=env, cwd=tmp_path, timeout=300,
            )
            text = out.read_text(encoding="utf-8")
            outputs.append(re.sub(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}", "<TS>", text))
        assert outputs[0] == outputs[1]

    def test_neighborhood_tiebreak_is_lexicographic(self):
        from metadarkmatter.core.classification.ani_matrix import ANIMatrix
        from metadarkmatter.core.novel_diversity.models import NovelCluster
        from metadarkmatter.core.novel_diversity.neighborhood import (
            PhylogeneticNeighborhoodAnalyzer,
        )

        genomes = ["G_d", "G_c", "G_b", "G_a"]
        ani = {g: {h: (100.0 if g == h else 80.0) for h in genomes} for g in genomes}
        matrix = ANIMatrix(ani)
        cluster = NovelCluster(
            cluster_id="NSP_001", taxonomic_call="Novel Species",
            nearest_genome="G_d", nearest_species="x", nearest_genus="X",
            nearest_family="F", novelty_band=5, read_count=5,
            mean_novelty_index=7.0, novelty_min=6.0, novelty_max=8.0,
            mean_placement_uncertainty=0.5, suggested_name="X sp.",
            confidence="High", phylogenetic_context="ctx", read_ids=["r1"],
        )
        genus_map = {"G_d": "Zeta", "G_c": "Beta", "G_b": "Alpha", "G_a": "Alpha"}
        nbr = PhylogeneticNeighborhoodAnalyzer(matrix, genus_map=genus_map).analyze([cluster])[0]
        assert nbr.neighborhood is not None
        genera = nbr.neighborhood.nearest_genera
        # Zeta (nearest genome) first; tied Alpha/Beta ordered by name, and
        # Alpha's representative is the lexicographically smallest genome.
        assert [g.genus for g in genera] == ["Zeta", "Alpha", "Beta"]
        assert genera[1].representative_genome == "G_a"

    def test_select_relevant_references_tiebreak(self):
        from metadarkmatter.core.novel_diversity.models import NovelCluster
        from metadarkmatter.visualization.report.components.extended_matrix_builder import (
            select_relevant_references,
        )

        all_genomes = [f"G{i:02d}" for i in range(12)]
        ani = {(a, b): 90.0 for a in all_genomes for b in all_genomes if a != b}

        def cl(cid: str, g: str) -> NovelCluster:
            return NovelCluster(
                cluster_id=cid, taxonomic_call="Novel Species", nearest_genome=g,
                nearest_species="x", nearest_genus="X", nearest_family="F",
                novelty_band=5, read_count=5, mean_novelty_index=7.0,
                novelty_min=6.0, novelty_max=8.0, mean_placement_uncertainty=0.5,
                suggested_name="X sp.", confidence="High",
                phylogenetic_context="ctx", read_ids=["r1"],
            )

        clusters = [cl("NSP_001", "G09"), cl("NSP_002", "G03")]
        selected = select_relevant_references(ani, clusters, all_genomes, max_refs=4)
        # Nearest genomes in sorted order (G03 first), then ties by name.
        assert selected == ["G00", "G01", "G03", "G09"]


# ---------------------------------------------------------------------------
# 3. HTML / JS injection from metadata names
# ---------------------------------------------------------------------------


EVIL = '"tul</script><img src=x onerror=alert(1)>'


class TestMetadataEscaping:
    def test_report_escapes_species_and_genus(self, tmp_path):
        species = dict.fromkeys(GENOMES, f"Francisella {EVIL}")
        genus = dict.fromkeys(GENOMES, "Allo<b>genus</b>")
        html_out = _render(
            _classifications(), tmp_path,
            ani_matrix=_ani_matrix(), genome_metadata=_metadata(species, genus),
        )
        assert "<img src=x" not in html_out
        assert "<b>genus</b>" not in html_out

    def test_json_for_script_escapes_script_close(self):
        import json

        from metadarkmatter.visualization.report.sections.phylogeny import (
            _json_for_script,
        )

        payload = {"label": "a</script><b>&"}
        out = _json_for_script(payload)
        assert "<" not in out and ">" not in out and "&" not in out
        assert json.loads(out) == payload

    def test_tree_tooltip_escapes_labels(self):
        from metadarkmatter.visualization.report.templates import PHYLOTREE_JS_TEMPLATE

        assert "escapeHtml(label)" in PHYLOTREE_JS_TEMPLATE
        assert "'<strong>' + label +" not in PHYLOTREE_JS_TEMPLATE

    def test_novel_cluster_table_escapes(self):
        from metadarkmatter.core.novel_diversity.models import NovelCluster
        from metadarkmatter.visualization.report.novel_section import (
            build_cluster_table_html,
            build_phylogenetic_context_html,
        )

        cluster = NovelCluster(
            cluster_id="NSP_001", taxonomic_call="Novel Species",
            nearest_genome="G1", nearest_species=EVIL, nearest_genus="<i>G</i>",
            nearest_family="F", novelty_band=5, read_count=5,
            mean_novelty_index=7.0, novelty_min=6.0, novelty_max=8.0,
            mean_placement_uncertainty=0.5, suggested_name=f"{EVIL} sp.",
            confidence="High", phylogenetic_context=f"Novel species within {EVIL}",
            read_ids=["r1"],
        )
        for out in (build_cluster_table_html([cluster]),
                    build_phylogenetic_context_html([cluster])):
            assert "<img" not in out
            assert "<i>G</i>" not in out

    def test_multi_sample_title_and_names_escaped(self):
        from metadarkmatter.visualization.report.multi_generator import (
            MultiSampleConfig,
            MultiSampleReportGenerator,
        )

        df = _classifications(with_off_target=False)
        gen = MultiSampleReportGenerator(
            {"<s>evil</s>": df, "ok": df},
            MultiSampleConfig(title="<script>alert(1)</script>"),
        )
        out = gen._build_html()
        assert "<script>alert(1)</script>" not in out
        assert "<s>evil</s>" not in out


# ---------------------------------------------------------------------------
# 4. Novel clustering: merge before size filter; chained merges
# ---------------------------------------------------------------------------


def _novel_df(novelties: list[float]) -> pl.DataFrame:
    n = len(novelties)
    return pl.DataFrame({
        "read_id": [f"r{i}" for i in range(n)],
        "best_match_genome": ["GCF_000001.1"] * n,
        "novelty_index": novelties,
        "top_hit_identity": [100.0 - x for x in novelties],
        "placement_uncertainty": [0.5] * n,
        "taxonomic_call": ["Novel Species"] * n,
    })


class TestClusterMerging:
    def test_sub_threshold_adjacent_bands_merge_before_filter(self):
        from metadarkmatter.core.novel_diversity import NovelDiversityAnalyzer

        # Two reads in band 5, two in band 10: neither alone reaches 3.
        df = _novel_df([9.0, 9.5, 10.2, 10.5])
        clusters = NovelDiversityAnalyzer(df, min_cluster_size=3).cluster_novel_reads()
        assert len(clusters) == 1
        assert clusters[0].read_count == 4

    def test_chained_band_merge(self):
        from metadarkmatter.core.novel_diversity import NovelDiversityAnalyzer

        # Bands 5, 10, 15 with near-contiguous novelty ranges.
        df = _novel_df([9.5, 10.0, 14.5, 15.5])
        clusters = NovelDiversityAnalyzer(df, min_cluster_size=1).cluster_novel_reads()
        assert len(clusters) == 1
        assert clusters[0].read_count == 4
        assert clusters[0].novelty_band == 5


# ---------------------------------------------------------------------------
# 5. Category coverage in summary grid and scatter
# ---------------------------------------------------------------------------


class TestCategoryCoverage:
    def test_summary_grid_includes_ambiguous_within_genus(self):
        from metadarkmatter.visualization.report.generator import ReportGenerator

        section = ReportGenerator(_classifications())._build_summary_section()
        assert "Ambiguous Within Genus" in section

    def test_scatter_includes_all_present_categories(self):
        from metadarkmatter.visualization.plots.scatter_2d import (
            NoveltyUncertaintyScatter,
        )

        fig = NoveltyUncertaintyScatter(
            _classifications(), enable_single_hit_toggle=False
        ).create_figure()
        names = " ".join(t.name for t in fig.data if t.name)
        assert "Ambiguous Within Genus" in names
        # Off-target rows have null uncertainty and cannot be placed.
        assert "Off-target" not in names


# ---------------------------------------------------------------------------
# 6. Tree units and label-based matrix indexing
# ---------------------------------------------------------------------------


def _ani_pd(order: list[str]) -> pd.DataFrame:
    return _ani_matrix().to_pandas().set_index("genome").loc[order]


class TestTreeConsistency:
    def test_branch_length_in_percent_units(self):
        from metadarkmatter.core.phylogeny.placement import _estimate_branch_length

        assert _estimate_branch_length(5.0, "novel_species") == pytest.approx(5.0)
        assert _estimate_branch_length(20.0, "novel_genus") == pytest.approx(20.0)

    def test_tree_independent_of_row_order(self):
        from io import StringIO

        from Bio import Phylo

        from metadarkmatter.core.phylogeny.tree_builder import ani_to_newick

        def dists(newick: str) -> dict[tuple[str, str], float]:
            tree = Phylo.read(StringIO(newick), "newick")
            return {
                (a, b): round(tree.distance(a, b), 6)
                for a in GENOMES for b in GENOMES if a < b
            }

        a = ani_to_newick(_ani_pd(GENOMES))
        b = ani_to_newick(_ani_pd(list(reversed(GENOMES))))
        assert a is not None and b is not None
        assert dists(a) == dists(b)
        # Distances are in percent (100 - ANI): 01-02 ANI 95.5 -> 4.5
        assert dists(a)[("GCF_000001.1", "GCF_000002.1")] == pytest.approx(4.5, abs=0.01)

    def test_heatmap_independent_of_row_order(self):
        from metadarkmatter.visualization.report.components.heatmap_builder import (
            ANI_HEATMAP_CONFIG,
            _build_similarity_heatmap,
        )

        shuffled = _ani_matrix()[[2, 0, 3, 1]]
        fig, _stats, _ok = _build_similarity_heatmap(shuffled, {}, ANI_HEATMAP_CONFIG)
        z = fig.data[0].z
        labels = list(fig.data[0].x)
        ref = _ani_matrix().to_pandas().set_index("genome")
        for i, a in enumerate(labels):
            for j, b in enumerate(labels):
                assert z[i][j] == pytest.approx(ref.loc[a, b])


# ---------------------------------------------------------------------------
# 7. D3 fallback element id
# ---------------------------------------------------------------------------


def test_d3_fallback_targets_existing_container():
    from metadarkmatter.visualization.report.templates import (
        PHYLOGENY_SECTION_TEMPLATE,
        PHYLOTREE_JS_TEMPLATE,
    )

    assert 'id="phylogeny-container"' in PHYLOGENY_SECTION_TEMPLATE
    assert "getElementById('phylogeny-tree')" not in PHYLOTREE_JS_TEMPLATE


# ---------------------------------------------------------------------------
# 8. Empty / all-Off-target CSV
# ---------------------------------------------------------------------------


_CSV_HEADER = (
    "read_id,best_match_genome,top_hit_identity,novelty_index,placement_uncertainty,"
    "num_ambiguous_hits,identity_gap,confidence_score,taxonomic_call,is_novel,"
    "p_known_species,posterior_entropy,inferred_uncertainty"
)


class TestDegenerateInputs:
    def test_empty_csv(self, tmp_path):
        csv = tmp_path / "empty.csv"
        csv.write_text(_CSV_HEADER + "\n")
        html_out = _render(pl.read_csv(csv), tmp_path, ani_matrix=_ani_matrix())
        assert "<html" in html_out.lower()

    def test_all_off_target_csv(self, tmp_path):
        csv = tmp_path / "ot.csv"
        row = "r{i},GCA_9.1,99.0,1.0,,0,,10.0,Off-target,false,,,"
        csv.write_text(_CSV_HEADER + "\n" + "\n".join(row.format(i=i) for i in range(3)) + "\n")
        df = pl.read_csv(csv)
        assert df.schema["placement_uncertainty"] == pl.String  # precondition
        html_out = _render(df, tmp_path, ani_matrix=_ani_matrix())
        assert "<html" in html_out.lower()
