"""Tests for genus / reference grouping of novel reads (NovelDiversityAnalyzer)."""

from __future__ import annotations

import json

import polars as pl
import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.main import app
from metadarkmatter.core.classification.ani_matrix import ANIMatrix
from metadarkmatter.core.metadata import GenomeMetadata
from metadarkmatter.core.novel_diversity.clustering import NovelDiversityAnalyzer

# Two genera. Genus Ka has three references; KA1/KA2 are 88% ANI to each
# other, KA3 is 78% to both (a distant lineage of the same genus). Genus Kb
# has one reference at 85% ANI to KA1, so unrestricted 80% single linkage
# would chain both genera together.
GENOMES = {
    "KA1": ("Ka", "Ka one"), "KA2": ("Ka", "Ka two"), "KA3": ("Ka", "Ka three"),
    "KB1": ("Kb", "Kb one"),
}
ANI = {
    "KA1": {"KA2": 88.0, "KA3": 78.0, "KB1": 85.0},
    "KA2": {"KA3": 78.0, "KB1": 82.0},
    "KA3": {"KB1": 75.0},
}

# org1 spreads over KA1 and KA2 at mixed novelty bands; org2 hits KB1;
# org3 hits the distant KA3. Two org1 reads are Novel Genus.
READS = (
    [("org1", "KA1", 8.0 + i * 0.7, "Novel Species") for i in range(8)]
    + [("org1", "KA2", 9.0 + i * 0.9, "Novel Species") for i in range(6)]
    + [("org1", "KA2", 21.0, "Novel Genus"), ("org1", "KA1", 22.0, "Novel Genus")]
    + [("org2", "KB1", 12.0 + i * 0.3, "Novel Species") for i in range(5)]
    + [("org3", "KA3", 23.0 + i * 0.2, "Novel Genus") for i in range(4)]
    + [("known", "KA1", 1.0, "Known Species")] * 3
)


def _df() -> pl.DataFrame:
    return pl.DataFrame({
        "read_id": [f"{o}_r{i}" for i, (o, _g, _n, _c) in enumerate(READS)],
        "best_match_genome": [g for _o, g, _n, _c in READS],
        "novelty_index": [n for _o, _g, n, _c in READS],
        "placement_uncertainty": [1.0] * len(READS),
        "taxonomic_call": [c for _o, _g, _n, c in READS],
    })


def _ani() -> ANIMatrix:
    full: dict[str, dict[str, float]] = {g: {g: 100.0} for g in GENOMES}
    for a, row in ANI.items():
        for b, v in row.items():
            full[a][b] = v
            full[b][a] = v
    return ANIMatrix(full)


@pytest.fixture
def metadata_file(tmp_path):
    path = tmp_path / "genome_metadata.tsv"
    lines = ["accession\tspecies\tgenus\tfamily"] + [
        f"{acc}\t{sp}\t{genus}\tFam" for acc, (genus, sp) in GENOMES.items()
    ]
    path.write_text("\n".join(lines) + "\n")
    return path


def _groups(clusters) -> list[set[str]]:
    return sorted(({r.split("_")[0] for r in c.read_ids} for c in clusters), key=sorted)


class TestGenusGrouping:
    def test_auto_uses_genus_and_splits_distant_lineage(self, metadata_file):
        an = NovelDiversityAnalyzer(
            _df(), metadata=GenomeMetadata.from_file(metadata_file), ani_matrix=_ani()
        )
        clusters = an.cluster_novel_reads()
        assert an.grouping == "genus"
        by_org = {frozenset({r.split("_")[0] for r in c.read_ids}): c for c in clusters}
        assert set(by_org) == {frozenset({"org1"}), frozenset({"org2"}), frozenset({"org3"})}
        org1 = by_org[frozenset({"org1"})]
        # Not split by call or band; call is the majority call.
        assert org1.read_count == 16
        assert org1.taxonomic_call == "Novel Species"
        assert (org1.novel_species_reads, org1.novel_genus_reads) == (14, 2)
        assert org1.contributing_genomes == ["KA1", "KA2"]
        assert org1.nearest_genome == "KA1"  # most reads (9 vs 7)
        assert org1.nearest_genus == "Ka"
        assert org1.grouping == "genus"
        org3 = by_org[frozenset({"org3"})]
        assert org3.taxonomic_call == "Novel Genus"
        assert org3.cluster_id.startswith("NGN_")

    def test_without_ani_matrix_whole_genus_is_one_group(self, metadata_file):
        an = NovelDiversityAnalyzer(_df(), metadata=GenomeMetadata.from_file(metadata_file))
        assert _groups(an.cluster_novel_reads()) == [{"org1", "org3"}, {"org2"}]

    def test_legacy_chains_genera(self, metadata_file):
        """The legacy 80% single linkage merges both genera into one neighbourhood."""
        an = NovelDiversityAnalyzer(
            _df(), metadata=GenomeMetadata.from_file(metadata_file), ani_matrix=_ani(),
            grouping="legacy",
        )
        clusters = an.cluster_novel_reads()
        assert an.grouping == "legacy"
        # KB1 (genus Kb) is pulled into the KA1 neighbourhood: org2 is mixed with org1.
        assert any({"org1", "org2"} <= g for g in _groups(clusters))
        assert all(c.novel_species_reads is None for c in clusters)

    def test_deterministic(self, metadata_file):
        md = GenomeMetadata.from_file(metadata_file)
        a = [c.model_dump() for c in NovelDiversityAnalyzer(_df(), metadata=md, ani_matrix=_ani()).cluster_novel_reads()]
        shuffled = _df().sample(fraction=1.0, shuffle=True, seed=3)
        b = [c.model_dump() for c in NovelDiversityAnalyzer(shuffled, metadata=md, ani_matrix=_ani()).cluster_novel_reads()]
        assert a == b


class TestReferenceGrouping:
    def test_auto_falls_back_without_metadata(self, caplog):
        with caplog.at_level("WARNING"):
            an = NovelDiversityAnalyzer(_df(), ani_matrix=_ani())
            clusters = an.cluster_novel_reads()
        assert an.grouping == "reference"
        assert "No genus information" in caplog.text
        assert sorted(c.nearest_genome for c in clusters) == ["KA1", "KA2", "KA3", "KB1"]

    def test_explicit_reference(self, metadata_file):
        an = NovelDiversityAnalyzer(
            _df(), metadata=GenomeMetadata.from_file(metadata_file), grouping="reference"
        )
        assert len(an.cluster_novel_reads()) == 4

    def test_invalid_grouping(self):
        with pytest.raises(ValueError, match="grouping"):
            NovelDiversityAnalyzer(_df(), grouping="species")


class TestCli:
    def _cls(self, tmp_path):
        path = tmp_path / "cls.csv"
        _df().write_csv(path)
        return path

    def test_summarize_genus_with_ani(self, tmp_path, metadata_file):
        ani = tmp_path / "ani.csv"
        names, m = list(GENOMES), _ani()
        ani.write_text(
            "genome," + ",".join(names) + "\n"
            + "".join(f"{a}," + ",".join(str(m.get_ani(a, b)) for b in names) + "\n" for a in names)
        )
        out = tmp_path / "s.json"
        res = CliRunner().invoke(app, [
            "report", "summarize", str(self._cls(tmp_path)), "-o", str(out), "-f", "json",
            "--metadata", str(metadata_file), "--ani", str(ani), "-q",
        ])
        assert res.exit_code == 0, res.output
        data = json.loads(out.read_text())
        assert len(data["clusters"]) == 3
        assert {c["grouping"] for c in data["clusters"]} == {"genus"}

    def test_summarize_rejects_bad_grouping(self, tmp_path):
        res = CliRunner().invoke(app, [
            "report", "summarize", str(self._cls(tmp_path)), "-o", str(tmp_path / "s.tsv"),
            "--grouping", "species",
        ])
        assert res.exit_code == 1
        assert "Invalid --grouping" in res.output

    def test_generate_rejects_bad_grouping(self, tmp_path):
        res = CliRunner().invoke(app, [
            "report", "generate", "--classifications", str(self._cls(tmp_path)),
            "--output", str(tmp_path / "r.html"), "--cluster-grouping", "species",
        ])
        assert res.exit_code == 2
        assert "Invalid --cluster-grouping" in res.output
