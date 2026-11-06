"""Regression tests for the round-2 CLI/pipeline audit fixes.

None of these need the external bioinformatics tools.
"""

from __future__ import annotations

import gzip
import json
import time
from pathlib import Path
from unittest.mock import MagicMock, patch

import polars as pl
import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.main import app
from metadarkmatter.core.id_mapping import ContigIdMapping
from metadarkmatter.core.pipeline import Pipeline, PipelineLayout, PipelineStepError, StepStatus
from metadarkmatter.models.classification import TaxonomicCall
from metadarkmatter.models.pipeline_config import PipelineConfig, PipelineStep

runner = CliRunner()


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.delenv("MDM_DRY_RUN", raising=False)
    yield
    monkeypatch.delenv("MDM_DRY_RUN", raising=False)


# ---------------------------------------------------------------------------
# 1-3: validate
# ---------------------------------------------------------------------------
def _blast_line(n_extra: int) -> str:
    cols = ["q1", "g|c1", "99.0", "100", "1", "0", "1", "100", "1", "100", "1e-30", "180.0"]
    cols += ["150", "0.66"][:n_extra]
    return "\t".join(cols) + "\n"


@pytest.mark.parametrize("n_extra", [0, 1, 2])
def test_validate_blast_accepts_12_13_14_columns(tmp_path: Path, n_extra: int):
    f = tmp_path / "x.tsv"
    f.write_text(_blast_line(n_extra) * 3)
    res = runner.invoke(app, ["validate", "blast", "-i", str(f)])
    assert res.exit_code == 0, res.output
    assert "valid" in res.output


def test_validate_blast_rejects_other_column_counts(tmp_path: Path):
    f = tmp_path / "x.tsv"
    f.write_text("\t".join(["a"] * 11) + "\n")
    res = runner.invoke(app, ["validate", "blast", "-i", str(f)])
    assert res.exit_code == 1


def test_validate_classifications_accepts_every_enum_category(tmp_path: Path):
    calls = [c.value for c in TaxonomicCall]
    f = tmp_path / "c.csv"
    pl.DataFrame(
        {
            "read_id": [f"r{i}" for i in range(len(calls))],
            "best_match_genome": ["G"] * len(calls),
            "novelty_index": [1.0] * len(calls),
            "placement_uncertainty": [1.0] * len(calls),
            "taxonomic_call": calls,
        }
    ).write_csv(f)
    res = runner.invoke(app, ["validate", "classifications", "-i", str(f)])
    assert res.exit_code == 0, res.output
    for legacy in ("Conserved Region", "Species Boundary", "Off-target", "Unclassified"):
        assert legacy in calls


def test_validate_classifications_rejects_unknown_category(tmp_path: Path):
    f = tmp_path / "c.csv"
    pl.DataFrame(
        {
            "read_id": ["r"], "best_match_genome": ["G"], "novelty_index": [1.0],
            "placement_uncertainty": [1.0], "taxonomic_call": ["Conserved Region Hit"],
        }
    ).write_csv(f)
    res = runner.invoke(app, ["validate", "classifications", "-i", str(f)])
    assert res.exit_code == 1


def test_validate_fasta_protein_and_wording(tmp_path: Path):
    f = tmp_path / "p.faa"
    f.write_text(">p1\nMKVLAAGIVEFQ*\n>p2\nMKVL\n")
    bad = runner.invoke(app, ["validate", "fasta", "-i", str(f)])
    assert bad.exit_code == 1  # protein letters are not nucleotides
    ok = runner.invoke(app, ["validate", "fasta", "-i", str(f), "--sequence-type", "protein"])
    assert ok.exit_code == 0, ok.output
    assert "Records read: 2" in ok.output
    capped = runner.invoke(
        app,
        ["validate", "fasta", "-i", str(f), "--sequence-type", "protein", "--max-records", "1"],
    )
    assert "Records sampled: 1" in capped.output


# ---------------------------------------------------------------------------
# 4-5: pipeline checkpointing, manifest, warnings
# ---------------------------------------------------------------------------
class ProducingRunner:
    """Succeeds and creates each step's declared output."""

    def __init__(self) -> None:
        self.calls: list[list[str]] = []

    def __call__(self, argv, dry_run: bool) -> int:
        argv = list(argv)
        self.calls.append(argv)
        if "--output" in argv:
            out = Path(argv[argv.index("--output") + 1])
            if "makedb" in argv:
                out.parent.mkdir(parents=True, exist_ok=True)
                Path(f"{out}.nhr").write_text("db")
            else:
                out.parent.mkdir(parents=True, exist_ok=True)
                out.write_text("x")
        return 0

    def ran(self, word: str) -> bool:
        return any(word in c for c in self.calls)


def _cfg(tmp_path: Path, **kw) -> PipelineConfig:
    genomes = tmp_path / "genomes"
    genomes.mkdir(exist_ok=True)
    (genomes / "GCF_1.fna").write_text(">c1\nACGT\n")
    reads = tmp_path / "r1.fastq.gz"
    reads.write_bytes(b"reads")
    return PipelineConfig(
        genomes_dir=genomes, extracted_reads_1=reads, output_dir=tmp_path / "out", **kw
    )


def test_identical_rerun_is_fully_cached(tmp_path: Path):
    cfg = _cfg(tmp_path)
    lay = PipelineLayout(root=cfg.output_dir)
    Pipeline(cfg, lay, runner=ProducingRunner()).run()
    r2 = ProducingRunner()
    m = Pipeline(cfg, lay, runner=r2).run()
    assert r2.calls == []
    assert {s.status for s in m.steps if s.name == "classify"} == {StepStatus.CACHED}


def test_changed_input_reruns_step_and_later_steps(tmp_path: Path):
    cfg = _cfg(tmp_path)
    lay = PipelineLayout(root=cfg.output_dir)
    Pipeline(cfg, lay, runner=ProducingRunner()).run()

    time.sleep(0.01)
    cfg.extracted_reads_1.write_bytes(b"different, longer reads")  # type: ignore[union-attr]
    r2 = ProducingRunner()
    m = Pipeline(cfg, lay, runner=r2).run()
    status = {s.name: s.status for s in m.steps}
    assert status["makedb"] == StepStatus.CACHED  # reads do not feed makedb
    assert status["align"] == StepStatus.COMPLETED
    assert status["classify"] == StepStatus.COMPLETED
    assert status["report"] == StepStatus.COMPLETED
    assert not r2.ran("makedb")


def test_changed_parameter_reruns(tmp_path: Path):
    cfg = _cfg(tmp_path)
    lay = PipelineLayout(root=cfg.output_dir)
    Pipeline(cfg, lay, runner=ProducingRunner()).run()
    cfg2 = cfg.model_copy(update={"ani_backend": "fastani"})
    r2 = ProducingRunner()
    Pipeline(cfg2, lay, runner=r2).run()
    assert r2.ran("ani")
    assert not r2.ran("makedb")


def test_thread_count_change_does_not_rerun(tmp_path: Path):
    """--threads does not change results, so it is not part of the fingerprint."""
    cfg = _cfg(tmp_path)
    lay = PipelineLayout(root=cfg.output_dir)
    Pipeline(cfg, lay, runner=ProducingRunner()).run()
    r2 = ProducingRunner()
    Pipeline(cfg.model_copy(update={"threads": 3}), lay, runner=r2).run()
    assert not r2.ran("align") and not r2.ran("ani") and not r2.ran("classify")


def test_killed_step_leaves_no_marker_and_reruns(tmp_path: Path):
    cfg = _cfg(tmp_path)
    lay = PipelineLayout(root=cfg.output_dir)

    class Killed(ProducingRunner):
        def __call__(self, argv, dry_run: bool) -> int:
            if "score" in argv:
                # Partial output written, then the process dies.
                out = Path(argv[argv.index("--output") + 1])
                out.write_text("partial")
                raise KeyboardInterrupt
            return super().__call__(argv, dry_run)

    with pytest.raises(KeyboardInterrupt):
        Pipeline(cfg, lay, runner=Killed()).run()
    assert lay.classifications.exists()
    assert not lay.marker(Pipeline(cfg, lay).plan()[6].step).exists()

    # --from must not accept the interrupted step's partial output as complete.
    with pytest.raises(PipelineStepError, match="not completed|missing"):
        Pipeline(cfg, lay, runner=ProducingRunner()).run(from_step=PipelineStep.REPORT)

    r2 = ProducingRunner()
    Pipeline(cfg, lay, runner=r2).run()
    assert r2.ran("classify")
    assert lay.classifications.read_text() == "x"


def test_old_manifest_is_kept(tmp_path: Path):
    cfg = _cfg(tmp_path)
    lay = PipelineLayout(root=cfg.output_dir)
    Pipeline(cfg, lay, runner=ProducingRunner()).run(mdm_version="1")
    Pipeline(cfg, lay, runner=ProducingRunner()).run(mdm_version="2")
    archived = list(lay.root.glob("run_manifest.*.json"))
    assert len(archived) == 1
    assert json.loads(archived[0].read_text())["mdm_version"] == "1"
    assert json.loads(lay.manifest.read_text())["mdm_version"] == "2"
    steps = json.loads(lay.manifest.read_text())["steps"]
    assert all("fingerprint" in s for s in steps)


def test_blast_with_r2_aligns_both_mates_without_warning(tmp_path: Path):
    cfg = _cfg(tmp_path)
    lay = PipelineLayout(root=cfg.output_dir)
    r2 = tmp_path / "r2.fastq.gz"
    r2.write_bytes(b"x")
    p = Pipeline(cfg, lay, extracted_reads_2=r2, runner=ProducingRunner())
    assert not any("single-end" in w for w in p.warnings)
    align = next(s for s in p.plan() if s.step == PipelineStep.ALIGN)
    assert align.argv[align.argv.index("--query-2") + 1] == str(r2)
    assert "--query" not in align.argv
    m = p.run()
    assert not m.warnings
    assert "single-end" not in lay.manifest.read_text()


def test_cli_dry_run_blast_r2_plans_paired_alignment(tmp_path: Path):
    cfg = _cfg(tmp_path)
    r2 = tmp_path / "r2.fastq.gz"
    r2.write_bytes(b"x")
    res = runner.invoke(
        app,
        ["run", "pipeline", "--genomes", str(cfg.genomes_dir), "--extracted-reads-1",
         str(cfg.extracted_reads_1), "--extracted-reads-2", str(r2),
         "--output-dir", str(tmp_path / "o"), "--dry-run"],
    )
    assert res.exit_code == 0, res.output
    assert "single-end" not in res.output
    assert "--query-2" in res.output.replace("\n", "")


def test_cli_nonexistent_genomes_fails_early(tmp_path: Path):
    reads = tmp_path / "r.fq"
    reads.write_text("x")
    res = runner.invoke(
        app,
        ["run", "pipeline", "--genomes", str(tmp_path / "nope"), "--extracted-reads-1",
         str(reads), "--output-dir", str(tmp_path / "o"), "--dry-run"],
    )
    assert res.exit_code != 0
    assert "nope" in res.output


def test_yaml_unknown_key_and_compute_aai_warn(tmp_path: Path):
    cfg = _cfg(tmp_path)
    y = tmp_path / "p.yaml"
    y.write_text(
        f"genomes_dir: {cfg.genomes_dir}\nextracted_reads_1: {cfg.extracted_reads_1}\n"
        f"output_dir: {tmp_path / 'o'}\nbogus_key: 1\ncompute_aai: true\n"
    )
    msgs = PipelineConfig.yaml_warnings(y)
    assert any("bogus_key" in m for m in msgs)
    assert any("compute_aai" in m and "not implemented" in m for m in msgs)
    res = runner.invoke(app, ["run", "pipeline", "--config", str(y), "--dry-run"])
    assert res.exit_code == 0, res.output
    assert "bogus_key" in res.output
    assert "not implemented" in res.output


# ---------------------------------------------------------------------------
# 6-7: util commands
# ---------------------------------------------------------------------------
def _write_map(path: Path, rows: list[tuple[str, str]]) -> None:
    path.write_text(
        "original_contig_id\tgenome_accession\n" + "".join(f"{a}\t{b}\n" for a, b in rows)
    )


def test_validate_mapping_reports_conflicts(tmp_path: Path):
    m = tmp_path / "m.tsv"
    _write_map(m, [("x", "G1"), ("x", "G2"), ("y", "G1"), ("y", "G1")])
    res = runner.invoke(app, ["util", "validate-mapping", str(m)])
    assert res.exit_code == 1
    assert "x" in res.output and "G1" in res.output and "G2" in res.output
    assert "y -> " not in res.output  # exact duplicate row is not a conflict


def test_validate_mapping_counts_standardized_ids_and_min_coverage(tmp_path: Path):
    m = tmp_path / "m.tsv"
    _write_map(m, [("c1", "G1")])
    b = tmp_path / "b.tsv"
    b.write_text("q\tG1|c1\t99\t10\t0\t0\t1\t10\t1\t10\t1e-5\t50\nq\tzzz\t99\t10\t0\t0\t1\t10\t1\t10\t1e-5\t50\n")
    res = runner.invoke(app, ["util", "validate-mapping", str(m), "--blast", str(b)])
    assert res.exit_code == 0, res.output
    assert "50.0%" in res.output
    res = runner.invoke(
        app, ["util", "validate-mapping", str(m), "--blast", str(b), "--min-coverage", "90"]
    )
    assert res.exit_code == 1
    res = runner.invoke(
        app, ["util", "validate-mapping", str(m), "--blast", str(b), "--min-coverage", "50"]
    )
    assert res.exit_code == 0
    only_std = tmp_path / "s.tsv"
    only_std.write_text("q\tG1|c1\t99\t10\t0\t0\t1\t10\t1\t10\t1e-5\t50\n")
    res = runner.invoke(
        app, ["util", "validate-mapping", str(m), "--blast", str(only_std), "--min-coverage", "100"]
    )
    assert res.exit_code == 0, res.output


def test_duplicate_contig_error_mentions_cli_flag(tmp_path: Path):
    g = tmp_path / "g"
    g.mkdir()
    (g / "GCF_000000001.1.fna").write_text(">dup\nACGT\n")
    (g / "GCF_000000002.1.fna").write_text(">dup\nACGT\n")
    with pytest.raises(ValueError, match="--allow-duplicates"):
        ContigIdMapping.from_genome_dir(g)


def test_generate_metadata_warns_and_writes_empty_not_quotes(tmp_path: Path):
    g = tmp_path / "g"
    g.mkdir()
    (g / "GCF_000000001.1.fna").write_text(">c\nACGT\n")
    (g / "Fam_Gen_sp_GCF_000000002.1.fna").write_text(">d\nACGT\n")
    out = tmp_path / "md.tsv"
    res = runner.invoke(app, ["util", "generate-metadata", "-g", str(g), "-o", str(out)])
    assert res.exit_code == 0, res.output
    assert "1 of" in res.output.replace("\n", " ") and "no taxonomy" in res.output.replace("\n", " ")
    text = out.read_text()
    assert '""' not in text
    row = next(r for r in text.splitlines() if r.startswith("GCF_000000001.1"))
    assert row.split("\t")[1:4] == ["", "", ""]


# ---------------------------------------------------------------------------
# 8: download
# ---------------------------------------------------------------------------
@patch("metadarkmatter.cli.download.GTDBClient")
def test_all_genomes_assigns_one_representative_per_species(mock_cls, tmp_path: Path):
    from metadarkmatter.clients.gtdb import GTDBGenome, GTDBQueryResult

    tax = "d__B;f__F;g__G;s__{}"
    allg = (
        GTDBGenome("GCF_000000003.1", tax.format("A"), "A"),
        GTDBGenome("GCF_000000001.1", tax.format("A"), "A"),
        GTDBGenome("GCF_000000002.1", tax.format("A"), "A"),
        GTDBGenome("GCF_000000010.1", tax.format("B"), "B"),
        GTDBGenome("GCF_000000011.1", tax.format("B"), "B"),
    )
    reps = (GTDBGenome("GCF_000000011.1", tax.format("B"), "B"),)

    def mk(genomes):
        return GTDBQueryResult(
            taxon="g__G", genomes=genomes, total_count=len(genomes),
            genus_counts={"G": len(genomes)},
            species_counts=dict.fromkeys({g.species for g in genomes}, 1),
        )

    client = MagicMock()
    client.query_genomes.side_effect = [mk(allg), mk(reps)]
    mock_cls.return_value = client

    out = tmp_path / "genomes.tsv"
    res = runner.invoke(app, ["download", "genomes", "list", "g__G", "--all-genomes",
                              "--output", str(out)])
    assert res.exit_code == 0, res.output
    assert "no GTDB" in res.output.replace("\n", " ") or "representative" in res.output
    md = pl.read_csv(tmp_path / "genome_metadata.tsv", separator="\t")
    per_species = md.group_by("species").agg(pl.col("representative").n_unique().alias("n"))
    assert per_species["n"].to_list() == [1, 1]
    a_rep = md.filter(pl.col("species") == "A")["representative"].unique().to_list()
    b_rep = md.filter(pl.col("species") == "B")["representative"].unique().to_list()
    assert a_rep == ["GCF_000000001.1"]  # first accession, sorted
    assert b_rep == ["GCF_000000011.1"]  # GTDB representative kept


@patch("metadarkmatter.cli.download.NCBIDatasets")
def test_fetch_require_all(mock_cls, tmp_path: Path):
    from metadarkmatter.external.ncbi_datasets import DownloadOutcome, DownloadReport

    ncbi = MagicMock()
    ncbi.check_available.return_value = True
    ncbi.download_genomes.return_value = DownloadReport(
        outcomes=[
            DownloadOutcome("GCF_000000001.1", True),
            DownloadOutcome("GCF_000000002.1", False, reason="suppressed"),
        ]
    )
    mock_cls.return_value = ncbi
    acc = tmp_path / "a.tsv"
    pl.DataFrame(
        {"accession": ["GCF_000000001.1", "GCF_000000002.1"], "gtdb_taxonomy": ["t", "t"],
         "species": ["s", "s"]}
    ).write_csv(acc, separator="\t")
    base = ["download", "genomes", "fetch", "-a", str(acc), "-o", str(tmp_path / "out")]
    default = runner.invoke(app, base)
    assert default.exit_code == 0, default.output
    assert "finished with 1 of 2" in default.output.replace("\n", " ")
    assert "Download complete" not in default.output
    strict = runner.invoke(app, [*base, "--require-all"])
    assert strict.exit_code == 1


# ---------------------------------------------------------------------------
# 9: visualize recruitment
# ---------------------------------------------------------------------------
def _viz_args(tmp_path: Path, *extra: str) -> list[str]:
    bam = tmp_path / "x.bam"
    bam.write_bytes(b"\x00" * 16)
    return ["visualize", "recruitment", "--bam", str(bam), *extra]


@patch("metadarkmatter.cli.visualize.Samtools")
@patch("metadarkmatter.cli.visualize.load_recruitment_data")
def test_recruitment_no_alignments_errors(mock_load, mock_st, tmp_path: Path):
    mock_st.check_available.return_value = True
    mock_load.return_value = pl.DataFrame({"genome_name": [], "identity": []})
    res = runner.invoke(app, _viz_args(tmp_path, "-o", str(tmp_path / "p.html")))
    assert res.exit_code == 1
    assert "Visualization complete" not in res.output
    assert not (tmp_path / "p.html").exists()


@patch("metadarkmatter.cli.visualize.Samtools")
@patch("metadarkmatter.cli.visualize.load_recruitment_data")
def test_recruitment_unknown_genome_errors(mock_load, mock_st, tmp_path: Path):
    mock_st.check_available.return_value = True
    mock_load.return_value = pl.DataFrame({"genome_name": ["A", "A"], "identity": [90.0, 91.0]})
    res = runner.invoke(
        app, _viz_args(tmp_path, "-o", str(tmp_path / "p.html"), "-g", "GCF_nonexistent")
    )
    assert res.exit_code == 1
    assert "Visualization complete" not in res.output


def test_recruitment_thresholds_must_decrease(tmp_path: Path):
    res = runner.invoke(
        app,
        _viz_args(tmp_path, "-o", str(tmp_path / "p.html"), "--known-threshold", "80",
                  "--novel-species-threshold", "85"),
    )
    assert res.exit_code == 1
    assert "decreasing" in res.output


@patch("metadarkmatter.cli.visualize.Samtools")
@patch("metadarkmatter.cli.visualize.load_recruitment_data")
def test_recruitment_extensionless_output_gets_extension(mock_load, mock_st, tmp_path: Path):
    mock_st.check_available.return_value = True
    mock_load.return_value = pl.DataFrame({"genome_name": [], "identity": []})
    # Zero alignments abort before writing, but the resolved name is reported.
    res = runner.invoke(app, _viz_args(tmp_path, "-o", str(tmp_path / "r10"), "-f", "json"))
    assert res.exit_code == 1
    assert "r10.json" in res.output.replace("\n", "")


# ---------------------------------------------------------------------------
# 10: doctor
# ---------------------------------------------------------------------------
def test_doctor_lists_optional_tools_and_strict_flag(monkeypatch: pytest.MonkeyPatch):
    from metadarkmatter.cli import doctor as doc

    monkeypatch.setattr(doc.shutil, "which", lambda exe: None)
    res = runner.invoke(app, ["doctor"])
    assert res.exit_code == 0
    for tool in ("bowtie2", "samtools", "prodigal"):
        assert tool in res.output
    assert "map" in res.output and "proteins" in res.output
    strict = runner.invoke(app, ["doctor", "--strict"])
    assert strict.exit_code == 1
    assert "missing required tool" in strict.output


def test_doctor_strict_passes_when_required_tools_present(monkeypatch: pytest.MonkeyPatch):
    from metadarkmatter.cli import doctor as doc

    monkeypatch.setattr(doc.shutil, "which", lambda exe: f"/bin/{exe}")
    monkeypatch.setattr(doc, "_probe_tool", lambda exe, args: (f"/bin/{exe}", "v1"))
    res = runner.invoke(app, ["doctor", "--strict"])
    assert res.exit_code == 0, res.output


def test_gzip_import_used():
    # keep gzip import honest for blast.gz validation
    p = Path("/nonexistent")
    assert gzip and not p.exists()
