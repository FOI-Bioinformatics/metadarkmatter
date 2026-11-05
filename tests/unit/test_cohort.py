"""Tests for ``mdm run cohort`` (sample sheet, orchestration, summary, CLI)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.main import app
from metadarkmatter.core.pipeline import (
    Cohort,
    CohortError,
    CohortLayout,
    CohortSample,
    read_samples_tsv,
)
from metadarkmatter.models.pipeline_config import Aligner, PipelineStep

runner = CliRunner()


@pytest.fixture(autouse=True)
def _clean_flags(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.delenv("MDM_DRY_RUN", raising=False)


def _touch(path: Path, text: str = "x") -> Path:
    if not path.exists():
        path.write_text(text)
    return path


def _sheet(tmp_path: Path, rows: list[str], header: str = "sample\treads_1\treads_2") -> Path:
    p = tmp_path / "samples.tsv"
    p.write_text("\n".join([header, *rows]) + "\n")
    return p


# ---------------------------------------------------------------- sample sheet


def test_sheet_parses_paired_and_extracted(tmp_path: Path):
    _touch(tmp_path / "a1.fq")
    _touch(tmp_path / "a2.fq")
    _touch(tmp_path / "b1.fq")
    sheet = tmp_path / "s.tsv"
    sheet.write_text(
        "sample\treads_1\treads_2\textracted_reads_1\textracted_reads_2\n"
        "# a comment\n"
        "A\ta1.fq\ta2.fq\t\t\n"
        "B\t\t\tb1.fq\t\n"
    )
    a, b = read_samples_tsv(sheet)
    assert (a.name, a.reads_1, a.reads_2) == ("A", tmp_path / "a1.fq", tmp_path / "a2.fq")
    assert b.extracted_reads_1 == tmp_path / "b1.fq" and b.reads_1 is None


@pytest.mark.parametrize(
    ("rows", "fragment"),
    [
        (["A\ta1.fq\t", "A\ta1.fq\t"], "duplicate sample name"),
        (["A\tmissing.fq\t"], "file not found"),
        (["bad name\ta1.fq\t"], "sample names may only"),
        (["x_classifications\ta1.fq\t"], "must not contain"),
        (["A\t\t"], "provide reads_1 or extracted_reads_1"),
    ],
)
def test_sheet_validation_errors(tmp_path: Path, rows: list[str], fragment: str):
    _touch(tmp_path / "a1.fq")
    with pytest.raises(CohortError, match=fragment):
        read_samples_tsv(_sheet(tmp_path, rows))


def test_sheet_requires_sample_column_and_known_columns(tmp_path: Path):
    with pytest.raises(CohortError, match="'sample' column"):
        read_samples_tsv(_sheet(tmp_path, ["x"], header="name"))
    with pytest.raises(CohortError, match="Unknown column"):
        read_samples_tsv(_sheet(tmp_path, ["A\tx"], header="sample\tfoo"))


def test_sheet_extracted_and_raw_mutually_exclusive(tmp_path: Path):
    _touch(tmp_path / "a1.fq")
    sheet = _sheet(tmp_path, ["A\ta1.fq\t\ta1.fq"], header="sample\treads_1\treads_2\textracted_reads_1")
    with pytest.raises(CohortError, match="not both"):
        read_samples_tsv(sheet)


# ---------------------------------------------------------------- orchestration


class FakeRunner:
    """Simulates the tools by creating each step's expected output."""

    def __init__(self, fail_sample: str | None = None) -> None:
        self.calls: list[list[str]] = []
        self.fail_sample = fail_sample

    def __call__(self, argv, dry_run: bool) -> int:
        argv = list(argv)
        self.calls.append(argv)
        key = argv[3:5]
        out = Path(argv[argv.index("--output") + 1]) if "--output" in argv else None
        if self.fail_sample and key == ["blast", "align"] and f"/{self.fail_sample}" in str(out):
            return 1
        if key == ["blast", "makedb"] and out is not None:
            out.parent.mkdir(parents=True, exist_ok=True)
            Path(f"{out}.nhr").write_text("db")
        elif key == ["ani", "compute"] and out is not None:
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text("ani")
        elif key == ["blast", "align"] and out is not None:
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text("aln")
        elif key == ["score", "classify"] and out is not None:
            out.parent.mkdir(parents=True, exist_ok=True)
            sample = out.parent.name
            calls = ["Known Species"] * 3 + (["Novel Species"] if sample == "A" else ["Novel Genus"] * 2)
            out.write_text("read_id,taxonomic_call\n" + "".join(f"r{i},{c}\n" for i, c in enumerate(calls)))
        elif key == ["report", "generate"] and out is not None:
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text("<html/>")
        elif key == ["report", "multi"] and out is not None:
            out.write_text("<html/>")
        return 0


def _cohort(tmp_path: Path, runner_: FakeRunner, **kw) -> tuple[Cohort, CohortLayout]:
    genomes = tmp_path / "genomes"
    genomes.mkdir(exist_ok=True)
    _touch(genomes / "GCF_1.fna", ">c\nACGT\n")
    samples = []
    for name in ("A", "B"):
        r1 = _touch(tmp_path / f"{name}_1.fq")
        r2 = _touch(tmp_path / f"{name}_2.fq") if name == "A" else None
        samples.append(CohortSample(name, extracted_reads_1=r1, extracted_reads_2=r2))
    layout = CohortLayout(root=tmp_path / "cohort")
    cohort = Cohort(
        samples, {"genomes_dir": genomes, "threads": 2}, layout, runner=runner_, **kw
    )
    return cohort, layout


def test_reference_built_once_and_samples_share_it(tmp_path: Path):
    fake = FakeRunner()
    cohort, layout = _cohort(tmp_path, fake)
    result = cohort.run(mdm_version="t")
    assert result.ok

    keys = [tuple(c[3:5]) for c in fake.calls]
    assert keys.count(("blast", "makedb")) == 1
    assert keys.count(("ani", "compute")) == 1
    assert keys.count(("blast", "align")) == 2
    assert keys[-1] == ("report", "multi")

    aligns = [c for c in fake.calls if c[3:5] == ["blast", "align"]]
    # sample A is paired, B single-end; both use the shared reference database
    assert "--query-2" in aligns[0] and "--query" in aligns[1]
    for c in aligns:
        db = c[c.index("--database") + 1]
        assert db == str(layout.reference_dir / "blastdb" / "pangenome")
    # file stems follow the sample name
    assert aligns[0][aligns[0].index("--output") + 1].endswith("samples/A/A.blast.tsv.gz")
    classify = [c for c in fake.calls if c[3:5] == ["score", "classify"]]
    for c in classify:
        assert c[c.index("--ani") + 1] == str(layout.reference_dir / "ani_matrix.csv")

    assert (layout.reference_dir / "run_manifest.json").exists()
    assert (layout.sample_dir("A") / "run_manifest.json").exists()
    assert layout.report.exists()
    assert (layout.multi_inputs_dir / "A_classifications.csv").resolve() == (
        layout.sample_dir("A") / "classifications.csv"
    ).resolve()


def test_summary_tsv_contents(tmp_path: Path):
    cohort, layout = _cohort(tmp_path, FakeRunner())
    cohort.run()
    lines = [ln.split("\t") for ln in layout.summary_tsv.read_text().splitlines()]
    assert lines[0] == [
        "sample", "status", "total_reads", "Known Species", "Novel Genus", "Novel Species",
        "novel_fraction", "error",
    ]
    assert lines[1] == ["A", "completed", "4", "3", "0", "1", "0.250000", ""]
    assert lines[2] == ["B", "completed", "5", "3", "2", "0", "0.400000", ""]


def test_failed_sample_does_not_stop_others(tmp_path: Path):
    fake = FakeRunner(fail_sample="A")
    cohort, layout = _cohort(tmp_path, fake)
    result = cohort.run()
    assert not result.ok
    assert [s.name for s in result.failed] == ["A"]
    rows = {ln.split("\t")[0]: ln.split("\t") for ln in layout.summary_tsv.read_text().splitlines()[1:]}
    assert rows["A"][1] == "failed" and "align" in rows["A"][-1]
    assert rows["B"][1] == "completed"
    # the cohort report is still built from the surviving sample
    assert not (layout.multi_inputs_dir / "A_classifications.csv").exists()
    assert (layout.multi_inputs_dir / "B_classifications.csv").is_symlink()
    manifest = json.loads(layout.manifest.read_text())
    assert [s["status"] for s in manifest["samples"]] == ["failed", "completed"]


def test_reference_failure_is_fatal(tmp_path: Path):
    class Boom(FakeRunner):
        def __call__(self, argv, dry_run):
            return 1 if list(argv)[3:5] == ["blast", "makedb"] else super().__call__(argv, dry_run)

    cohort, layout = _cohort(tmp_path, Boom())
    with pytest.raises(CohortError, match="reference"):
        cohort.run()
    assert not layout.samples_dir.exists()


def test_rerun_is_cached_per_sample(tmp_path: Path):
    fake = FakeRunner()
    cohort, _ = _cohort(tmp_path, fake)
    cohort.run()
    n_first = len(fake.calls)
    fake2 = FakeRunner()
    cohort2, _ = _cohort(tmp_path, fake2)
    cohort2.run()
    # only the multi-sample report step is re-issued
    assert [c[3:5] for c in fake2.calls] == [["report", "multi"]]
    assert n_first > 1


def test_dry_run_plans_everything_without_writing(tmp_path: Path):
    fake = FakeRunner()
    logs: list[str] = []
    cohort, layout = _cohort(tmp_path, fake, dry_run=True, log=logs.append)
    result = cohort.run()
    assert fake.calls == [] and not layout.root.exists()
    assert result.ok
    text = "\n".join(logs)
    assert "[reference] [dry-run] makedb" in text and "[reference] [dry-run] ani" in text
    for name in ("A", "B"):
        assert f"[{name}] [dry-run] align" in text and f"[{name}] [dry-run] classify" in text
    assert "report multi" in text


def test_from_step_applies_to_samples(tmp_path: Path):
    cohort, _ = _cohort(tmp_path, FakeRunner())
    cohort.run()
    fake = FakeRunner()
    again, _ = _cohort(tmp_path, fake)
    again.run(from_step=PipelineStep.CLASSIFY)
    keys = [c[3:5] for c in fake.calls]
    assert keys.count(["score", "classify"]) == 2 and ["blast", "align"] not in keys


def test_mmseqs2_cohort_uses_shared_mmseqs_db(tmp_path: Path):
    fake = FakeRunner()
    cohort, layout = _cohort(tmp_path, fake, aligner=Aligner.MMSEQS2, dry_run=True)
    cohort.run()
    assert cohort.reference_pipeline().layout.database(Aligner.MMSEQS2) == (
        layout.reference_dir / "mmseqs_db" / "pangenome"
    )


# ---------------------------------------------------------------- CLI


def test_cli_cohort_dry_run(tmp_path: Path):
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    (genomes / "GCF_1.fna").write_text(">c\nACGT\n")
    _touch(tmp_path / "a.fq")
    _touch(tmp_path / "b.fq")
    sheet = _sheet(tmp_path, ["A\ta.fq", "B\tb.fq"], header="sample\textracted_reads_1")
    out = tmp_path / "cohort"
    res = runner.invoke(
        app,
        ["run", "cohort", "--samples", str(sheet), "--genomes", str(genomes),
         "--output-dir", str(out), "--dry-run"],
    )
    assert res.exit_code == 0, res.output
    assert "planned" in res.output and "[A]" in res.output and "[B]" in res.output
    assert not out.exists()


def test_cli_cohort_invalid_sheet_exits_nonzero(tmp_path: Path):
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    sheet = _sheet(tmp_path, ["A\tnope.fq"], header="sample\textracted_reads_1")
    res = runner.invoke(
        app,
        ["run", "cohort", "--samples", str(sheet), "--genomes", str(genomes),
         "--output-dir", str(tmp_path / "o")],
    )
    assert res.exit_code != 0
    assert "file not found" in res.output.replace("\n", " ")


def test_cli_cohort_raw_reads_need_taxid_and_kraken(tmp_path: Path):
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    _touch(tmp_path / "a.fq")
    sheet = _sheet(tmp_path, ["A\ta.fq"], header="sample\treads_1")
    res = runner.invoke(
        app,
        ["run", "cohort", "--samples", str(sheet), "--genomes", str(genomes),
         "--output-dir", str(tmp_path / "o"), "--dry-run"],
    )
    assert res.exit_code != 0
    assert "family-taxid" in res.output.replace("\n", " ")
