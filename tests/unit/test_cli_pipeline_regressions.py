"""Regression tests for CLI / pipeline-orchestration bugs found in the 2026-11 audit.

Each test names the bug it guards. None of them require the external tools:
tool availability and subprocess calls are patched where needed.
"""

from __future__ import annotations

import gzip
import json
import subprocess
from pathlib import Path
from unittest.mock import patch

import polars as pl
import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.main import app
from metadarkmatter.core.pipeline import Pipeline, PipelineLayout, StepStatus
from metadarkmatter.core.pipeline.orchestrator import _db_present
from metadarkmatter.models.pipeline_config import Aligner, PipelineConfig, PipelineStep


@pytest.fixture(autouse=True)
def _clean_global_flags(monkeypatch: pytest.MonkeyPatch):
    """The global --dry-run sets MDM_DRY_RUN process-wide; isolate each test."""
    monkeypatch.delenv("MDM_DRY_RUN", raising=False)
    yield
    monkeypatch.delenv("MDM_DRY_RUN", raising=False)


class RecordingRunner:
    def __init__(self) -> None:
        self.calls: list[list[str]] = []

    def __call__(self, argv, dry_run: bool) -> int:
        self.calls.append(list(argv))
        return 0


def _genomes(tmp_path: Path) -> Path:
    g = tmp_path / "genomes"
    g.mkdir()
    (g / "GCF_000000001.1.fna").write_text(">c1\nACGT\n")
    return g


def _argv_after(argv: list[str], flag: str) -> str:
    return argv[argv.index(flag) + 1]


def _plan(cfg: PipelineConfig, **kw) -> dict[PipelineStep, list[str]]:
    pipeline = Pipeline(cfg, PipelineLayout(root=cfg.output_dir), **kw)
    return {s.step: s.argv for s in pipeline.plan()}


# ---------------------------------------------------------------------------
# Orchestrator: kraken2 paths, metadata, mmseqs2 query args, DB detection,
# --from, manifest timestamps
# ---------------------------------------------------------------------------


def _full_config(tmp_path: Path, **kw) -> PipelineConfig:
    reads = tmp_path / "s_R1.fastq.gz"
    reads.write_bytes(b"x")
    fields: dict = {
        "family": "f__X",
        "family_taxid": 34064,
        "reads_1": reads,
        "kraken_db": tmp_path / "kdb",
        "output_dir": tmp_path / "out",
    }
    fields.update(kw)
    return PipelineConfig(**fields)


def test_kraken_classify_gets_output_dir_and_extract_reads_real_kraken_file(tmp_path: Path):
    """Bug 2: --output for 'kraken2 classify' is a directory, not a .kraken path."""
    cfg = _full_config(tmp_path)
    layout = PipelineLayout(root=cfg.output_dir, family_taxid=cfg.family_taxid)
    plan = _plan(cfg)
    kc = plan[PipelineStep.KRAKEN_CLASSIFY]
    out_dir = Path(_argv_after(kc, "--output"))
    sample = _argv_after(kc, "--sample-name")
    assert out_dir == layout.root
    # classify writes {output}/{sample}.kraken; extract must read exactly that
    ke = plan[PipelineStep.KRAKEN_EXTRACT]
    assert Path(_argv_after(ke, "--kraken-output")) == out_dir / f"{sample}.kraken"


def test_extracted_read_names_match_kraken2_extract(tmp_path: Path):
    """Bug 1: layout must use the {sample}_taxid{taxid}_R1.fastq.gz name extract writes."""
    cfg = _full_config(tmp_path)
    plan = _plan(cfg)
    ke = plan[PipelineStep.KRAKEN_EXTRACT]
    sample = _argv_after(ke, "--sample-name")
    extraction_dir = Path(_argv_after(ke, "--output"))
    expected_r1 = extraction_dir / f"{sample}_taxid34064_R1.fastq.gz"
    align = plan[PipelineStep.ALIGN]
    assert Path(_argv_after(align, "--query")) == expected_r1


def test_metadata_passed_when_download_will_create_it(tmp_path: Path):
    """Bug 3: metadata is decided before download creates genome_metadata.tsv."""
    cfg = _full_config(tmp_path)
    plan = _plan(cfg)
    metadata = str(cfg.output_dir / "genome_metadata.tsv")
    for step in (PipelineStep.ANI, PipelineStep.CLASSIFY, PipelineStep.REPORT):
        assert _argv_after(plan[step], "--metadata") == metadata


def test_mmseqs2_single_end_raw_reads_uses_query(tmp_path: Path):
    """Bug 5: no --query-2 for single-end input."""
    cfg = _full_config(tmp_path, aligner=Aligner.MMSEQS2)
    align = _plan(cfg)[PipelineStep.ALIGN]
    assert "--query" in align
    assert "--query-1" not in align and "--query-2" not in align


def test_mmseqs2_paired_raw_reads_uses_query_1_and_2(tmp_path: Path):
    r2 = tmp_path / "s_R2.fastq.gz"
    r2.write_bytes(b"x")
    cfg = _full_config(tmp_path, aligner=Aligner.MMSEQS2, reads_2=r2)
    align = _plan(cfg)[PipelineStep.ALIGN]
    assert _argv_after(align, "--query-2").endswith("sample_taxid34064_R2.fastq.gz")


def test_mmseqs2_pre_extracted_single_and_paired(tmp_path: Path):
    reads = tmp_path / "r1.fq"
    reads.write_text("@r\nA\n+\nI\n")
    cfg = PipelineConfig(
        genomes_dir=_genomes(tmp_path), extracted_reads_1=reads,
        output_dir=tmp_path / "out", aligner=Aligner.MMSEQS2,
    )
    single = _plan(cfg)[PipelineStep.ALIGN]
    assert _argv_after(single, "--query") == str(reads)

    r2 = tmp_path / "r2.fq"
    paired = _plan(cfg, extracted_reads_2=r2)[PipelineStep.ALIGN]
    assert _argv_after(paired, "--query-1") == str(reads)
    assert _argv_after(paired, "--query-2") == str(r2)


def test_db_present_ignores_intermediate_files(tmp_path: Path):
    """Bug 11: concatenated FASTA / mapping TSV are not a finished database."""
    prefix = tmp_path / "blastdb" / "pangenome"
    prefix.parent.mkdir()
    (prefix.parent / "pangenome_pangenome.fasta").write_text(">a\nA\n")
    (prefix.parent / "pangenome_contig_mapping.tsv").write_text("x\n")
    assert not _db_present(prefix, Aligner.BLAST)
    assert not _db_present(prefix, Aligner.MMSEQS2)

    (prefix.parent / "pangenome.00.nhr").write_text("x")
    assert _db_present(prefix, Aligner.BLAST)
    (prefix.parent / "pangenome.dbtype").write_text("x")
    assert _db_present(prefix, Aligner.MMSEQS2)


def _compute_only(tmp_path: Path) -> tuple[PipelineConfig, PipelineLayout]:
    reads = tmp_path / "reads.fastq.gz"
    reads.write_bytes(b"x")
    cfg = PipelineConfig(
        genomes_dir=_genomes(tmp_path), extracted_reads_1=reads, output_dir=tmp_path / "out"
    )
    return cfg, PipelineLayout(root=cfg.output_dir)


def test_from_step_forces_rerun_of_selected_and_later_steps(tmp_path: Path):
    """Bug 13: --from must re-run the selected step even if outputs exist."""
    cfg, layout = _compute_only(tmp_path)
    db = layout.database(Aligner.BLAST)
    db.parent.mkdir(parents=True)
    Path(f"{db}.nhr").write_text("x")
    for p in (layout.alignment(Aligner.BLAST), layout.ani_matrix, layout.classifications,
              layout.report):
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("x")

    runner = RecordingRunner()
    manifest = Pipeline(cfg, layout, runner=runner).run(from_step=PipelineStep.CLASSIFY)
    statuses = {s.name: s.status for s in manifest.steps}
    assert statuses["ani"] == StepStatus.CACHED
    assert statuses["classify"] == StepStatus.COMPLETED
    assert statuses["report"] == StepStatus.COMPLETED
    assert len(runner.calls) == 2


def test_manifest_records_start_and_finish(tmp_path: Path):
    """Bug 17: started_at / finished_at were never set."""
    cfg, layout = _compute_only(tmp_path)
    Pipeline(cfg, layout, runner=RecordingRunner()).run()
    data = json.loads(layout.manifest.read_text())
    assert data["started_at"] and data["finished_at"]


# ---------------------------------------------------------------------------
# mdm run CLI: --config precedence, --extracted-reads-2
# ---------------------------------------------------------------------------


class _CapturingPipeline:
    captured: dict = {}

    def __init__(self, cfg, layout, **kw):
        type(self).captured = {"cfg": cfg, "layout": layout, **kw}

    def run(self, from_step=None, mdm_version=""):
        from metadarkmatter.core.pipeline import RunManifest

        return RunManifest(mdm_version="", output_dir="", aligner="", config={})


def test_run_config_yaml_values_not_overridden_by_cli_defaults(tmp_path: Path):
    """Bug 7: YAML aligner/threads/output_dir must be honoured."""
    reads = tmp_path / "r1.fq"
    reads.write_text("@r\nA\n+\nI\n")
    r2 = tmp_path / "r2.fq"
    r2.write_text("@r\nA\n+\nI\n")
    yaml_path = tmp_path / "run.yaml"
    yaml_path.write_text(
        f"genomes_dir: {_genomes(tmp_path)}\n"
        f"extracted_reads_1: {reads}\n"
        f"extracted_reads_2: {r2}\n"
        f"output_dir: {tmp_path / 'yaml_out'}\n"
        "aligner: mmseqs2\nthreads: 3\n"
    )
    with patch("metadarkmatter.cli.run.Pipeline", _CapturingPipeline):
        result = CliRunner().invoke(
            app, ["run", "pipeline", "--config", str(yaml_path), "--dry-run", "--quiet"]
        )
    assert result.exit_code == 0, result.output
    cfg = _CapturingPipeline.captured["cfg"]
    assert cfg.aligner == Aligner.MMSEQS2
    assert cfg.threads == 3
    assert cfg.output_dir == tmp_path / "yaml_out"
    assert _CapturingPipeline.captured["extracted_reads_2"] == r2

    # explicit CLI flag still wins
    with patch("metadarkmatter.cli.run.Pipeline", _CapturingPipeline):
        result = CliRunner().invoke(
            app,
            ["run", "pipeline", "--config", str(yaml_path), "--threads", "5", "--dry-run", "-q"],
        )
    assert result.exit_code == 0, result.output
    assert _CapturingPipeline.captured["cfg"].threads == 5


def test_run_without_config_requires_output_dir(tmp_path: Path):
    reads = tmp_path / "r1.fq"
    reads.write_text("@r\nA\n+\nI\n")
    result = CliRunner().invoke(
        app,
        ["run", "pipeline", "--genomes", str(_genomes(tmp_path)),
         "--extracted-reads-1", str(reads), "--dry-run"],
    )
    assert result.exit_code == 1
    assert "--output-dir" in result.output


# ---------------------------------------------------------------------------
# Classification runner / score CLI
# ---------------------------------------------------------------------------


def _blast_and_ani(tmp_path: Path) -> tuple[Path, Path]:
    genomes = ["GCF_000123456.1", "GCF_000789012.1"]
    pl.DataFrame(
        {"genome": genomes, genomes[0]: [100.0, 95.0], genomes[1]: [95.0, 100.0]}
    ).write_csv(tmp_path / "ani.csv")
    lines = [
        f"r{i}\t{g}|c1\t98.0\t150\t2\t0\t1\t150\t1\t150\t1e-50\t250\t150"
        for i in range(3)
        for g in genomes[:1]
    ]
    blast = tmp_path / "hits.tsv.gz"
    with gzip.open(blast, "wt") as f:
        f.write("\n".join(lines) + "\n")
    return blast, tmp_path / "ani.csv"


def test_streaming_with_id_mapping_is_an_error(tmp_path: Path):
    """Bug 6: --streaming silently ignored --id-mapping."""
    blast, ani = _blast_and_ani(tmp_path)
    mapping = tmp_path / "map.tsv"
    mapping.write_text("contig_id\tgenome_accession\nc1\tGCF_000123456.1\n")
    result = CliRunner().invoke(
        app,
        ["score", "classify", "--alignment", str(blast), "--ani", str(ani),
         "--id-mapping", str(mapping), "--streaming", "--output", str(tmp_path / "o.csv")],
    )
    assert result.exit_code != 0
    assert "--streaming" in result.output
    assert not (tmp_path / "o.csv").exists()


def test_global_dry_run_honoured_by_score_classify_and_counts_gz_lines(tmp_path: Path):
    """Bugs 8 + 17: global --dry-run ignored; .gz line count was of compressed bytes."""
    blast, ani = _blast_and_ani(tmp_path)
    out = tmp_path / "o.csv"
    result = CliRunner().invoke(
        app,
        ["--dry-run", "score", "classify", "--alignment", str(blast), "--ani", str(ani),
         "--output", str(out)],
    )
    assert result.exit_code == 0, result.output
    assert "DRY RUN" in result.output
    assert "Alignments: ~3" in result.output
    assert not out.exists()


def test_preset_keeps_explicit_flags_and_preset_only_fields():
    """Bug 18: --preset silently dropped CLI flags / preset-only fields."""
    from metadarkmatter.cli.score import THRESHOLD_PRESETS, _ClassifyOptions, build_classify_config

    opts: dict = {
        "alignment_mode": "nucleotide", "bitscore_threshold": 95.0,
        "min_alignment_length": 100, "min_alignment_fraction": 0.3, "max_evalue": 1e-10,
        "min_percent_identity": 70.0, "min_bitscore": 50.0, "min_read_length": 0,
        "min_query_coverage": 0.0, "coverage_weight_mode": "linear",
        "coverage_weight_strength": 0.5, "uncertainty_mode": "second",
        "single_hit_uncertainty_threshold": 10.0, "target_family": "f__Francisellaceae",
        "family_ratio_threshold": 0.9, "include_legacy_scores": True, "bayesian": False,
        "config_file": None, "preset": "literature-strict",
    }
    cfg, _ = build_classify_config(_ClassifyOptions(**opts))
    assert cfg.target_family == "f__Francisellaceae"
    assert cfg.max_evalue == 1e-10
    assert cfg.min_percent_identity == 70.0
    assert cfg.min_bitscore == 50.0
    assert cfg.family_ratio_threshold == 0.9
    assert cfg.include_legacy_scores is True

    opts["alignment_mode"] = "protein"
    cfg2, _ = build_classify_config(_ClassifyOptions(**opts))
    preset = THRESHOLD_PRESETS["literature-strict"]
    assert cfg2.alignment_mode == "protein"
    assert cfg2.model_dump(exclude={
        "alignment_mode", "target_family", "max_evalue", "min_percent_identity",
        "min_bitscore", "family_ratio_threshold", "include_legacy_scores",
    }) == preset.model_dump(exclude={
        "alignment_mode", "target_family", "max_evalue", "min_percent_identity",
        "min_bitscore", "family_ratio_threshold", "include_legacy_scores",
    })


def test_report_invalid_mode_message(tmp_path: Path):
    """Bug 17: bogus --report-mode printed an empty 'Error generating report:'."""
    csv = tmp_path / "c.csv"
    csv.write_text("read_id,taxonomic_call\nr1,Known Species\n")
    result = CliRunner().invoke(
        app,
        ["report", "generate", "--classifications", str(csv), "--output",
         str(tmp_path / "r.html"), "--report-mode", "bogus"],
    )
    assert result.exit_code == 2
    assert "Invalid --report-mode 'bogus'" in result.output
    assert "Error generating report" not in result.output


# ---------------------------------------------------------------------------
# blast CLI / wrapper
# ---------------------------------------------------------------------------


def test_blast_dry_run_global_flag_without_tools(tmp_path: Path):
    """Bug 8: global --dry-run ignored by blast makedb/align; missing tools blocked it."""
    q = tmp_path / "q.fa"
    q.write_text(">r\nACGT\n")
    with (
        patch("metadarkmatter.cli.blast._check_blastn_available",
              return_value=(False, "blastn (BLAST+)")),
        patch("metadarkmatter.cli.blast._check_makeblastdb_available",
              return_value=(False, "makeblastdb (BLAST+)")),
    ):
        r1 = CliRunner().invoke(
            app,
            ["--dry-run", "blast", "makedb", "--genomes", str(_genomes(tmp_path)),
             "--output", str(tmp_path / "db" / "p")],
        )
        r2 = CliRunner().invoke(
            app,
            ["--dry-run", "blast", "align", "--query", str(q), "--database",
             str(tmp_path / "db" / "p"), "--output", str(tmp_path / "x.tsv")],
        )
    assert r1.exit_code == 0, r1.output
    assert r2.exit_code == 0, r2.output
    assert not (tmp_path / "db").exists()
    assert not (tmp_path / "x.tsv").exists()


def _blast_align_dry_run(tmp_path: Path, *extra: str) -> str:
    q = tmp_path / "q.fa"
    q.write_text(">r\nACGT\n")
    db = tmp_path / "pan.v1"
    Path(f"{db}.nal").write_text("x")  # multi-volume alias, dotted prefix
    with patch("metadarkmatter.cli.blast._check_blastn_available", return_value=(True, "")), \
         patch("metadarkmatter.external.blast.BlastN.get_executable",
               return_value=Path("/usr/bin/blastn")):
        result = CliRunner().invoke(
            app,
            ["blast", "align", "--query", str(q), "--database", str(db), "--dry-run", *extra],
            env={"COLUMNS": "400"},
        )
    assert result.exit_code == 0, result.output
    return result.output


def test_blast_align_compress_follows_extension(tmp_path: Path):
    """Bug 12: '--output x.tsv' silently became x.tsv.gz."""
    plain = _blast_align_dry_run(tmp_path, "--output", str(tmp_path / "x.tsv"))
    assert "gzip" not in plain
    gz = _blast_align_dry_run(tmp_path, "--output", str(tmp_path / "x.tsv.gz"))
    assert "gzip" in gz


def test_blast_align_task_default_word_size(tmp_path: Path):
    """Coordinator E: dc-megablast needs word size 11/12, not the blastn default 7."""
    out = _blast_align_dry_run(
        tmp_path, "--output", str(tmp_path / "x.tsv"), "--task", "dc-megablast"
    )
    assert "-word_size 11" in out
    out = _blast_align_dry_run(tmp_path, "--output", str(tmp_path / "x.tsv"))
    assert "-word_size 7" in out


def test_blast_db_exists_accepts_volumes_and_dotted_prefix(tmp_path: Path):
    """Bug 11: only <db>.nhr was accepted, and with_suffix corrupted dotted prefixes."""
    from metadarkmatter.cli.blast import blast_db_exists

    prefix = tmp_path / "pangenome.v2"
    assert not blast_db_exists(prefix)
    Path(f"{prefix}.00.nhr").write_text("x")
    assert blast_db_exists(prefix)
    other = tmp_path / "db.v3"
    Path(f"{other}.nal").write_text("x")
    assert blast_db_exists(other)
    (tmp_path / "db.nhr").write_text("x")  # with_suffix('.nhr') of db.v4 would hit this
    assert not blast_db_exists(tmp_path / "db.v4")


def test_parse_seqids_disabled_for_long_ids(tmp_path: Path):
    """Coordinator D: makeblastdb -parse_seqids fails for IDs > 50 chars."""
    from metadarkmatter.cli.blast import _fasta_ids_fit_parse_seqids

    short = tmp_path / "s.fa"
    short.write_text(">GCF_000195955.2|NZ_CP000001.1 desc\nACGT\n")
    long = tmp_path / "l.fa"
    long.write_text(">GCA_000123456.1|NODE_12345_length_67890_cov_12.345678_extra\nACGT\n")
    assert _fasta_ids_fit_parse_seqids(short)
    assert not _fasta_ids_fit_parse_seqids(long)


def test_gzipped_fasta_query_not_treated_as_fastq(tmp_path: Path):
    """Bug 9: gz FASTA decompressed to a .fastq temp lost records in the fallback."""
    from metadarkmatter.cli.blast import _sniff_sequence_suffix
    from metadarkmatter.external.blast import BlastN

    fa_gz = tmp_path / "reads.gz"
    with gzip.open(fa_gz, "wt") as f:
        f.write(">r1\nACGT\nACGT\n>r2\nGGGG\n>r3\nTTTT\n")
    assert _sniff_sequence_suffix(fa_gz) == ".fasta"
    fq = tmp_path / "x.fa"
    fq.write_text("@r1\nACGT\n+\nIIII\n")
    assert _sniff_sequence_suffix(fq) == ".fastq"

    # Content wins over a misleading extension in the wrapper too.
    mislabeled = tmp_path / "tmp123.fastq"
    mislabeled.write_text(">r1\nACGT\n>r2\nGGGG\n")
    assert BlastN._is_fastq_file(mislabeled) is False
    assert BlastN._is_fastq_file(fq) is True


# ---------------------------------------------------------------------------
# mmseqs2 / blastx helpers
# ---------------------------------------------------------------------------


def test_paired_concatenation_makes_mate_ids_distinct(tmp_path: Path):
    """Bug 15 + coordinator C: identical mate IDs merged; gz detected per file;
    a missing trailing newline merged records."""
    from rich.console import Console

    from metadarkmatter.cli.mmseqs2 import _concatenate_paired_reads

    r1 = tmp_path / "R1.fastq.gz"
    with gzip.open(r1, "wt") as f:
        f.write("@SRR1.1 len=4\nACGT\n+\nIIII\n@SRR1.2\nACGT\n+\nIIII")  # no final newline
    r2 = tmp_path / "R2.fastq"  # plain, while R1 is gzipped
    r2.write_text("@SRR1.1 len=4\nTTTT\n+\nIIII\n@SRR1.2/2\nTTTT\n+\nIIII\n")
    out = tmp_path / "cat.fastq"
    n = _concatenate_paired_reads(r1, r2, out, Console(quiet=True), quiet=True)
    assert n == 4
    headers = [ln for i, ln in enumerate(out.read_text().splitlines()) if i % 4 == 0]
    assert headers == ["@SRR1.1/1 len=4", "@SRR1.2/1", "@SRR1.1/2 len=4", "@SRR1.2/2"]


def test_mmseqs2_translated_search_uses_nucleotide_query_db(tmp_path: Path):
    """Coordinator B: --search-type 2 built an amino-acid query DB (dbtype 1)."""
    from metadarkmatter.external.mmseqs2 import MMseqs2

    q = tmp_path / "q.fa"
    q.write_text(">r\nACGT\n")
    seen: list[list[str]] = []

    class _Stop(Exception):
        pass

    def fake_run(cmd, **kw):
        seen.append(cmd)
        raise _Stop

    with patch.object(MMseqs2, "get_executable", return_value=Path("/usr/bin/mmseqs")), \
         patch("metadarkmatter.external.mmseqs2.subprocess.run", side_effect=fake_run), \
         pytest.raises(_Stop):
        MMseqs2().search_multistep(
            query=q, database=tmp_path / "db", output=tmp_path / "o.tsv",
            tmp_dir=tmp_path / "tmp", search_type=2,
        )
    createdb = seen[0]
    assert createdb[1] == "createdb"
    assert createdb[createdb.index("--dbtype") + 1] == "2"


def test_blastx_protein_accession_from_faa_gz(tmp_path: Path):
    """Bug 14: .stem gave 'GCF_x.1.faa' for .faa.gz files; newline guard (C)."""
    from metadarkmatter.cli.blastx import _concatenate_proteins, _protein_file_accession

    assert _protein_file_accession(Path("GCF_000195955.2.faa.gz")) == "GCF_000195955.2"
    assert _protein_file_accession(Path("GCF_000195955.2_ASM1_protein.faa")) == "GCF_000195955.2"

    d = tmp_path / "p"
    d.mkdir()
    with gzip.open(d / "GCF_000000001.1.faa.gz", "wt") as f:
        f.write(">p1\nMKV")  # no trailing newline
    (d / "GCF_000000002.1.faa").write_text(">p2\nMAA\n")
    out = tmp_path / "all.faa"
    _concatenate_proteins(d, out, pattern="*.faa*")
    assert out.read_text().splitlines() == [
        ">GCF_000000001.1|p1", "MKV", ">GCF_000000002.1|p2", "MAA",
    ]


# ---------------------------------------------------------------------------
# ani compute
# ---------------------------------------------------------------------------


def _ani_inputs(tmp_path: Path) -> tuple[Path, Path]:
    g = tmp_path / "genomes"
    g.mkdir()
    for acc in ("GCF_000123456.1", "GCF_000789012.1", "GCA_000111222.1"):
        (g / f"{acc}.fna").write_text(">c\nACGT\n")
    md = tmp_path / "genome_metadata.tsv"
    md.write_text(
        "accession\tspecies\tgenus\tfamily\trepresentative\tgtdb_taxonomy\n"
        "GCF_000123456.1\ts1\tg\tF\tGCF_000123456.1\tx\n"
        "GCF_000789012.1\ts1\tg\tF\tGCF_000123456.1\tx\n"
        "GCA_000111222.1\ts2\tg\tF\tGCA_000111222.1\tx\n"
    )
    return g, md


def test_fastani_representatives_only_list(tmp_path: Path):
    """Bug 4: fastANI list file ignored --representatives-only."""
    from metadarkmatter.external import FastANI

    g, md = _ani_inputs(tmp_path)
    captured: dict[str, str] = {}

    def fake(self, query_list, reference_list, output, **kw):
        captured["list"] = Path(query_list).read_text()
        names = captured["list"].split()
        with open(output, "w") as f:
            for a in names:
                for b in names:
                    f.write(f"{a}\t{b}\t97.0\t10\t10\n")

        class R:
            command_string = "fastANI"
            stderr = ""

        return R()

    with patch.object(FastANI, "check_available", return_value=True), \
         patch.object(FastANI, "run_or_raise", fake):
        result = CliRunner().invoke(
            app,
            ["ani", "compute", "--genomes", str(g), "--output", str(tmp_path / "a.csv"),
             "--metadata", str(md), "--representatives-only", "--backend", "fastani"],
        )
    assert result.exit_code == 0, result.output
    listed = sorted(Path(p).name for p in captured["list"].split())
    assert listed == ["GCA_000111222.1.fna", "GCF_000123456.1.fna"]


def test_fastani_batch_failure_is_reported_not_traceback(tmp_path: Path):
    """Coordinator F: batched fastANI RuntimeError escaped as a traceback."""
    from metadarkmatter.external import FastANI

    g, _ = _ani_inputs(tmp_path)

    def boom(self, **kw):
        raise RuntimeError("fastANI batch failed (exit 1): bad genome")

    with patch.object(FastANI, "check_available", return_value=True), \
         patch.object(FastANI, "compute_parallel_batches", boom):
        result = CliRunner().invoke(
            app,
            ["ani", "compute", "--genomes", str(g), "--output", str(tmp_path / "a.csv"),
             "--backend", "fastani", "--fastani-batches", "2"],
        )
    assert result.exit_code == 1
    assert "fastANI failed" in result.output
    assert not isinstance(result.exception, RuntimeError)


# ---------------------------------------------------------------------------
# external/base timeouts
# ---------------------------------------------------------------------------


def test_no_default_timeout_and_timeout_is_execution_error():
    """Coordinator A: 1 h default timeout killed long runs with an uncaught error."""
    from metadarkmatter.external.base import ExternalTool, ToolExecutionError, ToolTimeoutError

    assert ExternalTool.DEFAULT_TIMEOUT is None
    err = ToolTimeoutError("blastn", 10.0, ["blastn", "-query", "x"])
    assert isinstance(err, ToolExecutionError)
    assert "timed out" in err.message


def test_run_passes_no_timeout_by_default(tmp_path: Path):
    from metadarkmatter.external.blast import MakeBlastDb

    seen: dict = {}

    def fake_run(cmd, **kw):
        seen.update(kw)
        return subprocess.CompletedProcess(cmd, 0, "", "")

    with patch.object(MakeBlastDb, "get_executable", return_value=Path("/usr/bin/makeblastdb")), \
         patch("metadarkmatter.external.base.subprocess.run", side_effect=fake_run):
        MakeBlastDb().run(input_fasta=tmp_path / "a.fa", output_db=tmp_path / "db")
    assert seen["timeout"] is None
