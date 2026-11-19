# NCBI Subspecies Labels Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Label genomes inside selected GTDB species with their NCBI species ("subspecies"), report the closest subspecies per Known Species read, and summarise subspecies counts in the HTML report.

**Architecture:** A new module `core/subspecies.py` holds all subspecies logic (name parsing, panel selection, per-read assignment, report summary). `download genomes list` builds the panel via NCBI Datasets and writes a `subspecies` metadata column. `VectorizedClassifier` collects the genomes tied at the top weighted bitscore and, given a subspecies map, adds two columns; because streaming calls `classify_file` per partition, both paths get the columns from one change. The report reads the columns plus metadata.

**Tech Stack:** Python 3.11, polars, typer, pydantic, pytest; NCBI Datasets CLI (`datasets`).

**Spec:** `docs/superpowers/specs/2026-11-19-ncbi-subspecies-design.md`

## Global Constraints

- Metadata without a `subspecies` column (or with it all empty) must produce byte-identical classification CSVs and reports to the current code.
- All panel genomes keep `species` = GTDB species name and `representative` = GTDB species representative.
- Only reads with `taxonomic_call == "Known Species"` get non-null `closest_subspecies` / `subspecies_candidates`.
- Unresolved label string: `"unresolved"` (exact). Candidates separator: `";"`, sorted ascending.
- Panel defaults: `per_subspecies=2`, `per_subspecies_large=3`, large threshold: more than 50 genomes.
- Ranking order: RefSeq category `"reference genome"` first; assembly level `Complete Genome` < `Chromosome` < `Scaffold` < `Contig` < other; `GCF_` before `GCA_`; contig N50 descending; accession ascending.
- Determinism: sort every `group_by`/`unique` output (repo convention, CLAUDE.md).
- Modest scientific language in docs, help text and report notes.
- Run tests with: `uv run --frozen pytest <path> -v`

## Review Focus

1. Organism names like "Brucella sp. 09RB8471", "Candidatus Liberibacter asiaticus str. X", "[Brucella] something", or a single word: only real binomials become subspecies; "sp." and malformed names are skipped. Tests in Task 2.
2. `datasets` not installed or returning a non-zero exit: `download genomes list --ncbi-subspecies` exits 1 with the install hint, and writes no partial metadata. Test in Task 4.
3. `--ncbi-subspecies` naming a species absent from the queried taxon (e.g. typo, or outside the family): exit 1 with a message naming it. Test in Task 4.
4. Tied top hits that include genomes without a subspecies (other GTDB species or unlabelled): those are ignored; a read with no labelled tied hit gets null. Test in Task 3.
5. Streaming and in-memory runs on the same alignment produce identical `closest_subspecies` / `subspecies_candidates`. Test in Task 5.

---

### Task 1: NCBI assembly summary lookup

**Files:**
- Modify: `src/metadarkmatter/external/ncbi_datasets.py` (add dataclass + function near the top-level helpers, and a method on `NCBIDatasets`)
- Test: `tests/unit/test_ncbi_datasets.py`

**Interfaces:**
- Produces:
  - `@dataclass(frozen=True) class AssemblySummary: accession: str; organism_name: str; assembly_level: str; refseq_category: str; contig_n50: int`
  - `parse_assembly_summary_jsonl(text: str) -> list[AssemblySummary]` reading `datasets summary genome accession --as-json-lines` output. JSON paths: `accession`, `organism.organism_name`, `assembly_info.assembly_level`, `assembly_info.refseq_category` (missing -> `""`), `assembly_stats.contig_n50` (missing -> `0`). Blank lines skipped.
  - `NCBIDatasets.fetch_assembly_summaries(self, accessions: list[str]) -> list[AssemblySummary]`: writes accessions to a temp file, runs `datasets summary genome accession --inputfile <tmp> --as-json-lines` via the base-class run mechanism, parses stdout. Raises `ToolNotFoundError` (`external/base.py:123`) when the executable is missing and `ToolExecutionError` (`external/base.py:138`) on non-zero exit.

- [ ] **Step 1: Write failing tests** `test_parse_assembly_summary_jsonl_fields`, `test_parse_assembly_summary_jsonl_missing_optional_fields` (no refseq_category, no contig_n50 -> `""`, `0`), `test_fetch_assembly_summaries_builds_command` (patch `subprocess.run`; assert argv contains `["summary", "genome", "accession", "--inputfile"]` and `"--as-json-lines"`).
- [ ] **Step 2: Run** `uv run --frozen pytest tests/unit/test_ncbi_datasets.py -v` - expect the three new tests to FAIL (ImportError / AttributeError).
- [ ] **Step 3: Implement** the dataclass, parser and method.
- [ ] **Step 4: Run** the same command - expect PASS.
- [ ] **Step 5: Commit** `feat(ncbi): fetch assembly summaries via datasets CLI`

### Task 2: Subspecies name parsing and panel selection

**Files:**
- Create: `src/metadarkmatter/core/subspecies.py`
- Test: `tests/unit/test_subspecies.py`

**Interfaces:**
- Consumes: `AssemblySummary` (Task 1).
- Produces:
  - `UNRESOLVED: Final = "unresolved"`
  - `ncbi_species_name(organism_name: str) -> str | None`: first two whitespace tokens after stripping a leading `"Candidatus "`; returns None when fewer than two tokens, the second token is `"sp."`, or the first token is not capitalised alphabetic (rejects `[Brucella]`).
  - `@dataclass(frozen=True) class PanelEntry: accession: str; gtdb_species: str; subspecies: str; ncbi_organism_name: str; assembly_level: str; refseq_category: str; contig_n50: int; rank: int; n_genomes_in_subspecies: int`
  - `select_subspecies_panel(summaries: list[AssemblySummary], gtdb_species: str, gtdb_representative: str, *, per_subspecies: int = 2, per_subspecies_large: int = 3, large_threshold: int = 50) -> list[PanelEntry]`: groups by `ncbi_species_name`, ranks per Global Constraints, keeps top N (large count when group size > `large_threshold`). The GTDB representative is always included (with its own subspecies, rank as ranked) even if outside the top N. Output sorted by (subspecies, rank).
  - `write_panel_tsv(entries: list[PanelEntry], path: Path) -> None`: columns in dataclass field order, tab-separated.

- [ ] **Step 1: Write failing tests**
  - `test_ncbi_species_name_cases`: `"Brucella abortus 2308"` -> `"Brucella abortus"`; `"Brucella sp. 09RB8471"` -> None; `"Candidatus Liberibacter asiaticus str. psy62"` -> `"Liberibacter asiaticus"`; `"[Brucella] x"` -> None; `"Brucella"` -> None.
  - `test_select_panel_ranking`: four abortus genomes (reference genome Scaffold GCF; Complete GCA N50 9e6; Complete GCF N50 2e6; Complete GCF N50 3e6) with `per_subspecies=2` -> accessions `[reference Scaffold, Complete GCF N50 3e6]`.
  - `test_select_panel_large_group_keeps_three`: 51 genomes of one species -> 3 entries; 50 -> 2.
  - `test_select_panel_always_keeps_representative`: representative ranked last still appears.
  - `test_select_panel_skips_sp`: `"Brucella sp."` genomes produce no entries.
  - `test_write_panel_tsv_header`: header equals the dataclass field names joined by tab.
- [ ] **Step 2: Run** `uv run --frozen pytest tests/unit/test_subspecies.py -v` - expect FAIL (module missing).
- [ ] **Step 3: Implement** the functions in `core/subspecies.py`.
- [ ] **Step 4: Run** - expect PASS.
- [ ] **Step 5: Commit** `feat(subspecies): NCBI species parsing and panel selection`

### Task 3: Metadata subspecies column and per-read assignment

**Files:**
- Modify: `src/metadarkmatter/models/genomes.py` (`GenomeAccession`, `AccessionList.from_tsv`, `AccessionList.to_metadata_tsv`)
- Modify: `src/metadarkmatter/core/metadata.py` (`GenomeMetadata.__init__`, new accessors, `join_classifications`)
- Modify: `src/metadarkmatter/core/subspecies.py`
- Test: `tests/unit/test_metadata.py`, `tests/unit/test_genome_models.py`, `tests/unit/test_subspecies.py`

**Interfaces:**
- Produces:
  - `GenomeAccession.subspecies: str = ""`
  - `AccessionList.to_metadata_tsv` appends a final `subspecies` column only when at least one accession has a non-empty subspecies; `from_tsv` reads it when present.
  - `GenomeMetadata.get_subspecies(accession: str) -> str | None`, `GenomeMetadata.has_subspecies -> bool` (any non-empty), `GenomeMetadata.subspecies_map() -> dict[str, str]` (non-empty only).
  - `join_classifications` adds `subspecies_reference` = subspecies of `best_match_genome` only when `has_subspecies` (name avoids clashing with the per-read `closest_subspecies`).
  - `assign_subspecies(df: pl.DataFrame, subspecies_map: dict[str, str]) -> pl.DataFrame` in `core/subspecies.py`. Input has `taxonomic_call` and `top_hit_genomes` (`list[str]`). Output drops `top_hit_genomes`, adds `closest_subspecies` (Utf8) and `subspecies_candidates` (Utf8). Rule: labels = sorted unique `subspecies_map[g]` for g in top_hit_genomes that are in the map; if call != "Known Species" or no labels -> both null; one label -> that label; several -> `UNRESOLVED`. Candidates = labels joined by `";"`.

- [ ] **Step 1: Write failing tests**
  - `test_metadata_tsv_without_subspecies_unchanged` (genome models): no subspecies set -> header is exactly `accession\tspecies\tgenus\tfamily\trepresentative\tgtdb_taxonomy`.
  - `test_metadata_tsv_with_subspecies_roundtrip`: written column round-trips through `from_tsv` and `GenomeMetadata.from_file`.
  - `test_genome_metadata_subspecies_accessors`: empty cells -> `get_subspecies` None, absent from `subspecies_map`, `has_subspecies` False when all empty.
  - `test_assign_subspecies_rules` with rows: Known Species, tied `[A(abortus)]` -> `"Brucella abortus"`, candidates `"Brucella abortus"`; Known Species, `[A(abortus), S(suis)]` -> `"unresolved"`, `"Brucella abortus;Brucella suis"`; Known Species, `[A(abortus), X(unlabelled)]` -> `"Brucella abortus"`; Known Species, `[X]` -> null/null; Novel Species, `[A]` -> null/null. Assert `top_hit_genomes` not in output columns.
- [ ] **Step 2: Run** `uv run --frozen pytest tests/unit/test_metadata.py tests/unit/test_genome_models.py tests/unit/test_subspecies.py -v` - expect new tests FAIL.
- [ ] **Step 3: Implement** model, metadata and `assign_subspecies` changes.
- [ ] **Step 4: Run** the same command plus `tests/unit/test_representative_mapping.py` - expect PASS.
- [ ] **Step 5: Commit** `feat(metadata): optional subspecies column and per-read assignment`

### Task 4: `download genomes list --ncbi-subspecies`

**Files:**
- Modify: `src/metadarkmatter/cli/download.py` (`list_genomes`)
- Test: `tests/unit/test_cli_blast_ani_download.py` (or new `tests/unit/test_cli_download_subspecies.py`)

**Interfaces:**
- Consumes: `GTDBClient.query_genomes`, `NCBIDatasets.fetch_assembly_summaries`, `select_subspecies_panel`, `ncbi_species_name`, `write_panel_tsv`, `GenomeAccession.subspecies`.
- Produces CLI options: `--ncbi-subspecies TEXT` (repeatable, list), `--per-subspecies INT` (default 2, min 1), `--per-subspecies-large INT` (default 3, min 1). Output file `subspecies_panel.tsv` next to `genome_metadata.tsv`.

Behaviour, per listed species `s` (after the main query, before writing files):
- `s` must start with `s__` and its name (without prefix) must be a species in the main result; otherwise exit 1 naming `s`.
- Query all genomes of `s` with `query_genomes(s, representatives_only=False)`; fetch summaries for them.
- reps-only mode: append the panel genomes not already in the list as `GenomeAccession` with the GTDB taxonomy of `s`, `is_representative=False`; set `subspecies` on every panel genome (including the representative).
- `--all-genomes` mode: set `subspecies = ncbi_species_name(organism_name) or ""` on every genome of `s` in the list; the panel TSV is still written for the record.
- Missing `datasets` or failed run: exit 1 with `NCBIDatasets.INSTALL_HINT`, before any file is written.

- [ ] **Step 1: Write failing tests** (patch `GTDBClient.query_genomes` and `NCBIDatasets.fetch_assembly_summaries`; use typer `CliRunner`):
  - `test_list_ncbi_subspecies_reps_only_appends_panel`: metadata has the representative plus panel genomes, all with `species` = GTDB species and `representative` = rep accession, and a `subspecies` column; `subspecies_panel.tsv` exists.
  - `test_list_ncbi_subspecies_all_genomes_labels_all`.
  - `test_list_ncbi_subspecies_unknown_species_exits` (Review Focus 3).
  - `test_list_ncbi_subspecies_datasets_missing_exits` (patch to raise `ToolNotFoundError`, and a second case raising `ToolExecutionError`; assert exit 1, install hint in output for the first, no `genome_metadata.tsv` written in either) (Review Focus 2).
- [ ] **Step 2: Run** the test file - expect FAIL.
- [ ] **Step 3: Implement** in `list_genomes`; keep the new code in a helper `_apply_ncbi_subspecies(...)` in `cli/download.py` so `list_genomes` grows by only the option declarations and one call.
- [ ] **Step 4: Run** - expect PASS; also run `tests/unit/test_cli_blast_ani_download.py` in full.
- [ ] **Step 5: Commit** `feat(download): --ncbi-subspecies panel selection`

### Task 5: Classifier and runner wiring

**Files:**
- Modify: `src/metadarkmatter/core/classification/classifiers/vectorized.py`
- Modify: `src/metadarkmatter/core/classification/runner.py`
- Test: `tests/unit/test_subspecies_classification.py` (new)

**Interfaces:**
- Consumes: `assign_subspecies`, `GenomeMetadata.subspecies_map`, `has_subspecies`.
- Produces: `VectorizedClassifier.__init__(..., subspecies_map: dict[str, str] | None = None)`. When non-empty:
  - Step 1 aggregation (the `best_hits` `group_by("qseqid")` near line 785) also collects `top_hit_genomes`: `genome_name` filtered to `weighted_bitscore == max`, unique, sorted.
  - `top_hit_genomes` is carried to the final column selection (near line 1177) and `assign_subspecies` is applied to the final frame, so `closest_subspecies`, `subspecies_candidates` are the last two columns.
  - Off-target rows (family validation path) get both columns null.
  - When `subspecies_map` is None or empty, no new columns and no extra aggregation.
- Runner: pass `subspecies_map=genome_metadata.subspecies_map()` when `genome_metadata` and `has_subspecies`, in both streaming and in-memory constructors; update the streaming note to say subspecies columns are still produced.

- [ ] **Step 1: Write failing tests** using a small synthetic alignment (3 reads, genomes A/S labelled, X unlabelled, ANI matrix with all three at 99.5%) via `VectorizedClassifier.classify_file`:
  - `test_classifier_adds_subspecies_columns`
  - `test_classifier_without_subspecies_map_has_no_columns`
  - `test_streaming_matches_in_memory_subspecies` (`stream_to_file` with `partition_size=2`, csv) (Review Focus 5)
- [ ] **Step 2: Run** `uv run --frozen pytest tests/unit/test_subspecies_classification.py -v` - expect FAIL.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** the new file plus `tests/unit/test_classifier_equivalence.py tests/unit/test_streaming_filter_consistency.py tests/integration/test_pipeline_determinism.py` - expect PASS.
- [ ] **Step 5: Commit** `feat(classify): closest_subspecies for Known Species reads`

### Task 6: Report subspecies table

**Files:**
- Modify: `src/metadarkmatter/core/subspecies.py`
- Modify: `src/metadarkmatter/visualization/report/templates.py`
- Modify: `src/metadarkmatter/visualization/report/sections/summary.py`
- Test: `tests/unit/test_report_generator.py` (or new `tests/unit/test_report_subspecies.py`)

**Interfaces:**
- Produces:
  - `summarize_subspecies(df: pl.DataFrame, metadata: GenomeMetadata | None) -> pl.DataFrame` with columns `gtdb_species, subspecies, reads, pct`. Uses rows with non-null `closest_subspecies`; GTDB species from the `species` column if present, else `metadata.get_species(best_match_genome)`, else `"Unknown species"`. `pct` = reads / Known Species reads of that GTDB species x 100. Sorted by gtdb_species asc, reads desc, subspecies asc. Empty frame when the column is absent.
  - `SUBSPECIES_TABLE_TEMPLATE`, `SUBSPECIES_ROW_TEMPLATE` in templates.py, reusing the `mini-table` class. Title: "Subspecies (NCBI species within GTDB species)". Header columns: GTDB species, Subspecies, Reads, "% of Known Species". Note text (exact): "Subspecies are assigned from the best-hit reference genomes of Known Species reads and are indicative for closely related taxa. Reads whose best hits span several subspecies are listed as unresolved."
  - `SummaryMixin._build_subspecies_table() -> str` ("" when no rows); rendered directly after the two-column Key Findings / Top Species row.

- [ ] **Step 1: Write failing tests**
  - `test_summarize_subspecies_counts_and_pct` (two GTDB species, one with an unresolved row; pct denominators per GTDB species).
  - `test_summarize_subspecies_uses_metadata_when_species_column_absent`.
  - `test_report_has_subspecies_table` (title and note text present) and `test_report_without_subspecies_unchanged` (title absent).
- [ ] **Step 2: Run** - expect FAIL.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** the new tests plus `tests/unit/test_report_determinism_v2.py` - expect PASS.
- [ ] **Step 5: Commit** `feat(report): subspecies table on Summary tab`

### Task 7: Pipeline options and documentation

**Files:**
- Modify: `src/metadarkmatter/models/pipeline_config.py` (fields `ncbi_subspecies: list[str] = []`, `per_subspecies: int = 2`, `per_subspecies_large: int = 3`)
- Modify: `src/metadarkmatter/cli/run.py` (matching `--ncbi-subspecies` repeatable, `--per-subspecies`, `--per-subspecies-large`; forwarded into the config dict)
- Modify: `src/metadarkmatter/core/pipeline/orchestrator.py` (append the flags to `download_argv` when `ncbi_subspecies` is non-empty)
- Modify: `CLAUDE.md` (metadata format block: optional `subspecies` column), `docs/REFERENCE.md` (new options and columns), `CHANGELOG.md` (Unreleased entry)
- Test: `tests/unit/test_pipeline_design.py`

- [ ] **Step 1: Write failing test** `test_download_argv_includes_ncbi_subspecies`: config with `ncbi_subspecies=["s__Brucella melitensis"]`, `per_subspecies=2` -> download step argv contains `["--ncbi-subspecies", "s__Brucella melitensis", "--per-subspecies", "2", "--per-subspecies-large", "3"]`; with empty list -> none of these flags.
- [ ] **Step 2: Run** `uv run --frozen pytest tests/unit/test_pipeline_design.py -v` - expect FAIL.
- [ ] **Step 3: Implement** config, CLI and orchestrator changes; update the three docs.
- [ ] **Step 4: Run** full suite `uv run --frozen pytest -q` and `uv run --frozen ruff check src tests` and `uv run --frozen mypy src` - expect all PASS.
- [ ] **Step 5: Commit** `feat(pipeline): forward NCBI subspecies options; docs`

### Task 8: End-to-end check on B_iv_01_M_S1

No code. Using the extracted, subsampled Brucella reads in `/Volumes/sekvens2/refbio_bacteria/metadarkmatter/B_iv_01_M_S1/`:

- [ ] **Step 1:** `mdm download genomes list g__Brucella --ncbi-subspecies "s__Brucella melitensis" --output refs_v2/genomes.tsv` - expect 44 representatives plus panel genomes (about 26), `refs_v2/subspecies_panel.tsv` present.
- [ ] **Step 2:** Fetch, makedb, align, ANI, classify with `--metadata refs_v2/genome_metadata.tsv`, report with `--metadata`.
- [ ] **Step 3:** Confirm `closest_subspecies` is non-null only on Known Species rows, and the report shows the subspecies table. Record the subspecies counts for the user.
