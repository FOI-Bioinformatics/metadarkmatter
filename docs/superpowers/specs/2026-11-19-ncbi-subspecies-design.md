# NCBI subspecies labels within a GTDB species

Date: 2026-11-19
Status: draft for approval (builder check done)

## Intent

GTDB merges several NCBI species into one GTDB species. The clearest case is
s__Brucella melitensis, which contains NCBI B. melitensis, B. abortus, B. suis,
B. canis, B. ovis, B. ceti, B. pinnipedialis, B. neotomae, B. microti,
B. inopinata, B. vulpis and B. amazoniensis (1,337 genomes in GTDB as of
2026-11-19). For Brucella samples the user wants the report to state the
closest NCBI species for reads placed in such a GTDB species. The NCBI name is
treated as a "subspecies" level below the GTDB species; the GTDB species call
and all ANI-based classification logic are unchanged.

## Goals

- Select, for a given GTDB species, a small set of reference genomes per NCBI
  species (the "subspecies panel") with a documented, deterministic rule.
- Record the NCBI species of each panel genome in genome_metadata.tsv as an
  optional `subspecies` column.
- Report a per-read `closest_subspecies` in the classification output, with an
  explicit "unresolved" value when the best hits do not discriminate.
- Show subspecies read counts in the HTML report, grouped under the GTDB
  species.
- Work for any GTDB species, not only Brucella.

## Non-goals

- No change to taxonomic_call, novelty, ambiguity or species-boundary logic.
- No separate re-alignment or competitive mapping step (approach B, rejected).
- No subspecies-level ANI thresholds or novelty calls.

## Anti-goals

- Panel genomes must not become separate GTDB species (approach C). All panel
  genomes keep the GTDB species name and the GTDB species representative in
  the `representative` column, so ANI lookups and ambiguity behave as before.
- The report must not present a subspecies label as more certain than the
  data allow: ties are reported as "unresolved", and the report states that
  subspecies assignment from short reads is indicative.

## Constraints

- Members of the classical Brucella group share >99% ANI. Single 150 bp reads
  will often hit several subspecies with identical scores.
- The GTDB API genomes-detail endpoint does not return NCBI taxonomy. NCBI
  organism names come from NCBI Datasets (`datasets summary genome accession`),
  which is already an external dependency candidate; if unavailable, the
  command errors with an install hint.
- genome_metadata.tsv files without the `subspecies` column must keep working
  unchanged.

## Approaches considered

- A. Metadata column plus selection command, subspecies taken from the
  best-hit genome with tie handling. Chosen: reuses the existing alignment,
  no extra compute, general.
- B. Second-stage competitive mapping of reads in the GTDB species against
  the panel. Rejected: a full extra step for modest gain at >99% ANI.
- C. Treat NCBI species as separate species in metadata. Rejected: breaks
  representative mapping and species-boundary logic.

## Design

### 1. Panel selection (download)

New options on `metadarkmatter download genomes list`:

- `--ncbi-subspecies TEXT` (repeatable): GTDB species (e.g.
  "s__Brucella melitensis") to expand into a subspecies panel.
- `--per-subspecies INT` (default 2): genomes kept per NCBI species.
- `--per-subspecies-large INT` (default 3) and threshold 50 genomes: NCBI
  species with more than 50 genomes in the GTDB species keep this many,
  to sample more of their diversity.

Procedure per listed GTDB species:

1. Query all genomes of the GTDB species (`--all-genomes` behaviour).
2. Fetch NCBI organism name, assembly level, RefSeq category and contig N50
   for those accessions via NCBI Datasets.
3. NCBI species = first two words of the organism name. Drop "Genus sp.".
4. Rank genomes within each NCBI species by: RefSeq category "reference
   genome" first; assembly level (Complete Genome, Chromosome, Scaffold,
   Contig); GCF before GCA; contig N50 descending; accession ascending.
5. Keep the top N. Always keep the GTDB species representative.
6. Append panel genomes to the accession list and genome_metadata.tsv with:
   species = GTDB species name, representative = GTDB representative,
   subspecies = NCBI species. Non-panel genomes get an empty subspecies.
7. Write `subspecies_panel.tsv` next to genome_metadata.tsv recording, per
   selected genome: accession, gtdb_species, ncbi_organism_name, subspecies,
   assembly level, RefSeq category, contig N50, rank, and the count of
   genomes per NCBI species in the GTDB species. This documents the selection.

`download genomes fetch` needs no change.

### 2. Metadata

`GenomeMetadata` reads an optional `subspecies` column (empty or null means
none), exposes `get_subspecies(accession)` and `has_subspecies`, and
`join_classifications` adds `subspecies` when present.

### 3. Classification

When metadata has subspecies, `score classify` adds two columns:

- `closest_subspecies` (Known Species reads only, else null): the subspecies shared by all hits whose
  weighted_bitscore equals the read's maximum; "unresolved" if those top hits
  carry more than one distinct subspecies; null if the best-hit genome has no
  subspecies.
- `subspecies_candidates`: semicolon-separated sorted list of distinct
  subspecies among those top hits (empty when null).

Only exact ties at the maximum weighted_bitscore count; no tolerance
parameter in this version. Columns are absent when metadata has no
subspecies, so existing outputs do not change.

### 4. Report

When `closest_subspecies` is present, the Summary tab gains a "Subspecies
(NCBI species within GTDB species)" table: one block per GTDB species with
reads, rows = subspecies including "unresolved", columns = reads and percent
of that GTDB species' Known Species reads, sorted by reads. Only Known
Species reads are counted. A one-line note states that subspecies
are assigned from best-hit genomes and are indicative for closely related
taxa.

### 5. Pipeline

`run pipeline` already passes `--metadata`; no change beyond forwarding
`--ncbi-subspecies` / `--per-subspecies` when the pipeline performs the
download step.

## Decisions after builder check

- With `--all-genomes`, every genome of a listed GTDB species gets its NCBI
  species as `subspecies` (the NCBI lookup covers all of them anyway); the
  panel rule only limits which genomes are added in the default reps-only
  mode. (Claude's call.)
- Among tied top hits, genomes without a subspecies label are ignored; null
  only when no tied hit is labelled. Reads whose best hits fall outside the
  expanded GTDB species are not labelled. (Claude's call.)
- `score classify --streaming` also produces `closest_subspecies` and
  `subspecies_candidates`, since Brucella-rich samples are the ones that need
  streaming; the report rebuilds the GTDB species grouping from metadata
  when the species column is absent. (Claude's call.)
- `run pipeline` gets `--ncbi-subspecies` / `--per-subspecies` CLI flags and
  matching config fields; its download step stays reps-only. (Claude's call.)
- Only reads with taxonomic_call == "Known Species" receive a subspecies
  label; all other calls get null in both subspecies columns. The report
  table counts Known Species reads only, and percentages use the Known
  Species reads of that GTDB species as denominator. (User decision.)

## Left to the builder

- Reuse `src/metadarkmatter/external/ncbi_datasets.py` for the NCBI lookup.
- Module placement of the NCBI Datasets client and the panel selection code.
- Caching of NCBI Datasets responses (follow the GTDB client cache pattern).
- Exact HTML/CSS of the report table (follow existing templates).
- Tests: selection ranking, metadata backward compatibility, tie handling,
  report rendering with and without subspecies.
