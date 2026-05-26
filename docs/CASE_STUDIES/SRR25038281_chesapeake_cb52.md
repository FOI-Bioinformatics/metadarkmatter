# Case study: Chesapeake Bay Francisellaceae strain CB52 (SRR25038281)

End-to-end run of the metadarkmatter pipeline against a published novel
microbial isolate. The result agrees with the paper's claim that CB52
is a novel Francisellaceae bacterium and demonstrates the package
working on real (non-synthetic) sequencing data.

## Sample

| Attribute | Value |
|---|---|
| SRA run | SRR25038281 |
| BioProject | PRJNA929474 (Chesapeake Bay environmental bacteria) |
| BioSample | SAMN34231785 |
| Sample name | W20_SRM_FM12 |
| Organism (declared) | Francisellaceae bacterium CB52 (taxon 3042159) |
| Source | Chesapeake Bay water, USNA, 2020-02-20 |
| Platform | Illumina NovaSeq 6000, paired-end |
| Read count | 3,829,373 spots (1,125 Mbp) |
| Submitter | United States Naval Academy |

## Pipeline invocation

The full v0.2.0 pipeline. Kraken2 step skipped because the sample is
already a Francisella isolate (no host or off-target contamination
worth filtering at the family level).

```bash
# Reference: 48 GTDB f__Francisellaceae representatives
mdm download genomes list "f__Francisellaceae" --output genomes.tsv
mdm download genomes fetch --accessions genomes.tsv --output-dir genomes/
mdm blast makedb --genomes genomes/ --output blastdb/pangenome
mdm ani compute --genomes genomes/ --output ani_matrix.csv --threads 8

# CB52 reads downloaded from ENA mirror
curl -O https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR250/081/SRR25038281/SRR25038281_1.fastq.gz
curl -O https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR250/081/SRR25038281/SRR25038281_2.fastq.gz

mdm blast align \
    --query SRR25038281_1.fastq.gz \
    --database blastdb/pangenome \
    --output cb52.blast.tsv.gz \
    --threads 8 \
    --task dc-megablast \
    --word-size 11

# Streaming mode because the alignment is ~10M hits and the default
# vectorized classifier needs ~18 GB RAM for this input.
mdm score classify \
    --alignment cb52.blast.tsv.gz \
    --ani ani_matrix.csv \
    --metadata genome_metadata.tsv \
    --output classifications.csv \
    --streaming --chunk-size 500000

mdm report generate \
    --classifications classifications.csv \
    --metadata genome_metadata.tsv \
    --ani ani_matrix.csv \
    --output report.html \
    --sample-name "Chesapeake Bay CB52 (SRR25038281)"
```

## Result

The R1 mate alone produced 3,334,828 classified reads in 4 seconds of
streaming-mode classifier wall time (10M alignments).

### Category distribution

| Taxonomic call | Reads | Fraction | Mean pident |
|---|---|---|---|
| Novel Species | 1,321,650 | 39.6% | 88.9% |
| Known Species | 1,119,438 | 33.6% | 96.6% |
| Ambiguous | 400,735 | 12.0% | 89.7% |
| Novel Genus | 251,008 | 7.5% | 78.0% |
| Conserved Region | 219,227 | 6.6% | 87.1% |
| Unclassified | 20,427 | 0.6% | 71.1% |
| Species Boundary | 2,343 | 0.07% | 88.9% |

The mean percent-identity per category tracks the expected ordering
(Known Species highest, Unclassified lowest), confirming the
classifier behaves as designed on real data.

### Closest known references

| Reference (GTDB rep) | Reads matched | Mean pident |
|---|---|---|
| Francisella sp038141505 (`GCA_038141505.1`) | 1,673,437 | 94.1% |
| Francisella sp047609165 (`GCA_047609165.1`) | 725,859 | 89.1% |
| Francisella sp965240355 (`GCA_965240355.1`) | 344,755 | 88.6% |

CB52's nearest GTDB reference, `Francisella sp038141505`, is itself
an uncharacterised Francisella species (a `sp...` placeholder name in
GTDB, no formal species description). CB52 matches it at ~94% mean
identity, just below the 95-96% species boundary. This is exactly
the regime metadarkmatter is designed to flag: a novel species that
sits close to but distinct from the nearest known reference.

### Verification against the paper

PRJNA929474 describes "novel bacterial strains and species from the
Chesapeake Bay." The sample title declares CB52 as
`Francisellaceae bacterium CB52` rather than a named species,
consistent with the submitting lab considering it a novel isolate.
The metadarkmatter classification (39.6% Novel Species reads,
mean ~94% identity to the closest Francisella GTDB representative)
agrees with that designation.

## Notes for reproducibility

- **Disk.** The uncompressed BLAST output reached ~10 GB before
  gzipping to 1.3 GB. Run with `--task dc-megablast --word-size 11`
  on a host with at least 15 GB free on the working partition.
- **Memory.** The vectorized classifier estimated 18 GB RAM for
  this alignment. `--streaming --chunk-size 500000` keeps memory
  bounded and finishes in seconds.
- **Backend.** `mdm ani compute` autoselected `fastANI` because it
  was on PATH. `skani` would also work; both produce the same
  category-level result for this sample.
- **R2 not aligned.** Only R1 was aligned because the BLAST wrapper
  expects single-end input. Adding `--query` support for paired-end
  reads is on the roadmap for a future release.

## Reproduction script

See `scripts/run_pipeline.sh` for the canonical orchestration of
these steps. Adapt the family, kraken2 step, and target-taxid for
other samples.

## Generated outputs

- `classifications.csv` (727 MB) — full per-read classification
- `report.html` (102 MB) — self-contained interactive report
