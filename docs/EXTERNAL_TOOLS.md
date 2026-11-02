# External Tools and Reproducibility

metadarkmatter shells out to a small set of bioinformatics tools. This
document records the **tested** versions and the conventions used to
keep results reproducible across lab machines.

## Python environment

The repository ships a `uv.lock` file pinning every transitive Python
dependency. To recreate the exact environment a result was produced in:

```bash
uv sync --frozen
```

(or, with pip: `pip install -r <(uv export --frozen --no-hashes)`).

The Python interpreter is constrained to `>=3.11`. The CI matrix
exercises 3.11 and 3.12; other versions are unsupported.

The lockfile is regenerated with `uv lock --upgrade` whenever a
runtime dependency is bumped intentionally. Do not edit `uv.lock`
by hand.

## External tools

Tools below are shelled out via `subprocess.run` from
`metadarkmatter.external.*` wrappers. Each row gives the tool, the
version we test against, and the install hint that the wrappers print
when a tool is missing. Run `metadarkmatter doctor` on a target
machine to verify the installed versions.

| Tool | Tested version | Purpose | Install hint |
|------|----------------|---------|--------------|
| `blastn` / `makeblastdb` | 2.16.0+ | competitive read recruitment (`blast makedb` / `blast align`) | `conda install -c bioconda blast` |
| `mmseqs` | 17-b804f or newer | fast alignment alternative to BLAST | `conda install -c bioconda mmseqs2` |
| `kraken2` | 2.1.3 or newer | per-read taxonomic pre-filter | `conda install -c bioconda kraken2` |
| `fastANI` | 1.34 or newer | ANI matrix construction | `conda install -c bioconda fastani` |
| `skani` | 0.3.1 or newer | faster ANI matrix alternative | `conda install -c bioconda skani` |
| `diamond` | 2.1.9 or newer | protein-level read alignment (`blastx makedb` / `blastx align`) and AAI | `conda install -c bioconda diamond` |
| `mashtree` | 1.4 or newer | optional alternative phylogeny | `conda install -c bioconda mashtree` |
| `mash` | 2.3 or newer | required by mashtree | `conda install -c bioconda mash` |

Further tools used by specific commands (no tested version recorded;
not reported by `doctor`):

| Tool | Used by | Install hint |
|------|---------|--------------|
| `extract_kraken_reads.py` (KrakenTools) | `kraken2 extract` | `conda install -c bioconda krakentools` |
| `bowtie2` / `bowtie2-build` | `map reads` | `conda install -c bioconda bowtie2` |
| `samtools` | `map reads`, BAM-based recruitment | `conda install -c bioconda samtools` |
| `prodigal` | `proteins predict` | `conda install -c bioconda prodigal` |
| `datasets` / `dataformat` (NCBI) | `download genomes fetch` | `conda install -c conda-forge ncbi-datasets-cli` |

Only the tools relevant to your workflow need to be installed. For
example, a nucleotide-only pipeline with BLAST and FastANI does not
require Diamond or MMseqs2.

The wrappers do not enforce minimum versions; a tool is only required
to be on `PATH`. `doctor` reports the version strings of `blastn`,
`blastx`, `makeblastdb`, `mmseqs`, `kraken2`, `fastANI`, `skani`,
`diamond`, `mashtree` and `mash` (`blastx` is listed by `doctor`, but
`blastx align` runs Diamond, not NCBI `blastx`). External tools run
without a default timeout; when a caller sets one explicitly and it
is exceeded, this is reported as a tool error.

## Reproducibility checklist

Before publishing a result we recommend the following steps:

1. Run `uv sync --frozen` to recreate the pinned Python environment.
2. Run `metadarkmatter doctor > doctor.txt` and archive the output
   alongside the result. This records Python version, package
   versions, and external-tool versions and paths.
3. Set `METADARKMATTER_SEED=<integer>` if you want to use a non-default
   seed; every randomness-using site in metadarkmatter (plot
   subsampling, adaptive GMM fits) reads this single variable, defaulting
   to 42 when unset.
4. For pipeline-level reproducibility, record the SHA of the
   metadarkmatter commit and the exact CLI invocation.

## Container image

A `Containerfile` at the repository root builds an internal-lab image
with Python 3.12, the locked Python environment, and the external
tools from the first table (plus `samtools`) pre-installed from
bioconda. The Containerfile specifies lower bounds (e.g. `blast>=2.16`,
`mmseqs2>=15`, `skani>=0.3`, `diamond>=2.1`) rather than exact versions,
so record `doctor` output from the built image. KrakenTools, Bowtie2,
Prodigal and the NCBI datasets CLI are not included. Build and run:

```bash
podman build -t metadarkmatter:dev .
# or: docker build -t metadarkmatter:dev .

podman run --rm -v "$PWD/data:/data" metadarkmatter:dev doctor
podman run --rm -v "$PWD/data:/data" metadarkmatter:dev \
    score classify --alignment /data/sample.blast.tsv.gz \
                   --ani /data/ani_matrix.csv \
                   --output /data/classifications.csv
```

The build step runs `metadarkmatter doctor` as a smoke test so a broken
image fails the build instead of failing at first use. The image is
not pushed to a public registry; pass a build tag suited to the local
infrastructure (e.g. an internal registry path).

## Out of scope

A public Docker image and PyPI release are not currently planned;
metadarkmatter is operated as an internal lab tool.
