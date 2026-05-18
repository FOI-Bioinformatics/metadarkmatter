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
| `blastn` / `blastx` / `makeblastdb` | 2.16.0+ | competitive read recruitment | `conda install -c bioconda blast` |
| `mmseqs` | 17-b804f or newer | fast alignment alternative to BLAST | `conda install -c bioconda mmseqs2` |
| `kraken2` | 2.1.3 or newer | per-read taxonomic pre-filter | `conda install -c bioconda kraken2` |
| `fastANI` | 1.34 or newer | ANI matrix construction | `conda install -c bioconda fastani` |
| `skani` | 0.3.1 or newer | faster ANI matrix alternative | `conda install -c bioconda skani` |
| `diamond` | 2.1.9 or newer | protein-level alignment for AAI | `conda install -c bioconda diamond` |
| `mashtree` | 1.4 or newer | optional alternative phylogeny | `conda install -c bioconda mashtree` |
| `mash` | 2.3 or newer | required by mashtree | `conda install -c bioconda mash` |

Only the tools relevant to your workflow need to be installed. For
example, a nucleotide-only pipeline with BLAST and FastANI does not
require Diamond or MMseqs2.

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
with Python 3.12, the locked Python environment, and all external
bioinformatics tools pre-installed from bioconda. Build and run:

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
