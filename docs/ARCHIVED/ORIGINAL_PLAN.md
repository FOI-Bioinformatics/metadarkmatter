# ARCHIVED: Original Implementation Plan

> **ARCHIVED DOCUMENT**: This describes the originally planned Snakemake workflow architecture that was NOT implemented. The current metadarkmatter package uses a CLI-based approach instead. See the main documentation for current usage.

---

The content below is preserved for historical reference only.

---

# Implementation Strategy: Hybrid Snakemake + Python

## Original Recommendation

The original plan recommended combining:
- **Snakemake** for workflow orchestration (Phases I & II)
- **Optimized Python package** for Phase III ANI-weighted placement algorithm

## What Was Actually Built

The project implemented a **CLI-based Python package** that provides individual commands for each workflow step:

- `metadarkmatter download` - Genome acquisition from GTDB/NCBI
- `metadarkmatter extract` - Kraken2 wrapper for read extraction
- `metadarkmatter blast` - BLASTN wrapper for competitive alignment
- `metadarkmatter map` - Bowtie2 wrapper for read mapping
- `metadarkmatter score` - Core ANI-weighted classification algorithm
- `metadarkmatter visualize` - Recruitment plot generation
- `metadarkmatter report` - HTML report generation

## Why Snakemake Was Not Implemented

The CLI-based approach provides:
1. Simpler installation (single pip install)
2. More flexibility for users to integrate with existing pipelines
3. Easier debugging of individual steps
4. No Snakemake dependency required

Users who need workflow automation can easily wrap the CLI commands with their preferred workflow manager (Snakemake, Nextflow, shell scripts).

---

## Original Snakemake Directory Structure (Not Built)

```
workflow/
├── Snakefile
├── rules/
│   ├── extraction.smk
│   ├── indexing.smk
│   ├── mapping.smk
│   └── ani_placement.smk
├── envs/
│   ├── kraken2.yaml
│   ├── blast.yaml
│   └── mapping.yaml
└── scripts/
    └── concatenate_genomes.py
```

This structure was never created.
