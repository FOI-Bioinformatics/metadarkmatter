"""
Configuration model for the end-to-end ``mdm run`` pipeline.

``PipelineConfig`` captures everything needed to drive the full workflow
(download -> kraken2 extract -> align -> ANI -> classify -> report) from a
single command. It can be built from a YAML file (``from_yaml``) and/or from
CLI flags, mirroring the precedence used by ``score classify``.

The pipeline reuses the existing ``mdm`` subcommands as steps, so this config
only needs the cross-cutting inputs and toggles; per-step tuning is left to the
individual commands' own defaults (or a future nested config).
"""

from __future__ import annotations

from enum import StrEnum
from pathlib import Path
from typing import Any, Self

import yaml
from pydantic import BaseModel, Field, model_validator


class Aligner(StrEnum):
    """Sequence aligner used for the alignment step."""

    BLAST = "blast"
    MMSEQS2 = "mmseqs2"


class PipelineStep(StrEnum):
    """Ordered pipeline steps. The order here is the execution order."""

    DOWNLOAD = "download"
    KRAKEN_CLASSIFY = "kraken-classify"
    KRAKEN_EXTRACT = "kraken-extract"
    MAKEDB = "makedb"
    ALIGN = "align"
    ANI = "ani"
    CLASSIFY = "classify"
    REPORT = "report"

    @classmethod
    def ordered(cls) -> list[PipelineStep]:
        return [
            cls.DOWNLOAD,
            cls.KRAKEN_CLASSIFY,
            cls.KRAKEN_EXTRACT,
            cls.MAKEDB,
            cls.ALIGN,
            cls.ANI,
            cls.CLASSIFY,
            cls.REPORT,
        ]


class PipelineConfig(BaseModel):
    """Inputs and toggles for a full ``mdm run`` pipeline.

    Steps are inferred from which inputs are supplied (e.g. a pre-populated
    ``genomes_dir`` skips download; a pre-extracted ``extracted_reads_1`` skips
    the kraken2 steps). Use ``report`` / ``compute_aai`` to toggle optional
    stages explicitly.
    """

    # Identification
    family: str | None = Field(
        default=None,
        description="GTDB family (e.g. 'f__Francisellaceae') to download references for.",
    )
    family_taxid: int | None = Field(
        default=None,
        description="NCBI taxid for the family; required for the kraken2 extraction step.",
    )

    # Inputs
    reads_1: Path | None = Field(default=None, description="Path to R1 reads (FASTQ[.gz]).")
    reads_2: Path | None = Field(default=None, description="Path to R2 reads for paired-end input.")
    kraken_db: Path | None = Field(default=None, description="Kraken2 database directory.")
    genomes_dir: Path | None = Field(
        default=None,
        description="Pre-downloaded reference genomes directory; when set, download is skipped.",
    )
    kraken_output: Path | None = Field(
        default=None,
        description="Pre-computed Kraken2 .kraken file; when set, kraken classify is skipped.",
    )
    extracted_reads_1: Path | None = Field(
        default=None,
        description="Pre-extracted family R1 reads; when set, both kraken steps are skipped.",
    )

    # Output + execution
    output_dir: Path = Field(description="Root output directory for the run.")
    aligner: Aligner = Field(default=Aligner.BLAST, description="Aligner for the alignment step.")
    threads: int = Field(default=8, ge=1, description="Threads passed to each tool step.")

    # Optional stage toggles
    compute_aai: bool = Field(default=False, description="Also compute an AAI matrix.")
    report: bool = Field(default=True, description="Generate the HTML report at the end.")

    model_config = {"frozen": True}

    @model_validator(mode="after")
    def _validate_minimal_inputs(self) -> Self:
        """Ensure the supplied inputs can produce at least an alignment.

        We require *either* pre-extracted reads, *or* the inputs needed to
        produce them (reads + a way to get/extract family reads). Reference
        genomes must be obtainable (a genomes_dir or a family to download).
        """
        if self.genomes_dir is None and self.family is None:
            msg = "Provide --genomes (pre-downloaded references) or --family (to download them)."
            raise ValueError(msg)

        has_reads_for_extraction = self.reads_1 is not None
        if self.extracted_reads_1 is None and not has_reads_for_extraction:
            msg = "Provide --reads-1 (raw reads) or --extracted-reads-1 (pre-extracted family reads)."
            raise ValueError(msg)

        # If we must run kraken2 extraction, we need a taxid and a source of the
        # .kraken file (either pre-computed or a kraken_db to compute it).
        needs_extraction = self.extracted_reads_1 is None
        if needs_extraction:
            if self.family_taxid is None:
                msg = "--family-taxid is required to extract family reads with Kraken2."
                raise ValueError(msg)
            if self.kraken_output is None and self.kraken_db is None:
                msg = "Provide --kraken-db (to classify) or --kraken-output (a precomputed .kraken)."
                raise ValueError(msg)
        return self

    @classmethod
    def from_yaml(cls, path: Path, **overrides: Any) -> PipelineConfig:
        """Load a pipeline config from YAML, applying non-None CLI overrides.

        YAML keys map directly to field names. ``overrides`` (typically CLI
        flags) take precedence when not None, matching the precedence used by
        ``score classify``.
        """
        raw = yaml.safe_load(path.read_text()) or {}
        if not isinstance(raw, dict):
            msg = f"Pipeline YAML must be a mapping, got {type(raw).__name__}"
            raise ValueError(msg)
        merged = {**raw, **{k: v for k, v in overrides.items() if v is not None}}
        return cls(**merged)
