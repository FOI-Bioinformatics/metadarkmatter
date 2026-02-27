"""
Wrappers for external bioinformatics tools.

Provides Python interfaces to BLAST, Bowtie2, Kraken2,
samtools, and other pipeline dependencies.
"""

from metadarkmatter.external.base import (
    ExternalTool,
    ToolExecutionError,
    ToolNotFoundError,
    ToolResult,
    ToolTimeoutError,
)
from metadarkmatter.external.blast import BlastN, MakeBlastDb
from metadarkmatter.external.bowtie2 import Bowtie2, Bowtie2Build, concatenate_genomes
from metadarkmatter.external.fastani import FastANI, create_genome_list_file
from metadarkmatter.external.kraken import (
    ExtractKrakenReads,
    Kraken2,
    KrakenReport,
)
from metadarkmatter.external.mashtree import Mashtree
from metadarkmatter.external.mmseqs2 import MMseqs2
from metadarkmatter.external.ncbi_datasets import (
    DownloadOutcome,
    DownloadReport,
    NCBIDatasets,
)
from metadarkmatter.external.prodigal import Prodigal
from metadarkmatter.external.samtools import Samtools
from metadarkmatter.external.skani import Skani

__all__ = [
    "BlastN",
    "Bowtie2",
    "DownloadOutcome",
    "DownloadReport",
    "Bowtie2Build",
    "ExternalTool",
    "ExtractKrakenReads",
    "FastANI",
    "Kraken2",
    "KrakenReport",
    "MakeBlastDb",
    "Mashtree",
    "MMseqs2",
    "NCBIDatasets",
    "Prodigal",
    "Samtools",
    "Skani",
    "ToolExecutionError",
    "ToolNotFoundError",
    "ToolResult",
    "ToolTimeoutError",
    "concatenate_genomes",
    "create_genome_list_file",
]
