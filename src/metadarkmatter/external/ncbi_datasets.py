"""
NCBI Datasets CLI wrapper for genome downloads.

Provides a Python interface to the NCBI datasets command-line tool
for downloading genome assemblies by accession.
"""

from __future__ import annotations

import logging
import re
import shutil
import time
import zipfile
from collections.abc import Callable
from pathlib import Path
from typing import ClassVar

from metadarkmatter.external.base import ExternalTool, ToolResult

# Rate limiting defaults for NCBI API
DEFAULT_BATCH_DELAY = 0.5  # seconds between batch requests

logger = logging.getLogger(__name__)

# Pattern for valid NCBI accession (GCF/GCA followed by underscore and digits)
_ACCESSION_PATTERN = re.compile(r"^GC[AF]_\d+\.\d+$")


class NCBIDatasets(ExternalTool):
    """Wrapper for NCBI datasets CLI tool.

    The NCBI datasets CLI is the official tool for downloading
    data from NCBI, including genome assemblies.

    Install via conda:
        conda install -c conda-forge ncbi-datasets-cli
    """

    TOOL_NAME: ClassVar[str] = "datasets"
    TOOL_ALIASES: ClassVar[tuple[str, ...]] = ("ncbi-datasets-cli",)
    INSTALL_HINT: ClassVar[str] = "conda install -c conda-forge ncbi-datasets-cli"

    def build_command(self, **kwargs: object) -> list[str]:
        """Build download command for genome accessions.

        Keyword Args:
            accessions: List of NCBI accessions to download
            output_file: Path for output zip file
            include_sequence: Include genome sequences (default True)
            include_protein: Include protein FASTA files (default False)

        Returns:
            Command as list of strings
        """
        accessions = kwargs.get("accessions", [])
        output_file = kwargs.get("output_file")
        include_sequence = kwargs.get("include_sequence", True)
        include_protein = kwargs.get("include_protein", False)

        if not accessions:
            msg = "accessions parameter is required"
            raise ValueError(msg)
        if not output_file:
            msg = "output_file parameter is required"
            raise ValueError(msg)

        cmd = [
            str(self.get_executable()),
            "download",
            "genome",
            "accession",
        ]

        # Add accessions
        if isinstance(accessions, list):
            cmd.extend(accessions)
        else:
            cmd.append(str(accessions))

        # Output file
        cmd.extend(["--filename", str(output_file)])

        # Build include list
        include_types = []
        if include_sequence:
            include_types.append("genome")
        if include_protein:
            include_types.append("protein")

        if include_types:
            cmd.extend(["--include", ",".join(include_types)])

        return cmd

    def build_download_command(
        self,
        accessions: list[str],
        output_file: Path,
        *,
        include_sequence: bool = True,
        include_protein: bool = False,
    ) -> list[str]:
        """Build command to download genomes by accession.

        Args:
            accessions: List of NCBI accessions
            output_file: Path for output zip file
            include_sequence: Include genome FASTA sequences
            include_protein: Include protein FASTA files

        Returns:
            Command as list of strings
        """
        return self.build_command(
            accessions=accessions,
            output_file=output_file,
            include_sequence=include_sequence,
            include_protein=include_protein,
        )

    def download_batch(
        self,
        accessions: list[str],
        output_dir: Path,
        *,
        dry_run: bool = False,
        timeout: float | None = 1800.0,  # 30 minutes default
    ) -> ToolResult:
        """Download a batch of genomes.

        Args:
            accessions: List of NCBI accessions
            output_dir: Output directory for downloaded files
            dry_run: If True, show command without executing
            timeout: Timeout in seconds for download

        Returns:
            ToolResult with command and output
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        output_zip = output_dir / "ncbi_dataset.zip"

        return self.run(
            accessions=accessions,
            output_file=output_zip,
            include_sequence=True,
            dry_run=dry_run,
            timeout=timeout,
        )

    def download_genomes(
        self,
        accessions: list[str],
        output_dir: Path,
        *,
        batch_size: int = 100,
        dry_run: bool = False,
        decompress: bool = True,
        skip_if_exists: bool = True,
        include_protein: bool = False,
        timeout_per_batch: float = 1800.0,
        batch_delay: float = DEFAULT_BATCH_DELAY,
        progress_callback: Callable[[int, int], None] | None = None,
    ) -> ToolResult:
        """Download genomes in batches and extract FASTA files.

        Downloads genomes from NCBI in batches, extracts the FASTA files,
        and optionally decompresses them. Includes rate limiting to avoid
        overwhelming the NCBI API.

        Args:
            accessions: List of NCBI accessions to download.
            output_dir: Output directory for genome FASTA files.
            batch_size: Number of accessions per batch request.
            dry_run: If True, show commands without executing.
            decompress: If True, decompress .gz files.
            skip_if_exists: Skip accessions with existing FASTA files.
            include_protein: If True, also download protein FASTA files (.faa).
            timeout_per_batch: Timeout per batch in seconds.
            batch_delay: Delay between batch requests in seconds (rate limiting).
            progress_callback: Optional callback(completed, total).

        Returns:
            ToolResult with combined output from all batches.
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        # Filter accessions if skip_if_exists
        if skip_if_exists:
            accessions_to_download = []
            for acc in accessions:
                # Check for both .fna and .fna.gz
                fasta_path = output_dir / f"{acc}.fna"
                fasta_gz_path = output_dir / f"{acc}.fna.gz"
                if not fasta_path.exists() and not fasta_gz_path.exists():
                    accessions_to_download.append(acc)
            accessions = accessions_to_download

        if not accessions:
            return ToolResult(
                command=("datasets", "download", "genome", "accession"),
                return_code=0,
                stdout="All genomes already exist, nothing to download.",
                stderr="",
                elapsed_seconds=0.0,
            )

        total = len(accessions)
        completed = 0
        all_stdout = []
        all_stderr = []
        total_elapsed = 0.0

        # Process in batches
        for i in range(0, total, batch_size):
            batch = accessions[i : i + batch_size]
            batch_zip = output_dir / f"ncbi_dataset_batch_{i}.zip"

            if dry_run:
                cmd = self.build_download_command(
                    batch, batch_zip, include_protein=include_protein
                )
                all_stdout.append(f"[dry-run] Would execute: {' '.join(cmd)}")
                completed += len(batch)
                if progress_callback:
                    progress_callback(completed, total)
                continue

            # Download batch
            result = self.run(
                accessions=batch,
                output_file=batch_zip,
                include_sequence=True,
                include_protein=include_protein,
                timeout=timeout_per_batch,
            )

            all_stdout.append(result.stdout)
            all_stderr.append(result.stderr)
            total_elapsed += result.elapsed_seconds

            if result.success and batch_zip.exists():
                # Extract FASTA files (nucleotide and protein if requested)
                self._extract_fastas(
                    batch_zip,
                    output_dir,
                    decompress=decompress,
                    include_protein=include_protein,
                )
                # Clean up zip file
                batch_zip.unlink()

            completed += len(batch)
            if progress_callback:
                progress_callback(completed, total)

            # Rate limiting: delay between batches (skip for last batch)
            if batch_delay > 0 and i + batch_size < total:
                time.sleep(batch_delay)

        return ToolResult(
            command=("datasets", "download", "genome", "accession", f"[{total} accessions]"),
            return_code=0,
            stdout="\n".join(all_stdout),
            stderr="\n".join(all_stderr),
            elapsed_seconds=total_elapsed,
        )

    def _validate_accession(self, accession: str) -> bool:
        """Validate NCBI accession format for security.

        Args:
            accession: Potential accession string to validate

        Returns:
            True if accession matches expected NCBI format
        """
        return bool(_ACCESSION_PATTERN.match(accession))

    def _safe_output_path(
        self,
        filename: str,
        output_dir: Path,
    ) -> Path | None:
        """Validate that output path stays within output directory.

        Prevents Zip Slip path traversal attacks by ensuring the resolved
        output path is within the intended output directory.

        Args:
            filename: Proposed output filename
            output_dir: Intended output directory

        Returns:
            Safe resolved path, or None if path would escape output directory
        """
        # Construct the target path
        target = (output_dir / filename).resolve()
        output_resolved = output_dir.resolve()

        # Ensure target is within output directory
        try:
            target.relative_to(output_resolved)
            return target
        except ValueError:
            # Path would escape output directory
            logger.warning(
                "Skipping potentially malicious path: %s (would escape %s)",
                filename,
                output_dir,
            )
            return None

    def _extract_fastas(
        self,
        zip_path: Path,
        output_dir: Path,
        *,
        decompress: bool = True,
        include_protein: bool = False,
    ) -> list[Path]:
        """Extract FASTA files from NCBI dataset zip.

        Args:
            zip_path: Path to downloaded zip file
            output_dir: Output directory for FASTA files
            decompress: If True, decompress .gz files
            include_protein: If True, also extract protein files (.faa)

        Returns:
            List of extracted FASTA file paths
        """
        extracted = []

        # Build list of extensions to extract
        nucleotide_exts = (".fna", ".fna.gz")
        protein_exts = (".faa", ".faa.gz")
        valid_exts = nucleotide_exts + (protein_exts if include_protein else ())

        with zipfile.ZipFile(zip_path, "r") as zf:
            for name in zf.namelist():
                # NCBI datasets puts FASTAs in ncbi_dataset/data/ACCESSION/
                if not name.endswith(valid_exts):
                    continue

                # Determine file type (nucleotide or protein)
                is_protein = name.endswith(protein_exts)

                # Extract accession from path
                parts = Path(name).parts
                # Find the accession directory (usually after "data")
                accession = None
                for j, part in enumerate(parts):
                    if part == "data" and j + 1 < len(parts):
                        accession = parts[j + 1]
                        break

                if not accession:
                    continue

                # Validate accession format (security check)
                if not self._validate_accession(accession):
                    logger.warning(
                        "Skipping file with invalid accession format: %s",
                        accession,
                    )
                    continue

                # Determine output filename based on file type
                if is_protein:
                    if name.endswith(".faa.gz"):
                        out_name = f"{accession}.faa.gz"
                    else:
                        out_name = f"{accession}.faa"
                else:
                    if name.endswith(".fna.gz"):
                        out_name = f"{accession}.fna.gz"
                    else:
                        out_name = f"{accession}.fna"

                # Validate output path (Zip Slip protection)
                out_path = self._safe_output_path(out_name, output_dir)
                if out_path is None:
                    continue

                # Extract file
                with zf.open(name) as src, out_path.open("wb") as dst:
                    shutil.copyfileobj(src, dst)

                # Decompress if requested
                if decompress and out_path.suffix == ".gz":
                    import gzip

                    decompressed_path = out_path.with_suffix("")
                    with (
                        gzip.open(out_path, "rb") as gz_in,
                        decompressed_path.open("wb") as out,
                    ):
                        shutil.copyfileobj(gz_in, out)
                    out_path.unlink()
                    out_path = decompressed_path

                extracted.append(out_path)

        return extracted


class NCBIDatazip(ExternalTool):
    """Wrapper for NCBI dataformat tool (for processing dataset files)."""

    TOOL_NAME: ClassVar[str] = "dataformat"
    TOOL_ALIASES: ClassVar[tuple[str, ...]] = ()
    INSTALL_HINT: ClassVar[str] = "conda install -c conda-forge ncbi-datasets-cli"

    def build_command(self, **kwargs: object) -> list[str]:
        """Build dataformat command."""
        subcommand = kwargs.get("subcommand", "tsv")
        package = kwargs.get("package", "genome")

        cmd = [
            str(self.get_executable()),
            subcommand,
            package,
        ]

        input_file = kwargs.get("input_file")
        if input_file:
            cmd.extend(["--package", str(input_file)])

        return cmd
