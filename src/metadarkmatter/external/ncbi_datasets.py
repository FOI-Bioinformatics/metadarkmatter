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
from dataclasses import dataclass, field
from pathlib import Path
from typing import ClassVar

from metadarkmatter.external.base import ExternalTool, ToolResult

# Rate limiting defaults for NCBI API
DEFAULT_BATCH_DELAY = 0.5  # seconds between batch requests

logger = logging.getLogger(__name__)

# Pattern for valid NCBI accession (GCF/GCA followed by underscore and digits)
_ACCESSION_PATTERN = re.compile(r"^GC[AF]_\d+\.\d+$")

# Pattern to split accession into base and version: GCF_000195955.2 -> (GCF_000195955, 2)
_ACCESSION_VERSION_PATTERN = re.compile(r"^(GC[AF]_\d+)\.(\d+)$")


@dataclass(frozen=True)
class DownloadOutcome:
    """Outcome of a single genome download attempt."""

    accession: str
    success: bool
    resolved_version: str | None = None
    reason: str = ""

    @property
    def was_version_recovered(self) -> bool:
        """True if downloaded via an alternative accession version."""
        return self.success and self.resolved_version is not None


@dataclass
class DownloadReport:
    """Structured report of genome download outcomes."""

    outcomes: list[DownloadOutcome] = field(default_factory=list)
    elapsed_seconds: float = 0.0

    @property
    def succeeded(self) -> list[DownloadOutcome]:
        return [o for o in self.outcomes if o.success]

    @property
    def failed(self) -> list[DownloadOutcome]:
        return [o for o in self.outcomes if not o.success]

    @property
    def recovered(self) -> list[DownloadOutcome]:
        """Accessions recovered via alternative version."""
        return [o for o in self.outcomes if o.was_version_recovered]

    def write_failures_tsv(self, path: Path) -> None:
        """Write failed accessions to a TSV file."""
        failed = self.failed
        if not failed:
            return
        lines = ["accession\treason"]
        for o in failed:
            lines.append(f"{o.accession}\t{o.reason}")
        path.write_text("\n".join(lines) + "\n")


def _parse_accession_version(accession: str) -> tuple[str, int] | None:
    """Split accession into base and version number.

    Args:
        accession: e.g. "GCF_000195955.2"

    Returns:
        (base, version) tuple, e.g. ("GCF_000195955", 2), or None.
    """
    m = _ACCESSION_VERSION_PATTERN.match(accession)
    if m:
        return m.group(1), int(m.group(2))
    return None


def _alternative_versions(accession: str, max_attempts: int = 3) -> list[str]:
    """Generate alternative accession versions to try.

    Interleaves higher and lower versions: N+1, N-1, N+2, N-2, ...
    Higher versions are tried first at each step (more likely for NCBI
    updates). Returns up to max_attempts alternatives.
    """
    parsed = _parse_accession_version(accession)
    if parsed is None:
        return []
    base, version = parsed
    candidates: list[str] = []
    for offset in range(1, max_attempts + 1):
        # Try higher version
        candidates.append(f"{base}.{version + offset}")
        # Try lower version
        v = version - offset
        if v >= 1:
            candidates.append(f"{base}.{v}")
    return candidates[:max_attempts]


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
        retry_versions: bool = True,
        max_version_attempts: int = 3,
        progress_callback: Callable[[int, int], None] | None = None,
    ) -> DownloadReport:
        """Download genomes in batches and extract FASTA files.

        Downloads genomes from NCBI in batches, extracts the FASTA files,
        and optionally decompresses them. Tracks per-accession outcomes and
        retries failed accessions with alternative version numbers.

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
            retry_versions: If True, retry failed accessions with alt versions.
            max_version_attempts: Max alternative versions to try per accession.
            progress_callback: Optional callback(completed, total).

        Returns:
            DownloadReport with per-accession outcomes.
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        # Filter accessions if skip_if_exists
        if skip_if_exists:
            accessions_to_download = []
            for acc in accessions:
                fasta_path = output_dir / f"{acc}.fna"
                fasta_gz_path = output_dir / f"{acc}.fna.gz"
                if not fasta_path.exists() and not fasta_gz_path.exists():
                    accessions_to_download.append(acc)
            accessions = accessions_to_download

        if not accessions:
            return DownloadReport(outcomes=[], elapsed_seconds=0.0)

        total = len(accessions)
        completed = 0
        total_elapsed = 0.0
        all_failed: list[str] = []
        outcomes: list[DownloadOutcome] = []

        # Process in batches
        for i in range(0, total, batch_size):
            batch = accessions[i : i + batch_size]
            batch_zip = output_dir / f"ncbi_dataset_batch_{i}.zip"

            if dry_run:
                for acc in batch:
                    outcomes.append(DownloadOutcome(accession=acc, success=True))
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

            total_elapsed += result.elapsed_seconds

            extracted_accessions: set[str] = set()
            if result.success and batch_zip.exists():
                _, extracted_accessions = self._extract_fastas(
                    batch_zip,
                    output_dir,
                    decompress=decompress,
                    include_protein=include_protein,
                )
                batch_zip.unlink()

            # Track per-accession outcomes for this batch
            batch_set = set(batch)
            succeeded_in_batch = batch_set & extracted_accessions
            failed_in_batch = batch_set - extracted_accessions

            for acc in succeeded_in_batch:
                outcomes.append(DownloadOutcome(accession=acc, success=True))
            all_failed.extend(sorted(failed_in_batch))

            completed += len(batch)
            if progress_callback:
                progress_callback(completed, total)

            if batch_delay > 0 and i + batch_size < total:
                time.sleep(batch_delay)

        # Retry failed accessions with alternative versions
        if retry_versions and all_failed and not dry_run:
            retry_outcomes = self._retry_with_alternative_versions(
                all_failed,
                output_dir,
                max_version_attempts=max_version_attempts,
                decompress=decompress,
                include_protein=include_protein,
                timeout=timeout_per_batch,
            )
            total_elapsed += sum(
                o.elapsed_seconds
                for o in [retry_outcomes]
                if hasattr(o, "elapsed_seconds")
            )
            outcomes.extend(retry_outcomes)
        else:
            # No retry — mark all failed as final failures
            for acc in all_failed:
                outcomes.append(DownloadOutcome(
                    accession=acc, success=False, reason="not found in NCBI",
                ))

        return DownloadReport(outcomes=outcomes, elapsed_seconds=total_elapsed)

    def _retry_with_alternative_versions(
        self,
        failed_accessions: list[str],
        output_dir: Path,
        *,
        max_version_attempts: int = 3,
        decompress: bool = True,
        include_protein: bool = False,
        timeout: float = 1800.0,
    ) -> list[DownloadOutcome]:
        """Retry failed accessions by trying alternative version numbers.

        For each failed accession (e.g. GCA_012345678.1), tries alternative
        versions (.2, .3, etc.) and saves under the original accession name.

        Args:
            failed_accessions: Accessions that failed initial download.
            output_dir: Output directory for FASTA files.
            max_version_attempts: Max alternative versions to try.
            decompress: If True, decompress .gz files.
            include_protein: If True, also extract protein files.
            timeout: Download timeout in seconds.

        Returns:
            List of outcomes for each failed accession.
        """
        if not failed_accessions:
            return []

        # Build candidate -> original mapping for all failed accessions
        # candidate_to_original: maps alt version -> original accession
        candidate_to_original: dict[str, str] = {}
        # Track which originals have no candidates (no version to parse)
        no_candidates: list[str] = []

        for acc in failed_accessions:
            candidates = _alternative_versions(acc, max_attempts=max_version_attempts)
            if not candidates:
                no_candidates.append(acc)
                continue
            for candidate in candidates:
                # First candidate mapping wins (don't overwrite)
                if candidate not in candidate_to_original:
                    candidate_to_original[candidate] = acc

        outcomes: list[DownloadOutcome] = []

        # Mark accessions with no parseable version as final failures
        for acc in no_candidates:
            outcomes.append(DownloadOutcome(
                accession=acc, success=False,
                reason="no version number to retry",
            ))

        if not candidate_to_original:
            return outcomes

        all_candidates = list(candidate_to_original.keys())
        logger.info(
            "Retrying %d failed accessions with %d alternative versions",
            len(failed_accessions) - len(no_candidates),
            len(all_candidates),
        )

        # Download all candidates in one batch
        batch_zip = output_dir / "ncbi_dataset_version_retry.zip"
        result = self.run(
            accessions=all_candidates,
            output_file=batch_zip,
            include_sequence=True,
            include_protein=include_protein,
            timeout=timeout,
        )

        recovered_originals: set[str] = set()
        if result.success and batch_zip.exists():
            # Build rename map: NCBI accession -> original accession
            rename_map = {
                cand: orig
                for cand, orig in candidate_to_original.items()
            }
            _, extracted_accessions = self._extract_fastas(
                batch_zip,
                output_dir,
                decompress=decompress,
                include_protein=include_protein,
                rename_map=rename_map,
            )
            batch_zip.unlink()

            # Determine which originals were recovered
            for extracted_acc in extracted_accessions:
                original = candidate_to_original.get(extracted_acc)
                if original and original not in recovered_originals:
                    recovered_originals.add(original)
                    outcomes.append(DownloadOutcome(
                        accession=original,
                        success=True,
                        resolved_version=extracted_acc,
                    ))
                    logger.info(
                        "Recovered %s via alternative version %s",
                        original, extracted_acc,
                    )
        elif batch_zip.exists():
            batch_zip.unlink(missing_ok=True)

        # Mark remaining failures
        all_original_failed = set(failed_accessions) - recovered_originals - set(no_candidates)
        for acc in sorted(all_original_failed):
            candidates = _alternative_versions(acc, max_attempts=max_version_attempts)
            tried_str = ", ".join([acc] + candidates)
            outcomes.append(DownloadOutcome(
                accession=acc, success=False,
                reason=f"not found in NCBI (tried: {tried_str})",
            ))

        return outcomes

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
        rename_map: dict[str, str] | None = None,
    ) -> tuple[list[Path], set[str]]:
        """Extract FASTA files from NCBI dataset zip.

        Args:
            zip_path: Path to downloaded zip file
            output_dir: Output directory for FASTA files
            decompress: If True, decompress .gz files
            include_protein: If True, also extract protein files (.faa)
            rename_map: Optional mapping of NCBI accession -> desired output
                accession. Used when downloading alternative versions to save
                files under the originally requested accession name.

        Returns:
            Tuple of (extracted file paths, set of accessions extracted).
            When rename_map is provided, the returned accession set contains
            the *original* (pre-rename) accessions from NCBI.
        """
        extracted = []
        extracted_accessions: set[str] = set()

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

                # Determine output accession (may be renamed for version recovery)
                output_accession = accession
                if rename_map and accession in rename_map:
                    output_accession = rename_map[accession]

                # Determine output filename based on file type
                if is_protein:
                    if name.endswith(".faa.gz"):
                        out_name = f"{output_accession}.faa.gz"
                    else:
                        out_name = f"{output_accession}.faa"
                else:
                    if name.endswith(".fna.gz"):
                        out_name = f"{output_accession}.fna.gz"
                    else:
                        out_name = f"{output_accession}.fna"

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
                extracted_accessions.add(accession)

        return extracted, extracted_accessions


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
