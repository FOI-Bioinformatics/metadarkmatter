"""
E2E test fixtures for metadarkmatter CLI testing.

Provides fixtures that combine test factories with CLI invocation
helpers for comprehensive end-to-end testing.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable
from unittest.mock import MagicMock

import pytest
from typer.testing import CliRunner

from metadarkmatter.cli.main import app
from tests.factories import E2ETestDataset

if TYPE_CHECKING:
    from click.testing import Result


# =============================================================================
# Fixture Directory Path
# =============================================================================

FIXTURES_DIR = Path(__file__).parent.parent / "fixtures"


# =============================================================================
# Core Fixtures
# =============================================================================


@pytest.fixture
def e2e_runner() -> CliRunner:
    """Provide a CLI runner for E2E tests."""
    return CliRunner()


@pytest.fixture
def e2e_temp_dir(tmp_path: Path) -> Path:
    """Provide a temporary directory for E2E test files."""
    return tmp_path


@pytest.fixture
def test_dataset(e2e_temp_dir: Path) -> E2ETestDataset:
    """Provide a seeded test dataset generator."""
    return E2ETestDataset(e2e_temp_dir, seed=42)


# =============================================================================
# Standard Dataset Fixtures
# =============================================================================


@pytest.fixture
def standard_dataset(test_dataset: E2ETestDataset) -> tuple[Path, Path]:
    """Create a standard mixed dataset for testing."""
    return test_dataset.create_standard_dataset()


@pytest.fixture
def batch_dataset(test_dataset: E2ETestDataset) -> tuple[Path, Path]:
    """Create a batch dataset (multiple BLAST files) for testing."""
    return test_dataset.create_batch_dataset(n_samples=3)


@pytest.fixture
def single_read_dataset(test_dataset: E2ETestDataset) -> tuple[Path, Path]:
    """Create a minimal single-read dataset for edge case testing."""
    return test_dataset.create_single_read_dataset()


@pytest.fixture
def boundary_dataset(test_dataset: E2ETestDataset) -> tuple[Path, Path]:
    """Create a dataset with reads at classification boundary values."""
    return test_dataset.create_boundary_novelty_dataset()


@pytest.fixture
def empty_blast_file(test_dataset: E2ETestDataset) -> Path:
    """Create an empty BLAST file for edge case testing."""
    return test_dataset.create_empty_blast_file()


# =============================================================================
# CLI Invocation Helpers
# =============================================================================


@pytest.fixture
def run_classify(
    e2e_runner: CliRunner,
    e2e_temp_dir: Path,
) -> Callable[..., Result]:
    """
    Fixture that returns a function to run the classify command.

    Usage:
        result = run_classify(
            blast=blast_path,
            ani=ani_path,
            output_name="result.csv",
            fast=True,
        )
    """
    def _run(
        blast: Path,
        ani: Path,
        output_name: str = "classifications.csv",
        summary_name: str | None = None,
        output_format: str = "csv",
        bitscore_threshold: float = 95.0,
        fast: bool = False,
        parallel: bool = False,
        streaming: bool = False,
        verbose: bool = False,
        extra_args: list[str] | None = None,
    ) -> Result:
        output_path = e2e_temp_dir / output_name

        args = [
            "score", "classify",
            "--blast", str(blast),
            "--ani", str(ani),
            "--output", str(output_path),
            "--format", output_format,
            "--bitscore-threshold", str(bitscore_threshold),
        ]

        if summary_name:
            args.extend(["--summary", str(e2e_temp_dir / summary_name)])

        if fast:
            args.append("--fast")
        if parallel:
            args.append("--parallel")
        if streaming:
            args.append("--streaming")
        if verbose:
            args.append("--verbose")

        if extra_args:
            args.extend(extra_args)

        return e2e_runner.invoke(app, args)

    return _run


@pytest.fixture
def run_batch(
    e2e_runner: CliRunner,
    e2e_temp_dir: Path,
) -> Callable[..., Result]:
    """
    Fixture that returns a function to run the batch command.

    Usage:
        result = run_batch(
            blast_dir=blast_dir,
            ani=ani_path,
            pattern="*.blast.tsv",
        )
    """
    def _run(
        blast_dir: Path,
        ani: Path,
        output_dir_name: str = "batch_output",
        pattern: str = "*.blast.tsv",
        output_format: str = "csv",
        bitscore_threshold: float = 95.0,
        fast: bool = False,
        parallel: bool = False,
        verbose: bool = False,
        extra_args: list[str] | None = None,
    ) -> Result:
        output_dir = e2e_temp_dir / output_dir_name

        args = [
            "score", "batch",
            "--blast-dir", str(blast_dir),
            "--ani", str(ani),
            "--output-dir", str(output_dir),
            "--pattern", pattern,
            "--format", output_format,
            "--bitscore-threshold", str(bitscore_threshold),
        ]

        if fast:
            args.append("--fast")
        if parallel:
            args.append("--parallel")
        if verbose:
            args.append("--verbose")

        if extra_args:
            args.extend(extra_args)

        return e2e_runner.invoke(app, args)

    return _run


# =============================================================================
# Output Path Helpers
# =============================================================================


@pytest.fixture
def output_paths(e2e_temp_dir: Path) -> dict[str, Path]:
    """
    Provide standard output paths for tests.

    Returns dict with keys:
    - classification_csv
    - classification_parquet
    - summary_json
    - batch_dir
    """
    return {
        "classification_csv": e2e_temp_dir / "classifications.csv",
        "classification_parquet": e2e_temp_dir / "classifications.parquet",
        "summary_json": e2e_temp_dir / "summary.json",
        "batch_dir": e2e_temp_dir / "batch_output",
    }


# =============================================================================
# Validation Fixtures
# =============================================================================


@pytest.fixture
def validate_output() -> Callable[[Path], Any]:
    """
    Fixture that returns a validation function for classification outputs.

    Combines file loading with comprehensive validation.
    """
    from tests.utils.assertions import validate_classification_output
    return validate_classification_output


@pytest.fixture
def validate_summary() -> Callable[[Path, Any], Any]:
    """
    Fixture that returns a validation function for summary outputs.
    """
    from tests.utils.assertions import validate_summary_output
    return validate_summary_output


# =============================================================================
# Download Command Fixtures
# =============================================================================


@pytest.fixture
def gtdb_francisella_fixture() -> dict[str, Any]:
    """Load saved GTDB API response for g__Francisella."""
    fixture_path = FIXTURES_DIR / "gtdb_francisella_response.json"
    with open(fixture_path) as f:
        return json.load(f)


@pytest.fixture
def gtdb_francisellaceae_fixture() -> dict[str, Any]:
    """Load saved GTDB API response for f__Francisellaceae."""
    fixture_path = FIXTURES_DIR / "gtdb_francisellaceae_response.json"
    with open(fixture_path) as f:
        return json.load(f)


@pytest.fixture
def mock_gtdb_response(monkeypatch: pytest.MonkeyPatch, gtdb_francisella_fixture: dict) -> dict:
    """Mock GTDB API with Francisella response for genus queries."""
    import httpx

    def mock_get(self: Any, endpoint: str, params: dict | None = None) -> MagicMock:
        """Mock httpx.Client.get to return fixture data."""
        response = MagicMock()
        response.status_code = 200
        response.json.return_value = gtdb_francisella_fixture
        response.raise_for_status = MagicMock()
        return response

    monkeypatch.setattr(httpx.Client, "get", mock_get)
    return gtdb_francisella_fixture


@pytest.fixture
def mock_gtdb_family_response(
    monkeypatch: pytest.MonkeyPatch, gtdb_francisellaceae_fixture: dict
) -> dict:
    """Mock GTDB API with Francisellaceae response for family queries."""
    import httpx

    def mock_get(self: Any, endpoint: str, params: dict | None = None) -> MagicMock:
        """Mock httpx.Client.get to return fixture data."""
        response = MagicMock()
        response.status_code = 200
        response.json.return_value = gtdb_francisellaceae_fixture
        response.raise_for_status = MagicMock()
        return response

    monkeypatch.setattr(httpx.Client, "get", mock_get)
    return gtdb_francisellaceae_fixture


@pytest.fixture
def mock_gtdb_404(monkeypatch: pytest.MonkeyPatch) -> None:
    """Mock GTDB API to return 404 for unknown taxa."""
    import httpx

    def mock_get(self: Any, endpoint: str, params: dict | None = None) -> MagicMock:
        """Mock httpx.Client.get to raise 404."""
        response = MagicMock()
        response.status_code = 404

        def raise_status() -> None:
            raise httpx.HTTPStatusError(
                "Not Found",
                request=MagicMock(),
                response=response,
            )

        response.raise_for_status = raise_status
        return response

    monkeypatch.setattr(httpx.Client, "get", mock_get)


@pytest.fixture
def sample_accession_tsv(tmp_path: Path) -> Path:
    """Create sample accession TSV file for fetch command testing."""
    tsv_content = """accession\tgtdb_taxonomy\tspecies
GCF_000833475.1\td__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Francisellales;f__Francisellaceae;g__Francisella;s__Francisella tularensis\tFrancisella tularensis
GCF_000156715.1\td__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Francisellales;f__Francisellaceae;g__Francisella;s__Francisella philomiragia\tFrancisella philomiragia
GCF_001885235.1\td__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Francisellales;f__Francisellaceae;g__Francisella;s__Francisella hispaniensis\tFrancisella hispaniensis
"""
    tsv_path = tmp_path / "francisella_accessions.tsv"
    tsv_path.write_text(tsv_content)
    return tsv_path


@pytest.fixture
def mock_ncbi_missing(monkeypatch: pytest.MonkeyPatch) -> None:
    """Mock NCBIDatasets as unavailable (executable not found)."""
    from metadarkmatter.external.ncbi_datasets import NCBIDatasets

    def mock_check_available(self: Any) -> bool:
        return False

    monkeypatch.setattr(NCBIDatasets, "check_available", mock_check_available)
    # Clear the executable cache
    NCBIDatasets._executable_cache.clear()


@pytest.fixture
def mock_ncbi_available(monkeypatch: pytest.MonkeyPatch) -> None:
    """Mock NCBIDatasets as available."""
    from metadarkmatter.external.ncbi_datasets import NCBIDatasets

    NCBIDatasets._executable_cache["datasets"] = Path("/usr/bin/datasets")


@pytest.fixture
def run_download_list(e2e_runner: CliRunner, e2e_temp_dir: Path) -> Callable[..., "Result"]:
    """Fixture to run download genomes list command."""

    def _run(
        taxon: str,
        output_name: str = "accessions.tsv",
        all_genomes: bool = False,
        include_metadata: bool = False,
        verbose: bool = False,
        quiet: bool = False,
        dry_run: bool = False,
        extra_args: list[str] | None = None,
    ) -> "Result":
        output_path = e2e_temp_dir / output_name

        args = [
            "download", "genomes", "list",
            taxon,
            "--output", str(output_path),
        ]

        if all_genomes:
            args.append("--all-genomes")
        if include_metadata:
            args.append("--include-metadata")
        if verbose:
            args.append("--verbose")
        if quiet:
            args.append("--quiet")
        if dry_run:
            args.append("--dry-run")

        if extra_args:
            args.extend(extra_args)

        return e2e_runner.invoke(app, args)

    return _run


@pytest.fixture
def run_download_fetch(e2e_runner: CliRunner, e2e_temp_dir: Path) -> Callable[..., "Result"]:
    """Fixture to run download genomes fetch command."""

    def _run(
        accession_file: Path,
        output_dir_name: str = "genomes",
        batch_size: int = 100,
        skip_if_exists: bool = True,  # Default is True in CLI
        redownload: bool = False,
        verbose: bool = False,
        quiet: bool = False,
        dry_run: bool = False,
        extra_args: list[str] | None = None,
    ) -> "Result":
        output_dir = e2e_temp_dir / output_dir_name

        args = [
            "download", "genomes", "fetch",
            "--accessions", str(accession_file),
            "--output-dir", str(output_dir),
            "--batch-size", str(batch_size),
        ]

        if redownload:
            args.append("--redownload")
        if verbose:
            args.append("--verbose")
        if quiet:
            args.append("--quiet")
        if dry_run:
            args.append("--dry-run")

        if extra_args:
            args.extend(extra_args)

        return e2e_runner.invoke(app, args)

    return _run
