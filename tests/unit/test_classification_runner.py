"""Tests for the CLI-free classification runner.

These exercise ``core.classification.runner`` directly (no Typer/Rich harness),
which is the whole point of the extraction: the load -> classify -> finalize
pipeline must be drivable in-process and unit-testable on its own.
"""

from __future__ import annotations

import builtins
from pathlib import Path
from unittest.mock import patch

import polars as pl
import pytest

from metadarkmatter.core.ani_placement import ANIMatrix
from metadarkmatter.core.classification.runner import (
    ClassificationRequest,
    finalize_classification,
    run_classification,
    validate_ani_genome_coverage,
)
from metadarkmatter.core.exceptions import ConfigurationError
from metadarkmatter.models.config import ScoringConfig

# The coverage validator extracts the accession from the standardized
# ``{accession}|{contig}`` FASTA-header format, so these tests use a BLAST
# fixture written in that form (the shared conftest fixture predates it).
_GENOMES = ["GCF_000123456.1", "GCF_000789012.1", "GCA_000111222.1"]


@pytest.fixture
def pipe_blast_file(temp_dir: Path) -> Path:
    """BLAST tabular file whose sseqids use the ``{accession}|{contig}`` format."""
    rows = []
    for read_idx in range(5):
        for g_idx, genome in enumerate(_GENOMES):
            rows.append(
                {
                    "qseqid": f"read_{read_idx:03d}",
                    "sseqid": f"{genome}|contig{g_idx}",
                    "pident": 98.5 - g_idx * 5,
                    "length": 150,
                    "mismatch": 2 + g_idx * 3,
                    "gapopen": 0,
                    "qstart": 1,
                    "qend": 150,
                    "sstart": 1000,
                    "send": 1150,
                    "evalue": 1e-50,
                    "bitscore": 250.0 - g_idx * 20,
                }
            )
    blast_path = temp_dir / "pipe.blast.tsv"
    pl.DataFrame(rows).write_csv(blast_path, separator="\t", include_header=False)
    return blast_path


def _request(blast_file: Path, temp_ani_file: Path, output: Path, **overrides):
    """Build a ClassificationRequest with sensible defaults for these tests."""
    fields: dict = {
        "alignment": blast_file,
        "ani": temp_ani_file,
        "output": output,
        "config": ScoringConfig(),
    }
    fields.update(overrides)
    return ClassificationRequest(**fields)


# ---------------------------------------------------------------------------
# validate_ani_genome_coverage
# ---------------------------------------------------------------------------


def test_validate_coverage_full(pipe_blast_file: Path, temp_ani_file: Path) -> None:
    """All BLAST genomes present in the ANI matrix => 100% coverage, no misses."""
    ani = ANIMatrix.from_file(temp_ani_file)
    matched, total, pct, missing = validate_ani_genome_coverage(pipe_blast_file, ani)

    assert total == 3
    assert matched == 3
    assert pct == pytest.approx(100.0)
    assert missing == set()


def test_validate_coverage_with_missing_genome(
    pipe_blast_file: Path, temp_dir: Path
) -> None:
    """A genome absent from the ANI matrix is reported as missing."""
    # ANI matrix covering only two of the three BLAST genomes.
    genomes = ["GCF_000123456.1", "GCF_000789012.1"]
    data: dict[str, list] = {"genome": genomes}
    data["GCF_000123456.1"] = [100.0, 95.5]
    data["GCF_000789012.1"] = [95.5, 100.0]
    ani_path = temp_dir / "partial.ani.csv"
    pl.DataFrame(data).write_csv(ani_path)

    ani = ANIMatrix.from_file(ani_path)
    matched, total, pct, missing = validate_ani_genome_coverage(pipe_blast_file, ani)

    assert total == 3
    assert matched == 2
    assert pct == pytest.approx(2 / 3 * 100.0)
    assert missing == {"GCA_000111222.1"}


def test_validate_coverage_representative_mapping_scopes_denominator(
    pipe_blast_file: Path, temp_ani_file: Path
) -> None:
    """With a representative mapping, only mapped (family) genomes are counted."""
    ani = ANIMatrix.from_file(temp_ani_file)
    # Only one BLAST genome is in the mapping; the denominator collapses to it.
    rep_map = {"GCF_000123456.1": "GCF_000123456.1"}
    matched, total, pct, missing = validate_ani_genome_coverage(
        pipe_blast_file, ani, representative_mapping=rep_map
    )

    assert total == 1
    assert matched == 1
    assert pct == pytest.approx(100.0)
    assert missing == set()


# ---------------------------------------------------------------------------
# finalize_classification
# ---------------------------------------------------------------------------


def test_finalize_none_returns_zero(temp_dir: Path) -> None:
    """A None DataFrame (e.g. streaming) yields zero and writes nothing."""
    output = temp_dir / "out.csv"
    assert finalize_classification(None, None, output, "csv") == 0
    assert not output.exists()


def test_finalize_writes_rows(temp_dir: Path) -> None:
    """A non-empty DataFrame is written and its row count returned."""
    output = temp_dir / "out.csv"
    df = pl.DataFrame({"read_id": ["r1", "r2"], "taxonomic_call": ["Known", "Novel"]})

    assert finalize_classification(df, None, output, "csv") == 2
    assert output.exists()
    assert len(pl.read_csv(output)) == 2


def test_finalize_empty_skips_write(temp_dir: Path) -> None:
    """An empty DataFrame returns zero and does not create the output file."""
    output = temp_dir / "out.csv"
    df = pl.DataFrame({"read_id": [], "taxonomic_call": []})

    assert finalize_classification(df, None, output, "csv") == 0
    assert not output.exists()


def test_finalize_writes_parquet(temp_dir: Path) -> None:
    """The requested parquet output format is honoured."""
    output = temp_dir / "out.parquet"
    df = pl.DataFrame({"read_id": ["r1"], "taxonomic_call": ["Known"]})

    assert finalize_classification(df, None, output, "parquet") == 1
    assert output.exists()
    assert len(pl.read_parquet(output)) == 1


def test_finalize_joins_metadata(temp_dir: Path) -> None:
    """When metadata is supplied, its join_classifications hook is applied."""
    from unittest.mock import MagicMock

    output = temp_dir / "out.csv"
    df = pl.DataFrame({"read_id": ["r1"], "taxonomic_call": ["Known"]})
    metadata = MagicMock()
    metadata.join_classifications.return_value = df

    finalize_classification(df, metadata, output, "csv")
    metadata.join_classifications.assert_called_once_with(df)


# ---------------------------------------------------------------------------
# run_classification (end-to-end, non-streaming)
# ---------------------------------------------------------------------------


def test_run_classification_basic(
    pipe_blast_file: Path, temp_ani_file: Path, temp_dir: Path
) -> None:
    """Non-streaming run classifies reads, writes output, and reports coverage."""
    output = temp_dir / "classifications.csv"
    result = run_classification(_request(pipe_blast_file, temp_ani_file, output))

    assert result.num_classified > 0
    assert result.classification_df is not None
    assert len(result.classification_df) == result.num_classified
    assert result.coverage is not None
    assert result.coverage[2] == pytest.approx(100.0)
    assert output.exists()


def test_run_classification_invokes_log_callback(
    pipe_blast_file: Path, temp_ani_file: Path, temp_dir: Path
) -> None:
    """Status text is delivered through the injected log callback, not stdout."""
    output = temp_dir / "classifications.csv"
    messages: list[str] = []

    run_classification(
        _request(pipe_blast_file, temp_ani_file, output), log=messages.append
    )

    joined = "\n".join(messages)
    assert "Loaded ANI matrix" in joined
    assert "Classified" in joined


def test_run_classification_streaming_writes_without_dataframe(
    pipe_blast_file: Path, temp_ani_file: Path, temp_dir: Path
) -> None:
    """Streaming mode writes incrementally and returns no in-memory DataFrame."""
    output = temp_dir / "streamed.csv"
    result = run_classification(
        _request(pipe_blast_file, temp_ani_file, output, streaming=True)
    )

    assert result.classification_df is None
    assert result.num_classified > 0
    assert output.exists()


def test_run_classification_rejects_genomes_and_id_mapping(
    pipe_blast_file: Path, temp_ani_file: Path, temp_dir: Path
) -> None:
    """Supplying both --genomes and --id-mapping is a configuration error."""
    output = temp_dir / "out.csv"
    request = _request(
        pipe_blast_file,
        temp_ani_file,
        output,
        genomes=temp_dir,
        id_mapping=temp_dir / "mapping.tsv",
    )

    with pytest.raises(ConfigurationError, match="both --genomes and --id-mapping"):
        run_classification(request)


def test_run_classification_adaptive_missing_dependency(
    pipe_blast_file: Path, temp_ani_file: Path, temp_dir: Path
) -> None:
    """A missing scikit-learn for adaptive thresholds raises a typed error."""
    output = temp_dir / "out.csv"
    request = _request(
        pipe_blast_file, temp_ani_file, output, adaptive_thresholds=True
    )

    real_import = builtins.__import__

    def fake_import(name: str, *args, **kwargs):
        if "adaptive" in name:
            raise ImportError("No module named 'sklearn'")
        return real_import(name, *args, **kwargs)

    with patch("builtins.__import__", side_effect=fake_import):
        with pytest.raises(ConfigurationError, match="scikit-learn"):
            run_classification(request)
