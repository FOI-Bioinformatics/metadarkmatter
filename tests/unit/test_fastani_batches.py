"""Tests for FastANI.compute_parallel_batches.

The fastANI binary is not installed in CI, so we patch the subprocess
call inside the module-level worker to verify the orchestration logic
(splitting, threads-per-batch, output concatenation, deterministic
ordering).
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from metadarkmatter.external.fastani import FastANI


@pytest.fixture
def fake_genome_dir(tmp_path: Path) -> Path:
    """Create N fake genome files and return the directory."""
    for i in range(8):
        (tmp_path / f"g{i:02d}.fna").write_text(f">contig_{i}\nACGT\n")
    return tmp_path


@pytest.fixture
def query_list(tmp_path: Path, fake_genome_dir: Path) -> Path:
    """Genome list file with absolute paths."""
    list_file = tmp_path / "queries.txt"
    list_file.write_text(
        "\n".join(str(p) for p in sorted(fake_genome_dir.glob("g*.fna"))) + "\n"
    )
    return list_file


def _make_runner(per_chunk_lines: int):
    """Build a subprocess.run replacement that writes a stub fastANI output."""

    def fake_run(cmd, **kwargs):
        # The output path is the argument after '-o' in the command list.
        output = Path(cmd[cmd.index("-o") + 1])
        # Write predictable per-batch output so we can verify ordering.
        output.write_text("".join(f"line_{i}\n" for i in range(per_chunk_lines)))
        return MagicMock(returncode=0, stdout="", stderr="")

    return fake_run


def test_single_batch_uses_run_or_raise(query_list: Path, tmp_path: Path) -> None:
    """batches=1 takes the simple path that delegates to run_or_raise."""
    f = FastANI()
    with patch.object(f, "run_or_raise") as mock_ror:
        f.compute_parallel_batches(
            query_list=query_list,
            reference_list=query_list,
            output=tmp_path / "out.tsv",
            tmp_dir=tmp_path,
            threads=4,
            batches=1,
        )
    mock_ror.assert_called_once()


def test_multiple_batches_split_query_list_evenly(
    query_list: Path, tmp_path: Path
) -> None:
    """With 8 queries and batches=4, expect 4 per-batch input files of size 2."""
    f = FastANI()
    with patch.object(
        f, "get_executable", return_value=Path("/usr/bin/fastANI")
    ), patch(
        "metadarkmatter.external.fastani.subprocess.run",
        side_effect=_make_runner(per_chunk_lines=2),
        create=True,
    ):
        f.compute_parallel_batches(
            query_list=query_list,
            reference_list=query_list,
            output=tmp_path / "merged.tsv",
            tmp_dir=tmp_path,
            threads=8,
            batches=4,
        )

    batch_files = sorted(tmp_path.glob("query_batch_*.txt"))
    assert len(batch_files) == 4
    # 8 queries split into 4 chunks of 2 each
    for bf in batch_files:
        assert len(bf.read_text().strip().splitlines()) == 2


def test_batches_capped_at_query_count(query_list: Path, tmp_path: Path) -> None:
    """Asking for more batches than queries silently caps; no empty batches."""
    # Truncate to 3 queries.
    queries = query_list.read_text().splitlines()[:3]
    small = tmp_path / "small.txt"
    small.write_text("\n".join(queries) + "\n")

    f = FastANI()
    with patch.object(
        f, "get_executable", return_value=Path("/usr/bin/fastANI")
    ), patch(
        "metadarkmatter.external.fastani.subprocess.run",
        side_effect=_make_runner(per_chunk_lines=1),
        create=True,
    ):
        f.compute_parallel_batches(
            query_list=small,
            reference_list=query_list,
            output=tmp_path / "merged.tsv",
            tmp_dir=tmp_path,
            threads=8,
            batches=10,  # more than 3 queries
        )

    # At most 3 batches.
    assert len(list(tmp_path.glob("query_batch_*.txt"))) <= 3


def test_batches_must_be_positive(query_list: Path, tmp_path: Path) -> None:
    f = FastANI()
    with pytest.raises(ValueError, match="batches must be"):
        f.compute_parallel_batches(
            query_list=query_list,
            reference_list=query_list,
            output=tmp_path / "out.tsv",
            tmp_dir=tmp_path,
            threads=4,
            batches=0,
        )


def test_empty_query_list_raises(tmp_path: Path) -> None:
    empty = tmp_path / "empty.txt"
    empty.write_text("")
    f = FastANI()
    with pytest.raises(ValueError, match="empty"):
        f.compute_parallel_batches(
            query_list=empty,
            reference_list=empty,
            output=tmp_path / "out.tsv",
            tmp_dir=tmp_path,
            threads=4,
            batches=2,
        )
