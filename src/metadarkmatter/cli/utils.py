"""
Shared CLI utilities for metadarkmatter commands.

Provides common functionality used across CLI modules.
"""

from __future__ import annotations

from collections.abc import Generator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

from rich.console import Console
from rich.progress import (
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
)


@contextmanager
def spinner_progress(
    description: str,
    console: Console | None = None,
    quiet: bool = False,
) -> Generator[Progress, None, None]:
    """Context manager for spinner-style progress display.

    Creates a standardized spinner progress bar used throughout the CLI.
    The spinner is suppressed when quiet mode is enabled.

    Args:
        description: Task description to display.
        console: Rich Console instance. If None and not quiet, creates one.
        quiet: If True, suppress the progress display entirely.

    Yields:
        Progress instance (even when quiet, for API consistency).

    Example:
        >>> with spinner_progress("Loading data...", console, quiet) as progress:
        ...     # Perform work
        ...     pass
    """
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        TimeElapsedColumn(),
        console=console if not quiet else None,
    ) as progress:
        progress.add_task(description=description, total=None)
        yield progress


def extract_sample_name(path: Path, suffixes: tuple[str, ...] = ()) -> str:
    """Extract sample name from a file path.

    Removes common suffixes from the file stem to get a clean sample name.

    Args:
        path: Path to the file.
        suffixes: Additional suffixes to remove from the stem.

    Returns:
        Cleaned sample name.

    Example:
        >>> from pathlib import Path
        >>> extract_sample_name(Path("sample_001_classifications.csv"))
        'sample_001'
        >>> extract_sample_name(Path("sample.blast.tsv.gz"), (".blast", ".tsv"))
        'sample'
    """
    name = path.stem
    # Handle double extensions like .tsv.gz
    if path.suffix == ".gz" and "." in name:
        name = Path(name).stem

    # Remove common suffixes
    default_suffixes = (
        "_classifications",
        "classifications_",
        ".kraken",
        ".blast",
        ".tsv",
    )
    for suffix in default_suffixes + suffixes:
        name = name.replace(suffix, "")

    return name


def extract_sample_name_from_reads(reads_path: Path) -> str:
    """Extract sample name from a FASTQ reads file path.

    Handles common FASTQ naming conventions for paired-end sequencing reads,
    removing read pair indicators (_R1, _R2, _1, _2) and file extensions.

    Args:
        reads_path: Path to a reads file (FASTQ, optionally gzipped).

    Returns:
        Cleaned sample name without read pair indicators or extensions.

    Example:
        >>> from pathlib import Path
        >>> extract_sample_name_from_reads(Path("sample_001_R1.fastq.gz"))
        'sample_001'
        >>> extract_sample_name_from_reads(Path("sample_001_1.fq.gz"))
        'sample_001'
        >>> extract_sample_name_from_reads(Path("sample.fastq"))
        'sample'
    """
    name = reads_path.name

    # Remove common paired-end suffixes (ordered by specificity)
    read_suffixes = [
        "_R1.fastq.gz", "_R1.fq.gz", "_1.fastq.gz", "_1.fq.gz",
        "_R2.fastq.gz", "_R2.fq.gz", "_2.fastq.gz", "_2.fq.gz",
        "_R1.fastq", "_R1.fq", "_1.fastq", "_1.fq",
        "_R2.fastq", "_R2.fq", "_2.fastq", "_2.fq",
        ".fastq.gz", ".fq.gz", ".fastq", ".fq",
    ]

    for suffix in read_suffixes:
        if name.endswith(suffix):
            return name[:-len(suffix)]

    # Fallback: split on first dot
    return name.split(".")[0]


class QuietConsole:
    """Console wrapper that suppresses output in quiet mode.

    This class wraps a Rich Console instance and conditionally suppresses
    print output when quiet mode is enabled. All other console methods
    are delegated to the wrapped instance.

    Example:
        >>> console = Console()
        >>> qc = QuietConsole(console, quiet=True)
        >>> qc.print("This won't be shown")  # Suppressed
        >>> qc.print("This will be shown", force=True)  # Not suppressed
    """

    def __init__(self, console: Console, quiet: bool = False):
        """Initialize QuietConsole wrapper.

        Args:
            console: Rich Console instance to wrap.
            quiet: If True, suppress print output.
        """
        self._console = console
        self._quiet = quiet

    @property
    def console(self) -> Console:
        """Access the underlying Rich Console instance.

        Use this when you need to print regardless of quiet mode,
        such as for tables or other structured output.

        Returns:
            The wrapped Console instance.
        """
        return self._console

    def print(self, *args: Any, **kwargs: Any) -> None:
        """Print to console unless quiet mode is enabled.

        Args:
            *args: Positional arguments passed to Console.print.
            **kwargs: Keyword arguments passed to Console.print.
        """
        if not self._quiet:
            self._console.print(*args, **kwargs)

    def __getattr__(self, name: str) -> Any:
        """Delegate attribute access to wrapped console.

        Args:
            name: Attribute name to access.

        Returns:
            The attribute from the wrapped Console instance.
        """
        return getattr(self._console, name)
