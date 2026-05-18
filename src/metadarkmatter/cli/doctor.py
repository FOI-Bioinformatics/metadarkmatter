"""
Doctor command: report the operational state of a metadarkmatter install.

Verifies the Python interpreter, package version, key Python dependencies,
and the discoverability and reported versions of the external bioinformatics
tools that metadarkmatter shells out to. Intended as a first stop when a
user reports a problem in an internal-lab context: the output is plain text
and meant to be pasted directly into a bug report.
"""

from __future__ import annotations

import platform
import shutil
import subprocess
import sys
from importlib import metadata as importlib_metadata

import typer
from rich.console import Console
from rich.table import Table

from metadarkmatter import __version__

app = typer.Typer(
    name="doctor",
    help="Report environment, dependencies, and external tool versions.",
    no_args_is_help=False,
)

console = Console()


# (executable name, single-shot version flag) for each external tool.
# Many of these tools print versions to stderr and exit non-zero, so we
# capture both streams and tolerate non-zero exit codes.
_EXTERNAL_TOOLS: tuple[tuple[str, list[str]], ...] = (
    ("blastn", ["-version"]),
    ("blastx", ["-version"]),
    ("makeblastdb", ["-version"]),
    ("mmseqs", ["version"]),
    ("kraken2", ["--version"]),
    ("fastANI", ["--version"]),
    ("skani", ["--version"]),
    ("diamond", ["version"]),
    ("mashtree", ["--version"]),
    ("mash", ["--version"]),
)

# Python packages reported when present; absence is not an error here.
_PYTHON_PACKAGES: tuple[str, ...] = (
    "polars",
    "pyarrow",
    "numpy",
    "scipy",
    "biopython",
    "pydantic",
    "pydantic-settings",
    "typer",
    "rich",
    "pyyaml",
    "httpx",
    "plotly",
    "scikit-learn",
    "matplotlib",
)


def _probe_tool(executable: str, version_args: list[str]) -> tuple[str, str]:
    """Return (path, version_string) for the given external tool.

    Returns ("", "not found") when the tool is not on PATH. The version string
    is the first non-empty line of combined stdout/stderr, truncated to keep
    the report scannable.
    """
    path = shutil.which(executable)
    if path is None:
        return ("", "not found")

    try:
        result = subprocess.run(  # noqa: S603 - args are a constant list
            [path, *version_args],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except (subprocess.TimeoutExpired, OSError) as exc:
        return (path, f"probe failed: {exc}")

    combined = (result.stdout or "") + (result.stderr or "")
    first_line = next((line.strip() for line in combined.splitlines() if line.strip()), "")
    if not first_line:
        first_line = f"exit {result.returncode}"
    if len(first_line) > 120:
        first_line = first_line[:117] + "..."
    return (path, first_line)


def _probe_package(name: str) -> str:
    try:
        return importlib_metadata.version(name)
    except importlib_metadata.PackageNotFoundError:
        return "not installed"


@app.callback(invoke_without_command=True)
def doctor() -> None:
    """Print a diagnostic report of the metadarkmatter environment."""
    console.print(f"[bold]metadarkmatter[/bold] {__version__}")
    console.print(f"Python {sys.version.split()[0]} ({sys.executable})")
    console.print(f"Platform: {platform.platform()}")
    console.print()

    pkg_table = Table(title="Python dependencies", show_header=True, header_style="bold")
    pkg_table.add_column("Package")
    pkg_table.add_column("Version")
    for pkg in _PYTHON_PACKAGES:
        pkg_table.add_row(pkg, _probe_package(pkg))
    console.print(pkg_table)
    console.print()

    tool_table = Table(title="External tools", show_header=True, header_style="bold")
    tool_table.add_column("Tool")
    tool_table.add_column("Path")
    tool_table.add_column("Version")
    missing: list[str] = []
    for exe, args in _EXTERNAL_TOOLS:
        path, version = _probe_tool(exe, args)
        if not path:
            missing.append(exe)
            tool_table.add_row(exe, "-", "[yellow]not found[/yellow]")
        else:
            tool_table.add_row(exe, path, version)
    console.print(tool_table)

    if missing:
        console.print(
            "\n[yellow]"
            f"{len(missing)} tool(s) not on PATH: {', '.join(missing)}. "
            "Install only the tools your workflow needs."
            "[/yellow]"
        )
