"""
Doctor command: report the operational state of a metadarkmatter install.

Verifies the Python interpreter, package version, key Python dependencies,
and the discoverability and reported versions of the external bioinformatics
tools that metadarkmatter shells out to. Intended as a first stop when a
user reports a problem in an internal-lab context: the output is plain text
and meant to be pasted directly into a bug report.
"""

from __future__ import annotations

import os
import platform
import sys
from importlib import metadata as importlib_metadata
from pathlib import Path

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
    # Needed only by specific subcommands (see _TOOL_NEEDED_BY).
    ("bowtie2", ["--version"]),
    ("bowtie2-build", ["--version"]),
    ("samtools", ["--version"]),
    ("prodigal", ["-v"]),
)

# Which command(s) need each tool that is not part of the core pipeline.
_TOOL_NEEDED_BY: dict[str, str] = {
    "bowtie2": "optional: 'map' (read mapping)",
    "bowtie2-build": "optional: 'map' (read mapping)",
    "samtools": "optional: 'map', 'visualize recruitment'",
    "prodigal": "optional: 'proteins' (gene prediction)",
}

# Tools the default pipeline cannot run without. Each entry is a group of
# alternatives: at least one member of every group must be on PATH.
_REQUIRED_TOOL_GROUPS: tuple[tuple[str, ...], ...] = (
    ("blastn",),
    ("makeblastdb",),
    ("kraken2",),
    ("skani", "fastANI"),
)

# Python packages the CLI cannot start or run its core commands without.
_REQUIRED_PACKAGES: tuple[str, ...] = (
    "polars",
    "numpy",
    "pydantic",
    "typer",
    "rich",
    "pyyaml",
)

# bioconda install hint per external tool (the conda package name differs from
# the executable name in a few cases, and several executables share a package).
_TOOL_INSTALL_HINTS: dict[str, str] = {
    "blastn": "conda install -c bioconda blast",
    "blastx": "conda install -c bioconda blast",
    "makeblastdb": "conda install -c bioconda blast",
    "mmseqs": "conda install -c bioconda mmseqs2",
    "kraken2": "conda install -c bioconda kraken2",
    "fastANI": "conda install -c bioconda fastani",
    "skani": "conda install -c bioconda skani",
    "diamond": "conda install -c bioconda diamond",
    "mashtree": "conda install -c bioconda mashtree",
    "mash": "conda install -c bioconda mash",
    "bowtie2": "conda install -c bioconda bowtie2",
    "bowtie2-build": "conda install -c bioconda bowtie2",
    "samtools": "conda install -c bioconda samtools",
    "prodigal": "conda install -c bioconda prodigal",
}

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
    """Return (path, version_string) for the given external tool (see core.tool_versions)."""
    from metadarkmatter.core.tool_versions import probe_tool

    return probe_tool(executable, version_args)


def _probe_package(name: str) -> str:
    try:
        return importlib_metadata.version(name)
    except importlib_metadata.PackageNotFoundError:
        return "not installed"


_ENV_VARS = (
    ("METADARKMATTER_SEED", "global random seed (default 42)"),
    ("METADARKMATTER_GTDB_CACHE_DIR", "GTDB on-disk cache (empty = disabled)"),
)

# Files a built Kraken2 database directory must contain.
_KRAKEN2_DB_FILES = ("hash.k2d", "opts.k2d", "taxo.k2d")


def _check_kraken_db(db_dir: Path) -> tuple[bool, str]:
    """Validate a Kraken2 database directory.

    Returns (ok, detail). A valid DB is a directory containing the three
    .k2d index files Kraken2 produces.
    """
    if not db_dir.exists():
        return (False, f"path does not exist: {db_dir}")
    if not db_dir.is_dir():
        return (False, f"not a directory: {db_dir}")
    missing = [f for f in _KRAKEN2_DB_FILES if not (db_dir / f).is_file()]
    if missing:
        return (False, f"missing index files: {', '.join(missing)}")
    return (True, f"valid Kraken2 DB ({', '.join(_KRAKEN2_DB_FILES)} present)")


def _check_gtdb_cache() -> tuple[str, str]:
    """Report the GTDB cache directory state from the environment."""
    raw = os.environ.get("METADARKMATTER_GTDB_CACHE_DIR")
    if raw is None:
        return ("default", "~/.cache/metadarkmatter/gtdb (used if present)")
    if raw == "":
        return ("disabled", "METADARKMATTER_GTDB_CACHE_DIR is empty")
    cache = Path(raw)
    if cache.is_dir():
        n = sum(1 for _ in cache.glob("*"))
        return ("ok", f"{cache} ({n} cached entries)")
    return ("missing", f"{cache} does not exist yet (created on first use)")


@app.callback(invoke_without_command=True)
def doctor(
    kraken_db: Path | None = typer.Option(
        None,
        "--kraken-db",
        help="Optionally validate a Kraken2 database directory.",
    ),
    strict: bool = typer.Option(
        False,
        "--strict",
        help=(
            "Exit non-zero if a required tool or package is missing "
            "(blastn, makeblastdb, kraken2, skani or fastANI; core Python "
            "packages) or a given --kraken-db is invalid."
        ),
    ),
) -> None:
    """Print a diagnostic report of the metadarkmatter environment."""
    console.print(f"[bold]metadarkmatter[/bold] {__version__}")
    console.print(f"Python {sys.version.split()[0]} ({sys.executable})")
    console.print(f"Platform: {platform.platform()}")
    console.print()

    env_table = Table(title="Environment variables", show_header=True, header_style="bold")
    env_table.add_column("Variable")
    env_table.add_column("Value")
    env_table.add_column("Purpose", overflow="fold")
    for var, purpose in _ENV_VARS:
        raw = os.environ.get(var)
        display = raw if raw is not None else "[dim]unset (default)[/dim]"
        env_table.add_row(var, display, purpose)
    console.print(env_table)
    console.print()

    pkg_table = Table(title="Python dependencies", show_header=True, header_style="bold")
    pkg_table.add_column("Package")
    pkg_table.add_column("Version")
    missing_packages: list[str] = []
    for pkg in _PYTHON_PACKAGES:
        version = _probe_package(pkg)
        if version == "not installed" and pkg in _REQUIRED_PACKAGES:
            missing_packages.append(pkg)
        pkg_table.add_row(pkg, version)
    console.print(pkg_table)
    console.print()

    tool_table = Table(title="External tools", show_header=True, header_style="bold")
    tool_table.add_column("Tool")
    tool_table.add_column("Path")
    tool_table.add_column("Version")
    tool_table.add_column("Needed by", overflow="fold")
    missing: list[str] = []
    present: set[str] = set()
    for exe, args in _EXTERNAL_TOOLS:
        path, version = _probe_tool(exe, args)
        if not path:
            missing.append(exe)
            tool_table.add_row(
                exe, "-", "[yellow]not found[/yellow]", _TOOL_NEEDED_BY.get(exe, "")
            )
        else:
            present.add(exe)
            tool_table.add_row(exe, path, version, _TOOL_NEEDED_BY.get(exe, ""))
    console.print(tool_table)

    if missing:
        console.print(
            "\n[yellow]"
            f"{len(missing)} tool(s) not on PATH: {', '.join(missing)}. "
            "Install only the tools your workflow needs:[/yellow]"
        )
        # De-duplicate install hints (blast covers blastn/blastx/makeblastdb).
        seen_hints: set[str] = set()
        for exe in missing:
            hint = _TOOL_INSTALL_HINTS.get(exe)
            if hint and hint not in seen_hints:
                seen_hints.add(hint)
                console.print(f"  [dim]{hint}[/dim]")

    # Database / cache checks.
    console.print()
    db_table = Table(title="Databases & caches", show_header=True, header_style="bold")
    db_table.add_column("Resource")
    db_table.add_column("State")
    db_table.add_column("Detail", overflow="fold")

    gtdb_state, gtdb_detail = _check_gtdb_cache()
    db_table.add_row("GTDB cache", gtdb_state, gtdb_detail)

    if kraken_db is not None:
        ok, detail = _check_kraken_db(kraken_db)
        db_table.add_row("Kraken2 DB", "ok" if ok else "[yellow]problem[/yellow]", detail)
    else:
        db_table.add_row(
            "Kraken2 DB", "[dim]not checked[/dim]", "pass --kraken-db PATH to validate"
        )
    console.print(db_table)

    # Strict mode: exit non-zero when required components are missing.
    problems: list[str] = []
    for group in _REQUIRED_TOOL_GROUPS:
        if not any(t in present for t in group):
            problems.append("missing required tool: " + " or ".join(group))
    problems.extend(f"missing required package: {p}" for p in missing_packages)
    if kraken_db is not None and not _check_kraken_db(kraken_db)[0]:
        problems.append(f"invalid Kraken2 DB: {kraken_db}")
    if problems:
        console.print()
        for problem in problems:
            console.print(f"[red]{problem}[/red]")
        if strict:
            raise typer.Exit(code=1)
        console.print("[dim]Re-run with --strict to exit non-zero on these.[/dim]")
