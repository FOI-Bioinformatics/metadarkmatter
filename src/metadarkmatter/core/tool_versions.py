"""Version probes for external tools (shared by ``mdm doctor`` and pipelines)."""

from __future__ import annotations

import shutil
import subprocess

# Arguments that make each tool print its version.
VERSION_ARGS: dict[str, list[str]] = {
    "blastn": ["-version"],
    "blastx": ["-version"],
    "makeblastdb": ["-version"],
    "mmseqs": ["version"],
    "kraken2": ["--version"],
    "fastANI": ["--version"],
    "skani": ["--version"],
    "diamond": ["version"],
    "mashtree": ["--version"],
    "mash": ["--version"],
    "bowtie2": ["--version"],
    "bowtie2-build": ["--version"],
    "samtools": ["--version"],
    "prodigal": ["-v"],
}


def probe_tool(executable: str, version_args: list[str] | None = None) -> tuple[str, str]:
    """Return (path, version_string) for an external tool.

    Returns ("", "not found") when the tool is not on PATH. The version string
    is the first non-empty line of combined stdout/stderr, truncated to keep
    reports scannable.
    """
    if version_args is None:
        version_args = VERSION_ARGS.get(executable, ["--version"])
    path = shutil.which(executable)
    if path is None:
        return ("", "not found")
    try:
        result = subprocess.run(
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
