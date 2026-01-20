"""
Base classes for wrapping external bioinformatics tools.

Provides a consistent interface for executing command-line tools
with proper error handling, timeout support, and dry-run capability.
"""

from __future__ import annotations

import logging
import re
import shutil
import subprocess
import time
from abc import ABC, abstractmethod
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar

from metadarkmatter.core.exceptions import MetadarkmatterError

logger = logging.getLogger(__name__)

# Pattern for safe path characters (alphanumeric, underscore, hyphen, dot, slash)
# This is conservative but covers typical bioinformatics file paths
_SAFE_PATH_PATTERN = re.compile(r"^[\w\-./]+$")


class UnsafePathError(MetadarkmatterError):
    """Raised when a file path contains potentially unsafe characters."""

    def __init__(self, path: Path, reason: str = ""):
        detail = f": {reason}" if reason else ""
        super().__init__(
            message=f"Unsafe path detected: {path}{detail}",
            suggestion=(
                "Ensure file paths contain only alphanumeric characters, "
                "underscores, hyphens, and periods. Avoid spaces and special characters."
            ),
        )
        self.path = path


def validate_path_safe(
    path: Path,
    *,
    must_exist: bool = False,
    resolve: bool = True,
) -> Path:
    """Validate that a path is safe for use in subprocess commands.

    This provides defense-in-depth for subprocess calls. While list-based
    subprocess execution prevents shell injection, special characters in
    paths could still cause issues with underlying bioinformatics tools.

    Args:
        path: Path to validate
        must_exist: If True, raise error if path doesn't exist
        resolve: If True, resolve the path to its absolute form

    Returns:
        The validated (and optionally resolved) path

    Raises:
        UnsafePathError: If path contains suspicious characters
        FileNotFoundError: If must_exist=True and path doesn't exist
    """
    if resolve:
        path = path.resolve()

    path_str = str(path)

    # Check for null bytes (could truncate path)
    if "\x00" in path_str:
        raise UnsafePathError(path, "contains null byte")

    # Check for path traversal patterns
    if ".." in path.parts:
        # Allow .. only if it resolves to a path that doesn't escape
        # We already resolved above, so this catches literal .. in resolved path
        # Log at warning level for security auditing
        logger.warning("Path contains path traversal pattern (..): %s", path)

    # Warn about unusual characters (but don't block - just log)
    if not _SAFE_PATH_PATTERN.match(path_str):
        logger.warning(
            "Path contains unusual characters (may cause issues): %s",
            path,
        )

    if must_exist and not path.exists():
        raise FileNotFoundError(f"Path does not exist: {path}")

    return path


def validate_paths_safe(
    *paths: Path,
    must_exist: bool = False,
    resolve: bool = True,
) -> tuple[Path, ...]:
    """Validate multiple paths for safety.

    Args:
        *paths: Paths to validate
        must_exist: If True, raise error if any path doesn't exist
        resolve: If True, resolve paths to absolute form

    Returns:
        Tuple of validated paths

    Raises:
        UnsafePathError: If any path contains suspicious characters
        FileNotFoundError: If must_exist=True and any path doesn't exist
    """
    return tuple(
        validate_path_safe(p, must_exist=must_exist, resolve=resolve)
        for p in paths
    )


class ToolNotFoundError(MetadarkmatterError):
    """Raised when a required external tool is not installed or not in PATH."""

    def __init__(self, tool_name: str, install_hint: str = ""):
        suggestion = f"Install {tool_name} and ensure it is in your PATH."
        if install_hint:
            suggestion = f"{suggestion}\n\nInstallation:\n  {install_hint}"

        super().__init__(
            message=f"Required tool '{tool_name}' not found in PATH",
            suggestion=suggestion,
        )
        self.tool_name = tool_name


class ToolExecutionError(MetadarkmatterError):
    """Raised when an external tool returns a non-zero exit code."""

    def __init__(
        self,
        tool_name: str,
        command: list[str],
        return_code: int,
        stderr: str,
    ):
        # Truncate long commands and stderr for display
        cmd_str = " ".join(command)
        if len(cmd_str) > 200:
            cmd_str = cmd_str[:200] + "..."

        stderr_display = stderr.strip()
        if len(stderr_display) > 500:
            stderr_display = stderr_display[:500] + "\n...[truncated]"

        super().__init__(
            message=(
                f"{tool_name} failed with exit code {return_code}\n\n"
                f"Command: {cmd_str}\n\n"
                f"Error output:\n{stderr_display}"
            ),
            suggestion=(
                "Check the command parameters and input files. "
                "Run with --verbose for detailed output."
            ),
        )
        self.tool_name = tool_name
        self.command = command
        self.return_code = return_code
        self.stderr = stderr


class ToolTimeoutError(MetadarkmatterError):
    """Raised when an external tool exceeds the specified timeout."""

    def __init__(self, tool_name: str, timeout_seconds: float, command: list[str]):
        cmd_str = " ".join(command)
        if len(cmd_str) > 200:
            cmd_str = cmd_str[:200] + "..."

        super().__init__(
            message=(
                f"{tool_name} timed out after {timeout_seconds:.0f} seconds\n\n"
                f"Command: {cmd_str}"
            ),
            suggestion=(
                "Increase the timeout or check if the tool is hanging. "
                "Consider processing smaller input files."
            ),
        )
        self.tool_name = tool_name
        self.timeout_seconds = timeout_seconds
        self.command = command


@dataclass(frozen=True)
class ToolResult:
    """Result from running an external tool.

    Attributes:
        command: The command that was executed.
        return_code: Exit code from the process.
        stdout: Standard output from the process.
        stderr: Standard error from the process.
        elapsed_seconds: Wall-clock time for execution.
    """

    command: tuple[str, ...]
    return_code: int
    stdout: str
    stderr: str
    elapsed_seconds: float

    @property
    def success(self) -> bool:
        """Return True if the tool exited with code 0."""
        return self.return_code == 0

    @property
    def command_string(self) -> str:
        """Return the command as a space-separated string."""
        return " ".join(self.command)

    def format_output(self, max_lines: int = 50) -> str:
        """Format stdout/stderr for display.

        Args:
            max_lines: Maximum number of lines to include.

        Returns:
            Formatted string with stdout and stderr.
        """
        parts = []

        if self.stdout.strip():
            stdout_lines = self.stdout.strip().split("\n")
            if len(stdout_lines) > max_lines:
                truncated = len(stdout_lines) - max_lines
                stdout_lines = [*stdout_lines[:max_lines], f"... ({truncated} more lines)"]
            parts.append("stdout:\n" + "\n".join(stdout_lines))

        if self.stderr.strip():
            stderr_lines = self.stderr.strip().split("\n")
            if len(stderr_lines) > max_lines:
                truncated = len(stderr_lines) - max_lines
                stderr_lines = [*stderr_lines[:max_lines], f"... ({truncated} more lines)"]
            parts.append("stderr:\n" + "\n".join(stderr_lines))

        return "\n\n".join(parts) if parts else "(no output)"


class ExternalTool(ABC):
    """Abstract base class for wrapping external command-line tools.

    Subclasses must define:
        TOOL_NAME: Primary executable name (e.g., "kraken2")
        build_command: Method to construct the command arguments

    Optional class attributes:
        TOOL_ALIASES: Alternative executable names to search
        INSTALL_HINT: Instructions for installing the tool

    Dependency injection:
        Use set_executable_resolver() to inject a custom resolver for testing.
        This allows tests to mock executable lookups without requiring tools.
    """

    TOOL_NAME: ClassVar[str]
    TOOL_ALIASES: ClassVar[tuple[str, ...]] = ()
    INSTALL_HINT: ClassVar[str] = ""

    _executable_cache: ClassVar[dict[str, Path | None]] = {}
    # Dependency injection for testability - defaults to shutil.which
    _executable_resolver: ClassVar[Callable[[str], str | None]] = staticmethod(
        shutil.which
    )

    @classmethod
    def check_available(cls) -> bool:
        """Check if the tool is installed and available in PATH.

        Returns:
            True if the tool can be found, False otherwise.
        """
        try:
            cls.get_executable()
            return True
        except ToolNotFoundError:
            return False

    @classmethod
    def get_executable(cls) -> Path:
        """Find the tool executable in PATH.

        Uses the configured resolver (default: shutil.which) to locate
        executables. Override with set_executable_resolver() for testing.

        Returns:
            Path to the executable.

        Raises:
            ToolNotFoundError: If the tool cannot be found.
        """
        # Check cache first
        if cls.TOOL_NAME in cls._executable_cache:
            cached = cls._executable_cache[cls.TOOL_NAME]
            if cached is not None:
                return cached
            raise ToolNotFoundError(cls.TOOL_NAME, cls.INSTALL_HINT)

        # Search for primary name and aliases using injected resolver
        names_to_try = (cls.TOOL_NAME, *cls.TOOL_ALIASES)

        for name in names_to_try:
            exe_path = cls._executable_resolver(name)
            if exe_path:
                path = Path(exe_path)
                cls._executable_cache[cls.TOOL_NAME] = path
                return path

        # Not found - cache the failure
        cls._executable_cache[cls.TOOL_NAME] = None
        raise ToolNotFoundError(cls.TOOL_NAME, cls.INSTALL_HINT)

    @classmethod
    def clear_cache(cls) -> None:
        """Clear the executable location cache."""
        cls._executable_cache.clear()

    @classmethod
    def set_executable_resolver(
        cls,
        resolver: Callable[[str], str | None],
    ) -> None:
        """Inject a custom executable resolver for testing.

        This allows tests to mock executable lookups without requiring
        the actual bioinformatics tools to be installed.

        Args:
            resolver: Function that takes a tool name and returns
                the path to the executable or None if not found.

        Example:
            # Mock resolver that pretends all tools exist
            def mock_resolver(name):
                return f"/usr/bin/{name}"

            ExternalTool.set_executable_resolver(mock_resolver)
            ExternalTool.clear_cache()  # Clear cached lookups

            # Run tests...

            ExternalTool.reset_executable_resolver()
        """
        cls._executable_resolver = staticmethod(resolver)
        cls.clear_cache()  # Clear cached lookups when resolver changes

    @classmethod
    def reset_executable_resolver(cls) -> None:
        """Reset the executable resolver to the default (shutil.which)."""
        cls._executable_resolver = staticmethod(shutil.which)
        cls.clear_cache()

    @abstractmethod
    def build_command(self, **kwargs: object) -> list[str]:
        """Build the command-line arguments for this tool.

        Returns:
            List of command-line arguments (including the executable).
        """
        ...

    def run(
        self,
        *,
        timeout: float | None = None,
        dry_run: bool = False,
        capture_output: bool = True,
        **kwargs: object,
    ) -> ToolResult:
        """Execute the tool with the specified arguments.

        Args:
            timeout: Maximum execution time in seconds (None for no limit).
            dry_run: If True, return command without execution.
            capture_output: If True, capture stdout/stderr.
            **kwargs: Arguments passed to build_command().

        Returns:
            ToolResult with command, exit code, and output.

        Raises:
            ToolNotFoundError: If the tool is not installed.
            ToolTimeoutError: If execution exceeds timeout.
        """
        command = self.build_command(**kwargs)
        command_tuple = tuple(command)

        if dry_run:
            return ToolResult(
                command=command_tuple,
                return_code=0,
                stdout="[dry-run] Command not executed",
                stderr="",
                elapsed_seconds=0.0,
            )

        start_time = time.perf_counter()

        try:
            if capture_output:
                result = subprocess.run(
                    command,
                    capture_output=True,
                    text=True,
                    timeout=timeout,
                )
                stdout = result.stdout
                stderr = result.stderr
            else:
                # Redirect to DEVNULL to prevent pipe buffer blocking
                # This is critical for tools like Diamond that write progress
                # to stderr - without DEVNULL, the buffer fills and blocks
                result = subprocess.run(
                    command,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    timeout=timeout,
                )
                stdout = ""
                stderr = ""

            elapsed = time.perf_counter() - start_time

            return ToolResult(
                command=command_tuple,
                return_code=result.returncode,
                stdout=stdout,
                stderr=stderr,
                elapsed_seconds=elapsed,
            )

        except subprocess.TimeoutExpired as e:
            raise ToolTimeoutError(
                self.TOOL_NAME,
                timeout or 0,
                command,
            ) from e

        except FileNotFoundError as e:
            # This can happen if the executable disappears between check and run
            raise ToolNotFoundError(self.TOOL_NAME, self.INSTALL_HINT) from e

    def run_or_raise(
        self,
        *,
        timeout: float | None = None,
        dry_run: bool = False,
        capture_output: bool = True,
        **kwargs: object,
    ) -> ToolResult:
        """Execute the tool and raise an exception on failure.

        Same as run() but raises ToolExecutionError if exit code is non-zero.

        Returns:
            ToolResult with command, exit code, and output.

        Raises:
            ToolNotFoundError: If the tool is not installed.
            ToolTimeoutError: If execution exceeds timeout.
            ToolExecutionError: If the tool returns non-zero exit code.
        """
        result = self.run(
            timeout=timeout,
            dry_run=dry_run,
            capture_output=capture_output,
            **kwargs,
        )

        if not result.success and not dry_run:
            raise ToolExecutionError(
                self.TOOL_NAME,
                list(result.command),
                result.return_code,
                result.stderr,
            )

        return result

    def get_version(self, *, timeout: float = 30.0) -> str | None:
        """Attempt to get the tool version.

        Returns:
            Version string if available, None otherwise.
        """
        # Common version flags to try
        version_flags = ["--version", "-version", "-v", "version"]

        for flag in version_flags:
            try:
                exe = str(self.get_executable())
                result = subprocess.run(
                    [exe, flag],
                    capture_output=True,
                    text=True,
                    timeout=timeout,
                )
                if result.returncode == 0 and result.stdout.strip():
                    # Return first non-empty line
                    for line in result.stdout.strip().split("\n"):
                        if line.strip():
                            return line.strip()
            except (subprocess.TimeoutExpired, FileNotFoundError, ToolNotFoundError):
                continue

        return None
