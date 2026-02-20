"""Unit tests for external tool base classes and path validation.

Covers UnsafePathError, validate_path_safe, validate_paths_safe,
ToolResult.format_output, ExternalTool lifecycle (check_available,
get_executable, run, run_or_raise, get_version), and dependency
injection via set_executable_resolver / reset_executable_resolver.

All tests use mocking so that no external bioinformatics tools are
required on the host system.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import ClassVar
from unittest.mock import MagicMock, patch

import pytest

from metadarkmatter.external.base import (
    ExternalTool,
    ToolExecutionError,
    ToolNotFoundError,
    ToolResult,
    ToolTimeoutError,
    UnsafePathError,
    validate_path_safe,
    validate_paths_safe,
)


# ---------------------------------------------------------------------------
# Concrete subclass of ExternalTool used throughout the test module
# ---------------------------------------------------------------------------


class _MockTool(ExternalTool):
    """Minimal concrete subclass of ExternalTool for testing."""

    TOOL_NAME: ClassVar[str] = "mocktool"
    TOOL_ALIASES: ClassVar[tuple[str, ...]] = ("mocktool-alt",)
    INSTALL_HINT: ClassVar[str] = "pip install mocktool"

    def build_command(self, **kwargs: object) -> list[str]:
        exe = str(self.get_executable())
        cmd = [exe, "run"]
        if "extra" in kwargs:
            cmd.append(str(kwargs["extra"]))
        return cmd


@pytest.fixture(autouse=True)
def _clear_caches():
    """Ensure the executable cache is empty before and after each test."""
    ExternalTool._executable_cache.clear()
    yield
    ExternalTool._executable_cache.clear()


# =========================================================================
# UnsafePathError
# =========================================================================


class TestUnsafePathError:
    """Tests for UnsafePathError exception formatting."""

    def test_without_reason(self):
        """Should include path in message when no reason is given."""
        err = UnsafePathError(Path("/bad/path"))
        assert "/bad/path" in str(err)
        assert "Unsafe path" in str(err)

    def test_with_reason(self):
        """Should include both path and reason in message."""
        err = UnsafePathError(Path("/bad/path"), reason="contains null byte")
        assert "/bad/path" in str(err)
        assert "contains null byte" in str(err)

    def test_stores_path_attribute(self):
        """Should store the offending path as an attribute."""
        p = Path("/some/path")
        err = UnsafePathError(p)
        assert err.path == p

    def test_includes_suggestion(self):
        """Should include a suggestion about safe characters."""
        err = UnsafePathError(Path("/bad"))
        assert "alphanumeric" in str(err).lower() or "Suggestion" in str(err)


# =========================================================================
# validate_path_safe
# =========================================================================


class TestValidatePathSafe:
    """Tests for the validate_path_safe utility function."""

    def test_simple_safe_path(self, tmp_path: Path):
        """Should accept a simple safe path and resolve it."""
        p = tmp_path / "data"
        result = validate_path_safe(p)
        assert result.is_absolute()

    def test_resolve_disabled(self, tmp_path: Path):
        """Should return path as-is when resolve=False."""
        p = Path("relative/path")
        result = validate_path_safe(p, resolve=False)
        assert result == p

    def test_null_byte_raises(self):
        """Should raise UnsafePathError for paths with null bytes."""
        p = Path("/some/path\x00evil")
        with pytest.raises(UnsafePathError, match="null byte"):
            validate_path_safe(p, resolve=False)

    def test_path_traversal_warning(self, caplog: pytest.LogCaptureFixture):
        """Should log a warning when the resolved path still contains '..'."""
        # Create a Path whose .parts include '..'.  We disable resolve so
        # that the literal '..' is preserved in parts.
        p = Path("/safe/../traversal")
        with caplog.at_level("WARNING"):
            validate_path_safe(p, resolve=False)
        assert "path traversal" in caplog.text.lower() or ".." in caplog.text

    def test_unusual_characters_warning(self, caplog: pytest.LogCaptureFixture):
        """Should log a warning for paths with spaces or special characters."""
        p = Path("/path with spaces/file.txt")
        with caplog.at_level("WARNING"):
            validate_path_safe(p, resolve=False)
        assert "unusual characters" in caplog.text.lower()

    def test_must_exist_raises_when_missing(self, tmp_path: Path):
        """Should raise FileNotFoundError when must_exist=True and path is absent."""
        p = tmp_path / "does_not_exist.txt"
        with pytest.raises(FileNotFoundError, match="does not exist"):
            validate_path_safe(p, must_exist=True)

    def test_must_exist_succeeds_when_present(self, tmp_path: Path):
        """Should succeed when must_exist=True and path does exist."""
        p = tmp_path / "exists.txt"
        p.write_text("content")
        result = validate_path_safe(p, must_exist=True)
        assert result.exists()


# =========================================================================
# validate_paths_safe
# =========================================================================


class TestValidatePathsSafe:
    """Tests for the validate_paths_safe batch helper."""

    def test_validates_multiple_paths(self, tmp_path: Path):
        """Should return a tuple of validated paths."""
        p1 = tmp_path / "a.txt"
        p2 = tmp_path / "b.txt"
        results = validate_paths_safe(p1, p2)
        assert isinstance(results, tuple)
        assert len(results) == 2
        for r in results:
            assert r.is_absolute()

    def test_raises_on_unsafe_path(self):
        """Should raise UnsafePathError if any path is unsafe."""
        good = Path("/safe/path")
        bad = Path("/has\x00null")
        with pytest.raises(UnsafePathError):
            validate_paths_safe(good, bad, resolve=False)

    def test_must_exist_propagated(self, tmp_path: Path):
        """Should propagate must_exist to each individual validation."""
        p = tmp_path / "missing.txt"
        with pytest.raises(FileNotFoundError):
            validate_paths_safe(p, must_exist=True)


# =========================================================================
# ToolResult.format_output
# =========================================================================


class TestToolResultFormatOutput:
    """Tests for ToolResult.format_output method."""

    def test_no_output(self):
        """Should return '(no output)' when both stdout and stderr are empty."""
        result = ToolResult(
            command=("cmd",),
            return_code=0,
            stdout="",
            stderr="",
            elapsed_seconds=0.0,
        )
        assert result.format_output() == "(no output)"

    def test_stdout_only(self):
        """Should format stdout when only stdout is present."""
        result = ToolResult(
            command=("cmd",),
            return_code=0,
            stdout="line1\nline2\n",
            stderr="",
            elapsed_seconds=0.0,
        )
        formatted = result.format_output()
        assert "stdout:" in formatted
        assert "line1" in formatted
        assert "line2" in formatted

    def test_stderr_only(self):
        """Should format stderr when only stderr is present."""
        result = ToolResult(
            command=("cmd",),
            return_code=1,
            stdout="",
            stderr="error info\n",
            elapsed_seconds=0.0,
        )
        formatted = result.format_output()
        assert "stderr:" in formatted
        assert "error info" in formatted

    def test_both_streams(self):
        """Should include both stdout and stderr when present."""
        result = ToolResult(
            command=("cmd",),
            return_code=0,
            stdout="output\n",
            stderr="warning\n",
            elapsed_seconds=0.0,
        )
        formatted = result.format_output()
        assert "stdout:" in formatted
        assert "stderr:" in formatted

    def test_truncation_stdout(self):
        """Should truncate stdout exceeding max_lines."""
        many_lines = "\n".join(f"line_{i}" for i in range(100))
        result = ToolResult(
            command=("cmd",),
            return_code=0,
            stdout=many_lines,
            stderr="",
            elapsed_seconds=0.0,
        )
        formatted = result.format_output(max_lines=10)
        assert "more lines" in formatted

    def test_truncation_stderr(self):
        """Should truncate stderr exceeding max_lines."""
        many_lines = "\n".join(f"err_{i}" for i in range(100))
        result = ToolResult(
            command=("cmd",),
            return_code=1,
            stdout="",
            stderr=many_lines,
            elapsed_seconds=0.0,
        )
        formatted = result.format_output(max_lines=5)
        assert "more lines" in formatted


# =========================================================================
# ToolExecutionError - long command truncation
# =========================================================================


class TestToolExecutionErrorTruncation:
    """Tests for long command / stderr truncation in ToolExecutionError."""

    def test_long_command_is_truncated(self):
        """Should truncate the command string when it exceeds 200 characters."""
        long_cmd = ["tool"] + [f"arg_{i:05d}" for i in range(100)]
        err = ToolExecutionError(
            tool_name="tool",
            command=long_cmd,
            return_code=1,
            stderr="oops",
        )
        assert "..." in str(err)


# =========================================================================
# ToolTimeoutError - long command truncation
# =========================================================================


class TestToolTimeoutErrorTruncation:
    """Tests for long command truncation in ToolTimeoutError."""

    def test_long_command_is_truncated(self):
        """Should truncate the command string when it exceeds 200 characters."""
        long_cmd = ["tool"] + [f"arg_{i:05d}" for i in range(100)]
        err = ToolTimeoutError(
            tool_name="tool",
            timeout_seconds=120.0,
            command=long_cmd,
        )
        assert "..." in str(err)


# =========================================================================
# ExternalTool.check_available
# =========================================================================


class TestCheckAvailable:
    """Tests for ExternalTool.check_available class method."""

    def test_returns_true_when_found(self):
        """Should return True when the executable is found."""
        _MockTool._executable_cache["mocktool"] = Path("/usr/bin/mocktool")
        assert _MockTool.check_available() is True

    def test_returns_false_when_not_found(self):
        """Should return False when no executable is found."""
        _MockTool._executable_resolver = staticmethod(lambda name: None)
        try:
            assert _MockTool.check_available() is False
        finally:
            _MockTool.reset_executable_resolver()


# =========================================================================
# ExternalTool.get_executable
# =========================================================================


class TestGetExecutable:
    """Tests for ExternalTool.get_executable class method."""

    def test_returns_cached_path(self):
        """Should return the cached path without calling the resolver."""
        _MockTool._executable_cache["mocktool"] = Path("/cached/mocktool")
        assert _MockTool.get_executable() == Path("/cached/mocktool")

    def test_raises_when_cached_none(self):
        """Should raise ToolNotFoundError when cached value is None."""
        _MockTool._executable_cache["mocktool"] = None
        with pytest.raises(ToolNotFoundError, match="mocktool"):
            _MockTool.get_executable()

    def test_finds_primary_name(self):
        """Should find the tool by its primary name via the resolver."""
        _MockTool._executable_resolver = staticmethod(
            lambda name: "/usr/bin/mocktool" if name == "mocktool" else None
        )
        try:
            result = _MockTool.get_executable()
            assert result == Path("/usr/bin/mocktool")
            # Verify it was cached
            assert _MockTool._executable_cache["mocktool"] == Path("/usr/bin/mocktool")
        finally:
            _MockTool.reset_executable_resolver()

    def test_finds_alias(self):
        """Should fall back to alias if primary name is not found."""
        _MockTool._executable_resolver = staticmethod(
            lambda name: "/usr/bin/mocktool-alt" if name == "mocktool-alt" else None
        )
        try:
            result = _MockTool.get_executable()
            assert result == Path("/usr/bin/mocktool-alt")
        finally:
            _MockTool.reset_executable_resolver()

    def test_caches_failure(self):
        """Should cache None and raise on subsequent calls."""
        _MockTool._executable_resolver = staticmethod(lambda name: None)
        try:
            with pytest.raises(ToolNotFoundError):
                _MockTool.get_executable()
            # Verify failure is cached
            assert _MockTool._executable_cache["mocktool"] is None
        finally:
            _MockTool.reset_executable_resolver()

    def test_install_hint_in_error(self):
        """Should include install hint in ToolNotFoundError."""
        _MockTool._executable_resolver = staticmethod(lambda name: None)
        try:
            with pytest.raises(ToolNotFoundError) as exc_info:
                _MockTool.get_executable()
            assert "pip install mocktool" in str(exc_info.value)
        finally:
            _MockTool.reset_executable_resolver()


# =========================================================================
# ExternalTool.clear_cache
# =========================================================================


class TestClearCache:
    """Tests for ExternalTool.clear_cache class method."""

    def test_clears_all_entries(self):
        """Should remove all entries from the executable cache."""
        _MockTool._executable_cache["mocktool"] = Path("/usr/bin/mocktool")
        _MockTool._executable_cache["other"] = Path("/usr/bin/other")
        _MockTool.clear_cache()
        assert len(_MockTool._executable_cache) == 0


# =========================================================================
# set_executable_resolver / reset_executable_resolver
# =========================================================================


class TestExecutableResolverInjection:
    """Tests for dependency injection of executable resolvers."""

    def test_set_resolver_clears_cache(self):
        """Should clear the cache when a new resolver is set."""
        _MockTool._executable_cache["mocktool"] = Path("/old/path")
        _MockTool.set_executable_resolver(lambda name: "/new/path")
        assert len(_MockTool._executable_cache) == 0

    def test_reset_restores_default(self):
        """Should restore shutil.which as the default resolver."""
        _MockTool.set_executable_resolver(lambda name: "/custom")
        _MockTool.reset_executable_resolver()
        # After reset, calling get_executable uses shutil.which.
        # Since 'mocktool' likely is not installed, this should raise.
        with pytest.raises(ToolNotFoundError):
            _MockTool.get_executable()

    def test_custom_resolver_is_used(self):
        """Should use the injected resolver for executable lookups."""
        _MockTool.set_executable_resolver(lambda name: f"/custom/{name}")
        try:
            path = _MockTool.get_executable()
            assert path == Path("/custom/mocktool")
        finally:
            _MockTool.reset_executable_resolver()


# =========================================================================
# ExternalTool.run
# =========================================================================


class TestExternalToolRun:
    """Tests for ExternalTool.run with mocked subprocess calls."""

    def _setup_tool(self) -> _MockTool:
        """Create a _MockTool with a cached executable path."""
        _MockTool._executable_cache["mocktool"] = Path("/usr/bin/mocktool")
        return _MockTool()

    def test_dry_run_returns_without_executing(self):
        """Should return a dry-run result without calling subprocess."""
        tool = self._setup_tool()
        result = tool.run(dry_run=True)
        assert result.success is True
        assert "[dry-run]" in result.stdout
        assert result.elapsed_seconds == 0.0
        assert "mocktool" in result.command_string

    @patch("metadarkmatter.external.base.subprocess.run")
    def test_capture_output_true(self, mock_run: MagicMock):
        """Should capture stdout and stderr when capture_output=True."""
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="captured output",
            stderr="captured warning",
        )
        tool = self._setup_tool()
        result = tool.run(capture_output=True)

        assert result.success is True
        assert result.stdout == "captured output"
        assert result.stderr == "captured warning"
        assert result.elapsed_seconds >= 0.0
        mock_run.assert_called_once()
        call_kwargs = mock_run.call_args
        assert call_kwargs.kwargs["capture_output"] is True
        assert call_kwargs.kwargs["text"] is True

    @patch("metadarkmatter.external.base.subprocess.run")
    def test_capture_output_false(self, mock_run: MagicMock):
        """Should redirect to DEVNULL when capture_output=False."""
        mock_run.return_value = MagicMock(returncode=0)
        tool = self._setup_tool()
        result = tool.run(capture_output=False)

        assert result.success is True
        assert result.stdout == ""
        assert result.stderr == ""
        call_kwargs = mock_run.call_args
        assert call_kwargs.kwargs["stdout"] == subprocess.DEVNULL
        assert call_kwargs.kwargs["stderr"] == subprocess.DEVNULL

    @patch("metadarkmatter.external.base.subprocess.run")
    def test_nonzero_exit_code(self, mock_run: MagicMock):
        """Should return failure result for non-zero exit code."""
        mock_run.return_value = MagicMock(
            returncode=1,
            stdout="",
            stderr="error occurred",
        )
        tool = self._setup_tool()
        result = tool.run()

        assert result.success is False
        assert result.return_code == 1
        assert result.stderr == "error occurred"

    @patch("metadarkmatter.external.base.subprocess.run")
    def test_timeout_raises_tool_timeout_error(self, mock_run: MagicMock):
        """Should raise ToolTimeoutError when subprocess times out."""
        mock_run.side_effect = subprocess.TimeoutExpired(
            cmd=["mocktool", "run"], timeout=30
        )
        tool = self._setup_tool()
        with pytest.raises(ToolTimeoutError, match="mocktool"):
            tool.run(timeout=30.0)

    @patch("metadarkmatter.external.base.subprocess.run")
    def test_file_not_found_raises_tool_not_found(self, mock_run: MagicMock):
        """Should raise ToolNotFoundError when executable disappears at runtime."""
        mock_run.side_effect = FileNotFoundError("No such file")
        tool = self._setup_tool()
        with pytest.raises(ToolNotFoundError, match="mocktool"):
            tool.run()

    @patch("metadarkmatter.external.base.subprocess.run")
    def test_timeout_parameter_forwarded(self, mock_run: MagicMock):
        """Should forward the timeout parameter to subprocess.run."""
        mock_run.return_value = MagicMock(
            returncode=0, stdout="", stderr=""
        )
        tool = self._setup_tool()
        tool.run(timeout=60.0)

        call_kwargs = mock_run.call_args
        assert call_kwargs.kwargs["timeout"] == 60.0

    @patch("metadarkmatter.external.base.subprocess.run")
    def test_extra_kwargs_passed_to_build_command(self, mock_run: MagicMock):
        """Should forward extra kwargs through to build_command."""
        mock_run.return_value = MagicMock(
            returncode=0, stdout="", stderr=""
        )
        tool = self._setup_tool()
        result = tool.run(extra="--verbose")
        assert "--verbose" in result.command


# =========================================================================
# ExternalTool.run_or_raise
# =========================================================================


class TestExternalToolRunOrRaise:
    """Tests for ExternalTool.run_or_raise."""

    def _setup_tool(self) -> _MockTool:
        _MockTool._executable_cache["mocktool"] = Path("/usr/bin/mocktool")
        return _MockTool()

    @patch("metadarkmatter.external.base.subprocess.run")
    def test_success_returns_result(self, mock_run: MagicMock):
        """Should return the result on success."""
        mock_run.return_value = MagicMock(
            returncode=0, stdout="ok", stderr=""
        )
        tool = self._setup_tool()
        result = tool.run_or_raise()
        assert result.success is True

    @patch("metadarkmatter.external.base.subprocess.run")
    def test_failure_raises_execution_error(self, mock_run: MagicMock):
        """Should raise ToolExecutionError on non-zero exit code."""
        mock_run.return_value = MagicMock(
            returncode=2, stdout="", stderr="failed"
        )
        tool = self._setup_tool()
        with pytest.raises(ToolExecutionError) as exc_info:
            tool.run_or_raise()
        assert exc_info.value.return_code == 2
        assert exc_info.value.tool_name == "mocktool"

    def test_dry_run_does_not_raise(self):
        """Should not raise even though exit code check is present in dry_run."""
        tool = self._setup_tool()
        result = tool.run_or_raise(dry_run=True)
        assert result.success is True

    @patch("metadarkmatter.external.base.subprocess.run")
    def test_timeout_propagated(self, mock_run: MagicMock):
        """Should propagate timeout to the underlying run call."""
        mock_run.side_effect = subprocess.TimeoutExpired(
            cmd=["mocktool"], timeout=10
        )
        tool = self._setup_tool()
        with pytest.raises(ToolTimeoutError):
            tool.run_or_raise(timeout=10.0)


# =========================================================================
# ExternalTool.get_version
# =========================================================================


class TestGetVersion:
    """Tests for ExternalTool.get_version."""

    def _setup_tool(self) -> _MockTool:
        _MockTool._executable_cache["mocktool"] = Path("/usr/bin/mocktool")
        return _MockTool()

    @patch("metadarkmatter.external.base.subprocess.run")
    def test_returns_first_nonempty_line(self, mock_run: MagicMock):
        """Should return the first non-empty line of --version output."""
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="mocktool v1.2.3\nCopyright 2024\n",
        )
        tool = self._setup_tool()
        version = tool.get_version()
        assert version == "mocktool v1.2.3"

    @patch("metadarkmatter.external.base.subprocess.run")
    def test_tries_multiple_flags(self, mock_run: MagicMock):
        """Should try multiple version flags until one succeeds."""
        def side_effect(cmd, **kwargs):
            flag = cmd[1] if len(cmd) > 1 else ""
            if flag == "-v":
                return MagicMock(returncode=0, stdout="v2.0.0\n")
            return MagicMock(returncode=1, stdout="")

        mock_run.side_effect = side_effect
        tool = self._setup_tool()
        version = tool.get_version()
        assert version == "v2.0.0"

    @patch("metadarkmatter.external.base.subprocess.run")
    def test_returns_none_when_all_fail(self, mock_run: MagicMock):
        """Should return None when no version flag produces output."""
        mock_run.return_value = MagicMock(returncode=1, stdout="")
        tool = self._setup_tool()
        version = tool.get_version()
        assert version is None

    @patch("metadarkmatter.external.base.subprocess.run")
    def test_handles_timeout_gracefully(self, mock_run: MagicMock):
        """Should return None when version commands time out."""
        mock_run.side_effect = subprocess.TimeoutExpired(
            cmd=["mocktool"], timeout=30
        )
        tool = self._setup_tool()
        version = tool.get_version()
        assert version is None

    @patch("metadarkmatter.external.base.subprocess.run")
    def test_handles_file_not_found_gracefully(self, mock_run: MagicMock):
        """Should return None when executable is missing at runtime."""
        mock_run.side_effect = FileNotFoundError("gone")
        tool = self._setup_tool()
        version = tool.get_version()
        assert version is None

    def test_handles_tool_not_found_gracefully(self):
        """Should return None when get_executable raises ToolNotFoundError."""
        _MockTool._executable_resolver = staticmethod(lambda name: None)
        _MockTool._executable_cache.clear()
        try:
            tool = _MockTool()
            version = tool.get_version()
            assert version is None
        finally:
            _MockTool.reset_executable_resolver()

    @patch("metadarkmatter.external.base.subprocess.run")
    def test_skips_empty_stdout(self, mock_run: MagicMock):
        """Should skip flags whose stdout is empty even with returncode 0."""
        call_count = 0

        def side_effect(cmd, **kwargs):
            nonlocal call_count
            call_count += 1
            if call_count <= 2:
                # First two flags return empty stdout with success
                return MagicMock(returncode=0, stdout="   \n")
            return MagicMock(returncode=0, stdout="version 3.0\n")

        mock_run.side_effect = side_effect
        tool = self._setup_tool()
        version = tool.get_version()
        assert version == "version 3.0"
