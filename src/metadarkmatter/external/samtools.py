"""
Samtools wrapper class.

Provides Python interface for common samtools operations:
- view: Convert between SAM/BAM formats
- sort: Sort alignments by coordinate or read name
- index: Create BAM index (.bai) for random access
"""

from __future__ import annotations

from pathlib import Path

from metadarkmatter.external.base import ExternalTool, ToolResult


class Samtools(ExternalTool):
    """Wrapper for samtools utilities.

    Provides methods for common BAM file operations without
    requiring pysam as a dependency.

    Example:
        >>> samtools = Samtools()
        >>> # Convert SAM to BAM
        >>> samtools.view(
        ...     input_file=Path("aligned.sam"),
        ...     output_file=Path("aligned.bam"),
        ...     output_bam=True,
        ... )
        >>> # Sort BAM
        >>> samtools.sort(
        ...     input_file=Path("aligned.bam"),
        ...     output_file=Path("aligned.sorted.bam"),
        ... )
        >>> # Index BAM
        >>> samtools.index(bam_file=Path("aligned.sorted.bam"))
    """

    TOOL_NAME = "samtools"
    INSTALL_HINT = "conda install -c bioconda samtools"

    def build_command(self, **kwargs: object) -> list[str]:
        """Build samtools command.

        This base method is not used directly. Use view(), sort(),
        index() methods instead.
        """
        # Subcommand-specific building is done in individual methods
        raise NotImplementedError(
            "Use view(), sort(), or index() methods instead"
        )

    def _build_view_command(
        self,
        *,
        input_file: Path,
        output_file: Path,
        output_bam: bool = True,
        threads: int = 4,
        include_header: bool = True,
        min_mapq: int | None = None,
        flags_required: int | None = None,
        flags_excluded: int | None = None,
        regions: list[str] | None = None,
    ) -> list[str]:
        """Build samtools view command."""
        exe = str(self.get_executable())
        cmd = [exe, "view"]

        # Output format
        if output_bam:
            cmd.append("-b")

        # Include header
        if include_header:
            cmd.append("-h")

        # Threading
        cmd.extend(["--threads", str(threads)])

        # Quality filter
        if min_mapq is not None:
            cmd.extend(["-q", str(min_mapq)])

        # Flag filters
        if flags_required is not None:
            cmd.extend(["-f", str(flags_required)])

        if flags_excluded is not None:
            cmd.extend(["-F", str(flags_excluded)])

        # Output file
        cmd.extend(["-o", str(output_file)])

        # Input file
        cmd.append(str(input_file))

        # Optional regions
        if regions:
            cmd.extend(regions)

        return cmd

    def view(
        self,
        *,
        input_file: Path,
        output_file: Path,
        output_bam: bool = True,
        threads: int = 4,
        include_header: bool = True,
        min_mapq: int | None = None,
        flags_required: int | None = None,
        flags_excluded: int | None = None,
        regions: list[str] | None = None,
        timeout: float | None = None,
        dry_run: bool = False,
    ) -> ToolResult:
        """Convert between SAM/BAM formats with optional filtering.

        Args:
            input_file: Input SAM/BAM file.
            output_file: Output file path.
            output_bam: Output BAM format (vs SAM).
            threads: Number of threads.
            include_header: Include header in output.
            min_mapq: Minimum mapping quality to include.
            flags_required: SAM flags that must be set.
            flags_excluded: SAM flags that must not be set.
            regions: Genomic regions to extract.
            timeout: Maximum execution time.
            dry_run: Return command without execution.

        Returns:
            ToolResult with command and output.
        """
        cmd = self._build_view_command(
            input_file=input_file,
            output_file=output_file,
            output_bam=output_bam,
            threads=threads,
            include_header=include_header,
            min_mapq=min_mapq,
            flags_required=flags_required,
            flags_excluded=flags_excluded,
            regions=regions,
        )

        if dry_run:
            return ToolResult(
                command=tuple(cmd),
                return_code=0,
                stdout="[dry-run] Command not executed",
                stderr="",
                elapsed_seconds=0.0,
            )

        return self._run_command(cmd, timeout=timeout)

    def _build_sort_command(
        self,
        *,
        input_file: Path,
        output_file: Path,
        threads: int = 4,
        memory_per_thread: str = "768M",
        by_name: bool = False,
        temp_prefix: Path | None = None,
    ) -> list[str]:
        """Build samtools sort command."""
        exe = str(self.get_executable())
        cmd = [exe, "sort"]

        # Threading
        cmd.extend(["--threads", str(threads)])

        # Memory per thread
        cmd.extend(["-m", memory_per_thread])

        # Sort order
        if by_name:
            cmd.append("-n")

        # Temp file prefix
        if temp_prefix is not None:
            cmd.extend(["-T", str(temp_prefix)])

        # Output
        cmd.extend(["-o", str(output_file)])

        # Input
        cmd.append(str(input_file))

        return cmd

    def sort(
        self,
        *,
        input_file: Path,
        output_file: Path,
        threads: int = 4,
        memory_per_thread: str = "768M",
        by_name: bool = False,
        temp_prefix: Path | None = None,
        timeout: float | None = None,
        dry_run: bool = False,
    ) -> ToolResult:
        """Sort alignments by coordinate or read name.

        Args:
            input_file: Input BAM/SAM file.
            output_file: Output sorted BAM file.
            threads: Number of threads.
            memory_per_thread: Memory per sorting thread.
            by_name: Sort by read name (vs coordinate).
            temp_prefix: Prefix for temporary files.
            timeout: Maximum execution time.
            dry_run: Return command without execution.

        Returns:
            ToolResult with command and output.
        """
        cmd = self._build_sort_command(
            input_file=input_file,
            output_file=output_file,
            threads=threads,
            memory_per_thread=memory_per_thread,
            by_name=by_name,
            temp_prefix=temp_prefix,
        )

        if dry_run:
            return ToolResult(
                command=tuple(cmd),
                return_code=0,
                stdout="[dry-run] Command not executed",
                stderr="",
                elapsed_seconds=0.0,
            )

        return self._run_command(cmd, timeout=timeout)

    def _build_index_command(
        self,
        *,
        bam_file: Path,
        threads: int = 4,
        csi: bool = False,
    ) -> list[str]:
        """Build samtools index command."""
        exe = str(self.get_executable())
        cmd = [exe, "index"]

        # Threading
        cmd.extend(["--threads", str(threads)])

        # Index type
        if csi:
            cmd.append("-c")

        # Input BAM
        cmd.append(str(bam_file))

        return cmd

    def index(
        self,
        *,
        bam_file: Path,
        threads: int = 4,
        csi: bool = False,
        timeout: float | None = None,
        dry_run: bool = False,
    ) -> ToolResult:
        """Create BAM index for random access.

        Args:
            bam_file: Input sorted BAM file.
            threads: Number of threads.
            csi: Create CSI index (vs BAI, for large chromosomes).
            timeout: Maximum execution time.
            dry_run: Return command without execution.

        Returns:
            ToolResult with command and output.
        """
        cmd = self._build_index_command(
            bam_file=bam_file,
            threads=threads,
            csi=csi,
        )

        if dry_run:
            return ToolResult(
                command=tuple(cmd),
                return_code=0,
                stdout="[dry-run] Command not executed",
                stderr="",
                elapsed_seconds=0.0,
            )

        return self._run_command(cmd, timeout=timeout)

    def _run_command(
        self,
        cmd: list[str],
        timeout: float | None = None,
    ) -> ToolResult:
        """Execute a samtools command."""
        import subprocess
        import time

        from metadarkmatter.external.base import (
            ToolExecutionError,
            ToolTimeoutError,
        )

        start_time = time.perf_counter()

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
            elapsed = time.perf_counter() - start_time

            tool_result = ToolResult(
                command=tuple(cmd),
                return_code=result.returncode,
                stdout=result.stdout,
                stderr=result.stderr,
                elapsed_seconds=elapsed,
            )

            if not tool_result.success:
                raise ToolExecutionError(
                    self.TOOL_NAME,
                    cmd,
                    result.returncode,
                    result.stderr,
                )

            return tool_result

        except subprocess.TimeoutExpired as e:
            raise ToolTimeoutError(
                self.TOOL_NAME,
                timeout or 0,
                cmd,
            ) from e

    def sam_to_sorted_bam(
        self,
        *,
        input_sam: Path,
        output_bam: Path,
        threads: int = 4,
        cleanup_sam: bool = False,
        timeout: float | None = None,
        dry_run: bool = False,
    ) -> tuple[ToolResult, ToolResult, ToolResult]:
        """Convert SAM to sorted, indexed BAM in one call.

        Convenience method that runs view -> sort -> index pipeline.

        Args:
            input_sam: Input SAM file.
            output_bam: Output sorted BAM file.
            threads: Number of threads for each step.
            cleanup_sam: Delete input SAM after conversion.
            timeout: Maximum time per step.
            dry_run: Return commands without execution.

        Returns:
            Tuple of (view_result, sort_result, index_result).
        """
        # Intermediate unsorted BAM
        unsorted_bam = output_bam.with_suffix(".unsorted.bam")

        # SAM -> BAM
        view_result = self.view(
            input_file=input_sam,
            output_file=unsorted_bam,
            output_bam=True,
            threads=threads,
            timeout=timeout,
            dry_run=dry_run,
        )

        # Sort BAM
        sort_result = self.sort(
            input_file=unsorted_bam,
            output_file=output_bam,
            threads=threads,
            timeout=timeout,
            dry_run=dry_run,
        )

        # Index BAM
        index_result = self.index(
            bam_file=output_bam,
            threads=threads,
            timeout=timeout,
            dry_run=dry_run,
        )

        # Cleanup intermediate files
        if not dry_run:
            if unsorted_bam.exists():
                unsorted_bam.unlink()
            if cleanup_sam and input_sam.exists():
                input_sam.unlink()

        return view_result, sort_result, index_result
