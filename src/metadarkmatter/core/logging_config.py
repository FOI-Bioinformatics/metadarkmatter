"""Centralized logging configuration for metadarkmatter."""

from __future__ import annotations

import json
import logging
import sys
from datetime import datetime, timezone


class _JSONFormatter(logging.Formatter):
    """Emit each log record as a single JSON object."""

    def format(self, record: logging.LogRecord) -> str:
        return json.dumps(
            {
                "timestamp": datetime.fromtimestamp(
                    record.created, tz=timezone.utc
                ).isoformat(),
                "level": record.levelname,
                "logger": record.name,
                "message": record.getMessage(),
            }
        )


_TEXT_FORMAT = "%(asctime)s %(levelname)-8s %(name)s: %(message)s"


def setup_logging(
    level: str = "WARNING",
    log_format: str = "text",
    log_file: str | None = None,
) -> None:
    """Configure the root logger for the application.

    Args:
        level: Logging level name (DEBUG, INFO, WARNING, ERROR).
        log_format: Output format - "text" for human-readable, "json" for
            machine-parseable structured output.
        log_file: Optional path to a log file. Logs are written to both
            stderr and the file when provided.
    """
    root = logging.getLogger()
    root.setLevel(getattr(logging, level.upper(), logging.WARNING))
    root.handlers.clear()

    formatter: logging.Formatter
    if log_format == "json":
        formatter = _JSONFormatter()
    else:
        formatter = logging.Formatter(_TEXT_FORMAT)

    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setFormatter(formatter)
    root.addHandler(stderr_handler)

    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        root.addHandler(file_handler)
