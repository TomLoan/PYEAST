# The Software is copyright (c) Commonwealth Scientific and Industrial Research Organisation (CSIRO) 2023-2024.
"""
This is a config file (written in python) used for runtime configuration.
"""

# see https://loguru.readthedocs.io/en/stable/api/type_hints.html#module-autodoc_stub_file.loguru
from __future__ import (
    annotations,
)

import os
import sys
from pathlib import Path

import loguru
from loguru import logger

out_dir = "./output"

# Get entry point file name as default log name
default_log_name = Path(sys.argv[0]).stem
default_log_name = "log" if default_log_name == "" else default_log_name


def get_logger(log_name: str = default_log_name, log_dir: str = ".") -> loguru.Logger:
    """Return a configured loguru logger.

    Call this once from entrypoints to set up a new logger.
    In non-entrypoint modules, just use `from loguru import logger` directly.

    To set the log level, use the `LOGURU_LEVEL` environment variable before or during runtime. E.g. `os.environ["LOGURU_LEVEL"] = "INFO"`
    Available levels are `TRACE`, `DEBUG`, `INFO`, `SUCCESS`, `WARNING`, `ERROR`, and `CRITICAL`. Default is `INFO`.

    Log file will be written to `f"{log_dir}/{log_name}.log"`

    See https://github.com/Delgan/loguru#suitable-for-scripts-and-libraries
    From loguru import Record, RecordFile # See these classes for all the available format strings

    Parameters:
        log_name (str): Name of the log. Corresponding log file will be called {log_name}.log in the .
        log_dir (str): Directory to write the log file to. Default is the current working directory.
    Returns:
        Logger: A configured loguru logger.
    """

    # set global log level via env var.  Set to INFO if not already set.
    if os.getenv("LOGURU_LEVEL") is None:
        os.environ["LOGURU_LEVEL"] = "INFO"

    format_str = (
        "<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | <level>{level: <8}</level> : <level>{message}</level> "
        '(<cyan>{name}:{thread.name}:pid-{process}</cyan> "<cyan>{file.path}</cyan>:<cyan>{line}</cyan>")'
    )
    log_config = {
        "handlers": [
            {
                "sink": sys.stdout,
                "level": "DEBUG",
                "diagnose": False,
                "format": format_str,
            },
            {
                "sink": out_dir + f"/log{os.getpid()}.log",
                "enqueue": True,
                "mode": "a+",
                "level": "DEBUG",
                "colorize": False,
                "serialize": True,
                "diagnose": False,
                "rotation": "10 MB",
                "compression": "zip",
            },
        ]
    }
    logger.configure(**log_config)
    return logger
