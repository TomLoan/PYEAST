"""Path resolution utilities for PYEAST.

This module provides helper functions for resolving data and output paths
in a way that works both in development mode (git checkout) and when
installed as a package.
"""

from pathlib import Path
from pyeast.config import get_config


def get_data_path(subdirectory: str = "") -> Path:
    """Get path to data directory or subdirectory.

    Args:
        subdirectory: Optional subdirectory within data directory

    Returns:
        Resolved absolute path

    Examples:
        >>> get_data_path()  # Base data directory
        >>> get_data_path("component libraries")
        >>> get_data_path("templates")
    """
    config = get_config()
    base = config.data_dir
    return base / subdirectory if subdirectory else base


def get_output_path(subdirectory: str = "") -> Path:
    """Get path to output directory or subdirectory.

    Args:
        subdirectory: Optional subdirectory within output directory

    Returns:
        Resolved absolute path

    Examples:
        >>> get_output_path()  # Base output directory
        >>> get_output_path("project1")
    """
    config = get_config()
    base = config.output_dir
    return base / subdirectory if subdirectory else base


def get_component_libraries_path(private: bool = False) -> Path:
    """Get path to component libraries directory.

    Args:
        private: If True, return path to private component libraries

    Returns:
        Path to component libraries directory
    """
    if private:
        return get_data_path("private/component libraries")
    return get_data_path("component libraries")


def get_integration_sites_path(private: bool = False) -> Path:
    """Get path to integration sites directory.

    Args:
        private: If True, return path to private integration sites

    Returns:
        Path to integration sites directory
    """
    if private:
        return get_data_path("private/integration sites")
    return get_data_path("integration sites")


def get_primers_path(private: bool = False) -> Path:
    """Get path to primers directory.

    Args:
        private: If True, return path to private primers

    Returns:
        Path to primers directory
    """
    if private:
        return get_data_path("private/primers")
    return get_data_path("primers")


def get_templates_path(private: bool = False) -> Path:
    """Get path to templates directory.

    Args:
        private: If True, return path to private templates

    Returns:
        Path to templates directory
    """
    if private:
        return get_data_path("private/templates")
    return get_data_path("templates")


def is_dev_mode() -> bool:
    """Check if running in development mode (git checkout).

    Returns:
        True if running from git checkout, False otherwise
    """
    config = get_config()
    return config.is_dev_mode


def ensure_output_dir_exists(subdirectory: str = "") -> Path:
    """Ensure output directory exists, creating it if necessary.

    Args:
        subdirectory: Optional subdirectory within output directory

    Returns:
        Path to the output directory
    """
    output_path = get_output_path(subdirectory)
    output_path.mkdir(parents=True, exist_ok=True)
    return output_path


def get_private_equivalent(public_path: Path) -> Path:
    """Convert a public data path to its private equivalent.

    Given a path like: /data/component libraries/YeastToolKit
    Returns: /data/private/component libraries/YeastToolKit

    Args:
        public_path: Path within public data directory

    Returns:
        Equivalent path in private data directory
    """
    config = get_config()
    data_dir = config.data_dir

    # Try to get the relative path from the data directory
    try:
        relative_path = public_path.relative_to(data_dir)
        return get_data_path("private") / relative_path
    except ValueError:
        # If public_path is not relative to data_dir, just prepend private/
        return get_data_path("private") / public_path.name
