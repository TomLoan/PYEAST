"""Path resolution utilities for PYEAST.

This module provides helper functions for resolving data and output paths.
"""

from pathlib import Path

from pyeast.config import get_config


def get_data_path(subdirectory: str = "") -> Path:
    """Get path to data directory or subdirectory.

    Args:
        subdirectory: Optional subdirectory within data directory

    Returns:
        Resolved absolute path

    Raises:
        ValueError: If subdirectory contains path traversal attempts or escapes base directory

    Examples:
        >>> get_data_path()  # Base data directory
        >>> get_data_path("component_libraries")
        >>> get_data_path("templates")
    """
    config = get_config()
    base = config.data_dir

    if subdirectory:
        # Validate subdirectory doesn't contain path traversal attempts
        if ".." in subdirectory or subdirectory.startswith(("/", "\\")):
            raise ValueError(f"Invalid subdirectory path: {subdirectory}")

    result = base / subdirectory if subdirectory else base

    # Ensure result stays within base directory (catches symlinks too)
    try:
        result.resolve().relative_to(base.resolve())
    except ValueError:
        raise ValueError(f"Path {result} escapes data directory {base}")

    return result


def get_output_path(subdirectory: str = "") -> Path:
    """Get path to output directory or subdirectory.

    Args:
        subdirectory: Optional subdirectory within output directory

    Returns:
        Resolved absolute path

    Raises:
        ValueError: If subdirectory contains path traversal attempts or escapes base directory

    Examples:
        >>> get_output_path()  # Base output directory
        >>> get_output_path("project1")
    """
    config = get_config()
    base = config.output_dir

    if subdirectory:
        # Validate subdirectory doesn't contain path traversal attempts
        if ".." in subdirectory or subdirectory.startswith(("/", "\\")):
            raise ValueError(f"Invalid subdirectory path: {subdirectory}")

    result = base / subdirectory if subdirectory else base

    # Ensure result stays within base directory (catches symlinks too)
    try:
        result.resolve().relative_to(base.resolve())
    except ValueError:
        raise ValueError(f"Path {result} escapes output directory {base}")

    return result


def get_component_libraries_path(private: bool = False) -> Path:
    """Get path to component libraries directory.

    Args:
        private: If True, return path to private component libraries

    Returns:
        Path to component libraries directory
    """
    if private:
        return get_data_path("private/component_libraries")
    return get_data_path("component_libraries")


def get_integration_sites_path(private: bool = False) -> Path:
    """Get path to integration sites directory.

    Args:
        private: If True, return path to private integration sites

    Returns:
        Path to integration sites directory
    """
    if private:
        return get_data_path("private/integration_sites")
    return get_data_path("integration_sites")


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


def ensure_output_dir_exists(subdirectory: str = "") -> Path:
    """Ensure output directory exists, creating it if necessary.

    Args:
        subdirectory: Optional subdirectory within output directory

    Returns:
        Path to the output directory

    Raises:
        ValueError: If path exists but is not a directory
    """
    output_path = get_output_path(subdirectory)

    if output_path.exists() and not output_path.is_dir():
        raise ValueError(f"{output_path} exists but is not a directory")

    output_path.mkdir(parents=True, exist_ok=True)
    return output_path


def get_private_equivalent(public_path: Path) -> Path:
    """Convert a public data path to its private equivalent.

    Given a path like: /data/component_libraries/YeastToolKit
    Returns: /data/private/component_libraries/YeastToolKit

    Args:
        public_path: Path within public data directory

    Returns:
        Equivalent path in private data directory

    Raises:
        ValueError: If public_path is not within the data directory
    """
    config = get_config()
    data_dir = config.data_dir

    # Try to get the relative path from the data directory
    try:
        relative_path = public_path.relative_to(data_dir)
        return get_data_path("private") / relative_path
    except ValueError:
        # Path is not within data directory - this is an error
        raise ValueError(f"Path {public_path} is not within data directory {data_dir}")
