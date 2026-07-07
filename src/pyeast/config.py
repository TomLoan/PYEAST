"""Configuration management for PYEAST.

This module handles configuration loading and data path resolution.
Configuration priority (highest to lowest):
1. Environment variables (PYEAST_DATA_DIR, PYEAST_OUTPUT_DIR)
2. Config file (~/PYEAST/config.yaml)
3. Default locations (~/PYEAST/data/)
"""

import os
from pathlib import Path
from typing import Optional

import yaml


class PyeastConfig:
    """Configuration manager for PYEAST."""

    def __init__(self):
        """Initialize configuration by resolving paths in priority order."""
        self.data_dir = self._resolve_data_dir()
        self.output_dir = self._resolve_output_dir()
        self.preferred_templates = self._resolve_preferred_templates()

    def _resolve_preferred_templates(self) -> list:
        """Resolve the preferred-templates list from the config file.

        When a part is found in several templates, these are preferred (matched against the
        template record name inside the .gb, not the filename or part name). Empty by default,
        in which case the shortest-matching-template sorter alone decides.
        """
        config_data = self._load_config_file()
        if config_data and 'preferred_templates' in config_data:
            value = config_data['preferred_templates']
            if isinstance(value, list):
                return value
        return []

    def _resolve_data_dir(self) -> Path:
        """Resolve data directory location.

        Priority:
        1. PYEAST_DATA_DIR environment variable
        2. data_dir in config file
        3. ~/PYEAST/data/ as default
        """
        # Priority 1: Environment variable
        env_data_dir = os.getenv('PYEAST_DATA_DIR')
        if env_data_dir:
            path = Path(env_data_dir)
            if path.exists():
                return path.resolve()

        # Priority 2: Config file
        config_data = self._load_config_file()
        if config_data and 'data_dir' in config_data:
            path = Path(config_data['data_dir'])
            if path.exists():
                return path.resolve()

        # Priority 3: Default user directory. Fallback for legacy installs.
        # Not intended for normal use
        return (Path.home() / "PYEAST" / "data").resolve()

    def _resolve_output_dir(self) -> Path:
        """Resolve output directory location.

        Priority:
        1. PYEAST_OUTPUT_DIR environment variable
        2. output_dir in config file
        3. ~/PYEAST/output/ as default
        """
        # Priority 1: Environment variable
        env_output_dir = os.getenv('PYEAST_OUTPUT_DIR')
        if env_output_dir:
            return Path(env_output_dir).resolve()

        # Priority 2: Config file
        config_data = self._load_config_file()
        if config_data and 'output_dir' in config_data:
            return Path(config_data['output_dir']).resolve()

        # Priority 3: Default user directory  
        return (Path.home() / "PYEAST" / "output").resolve()

    def _load_config_file(self) -> Optional[dict]:
        """Load configuration from ~/PYEAST/config.yaml if it exists.

        Falls back to the legacy ~/.pyeast/config.yaml location (read-only) so
        installs created before the layout change keep working until the user
        re-runs `pyeast init`.
        """
        config_file = Path.home() / "PYEAST" / "config.yaml"
        if not config_file.exists():
            # Legacy migration fallback (read-only; never written to again)
            legacy_file = Path.home() / ".pyeast" / "config.yaml"
            if legacy_file.exists():
                config_file = legacy_file
        if config_file.exists():
            try:
                with open(config_file) as f:
                    return yaml.safe_load(f)
            except (OSError, yaml.YAMLError, PermissionError) as e:
                # Config file exists but couldn't be loaded - warn user
                import warnings
                warnings.warn(f"Could not load config file {config_file}: {e}")
                pass
        return None

# Singleton instance
_config: Optional[PyeastConfig] = None


def get_config() -> PyeastConfig:
    """Get or create the global config instance."""
    global _config
    if _config is None:
        _config = PyeastConfig()
    return _config


def reset_config():
    """Reset the global config instance (useful for testing)."""
    global _config
    _config = None
