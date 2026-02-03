"""Configuration management for PYEAST.

This module handles configuration loading and data path resolution.
Configuration priority (highest to lowest):
1. Environment variables (PYEAST_DATA_DIR, PYEAST_OUTPUT_DIR)
2. Config file (~/.pyeast/config.yaml)
3. Dev mode detection (./data/ if running from git checkout)
4. Default locations (~/.pyeast/data/)
"""

import os
import yaml
from pathlib import Path
from typing import Optional


class PyeastConfig:
    """Configuration manager for PYEAST."""

    def __init__(self):
        """Initialize configuration by resolving paths in priority order."""
        self.data_dir = self._resolve_data_dir()
        self.output_dir = self._resolve_output_dir()

    def _resolve_data_dir(self) -> Path:
        """Resolve data directory location.

        Priority:
        1. PYEAST_DATA_DIR environment variable
        2. data_dir in config file
        3. ./data/ if in dev mode (git checkout)
        4. ~/.pyeast/data/ as default
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

        # Priority 3: Dev mode (git checkout with ./data/)
        if self._is_dev_mode():
            return (Path.cwd() / "data").resolve()

        # Priority 4: Default user directory
        return (Path.home() / ".pyeast" / "data").resolve()

    def _resolve_output_dir(self) -> Path:
        """Resolve output directory location.

        Priority:
        1. PYEAST_OUTPUT_DIR environment variable
        2. output_dir in config file
        3. ./output/ if in dev mode
        4. ~/.pyeast/output/ as default
        """
        # Priority 1: Environment variable
        env_output_dir = os.getenv('PYEAST_OUTPUT_DIR')
        if env_output_dir:
            return Path(env_output_dir).resolve()

        # Priority 2: Config file
        config_data = self._load_config_file()
        if config_data and 'output_dir' in config_data:
            return Path(config_data['output_dir']).resolve()

        # Priority 3: Dev mode
        if self._is_dev_mode():
            return (Path.cwd() / "output").resolve()

        # Priority 4: Default user directory
        return (Path.home() / ".pyeast" / "output").resolve()

    def _load_config_file(self) -> Optional[dict]:
        """Load configuration from ~/.pyeast/config.yaml if it exists."""
        config_file = Path.home() / ".pyeast" / "config.yaml"
        if config_file.exists():
            try:
                with open(config_file, 'r') as f:
                    return yaml.safe_load(f)
            except Exception:
                # Silently ignore config file errors
                pass
        return None

    def _is_dev_mode(self) -> bool:
        """Check if running from a git checkout (dev mode).

        Returns True if:
        - Current directory contains ./data/ directory
        - Current directory contains .git/ directory or file
        """
        cwd = Path.cwd()
        has_data = (cwd / "data").exists()
        has_git = (cwd / ".git").exists()
        return has_data and has_git

    @property
    def is_dev_mode(self) -> bool:
        """Public property to check if in dev mode."""
        return self._is_dev_mode()


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
