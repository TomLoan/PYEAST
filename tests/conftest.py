"""Pytest configuration and fixtures for PYEAST tests.

This module configures the test environment to use test fixtures
instead of live data directories.
"""

import pytest
from pathlib import Path
from unittest.mock import patch

import sys
import os

# Add src directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


@pytest.fixture(scope="session")
def test_data_dir():
    """Get the path to test fixtures directory."""
    return Path(__file__).parent / "fixtures"


@pytest.fixture(scope="function", autouse=True)
def mock_data_paths(request, test_data_dir, monkeypatch):
    """Mock all data path functions to use test fixtures.

    This fixture automatically applies to all tests and ensures
    that data paths resolve to the test fixtures directory.
    Tests marked with @pytest.mark.realdata are exempt and use the
    real configured data directory instead.
    """
    if request.node.get_closest_marker('realdata'):
        yield
        return

    # Mock environment variable
    monkeypatch.setenv('PYEAST_DATA_DIR', str(test_data_dir))

    # Mock the config to always use test data
    from pyeast.config import get_config, reset_config

    # Reset config to pick up the mocked environment variable
    reset_config()

    # Verify it worked
    config = get_config()
    assert config.data_dir == test_data_dir, f"Config data_dir is {config.data_dir}, expected {test_data_dir}"

    yield

    # Clean up after test
    reset_config()


@pytest.fixture
def temp_output_dir(tmp_path):
    """Provide a temporary output directory for tests."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    return output_dir


@pytest.fixture(scope="function")
def mock_output_path(temp_output_dir, monkeypatch):
    """Mock output path to use temporary directory."""
    monkeypatch.setenv('PYEAST_OUTPUT_DIR', str(temp_output_dir))

    from pyeast.config import reset_config
    reset_config()

    yield temp_output_dir

    reset_config()
