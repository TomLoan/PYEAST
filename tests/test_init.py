"""Tests for pyeast init command and configuration system."""

import subprocess
import pytest
import yaml
from pathlib import Path
from click.testing import CliRunner
from unittest.mock import patch, MagicMock

from pyeast.cli.main import cli, DATA_REPO_URL, _DATA_REPO_CLONE_DIR
from pyeast.config import PyeastConfig, get_config, reset_config


@pytest.fixture
def runner():
    """Provide Click CLI test runner."""
    return CliRunner()


@pytest.fixture
def isolated_config(tmp_path, monkeypatch):
    """Create isolated config environment for testing.

    This fixture disables the autouse mock_data_paths fixture
    so we can test config behavior in isolation.
    """
    # Create temporary home directory
    fake_home = tmp_path / "home"
    fake_home.mkdir()
    monkeypatch.setattr(Path, "home", lambda: fake_home)

    # Create temporary data directory
    fake_data = tmp_path / "data"
    fake_data.mkdir()

    # Create a fake CWD that's NOT in dev mode
    fake_cwd = tmp_path / "fake_cwd"
    fake_cwd.mkdir()

    # Mock Path.cwd() to return non-dev directory
    monkeypatch.setattr(Path, "cwd", lambda: fake_cwd)

    # Clear any environment variables
    monkeypatch.delenv('PYEAST_DATA_DIR', raising=False)
    monkeypatch.delenv('PYEAST_OUTPUT_DIR', raising=False)

    # Reset config
    reset_config()

    yield {
        'home': fake_home,
        'data': fake_data,
        'tmp': tmp_path,
        'cwd': fake_cwd
    }

    # Clean up
    reset_config()


class TestInitCommand:
    """Test the pyeast init command."""

    def test_init_with_data_dir(self, runner, isolated_config):
        """Test init command with --data-dir option."""
        data_dir = isolated_config['data']

        result = runner.invoke(cli, ['init', '--data-dir', str(data_dir)])

        assert result.exit_code == 0
        assert "Configured PYEAST to use data at" in result.output

        config_file = isolated_config['home'] / ".pyeast" / "config.yaml"
        assert config_file.exists()

        with open(config_file, 'r') as f:
            config_data = yaml.safe_load(f)

        assert 'data_dir' in config_data
        assert Path(config_data['data_dir']) == data_dir.resolve()

    def test_init_with_data_and_output_dir(self, runner, isolated_config):
        """Test init command with both --data-dir and --output-dir."""
        data_dir = isolated_config['data']
        output_dir = isolated_config['tmp'] / "output"
        output_dir.mkdir()

        result = runner.invoke(cli, [
            'init',
            '--data-dir', str(data_dir),
            '--output-dir', str(output_dir)
        ])

        assert result.exit_code == 0
        assert "Configured PYEAST to use data at" in result.output
        assert "Output directory set to" in result.output

        config_file = isolated_config['home'] / ".pyeast" / "config.yaml"
        with open(config_file, 'r') as f:
            config_data = yaml.safe_load(f)

        assert 'data_dir' in config_data
        assert 'output_dir' in config_data
        assert Path(config_data['data_dir']) == data_dir.resolve()
        assert Path(config_data['output_dir']) == output_dir.resolve()

    def test_init_with_nonexistent_data_dir(self, runner, isolated_config):
        """Test init command with non-existent data directory."""
        bad_dir = isolated_config['tmp'] / "nonexistent"

        result = runner.invoke(cli, ['init', '--data-dir', str(bad_dir)])

        assert result.exit_code != 0
        assert "Directory not found" in result.output

    def test_init_with_nonexistent_output_dir(self, runner, isolated_config):
        """Test init command with non-existent output directory (should warn but succeed)."""
        data_dir = isolated_config['data']
        bad_output = isolated_config['tmp'] / "nonexistent_output"

        result = runner.invoke(cli, [
            'init',
            '--data-dir', str(data_dir),
            '--output-dir', str(bad_output)
        ])

        assert result.exit_code == 0
        assert "Output directory does not exist" in result.output
        assert "It will be created when needed" in result.output

    def test_init_already_configured_no_reconfigure(self, runner, isolated_config):
        """Test init with no args when already configured — user declines reconfigure."""
        data_dir = isolated_config['data']

        # Run init once to configure
        runner.invoke(cli, ['init', '--data-dir', str(data_dir)])
        reset_config()

        # Run init again without args; user answers 'n' to reconfigure prompt
        result = runner.invoke(cli, ['init'], input='n\n')

        assert result.exit_code == 0
        assert "Data directory:" in result.output
        assert "Would you like to reconfigure?" in result.output

    def test_init_already_configured_reconfigure_manual_path(self, runner, isolated_config):
        """Test reconfiguring with a new manual path when already configured."""
        data_dir = isolated_config['data']
        new_data_dir = isolated_config['tmp'] / "new_data"
        new_data_dir.mkdir()

        # Initial configuration
        runner.invoke(cli, ['init', '--data-dir', str(data_dir)])
        reset_config()

        # Reconfigure: choose option 1 (manual path)
        result = runner.invoke(cli, ['init'], input=f'y\n1\n{new_data_dir}\n')

        assert result.exit_code == 0
        assert "Configured PYEAST to use data at" in result.output

        config_file = isolated_config['home'] / ".pyeast" / "config.yaml"
        with open(config_file, 'r') as f:
            config_data = yaml.safe_load(f)
        assert Path(config_data['data_dir']) == new_data_dir.resolve()

    def test_init_no_args_clones_repo(self, runner, isolated_config):
        """Test init with no args and nothing configured — should clone data repo."""
        clone_dir = isolated_config['home'] / ".pyeast" / "data-repo"

        mock_result = MagicMock()
        mock_result.returncode = 0

        with patch('pyeast.cli.main._DATA_REPO_CLONE_DIR', clone_dir), \
             patch('subprocess.run', return_value=mock_result) as mock_run:
            result = runner.invoke(cli, ['init'])

        assert result.exit_code == 0, result.output
        mock_run.assert_called_once_with(
            ["git", "clone", DATA_REPO_URL, str(clone_dir)],
            capture_output=True, text=True
        )
        assert "Cloning PYEAST data repository" in result.output

        config_file = isolated_config['home'] / ".pyeast" / "config.yaml"
        assert config_file.exists()
        with open(config_file, 'r') as f:
            config_data = yaml.safe_load(f)
        assert Path(config_data['data_dir']) == clone_dir

    def test_init_no_args_clone_fails(self, runner, isolated_config):
        """Test init with no args when git clone fails."""
        clone_dir = isolated_config['home'] / ".pyeast" / "data-repo"

        mock_result = MagicMock()
        mock_result.returncode = 1
        mock_result.stderr = "fatal: repository not found"

        with patch('pyeast.cli.main._DATA_REPO_CLONE_DIR', clone_dir), \
             patch('subprocess.run', return_value=mock_result):
            result = runner.invoke(cli, ['init'])

        assert result.exit_code != 0
        assert "Clone failed" in result.output

    def test_init_no_args_already_cloned(self, runner, isolated_config):
        """Test init with no args when clone dir already exists — should not re-clone."""
        clone_dir = isolated_config['home'] / ".pyeast" / "data-repo"
        clone_dir.mkdir(parents=True)

        with patch('pyeast.cli.main._DATA_REPO_CLONE_DIR', clone_dir), \
             patch('subprocess.run') as mock_run:
            result = runner.invoke(cli, ['init'])

        assert result.exit_code == 0, result.output
        mock_run.assert_not_called()
        assert "already cloned" in result.output

        config_file = isolated_config['home'] / ".pyeast" / "config.yaml"
        with open(config_file, 'r') as f:
            config_data = yaml.safe_load(f)
        assert Path(config_data['data_dir']) == clone_dir

    def test_init_already_configured_reconfigure_clone(self, runner, isolated_config):
        """Test reconfiguring via clone when already configured with a different path."""
        data_dir = isolated_config['data']
        clone_dir = isolated_config['home'] / ".pyeast" / "data-repo"

        # Initial configuration
        runner.invoke(cli, ['init', '--data-dir', str(data_dir)])
        reset_config()

        mock_result = MagicMock()
        mock_result.returncode = 0

        with patch('pyeast.cli.main._DATA_REPO_CLONE_DIR', clone_dir), \
             patch('subprocess.run', return_value=mock_result) as mock_run:
            # Reconfigure: choose option 2 (clone)
            result = runner.invoke(cli, ['init'], input='y\n2\n')

        assert result.exit_code == 0, result.output
        mock_run.assert_called_once_with(
            ["git", "clone", DATA_REPO_URL, str(clone_dir)],
            capture_output=True, text=True
        )


class TestConfigPriority:
    """Test configuration priority system."""

    def test_env_var_takes_priority(self, isolated_config, monkeypatch):
        """Test that environment variable has highest priority."""
        config_file = isolated_config['home'] / ".pyeast" / "config.yaml"
        config_file.parent.mkdir(parents=True, exist_ok=True)

        config_data_dir = isolated_config['tmp'] / "config_data"
        config_data_dir.mkdir()

        with open(config_file, 'w') as f:
            yaml.dump({'data_dir': str(config_data_dir)}, f)

        env_data_dir = isolated_config['tmp'] / "env_data"
        env_data_dir.mkdir()
        monkeypatch.setenv('PYEAST_DATA_DIR', str(env_data_dir))

        reset_config()
        config = get_config()

        assert config.data_dir == env_data_dir.resolve()

    def test_config_file_used_when_no_env_var(self, isolated_config, monkeypatch):
        """Test that config file is used when no environment variable."""
        monkeypatch.delenv('PYEAST_DATA_DIR', raising=False)

        config_file = isolated_config['home'] / ".pyeast" / "config.yaml"
        config_file.parent.mkdir(parents=True, exist_ok=True)

        config_data_dir = isolated_config['tmp'] / "config_data"
        config_data_dir.mkdir()

        with open(config_file, 'w') as f:
            yaml.dump({'data_dir': str(config_data_dir)}, f)

        reset_config()
        config = get_config()

        assert config.data_dir == config_data_dir.resolve()

    def test_default_path_when_nothing_configured(self, isolated_config, monkeypatch):
        """Test default path is used when nothing else configured."""
        monkeypatch.delenv('PYEAST_DATA_DIR', raising=False)

        reset_config()
        config = get_config()

        expected = (isolated_config['home'] / "PYEAST" / "data").resolve()
        assert config.data_dir == expected


class TestConfigOutputDir:
    """Test output directory configuration."""

    def test_output_dir_env_var(self, isolated_config, monkeypatch):
        """Test PYEAST_OUTPUT_DIR environment variable."""
        output_dir = isolated_config['tmp'] / "my_output"
        output_dir.mkdir()

        monkeypatch.setenv('PYEAST_OUTPUT_DIR', str(output_dir))

        reset_config()
        config = get_config()

        assert config.output_dir == output_dir.resolve()

    def test_output_dir_config_file(self, isolated_config, monkeypatch):
        """Test output_dir from config file."""
        monkeypatch.delenv('PYEAST_OUTPUT_DIR', raising=False)

        config_file = isolated_config['home'] / ".pyeast" / "config.yaml"
        config_file.parent.mkdir(parents=True, exist_ok=True)

        data_dir = isolated_config['data']
        output_dir = isolated_config['tmp'] / "config_output"
        output_dir.mkdir()

        with open(config_file, 'w') as f:
            yaml.dump({
                'data_dir': str(data_dir),
                'output_dir': str(output_dir)
            }, f)

        reset_config()
        config = get_config()

        assert config.output_dir == output_dir.resolve()

    def test_output_dir_default(self, isolated_config, monkeypatch):
        """Test default output directory."""
        monkeypatch.delenv('PYEAST_OUTPUT_DIR', raising=False)

        reset_config()
        config = get_config()

        expected = (isolated_config['home'] / "PYEAST" / "output").resolve()
        assert config.output_dir == expected


class TestConfigFileHandling:
    """Test config file reading and error handling."""

    def test_invalid_yaml_emits_warning(self, isolated_config):
        """Test that invalid YAML in config file emits warning and falls back to defaults."""
        config_file = isolated_config['home'] / ".pyeast" / "config.yaml"
        config_file.parent.mkdir(parents=True, exist_ok=True)

        with open(config_file, 'w') as f:
            f.write("this is not: valid: yaml: [[[")

        with pytest.warns(UserWarning, match="Could not load config file"):
            reset_config()
            config = get_config()

        expected = (isolated_config['home'] / "PYEAST" / "data").resolve()
        assert config.data_dir == expected

    def test_missing_config_file_ok(self, isolated_config):
        """Test that missing config file is OK."""
        config_file = isolated_config['home'] / ".pyeast" / "config.yaml"
        assert not config_file.exists()

        reset_config()
        config = get_config()

        expected = (isolated_config['home'] / "PYEAST" / "data").resolve()
        assert config.data_dir == expected

    def test_partial_config_file(self, isolated_config):
        """Test config file with only some fields."""
        config_file = isolated_config['home'] / ".pyeast" / "config.yaml"
        config_file.parent.mkdir(parents=True, exist_ok=True)

        data_dir = isolated_config['data']

        with open(config_file, 'w') as f:
            yaml.dump({'data_dir': str(data_dir)}, f)

        reset_config()
        config = get_config()

        assert config.data_dir == data_dir.resolve()

        expected_output = (isolated_config['home'] / "PYEAST" / "output").resolve()
        assert config.output_dir == expected_output
