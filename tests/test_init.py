"""Tests for pyeast init command and configuration system."""

import pytest
import yaml
from pathlib import Path
from click.testing import CliRunner
from unittest.mock import patch

from pyeast.cli.main import cli
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

        # Run init command
        result = runner.invoke(cli, ['init', '--data-dir', str(data_dir)])

        assert result.exit_code == 0
        assert "Configured PYEAST to use data at" in result.output

        # Check config file was created
        config_file = isolated_config['home'] / ".pyeast" / "config.yaml"
        assert config_file.exists()

        # Verify config content
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

        # Check config file
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

    def test_init_already_configured(self, runner, isolated_config):
        """Test init command when already configured."""
        data_dir = isolated_config['data']

        # Run init once
        runner.invoke(cli, ['init', '--data-dir', str(data_dir)])

        # Reset config to re-read from file
        reset_config()

        # Run init again without arguments
        result = runner.invoke(cli, ['init'])

        assert result.exit_code == 0
        assert "Data directory already configured" in result.output

    def test_init_no_args_shows_help(self, runner, isolated_config):
        """Test init command without arguments shows helpful message."""
        result = runner.invoke(cli, ['init'])

        assert result.exit_code == 0
        assert "PYEAST Data Setup" in result.output
        assert "Options:" in result.output
        assert "pyeast init --data-dir" in result.output

    def test_init_in_dev_mode(self, runner, isolated_config, monkeypatch):
        """Test init command when in dev mode."""
        # Create ./data/ and .git in the fake cwd
        dev_cwd = isolated_config['cwd']
        (dev_cwd / "data").mkdir()
        (dev_cwd / ".git").mkdir()

        # Reset config so it picks up the new dev mode
        reset_config()

        result = runner.invoke(cli, ['init'])

        assert result.exit_code == 0
        assert "Running in dev mode" in result.output


class TestConfigPriority:
    """Test configuration priority system."""

    def test_env_var_takes_priority(self, isolated_config, monkeypatch):
        """Test that environment variable has highest priority."""
        # Set up config file
        config_file = isolated_config['home'] / ".pyeast" / "config.yaml"
        config_file.parent.mkdir(parents=True, exist_ok=True)

        config_data_dir = isolated_config['tmp'] / "config_data"
        config_data_dir.mkdir()

        with open(config_file, 'w') as f:
            yaml.dump({'data_dir': str(config_data_dir)}, f)

        # Set environment variable to different path
        env_data_dir = isolated_config['tmp'] / "env_data"
        env_data_dir.mkdir()
        monkeypatch.setenv('PYEAST_DATA_DIR', str(env_data_dir))

        # Reset and get config
        reset_config()
        config = get_config()

        # Environment variable should win
        assert config.data_dir == env_data_dir.resolve()

    def test_config_file_used_when_no_env_var(self, isolated_config, monkeypatch):
        """Test that config file is used when no environment variable."""
        # Ensure no env var
        monkeypatch.delenv('PYEAST_DATA_DIR', raising=False)

        # Create config file
        config_file = isolated_config['home'] / ".pyeast" / "config.yaml"
        config_file.parent.mkdir(parents=True, exist_ok=True)

        config_data_dir = isolated_config['tmp'] / "config_data"
        config_data_dir.mkdir()

        with open(config_file, 'w') as f:
            yaml.dump({'data_dir': str(config_data_dir)}, f)

        # Reset and get config
        reset_config()
        config = get_config()

        # Config file should be used
        assert config.data_dir == config_data_dir.resolve()

    @patch('pyeast.config.PyeastConfig._is_dev_mode')
    def test_dev_mode_used_when_no_config(self, mock_dev_mode, isolated_config, monkeypatch):
        """Test that dev mode is detected when no other config."""
        # Ensure no env var or config file
        monkeypatch.delenv('PYEAST_DATA_DIR', raising=False)

        # Mock dev mode detection
        mock_dev_mode.return_value = True

        # Reset and get config
        reset_config()
        config = get_config()

        # Should use ./data/ from current directory
        expected = (Path.cwd() / "data").resolve()
        assert config.data_dir == expected

    def test_default_path_when_nothing_configured(self, isolated_config, monkeypatch):
        """Test default path is used when nothing else configured."""
        # Ensure no env var or config file
        monkeypatch.delenv('PYEAST_DATA_DIR', raising=False)

        # Ensure we're not in dev mode (no ./data/ or .git)
        # This is implicitly true in tmp_path

        # Reset and get config
        reset_config()
        config = get_config()

        # Should use default: ~/PYEAST/data/
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

        # Create config file with output_dir
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

        # Default should be ~/PYEAST/output/
        expected = (isolated_config['home'] / "PYEAST" / "output").resolve()
        assert config.output_dir == expected


class TestDevModeDetection:
    """Test development mode detection."""

    def test_dev_mode_detected_with_data_and_git(self, monkeypatch):
        """Test dev mode is detected when ./data/ and .git exist."""
        # Clear env var set by autouse fixture
        monkeypatch.delenv('PYEAST_DATA_DIR', raising=False)
        monkeypatch.delenv('PYEAST_OUTPUT_DIR', raising=False)

        # We're already IN a dev mode directory (PYEAST Local has ./data/ and .git/)
        # So just test that config detects it correctly
        reset_config()
        config = PyeastConfig()

        # The actual PYEAST Local directory has both ./data/ and .git/
        # so dev mode should be detected
        assert config.is_dev_mode is True

        # And data_dir should point to ./data/
        assert config.data_dir.name == "data"
        assert config.data_dir.exists()

    def test_dev_mode_not_detected_without_git(self, tmp_path, monkeypatch):
        """Test dev mode is not detected without .git."""
        # Create only ./data/
        data_dir = tmp_path / "data"
        data_dir.mkdir()

        # Mock cwd
        monkeypatch.setattr(Path, "cwd", lambda: tmp_path)

        reset_config()
        config = PyeastConfig()

        assert config.is_dev_mode is False

    def test_dev_mode_not_detected_without_data(self, tmp_path, monkeypatch):
        """Test dev mode is not detected without ./data/."""
        # Create only .git
        git_dir = tmp_path / ".git"
        git_dir.mkdir()

        # Mock cwd
        monkeypatch.setattr(Path, "cwd", lambda: tmp_path)

        reset_config()
        config = PyeastConfig()

        assert config.is_dev_mode is False


class TestConfigFileHandling:
    """Test config file reading and error handling."""

    def test_invalid_yaml_ignored(self, isolated_config):
        """Test that invalid YAML in config file is silently ignored."""
        config_file = isolated_config['home'] / ".pyeast" / "config.yaml"
        config_file.parent.mkdir(parents=True, exist_ok=True)

        # Write invalid YAML
        with open(config_file, 'w') as f:
            f.write("this is not: valid: yaml: [[[")

        # Should not crash, should use default
        reset_config()
        config = get_config()

        # Should use default path
        expected = (isolated_config['home'] / "PYEAST" / "data").resolve()
        assert config.data_dir == expected

    def test_missing_config_file_ok(self, isolated_config):
        """Test that missing config file is OK."""
        # Ensure config file doesn't exist
        config_file = isolated_config['home'] / ".pyeast" / "config.yaml"
        assert not config_file.exists()

        # Should not crash
        reset_config()
        config = get_config()

        # Should use default
        expected = (isolated_config['home'] / "PYEAST" / "data").resolve()
        assert config.data_dir == expected

    def test_partial_config_file(self, isolated_config):
        """Test config file with only some fields."""
        config_file = isolated_config['home'] / ".pyeast" / "config.yaml"
        config_file.parent.mkdir(parents=True, exist_ok=True)

        data_dir = isolated_config['data']

        # Write config with only data_dir
        with open(config_file, 'w') as f:
            yaml.dump({'data_dir': str(data_dir)}, f)

        reset_config()
        config = get_config()

        # data_dir should come from config
        assert config.data_dir == data_dir.resolve()

        # output_dir should use default
        expected_output = (isolated_config['home'] / "PYEAST" / "output").resolve()
        assert config.output_dir == expected_output
