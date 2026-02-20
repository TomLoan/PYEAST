"""Security tests for path validation in path_utils.

This module tests that path traversal attacks and other security issues
are properly prevented by the path utility functions.
"""


import pytest

from pyeast.config import reset_config
from pyeast.utils.path_utils import get_data_path, get_output_path, get_private_equivalent


class TestPathTraversalPrevention:
    """Test that path traversal attacks are blocked."""

    def test_get_data_path_blocks_parent_traversal(self):
        """Test that .. is blocked in get_data_path."""
        with pytest.raises(ValueError, match="Invalid subdirectory path"):
            get_data_path("../../../etc/passwd")

    def test_get_data_path_blocks_absolute_unix_path(self):
        """Test that absolute Unix paths are blocked."""
        with pytest.raises(ValueError, match="Invalid subdirectory path"):
            get_data_path("/etc/passwd")

    def test_get_data_path_blocks_absolute_windows_path(self):
        """Test that absolute Windows paths are blocked."""
        with pytest.raises(ValueError, match="Invalid subdirectory path"):
            get_data_path("\\Windows\\System32")

    def test_get_data_path_blocks_hidden_traversal(self):
        """Test that .. hidden in subdirectory is blocked."""
        with pytest.raises(ValueError, match="Invalid subdirectory path"):
            get_data_path("component libraries/../../evil")

    def test_get_output_path_blocks_parent_traversal(self):
        """Test that .. is blocked in get_output_path."""
        with pytest.raises(ValueError, match="Invalid subdirectory path"):
            get_output_path("../../../home/victim/.bashrc")

    def test_get_output_path_blocks_absolute_path(self):
        """Test that absolute paths are blocked in get_output_path."""
        with pytest.raises(ValueError, match="Invalid subdirectory path"):
            get_output_path("/tmp/evil")


class TestGetPrivateEquivalent:
    """Test get_private_equivalent security."""

    def test_private_equivalent_raises_on_external_path(self, tmp_path):
        """Test that paths outside data dir raise ValueError."""
        external_path = tmp_path / "external" / "malicious.fasta"

        with pytest.raises(ValueError, match="is not within data directory"):
            get_private_equivalent(external_path)

    def test_private_equivalent_works_for_valid_paths(self, test_data_dir):
        """Test that valid paths within data dir work correctly."""
        # This should work - path is within data dir
        public_path = test_data_dir / "component_libraries" / "YTK" / "promoter.fasta"
        private_path = get_private_equivalent(public_path)

        # Should be in private/component libraries/YTK/
        assert "private" in str(private_path)
        assert "promoter.fasta" in str(private_path)


class TestValidPathsStillWork:
    """Test that legitimate use cases are not broken by security fixes."""

    def test_get_data_path_normal_subdirectory(self):
        """Test that normal subdirectories work."""
        path = get_data_path("component_libraries")
        assert path.exists()
        assert "component_libraries" in str(path)

    def test_get_data_path_nested_subdirectory(self):
        """Test that nested subdirectories work."""
        path = get_data_path("component_libraries/YTK")
        # Path may not exist but should be constructed correctly
        assert "component_libraries" in str(path)
        assert "YTK" in str(path)

    def test_get_data_path_with_underscores(self):
        """Test that subdirectories with underscores work."""
        path = get_data_path("integration_sites")
        # Should construct path correctly
        assert "integration_sites" in str(path)

    def test_get_output_path_normal_use(self):
        """Test that normal output path usage works."""
        path = get_output_path("my_experiment")
        assert "my_experiment" in str(path)

    def test_get_data_path_empty_string(self):
        """Test that empty string returns base directory."""
        path = get_data_path("")
        assert path.exists()
        assert path.is_dir()

    def test_get_output_path_empty_string(self):
        """Test that empty string returns base output directory."""
        path = get_output_path("")
        # Output dir may not exist yet, but should be a valid path
        assert path.is_absolute()


class TestSymlinkAttackPrevention:
    """Test that symlink attacks are prevented."""

    def test_symlink_escape_blocked(self, tmp_path, monkeypatch):
        """Test that symlinks pointing outside base are blocked."""
        # Create a temporary data directory
        fake_data = tmp_path / "data"
        fake_data.mkdir()

        # Create a symlink pointing outside
        evil_link = fake_data / "evil_link"
        target = tmp_path / "outside"
        target.mkdir()

        try:
            evil_link.symlink_to(target)
        except OSError:
            # Windows may not allow symlinks without admin
            pytest.skip("Symlinks not available on this system")

        # Mock the data directory
        monkeypatch.setenv('PYEAST_DATA_DIR', str(fake_data))
        reset_config()

        # Try to access via symlink - should be blocked
        with pytest.raises(ValueError, match="escapes data directory"):
            get_data_path("evil_link")

        # Clean up
        reset_config()
