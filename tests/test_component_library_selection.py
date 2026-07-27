"""Tests for component library discovery and selection.

Covers library_has_plasmids and the require_plasmids filter on the interactive
component directory picker, which keeps the Golden Gate designer from offering
libraries it cannot assemble from.
"""

import click
import pytest

from pyeast.cli.main import get_component_dir
from pyeast.config import reset_config
from pyeast.utils.path_utils import library_has_plasmids


@pytest.fixture
def library_tree(tmp_path, monkeypatch, mock_data_paths):
    """Build a data dir covering every public/private plasmid combination.

    Layout:
        with_plasmids       public plasmids folder holding a .gb
        empty_plasmids      public plasmids folder with no GenBank files
        no_plasmids         no plasmids folder at all
        private_plasmids    plasmids only under private/
    """
    public = tmp_path / "component_libraries"
    private = tmp_path / "private" / "component_libraries"

    (public / "with_plasmids" / "plasmids").mkdir(parents=True)
    (public / "with_plasmids" / "plasmids" / "p1.gb").write_text("dummy")

    (public / "empty_plasmids" / "plasmids").mkdir(parents=True)

    (public / "no_plasmids").mkdir(parents=True)

    (public / "private_plasmids").mkdir(parents=True)
    (private / "private_plasmids" / "plasmids").mkdir(parents=True)
    (private / "private_plasmids" / "plasmids" / "p2.gbk").write_text("dummy")

    monkeypatch.setenv('PYEAST_DATA_DIR', str(tmp_path))
    reset_config()

    yield tmp_path

    reset_config()


class TestLibraryHasPlasmids:
    """Test the gg-viability check for a component library."""

    def test_public_plasmids_with_genbank(self, library_tree):
        assert library_has_plasmids("with_plasmids") is True

    def test_empty_plasmids_folder_rejected(self, library_tree):
        assert library_has_plasmids("empty_plasmids") is False

    def test_missing_plasmids_folder_rejected(self, library_tree):
        assert library_has_plasmids("no_plasmids") is False

    def test_private_only_plasmids_accepted(self, library_tree):
        assert library_has_plasmids("private_plasmids") is True

    def test_unknown_library_rejected(self, library_tree):
        assert library_has_plasmids("does_not_exist") is False


class TestGetComponentDirFiltering:
    """Test the require_plasmids filter on the interactive picker."""

    @staticmethod
    def _answer(monkeypatch, *responses):
        """Feed canned answers to the prompt, one per call.

        PromptSession is replaced wholesale because constructing the real one
        needs a console, which pytest does not provide.
        """
        answers = iter(responses)

        class FakePromptSession:
            def prompt(self, *args, **kwargs):
                return next(answers)

        monkeypatch.setattr("pyeast.cli.main.PromptSession", FakePromptSession)

    def test_plasmid_free_library_rejected_when_required(self, library_tree, monkeypatch):
        """A library without plasmids cannot be selected, even if typed exactly."""
        self._answer(monkeypatch, "no_plasmids", "with_plasmids")

        result = get_component_dir(require_plasmids=True)

        assert result == library_tree / "component_libraries" / "with_plasmids"

    def test_plasmid_free_library_allowed_by_default(self, library_tree, monkeypatch):
        """Other designers still see every library."""
        self._answer(monkeypatch, "no_plasmids")

        result = get_component_dir()

        assert result == library_tree / "component_libraries" / "no_plasmids"

    def test_private_only_library_returns_public_path(self, library_tree, monkeypatch):
        """The public path is returned so load_sequences can resolve both sides."""
        self._answer(monkeypatch, "private_plasmids")

        result = get_component_dir(require_plasmids=True)

        assert result == library_tree / "component_libraries" / "private_plasmids"

    def test_selection_is_case_insensitive(self, library_tree, monkeypatch):
        """Matching agrees with the case-insensitive completer."""
        self._answer(monkeypatch, "WITH_PLASMIDS")

        result = get_component_dir(require_plasmids=True)

        assert result == library_tree / "component_libraries" / "with_plasmids"

    def test_aborts_when_no_library_has_plasmids(self, tmp_path, monkeypatch, mock_data_paths):
        """Golden Gate aborts up front rather than failing mid-design."""
        (tmp_path / "component_libraries" / "parts_only").mkdir(parents=True)
        monkeypatch.setenv('PYEAST_DATA_DIR', str(tmp_path))
        reset_config()

        try:
            with pytest.raises(click.Abort):
                get_component_dir(require_plasmids=True)
        finally:
            reset_config()
