# The Software is copyright (c) Commonwealth Scientific and Industrial Research Organisation (CSIRO) 2023-2024.
import argparse
import fnmatch
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from pprint import pformat
from stat import S_IWRITE
from typing import Any, Callable

import pathspec
from dotenv import load_dotenv
from loguru import logger

# File patterns to search
file_patterns: list[str] = ["*.py", "*.yml", "*.yaml", "*.toml", "*.md", "*.txt", "Dockerfile", "*.txt", "*.rst"]


def get_substitution_map() -> dict[str, str]:
    """Initialise the project by loading the environment variables etc
    :return: dictionary of find-replace regexes
    """
    load_dotenv(verbose=True)

    logger.info(f"Initialising project with env: \n{pformat(dict({k: v for k, v in os.environ.items() if "NEW_" in k}))}")

    project_name: str = os.getenv("NEW_PROJECT_NAME").replace("_", "-")  # type: ignore
    project_slug: str = os.getenv("NEW_PROJECT_NAME").replace("-", "_")  # type: ignore
    project_author: str = os.getenv("NEW_PROJECT_AUTHOR")  # type: ignore
    project_email: str = os.getenv("NEW_PROJECT_EMAIL")  # type: ignore

    # Dictionary of find-replace regexes, set from environment variables
    substitutions: dict[str, str] = {
        r"python-template-uv": project_name,
        r"python_template_uv": project_slug,
        "author_name_here": project_author,
        "author_name_here@domain.com": project_email,
    }

    if not os.getenv("NEW_PROJECT_NAME") or not os.getenv("NEW_PROJECT_AUTHOR") or not os.getenv("NEW_PROJECT_EMAIL"):
        logger.error(
            f"Environment variables not set. Expecting NEW_PROJECT_NAME, NEW_PROJECT_AUTHOR, and NEW_PROJECT_EMAIL, but got "
            f"{os.getenv("NEW_PROJECT_NAME")}, {os.getenv("NEW_PROJECT_AUTHOR")}, {os.getenv("NEW_PROJECT_EMAIL")}. See Readme for details."
        )
        sys.exit(1)
    return substitutions


def load_gitignore(root_dir: str) -> pathspec.PathSpec:
    """Load the .gitignore file from the root directory."""
    gitignore_path = Path(root_dir) / ".gitignore"
    if gitignore_path.exists():
        with gitignore_path.open("r") as f:
            return pathspec.PathSpec.from_lines("gitwildmatch", f)
    else:
        raise FileNotFoundError(f"No .gitignore file found in {root_dir}")


def replace_in_file(file_path: Path, substitutions: dict[str, str], dry_run: bool, verbose: bool = False) -> None:
    """Replace all occurrences of find with replace in file."""
    if verbose:
        logger.debug(f"Processing {file_path}")

    content: str = file_path.read_text(encoding="utf-8", errors="ignore")
    for find, replace in substitutions.items():
        new_content, n_replaced = re.subn(find, replace, content)
        if (verbose or dry_run) and new_content != content:
            logger.info(f"Replaced {n_replaced} occurrences of {find} with {replace} in {file_path}")
        content = new_content
    if not dry_run:
        file_path.write_text(content, encoding="utf-8", errors="ignore")


def replace_str(orig_str: str, substitutions: dict[str, str]) -> str:
    """Replace all occurrences of find with replace in name."""
    for find, replace in substitutions.items():
        orig_str = re.sub(find, replace, orig_str)
    return orig_str


def process_directory(
    root_dir: str, subs_map: dict[str, str], file_patterns: list[str], dry_run: bool, undo: bool, untrack: bool, verbose: bool = False
) -> None:
    """Recursively process all files in the directory, applying the substitutions."""
    if undo:
        subs_map = {v: k for k, v in subs_map.items()}

    root_path: Path = Path(root_dir)
    gitignore_spec = load_gitignore(root_dir)
    all_files = list(root_path.glob("**/*"))
    matches = [file for file in all_files if not gitignore_spec.match_file(file)]

    if verbose:
        logger.debug(f"Using .gitignore, found {len(matches)} of {len(all_files)} files to process in {root_path.absolute()}")

    # Rename *inside* files first, otherwise we might rename a file or dir before we rename its contents
    for file_path in matches:
        if file_path.is_file() and any(fnmatch.fnmatch(file_path.name, pattern) for pattern in file_patterns):
            replace_in_file(file_path, subs_map, dry_run, verbose)

    # THEN search again to find and rename all the dirs and filenames
    for file_path in matches:
        new_name: str = replace_str(file_path.name, subs_map)
        if new_name != file_path.name:
            logger.info(f"File: {file_path}: Rename: {file_path.name} -> {new_name}")
            if not dry_run:
                file_path.rename(file_path.with_name(new_name))

    if untrack:
        shutil.rmtree(Path(".git"), onexc=delete_protected_file)

        subprocess.run(["git", "init"])

        # also clean up template files
        Path("src/initialise_project.py").unlink()
        Path("tests/test_template.py").unlink()

        logger.info("Cleaned up and re-initialized git.")


def delete_protected_file(func: Callable[..., Any], path: str, exc: BaseException) -> None:
    """Deletes protected files by setting their attributes to writable first."""
    if os.path.exists(path):
        os.chmod(path, S_IWRITE)
        os.unlink(path)


if __name__ == "__main__":
    get_substitution_map()
    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        description="Perform string substitution based on a dictionary of find-replace regexes."
    )
    parser.add_argument(
        "--execute",
        default=False,
        action="store_true",
        help="Execute the find-replace operation. For safety, this is disabled unless explicitly set.",
    )
    parser.add_argument("--undo", action="store_true", help="Reverse the find-replace map and re-run to attempt to undo changes.")
    parser.add_argument("--untrack", action="store_true", help="delete .git directory and re-initialize git")
    parser.add_argument("--verbose", action="store_true", help="Show verbose output.")
    args: argparse.Namespace = parser.parse_args()

    if not args.execute:
        logger.warning("Running in dry-run mode. No changes will be made.")

    root_directory: str = "."  # Set the root directory to start the search
    substitutions = get_substitution_map()
    process_directory(root_directory, substitutions, file_patterns, not args.execute, args.undo, args.untrack, args.verbose)

    logger.success("Template is now ready for use. Happy coding!")
