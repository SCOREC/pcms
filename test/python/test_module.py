"""Test module-level attributes of the pcms Python module."""

import re
import shutil
import subprocess
from pathlib import Path

import pytest

import pcms

REPO_ROOT = Path(__file__).resolve().parents[2]

# X.Y.Z at a release tag, X.Y.Z.devN+g<hash>[.dYYYYMMDD] otherwise
VERSION = re.compile(
    r"(?P<release>\d+\.\d+\.\d+)"
    r"(\.dev(?P<distance>\d+)\+g(?P<hash>[0-9a-f]+)(\.d\d{8})?)?")


def _git(*args):
    return subprocess.run(["git", *args], cwd=REPO_ROOT, capture_output=True,
                          text=True, check=True).stdout.strip()


def _project_version():
    cmake = (REPO_ROOT / "CMakeLists.txt").read_text()
    project = re.search(r"project\(\s*pcms\b.*?VERSION\s+(\S+)", cmake,
                        re.S | re.I)
    assert project is not None
    return project.group(1)


def test_version_is_pep440():
    assert VERSION.fullmatch(pcms.__version__), pcms.__version__


def test_git_revision_is_full_hash_or_unknown():
    revision = pcms.__git_revision__
    assert revision == "unknown" or re.fullmatch(r"[0-9a-f]{40}", revision)


def test_dev_version_names_the_git_revision():
    match = VERSION.fullmatch(pcms.__version__)
    if match["hash"] is None:
        pytest.skip("release build")
    assert pcms.__git_revision__.startswith(match["hash"])


def test_unknown_revision_uses_project_version():
    if pcms.__git_revision__ != "unknown":
        pytest.skip("built from a git checkout")
    assert pcms.__version__ == _project_version()


@pytest.mark.skipif(shutil.which("git") is None
                    or not (REPO_ROOT / ".git").exists(),
                    reason="needs the git checkout pcms was built from")
def test_version_matches_git_describe():
    if pcms.__git_revision__ != _git("rev-parse", "HEAD"):
        pytest.skip("pcms was built from a different commit")
    try:
        describe = _git("describe", "--tags", "--long", "--match", "v[0-9]*")
    except subprocess.CalledProcessError:
        # No release tag is reachable (e.g. a shallow clone), so the version
        # counts from the project() version with an unknown distance.
        tag, distance = "v" + _project_version(), None
    else:
        tag, distance, _ = describe.rsplit("-", 2)
    major, minor, patch = (int(x) for x in tag[1:].split("."))
    match = VERSION.fullmatch(pcms.__version__)
    if distance == "0" and match["hash"] is None:
        assert match["release"] == f"{major}.{minor}.{patch}"
    else:
        assert match["release"] == f"{major}.{minor}.{patch + 1}"
        assert match["distance"] == (distance or "0")
