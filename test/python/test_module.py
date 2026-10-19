"""Test module-level attributes of the pcms Python module."""

import re
from pathlib import Path

import pcms

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_version_matches_project_version():
    cmake = (REPO_ROOT / "CMakeLists.txt").read_text()
    project = re.search(r"project\(\s*pcms\b.*?VERSION\s+(\S+)", cmake,
                        re.S | re.I)
    assert project is not None
    assert pcms.__version__ == project.group(1)
