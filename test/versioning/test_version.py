"""Version generation and package metadata tests without PCMS dependencies."""

from email.parser import Parser
from pathlib import Path
import re
import shutil
import subprocess

import pytest
from scikit_build_core.build import prepare_metadata_for_build_wheel

ROOT = Path(__file__).resolve().parents[2]
PROJECT_VERSION = re.search(
    r"project\(\s*pcms\b.*?VERSION\s+(\S+)",
    (ROOT / "CMakeLists.txt").read_text(), re.S | re.I).group(1)
CMAKE = shutil.which("cmake")


def git(source, *args):
    return subprocess.check_output(
        ["git", "-C", str(source), *args], text=True).strip()


@pytest.fixture
def source(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    shutil.copy(ROOT / "pyproject.toml", source)
    (source / "tracked").write_text("initial\n")
    return source


def init_git(source):
    git(source, "init")
    git(source, "config", "user.email", "test@example.invalid")
    git(source, "config", "user.name", "Version tests")
    git(source, "add", ".")
    git(source, "commit", "-m", "initial")


def header(source, tmp_path, override=""):
    output = tmp_path / "version.hpp"
    subprocess.run([
        CMAKE, f"-DSOURCE_DIR={source}",
        f"-DGIT_EXECUTABLE={shutil.which('git')}",
        f"-DPROJECT_VERSION={PROJECT_VERSION}",
        f"-DVERSION_OVERRIDE={override}",
        f"-DINPUT={ROOT / 'src/pcms/version.hpp.in'}",
        f"-DOUTPUT={output}", "-P", str(ROOT / "cmake/PcmsVersion.cmake"),
    ], check=True)
    values = re.findall(r'return "([^"]+)";', output.read_text())
    return tuple(values)


def metadata(source, tmp_path, monkeypatch):
    monkeypatch.chdir(source)
    output = tmp_path / "metadata"
    output.mkdir(exist_ok=True)
    dist_info = prepare_metadata_for_build_wheel(str(output))
    return Parser().parsestr((output / dist_info / "METADATA").read_text())["Version"]


def test_archive_fallback(source, tmp_path, monkeypatch):
    assert metadata(source, tmp_path, monkeypatch) == PROJECT_VERSION
    assert header(source, tmp_path) == (PROJECT_VERSION, "unknown")


def test_sdist_keeps_development_version(source, tmp_path, monkeypatch):
    preserved = "0.3.1.dev511+g8a9fa572"
    (source / "PKG-INFO").write_text(
        f"Metadata-Version: 2.1\nName: pcms\nVersion: {preserved}\n")
    package_version = metadata(source, tmp_path, monkeypatch)
    assert package_version == preserved
    assert header(source, tmp_path, package_version) == (preserved, "unknown")


@pytest.mark.parametrize("state", ["release", "development", "dirty", "unrelated-tag"])
def test_git_versions(source, tmp_path, monkeypatch, state):
    init_git(source)
    if state == "unrelated-tag":
        git(source, "tag", "-a", "checkpoint", "-m", "unrelated")
    else:
        git(source, "tag", "v0.3.0")
    if state == "development":
        (source / "tracked").write_text("next\n")
        git(source, "commit", "-am", "next")
    elif state in ("dirty", "unrelated-tag"):
        (source / "tracked").write_text("dirty\n")
    generated, revision = header(source, tmp_path)
    assert revision == git(source, "rev-parse", "HEAD")
    if state == "release":
        assert generated == "0.3.0"
    else:
        distance = "1" if state == "development" else "0"
        assert re.fullmatch(
            rf"0\.3\.1\.dev{distance}\+g{revision[:7]}[0-9a-f]*"
            + (r"\.d\d{8}" if state != "development" else ""), generated)
    if state != "unrelated-tag":
        assert metadata(source, tmp_path, monkeypatch) == generated


def test_archive_inside_other_repository(source, tmp_path, monkeypatch):
    init_git(tmp_path)
    git(tmp_path, "tag", "v9.9.9")
    assert metadata(source, tmp_path, monkeypatch) == PROJECT_VERSION
    assert header(source, tmp_path) == (PROJECT_VERSION, "unknown")
