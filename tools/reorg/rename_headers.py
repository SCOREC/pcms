#!/usr/bin/env python3
"""Rename C++ headers from .h to .hpp and update every reference.

Usage (from the repository root):
    tools/reorg/rename_headers.py <dir> [<dir> ...]

Each <dir> is a path relative to the repository root, e.g. src/pcms/field.
Every ``*.h`` and ``*.h.in`` file below it is renamed with ``git mv`` to
``*.hpp`` / ``*.hpp.in``. Then every reference to a renamed header is
rewritten:

* ``#include`` / ``%include`` directives in C/C++/SWIG sources, whether the
  include is written relative to the including file or relative to ``src/``
  (``pcms/...``);
* file names listed in ``CMakeLists.txt`` / ``*.cmake`` files, relative to the
  directory of that CMake file, its binary-dir equivalent, or ``src/``.

The C API headers in ``src/pcms/capi`` are never renamed. The script is
deterministic: running it twice is a no-op.
"""

import os
import re
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SRC = os.path.join(ROOT, "src")
EXCLUDE_DIRS = [os.path.join(SRC, "pcms", "capi")]
CODE_EXT = (".h", ".hpp", ".c", ".cpp", ".cc", ".i", ".h.in", ".hpp.in")
CMAKE_NAMES = ("CMakeLists.txt",)
CMAKE_EXT = (".cmake", ".cmake.in")
SKIP_TREES = (".git", "build")

INCLUDE_RE = re.compile(r'^(\s*[#%]\s*include\s*)([<"])([^>"]+)([>"])', re.M)
CMAKE_TOKEN_RE = re.compile(r'(?<![\w./${}-])([\w./${}-]+\.h(?:\.in)?)(?![\w.])')


def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, check=True,
                          capture_output=True, text=True).stdout


def tracked_files():
    return [f for f in git("ls-files").splitlines()
            if not f.startswith(SKIP_TREES)]


def new_name(path):
    if path.endswith(".h.in"):
        return path[: -len(".h.in")] + ".hpp.in"
    assert path.endswith(".h")
    return path + "pp"


def collect(dirs, files):
    renames = {}
    for d in dirs:
        absd = os.path.abspath(os.path.join(ROOT, d))
        for f in files:
            absf = os.path.join(ROOT, f)
            if not absf.startswith(absd + os.sep):
                continue
            if any(absf.startswith(e + os.sep) for e in EXCLUDE_DIRS):
                continue
            if f.endswith((".h", ".h.in")):
                renames[absf] = os.path.join(ROOT, new_name(f))
    return renames


def rewrite_includes(path, text, renamed):
    def repl(m):
        inc = m.group(3)
        if not inc.endswith(".h"):
            return m.group(0)
        cands = [os.path.normpath(os.path.join(os.path.dirname(path), inc)),
                 os.path.normpath(os.path.join(SRC, inc))]
        # generated headers are included as "pcms/configuration.h" but live
        # in src/ as *.h.in
        if any(c in renamed or c + ".in" in renamed for c in cands):
            return f"{m.group(1)}{m.group(2)}{inc}pp{m.group(4)}"
        return m.group(0)
    return INCLUDE_RE.sub(repl, text)


def rewrite_cmake(path, text, renamed):
    base = os.path.dirname(path)

    def repl(m):
        tok = m.group(1)
        rel = re.sub(r"^\$\{(CMAKE_CURRENT_SOURCE_DIR|CMAKE_CURRENT_LIST_DIR|"
                     r"CMAKE_CURRENT_BINARY_DIR|CMAKE_SOURCE_DIR|"
                     r"PROJECT_SOURCE_DIR)\}/", "", tok)
        if "${" in rel:
            return tok
        cands = [os.path.normpath(os.path.join(base, rel)),
                 os.path.normpath(os.path.join(ROOT, rel)),
                 os.path.normpath(os.path.join(SRC, rel))]
        for c in cands:
            if c in renamed or c + ".in" in renamed:
                return new_name(tok)
        return tok
    return CMAKE_TOKEN_RE.sub(repl, text)


def main(dirs):
    files = tracked_files()
    renamed = collect(dirs, files)
    if not renamed:
        print("nothing to rename")
        return
    for old, new in sorted(renamed.items()):
        git("mv", os.path.relpath(old, ROOT), os.path.relpath(new, ROOT))
    # references are resolved against the *old* paths
    files = tracked_files()
    inv = {v: k for k, v in renamed.items()}
    changed = 0
    for f in files:
        absf = os.path.join(ROOT, f)
        name = os.path.basename(f)
        is_code = f.endswith(CODE_EXT)
        is_cmake = name in CMAKE_NAMES or f.endswith(CMAKE_EXT)
        if not (is_code or is_cmake):
            continue
        # a renamed file resolves its relative includes from its old location
        where = inv.get(absf, absf)
        with open(absf, encoding="utf-8") as fh:
            text = fh.read()
        new = rewrite_includes(where, text, renamed) if is_code else \
            rewrite_cmake(where, text, renamed)
        if new != text:
            with open(absf, "w", encoding="utf-8") as fh:
                fh.write(new)
            changed += 1
    print(f"renamed {len(renamed)} headers, updated {changed} files")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
