#!/usr/bin/env python3
"""Rewrite relative includes of PCMS headers as pcms/-rooted includes.

Usage (from the repository root):
    tools/reorg/absolute_includes.py <dir> [<dir> ...]

For every C/C++ source or header below each <dir>, each *quoted* include
whose path, resolved relative to the including file's directory, names a
tracked file under src/pcms (or a generated header there) is rewritten to
"pcms/<path>", keeping the quote style. Includes that are already
pcms/-rooted, angle-bracket includes, and SWIG %include directives are left
unchanged (the SWIG interface files are processed by `swig -fortran` from
their own directory). Running it twice is a no-op.
"""

import os
import re
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SRC_PCMS = os.path.join(ROOT, "src", "pcms")
CODE = re.compile(r"\.(h|hpp|c|cpp|cc)(\.in)?$")
INCLUDE = re.compile(r'^(\s*#\s*include\s*)"([^"]+)"', re.M)


def tracked():
    out = subprocess.run(["git", "ls-files"], cwd=ROOT, check=True,
                         capture_output=True, text=True).stdout.split()
    return out


def main(dirs):
    files = tracked()
    exist = {os.path.join(ROOT, f) for f in files}
    # generated headers live as *.in in the source tree
    exist |= {p[:-3] for p in exist if p.endswith(".hpp.in")}
    changed = rewrites = 0
    for f in files:
        absf = os.path.join(ROOT, f)
        if not CODE.search(f):
            continue
        if not any(absf.startswith(os.path.abspath(d) + os.sep) for d in dirs):
            continue

        def repl(m):
            nonlocal rewrites
            inc = m.group(2)
            if inc.startswith("pcms/"):
                return m.group(0)
            target = os.path.normpath(os.path.join(os.path.dirname(absf), inc))
            if target in exist and target.startswith(SRC_PCMS + os.sep):
                rewrites += 1
                rel = os.path.relpath(target, os.path.dirname(SRC_PCMS))
                return f'{m.group(1)}"{rel}"'
            return m.group(0)

        with open(absf, encoding="utf-8") as fh:
            text = fh.read()
        new = INCLUDE.sub(repl, text)
        if new != text:
            with open(absf, "w", encoding="utf-8") as fh:
                fh.write(new)
            changed += 1
    print(f"rewrote {rewrites} includes in {changed} files")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main([os.path.join(ROOT, d) for d in sys.argv[1:]])
