#!/usr/bin/env python3
"""Move files with git mv and update every reference to them.

Usage (from the repository root):
    tools/reorg/move_files.py <mapping-file>
    tools/reorg/move_files.py <old> <new> [<old> <new> ...]

A mapping file has one "<old path> <new path>" pair per line (paths relative
to the repository root, '#' starts a comment). A directory can be given as a
pair; every tracked file below it is moved.

After moving, references are updated so the tree keeps building:

* #include directives written pcms/-rooted, or (quoted only) relative to
  the including file, are rewritten to the new pcms/-rooted path when the
  moved file lives under include/ or src/ (relative includes are kept
  relative when the moved file is outside both);
* paths in CMakeLists.txt / *.cmake files that name a moved file (relative to
  that CMake file's directory, the repository root, or src/) are rewritten
  relative to the same base.

Only paths change; file contents are otherwise untouched. CMake entries keep
belonging to the same target; moving them between targets is a separate,
hand-written step. Running the script again with the same mapping is a
no-op.
"""

import os
import re
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SRC = os.path.join(ROOT, "src")
INCLUDE_DIR = os.path.join(ROOT, "include")
INCLUDE_ROOTS = (INCLUDE_DIR, SRC)  # pcms/... resolves against these
CODE = re.compile(r"\.(h|hpp|c|cpp|cc|i)(\.in)?$")
INCLUDE = re.compile(r'^(\s*#\s*include\s*)([<"])([^>"]+)([>"])', re.M)
CMAKE_TOKEN = re.compile(r'(?<![\w./${}-])([\w./${}-]+\.(?:hpp|h|cpp|c|hpp\.in|h\.in))(?![\w.])')
CMAKE_PREFIX = re.compile(r"^\$\{(CMAKE_CURRENT_SOURCE_DIR|CMAKE_CURRENT_LIST_DIR|"
                          r"CMAKE_SOURCE_DIR|PROJECT_SOURCE_DIR)\}/")


def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, check=True,
                          capture_output=True, text=True).stdout


def tracked():
    return git("ls-files").split()


def read_pairs(argv):
    if len(argv) == 1:
        pairs = []
        with open(argv[0]) as fh:
            for line in fh:
                line = line.split("#", 1)[0].split()
                if line:
                    pairs.append((line[0], line[1]))
        return pairs
    if len(argv) % 2:
        sys.exit("odd number of paths")
    return list(zip(argv[::2], argv[1::2]))


def expand(pairs, files):
    moves = {}
    for old, new in pairs:
        old, new = old.rstrip("/"), new.rstrip("/")
        if old in files:
            moves[old] = new
            continue
        below = [f for f in files if f.startswith(old + "/")]
        if not below:
            if new in files or any(f.startswith(new + "/") for f in files):
                continue  # already moved
            sys.exit(f"nothing tracked at {old}")
        for f in below:
            moves[f] = new + f[len(old):]
    return moves


def include_rel(path):
    """pcms/... include path of a file below include/ or src/, else None."""
    for root in INCLUDE_ROOTS:
        if path.startswith(root + os.sep):
            return os.path.relpath(path, root)
    return None


def main(argv):
    files = set(tracked())
    moves = expand(read_pairs(argv), files)
    moves = {o: n for o, n in moves.items() if o != n}
    if not moves:
        print("nothing to move")
        return
    for old, new in sorted(moves.items()):
        os.makedirs(os.path.join(ROOT, os.path.dirname(new)), exist_ok=True)
        git("mv", old, new)
    absmoves = {os.path.join(ROOT, o): os.path.join(ROOT, n)
                for o, n in moves.items()}
    inv = {v: k for k, v in absmoves.items()}

    def target_of(cands):
        for c in cands:
            if c in absmoves:
                return absmoves[c]
            if c + ".in" in absmoves:  # generated header
                return absmoves[c + ".in"][:-3]
        return None

    changed = 0
    for f in tracked():
        absf = os.path.join(ROOT, f)
        where = inv.get(absf, absf)  # resolve relative paths from old place
        name = os.path.basename(f)
        is_cmake = name == "CMakeLists.txt" or f.endswith((".cmake", ".cmake.in"))
        if not (CODE.search(f) or is_cmake):
            continue
        with open(absf, encoding="utf-8") as fh:
            text = fh.read()
        if is_cmake:
            base = os.path.dirname(where)
            newbase = os.path.dirname(absf)

            def repl(m):
                tok = m.group(1)
                pm = CMAKE_PREFIX.match(tok)
                rel = tok[pm.end():] if pm else tok
                if "${" in rel:
                    return tok
                for b in (base, ROOT, SRC, INCLUDE_DIR):
                    t = target_of([os.path.normpath(os.path.join(b, rel))])
                    if t:
                        nb = newbase if b == base else b
                        out = os.path.relpath(t, nb)
                        return (pm.group(0) if pm else "") + out
                return tok
            new = CMAKE_TOKEN.sub(repl, text)
        else:
            def repl(m):
                inc = m.group(3)
                cands = [os.path.normpath(os.path.join(r, inc))
                         for r in INCLUDE_ROOTS]
                if m.group(2) == '"':
                    # only quoted includes are looked up next to the including
                    # file; <...> names a library header (e.g. <mfem.hpp>)
                    cands.insert(0, os.path.normpath(
                        os.path.join(os.path.dirname(where), inc)))
                t = target_of(cands)
                if not t:
                    return m.group(0)
                out = include_rel(t)
                if out is None:
                    out = os.path.relpath(t, os.path.dirname(absf))
                return f"{m.group(1)}{m.group(2)}{out}{m.group(4)}"
            new = INCLUDE.sub(repl, text)
        if new != text:
            with open(absf, "w", encoding="utf-8") as fh:
                fh.write(new)
            changed += 1
    print(f"moved {len(moves)} files, updated {changed} files")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
