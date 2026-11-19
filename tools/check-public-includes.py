#!/usr/bin/env python3
"""Check the includes and the detail namespaces of PCMS's public headers.

Usage: tools/check-public-includes.py [<repository root>]

Every header below include/pcms/ may include the C++ standard library,
Kokkos, MPI and other PCMS public headers. Headers below
include/pcms/backends/<name>/ may in addition include the libraries of that
backend. Headers below include/pcms/coupler/ may include redev (the coupler's
transport). Any other include of a backend
library (Omega_h, MFEM, MeshFields, PETSc, redev, ADIOS2) is reported, and the
script exits with status 1. So is an include of a PCMS header that is not
public (not below include/pcms/ and not generated), e.g. the private
pcms/discretization/discretization.hpp.

Installed code that is not public lives in headers below a detail/ folder, in
a detail namespace: a header outside detail/ may not open a detail namespace,
and a header below detail/ may declare nothing outside one.
"""

import os
import re
import sys

BACKEND_LIBS = {
    "Omega_h": re.compile(r"^(Omega_h|r3d)"),
    "MFEM": re.compile(r"^mfem"),
    "MeshFields": re.compile(r"^(MeshField|meshfields)"),
    "PETSc": re.compile(r"^petsc"),
    "redev": re.compile(r"^redev"),
    "ADIOS2": re.compile(r"^adios2"),
}

# which backend libraries a directory below include/pcms may use
ALLOWED = {
    "backends/omega_h/": {"Omega_h"},
    "backends/mfem/": {"MFEM"},
    # the XGC reverse classification can be built from an Omega_h mesh
    "backends/xgc/": {"Omega_h"},
    "coupler/": {"redev"},
}

# generated into the build tree and installed next to the public headers
GENERATED = {"pcms/configuration.hpp", "pcms/version.hpp"}

INCLUDE = re.compile(r'^\s*#\s*include\s*[<"]([^>"]+)[>"]', re.M)

# comments, string and character literals, and preprocessor lines (with their
# continuations), which declare nothing
NOISE = re.compile(r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"|'
                   r"'(?:\\.|[^'\\])*'|^[ \t]*#(?:[^\n]*\\\n)*[^\n]*",
                   re.S | re.M)
NAMESPACE = re.compile(r'namespace\s+([\w:]*)\s*\{')


def detail_violations(rel, text):
    """The detail-namespace rule for one header (see the module comment)."""
    code = NOISE.sub(' ', text)
    in_detail_dir = '/detail/' in '/' + rel
    if not in_detail_dir:
        if re.search(r'namespace\s+(\w+::)*detail\b', code):
            return [f"include/pcms/{rel}: opens a detail namespace outside a "
                    "detail/ folder"]
        return []
    # the stack of open braces: the namespace path they open, or None for any
    # other brace
    stack = []
    i = 0
    while i < len(code):
        m = NAMESPACE.match(code, i)
        if m:
            stack.append(m.group(1).split('::'))
            i = m.end()
            continue
        c = code[i]
        if c == '{':
            stack.append(None)
        elif c == '}':
            if stack:
                stack.pop()
        elif not c.isspace() and c != ';':
            inside = any(s is None for s in stack) or any(
                'detail' in s for s in stack if s is not None)
            if not inside:
                return [f"include/pcms/{rel}: declares something outside a "
                        "detail namespace"]
        i += 1
    return []


def check(root):
    base = os.path.join(root, "include", "pcms")
    bad = []
    for dirpath, _, files in os.walk(base):
        for name in files:
            if not name.endswith((".hpp", ".h")):
                continue
            path = os.path.join(dirpath, name)
            rel = os.path.relpath(path, base).replace(os.sep, "/")
            allowed = set()
            for prefix, libs in ALLOWED.items():
                if rel.startswith(prefix):
                    allowed |= libs
            with open(path, encoding="utf-8") as fh:
                text = fh.read()
            bad += detail_violations(rel, text)
            for inc in INCLUDE.findall(text):
                if inc.startswith("pcms/"):
                    # only public (installed or generated) PCMS headers
                    if (not os.path.exists(os.path.join(base, inc[len("pcms/"):]))
                            and inc not in GENERATED):
                        bad.append(f"include/pcms/{rel}: includes <{inc}> "
                                   "(not a public header)")
                    continue
                for lib, pattern in BACKEND_LIBS.items():
                    if pattern.match(inc) and lib not in allowed:
                        bad.append(f"include/pcms/{rel}: includes <{inc}> ({lib})")
    return bad


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(__file__), "..")
    bad = check(os.path.abspath(root))
    for line in bad:
        print(line)
    if bad:
        print(f"{len(bad)} problem(s) in public headers")
        return 1
    print("public headers: no disallowed includes; detail namespaces only "
          "in detail/ headers")
    return 0


if __name__ == "__main__":
    sys.exit(main())
