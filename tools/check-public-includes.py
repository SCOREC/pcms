#!/usr/bin/env python3
"""Reject backend-library includes in PCMS's public headers.

Usage: tools/check-public-includes.py [<repository root>]

Every header below include/pcms/ may include the C++ standard library,
Kokkos, MPI and other PCMS public headers. Headers below
include/pcms/backends/<name>/ may in addition include the libraries of that
backend. Headers below include/pcms/coupler/ may include redev (the coupler's
transport, see tools/reorg/decisions.md). Any other include of a backend
library (Omega_h, MFEM, MeshFields, PETSc, redev, ADIOS2) is reported, and the
script exits with status 1.
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

INCLUDE = re.compile(r'^\s*#\s*include\s*[<"]([^>"]+)[>"]', re.M)


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
            for inc in INCLUDE.findall(text):
                if inc.startswith("pcms/"):
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
        print(f"{len(bad)} backend-library include(s) in public headers")
        return 1
    print("public headers: no disallowed backend-library includes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
