#!/usr/bin/env python3
"""Point callers of the Omega_h-typed helpers at the Omega_h backend.

Usage (from the repository root):
    tools/reorg/omega_h_boundary_calls.py

* pcms::MakeConstRank2View(<omega_h array>, <dim>) (the two-argument overload)
  -> pcms::backend::omega_h::MakeConstRank2View(...)
* pcms::CreateUniformGridFromMesh / unqualified uses in namespace pcms
  -> pcms::backend::omega_h::CreateUniformGridFromMesh
and adds #include "pcms/backends/omega_h/omega_h.hpp" to each changed file
that does not include it yet. Running it twice is a no-op.
"""
import re
import subprocess

files = subprocess.run(["git", "ls-files", "src", "test", "examples", "tools"],
                       check=True, capture_output=True, text=True).stdout.split()
OH = "pcms::backend::omega_h::"
two_arg = re.compile(r"(?<![\w:])pcms::MakeConstRank2View\((\s*[^,()]+(?:\([^()]*\))?\s*),(\s*\d+\s*)\)")
grid = re.compile(r"(?<![\w:])(?:pcms::)?CreateUniformGridFromMesh\b")
include = '#include "pcms/backends/omega_h/omega_h.hpp"\n'
for f in files:
    if not re.search(r"\.(cpp|hpp|h)$", f) or f.startswith("include/") \
            or f == "tools/reorg/omega_h_boundary_calls.py":
        continue
    text = open(f).read()
    new = two_arg.sub(lambda m: f"{OH}MakeConstRank2View({m.group(1)},{m.group(2)})", text)
    if f != "src/pcms/backends/omega_h/uniform_grid_binary_field.hpp":
        new = grid.sub(f"{OH}CreateUniformGridFromMesh", new)
    else:  # already in namespace pcms::backend::omega_h; only needs the header
        pass
    needs = new != text or (f == "src/pcms/backends/omega_h/uniform_grid_binary_field.hpp")
    if needs and include.strip() not in new:
        lines = new.split("\n")
        last = max(i for i, l in enumerate(lines) if l.startswith("#include"))
        lines.insert(last + 1, include.rstrip("\n"))
        new = "\n".join(lines)
    new = new.replace(OH + OH, OH)
    if new != text:
        open(f, "w").write(new)
        print("updated", f)
