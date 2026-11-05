#!/usr/bin/env python3
"""Give every header a path-derived include guard.

Usage (from the repository root):
    tools/reorg/normalize_guards.py

For each tracked header (*.h, *.hpp, *.h.in, *.hpp.in) whose first directives
are an `#ifndef X` / `#define X` guard closed by the file's last `#endif`, the
guard macro is replaced by a name derived from the header's include path:

  include/pcms/field/field.hpp           -> PCMS_FIELD_FIELD_HPP
  src/pcms/backends/omega_h/adj_search.hpp -> PCMS_BACKENDS_OMEGA_H_ADJ_SEARCH_HPP
  include/pcms/capi/client.h             -> PCMS_CAPI_CLIENT_H
  src/pcms/configuration.hpp.in          -> PCMS_CONFIGURATION_HPP
  test/field_test_utils.hpp              -> PCMS_TEST_FIELD_TEST_UTILS_HPP

The trailing comment of the closing #endif is updated as well (and added if
missing). Guards are never replaced by #pragma once. Running the script
twice is a no-op.
"""

import re
import subprocess
import sys


def guard_name(path):
    p = re.sub(r"\.in$", "", path)
    if p.startswith(("include/pcms/", "src/pcms/")):
        rel = p.split("/pcms/", 1)[1]
        prefix = "PCMS_"
    elif p.startswith("test/"):
        rel = p[len("test/"):]
        prefix = "PCMS_TEST_"
    else:
        rel = p
        prefix = "PCMS_"
    return prefix + re.sub(r"[^A-Za-z0-9]", "_", rel).upper()


def main():
    files = subprocess.run(["git", "ls-files"], check=True, capture_output=True,
                           text=True).stdout.split()
    headers = [f for f in files if re.search(r"\.(h|hpp)(\.in)?$", f)
               and not f.startswith(("tools/", "src/pcms/external/"))]
    names, changed, skipped = {}, 0, []
    for f in headers:
        new_guard = guard_name(f)
        if new_guard in names:
            sys.exit(f"guard clash: {f} and {names[new_guard]}")
        names[new_guard] = f
        text = open(f, encoding="utf-8").read()
        m = re.search(r"^([ \t]*#[ \t]*ifndef[ \t]+)(\w+)([ \t]*\n[ \t]*#[ \t]*define[ \t]+)\2\b",
                      text, re.M)
        if not m:
            skipped.append(f)
            continue
        old = m.group(2)
        endifs = list(re.finditer(r"^[ \t]*#[ \t]*endif\b[^\n]*", text, re.M))
        if not endifs or text[endifs[-1].end():].strip():
            skipped.append(f)
            continue
        new = (text[:m.start()] + m.group(1) + new_guard + m.group(3) + new_guard
               + text[m.end():])
        # recompute the last #endif in the edited text
        last = list(re.finditer(r"^[ \t]*#[ \t]*endif\b[^\n]*", new, re.M))[-1]
        new = new[:last.start()] + f"#endif // {new_guard}" + new[last.end():]
        if new != text:
            open(f, "w", encoding="utf-8").write(new)
            changed += 1
        del old
    print(f"updated {changed} of {len(headers)} headers")
    for f in skipped:
        print(f"SKIPPED (no standard guard): {f}", file=sys.stderr)
    return 1 if skipped else 0


if __name__ == "__main__":
    sys.exit(main())
