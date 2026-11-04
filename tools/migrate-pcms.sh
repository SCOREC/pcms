#!/bin/bash
# Migrate a downstream source tree from PCMS 0.3 to PCMS 0.4.
#
# Usage: tools/migrate-pcms.sh [--dry-run] [--in-tree] <directory or file>...
#
# Applies tools/reorg/header-map.txt and tools/reorg/symbol-map.txt:
#   * rewrites #include directives of renamed public headers (C, C++, Fortran
#     preprocessor files);
#   * renames moved public C++ names, renamed macros, CMake options and the
#     pcms::core target;
#   * migrates Python code with tools/reorg/python_api_names.py;
#   * reports, with file and line, every use of a header or name that is no
#     longer public or whose API changed shape ("manual" entries), with the
#     replacement to use. Those need a manual edit; see MIGRATION.md.
# Files are modified in place (back up or commit first). With --dry-run
# nothing is written; the planned changes are listed. --in-tree is for code
# inside PCMS (e.g. a pull request ported across the reorganization): includes
# of private headers and qualified uses of internal names are rewritten to
# their new locations instead of being reported.
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
dry=0
in_tree=0
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) dry=1; shift ;;
    --in-tree) in_tree=1; shift ;;
    *) break ;;
  esac
done
if [ $# -lt 1 ]; then
  sed -n '2,21p' "$0"
  exit 2
fi
for t in "$@"; do
  if [ ! -e "$t" ]; then
    echo "migrate-pcms.sh: $t does not exist" >&2
    exit 2
  fi
done

PCMS_MIGRATE_DRY=$dry PCMS_MIGRATE_IN_TREE=$in_tree python3 - "$here/reorg/header-map.txt" \
  "$here/reorg/symbol-map.txt" "$@" <<'PY'
import os
import re
import sys

header_map, symbol_map = sys.argv[1:3]
roots = sys.argv[3:]
dry = os.environ.get("PCMS_MIGRATE_DRY") == "1"
in_tree = os.environ.get("PCMS_MIGRATE_IN_TREE") == "1"

headers, private_headers = {}, {}
for line in open(header_map):
    line = line.split("#", 1)[0].split()
    if len(line) >= 2:
        if len(line) == 3 and line[2] == "private":
            private_headers[line[0]] = line[1]
        else:
            headers[line[0]] = line[1]

renames, cmake, private_syms, manual = [], [], [], []
for line in open(symbol_map):
    if line.startswith("#") or not line.strip():
        continue
    if line.startswith("manual"):
        old, new = line[len("manual"):].split("->", 1)
        manual.append((old.strip(), new.strip()))
        continue
    kind, old, new = line.split()
    if kind in ("cpp", "macro"):
        renames.append((old, new))
    elif kind in ("cmake-option", "cmake-target"):
        cmake.append((old, new))
    elif kind == "cpp-private":
        private_syms.append((old, new))

if in_tree:
    # inside PCMS, private headers and internal names are still usable
    headers.update({k: v for k, v in private_headers.items() if v != "-"})
    private_headers = {k: v for k, v in private_headers.items() if v == "-"}
    renames += private_syms
    private_syms = []

CODE = re.compile(r"\.(h|hh|hpp|hxx|c|cc|cpp|cxx|cu|F90|f90|F|inc)$")
CMAKE = re.compile(r"(CMakeLists\.txt|\.cmake)$")
INCLUDE = re.compile(r'^(\s*#\s*include\s*)([<"])([^>"]+)([>"])', re.M)


def word(old):
    # old may contain ::, so build the boundaries by hand
    return re.compile(r"(?<![\w:])" + re.escape(old) + r"(?![\w])")


def walk(roots):
    for root in roots:
        if os.path.isfile(root):
            yield root
            continue
        for dirpath, dirnames, files in os.walk(root):
            dirnames[:] = [d for d in dirnames if d not in (".git", "build")]
            for name in files:
                yield os.path.join(dirpath, name)


changed, reports = 0, []
for path in walk(roots):
    name = os.path.basename(path)
    is_code, is_cmake = CODE.search(name), CMAKE.search(name)
    if not (is_code or is_cmake):
        continue
    try:
        text = open(path, encoding="utf-8").read()
    except UnicodeDecodeError:
        continue
    new = text
    if is_code:
        def inc(m):
            path_ = m.group(3)
            if path_ in headers:
                return f"{m.group(1)}{m.group(2)}{headers[path_]}{m.group(4)}"
            return m.group(0)
        new = INCLUDE.sub(inc, new)
        for old, rep in renames:
            new = word(old).sub(rep, new)
    if is_cmake:
        for old, rep in cmake:
            new = word(old).sub(rep, new)
    # reports
    for n, line in enumerate(new.split("\n"), 1):
        m = INCLUDE.match(line)
        if m and m.group(3) in private_headers:
            reports.append(f"{path}:{n}: includes {m.group(3)}, which is no "
                           "longer installed (see MIGRATION.md)")
        for old, rep in private_syms:
            if word(old).search(line):
                reports.append(f"{path}:{n}: uses {old}, no longer public "
                               f"(now internal {rep}); see MIGRATION.md")
        for old, rep in manual:
            key = re.split(r"[ (]", old)[0]
            if key.startswith(("pcms::", "redev::")) and word(key).search(line):
                reports.append(f"{path}:{n}: {old} -> {rep}")
    if new != text:
        changed += 1
        print(("would update " if dry else "updated ") + path)
        if not dry:
            open(path, "w", encoding="utf-8").write(new)

print(f"{'would update' if dry else 'updated'} {changed} C/C++/CMake files")
if reports:
    print(f"\n{len(reports)} place(s) need a manual change:")
    for r in reports:
        print("  " + r)
PY

# Python code
mapfile -t pyfiles < <(find "$@" -name '*.py' -not -path '*/.git/*' \
  -exec grep -lE 'LagrangeFunctionSpace|PolynomialReconstructionFunctionSpace|OmegaHConservativeProjection|OmegaHControlVariateProjection|create_uniform_grid_binary_field|create_uniform_grid_from_mesh|compute_entity_centroids|map_entity_field_to_vertices_|\.mesh\(\)' {} +)
if [ ${#pyfiles[@]} -gt 0 ]; then
  if [ $dry = 1 ]; then
    printf 'would update (Python) %s\n' "${pyfiles[@]}"
  else
    python3 "$here/reorg/python_api_names.py" "${pyfiles[@]}"
  fi
fi
