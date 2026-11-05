#!/usr/bin/env python3
"""Move a backend's code into namespace pcms::backend::<name>.

Usage (from the repository root):
    tools/reorg/backend_namespace.py <name> <symbols-file>

<symbols-file> lists, one per line, the namespace-scope names declared by the
backend's files in namespace pcms ('#' starts a comment). The backend's files
are every tracked C++ file below src/pcms/backends/<name>/.

In the backend's own files:
* ``namespace pcms`` (opening line and closing comment) becomes
  ``namespace pcms::backend::<name>``; nested ``pcms::detail`` etc. are left;
* ``pcms::<Symbol>`` becomes ``pcms::backend::<name>::<Symbol>`` (e.g. in
  std::hash specializations);
* a compact ``namespace pcms::detail`` block becomes
  ``namespace pcms::backend::<name>::detail`` (nested under the new
  namespace); uses in other files of the names declared in the moved detail
  namespaces (``pcms::detail::X`` or ``detail::X``) are qualified with the
  new namespace, unless another file also declares X in a detail namespace;
* if the backend declares its own nested ``detail`` namespace, references
  to ``detail::X`` for names it does not declare there are written
  ``pcms::detail::X`` (inside pcms::backend::<name>, ``detail::`` would
  otherwise find only the backend's own detail namespace);
* uses of a library namespace that a backend namespace shadows are written
  with a leading ``::`` (``mfem``, which pcms::backend::mfem hides from
  code inside pcms::backend); other library namespaces (Kokkos, Omega_h,
  ...) are left as written.

In every other tracked C++ file, each use of a listed symbol that is not in
a comment or string literal and is not a member access is rewritten to
``pcms::backend::<name>::<Symbol>``; an existing ``pcms::`` qualifier is
replaced. Running the script twice is a no-op.
"""

import os
import re
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
CODE = re.compile(r"\.(h|hpp|cpp|cc|c)(\.in)?$")
# library namespaces shadowed by a pcms::backend::<name> namespace
EXTERNAL = ("mfem",)


def tracked():
    return subprocess.run(["git", "ls-files"], cwd=ROOT, check=True,
                          capture_output=True, text=True).stdout.split()


def code_spans(text):
    """Yield (start, end) of text that is code (not comment/string/char)."""
    i, n, start = 0, len(text), 0
    while i < n:
        c = text[i]
        two = text[i:i + 2]
        if two == "//" or two == "/*" or c in "\"'":
            if start < i:
                yield start, i
            if two == "//":
                j = text.find("\n", i)
                i = n if j < 0 else j
            elif two == "/*":
                j = text.find("*/", i + 2)
                i = n if j < 0 else j + 2
            elif c == '"' and text[max(0, i - 1)] == "R":  # raw string R"(...)"
                m = re.match(r'"([^(]*)\(', text[i:])
                end = text.find(")" + m.group(1) + '"', i) if m else -1
                i = n if end < 0 else end + len(m.group(1)) + 2
            else:
                j = i + 1
                while j < n and text[j] != c:
                    j += 2 if text[j] == "\\" else 1
                i = j + 1
            start = i
        elif c == "#" and (i == 0 or text[i - 1] == "\n"):
            # preprocessor line: keep #include paths untouched, but code in
            # #if/#define lines is still code
            j = text.find("\n", i)
            j = n if j < 0 else j
            if re.match(r"#\s*include", text[i:j]):
                if start < i:
                    yield start, i
                start = j
            i = j
        else:
            i += 1
    if start < n:
        yield start, n


def transform_code(text, fn):
    out, last = [], 0
    for s, e in code_spans(text):
        out.append(text[last:s])
        out.append(fn(text[s:e]))
        last = e
    out.append(text[last:])
    return "".join(out)


DECL = re.compile(r"^\s*(?:template\s*<.*>\s*)?(?:struct|class|using|enum class|enum)\s+(\w+)"
                  r"|^\s*(?:(?:inline|static|constexpr|[A-Z][A-Z0-9_]+|\[\[\w+\]\])\s+)*"
                  r"[\w:<>,]+[\s&*]+(\w+)\s*\(")


def collect_detail_names(paths):
    """Names declared directly in `namespace detail` blocks (nested in
    namespace pcms...) or compact `namespace pcms::detail` blocks of paths.
    Returns (all names, names of compact blocks)."""
    names, compact_names = set(), set()
    for f in paths:
        try:
            fh = open(os.path.join(ROOT, f), encoding="utf-8")
        except OSError:
            continue
        with fh:
            depth = None  # brace depth inside a detail namespace
            compact = False
            for line in fh:
                code = line.split("//", 1)[0]
                if depth is None:
                    if re.match(r"\s*namespace detail\b", code):
                        depth, compact = 0, False
                    elif re.match(r"namespace pcms::detail\b", code):
                        depth, compact = 0, True
                    else:
                        continue
                elif depth == 1:
                    m = DECL.match(code)
                    name = (m.group(1) or m.group(2)) if m else None
                    if name is None and re.match(r"\S", code) and "(" in code:
                        # unindented declaration with a complex return type:
                        # the name is the identifier before the first '('
                        f2 = re.search(r"(\w+)\s*\(", code)
                        if f2 and f2.group(1) not in ("if", "for", "while",
                                                      "switch", "return"):
                            name = f2.group(1)
                    if name:
                        names.add(name)
                        if compact:
                            compact_names.add(name)
                depth += code.count("{") - code.count("}")
                if depth <= 0 and "}" in code:
                    depth = None
    return names, compact_names


def read_list(path):
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if line:
                out.append(line)
    return out


def main(name, symfile):
    syms = read_list(symfile)
    dirs = (f"src/pcms/backends/{name}/", f"include/pcms/backends/{name}/")
    own = [f for f in tracked() if CODE.search(f) and f.startswith(dirs)]
    if not own:
        sys.exit(f"no files below {dirs[0]}")
    move_namespace(f"pcms::backend::{name}", own, syms, external_prefix=True)


def move_namespace(ns, own, syms, external_prefix):
    """Move the code of the files `own` from namespace pcms into ns and
    qualify uses of `syms` elsewhere (see the module docstring)."""
    files = [f for f in tracked() if CODE.search(f)]
    # names declared directly at the scope of a backend's own `namespace
    # detail` (not names merely used inside it)
    own_detail, compact_detail = collect_detail_names(own)
    # names a file outside `own` declares in a detail namespace: uses of
    # those names cannot be attributed to the moved code
    others = [f for f in files if f not in own]
    elsewhere, _ = collect_detail_names(others)
    moved_detail = own_detail - elsewhere
    if os.environ.get("REORG_DEBUG"):
        print("own detail names:", sorted(own_detail))
    detail_ref = re.compile(r"(?<![\w:])detail::(\w+)")
    compact_ref = re.compile(r"(?<![\w:])(?:::)?(?:pcms::)?detail::("
                             + "|".join(map(re.escape, sorted(moved_detail)))
                             + r")\b") if moved_detail else None
    sym_alt = "|".join(map(re.escape, sorted(syms, key=len, reverse=True)))
    qualified = re.compile(r"(?<![\w:])(?:::)?pcms::(" + sym_alt + r")\b")
    ext = re.compile(r"(?<![\w:])(" + "|".join(EXTERNAL) + r")::")
    changed = 0
    for f in files:
        path = os.path.join(ROOT, f)
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
        if f in own:
            new = re.sub(r"^namespace pcms([ \t]*\{|[ \t]*$)",
                         lambda m: f"namespace {ns}{m.group(1)}", text,
                         flags=re.M)
            new = re.sub(r"(\}[ \t]*//[ \t]*namespace) pcms[ \t]*$",
                         lambda m: f"{m.group(1)} {ns}", new, flags=re.M)
            # compact nested detail namespace: pcms::detail -> ns::detail
            new = re.sub(r"^namespace pcms::detail([ \t]*\{|[ \t]*$)",
                         lambda m: f"namespace {ns}::detail{m.group(1)}", new,
                         flags=re.M)
            new = re.sub(r"(\}[ \t]*//[ \t]*namespace) pcms::detail[ \t]*$",
                         lambda m: f"{m.group(1)} {ns}::detail", new, flags=re.M)

            def own_fn(code):
                code = qualified.sub(lambda m: f"{ns}::{m.group(1)}", code)
                if moved_detail:
                    code = compact_ref.sub(
                        lambda m: f"{ns}::detail::{m.group(1)}"
                        if m.group(0).lstrip(":").startswith("pcms::")
                        else m.group(0), code)
                if own_detail:
                    code = detail_ref.sub(
                        lambda m: m.group(0) if m.group(1) in own_detail
                        else f"pcms::detail::{m.group(1)}", code)
                if external_prefix:
                    code = ext.sub(lambda m: f"::{m.group(1)}::", code)
                return code
            new = transform_code(new, own_fn)
        else:
            def other_fn(code):
                code = qualified.sub(lambda m: f"{ns}::{m.group(1)}", code)
                if moved_detail:
                    code = compact_ref.sub(
                        lambda m: f"{ns}::detail::{m.group(1)}", code)
                # bare uses (not already qualified, not member access)
                code = re.sub(r"(?<![\w:.])(?<!->)(" + sym_alt + r")\b",
                              lambda m: f"{ns}::{m.group(1)}", code)
                return code
            new = transform_code(text, other_fn)
            # undo double qualification produced on already-rewritten text
            new = new.replace(f"{ns}::{ns}::", f"{ns}::")
        if new != text:
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(new)
            changed += 1
    print(f"updated {changed} files ({len(own)} files moved to {ns})")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
