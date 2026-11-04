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


def main(name, symfile):
    syms = []
    with open(symfile) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if line:
                syms.append(line)
    ns = f"pcms::backend::{name}"
    bdir = f"src/pcms/backends/{name}/"
    files = [f for f in tracked() if CODE.search(f)]
    own = [f for f in files if f.startswith(bdir)]
    if not own:
        sys.exit(f"no files below {bdir}")
    # names declared directly at the scope of a backend's own `namespace
    # detail` (not names merely used inside it)
    own_detail = set()
    decl = re.compile(r"^\s*(?:template\s*<.*>\s*)?(?:struct|class|using|enum class|enum)\s+(\w+)"
                      r"|^\s*(?:(?:inline|static|constexpr|KOKKOS_\w+|\[\[\w+\]\])\s+)*"
                      r"[\w:<>,]+[\s&*]+(\w+)\s*\(")
    for f in own:
        with open(os.path.join(ROOT, f), encoding="utf-8") as fh:
            depth = None  # brace depth inside a detail namespace
            for line in fh:
                code = line.split("//", 1)[0]
                if depth is None:
                    if re.match(r"\s*namespace detail\b", code):
                        depth = 0
                    else:
                        continue
                elif depth == 1:
                    m = decl.match(code)
                    if m:
                        own_detail.add(m.group(1) or m.group(2))
                    elif re.match(r"\S", code) and "(" in code:
                        # unindented declaration with a complex return type:
                        # the name is the identifier before the first '('
                        f = re.search(r"(\w+)\s*\(", code)
                        if f and f.group(1) not in ("if", "for", "while",
                                                    "switch", "return"):
                            own_detail.add(f.group(1))
                depth += code.count("{") - code.count("}")
                if depth <= 0 and "}" in code:
                    depth = None
    if os.environ.get("REORG_DEBUG"):
        print("own detail names:", sorted(own_detail))
    detail_ref = re.compile(r"(?<![\w:])detail::(\w+)")
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

            def own_fn(code):
                code = qualified.sub(lambda m: f"{ns}::{m.group(1)}", code)
                if own_detail:
                    code = detail_ref.sub(
                        lambda m: m.group(0) if m.group(1) in own_detail
                        else f"pcms::detail::{m.group(1)}", code)
                return ext.sub(lambda m: f"::{m.group(1)}::", code)
            new = transform_code(new, own_fn)
        else:
            def other_fn(code):
                code = qualified.sub(lambda m: f"{ns}::{m.group(1)}", code)
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
    print(f"updated {changed} files ({len(own)} backend files)")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
