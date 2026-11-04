#!/usr/bin/env python3
"""Move a sub-library's internal code into namespace pcms::<lib>.

Usage (from the repository root):
    tools/reorg/sublib_namespace.py <lib> <symbols-file> <files-file>

<files-file> lists the files whose code moves (the private headers of the
sub-library and the sources that implement them; '#' starts a comment);
<symbols-file> lists the namespace-scope names they declare that are used
elsewhere. Same rules as backend_namespace.py: `namespace pcms` becomes
`namespace pcms::<lib>` in the listed files, references to pcms::detail
names from a file set that declares its own nested detail namespace are
written pcms::detail::X, and uses of the listed symbols in every other file
become pcms::<lib>::X. Public API types are not listed and stay in pcms.
External library namespaces are left as they are (no sub-library namespace
shadows one). Running the script twice is a no-op.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import backend_namespace as bn  # noqa: E402

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    lib, symfile, filesfile = sys.argv[1:]
    bn.move_namespace(f"pcms::{lib}", bn.read_list(filesfile),
                      bn.read_list(symfile), external_prefix=False)
