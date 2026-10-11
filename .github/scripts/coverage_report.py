#!/usr/bin/env python3
"""Parse Cobertura XML coverage reports and render a Markdown summary.

Usage:
    coverage_report.py [LABEL=]PATH [[LABEL=]PATH ...]

Each positional argument is the path to a Cobertura XML file, optionally
prefixed with a label (e.g. ``Python=coverage.xml``). Missing files are
skipped. The report is written to stdout and starts with a marker comment so
that a bot can find and update a single PR comment.
"""

from __future__ import annotations

import argparse
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

MARKER = "<!-- pcms-coverage-report -->"
MAX_FILES_PER_SECTION = 200


def _to_float(value):
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _to_int(value):
    number = _to_float(value)
    return int(number) if number is not None else None


def _fmt_counts(covered, valid):
    """Format a percentage with absolute counts, or an em dash if unknown."""
    if valid is None or valid <= 0:
        return "\u2014"
    covered = covered or 0
    return f"{100.0 * covered / valid:.1f}% ({covered:,}/{valid:,})"


def _fmt_rate(rate, covered, valid):
    """Format a 0..1 rate, preferring absolute counts when available."""
    if rate is None:
        return "\u2014"
    if valid is not None and valid > 0 and covered is not None:
        return _fmt_counts(covered, valid)
    return f"{100.0 * rate:.1f}%"


def _parse(path):
    root = ET.parse(path).getroot()

    def counts(prefix):
        return (
            _to_int(root.get(f"{prefix}-covered")),
            _to_int(root.get(f"{prefix}-valid")),
        )

    lines_covered, lines_valid = counts("lines")
    branches_covered, branches_valid = counts("branches")

    classes = []
    for cls in root.iter("class"):
        filename = cls.get("filename") or cls.get("name") or "?"
        classes.append(
            {
                "filename": filename,
                "line_rate": _to_float(cls.get("line-rate")),
                "branch_rate": _to_float(cls.get("branch-rate")),
                "line_covered": _to_int(cls.get("lines-covered")),
                "line_valid": _to_int(cls.get("lines-valid")),
                "branch_covered": _to_int(cls.get("branches-covered")),
                "branch_valid": _to_int(cls.get("branches-valid")),
            }
        )

    return {
        "lines_covered": lines_covered,
        "lines_valid": lines_valid,
        "branches_covered": branches_covered,
        "branches_valid": branches_valid,
        "has_branch": branches_valid is not None and branches_valid > 0,
        "classes": classes,
    }


def _label_for(path):
    name = path.name.lower()
    if "cpp" in name or "c++" in name:
        return "C++"
    if name.startswith("coverage"):
        return "Python"
    return path.stem


def main(argv=None):
    parser = argparse.ArgumentParser(description="Render Cobertura XML as Markdown.")
    parser.add_argument("inputs", nargs="*", help="[LABEL=]PATH to a Cobertura XML file")
    args = parser.parse_args(argv)

    reports = []
    for spec in args.inputs:
        if "=" in spec:
            label, raw_path = spec.split("=", 1)
        else:
            label, raw_path = None, spec
        path = Path(raw_path)
        if not path.is_file():
            print(f"warning: skipping missing coverage file: {path}", file=sys.stderr)
            continue
        label = label or _label_for(path)
        reports.append((label, _parse(path)))

    lines = [MARKER, "", "## Code Coverage Report", ""]

    if not reports:
        lines.append("_No coverage data was found for this run._")
        print("\n".join(lines))
        return 0

    lines.append("| Source | Line Coverage | Branch Coverage |")
    lines.append("|---|---|---|")
    for label, report in reports:
        lines.append(
            f"| {label} | {_fmt_counts(report['lines_covered'], report['lines_valid'])} "
            f"| {_fmt_counts(report['branches_covered'], report['branches_valid'])} |"
        )
    lines.append("")

    for label, report in reports:
        classes = sorted(
            report["classes"],
            key=lambda c: (c["line_rate"] if c["line_rate"] is not None else 1.0, c["filename"]),
        )
        if not classes:
            continue
        lines.append("<details>")
        lines.append(f"<summary>{label}: {len(classes)} files</summary>")
        lines.append("")
        lines.append("| File | Line Coverage | Branch Coverage |")
        lines.append("|---|---|---|")
        for cls in classes[:MAX_FILES_PER_SECTION]:
            branch_cell = (
                _fmt_rate(cls["branch_rate"], cls["branch_covered"], cls["branch_valid"])
                if report["has_branch"]
                else "\u2014"
            )
            lines.append(
                f"| `{cls['filename']}` | {_fmt_rate(cls['line_rate'], cls['line_covered'], cls['line_valid'])} "
                f"| {branch_cell} |"
            )
        if len(classes) > MAX_FILES_PER_SECTION:
            lines.append(f"| _... {len(classes) - MAX_FILES_PER_SECTION} more files omitted_ | | |")
        lines.append("")
        lines.append("</details>")
        lines.append("")

    print("\n".join(lines).rstrip())
    return 0


if __name__ == "__main__":
    sys.exit(main())
