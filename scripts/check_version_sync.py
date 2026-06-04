#!/usr/bin/env python3
"""Assert the project version is identical across CMakeLists.txt, python/pyproject.toml,
and recipes/bioconda/meta.yaml.

A version bump must move all three together — otherwise the binary, the PyPI package, and
the bioconda recipe drift apart (and a stale recipe ships the wrong source). Wired into
CI as a fast lint. Exit non-zero on mismatch / unparseable version.
"""
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
SOURCES = {
    "CMakeLists.txt": r"project\([^)]*VERSION\s+([0-9]+\.[0-9]+\.[0-9]+)",
    "python/pyproject.toml": r'(?m)^\s*version\s*=\s*"([0-9]+\.[0-9]+\.[0-9]+)"',
    "recipes/bioconda/meta.yaml": r'set\s+version\s*=\s*"([0-9]+\.[0-9]+\.[0-9]+)"',
}


def main():
    versions = {}
    for path, pat in SOURCES.items():
        text = (ROOT / path).read_text()
        m = re.search(pat, text)
        versions[path] = m.group(1) if m else None
    for path, v in versions.items():
        print(f"  {v or '??':8s}  {path}")
    missing = [p for p, v in versions.items() if v is None]
    if missing:
        print(f"ERROR: could not parse a version from: {missing}", file=sys.stderr)
        return 1
    uniq = set(versions.values())
    if len(uniq) != 1:
        print(f"ERROR: project version is NOT in sync across sources: {versions}", file=sys.stderr)
        return 1
    print(f"OK: all sources at version {uniq.pop()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
