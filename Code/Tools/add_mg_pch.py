#!/usr/bin/env python3
"""
Add per-file <PrecompiledHeader>Create</PrecompiledHeader> entries for the
Debug_MG / Release_MG configurations.

add_mg_configs.py clones the conditioned ItemDefinitionGroups, but the "which
.cpp creates the PCH" setting lives as a per-file override inside a <ClCompile>
item, which that script does not touch. Without this, an _MG build fails with
C1083 "Cannot open precompiled header file".

Idempotent.
"""

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
CODE = REPO / "Code"

PAIRS = [("Debug", "Debug_MG"), ("Release", "Release_MG")]


def process(path):
    text = path.read_text(encoding="utf-8-sig")
    original = text

    for base, mg in PAIRS:
        base_line = (
            "<PrecompiledHeader Condition=\"'$(Configuration)|$(Platform)'=='"
            + base
            + "|x64'\">Create</PrecompiledHeader>"
        )
        mg_line = (
            "<PrecompiledHeader Condition=\"'$(Configuration)|$(Platform)'=='"
            + mg
            + "|x64'\">Create</PrecompiledHeader>"
        )
        if mg_line in text:
            continue
        if base_line not in text:
            continue
        # Preserve the leading indentation of the matched line.
        pattern = re.compile(r"([ \t]*)" + re.escape(base_line))
        m = pattern.search(text)
        indent = m.group(1) if m else "      "
        text = text.replace(base_line, base_line + "\n" + indent + mg_line, 1)

    if text == original:
        return f"UNCHANGED {path.name}"
    path.write_text(text, encoding="utf-8")
    return f"OK {path.name}"


def main():
    targets = sorted(CODE.rglob("*.vcxproj"))
    if not targets:
        print("no vcxproj found", file=sys.stderr)
        return 1
    for p in targets:
        print(process(p))
    return 0


if __name__ == "__main__":
    sys.exit(main())
