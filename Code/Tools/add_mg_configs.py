#!/usr/bin/env python3
"""
Add Debug_MG / Release_MG x64 configurations to the CudaLatticeGauge VS projects.

These configurations carry the multi-GPU build: they define _CLG_MULTI_GPU=1 and
reference MS-MPI. The stock Debug/Release configurations are left untouched and
never reference MPI, so a machine without MS-MPI installed can still open the
solution and build the single-GPU configurations (see Docs/MultiGPU-Plan.md 9.1).

Idempotent: running twice does not duplicate anything.
"""

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]   # Tools -> Code -> repo root
CODE = REPO / "Code"

# Release_MG derives from Release, Debug_MG from Debug.
PAIRS = [("Debug", "Debug_MG"), ("Release", "Release_MG")]

MPI_INC = "$(MSMPI_INC);"
MPI_LIB_DIR = "$(MSMPI_LIB64);"
MPI_LIB = "msmpi.lib;"


def clone_project_configurations(text):
    """Add <ProjectConfiguration> entries for the _MG configs."""
    for base, mg in PAIRS:
        if f'ProjectConfiguration Include="{mg}|x64"' in text:
            continue
        block = (
            f'    <ProjectConfiguration Include="{base}|x64">\n'
            f"      <Configuration>{base}</Configuration>\n"
            f"      <Platform>x64</Platform>\n"
            f"    </ProjectConfiguration>\n"
        )
        if block not in text:
            return None, f"could not find ProjectConfiguration block for {base}|x64"
        new_block = block + (
            f'    <ProjectConfiguration Include="{mg}|x64">\n'
            f"      <Configuration>{mg}</Configuration>\n"
            f"      <Platform>x64</Platform>\n"
            f"    </ProjectConfiguration>\n"
        )
        text = text.replace(block, new_block, 1)
    return text, None


def clone_conditioned_groups(text):
    """
    Duplicate every group conditioned on '<base>|x64' into one conditioned on
    '<mg>|x64'. Covers PropertyGroup / ImportGroup / ItemDefinitionGroup.
    """
    for base, mg in PAIRS:
        base_cond = f"'$(Configuration)|$(Platform)'=='{base}|x64'"
        mg_cond = f"'$(Configuration)|$(Platform)'=='{mg}|x64'"
        if mg_cond in text:
            continue

        # Match whole top-level groups carrying the base condition.
        pattern = re.compile(
            r'([ \t]*<(PropertyGroup|ImportGroup|ItemDefinitionGroup)\b[^>]*'
            + re.escape(base_cond)
            + r'[^>]*>.*?</\2>\n)',
            re.DOTALL,
        )

        out = []
        pos = 0
        for m in pattern.finditer(text):
            group = m.group(1)
            clone = group.replace(base_cond, mg_cond)
            clone = inject_mpi(clone)
            out.append(text[pos:m.end()])
            out.append(clone)
            pos = m.end()
        if not out:
            return None, f"no groups found for {base}|x64"
        out.append(text[pos:])
        text = "".join(out)
    return text, None


def inject_mpi(group):
    """Add MS-MPI include/lib and the _CLG_MULTI_GPU define to an _MG group."""
    # C++ preprocessor defines
    group = re.sub(
        r"(<PreprocessorDefinitions>)",
        r"\1_CLG_MULTI_GPU=1;",
        group,
    )
    # CUDA defines (.cu files need the macro too -- easy to miss)
    group = re.sub(
        r"(<Defines>)",
        r"\1_CLG_MULTI_GPU=1;",
        group,
    )
    # Header search paths, for both cl and nvcc
    group = re.sub(
        r"(<AdditionalIncludeDirectories>)",
        r"\g<1>" + MPI_INC,
        group,
    )
    # Link inputs
    group = re.sub(
        r"(<AdditionalLibraryDirectories>)",
        r"\g<1>" + MPI_LIB_DIR,
        group,
    )
    group = re.sub(
        r"(<AdditionalDependencies>)",
        r"\g<1>" + MPI_LIB,
        group,
    )
    return group


def process(path):
    text = path.read_text(encoding="utf-8-sig")
    original = text

    text, err = clone_project_configurations(text)
    if err:
        return f"SKIP {path.name}: {err}"
    text, err = clone_conditioned_groups(text)
    if err:
        return f"SKIP {path.name}: {err}"

    if text == original:
        return f"UNCHANGED {path.name} (already has _MG)"
    path.write_text(text, encoding="utf-8")
    return f"OK {path.name}"


def main():
    targets = sorted(CODE.rglob("*.vcxproj"))
    if not targets:
        print("no vcxproj found", file=sys.stderr)
        return 1
    only = sys.argv[1:]
    for p in targets:
        if only and p.stem not in only:
            continue
        print(process(p))
    return 0


if __name__ == "__main__":
    sys.exit(main())
