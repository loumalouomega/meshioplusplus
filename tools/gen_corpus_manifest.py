"""Write ``tests/python/meshes/MANIFEST.json``, the list of every reference fixture.

Each fixture under ``tests/python/meshes/`` and ``tests/python/input/`` gets a row
with its size, its SHA-256, the README that documents its directory, and its
``producer``, ``licence`` and ``source`` provenance metadata (roadmap 3.20). A
fixture that is a Git LFS pointer (a checkout without ``git lfs pull``) is
recorded with the ``oid`` and ``size`` the pointer carries, which are the hash
and size of the real file, so the manifest is the same with or without the LFS
content::

    python tools/gen_corpus_manifest.py            # rewrite the manifest
    python tools/gen_corpus_manifest.py --check    # exit 1 if it is stale

Every fixture directory must be registered in :data:`DIRECTORY_METADATA` with a
non-empty ``producer``, ``licence`` and ``source``; adding a fixture under an
unregistered directory fails the generator and ``test_corpus_manifest.py``.
Every directory marked ``third_party=True`` must also be attributed in
``CITATION.cff``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import re
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]
ROOT = REPO / "tests" / "python"
SECTIONS = ("meshes", "input")
MANIFEST = ROOT / "meshes" / "MANIFEST.json"
_POINTER_RE = re.compile(
    rb"^version https://git-lfs\.github\.com/spec/v1\noid sha256:([0-9a-f]{64})\nsize (\d+)\n"
)
_SKIP_NAMES = {"MANIFEST.json", ".gitattributes", ".DS_Store"}

# Longest-prefix provenance metadata per fixture directory (roadmap 3.20).
# Each entry maps a repository-relative directory prefix to:
#   producer:    tool and version (or generator script) that wrote the fixtures
#   licence:     SPDX identifier or license name governing the fixtures
#   source:      upstream repository URL, generator script, or "meshio++ test suite"
#   third_party: True if redistributed from an external project (checked in CITATION.cff)
DIRECTORY_METADATA: dict[str, dict[str, object]] = {
    "tests/python/input/mdpa": {
        "producer": "meshio++ / Kratos Multiphysics MDPA test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/abaqus": {
        "producer": "Abaqus / meshio++ test fixtures",
        "licence": "MIT",
        "source": "tests/python/meshes/abaqus/README.md",
        "third_party": False,
    },
    "tests/python/meshes/abaqus_fil": {
        "producer": "tools/gen_abaqus_fil_fixtures.py",
        "licence": "MIT",
        "source": "tools/gen_abaqus_fil_fixtures.py",
        "third_party": False,
    },
    "tests/python/meshes/abaqus_fil/bertoldi": {
        "producer": "Abaqus 6.23",
        "licence": "MIT",
        "source": "https://github.com/bertoldi-collab/inverse-design-textile-metamaterials",
        "third_party": True,
    },
    "tests/python/meshes/abaqus_fil/cjekel": {
        "producer": "Abaqus 6.14",
        "licence": "MIT",
        "source": "https://github.com/cjekel/Introduction-to-Python-Numerical-Analysis-for-Engineers-and-Scientist",
        "third_party": True,
    },
    "tests/python/meshes/abaqus_fil/pybaqus": {
        "producer": "Abaqus 2023 (6.23-1)",
        "licence": "MIT",
        "source": "https://github.com/cristobaltapia/pybaqus",
        "third_party": True,
    },
    "tests/python/meshes/ansys": {
        "producer": "Ansys MAPDL / Workbench / HyperMesh (via mapdl-archive)",
        "licence": "MIT",
        "source": "https://github.com/akaszynski/mapdl-archive",
        "third_party": True,
    },
    "tests/python/meshes/ansys/rst": {
        "producer": "Ansys MAPDL (via pymapdl-reader)",
        "licence": "MIT",
        "source": "https://github.com/ansys/pymapdl-reader",
        "third_party": True,
    },
    "tests/python/meshes/cad_json": {
        "producer": "tools/gen_cad_json_fixtures.py",
        "licence": "MIT",
        "source": "tools/gen_cad_json_fixtures.py",
        "third_party": False,
    },
    "tests/python/meshes/cgns": {
        "producer": "cgnslib",
        "licence": "MIT",
        "source": "tests/python/meshes/cgns/README.md",
        "third_party": False,
    },
    "tests/python/meshes/code_aster": {
        "producer": "Code_Aster / meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/comsol": {
        "producer": "COMSOL Multiphysics / meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/elmer": {
        "producer": "ElmerGrid / meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/ensight": {
        "producer": "EnSight / meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/exodus": {
        "producer": "PeriLab (DLR) / Exodus II",
        "licence": "BSD-3-Clause",
        "source": "tests/python/meshes/exodus/README.md",
        "third_party": True,
    },
    "tests/python/meshes/febio": {
        "producer": "FEBio / meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/femap": {
        "producer": "Femap / EMSolution / MYSTRAN (via tools/gen_femap_fixtures.py)",
        "licence": "MIT",
        "source": "https://github.com/FrontISTR/FrontISTR",
        "third_party": True,
    },
    "tests/python/meshes/flac3d": {
        "producer": "FLAC3D / meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/flux": {
        "producer": "Altair Flux / meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/frd": {
        "producer": "CalculiX ccx / meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/gid": {
        "producer": "GiD / meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/gmsh": {
        "producer": "Gmsh / meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/libmesh": {
        "producer": "tools/gen_libmesh_fixtures.py (libMesh-1.3.0 layout)",
        "licence": "MIT",
        "source": "tools/gen_libmesh_fixtures.py",
        "third_party": False,
    },
    "tests/python/meshes/lsdyna": {
        "producer": "LS-DYNA / meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/lsdyna_binout": {
        "producer": "LS-DYNA / lasso-python LSDA writer",
        "licence": "MIT",
        "source": "https://github.com/ansys/example-data",
        "third_party": True,
    },
    "tests/python/meshes/lsdyna_d3plot": {
        "producer": "tools/gen_d3plot_reference.py (lasso-python 2.0.4)",
        "licence": "MIT",
        "source": "tools/gen_d3plot_reference.py",
        "third_party": False,
    },
    "tests/python/meshes/lsdyna_d3plot/dyna": {
        "producer": "LS-DYNA (via ansys/example-data)",
        "licence": "MIT",
        "source": "https://github.com/ansys/example-data",
        "third_party": True,
    },
    "tests/python/meshes/lsdyna_d3plot/lasso": {
        "producer": "LS-DYNA (via lasso-python)",
        "licence": "BSD-3-Clause",
        "source": "https://github.com/open-lasso-python/lasso-python",
        "third_party": True,
    },
    "tests/python/meshes/marc": {
        "producer": "tools/gen_marc_fixtures.py",
        "licence": "MIT",
        "source": "tools/gen_marc_fixtures.py",
        "third_party": False,
    },
    "tests/python/meshes/mdpa": {
        "producer": "Kratos Multiphysics / meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/med": {
        "producer": "Salome 9.2.2 / Code_Aster 13.6",
        "licence": "MIT",
        "source": "tests/python/meshes/med/README.md",
        "third_party": False,
    },
    "tests/python/meshes/medit": {
        "producer": "Medit / meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/mfem": {
        "producer": "MFEM / PyMFEM 4.10",
        "licence": "BSD-3-Clause",
        "source": "https://github.com/mfem/mfem",
        "third_party": True,
    },
    "tests/python/meshes/msh": {
        "producer": "Gmsh 4.2.2",
        "licence": "MIT",
        "source": "tests/python/meshes/msh/README.md",
        "third_party": False,
    },
    "tests/python/meshes/nastran": {
        "producer": "HyperMesh 2017.3 / OptiStruct (pyNastran)",
        "licence": "BSD-3-Clause",
        "source": "https://github.com/SteveDoyle2/pyNastran",
        "third_party": True,
    },
    "tests/python/meshes/nastran_h5": {
        "producer": "MSC Nastran 2020 (msc20200)",
        "licence": "BSD-3-Clause",
        "source": "https://github.com/SteveDoyle2/pyNastran",
        "third_party": True,
    },
    "tests/python/meshes/nastran_op2": {
        "producer": "MSC Nastran / Simcenter NX Nastran",
        "licence": "BSD-3-Clause",
        "source": "https://github.com/SteveDoyle2/pyNastran",
        "third_party": True,
    },
    "tests/python/meshes/netgen": {
        "producer": "Netgen / meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/neuroglancer": {
        "producer": "Neuroglancer / meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/obj": {
        "producer": "meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/off": {
        "producer": "meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/patran": {
        "producer": "meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/patran/real": {
        "producer": "P3/PATRAN 3.0 / PATRAN 2.5 (WARP3D / Tahoe)",
        "licence": "NCSA AND BSD-3-Clause",
        "source": "https://github.com/rhdodds/warp3d",
        "third_party": True,
    },
    "tests/python/meshes/pcd": {
        "producer": "Point Cloud Library (PCL) test corpus",
        "licence": "BSD-3-Clause",
        "source": "https://github.com/PointCloudLibrary/pcl",
        "third_party": True,
    },
    "tests/python/meshes/ply": {
        "producer": "meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/radioss": {
        "producer": "tools/gen_radioss_fixtures.py / Gmsh 4.15",
        "licence": "MIT",
        "source": "tools/gen_radioss_fixtures.py",
        "third_party": False,
    },
    "tests/python/meshes/radioss_anim": {
        "producer": "OpenRadioss / meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/radioss_th/column": {
        "producer": "OpenRadioss latest-20260728",
        "licence": "MIT",
        "source": "tests/python/meshes/radioss_th/column/README.md",
        "third_party": False,
    },
    "tests/python/meshes/su2": {
        "producer": "SU2 / mixgrid_cube_v5",
        "licence": "MIT",
        "source": "https://su2code.github.io/docs_v7/Mesh-File/",
        "third_party": False,
    },
    "tests/python/meshes/szplt": {
        "producer": "tools/gen_szplt_fixtures.cpp (TecIO 2018.3, SZL v105)",
        "licence": "MIT",
        "source": "tools/gen_szplt_fixtures.cpp",
        "third_party": False,
    },
    "tests/python/meshes/tecplot": {
        "producer": "Tecplot / meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/tecplot/plt": {
        "producer": "tools/gen_tecplot_plt_fixtures.py (TecIO #!TDV112)",
        "licence": "MIT",
        "source": "tools/gen_tecplot_plt_fixtures.py",
        "third_party": False,
    },
    "tests/python/meshes/tecplot/visit": {
        "producer": "Tecplot #!TDV71-#!TDV108 (via VisIt test data)",
        "licence": "BSD-3-Clause",
        "source": "https://github.com/visit-dav/visit",
        "third_party": True,
    },
    "tests/python/meshes/tetgen": {
        "producer": "TetGen / meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/ugrid": {
        "producer": "meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/unv": {
        "producer": "I-deas UNV / meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/vtk": {
        "producer": "VTK / meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/vtu": {
        "producer": "VTK / meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/vtx": {
        "producer": "tools/gen_vtx_fixtures.py (DOLFINx 0.11.0, ADIOS2 2.12.1)",
        "licence": "MIT",
        "source": "tools/gen_vtx_fixtures.py",
        "third_party": False,
    },
    "tests/python/meshes/wkt": {
        "producer": "meshio++ test fixtures",
        "licence": "MIT",
        "source": "meshio++ test suite",
        "third_party": False,
    },
    "tests/python/meshes/z88": {
        "producer": "tools/gen_z88_fixtures.py / Z88OS V15 Z88R",
        "licence": "MIT",
        "source": "tools/gen_z88_fixtures.py",
        "third_party": False,
    },
}

_SORTED_PREFIXES = sorted(DIRECTORY_METADATA, key=len, reverse=True)


def metadata_for(rel: str) -> dict[str, object]:
    """Return provenance metadata for repository-relative fixture path ``rel``."""
    for prefix in _SORTED_PREFIXES:
        if rel == prefix or rel.startswith(prefix + "/"):
            return DIRECTORY_METADATA[prefix]
    raise ValueError(
        f"Fixture {rel!r} has no DIRECTORY_METADATA entry in tools/gen_corpus_manifest.py"
    )


def third_party_directories() -> list[str]:
    """Return every fixture directory marked ``third_party=True``."""
    return sorted(
        prefix
        for prefix, meta in DIRECTORY_METADATA.items()
        if bool(meta.get("third_party"))
    )


def fixture_paths() -> list[str]:
    """Every fixture file as a sorted list of repository-relative POSIX paths."""
    found = []
    for section in SECTIONS:
        base = ROOT / section
        if not base.is_dir():
            continue
        for path in base.rglob("*"):
            if not path.is_file() or "__pycache__" in path.parts:
                continue
            if path.name in _SKIP_NAMES or path.name.upper().startswith("README"):
                continue
            found.append(path.relative_to(REPO).as_posix())
    return sorted(found)


def _digest(path: pathlib.Path) -> tuple[str, int, bool]:
    data = path.read_bytes()
    match = _POINTER_RE.match(data)
    if match:
        return match.group(1).decode(), int(match.group(2)), True
    return hashlib.sha256(data).hexdigest(), len(data), False


def _readme_for(path: str) -> str | None:
    """The README closest above the fixture, within its section, or ``None``."""
    section_root = ROOT / pathlib.PurePosixPath(path).parts[2]
    folder = (REPO / path).parent
    while True:
        for name in ("README.md", "README"):
            candidate = folder / name
            if candidate.is_file():
                return candidate.relative_to(REPO).as_posix()
        if folder == section_root or folder == REPO:
            return None
        folder = folder.parent


def build() -> dict[str, object]:
    rows = []
    for rel in fixture_paths():
        sha, size, _pointer = _digest(REPO / rel)
        meta = metadata_for(rel)
        rows.append(
            {
                "path": rel,
                "bytes": size,
                "sha256": sha,
                "readme": _readme_for(rel),
                "producer": meta["producer"],
                "licence": meta["licence"],
                "source": meta["source"],
            }
        )
    return {"schema_version": 2, "fixtures": rows}


def render(doc: dict[str, object]) -> str:
    return json.dumps(doc, indent=1) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--check", action="store_true", help="exit 1 if the manifest is stale"
    )
    args = parser.parse_args(argv)
    text = render(build())
    if args.check:
        current = MANIFEST.read_text(encoding="utf-8") if MANIFEST.is_file() else ""
        if current != text:
            print(
                f"{MANIFEST.relative_to(REPO)} is stale: run tools/gen_corpus_manifest.py"
            )
            return 1
        return 0
    MANIFEST.write_text(text, encoding="utf-8")
    print(
        f"wrote {MANIFEST.relative_to(REPO)} ({len(json.loads(text)['fixtures'])} fixtures)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
