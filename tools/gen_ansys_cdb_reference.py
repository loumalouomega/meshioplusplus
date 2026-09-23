"""Freeze mapdl-archive's reading of the Ansys ``.cdb`` fixtures.

Writes ``tests/python/meshes/ansys/mapdl_archive_reference.npz``: for each
fixture, the VTK grid mapdl-archive builds (points, cell types, offsets and
connectivity) and the sizes of its node and element components.
``test_ansysInp.py`` compares meshio++'s cells against it as sets of node
coordinates, so the node order inside a cell does not matter and the check runs
without mapdl-archive installed.

mapdl-archive (MIT) is not a dependency. Run this in a throw-away environment:

    uv venv --python 3.12 /tmp/mapdl && uv pip install --python /tmp/mapdl/bin/python \
        mapdl-archive==0.4.2 pyvista==0.49.0
    /tmp/mapdl/bin/python tools/gen_ansys_cdb_reference.py
"""

import pathlib

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
FIXTURES = HERE.parent / "tests" / "python" / "meshes" / "ansys"


def main():
    import mapdl_archive

    out = {}
    for path in sorted(FIXTURES.glob("*.cdb")):
        archive = mapdl_archive.Archive(str(path))
        grid = archive.grid
        key = path.stem
        out[f"{key}/points"] = np.asarray(grid.points, dtype=np.float64)
        out[f"{key}/celltypes"] = np.asarray(grid.celltypes, dtype=np.uint8)
        out[f"{key}/offsets"] = np.asarray(grid.cell_offsets, dtype=np.int64)
        out[f"{key}/connectivity"] = np.asarray(grid.cell_connectivity, dtype=np.int64)
        for kind, comps in (
            ("node", archive.node_components),
            ("elem", archive.element_components),
        ):
            names = sorted(comps)
            out[f"{key}/{kind}_components"] = np.array(names, dtype=str)
            out[f"{key}/{kind}_sizes"] = np.array(
                [len(comps[n]) for n in names], dtype=np.int64
            )
    target = FIXTURES / "mapdl_archive_reference.npz"
    np.savez_compressed(target, **out)
    print(f"wrote {len(list(FIXTURES.glob('*.cdb')))} fixtures to {target}")


if __name__ == "__main__":
    main()
