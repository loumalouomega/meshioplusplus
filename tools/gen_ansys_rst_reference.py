"""Freeze pymapdl-reader's reading of the Ansys ``.rst``/``.rth`` fixtures.

Writes ``tests/python/meshes/ansys/rst/pymapdl_reference.npz``: for each fixture,
the result set times, and for every set the node numbers and the nodal DOF
solution pymapdl-reader returns (with the DOF labels), plus the cells of its VTK
grid (points, cell types, offsets, connectivity). ``test_ansys_rst.py`` compares
meshio++'s arrays against it, so the cross-check runs without pymapdl-reader
installed.

pymapdl-reader (MIT) is not a dependency. Run this in a throw-away environment:

    uv venv --python 3.12 /tmp/pymapdl && uv pip install --python /tmp/pymapdl/bin/python \
        "ansys-mapdl-reader @ git+https://github.com/ansys/pymapdl-reader@488ad8b" \
        pyvista==0.49.0
    cd /tmp && /tmp/pymapdl/bin/python /path/to/meshioplusplus/tools/gen_ansys_rst_reference.py

(The frozen file came from that commit, 0.56.dev0. Run it outside a source
checkout of pymapdl-reader, whose package directory would shadow the install.)

pymapdl-reader sizes a record's values by its length in 4-byte words, so on a
result set that holds only some nodes it reads past the record; none of these
fixtures has such a set (see doc/formats/ansys_rst.md).
"""

import pathlib

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
FIXTURES = HERE.parent / "tests" / "python" / "meshes" / "ansys" / "rst"


def main():
    from ansys.mapdl import reader as pymapdl_reader

    out = {}
    paths = sorted(FIXTURES.glob("*.rst")) + sorted(FIXTURES.glob("*.rth"))
    for path in paths:
        rst = pymapdl_reader.read_binary(str(path))
        key = path.name
        out[f"{key}/times"] = np.asarray(rst.time_values, dtype=np.float64)
        for s in range(rst.n_results):
            nnum, values = rst.nodal_solution(s)
            out[f"{key}/{s}/nnum"] = np.asarray(nnum, dtype=np.int64)
            out[f"{key}/{s}/values"] = np.asarray(values, dtype=np.float64)
            out[f"{key}/{s}/dofs"] = np.array(rst.result_dof(s), dtype=str)
        grid = rst.quadgrid
        out[f"{key}/points"] = np.asarray(grid.points, dtype=np.float64)
        out[f"{key}/celltypes"] = np.asarray(grid.celltypes, dtype=np.uint8)
        out[f"{key}/offsets"] = np.asarray(grid.cell_offsets, dtype=np.int64)
        out[f"{key}/connectivity"] = np.asarray(grid.cell_connectivity, dtype=np.int64)
    target = FIXTURES / "pymapdl_reference.npz"
    np.savez_compressed(target, **out)
    print(f"wrote {len(paths)} fixtures to {target}")


if __name__ == "__main__":
    main()
