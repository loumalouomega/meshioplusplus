#!/usr/bin/env python3
"""Regenerate ``tests/python/meshes/med/hexa27_penta18.med`` with the MED library.

The file holds one unit ``MED_HEXA27`` and one unit ``MED_PENTA18``, written by
MED-fichier itself through its Python bindings (``med``, packaged with the MED
system library), so the reader is checked against a file it did not write. The
nodes are the meshio++ reference elements put into MED's node order through
the ``med`` tables of ``meshioplusplus._node_order``; the test reads them back
and expects the reference elements again.

MED-fichier 5.0 wrote the committed file. Its ``INFOS_GENERALES`` version is
then set to 4.1.0, because the reader refuses a newer major and
``MEDfileVersionOpen`` cannot create a 4.1 file with that build; the groups and
datasets used here (``NOE/COO``, ``MAI/H27/NOD``, ``MAI/P18/NOD``) are laid out
the same in both.

``med`` is not installable with pip. Run this with a Python that has both
``med`` and ``h5py``, e.g. the project environment with the system packages on
the path (only ``_node_order.py`` is loaded from ``src/python``, by path)::

    PYTHONPATH=/usr/lib/python3.14/site-packages .venv/bin/python \
        tools/gen_med_quadratic_fixture.py
"""

import importlib.util
import pathlib

import h5py
import numpy as np
from med.medenum import (
    MED_ACC_CREAT,
    MED_CARTESIAN,
    MED_CELL,
    MED_FULL_INTERLACE,
    MED_HEXA27,
    MED_NO_DT,
    MED_NO_IT,
    MED_NODAL,
    MED_PENTA18,
    MED_SORT_DTIT,
    MED_UNSTRUCTURED_MESH,
)
from med.medfile import MEDfileClose, MEDfileOpen
from med.medmesh import (
    MEDFLOAT,
    MEDINT,
    MEDmeshCr,
    MEDmeshElementConnectivityWr,
    MEDmeshNodeCoordinateWr,
)

ROOT = pathlib.Path(__file__).resolve().parent.parent

# The _node_order module has no package-level dependencies beyond the standard
# library, so load it directly rather than importing the compiled package.
_spec = importlib.util.spec_from_file_location(
    "_node_order", ROOT / "src" / "python" / "meshioplusplus" / "_node_order.py"
)
_node_order = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_node_order)

OUT = ROOT / "tests" / "python" / "meshes" / "med" / "hexa27_penta18.med"

# fmt: off
HEX = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], dtype=float)
HEX_E = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
         (0, 4), (1, 5), (2, 6), (3, 7)]
HEX_F = [(0, 4, 7, 3), (1, 2, 6, 5), (0, 1, 5, 4), (3, 7, 6, 2), (0, 3, 2, 1), (4, 5, 6, 7)]
WED = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1]], dtype=float)
WED_E = [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)]
WED_F = [(0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]
# fmt: on


def reference(corners, edges, faces, body):
    parts = [corners]
    parts += [(corners[a] + corners[b])[None] / 2 for a, b in edges]
    parts += [corners[list(f)].mean(0)[None] for f in faces]
    if body:
        parts.append(corners.mean(0)[None])
    return np.vstack(parts)


def main():
    hexa = reference(HEX, HEX_E, HEX_F, True)
    penta = reference(WED, WED_E, WED_F, False) + [2.0, 0.0, 0.0]
    points = np.vstack([hexa, penta])
    # meshio[k] = med[to_meshio[k]], so med = meshio[from_meshio].
    hexa_conn = 1 + np.array(_node_order.node_order("med", "hexahedron27").from_meshio)
    penta_conn = 1 + 27 + np.array(_node_order.node_order("med", "wedge18").from_meshio)

    OUT.unlink(missing_ok=True)
    fid = MEDfileOpen(str(OUT), MED_ACC_CREAT)
    MEDmeshCr(
        fid,
        "mesh",
        3,
        3,
        MED_UNSTRUCTURED_MESH,
        "",
        "",
        MED_SORT_DTIT,
        MED_CARTESIAN,
        "XYZ",
        "",
    )
    MEDmeshNodeCoordinateWr(
        fid,
        "mesh",
        MED_NO_DT,
        MED_NO_IT,
        0.0,
        MED_FULL_INTERLACE,
        len(points),
        MEDFLOAT(points.ravel().tolist()),
    )
    for geo, conn in ((MED_HEXA27, hexa_conn), (MED_PENTA18, penta_conn)):
        MEDmeshElementConnectivityWr(
            fid,
            "mesh",
            MED_NO_DT,
            MED_NO_IT,
            0.0,
            MED_CELL,
            geo,
            MED_NODAL,
            MED_FULL_INTERLACE,
            1,
            MEDINT(conn.tolist()),
        )
    MEDfileClose(fid)
    with h5py.File(OUT, "r+") as f:
        info = f["INFOS_GENERALES"].attrs
        for key, value in (("MAJ", 4), ("MIN", 1), ("REL", 0)):
            info[key] = np.int32(value)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
