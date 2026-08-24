"""
I/O for the FLUX ``.pf3`` mesh format, following FEconv
<https://github.com/victorsndvg/FEconv>.

ASCII with French keyword headers.  Each element is a 12-integer record
(the 7th field is the type descriptor, the 8th the node count) followed by its
1-based connectivity; node coordinates live under ``COORDONNEES DES NOEUDS``.
Per-element region references are exposed as ``cell_data["pf3:ref"]``.

Node ordering within an element uses meshio's convention directly; this
round-trips losslessly but is not guaranteed identical to FLUX's internal
ordering for every element type.
"""

import numpy as np

from .._common import warn
from .._exceptions import ReadError
from .._files import open_file
from .._mesh import CellBlock, Mesh
from .._provenance import TAG as _PROVENANCE_TAG

__all__ = ["read", "write"]

# PF3 type descriptor (field 7) -> meshio type
_desc3_to_meshio = {
    2: "vertex",
    3: "line",
    4: "line3",
    5: "triangle",
    6: "triangle6",
    7: "quad",
    8: "quad8",
    10: "tetra",
    11: "tetra10",
    12: "wedge",
    13: "wedge15",
    15: "hexahedron",
    16: "hexahedron20",
    17: "pyramid",
}
# meshio type -> (desc1, desc2, desc3)
_meshio_to_desc = {
    "vertex": (1, 1, 2),
    "line": (2, 2, 3),
    "line3": (2, 3, 4),
    "triangle": (3, 7, 5),
    "triangle6": (3, 7, 6),
    "quad": (4, 202, 7),
    "quad8": (4, 303, 8),
    "tetra": (5, 4, 10),
    "tetra10": (5, 15, 11),
    "wedge": (6, 207, 12),
    "wedge15": (6, 307, 13),
    "hexahedron": (7, 2202, 15),
    "hexahedron20": (7, 3303, 16),
    "pyramid": (8, 4202, 17),
}
_dim_of = {  # topological dim -> header line index bucket
    0: "point",
    1: "edge",
    2: "sur",
    3: "vol",
}


def _header_value(lines, predicate):
    for line in lines:
        if predicate(line):
            return int(line.split()[0])
    raise ReadError("pf3: missing header field")


def read(filename):
    with open_file(filename, "r") as f:
        lines = f.read().splitlines()

    dim = _header_value(lines, lambda L: "NOMBRE DE DIMENSIONS" in L)
    nel = _header_value(
        lines,
        lambda L: "D'ELEMENTS" in L
        and not any(
            s in L
            for s in ("VOLUMIQUES", "SURFACIQUES", "LINEIQUES", "PONCTUELS", "MACRO")
        ),
    )
    nnod = _header_value(
        lines, lambda L: "NOMBRE DE POINTS" in L and "INTEGRATION" not in L
    )

    di = next(i for i, L in enumerate(lines) if "DESCRIPTEUR DE TOPOLOGIE" in L)
    ci = next(i for i, L in enumerate(lines) if "COORDONNEES DES NOEUDS" in L)

    etok = " ".join(lines[di + 1 : ci]).split()
    pos = 0
    groups = {}
    refs = {}
    for _ in range(nel):
        hdr = etok[pos : pos + 12]
        pos += 12
        ref = int(hdr[3])
        desc3 = int(hdr[6])
        lnn = int(hdr[7])
        nodes = [int(etok[pos + j]) - 1 for j in range(lnn)]
        pos += lnn
        if desc3 not in _desc3_to_meshio:
            raise ReadError(f"pf3: unknown element descriptor {desc3}")
        mtype = _desc3_to_meshio[desc3]
        groups.setdefault(mtype, []).append(nodes)
        refs.setdefault(mtype, []).append(ref)

    ctok = " ".join(lines[ci + 1 :]).split()
    points = np.empty((nnod, dim))
    p = 0
    for i in range(nnod):
        p += 1  # node index
        for j in range(dim):
            points[i, j] = float(ctok[p])
            p += 1

    cells = []
    cell_data = {"pf3:ref": []}
    for mtype, conn in groups.items():
        cells.append(CellBlock(mtype, np.array(conn, dtype=int)))
        cell_data["pf3:ref"].append(np.array(refs[mtype], dtype=int))

    return Mesh(points, cells, cell_data=cell_data if cells else {})


def write(filename, mesh):
    dim = mesh.points.shape[1]
    from .._mesh import topological_dimension

    counts = {"vol": 0, "sur": 0, "edge": 0, "point": 0}
    blocks = []
    for k, cb in enumerate(mesh.cells):
        if cb.type not in _meshio_to_desc:
            warn(f"pf3 does not support '{cb.type}' cells. Skipping.")
            continue
        td = topological_dimension[cb.type]
        counts[_dim_of[td]] += len(cb.data)
        blocks.append((k, cb))
    nel = sum(len(cb.data) for _, cb in blocks)

    ref_data = mesh.cell_data.get("pf3:ref")

    with open_file(filename, "w") as f:
        f.write(f" {_PROVENANCE_TAG}\n")
        f.write(f"{dim:8d}           NOMBRE DE DIMENSIONS DU DECOUPAGE\n")
        f.write(f"{nel:8d}           NOMBRE  D'ELEMENTS\n")
        f.write(f"{counts['vol']:8d}           NOMBRE  D'ELEMENTS VOLUMIQUES\n")
        f.write(f"{counts['sur']:8d}           NOMBRE  D'ELEMENTS SURFACIQUES\n")
        f.write(f"{counts['edge']:8d}           NOMBRE  D'ELEMENTS LINEIQUES\n")
        f.write(f"{counts['point']:8d}           NOMBRE  D'ELEMENTS PONCTUELS\n")
        f.write(f"{0:8d}           NOMBRE DE MACRO-ELEMENTS\n")
        f.write(f"{len(mesh.points):8d}           NOMBRE DE POINTS\n")
        f.write(f"{1:8d}           NOMBRE DE REGIONS\n")
        f.write(f"{0:8d}           NOMBRE DE REGIONS VOLUMIQUES\n")
        f.write(f"{0:8d}           NOMBRE DE REGIONS SURFACIQUES\n")
        f.write(f"{0:8d}           NOMBRE DE REGIONS LINEIQUES\n")
        f.write(f"{0:8d}           NOMBRE DE REGIONS PONCTUELLES\n")
        f.write(f"{0:8d}           NOMBRE DE REGIONS MACRO-ELEMENTAIRES\n")
        f.write(f"{20:8d}           NOMBRE DE NOEUDS DANS 1 ELEMENT (MAX)\n")
        f.write(f"{20:8d}           NOMBRE DE POINTS D'INTEGRATION / ELEMENT (MAX)\n")
        f.write(" NOMS DES REGIONS\n")
        f.write(" DESCRIPTEUR DE TOPOLOGIE DES ELEMENTS\n")

        eid = 0
        for k, cb in blocks:
            desc1, desc2, desc3 = _meshio_to_desc[cb.type]
            lnn = cb.data.shape[1]
            block_refs = None
            if ref_data is not None and k < len(ref_data):
                block_refs = ref_data[k]
            for r, row in enumerate(cb.data):
                eid += 1
                ref = int(block_refs[r]) if block_refs is not None else 0
                f.write(
                    f"{eid:8d}{desc1:8d}{desc2:8d}{ref:8d}{lnn:8d}{0:8d}"
                    f"{desc3:8d}{lnn:8d}{0:8d}{0:8d}{0:8d}{0:8d}\n"
                )
                f.write(" ".join(f"{v + 1:8d}" for v in row) + "\n")

        f.write(" COORDONNEES DES NOEUDS\n")
        for i, pt in enumerate(mesh.points):
            f.write(f"{i + 1:8d} " + " ".join(repr(float(x)) for x in pt) + "\n")
