"""
I/O for the Modulef Formatted Mesh (MFM) format used by FEconv
<https://github.com/victorsndvg/FEconv>.

MFM is an ASCII, single-element-type ("non-hybrid") mesh: a header of eight
integers followed by the connectivity, per-element reference arrays, the vertex
coordinates and a per-element subdomain array.  Only vertex coordinates are
stored, so higher-order (P2) elements would be straight-sided; meshio therefore
supports the linear element types losslessly and rejects the rest.
"""

import numpy as np

from .._exceptions import ReadError, WriteError
from .._files import open_file
from .._mesh import CellBlock, Mesh

__all__ = ["read", "write"]

# meshio linear type -> (lnv local vertices, lne local edges, lnf local faces)
_meshio_topology = {
    "line": (2, 1, 0),
    "triangle": (3, 3, 1),
    "quad": (4, 4, 1),
    "tetra": (4, 6, 4),
    "hexahedron": (8, 12, 6),
    "wedge": (6, 9, 5),
}
_num_nodes = {t: v for t, (v, _, _) in _meshio_topology.items()}


def _type_from_dims(lnv, lne, lnf, lnn):
    for t, (v, e, f) in _meshio_topology.items():
        if (v, e, f) == (lnv, lne, lnf) and _num_nodes[t] == lnn:
            return t
    raise ReadError(
        f"MFM: unsupported element (lnv={lnv}, lne={lne}, lnf={lnf}, lnn={lnn}). "
        "meshio++'s MFM only handles linear elements."
    )


def read(filename):
    with open_file(filename, "r") as f:
        # First non-empty line is the header of (up to) eight integers.
        line = f.readline()
        while line and line.strip() == "":
            line = f.readline()
        header = [int(v) for v in line.split()]
        if len(header) < 8:
            raise ReadError("MFM: expected a header of 8 integers.")
        nel, nnod, nver, dim, lnn, lnv, lne, lnf = header[:8]
        tokens = f.read().split()

    cell_type = _type_from_dims(lnv, lne, lnf, lnn)
    if lnn != lnv or nnod != nver:
        raise ReadError("MFM: only linear (P1) elements are supported.")

    pos = 0

    def take(n):
        nonlocal pos
        chunk = tokens[pos : pos + n]
        pos += n
        return chunk

    # Connectivity (vertex array `mm`), element-major.
    mm = np.array(take(lnv * nel), dtype=int).reshape(nel, lnv)
    # Reference arrays: nrc (faces, dim==3), nra (edges, dim>=2), nrv (vertices).
    if dim == 3:
        take(lnf * nel)  # nrc, discarded
    if dim >= 2:
        take(lne * nel)  # nra, discarded
    take(lnv * nel)  # nrv, discarded
    # Vertex coordinates, vertex-major.
    z = np.array(take(dim * nver), dtype=float).reshape(nver, dim)
    # Per-element subdomain / material reference.
    nsd = np.array(take(nel), dtype=int)

    cells = [CellBlock(cell_type, mm - 1)]
    cell_data = {"mfm:ref": [nsd]}
    return Mesh(z, cells, cell_data=cell_data)


def write(filename, mesh, float_fmt=".16e"):
    # MFM is single-type; collapse (only) same-type blocks.
    types = {c.type for c in mesh.cells}
    if len(types) != 1:
        raise WriteError(
            "MFM can only write a single element type, got "
            f"{', '.join(sorted(types)) or 'none'}."
        )
    cell_type = mesh.cells[0].type
    if cell_type not in _meshio_topology:
        raise WriteError(f"MFM does not support '{cell_type}' cells.")

    data = np.concatenate([c.data for c in mesh.cells])
    nel = data.shape[0]
    lnv, lne, lnf = _meshio_topology[cell_type]
    lnn = lnv
    points = mesh.points
    nver = nnod = points.shape[0]
    dim = points.shape[1]

    # Per-element reference (subdomain); default to 1, or use mfm:ref.
    if "mfm:ref" in mesh.cell_data:
        nsd = np.concatenate(mesh.cell_data["mfm:ref"]).astype(int)
    else:
        nsd = np.ones(nel, dtype=int)

    with open_file(filename, "w") as f:
        f.write(f"{nel} {nnod} {nver} {dim} {lnn} {lnv} {lne} {lnf}\n")
        # mm (vertex connectivity, 1-based)
        np.savetxt(f, (data + 1).reshape(nel, lnv), fmt="%d")
        # reference arrays (zeros): nrc (dim==3), nra (dim>=2), nrv
        if dim == 3:
            np.savetxt(f, np.zeros((nel, lnf), dtype=int), fmt="%d")
        if dim >= 2:
            np.savetxt(f, np.zeros((nel, lne), dtype=int), fmt="%d")
        np.savetxt(f, np.zeros((nel, lnv), dtype=int), fmt="%d")
        # vertex coordinates
        np.savetxt(f, points, fmt="%" + float_fmt)
        # subdomain array
        np.savetxt(f, nsd.reshape(nel, 1), fmt="%d")
