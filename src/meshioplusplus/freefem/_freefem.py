"""
I/O for the FreeFem++ ``.msh`` mesh format (as handled by FEconv
<https://github.com/victorsndvg/FEconv>).

ASCII.  2D: header ``nver ntri nedge``, then vertices ``x y ref``, triangles
``a b c ref`` and boundary edges ``a b ref``.  3D: header ``nver ntet ntri``,
vertices ``x y z ref``, tetrahedra ``a b c d ref`` and boundary triangles
``a b c ref``.  All 1-based; every entity carries an integer reference/label,
exposed as ``point_data``/``cell_data`` ``"freefem:ref"``.
"""

import numpy as np

from .._common import warn
from .._exceptions import ReadError, WriteError
from .._files import open_file
from .._mesh import CellBlock, Mesh

__all__ = ["read", "write"]


def _nonblank_lines(f):
    for line in f:
        s = line.split()
        if s:
            yield s


def read(filename):
    with open_file(filename, "r") as f:
        it = _nonblank_lines(f)
        try:
            header = next(it)
        except StopIteration:
            raise ReadError("FreeFem: empty file")
        if len(header) != 3:
            raise ReadError("FreeFem: expected a 3-integer header")
        nver, n_el1, n_el2 = (int(v) for v in header)

        # Vertices; the first one fixes the spatial dimension.
        verts = []
        first = next(it)
        dim = len(first) - 1
        if dim not in (2, 3):
            raise ReadError(f"FreeFem: unexpected vertex dimension {dim}")
        verts.append(first)
        for _ in range(nver - 1):
            verts.append(next(it))

        points = np.array([[float(x) for x in row[:dim]] for row in verts])
        point_ref = np.array([int(row[dim]) for row in verts])

        el1_type, lnv1 = ("triangle", 3) if dim == 2 else ("tetra", 4)
        el2_type, lnv2 = ("line", 2) if dim == 2 else ("triangle", 3)

        def read_block(n, lnv):
            conn = np.empty((n, lnv), dtype=int)
            ref = np.empty(n, dtype=int)
            for k in range(n):
                row = next(it)
                conn[k] = [int(v) for v in row[:lnv]]
                ref[k] = int(row[lnv])
            return conn - 1, ref

        el1_conn, el1_ref = read_block(n_el1, lnv1)
        el2_conn, el2_ref = read_block(n_el2, lnv2)

    cells = []
    cell_ref = []
    if n_el1 > 0:
        cells.append(CellBlock(el1_type, el1_conn))
        cell_ref.append(el1_ref)
    if n_el2 > 0:
        cells.append(CellBlock(el2_type, el2_conn))
        cell_ref.append(el2_ref)

    return Mesh(
        points,
        cells,
        point_data={"freefem:ref": point_ref},
        cell_data={"freefem:ref": cell_ref} if cell_ref else {},
    )


def write(filename, mesh):
    dim = mesh.points.shape[1]
    if dim not in (2, 3):
        raise WriteError(f"FreeFem can only write 2D or 3D meshes, got dim={dim}.")

    el1_type = "triangle" if dim == 2 else "tetra"
    el2_type = "line" if dim == 2 else "triangle"

    ref_blocks = mesh.cell_data.get("freefem:ref", [None] * len(mesh.cells))

    def gather(target_type):
        conns, refs = [], []
        for cb, rb in zip(mesh.cells, ref_blocks):
            if cb.type == target_type:
                conns.append(cb.data)
                refs.append(rb if rb is not None else np.zeros(len(cb.data), dtype=int))
        if conns:
            return np.concatenate(conns), np.concatenate(refs)
        return np.empty((0, 0), dtype=int), np.empty(0, dtype=int)

    el1_conn, el1_ref = gather(el1_type)
    el2_conn, el2_ref = gather(el2_type)

    skipped = {c.type for c in mesh.cells if c.type not in (el1_type, el2_type)}
    if skipped:
        warn(
            f"FreeFem ({dim}D) only supports {el1_type}/{el2_type}. Skipping {skipped}."
        )

    point_ref = mesh.point_data.get(
        "freefem:ref", np.zeros(len(mesh.points), dtype=int)
    )

    with open_file(filename, "w") as f:
        f.write(f"{len(mesh.points)} {len(el1_conn)} {len(el2_conn)}\n")
        for pt, r in zip(mesh.points, point_ref):
            f.write(" ".join(repr(float(x)) for x in pt) + f" {int(r)}\n")
        for conn, ref in ((el1_conn, el1_ref), (el2_conn, el2_ref)):
            for row, r in zip(conn, ref):
                f.write(" ".join(str(v + 1) for v in row) + f" {int(r)}\n")
