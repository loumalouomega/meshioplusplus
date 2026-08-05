"""Regular hexahedron lattices (``meshioplusplus.grid``).

The library's first mesh *generator*: every other operation transforms a mesh you
already have. It exists because :func:`meshioplusplus.voxelize` needs a lattice
anyway, and a lattice is a useful thing to be able to ask for on its own -- as a
background grid, as a test fixture that needs no file, or as the domain of a
synthetic field.

Points run **x fastest, then y, then z**, and cells run in the same order::

    vid(i, j, k) = (k * (ny + 1) + j) * (nx + 1) + i
    cid(i, j, k) = (k *  ny      + j) *  nx      + i

with each cell's eight nodes in the meshio/VTK ``hexahedron`` winding. That
ordering is inherited rather than chosen: it is what this repository's existing
C++/Python byte-identity fixtures already agree on, and it is VTK ImageData's own
layout, so ``cell_data`` reshapes to ``arr[z, y, x]``.

Coordinates are ``origin + index * spacing``, never accumulated -- an accumulating
loop makes the last plane depend on the cell count in the low bits, and would make
this module's implementation observable in the output.

This reference implementation is a byte-for-byte twin of the C++ core, pinned by
``tests/python/test_voxelize.py::test_grid_cpp_matches_python``.

Public API:

* :func:`grid` -- build a regular hexahedron lattice.
"""

from __future__ import annotations

import numpy as np

from ._mesh import Mesh

__all__ = ["grid"]

#: Default cell budget, a little above 256^3. A 512^3 grid is 134 million cells
#: and roughly 11.8 GB of points and connectivity, which is asked for by typo far
#: more often than on purpose.
DEFAULT_MAX_CELLS = 20000000


def _validate(dims, origin, spacing, max_cells):
    dims = np.asarray(dims, dtype=np.int64).reshape(-1)
    if dims.size != 3:
        raise ValueError("meshio++: grid: dims must be three cell counts")
    origin = np.asarray(origin, dtype=np.float64).reshape(-1)
    spacing = np.asarray(spacing, dtype=np.float64).reshape(-1)
    if origin.size != 3 or spacing.size != 3:
        raise ValueError(
            "meshio++: grid: origin and spacing must each have three components"
        )
    for k in range(3):
        if dims[k] < 0:
            raise ValueError(
                f"meshio++: grid: cell counts must not be negative, got {int(dims[k])} "
                f"on axis {k}"
            )
        if dims[k] > 0 and not spacing[k] > 0.0:
            raise ValueError(
                "meshio++: grid: spacing must be positive on every axis with cells, got "
                f"{float(spacing[k])} on axis {k}"
            )
    cells = int(dims[0]) * int(dims[1]) * int(dims[2]) if all(dims > 0) else 0
    if max_cells > 0 and cells > max_cells:
        raise ValueError(
            f"meshio++: grid: the requested grid has {cells} cells, above the limit of "
            f"{max_cells} (raise max_cells, coarsen the resolution, or use a band)"
        )
    return dims, origin, spacing


def _lattice_py(dims, origin, spacing):
    """The lattice itself, as points and connectivity."""
    nx, ny, nz = (int(d) for d in dims)
    if nx <= 0 or ny <= 0 or nz <= 0:
        return np.zeros((0, 3), dtype=np.float64), None

    px, py, pz = nx + 1, ny + 1, nz + 1
    # origin + index * spacing, evaluated per point. np.arange is deliberately
    # not used: it is free to compute its values by accumulation.
    i = np.arange(px, dtype=np.float64)
    j = np.arange(py, dtype=np.float64)
    k = np.arange(pz, dtype=np.float64)
    xs = origin[0] + i * spacing[0]
    ys = origin[1] + j * spacing[1]
    zs = origin[2] + k * spacing[2]
    # x fastest, then y, then z.
    gz, gy, gx = np.meshgrid(zs, ys, xs, indexing="ij")
    points = np.stack([gx.reshape(-1), gy.reshape(-1), gz.reshape(-1)], axis=1)

    ci, cj, ck = np.meshgrid(
        np.arange(nx, dtype=np.int64),
        np.arange(ny, dtype=np.int64),
        np.arange(nz, dtype=np.int64),
        indexing="ij",
    )
    # The cell loop is k outer, then j, then i, so transpose the ij-indexed grids
    # into that order before flattening.
    ci = ci.transpose(2, 1, 0).reshape(-1)
    cj = cj.transpose(2, 1, 0).reshape(-1)
    ck = ck.transpose(2, 1, 0).reshape(-1)
    base = (ck * py + cj) * px + ci
    top = base + px * py
    conn = np.stack(
        [
            base,
            base + 1,
            base + px + 1,
            base + px,
            top,
            top + 1,
            top + px + 1,
            top + px,
        ],
        axis=1,
    ).astype(np.int64)
    return points, conn


def _grid_py(dims, origin, spacing):
    points, conn = _lattice_py(dims, origin, spacing)
    cells = [] if conn is None else [("hexahedron", conn)]
    return Mesh(points, cells)


def lattice_from_mesh(mesh):
    """Recover the lattice a mesh *is*, from its geometry alone.

    The numpy twin of ``detail::lattice_from_mesh``, and the reason
    :mod:`meshioplusplus.vti` can write a grid at all: ImageData stores its
    geometry as ``Origin``/``Spacing``/``WholeExtent`` rather than as points, so
    a writer has to derive those three from the mesh.

    The recovery is exact rather than a fit. :func:`grid` writes every coordinate
    as ``origin + index * spacing`` evaluated independently, so all points on one
    plane carry bit-identical values and an exact ``np.unique`` recovers the
    planes with no tolerance. Only the *spacing* needs one -- consecutive gaps of
    ``origin + i*h`` differ in the last bits -- so uniformity is checked to a
    relative ``1e-9`` and the reported spacing is ``(last - first) / n``.

    Returns
    -------
    tuple or None
        ``(dims, origin, spacing)`` as three length-3 numpy arrays, or ``None``
        when the mesh is not a dense lattice. A *partial* lattice (``voxelize``'s
        ``surface``/``inside`` fills, or an octree) is deliberately not
        recoverable: its points do not tile the box, and inventing the missing
        ones would be a different mesh.
    """
    if len(mesh.cells) != 1:
        return None
    block = mesh.cells[0]
    if block.type != "hexahedron" or not isinstance(block.data, np.ndarray):
        return None
    conn = np.asarray(block.data)
    if conn.ndim != 2 or conn.shape[1] != 8:
        return None
    pts = np.asarray(mesh.points, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] < 3 or pts.shape[0] == 0:
        return None

    planes = [np.unique(pts[:, d]) for d in range(3)]
    if any(p.size < 2 for p in planes):
        return None
    if planes[0].size * planes[1].size * planes[2].size != pts.shape[0]:
        return None
    dims = np.array([p.size - 1 for p in planes], dtype=np.int64)
    if int(dims[0]) * int(dims[1]) * int(dims[2]) != conn.shape[0]:
        return None

    origin = np.array([p[0] for p in planes], dtype=np.float64)
    spacing = np.empty(3, dtype=np.float64)
    for d in range(3):
        n = int(dims[d])
        spacing[d] = (planes[d][-1] - planes[d][0]) / float(n)
        if not spacing[d] > 0.0:
            return None
        gaps = planes[d][1:] - planes[d][:-1]
        if np.any(np.abs(gaps - spacing[d]) > 1.0e-9 * spacing[d]):
            return None

    px, py = int(dims[0]) + 1, int(dims[1]) + 1
    # The points must be in the lattice's own x-fastest order, not merely occupy
    # its plane positions -- a permuted grid has identical plane sets and is a
    # different mesh. Exact comparison: `planes[d]` holds the very values the
    # points carry.
    g = np.arange(pts.shape[0], dtype=np.int64)
    idx = (g % px, (g // px) % py, g // (px * py))
    for d in range(3):
        if not np.array_equal(pts[:, d], planes[d][idx[d]]):
            return None

    nx, ny, nz = (int(v) for v in dims)
    c = np.arange(nx * ny * nz, dtype=np.int64)
    ci, cj, ck = c % nx, (c // nx) % ny, c // (nx * ny)
    base = (ck * py + cj) * px + ci
    top = base + px * py
    expect = np.stack(
        [
            base,
            base + 1,
            base + px + 1,
            base + px,
            top,
            top + 1,
            top + px + 1,
            top + px,
        ],
        axis=1,
    )
    if not np.array_equal(conn.astype(np.int64), expect):
        return None
    return dims, origin, spacing


def grid(
    dims, origin=(0.0, 0.0, 0.0), spacing=(1.0, 1.0, 1.0), max_cells=DEFAULT_MAX_CELLS
):
    """Build a regular hexahedron lattice.

    Parameters
    ----------
    dims :
        Cell counts ``(nx, ny, nz)``. A zero on any axis gives an empty mesh --
        a legal, if useless, request rather than an error.
    origin :
        The lo corner.
    spacing :
        Cell size per axis.
    max_cells :
        Refuse by name above this many cells. See :data:`DEFAULT_MAX_CELLS`.

    Returns
    -------
    Mesh
        One ``hexahedron`` cell block over a shared corner lattice.
    """
    dims, origin, spacing = _validate(dims, origin, spacing, max_cells)

    try:
        from . import _core

        return _core.grid(
            [int(d) for d in dims],
            [float(v) for v in origin],
            [float(v) for v in spacing],
            int(max_cells),
        )
    except (ValueError, TypeError):
        # A genuine user error must not fall through to the numpy path, which
        # would report it differently or not at all.
        raise
    except Exception:
        return _grid_py(dims, origin, spacing)
