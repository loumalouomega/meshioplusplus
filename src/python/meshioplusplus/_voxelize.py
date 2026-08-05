"""Voxelization (``meshioplusplus.voxelize``).

Build a regular grid around a mesh and, optionally, keep only the cells its
surface passes through or encloses.

The output is an **ordinary** :class:`~meshioplusplus._mesh.Mesh` -- one
``hexahedron`` cell block over a shared corner lattice -- and that is the whole
design. Once it is a mesh, every writer already handles it, ``view``/``screenshot``
already render it, ``crop``/``split``/``data calc`` already operate on it,
``--color-by`` already colours it and ``isosurface`` already contours whatever
field it carries. None of that needed new code.

Three fill rules:

``all``
    Every cell of the (optionally padded) bounding box. The mesh contributes only
    its extent, so this is a background grid rather than a shape.
``surface``
    Only cells a surface triangle passes through, by exact triangle/box overlap.
    Needs no closed surface, so it works on an open sheet.
``inside``
    Only cells whose centre is inside the surface. Needs a sign, and therefore a
    surface closed enough for the chosen convention.

This reference implementation is a byte-for-byte twin of the C++ core, pinned by
``tests/python/test_voxelize.py::test_cpp_matches_python``.

Public API:

* :func:`voxelize` -- build a grid around a mesh.
"""

from __future__ import annotations

import numpy as np

from ._grid import DEFAULT_MAX_CELLS, _lattice_py
from ._mesh import Mesh
from ._sdf import _sample_py, _soup

__all__ = ["voxelize"]

_PREFIX = "meshio++: voxelize: "
#: The same string, under a name the shared resolver can shadow safely.
_VOX_PREFIX = _PREFIX
_FILLS = ("all", "surface", "inside")


def _tri_box_overlap(centre, half, a, b, c):
    """Exact triangle / axis-aligned-box overlap by the separating axis theorem.

    Transcribed from ``detail/tri_box.hpp``: thirteen axes, in the same order.
    A bounding-box test is not this -- a long diagonal triangle overlaps far more
    boxes than it enters.
    """
    v = np.array([a - centre, b - centre, c - centre])
    e = np.array([v[1] - v[0], v[2] - v[1], v[0] - v[2]])

    def separated(p, rad):
        return p.min() > rad or p.max() < -rad

    for k in range(3):
        ex, ey, ez = abs(e[k][0]), abs(e[k][1]), abs(e[k][2])
        p = -e[k][2] * v[:, 1] + e[k][1] * v[:, 2]
        if separated(p, ez * half[1] + ey * half[2]):
            return False
        p = e[k][2] * v[:, 0] - e[k][0] * v[:, 2]
        if separated(p, ez * half[0] + ex * half[2]):
            return False
        p = -e[k][1] * v[:, 0] + e[k][0] * v[:, 1]
        if separated(p, ey * half[0] + ex * half[1]):
            return False

    for k in range(3):
        if separated(v[:, k], half[k]):
            return False

    normal = np.cross(e[0], e[1])
    vmin = np.where(normal > 0.0, -half - v[0], half - v[0])
    vmax = np.where(normal > 0.0, half - v[0], -half - v[0])
    if np.sum(normal * vmin) > 0.0:
        return False
    return np.sum(normal * vmax) >= 0.0


def _resolve_bounds(mesh, bounds, padding, padding_relative, _PREFIX=None):
    _PREFIX = _VOX_PREFIX if _PREFIX is None else _PREFIX
    if bounds is not None:
        b = np.asarray(bounds, dtype=np.float64).reshape(-1)
        if b.size != 6:
            raise ValueError(f"{_PREFIX}bounds must be (xlo, ylo, zlo, xhi, yhi, zhi)")
        lo, hi = b[:3].copy(), b[3:].copy()
        for k in range(3):
            if not hi[k] >= lo[k]:
                raise ValueError(
                    f"{_PREFIX}bounds are inverted on axis {k} "
                    f"(lo {lo[k]} > hi {hi[k]})"
                )
    else:
        pts = np.asarray(mesh.points, dtype=np.float64)
        if len(pts) == 0:
            raise ValueError(
                f"{_PREFIX}the mesh has no points, so it has no bounding box to cover "
                "(pass explicit bounds)"
            )
        if pts.shape[1] == 2:
            pts = np.column_stack([pts, np.zeros(len(pts))])
        lo, hi = pts.min(axis=0), pts.max(axis=0)

    diag = float(np.sqrt(np.sum((hi - lo) ** 2)))
    pad = padding + padding_relative * diag
    if pad < 0.0:
        raise ValueError(f"{_PREFIX}padding is negative")
    return lo - pad, hi + pad


def _resolve_lattice(
    mesh,
    resolution,
    cell_size,
    bounds,
    padding,
    padding_relative,
    max_cells,
    _PREFIX=None,
):
    """The numpy twin of ``detail::lattice_resolve``.

    The prefix is a parameter for the same reason it is in C++: ``compute_sdf``
    resolves the identical six fields the identical way, and its errors must name
    ``sdf`` rather than ``voxelize``.
    """
    _PREFIX = _VOX_PREFIX if _PREFIX is None else _PREFIX
    if (resolution is None) == (cell_size is None):
        raise ValueError(f"{_PREFIX}give exactly one of resolution and cell_size")
    lo, hi = _resolve_bounds(mesh, bounds, padding, padding_relative, _PREFIX)

    if resolution is not None:
        dims = np.asarray(resolution, dtype=np.int64).reshape(-1)
        if dims.size != 3:
            raise ValueError(f"{_PREFIX}resolution must be three cell counts")
        for k in range(3):
            if dims[k] <= 0:
                raise ValueError(
                    f"{_PREFIX}resolution must be positive on every axis, got "
                    f"{int(dims[k])} on axis {k}"
                )
        spacing = np.array(
            [(hi[k] - lo[k]) / float(dims[k]) for k in range(3)], dtype=np.float64
        )
    else:
        cell = float(cell_size)
        if not cell > 0.0:
            raise ValueError(f"{_PREFIX}cell_size must be positive")
        dims = np.zeros(3, dtype=np.int64)
        for k in range(3):
            extent = hi[k] - lo[k]
            dims[k] = int(np.ceil(extent / cell)) if extent > 0.0 else 0
            if dims[k] <= 0:
                raise ValueError(
                    f"{_PREFIX}the bounding box is degenerate on axis {k}, so a cell size "
                    "cannot fill it (pass an explicit resolution or bounds)"
                )
        spacing = np.full(3, cell, dtype=np.float64)

    cells = int(dims[0]) * int(dims[1]) * int(dims[2])
    if max_cells > 0 and cells > max_cells:
        raise ValueError(
            f"{_PREFIX}the requested grid has {cells} cells, above the limit of {max_cells} "
            "(raise max_cells, coarsen the resolution, or use a band)"
        )
    return lo, spacing, dims


def _voxelize_py(
    mesh,
    resolution,
    cell_size,
    bounds,
    padding,
    padding_relative,
    fill,
    attach_occupancy,
    max_cells,
    sign,
):
    origin, spacing, dims = _resolve_lattice(
        mesh, resolution, cell_size, bounds, padding, padding_relative, max_cells
    )
    nx, ny, nz = (int(d) for d in dims)
    ncells = nx * ny * nz

    occupied = None
    if fill == "surface":
        occupied = np.zeros(ncells, dtype=bool)
        _points, _verts, corners, _source = _soup(mesh)
        half = spacing * 0.5
        for t in range(len(corners)):
            a, b, c = corners[t]
            tlo = np.minimum(np.minimum(a, b), c)
            thi = np.maximum(np.maximum(a, b), c)
            lo_i = np.maximum(np.floor((tlo - origin) / spacing), 0).astype(np.int64)
            hi_i = np.minimum(np.floor((thi - origin) / spacing), dims - 1).astype(
                np.int64
            )
            if np.any(hi_i < lo_i):
                continue
            for k in range(lo_i[2], hi_i[2] + 1):
                for j in range(lo_i[1], hi_i[1] + 1):
                    for i in range(lo_i[0], hi_i[0] + 1):
                        cid = (k * ny + j) * nx + i
                        if occupied[cid]:
                            continue
                        centre = origin + (np.array([i, j, k]) + 0.5) * spacing
                        if _tri_box_overlap(centre, half, a, b, c):
                            occupied[cid] = True
    elif fill == "inside":
        if sign == "unsigned":
            raise ValueError(
                f"{_PREFIX}fill 'inside' needs a sign, but sign='unsigned' was given"
            )
        ii, jj, kk = np.meshgrid(
            np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij"
        )
        idx = np.stack(
            [
                ii.transpose(2, 1, 0).reshape(-1),
                jj.transpose(2, 1, 0).reshape(-1),
                kk.transpose(2, 1, 0).reshape(-1),
            ],
            axis=1,
        ).astype(np.float64)
        centres = origin + (idx + 0.5) * spacing
        dist, _cells, _band = _sample_py(mesh, centres, sign, "angle", 0.0)
        occupied = dist < 0.0

    points, conn = _lattice_py(dims, origin, spacing)
    if occupied is None:
        out = Mesh(points, [] if conn is None else [("hexahedron", conn)])
        kept = ncells
    else:
        keep = np.flatnonzero(occupied)
        kept = int(keep.size)
        sub = (
            conn[keep]
            if conn is not None and kept
            else np.zeros((0, 8), dtype=np.int64)
        )
        # Compact the referenced points in ascending order -- surface.cpp's
        # used/remap pattern, so the output does not depend on traversal.
        used = np.unique(sub) if kept else np.zeros(0, dtype=np.int64)
        remap = np.full(len(points), -1, dtype=np.int64)
        remap[used] = np.arange(len(used), dtype=np.int64)
        out = Mesh(
            points[used] if kept else np.zeros((0, 3)),
            [("hexahedron", remap[sub])] if kept else [],
        )
    if attach_occupancy and kept > 0:
        out.cell_data["voxel:occupancy"] = [np.ones(kept, dtype=np.int64)]

    return out, {
        "dims": [int(d) for d in dims],
        "origin": [float(v) for v in origin],
        "spacing": [float(v) for v in spacing],
        "num_occupied": kept,
    }


def voxelize(
    mesh,
    resolution=None,
    cell_size=None,
    bounds=None,
    padding=0.0,
    padding_relative=0.0,
    fill="all",
    attach_occupancy=False,
    max_cells=DEFAULT_MAX_CELLS,
    sign="pseudonormal",
    watertight_check="warn",
    return_report=False,
):
    """Build a regular grid around ``mesh``.

    Exactly one of ``resolution`` and ``cell_size`` must be given: defaulting one
    would silently pick a grid the caller did not choose, and the cost is cubic in
    that choice.
    """
    if fill not in _FILLS:
        raise ValueError(
            f"{_PREFIX}unknown fill '{fill}' (expected one of: {', '.join(_FILLS)})"
        )

    out = None
    report = None
    try:
        from . import _core

        res = _core.voxelize(
            mesh,
            None if resolution is None else [int(v) for v in resolution],
            None if cell_size is None else float(cell_size),
            None if bounds is None else [float(v) for v in bounds],
            float(padding),
            float(padding_relative),
            fill,
            bool(attach_occupancy),
            int(max_cells),
            sign,
            watertight_check,
        )
        out = res["mesh"]
        report = {
            "dims": res["dims"],
            "origin": res["origin"],
            "spacing": res["spacing"],
            "num_occupied": res["num_occupied"],
        }
    except (ValueError, TypeError):
        raise
    except Exception:
        out = None

    if out is None:
        out, report = _voxelize_py(
            mesh,
            resolution,
            cell_size,
            bounds,
            padding,
            padding_relative,
            fill,
            attach_occupancy,
            max_cells,
            sign,
        )

    return (out, report) if return_report else out
