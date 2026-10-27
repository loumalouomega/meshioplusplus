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

from . import _planar_skin
from ._fallback import core_op_declined
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


def _is_pixels(mesh, resolution, cell_size):
    """2-D: a resolution whose third entry is zero, or a cell size on a polyline."""
    if resolution is not None:
        r = np.asarray(resolution, dtype=np.int64).reshape(-1)
        return r.size == 3 and int(r[2]) == 0
    return cell_size is not None and _planar_skin.polyline_only(mesh)


def _pixels_py(
    mesh,
    resolution,
    cell_size,
    bounds,
    padding,
    padding_relative,
    fill,
    attach_occupancy,
    max_cells,
    plane="auto",
):
    """The numpy twin of ``vox_voxelize_pixels``: ``quad`` pixels around a planar skin."""
    axes = _planar_skin.frame(mesh, plane)
    prims, seg_ids, has_lines = _planar_skin.collect(mesh, _PREFIX, axes)
    if fill == "inside" and has_lines:
        _planar_skin.require_closed(seg_ids, _PREFIX)

    lo0, hi0 = _resolve_bounds(mesh, bounds, 0.0, 0.0)
    lo0, hi0 = (
        _planar_skin.permute(lo0, axes).copy(),
        _planar_skin.permute(hi0, axes).copy(),
    )
    diag = float(np.sqrt(np.sum((hi0 - lo0) ** 2)))
    pad = padding + padding_relative * diag
    if pad < 0.0:
        raise ValueError(f"{_PREFIX}padding is negative")
    # A planar skin keeps its plane: padding applies to the two in-plane axes.
    lo, hi = lo0.copy(), hi0.copy()
    for k in range(2):
        lo[k] -= pad
        hi[k] += pad

    spacing = np.zeros(3)
    dims = [0, 0]
    if resolution is not None:
        if cell_size is not None:
            raise ValueError(f"{_PREFIX}give exactly one of resolution and cell_size")
        r = np.asarray(resolution, dtype=np.int64).reshape(-1)
        for k in range(2):
            dims[k] = int(r[k])
            if dims[k] <= 0:
                raise ValueError(
                    f"{_PREFIX}resolution must be positive on every axis, got "
                    f"{dims[k]} on axis {k}"
                )
            spacing[k] = (hi[k] - lo[k]) / float(dims[k])
    else:
        cell = float(cell_size)
        if not cell > 0.0:
            raise ValueError(f"{_PREFIX}cell_size must be positive")
        for k in range(2):
            extent = hi[k] - lo[k]
            dims[k] = int(np.ceil(extent / cell)) if extent > 0.0 else 0
            if dims[k] <= 0:
                raise ValueError(
                    f"{_PREFIX}the bounding box is degenerate on axis {k}, so a cell "
                    "size cannot fill it (pass an explicit resolution or bounds)"
                )
            spacing[k] = cell
    nx, ny = dims
    ncells = nx * ny
    if max_cells > 0 and ncells > max_cells:
        raise ValueError(
            f"{_PREFIX}the requested grid has {ncells} cells, above the limit of "
            f"{max_cells} (raise max_cells, coarsen the resolution, or use a band)"
        )
    origin = lo.copy()

    occupied = None
    if fill == "surface":
        occupied = np.zeros(ncells, dtype=bool)
        half = np.array([spacing[0] * 0.5, spacing[1] * 0.5, 0.0])
        for a, b, _c in prims:
            ilo, ihi, skip = [0, 0], [0, 0], False
            for k in range(2):
                tlo, thi = min(a[k], b[k]), max(a[k], b[k])
                flo = np.floor((tlo - origin[k]) / spacing[k])
                fhi = np.floor((thi - origin[k]) / spacing[k])
                flo = 0.0 if flo < 0.0 else flo
                last = float(dims[k] - 1)
                fhi = last if fhi > last else fhi
                if fhi < flo:
                    skip = True
                    break
                ilo[k], ihi[k] = int(flo), int(fhi)
            if skip:
                continue
            for j in range(ilo[1], ihi[1] + 1):
                for i in range(ilo[0], ihi[0] + 1):
                    cid = j * nx + i
                    if occupied[cid]:
                        continue
                    centre = np.array(
                        [
                            origin[0] + (float(i) + 0.5) * spacing[0],
                            origin[1] + (float(j) + 0.5) * spacing[1],
                            origin[2],
                        ]
                    )
                    if _planar_skin.seg_box(centre, half, a, b):
                        occupied[cid] = True
    elif fill == "inside":
        occupied = np.zeros(ncells, dtype=bool)
        for c in range(ncells):
            centre = np.array(
                [
                    origin[0] + (float(c % nx) + 0.5) * spacing[0],
                    origin[1] + (float(c // nx) + 0.5) * spacing[1],
                    origin[2],
                ]
            )
            occupied[c] = _planar_skin.signed_distance(prims, centre, False) < 0.0

    px = nx + 1
    cells = (
        np.arange(ncells, dtype=np.int64)
        if occupied is None
        else np.flatnonzero(occupied)
    )
    kept = int(cells.size)
    base = (cells // nx) * px + cells % nx
    nodes = np.stack([base, base + 1, base + px + 1, base + px], axis=1)
    used = np.unique(nodes) if kept else np.zeros(0, dtype=np.int64)
    remap = np.full((nx + 1) * (ny + 1), -1, dtype=np.int64)
    remap[used] = np.arange(len(used), dtype=np.int64)
    pts = np.zeros((len(used), 3))
    pts[:, 0] = origin[0] + (used % px).astype(np.float64) * spacing[0]
    pts[:, 1] = origin[1] + (used // px).astype(np.float64) * spacing[1]
    pts[:, 2] = origin[2]
    pts = _planar_skin.unpermute(pts, axes)
    out = Mesh(pts, [("quad", remap[nodes])] if kept else [])
    if attach_occupancy and kept > 0:
        out.cell_data["voxel:occupancy"] = [np.ones(kept, dtype=np.int64)]
    return out, {
        "dims": [nx, ny, 1],
        "origin": [float(v) for v in origin],
        "spacing": [float(v) for v in spacing],
        "num_occupied": kept,
    }


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
    plane="auto",
):
    if _is_pixels(mesh, resolution, cell_size):
        return _pixels_py(
            mesh,
            resolution,
            cell_size,
            bounds,
            padding,
            padding_relative,
            fill,
            attach_occupancy,
            max_cells,
            plane,
        )
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
    plane="auto",
    return_report=False,
):
    """Build a regular grid around ``mesh``.

    Exactly one of ``resolution`` and ``cell_size`` must be given: defaulting one
    would silently pick a grid the caller did not choose, and the cost is cubic in
    that choice.

    A two-entry ``resolution`` selects 2-D: ``quad`` pixels around a planar skin (a
    closed ``line`` polyline, or the boundary of a 2-D surface mesh). A
    ``cell_size`` on a mesh of nothing but ``line`` cells selects it too. ``fill``
    ``"inside"`` then uses the winding number of the polyline and refuses one that
    does not close. The report's ``dims`` is ``[nx, ny, 1]``. ``plane`` (``"auto"``,
    ``"xy"``, ``"xz"`` or ``"yz"``) chooses the skin's plane; ``"auto"`` takes the
    axis of zero extent (``z``, then ``y``, then ``x``). The 3-D variant ignores it,
    as ``octree`` does.
    """
    if fill not in _FILLS:
        raise ValueError(
            f"{_PREFIX}unknown fill '{fill}' (expected one of: {', '.join(_FILLS)})"
        )
    if plane not in _planar_skin.PLANES:
        raise ValueError(
            f"meshio++: unknown plane '{plane}' "
            f"(expected one of: {', '.join(_planar_skin.PLANES)})"
        )

    if resolution is not None and len(resolution) == 2:
        resolution = [int(resolution[0]), int(resolution[1]), 0]  # 2-D pixels

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
            plane,
        )
        out = res["mesh"]
        report = {
            "dims": res["dims"],
            "origin": res["origin"],
            "spacing": res["spacing"],
            "num_occupied": res["num_occupied"],
        }
    except Exception as exc:
        if not core_op_declined(exc, "voxelize"):
            raise
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
            plane,
        )

    return (out, report) if return_report else out
