"""Spatial subsetting of a mesh (bounding box or half-space).

A dependency-free mesh *operation* (not a file format): it extracts the part of
a mesh inside an axis-aligned bounding box or a half-space, pruning unused points
and remapping connectivity and all data. The C++ core (``_core.crop_bbox`` /
``_core.crop_plane``) does the work; this module is the thin shim (try C++, fall
back to a pure-numpy reference) plus the ``point_sets``/``cell_sets`` remapping
the core cannot see across the Python boundary.

Public API:

* :func:`crop` — crop a mesh by ``bbox=`` or ``plane=``.
"""

from __future__ import annotations

import numpy as np

from ._mesh import Mesh


# --------------------------------------------------------------------------- #
# pure-numpy reference (rectangular blocks only; ragged is C++-core only)      #
# --------------------------------------------------------------------------- #
def _points_3d(mesh):
    p = np.asarray(mesh.points, dtype=np.float64)
    if p.ndim == 1:
        p = p.reshape(len(mesh.points), -1)
    out = np.zeros((p.shape[0], 3), dtype=np.float64)
    out[:, : min(p.shape[1], 3)] = p[:, : min(p.shape[1], 3)]
    return out


def _crop_py(mesh, mask, mode, record_ids):
    n = len(mesh.points)
    all_mode = mode == "all"

    kept_per_block = []
    used = np.zeros(n, dtype=bool)
    for cb in mesh.cells:
        data = np.asarray(cb.data)
        if data.ndim != 2:
            raise NotImplementedError(
                "crop numpy fallback: ragged/polyhedron blocks are C++-core only"
            )
        inside = mask[data]
        keep = inside.all(axis=1) if all_mode else inside.any(axis=1)
        kept = np.nonzero(keep)[0]
        kept_per_block.append(kept)
        if len(kept):
            used[data[kept].reshape(-1)] = True

    new_point = np.full(n, -1, dtype=np.int64)
    point_src = np.nonzero(used)[0]
    new_point[point_src] = np.arange(len(point_src), dtype=np.int64)

    pts = np.asarray(mesh.points)
    out_points = pts[point_src] if len(point_src) else pts[:0]

    new_cells = []
    cell_maps = []
    for cb, kept in zip(mesh.cells, kept_per_block):
        data = np.asarray(cb.data)
        conn = new_point[data[kept]] if len(kept) else data[:0].astype(np.int64)
        new_cells.append((cb.type, conn.astype(np.int64)))
        cm = np.full(len(cb.data), -1, dtype=np.int64)
        if len(kept):
            cm[kept] = np.arange(len(kept), dtype=np.int64)
        cell_maps.append(cm)

    point_data = {}
    for key, value in mesh.point_data.items():
        value = np.asarray(value)
        if value.shape[:1] == (n,):
            point_data[key] = value[point_src]
        else:
            point_data[key] = value.copy()
    if record_ids:
        point_data["crop:original_point_id"] = point_src.astype(np.int64)

    cell_data = {}
    for key, blk in mesh.cell_data.items():
        out_blk = []
        for b, value in enumerate(blk):
            if value is None:
                out_blk.append(None)
                continue
            value = np.asarray(value)
            kept = kept_per_block[b]
            if value.shape[:1] == (len(mesh.cells[b].data),):
                out_blk.append(value[kept] if len(kept) else value[:0])
            else:
                out_blk.append(value.copy())
        cell_data[key] = out_blk
    if record_ids:
        cell_data["crop:original_cell_id"] = [
            kept.astype(np.int64) for kept in kept_per_block
        ]

    field_data = {
        k: (v.copy() if isinstance(v, np.ndarray) else v)
        for k, v in mesh.field_data.items()
    }
    out = Mesh(
        out_points,
        new_cells,
        point_data=point_data,
        cell_data=cell_data,
        field_data=field_data,
    )
    return out, new_point, cell_maps


# --------------------------------------------------------------------------- #
# set remapping                                                                #
# --------------------------------------------------------------------------- #
def _remap_sets(src, out, point_map, cell_maps):
    point_sets = getattr(src, "point_sets", None)
    if point_sets:
        pm = np.asarray(point_map, dtype=np.int64)
        new_ps = {}
        for name, idx in point_sets.items():
            mapped = pm[np.asarray(idx, dtype=np.int64)]
            new_ps[name] = mapped[mapped >= 0]
        out.point_sets = new_ps
    cell_sets = getattr(src, "cell_sets", None)
    if cell_sets and cell_maps is not None:
        new_cs = {}
        for name, blocks in cell_sets.items():
            remapped = []
            for b, idx in enumerate(blocks):
                if idx is None or b >= len(cell_maps):
                    remapped.append(None)
                    continue
                mapped = np.asarray(cell_maps[b], dtype=np.int64)[
                    np.asarray(idx, dtype=np.int64)
                ]
                remapped.append(mapped[mapped >= 0])
            new_cs[name] = remapped
        out.cell_sets = new_cs


# --------------------------------------------------------------------------- #
# public API                                                                   #
# --------------------------------------------------------------------------- #
def crop(
    mesh,
    *,
    bbox=None,
    plane=None,
    mode: str = "all",
    record_ids: bool = False,
) -> Mesh:
    """Extract the part of a mesh inside a bounding box or a half-space.

    Parameters
    ----------
    mesh :
        the mesh to crop (unmodified).
    bbox :
        either ``(xmin, ymin, zmin, xmax, ymax, zmax)`` or ``(lo3, hi3)``.
    plane :
        ``(point3, normal3)`` — keeps the side where ``(p - point) . normal >= 0``.
    mode :
        ``"all"`` (default) keeps a cell only if every node is inside; ``"any"``
        keeps it if any node is inside.
    record_ids :
        when true, attach ``crop:original_point_id`` point_data and
        ``crop:original_cell_id`` cell_data of the source indices.

    Returns
    -------
    Mesh
        the cropped mesh with unused points pruned. ``point_sets``/``cell_sets``
        are remapped (entries outside the kept region are dropped).

    Raises
    ------
    ValueError
        if neither/both of ``bbox`` and ``plane`` are given, or ``mode`` is
        invalid.
    """
    if (bbox is None) == (plane is None):
        raise ValueError("crop: give exactly one of bbox= or plane=")
    if mode not in ("all", "any"):
        raise ValueError(f"crop: unknown mode '{mode}' (expected 'all' or 'any')")

    if bbox is not None:
        arr = np.asarray(bbox, dtype=np.float64).reshape(-1)
        if arr.size != 6:  # accepts 6 flat numbers or (lo3, hi3), both flatten to 6
            raise ValueError("crop: bbox must be 6 numbers or (lo3, hi3)")
        lo = [float(x) for x in arr[:3]]
        hi = [float(x) for x in arr[3:]]
    else:
        point, normal = plane
        point = [float(x) for x in np.asarray(point, dtype=np.float64).reshape(3)]
        normal = [float(x) for x in np.asarray(normal, dtype=np.float64).reshape(3)]

    out = None
    point_map = None
    cell_maps = None
    try:
        from . import _core

        if bbox is not None:
            res = _core.crop_bbox(mesh, lo, hi, mode, bool(record_ids))
        else:
            res = _core.crop_plane(mesh, point, normal, mode, bool(record_ids))
        out = res["mesh"]
        point_map = np.asarray(res["point_map"])
        cell_maps = [np.asarray(a) for a in res["cell_maps"]]
    except Exception:
        out = None
    used_cpp = out is not None
    if out is None:
        p3 = _points_3d(mesh)
        if bbox is not None:
            lo3 = np.asarray(lo, dtype=np.float64)
            hi3 = np.asarray(hi, dtype=np.float64)
            mask = np.all((p3 >= lo3) & (p3 <= hi3), axis=1)
        else:
            mask = (p3 - np.asarray(point)) @ np.asarray(normal) >= 0.0
        out, point_map, cell_maps = _crop_py(mesh, mask, mode, record_ids)

    # The C++ core carries named regions across itself (and therefore
    # `point_sets`/`cell_sets`, which are views over them), so remapping again
    # here would apply the maps twice. Only the numpy fallback needs it.
    if not used_cpp:
        _remap_sets(mesh, out, point_map, cell_maps)
    return out
