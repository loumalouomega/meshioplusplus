"""Subsetting a mesh: by bounding box, by half-space, or by a data predicate.

A dependency-free mesh *operation* (not a file format): it extracts the part of
a mesh inside an axis-aligned bounding box or a half-space, or the cells whose
value in one of its own ``cell_data`` arrays satisfies a comparison, pruning
unused points and remapping connectivity and all data. The C++ core
(``_core.crop_bbox`` / ``_core.crop_plane`` / ``_core.crop_predicate``) does the
work; this module is the thin shim (try C++, fall back to a pure-numpy
reference) plus the ``point_sets``/``cell_sets`` remapping the core cannot see
across the Python boundary.

The predicate form is deliberately **general rather than
inside/outside-a-surface specific** -- inside/outside then composes::

    f = mp.distance_to_surface(mesh, skin, location="center")
    inside = mp.crop(f, where=("sdf:distance", "<", 0.0))

and the same one mode also crops by ``quality:*``, by a material id, by
``partition:part``, or by anything ``data_calc`` can produce.

``mode`` is **not** accepted with ``where=``, and its absence is the honest
answer rather than an omission: ``bbox``/``plane`` test *points* and then need a
rule for reducing a cell's several nodes to one verdict, whereas a ``cell_data``
predicate is already one value per cell and has nothing to reduce.

Public API:

* :func:`crop` — crop a mesh by ``bbox=``, ``plane=`` or ``where=``.
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


#: The comparison vocabulary, shared with ``refine``'s predicate selector. The
#: two operations evaluate it through the same C++ function, so a second
#: transcription here would be the thing that drifts.
_COMPARES = ("<", "<=", ">", ">=", "==", "!=")


def _compare_py(values, compare, value):
    """``refine_compare_value``, vectorized. A non-finite value never matches."""
    values = np.asarray(values, dtype=np.float64).reshape(-1)
    finite = np.isfinite(values)
    if compare == "<":
        hit = values < value
    elif compare == "<=":
        hit = values <= value
    elif compare == ">":
        hit = values > value
    elif compare == ">=":
        hit = values >= value
    elif compare == "==":
        hit = values == value
    else:
        hit = values != value
    return finite & hit


def _crop_predicate_py(mesh, array, compare, value, record_ids):
    if array not in mesh.cell_data:
        if array in mesh.point_data:
            raise ValueError(
                f"crop: '{array}' is point_data; a crop predicate is per cell, so "
                "convert it first (`data to-cell`)"
            )
        available = ", ".join(sorted(mesh.cell_data)) or "(none)"
        raise ValueError(
            f"crop: no cell_data array named '{array}' (available: {available})"
        )
    blocks = mesh.cell_data[array]
    if len(blocks) != len(mesh.cells):
        raise ValueError(f"crop: cell_data '{array}' does not cover every cell block")

    kept_per_block = []
    for cb, blk in zip(mesh.cells, blocks):
        arr = np.asarray(blk)
        nc = len(cb.data)
        if arr.shape[0] != nc:
            raise ValueError(
                f"crop: cell_data '{array}' has {arr.shape[0]} rows on a block of "
                f"{nc} cells"
            )
        if nc and arr.ndim > 1 and int(np.prod(arr.shape[1:])) != 1:
            raise ValueError(
                f"crop: cell_data '{array}' must be scalar (one value per cell) to be "
                "used as a crop predicate"
            )
        kept_per_block.append(np.nonzero(_compare_py(arr, compare, value))[0])
    return _subset_py(mesh, kept_per_block, record_ids)


def _crop_py(mesh, mask, mode, record_ids):
    all_mode = mode == "all"
    kept_per_block = []
    for cb in mesh.cells:
        data = np.asarray(cb.data)
        if data.ndim != 2:
            raise NotImplementedError(
                "crop numpy fallback: ragged/polyhedron blocks are C++-core only"
            )
        inside = mask[data]
        keep = inside.all(axis=1) if all_mode else inside.any(axis=1)
        kept_per_block.append(np.nonzero(keep)[0])
    return _subset_py(mesh, kept_per_block, record_ids)


def _subset_py(mesh, kept_per_block, record_ids):
    """Prune points and remap everything, given the kept cells per block.

    The numpy counterpart of ``crop_finish``, which all three crop modes funnel
    through for the same reason: the pruning and remapping must not depend on how
    the kept set was chosen.
    """
    n = len(mesh.points)
    used = np.zeros(n, dtype=bool)
    for cb, kept in zip(mesh.cells, kept_per_block):
        data = np.asarray(cb.data)
        if data.ndim != 2:
            raise NotImplementedError(
                "crop numpy fallback: ragged/polyhedron blocks are C++-core only"
            )
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
    where=None,
    mode: str = "all",
    record_ids: bool = False,
) -> Mesh:
    """Extract part of a mesh, by geometry or by a ``cell_data`` predicate.

    Parameters
    ----------
    mesh :
        the mesh to crop (unmodified).
    bbox :
        either ``(xmin, ymin, zmin, xmax, ymax, zmax)`` or ``(lo3, hi3)``.
    plane :
        ``(point3, normal3)`` — keeps the side where ``(p - point) . normal >= 0``.
    where :
        ``(array_name, comparison, value)`` — keeps cells whose value in the
        scalar ``cell_data`` array satisfies the comparison, one of ``"<"``,
        ``"<="``, ``">"``, ``">="``, ``"=="``, ``"!="``. **A non-finite cell
        value never matches**, whatever the comparison: ``compute_quality``
        reports NaN where a metric does not apply, and predicating over
        ``quality:*`` on a mixed mesh is the headline use case.
    mode :
        ``"all"`` (default) keeps a cell only if every node is inside; ``"any"``
        keeps it if any node is inside. Applies to ``bbox``/``plane`` only;
        passing a non-default value with ``where=`` is an error rather than
        silently ignored, since it would mean nothing.
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
        if not exactly one of ``bbox``, ``plane`` and ``where`` is given, if
        ``mode`` is invalid or is set alongside ``where``, or if ``where`` names
        an array that is not scalar ``cell_data`` covering every block.
    """
    given = [bbox is not None, plane is not None, where is not None]
    if sum(given) != 1:
        raise ValueError("crop: give exactly one of bbox=, plane= or where=")
    if mode not in ("all", "any"):
        raise ValueError(f"crop: unknown mode '{mode}' (expected 'all' or 'any')")
    if where is not None and mode != "all":
        raise ValueError(
            "crop: mode= applies to bbox= and plane=, which test points; a where= "
            "predicate is already one value per cell and has nothing to reduce"
        )

    if where is not None:
        try:
            array, compare, value = where
        except (TypeError, ValueError):
            raise ValueError(
                "crop: where= must be (array_name, comparison, value)"
            ) from None
        array = str(array)
        compare = str(compare)
        if compare not in _COMPARES:
            raise ValueError(
                f"crop: unknown comparison '{compare}' "
                f"(expected one of: {', '.join(_COMPARES)})"
            )
        value = float(value)
    elif bbox is not None:
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
        elif plane is not None:
            res = _core.crop_plane(mesh, point, normal, mode, bool(record_ids))
        else:
            res = _core.crop_predicate(mesh, array, compare, value, bool(record_ids))
        out = res["mesh"]
        point_map = np.asarray(res["point_map"])
        cell_maps = [np.asarray(a) for a in res["cell_maps"]]
    except Exception:
        out = None
    used_cpp = out is not None
    if out is None:
        if where is not None:
            out, point_map, cell_maps = _crop_predicate_py(
                mesh, array, compare, value, record_ids
            )
        else:
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
