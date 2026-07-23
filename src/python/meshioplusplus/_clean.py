"""Mesh cleanup (weld + prune + de-dup, one pass).

A dependency-free mesh *operation* (not a file format): it welds coincident
points (spatial-hash bucket grid, never O(N^2)), drops degenerate and
exact-duplicate cells, and removes orphaned points, remapping connectivity and
data. The C++ core (``_core.clean``) does the work; this module is the thin shim
(try C++, fall back to a pure-numpy reference) plus the ``point_sets``/
``cell_sets`` remapping the core cannot see across the Python boundary.

With no options changed, the default runs remove-orphans + drop-degenerate +
drop-duplicate-cells (welding is opt-in via ``weld=True``).

Public API:

* :func:`clean` — clean a mesh, optionally returning a report of what changed.
"""

from __future__ import annotations

import numpy as np

from ._merge import _weld_map
from ._mesh import Mesh


# --------------------------------------------------------------------------- #
# pure-numpy reference (rectangular blocks only; ragged is C++-core only)      #
# --------------------------------------------------------------------------- #
def _points_2d(mesh):
    p = np.asarray(mesh.points, dtype=np.float64)
    if p.ndim == 1:
        p = p.reshape(len(mesh.points), -1)
    return p


def _corner_count(cell_type, npc):
    linear = {
        "vertex": 1,
        "line": 2,
        "triangle": 3,
        "quad": 4,
        "tetra": 4,
        "hexahedron": 8,
        "wedge": 6,
        "pyramid": 5,
    }
    # quadratic variants share the linear parent's corner count
    for base, cc in linear.items():
        if cell_type == base or cell_type.startswith(base):
            return cc
    return npc


def _clean_py(mesh, weld, atol, remove_orphans, drop_degenerate, drop_duplicate_cells):
    pts = _points_2d(mesh)
    n = pts.shape[0]
    dim = pts.shape[1] if pts.ndim == 2 else 0

    if weld and n and dim and atol > 0:
        remap, rep_global = _weld_map(pts, dim, atol)
        W = len(rep_global)
    else:
        remap = np.arange(n, dtype=np.int64)
        rep_global = np.arange(n, dtype=np.int64)
        W = n
    points_welded = n - W

    rep_used = np.zeros(W, dtype=bool)
    report = {"degenerate": 0, "duplicate": 0}

    blocks = []  # (type, kept_conn_repspace [k,npc], kept_cells [k])
    for cb in mesh.cells:
        data = np.asarray(cb.data)
        if data.ndim != 2:
            raise NotImplementedError(
                "clean numpy fallback: ragged/polyhedron blocks are C++-core only"
            )
        npc = data.shape[1]
        cc = min(_corner_count(cb.type, npc), npc)
        seen = set()
        kept_rows = []
        kept_cells = []
        for c in range(data.shape[0]):
            row = remap[data[c].astype(np.int64)]
            degenerate = False
            if drop_degenerate:
                corners = row[:cc]
                if len(np.unique(corners)) < cc:
                    degenerate = True
            if degenerate:
                report["degenerate"] += 1
                continue
            if drop_duplicate_cells:
                key = tuple(sorted(int(x) for x in row))
                if key in seen:
                    report["duplicate"] += 1
                    continue
                seen.add(key)
            kept_rows.append(row)
            kept_cells.append(c)
            rep_used[row] = True
        conn = (
            np.array(kept_rows, dtype=np.int64)
            if kept_rows
            else np.empty((0, npc), dtype=np.int64)
        )
        blocks.append((cb.type, conn, np.array(kept_cells, dtype=np.int64)))

    # orphan compaction
    new_point = np.full(W, -1, dtype=np.int64)
    if remove_orphans:
        kept_reps = np.nonzero(rep_used)[0]
    else:
        kept_reps = np.arange(W, dtype=np.int64)
    new_point[kept_reps] = np.arange(len(kept_reps), dtype=np.int64)
    P = len(kept_reps)
    points_removed_orphan = W - P

    point_src = rep_global[kept_reps]
    out_points = pts[point_src] if P else np.empty((0, dim), dtype=pts.dtype)
    out_points = out_points.astype(pts.dtype)

    new_cells = []
    cell_maps = []
    for cb, (t, conn, kept) in zip(mesh.cells, blocks):
        remapped = new_point[conn] if len(conn) else conn
        new_cells.append((t, remapped.astype(np.int64)))
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

    cell_data = {}
    for key, blk in mesh.cell_data.items():
        out_blk = []
        for b, value in enumerate(blk):
            if value is None:
                out_blk.append(None)
                continue
            value = np.asarray(value)
            kept = blocks[b][2]
            if value.shape[:1] == (len(mesh.cells[b].data),):
                out_blk.append(value[kept] if len(kept) else value[:0])
            else:
                out_blk.append(value.copy())
        cell_data[key] = out_blk

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
    point_map = new_point[remap]
    rep = {
        "points_welded": int(points_welded),
        "points_removed_orphan": int(points_removed_orphan),
        "cells_dropped_degenerate": int(report["degenerate"]),
        "cells_dropped_duplicate": int(report["duplicate"]),
    }
    return out, point_map, cell_maps, rep


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
def clean(
    mesh,
    *,
    weld: bool = False,
    atol: float = 1e-8,
    remove_orphans: bool = True,
    drop_degenerate: bool = True,
    drop_duplicate_cells: bool = True,
    return_report: bool = False,
):
    """Clean a mesh in one pass (weld / prune / de-dup).

    Parameters
    ----------
    mesh :
        the mesh to clean (unmodified).
    weld :
        fuse coincident points within ``atol`` (keep-first). Off by default.
    atol :
        absolute coincidence tolerance for welding.
    remove_orphans :
        drop points referenced by no surviving cell (default true).
    drop_degenerate :
        drop degenerate cells — a repeated corner node or a near-zero measure
        (default true).
    drop_duplicate_cells :
        drop exact-duplicate cells, keep-first per block (default true).
    return_report :
        when true, also return a dict of removal counts.

    Returns
    -------
    Mesh, or (Mesh, report)
        the cleaned mesh. ``point_sets``/``cell_sets`` indices are remapped
        (entries pointing at removed points/cells are dropped). ``report`` has
        ``points_welded``/``points_removed_orphan``/``cells_dropped_degenerate``/
        ``cells_dropped_duplicate``.
    """
    out = None
    point_map = None
    cell_maps = None
    report = None
    try:
        from . import _core

        res = _core.clean(
            mesh,
            bool(weld),
            float(atol),
            bool(remove_orphans),
            bool(drop_degenerate),
            bool(drop_duplicate_cells),
        )
        out = res["mesh"]
        point_map = np.asarray(res["point_map"])
        cell_maps = [np.asarray(a) for a in res["cell_maps"]]
        report = {
            "points_welded": int(res["points_welded"]),
            "points_removed_orphan": int(res["points_removed_orphan"]),
            "cells_dropped_degenerate": int(res["cells_dropped_degenerate"]),
            "cells_dropped_duplicate": int(res["cells_dropped_duplicate"]),
        }
    except Exception:
        out = None
    used_cpp = out is not None
    if out is None:
        out, point_map, cell_maps, report = _clean_py(
            mesh, weld, atol, remove_orphans, drop_degenerate, drop_duplicate_cells
        )

    # The C++ core carries named regions across itself (and therefore
    # `point_sets`/`cell_sets`, which are views over them), so remapping again
    # here would apply the maps twice. Only the numpy fallback needs it.
    if not used_cpp:
        _remap_sets(mesh, out, point_map, cell_maps)

    if return_report:
        return out, report
    return out
