"""Mesh merge / weld (combine two or more meshes into one).

A dependency-free mesh *operation* (not a file format): it concatenates an
ordered list of meshes — appending points, offsetting each input's connectivity,
merging cell blocks by type, concatenating ``point_data``/``cell_data`` per the
``data_policy``, and tagging provenance — and can optionally *weld* coincident
nodes within an absolute tolerance via a spatial-hash bucket grid (never
O(N^2)). The C++ core (``_core.merge``) does the work; this module is the thin
shim (try C++, fall back to a pure-numpy reference) plus the
``point_sets``/``cell_sets`` namespacing and remapping the core cannot see
across the Python boundary.

Semantics:

* **Concatenation** (default): indices stay valid via per-input point offsets;
  same-type cell blocks are merged in input order; a per-cell ``source_mesh_id``
  cell_data records each cell's originating input (``source_tag``, on by
  default).
* **Data policy** (``data_policy``): a ``point_data``/``cell_data`` key present
  in only some inputs is either dropped (``"intersection"``, the default) or
  filled with NaN in the missing rows (``"fill"``).
* **Welding** (``weld=True``): points within ``atol`` are merged; the surviving
  point keeps the value of the *first* contributing point (keep-first). With
  ``drop_duplicate_cells`` cells that become identical after welding are removed.
* **Name collisions**: overlapping ``point_sets``/``cell_sets`` names and
  ``field_data`` keys from different inputs are namespaced by source id
  (``"0:name"``, ``"1:name"``); unique names are kept as-is.

Public API:

* :func:`merge` — merge a sequence of meshes into a new :class:`Mesh`.
"""

from __future__ import annotations

import numpy as np

from ._mesh import Mesh


# --------------------------------------------------------------------------- #
# public API                                                                   #
# --------------------------------------------------------------------------- #
def merge(
    meshes,
    *,
    weld: bool = False,
    atol: float = 1e-8,
    source_tag: bool = True,
    data_policy: str = "intersection",
    drop_duplicate_cells: bool = False,
) -> Mesh:
    """Merge two or more meshes into one.

    Parameters
    ----------
    meshes :
        an ordered sequence of two or more :class:`Mesh` objects.
    weld :
        when true, merge points closer than ``atol`` across all inputs.
    atol :
        absolute coincidence tolerance for welding (default ``1e-8``).
    source_tag :
        when true (default), attach an int64 ``source_mesh_id`` cell_data
        recording each output cell's originating input index.
    data_policy :
        ``"intersection"`` (default) keeps only data keys present in every input;
        ``"fill"`` keeps every key, filling missing rows with NaN (as float64).
    drop_duplicate_cells :
        when true (and ``weld``), drop cells that become identical after welding.

    Returns
    -------
    Mesh
        the combined mesh. ``point_sets``/``cell_sets`` indices are remapped and
        namespaced by source id on collision.

    Raises
    ------
    ValueError
        if ``meshes`` is empty or ``data_policy`` is not
        ``"intersection"``/``"fill"``.
    """
    meshes = list(meshes)
    if len(meshes) < 1:
        raise ValueError("merge: at least one input mesh is required")
    if data_policy not in ("intersection", "fill"):
        raise ValueError(
            f"merge: unknown data_policy '{data_policy}' (expected 'intersection' or 'fill')"
        )

    out = None
    point_maps = None
    cell_maps = None
    try:
        from . import _core

        res = _core.merge(
            meshes,
            bool(weld),
            float(atol),
            bool(source_tag),
            str(data_policy),
            bool(drop_duplicate_cells),
        )
        out = res["mesh"]
        point_maps = [np.asarray(a) for a in res["point_maps"]]
        cell_maps = [np.asarray(a) for a in res["cell_maps"]]
    except ValueError:
        raise
    except Exception:
        out = None
    used_cpp = out is not None
    if out is None:
        out, point_maps, cell_maps = _merge_py(
            meshes, weld, atol, source_tag, data_policy, drop_duplicate_cells
        )

    # The C++ core remaps and namespaces regions itself (and therefore
    # `point_sets`/`cell_sets`, which are views over them); doing it again here
    # would apply the maps twice. Only the numpy fallback needs it.
    if not used_cpp:
        _merge_sets(meshes, out, point_maps, cell_maps)
    return out


# --------------------------------------------------------------------------- #
# set namespacing + remapping (Python-side side-channels the core cannot see)  #
# --------------------------------------------------------------------------- #
def _merge_sets(meshes, out, point_maps, cell_maps):
    # point_sets: remap indices through the per-input point map; namespace on
    # collision across inputs.
    ps_counts = {}
    for m in meshes:
        for name in getattr(m, "point_sets", {}) or {}:
            ps_counts[name] = ps_counts.get(name, 0) + 1
    new_ps = {}
    for i, m in enumerate(meshes):
        for name, idx in (getattr(m, "point_sets", {}) or {}).items():
            key = f"{i}:{name}" if ps_counts[name] > 1 else name
            mapped = np.asarray(point_maps[i])[np.asarray(idx, dtype=np.int64)]
            new_ps[key] = np.asarray(mapped, dtype=np.int64)
    if new_ps:
        out.point_sets = new_ps

    # cell_sets: convert per-block-local -> input-global -> output-global (via
    # the cell map) -> per-output-block-local; namespace on collision.
    cs_counts = {}
    for m in meshes:
        for name in getattr(m, "cell_sets", {}) or {}:
            cs_counts[name] = cs_counts.get(name, 0) + 1
    if not cs_counts:
        return
    out_off = np.cumsum([0] + [len(cb.data) for cb in out.cells])

    def to_block_local(gcell):
        j = int(np.searchsorted(out_off, gcell, side="right") - 1)
        return j, int(gcell - out_off[j])

    new_cs = {}
    for i, m in enumerate(meshes):
        cs = getattr(m, "cell_sets", {}) or {}
        if not cs:
            continue
        in_off = np.cumsum([0] + [len(cb.data) for cb in m.cells])
        cmap = np.asarray(cell_maps[i])
        for name, blocks in cs.items():
            key = f"{i}:{name}" if cs_counts[name] > 1 else name
            acc = [[] for _ in out.cells]
            for b, idx in enumerate(blocks):
                if idx is None:
                    continue
                gbase = int(in_off[b])
                for local in np.asarray(idx, dtype=np.int64):
                    og = int(cmap[gbase + int(local)])
                    if og < 0:
                        continue
                    j, jl = to_block_local(og)
                    acc[j].append(jl)
            new_cs[key] = [np.array(sorted(a), dtype=np.int64) for a in acc]
    if new_cs:
        out.cell_sets = new_cs


# --------------------------------------------------------------------------- #
# pure-numpy reference (used when the C++ core is unavailable)                 #
# --------------------------------------------------------------------------- #
def _weld_map(pts, dim, atol):
    """Spatial-hash weld: return (remap global->new, rep_global new->global).

    Mirrors the C++ core: quantize to atol-sized buckets, iterate points in
    global order, and weld a point onto the first earlier representative within
    ``atol`` found in the 3x3x3 bucket neighbourhood (keep-first).
    """
    total = pts.shape[0]
    ddim = min(dim, 3)
    keys = (
        np.floor(pts[:, :ddim] / atol).astype(np.int64)
        if ddim > 0
        else np.zeros((total, 0), dtype=np.int64)
    )
    grid = {}
    remap = np.empty(total, dtype=np.int64)
    rep_global = []
    atol2 = atol * atol
    for g in range(total):
        k = tuple(int(keys[g, d]) for d in range(ddim))
        kk = k + (0,) * (3 - len(k))
        found = -1
        for dz in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dx in (-1, 0, 1):
                    lst = grid.get((kk[0] + dx, kk[1] + dy, kk[2] + dz))
                    if not lst:
                        continue
                    for rep in lst:
                        rg = rep_global[rep]
                        d2 = float(np.sum((pts[g, :ddim] - pts[rg, :ddim]) ** 2))
                        if d2 <= atol2:
                            found = rep
                            break
                    if found >= 0:
                        break
                if found >= 0:
                    break
            if found >= 0:
                break
        if found >= 0:
            remap[g] = found
        else:
            nidx = len(rep_global)
            remap[g] = nidx
            rep_global.append(g)
            grid.setdefault(kk, []).append(nidx)
    return remap, np.array(rep_global, dtype=np.int64)


def _points_2d(mesh):
    p = np.asarray(mesh.points, dtype=np.float64)
    if p.ndim == 1:
        p = p.reshape(len(mesh.points), -1)
    return p


def _merge_py(meshes, weld, atol, source_tag, data_policy, drop_duplicate_cells):
    nps = [len(m.points) for m in meshes]
    poff = np.cumsum([0] + nps)
    total = int(poff[-1])
    dim = 0
    for m in meshes:
        p = _points_2d(m)
        dim = max(dim, p.shape[1] if p.size else 0)

    concat = np.zeros((total, dim), dtype=np.float64)
    for i, m in enumerate(meshes):
        if nps[i] == 0:
            continue
        p = _points_2d(m)
        pdim = p.shape[1] if p.size else 0
        if pdim:
            concat[poff[i] : poff[i + 1], :pdim] = p[:, :pdim]

    do_weld = bool(weld) and total > 0 and dim > 0 and atol > 0
    if do_weld:
        remap, rep_global = _weld_map(concat, dim, atol)
        out_points = concat[rep_global]
    else:
        remap = np.arange(total, dtype=np.int64)
        rep_global = np.arange(total, dtype=np.int64)
        out_points = concat
    point_maps = [remap[poff[i] : poff[i + 1]].copy() for i in range(len(meshes))]

    # group cell blocks by type (first-seen order); record contributions.
    out_types = []
    type_index = {}
    contribs = []  # per out block: list of (mesh, block, ncells, start)
    npc_list = []
    in_off = [np.cumsum([0] + [len(cb.data) for cb in m.cells]) for m in meshes]
    for i, m in enumerate(meshes):
        for b, cb in enumerate(m.cells):
            data = cb.data
            if not (isinstance(data, np.ndarray) and data.ndim == 2):
                raise NotImplementedError(
                    "merge numpy fallback: ragged/polyhedron cell blocks are "
                    "only supported by the C++ core"
                )
            t = cb.type
            if t not in type_index:
                type_index[t] = len(out_types)
                out_types.append(t)
                contribs.append([])
                npc_list.append(data.shape[1])
            oi = type_index[t]
            start = sum(c[2] for c in contribs[oi])
            contribs[oi].append((i, b, data.shape[0], start))

    # build remapped connectivity + source ids per out block (pre-dedup).
    built_conn = []
    built_src = []
    for oi, _ in enumerate(out_types):
        npc = npc_list[oi]
        pre = sum(c[2] for c in contribs[oi])
        conn = np.empty((pre, npc), dtype=np.int64)
        src = np.empty(pre, dtype=np.int64)
        for i, b, ncells, start in contribs[oi]:
            data = np.asarray(meshes[i].cells[b].data, dtype=np.int64)
            conn[start : start + ncells] = remap[poff[i] + data]
            src[start : start + ncells] = i
        built_conn.append(conn)
        built_src.append(src)

    # dedup (keep-first) if requested.
    finalpos = []
    final_to_pre = []
    for oi in range(len(out_types)):
        conn = built_conn[oi]
        pre = conn.shape[0]
        if weld and drop_duplicate_cells:
            seen = {}
            fp = np.full(pre, -1, dtype=np.int64)
            ftp = []
            for p in range(pre):
                key = tuple(sorted(int(x) for x in conn[p]))
                if key in seen:
                    fp[p] = -1
                else:
                    fp[p] = len(ftp)
                    seen[key] = fp[p]
                    ftp.append(p)
            finalpos.append(fp)
            final_to_pre.append(np.array(ftp, dtype=np.int64))
        else:
            finalpos.append(np.arange(pre, dtype=np.int64))
            final_to_pre.append(np.arange(pre, dtype=np.int64))

    new_cells = []
    out_block_sizes = []
    for oi, t in enumerate(out_types):
        ftp = final_to_pre[oi]
        new_cells.append((t, built_conn[oi][ftp]))
        out_block_sizes.append(len(ftp))
    out_off = np.cumsum([0] + out_block_sizes)

    point_data = _merge_point_data(
        meshes, nps, poff, total, rep_global, do_weld, data_policy
    )
    cell_data = _merge_cell_data(
        meshes, out_types, contribs, built_conn, final_to_pre, data_policy
    )
    if source_tag and out_types:
        cell_data["source_mesh_id"] = [
            built_src[oi][final_to_pre[oi]].astype(np.int64)
            for oi in range(len(out_types))
        ]

    # field_data: namespace on collision.
    field_data = {}
    fcounts = {}
    for m in meshes:
        for name in m.field_data:
            fcounts[name] = fcounts.get(name, 0) + 1
    for i, m in enumerate(meshes):
        for name, val in m.field_data.items():
            key = f"{i}:{name}" if fcounts[name] > 1 else name
            field_data[key] = val.copy() if isinstance(val, np.ndarray) else val

    # per-input cell maps (input global cell -> output global cell, -1 dropped).
    cell_maps = [
        np.full(int(in_off[i][-1]), -1, dtype=np.int64) for i in range(len(meshes))
    ]
    for oi in range(len(out_types)):
        fp = finalpos[oi]
        for i, b, ncells, start in contribs[oi]:
            gbase = int(in_off[i][b])
            for c in range(ncells):
                fl = int(fp[start + c])
                cell_maps[i][gbase + c] = -1 if fl < 0 else int(out_off[oi]) + fl

    out = Mesh(
        out_points,
        new_cells,
        point_data=point_data,
        cell_data=cell_data,
        field_data=field_data,
    )
    return out, point_maps, cell_maps


def _resolve_dtype(present, any_missing):
    if any_missing:
        return np.float64
    dts = {a.dtype for a in present}
    return present[0].dtype if len(dts) == 1 else np.float64


def _merge_point_data(meshes, nps, poff, total, rep_global, do_weld, data_policy):
    point_data = {}
    names = set()
    for m in meshes:
        names |= set(m.point_data.keys())
    for name in sorted(names):
        arrays = []
        present = []
        any_missing = False
        cols = None
        ok = True
        for i, m in enumerate(meshes):
            a = m.point_data.get(name)
            if a is not None and len(a) == nps[i]:
                a = np.asarray(a)
                arrays.append(a)
                present.append(a)
                c = a.shape[1] if a.ndim == 2 else 1
                if cols is None:
                    cols = c
                elif c != cols:
                    ok = False
            else:
                arrays.append(None)
                any_missing = True
                if a is not None:
                    ok = False
        if cols is None or not ok:
            continue
        if any_missing and data_policy == "intersection":
            continue
        dt = _resolve_dtype(present, any_missing)
        col_is_1d = all(a.ndim == 1 for a in present)
        shape = (total,) if col_is_1d else (total, cols)
        concat = np.full(shape, np.nan if any_missing else 0, dtype=dt)
        for i, a in enumerate(arrays):
            sl = slice(poff[i], poff[i + 1])
            if a is None:
                concat[sl] = np.nan
            else:
                concat[sl] = a if col_is_1d else a.reshape(nps[i], cols)
        point_data[name] = concat[rep_global] if do_weld else concat
    return point_data


def _merge_cell_data(
    meshes, out_types, contribs, built_conn, final_to_pre, data_policy
):
    cell_data = {}
    names = set()
    for m in meshes:
        names |= set(m.cell_data.keys())
    for name in sorted(names):
        any_missing = False
        cols = None
        ok = True
        present = []
        for oi in range(len(out_types)):
            for i, b, ncells, start in contribs[oi]:
                cd = meshes[i].cell_data.get(name)
                if (
                    cd is not None
                    and b < len(cd)
                    and cd[b] is not None
                    and len(cd[b]) == ncells
                ):
                    a = np.asarray(cd[b])
                    present.append(a)
                    c = a.shape[1] if a.ndim == 2 else 1
                    if cols is None:
                        cols = c
                    elif c != cols:
                        ok = False
                else:
                    any_missing = True
        if cols is None or not ok:
            continue
        if any_missing and data_policy == "intersection":
            continue
        dt = _resolve_dtype(present, any_missing)
        col_is_1d = all(a.ndim == 1 for a in present)
        blocks = []
        for oi in range(len(out_types)):
            pre = built_conn[oi].shape[0]
            shape = (pre,) if col_is_1d else (pre, cols)
            arr = np.full(shape, np.nan if any_missing else 0, dtype=dt)
            for i, b, ncells, start in contribs[oi]:
                cd = meshes[i].cell_data.get(name)
                if (
                    cd is not None
                    and b < len(cd)
                    and cd[b] is not None
                    and len(cd[b]) == ncells
                ):
                    a = np.asarray(cd[b])
                    arr[start : start + ncells] = (
                        a if col_is_1d else a.reshape(ncells, cols)
                    )
                else:
                    arr[start : start + ncells] = np.nan
            blocks.append(arr[final_to_pre[oi]])
        cell_data[name] = blocks
    return cell_data
