"""Partition a mesh into several submeshes (by type / region / component).

A dependency-free mesh *operation* (not a file format): it splits one mesh into
several, each pruned to its own points with connectivity and data remapped. The
C++ core (``_core.split``) does the work by cell type, connected component, or an
integer ``cell_data`` tag; this shim adds a Python-only ``by="region"`` path that
partitions by named ``cell_sets`` when the mesh carries them, plus the
``point_sets``/``cell_sets`` remapping the core cannot see across the boundary.

Public API:

* :func:`split` — split a mesh; returns an ordered ``dict`` ``{key: Mesh}``.
"""

from __future__ import annotations

import numpy as np

from ._mesh import Mesh


# --------------------------------------------------------------------------- #
# numpy submesh builder (rectangular blocks; drops now-empty blocks)          #
# --------------------------------------------------------------------------- #
def _build_subset_py(mesh, kept_per_block):
    n = len(mesh.points)
    used = np.zeros(n, dtype=bool)
    for cb, kept in zip(mesh.cells, kept_per_block):
        data = np.asarray(cb.data)
        if data.ndim != 2:
            raise NotImplementedError(
                "split numpy fallback: ragged/polyhedron blocks are C++-core only"
            )
        if len(kept):
            used[data[np.asarray(kept, dtype=np.int64)].reshape(-1)] = True
    new_point = np.full(n, -1, dtype=np.int64)
    point_src = np.nonzero(used)[0]
    new_point[point_src] = np.arange(len(point_src), dtype=np.int64)

    pts = np.asarray(mesh.points)
    out_points = pts[point_src] if len(point_src) else pts[:0]

    new_cells = []
    cell_maps = []
    for cb, kept in zip(mesh.cells, kept_per_block):
        kept = np.asarray(kept, dtype=np.int64)
        cm = np.full(len(cb.data), -1, dtype=np.int64)
        if len(kept):
            cm[kept] = np.arange(len(kept), dtype=np.int64)
            conn = new_point[np.asarray(cb.data)[kept]]
            new_cells.append((cb.type, conn.astype(np.int64)))
        cell_maps.append(cm)

    point_data = {}
    for key, value in mesh.point_data.items():
        value = np.asarray(value)
        point_data[key] = value[point_src] if value.shape[:1] == (n,) else value.copy()

    # cell_data: only for emitted (non-empty) blocks, in emit order
    cell_data = {}
    for key, blk in mesh.cell_data.items():
        out_blk = []
        for b, (cb, kept) in enumerate(zip(mesh.cells, kept_per_block)):
            if not len(kept):
                continue
            value = blk[b] if b < len(blk) else None
            if value is None:
                out_blk.append(None)
            else:
                value = np.asarray(value)
                out_blk.append(value[np.asarray(kept, dtype=np.int64)])
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
    return out, new_point, cell_maps


def _groups_by_type(mesh):
    groups = {}
    order = []
    for b, cb in enumerate(mesh.cells):
        if cb.type not in groups:
            groups[cb.type] = [[] for _ in mesh.cells]
            order.append(cb.type)
        groups[cb.type][b] = list(range(len(cb.data)))
    return [(k, groups[k]) for k in order]


def _groups_by_component(mesh):
    nblocks = len(mesh.cells)
    base = np.cumsum([0] + [len(cb.data) for cb in mesh.cells])
    total = int(base[-1])
    parent = np.arange(total, dtype=np.int64)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    node_first = np.full(len(mesh.points), -1, dtype=np.int64)
    for b, cb in enumerate(mesh.cells):
        data = np.asarray(cb.data)
        for c in range(data.shape[0]):
            gc = int(base[b] + c)
            for node in data[c]:
                node = int(node)
                if node_first[node] < 0:
                    node_first[node] = gc
                else:
                    union(int(node_first[node]), gc)

    root_index = {}
    order = []
    gc_group = np.zeros(total, dtype=np.int64)
    for gc in range(total):
        r = find(gc)
        if r not in root_index:
            root_index[r] = len(order)
            order.append([[] for _ in range(nblocks)])
        gc_group[gc] = root_index[r]

    for b in range(nblocks):
        for c in range(int(base[b + 1] - base[b])):
            order[int(gc_group[int(base[b] + c)])][b].append(c)
    return [(str(i), g) for i, g in enumerate(order)]


def _groups_by_tag(mesh, tag):
    nblocks = len(mesh.cells)
    name = tag
    if not name:
        for k, blk in mesh.cell_data.items():
            if len(blk) == nblocks and all(
                b is not None
                and np.asarray(b).dtype.kind in "iu"
                and len(b) == len(cb.data)
                for b, cb in zip(blk, mesh.cells)
            ):
                name = k
                break
    if not name or name not in mesh.cell_data:
        raise ValueError(
            "split: no integer cell_data tag found (pass an explicit name)"
        )
    blk = mesh.cell_data[name]
    groups = {}
    order = []
    for b, cb in enumerate(mesh.cells):
        vals = np.asarray(blk[b])
        for c in range(len(cb.data)):
            v = int(vals[c])
            if v not in groups:
                groups[v] = [[] for _ in mesh.cells]
                order.append(v)
            groups[v][b].append(c)
    return [(str(v), groups[v]) for v in order]


def _groups_by_region(mesh):
    """Pure-Python twin of the C++ ``split_by_region``: one group per named
    Cell region, independent of the others (a cell may land in several).

    Fallback only -- the C++ path (criterion ``"regions"``) is tried first;
    this runs when ``_core`` is unavailable.
    """
    from ._regions import block_bases

    bases = block_bases(mesh.cells)
    out = []
    for r in getattr(mesh, "regions", []):
        if r.kind != "cell":
            continue
        kept = [[] for _ in mesh.cells]
        for g in np.asarray(r.entries).reshape(-1):
            b = int(np.searchsorted(bases, g, side="right")) - 1
            if 0 <= b < len(mesh.cells):
                kept[b].append(int(g - bases[b]))
        out.append((r.name, kept))
    return out


def _split_py(mesh, by, tag):
    if by == "type":
        grouped = _groups_by_type(mesh)
    elif by == "component":
        grouped = _groups_by_component(mesh)
    elif by in ("region", "tag"):
        grouped = _groups_by_tag(mesh, tag)
    elif by == "regions":
        grouped = _groups_by_region(mesh)
    else:
        raise ValueError(f"split: unknown criterion '{by}'")
    pieces = []
    for key, kept in grouped:
        out, pm, cm = _build_subset_py(mesh, kept)
        pieces.append((key, out, pm, cm))
    return pieces


# --------------------------------------------------------------------------- #
# set remapping (per piece)                                                    #
# --------------------------------------------------------------------------- #
def _cell_set_is_remappable(mesh, blocks):
    """Whether a ``cell_sets`` entry is per-block cell-index lists.

    Formats can stash other per-block side data in ``cell_sets`` — gmsh's
    ``gmsh:bounding_entities`` holds entity tags (negative values included),
    not cell indices. Split reorders/drops blocks, so such entries cannot be
    carried meaningfully and are skipped rather than crashed on.
    """
    for b, cb in enumerate(mesh.cells):
        if b >= len(blocks) or blocks[b] is None:
            continue
        idx = np.asarray(blocks[b]).reshape(-1)
        if len(idx) == 0:
            continue
        if idx.dtype.kind not in "iu":
            return False
        if idx.min() < 0 or idx.max() >= len(cb.data):
            return False
    return True


def _remap_piece_sets(src, out, point_map, cell_maps):
    point_sets = getattr(src, "point_sets", None)
    if point_sets:
        pm = np.asarray(point_map, dtype=np.int64)
        new_ps = {}
        for name, idx in point_sets.items():
            idx = np.asarray(idx)
            if idx.dtype.kind not in "iu" or (
                len(idx) and (idx.min() < 0 or idx.max() >= len(pm))
            ):
                continue  # not point indices -- cannot be remapped
            mapped = pm[idx.astype(np.int64)]
            new_ps[name] = mapped[mapped >= 0]
        if new_ps:
            out.point_sets = new_ps
    cell_sets = getattr(src, "cell_sets", None)
    if cell_sets and cell_maps is not None:
        emit_blocks = [
            b for b in range(len(cell_maps)) if np.any(np.asarray(cell_maps[b]) >= 0)
        ]
        new_cs = {}
        for name, blocks in cell_sets.items():
            if not _cell_set_is_remappable(src, blocks):
                continue  # side data, not cell indices -- cannot be carried
            out_blocks = []
            for b in emit_blocks:
                if b < len(blocks) and blocks[b] is not None:
                    mapped = np.asarray(cell_maps[b], dtype=np.int64)[
                        np.asarray(blocks[b], dtype=np.int64).reshape(-1)
                    ]
                    out_blocks.append(mapped[mapped >= 0])
                else:
                    out_blocks.append(np.empty(0, dtype=np.int64))
            new_cs[name] = out_blocks
        if new_cs:
            out.cell_sets = new_cs


def _split_by_cell_sets(mesh):
    """Python-only region split by named cell_sets (one piece per set)."""
    result = {}
    for name, blocks in mesh.cell_sets.items():
        kept = []
        for b, cb in enumerate(mesh.cells):
            if b < len(blocks) and blocks[b] is not None:
                kept.append([int(x) for x in np.asarray(blocks[b], dtype=np.int64)])
            else:
                kept.append([])
        out, pm, cm = _build_subset_py(mesh, kept)
        _remap_piece_sets(mesh, out, pm, cm)
        result[name] = out
    return result


# --------------------------------------------------------------------------- #
# public API                                                                   #
# --------------------------------------------------------------------------- #
def split(mesh, by: str = "type", tag: str | None = None) -> dict:
    """Split a mesh into pieces.

    Parameters
    ----------
    mesh :
        the mesh to split.
    by :
        ``"type"`` (one piece per cell type), ``"component"`` (one per connected
        component), ``"region"`` (by named ``cell_sets`` if present, else by the
        integer ``cell_data`` tag ``tag``), or ``"regions"`` (plural -- one piece
        per named **Cell** region, cross-binding-capable since it runs in the
        C++ core rather than this shim; unlike every other criterion it is not
        a partition, since regions may overlap: a cell can land in zero, one or
        several pieces).
    tag :
        for ``by="region"`` without ``cell_sets``, the integer ``cell_data`` name
        to split on (default: the first integer cell_data). Unused by
        ``by="regions"``.

    Returns
    -------
    dict
        an ordered ``{key: Mesh}`` mapping. Keys are type names / component
        indices / set names / tag values / region names. ``point_sets``/
        ``cell_sets`` are remapped into each piece.

    Raises
    ------
    ValueError
        for an unknown ``by`` or when ``by="region"`` finds no tag/sets.
    """
    if by == "region" and getattr(mesh, "cell_sets", None):
        return _split_by_cell_sets(mesh)

    by_core = "tag" if by == "region" else by
    if by_core not in ("type", "component", "tag", "regions"):
        raise ValueError(
            f"split: unknown criterion '{by}' "
            "(expected 'type', 'region', 'regions', or 'component')"
        )

    pieces = None
    try:
        from . import _core

        raw = _core.split(mesh, by_core, tag or "")
        pieces = [
            (
                p["key"],
                p["mesh"],
                np.asarray(p["point_map"]),
                [np.asarray(c) for c in p["cell_maps"]],
            )
            for p in raw
        ]
    except ValueError:
        raise
    except Exception:
        pieces = None
    used_cpp = pieces is not None
    if pieces is None:
        pieces = _split_py(mesh, by_core, tag or "")

    result = {}
    for key, out, pm, cm in pieces:
        # The C++ core carries named regions onto each piece itself (and
        # therefore `point_sets`/`cell_sets`, which are views over them);
        # remapping again here would apply the maps twice.
        if not used_cpp:
            _remap_piece_sets(mesh, out, pm, cm)
        result[key] = out
    return result
