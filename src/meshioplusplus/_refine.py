"""Uniform mesh refinement: subdivide every cell into same-type children.

A dependency-free mesh *operation* (not a file format): it increases a mesh's
resolution while leaving the object it describes intact. The C++ core
(``_core.refine``) does the work; this module is the thin shim (try C++, fall
back to a pure-numpy reference) plus the ``point_sets``/``cell_sets`` remapping
the core cannot see across the Python boundary.

One level applies a fixed template per cell type: ``line`` -> 2, ``triangle``
-> 4, ``quad`` -> 4, ``tetra`` -> 8, ``wedge`` -> 8, ``hexahedron`` -> 8. New
nodes sit at the midpoints of the parent's edges, quad faces and (hexahedron
only) body, and carry the mean of that entity's corner values for every
``point_data`` array. Mid-edge and quad-face-centre nodes are **shared** between
every cell touching the entity, so the refined mesh has no hanging nodes.

Higher-order cells, pyramids, and ragged polygon/polyhedron blocks have no
same-type subdivision and raise :class:`ValueError`.

Public API:

* :func:`refine` — uniformly refine a mesh.
"""

from __future__ import annotations

import numpy as np

from ._mesh import Mesh

# --------------------------------------------------------------------------- #
# templates (mirroring cpp/src/operations/refine.cpp)                          #
# --------------------------------------------------------------------------- #
# Local node ids address one flat space per parent cell: corners first, then one
# node per edge (in _EDGES order), then one per quad face (in _QUAD_FACES
# order), then the body node when the type has one. That layout coincides with
# each type's own full-Lagrange numbering (line3 / triangle6 / quad9 / tetra10 /
# wedge18 / hexahedron27).
_EDGES = {
    "line": [(0, 1)],
    "triangle": [(0, 1), (1, 2), (2, 0)],
    "quad": [(0, 1), (1, 2), (2, 3), (3, 0)],
    "tetra": [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)],
    "hexahedron": [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 4),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ],
    "wedge": [
        (0, 1),
        (1, 2),
        (2, 0),
        (3, 4),
        (4, 5),
        (5, 3),
        (0, 3),
        (1, 4),
        (2, 5),
    ],
}

_QUAD_FACES = {
    "line": [],
    "triangle": [],
    "tetra": [],
    "quad": [(0, 1, 2, 3)],
    "hexahedron": [
        (0, 1, 5, 4),
        (1, 2, 6, 5),
        (2, 3, 7, 6),
        (3, 0, 4, 7),
        (0, 1, 2, 3),
        (4, 5, 6, 7),
    ],
    "wedge": [(0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)],
}

_HAS_BODY = {"hexahedron"}

_NUM_CORNERS = {
    "line": 2,
    "triangle": 3,
    "quad": 4,
    "tetra": 4,
    "wedge": 6,
    "hexahedron": 8,
}

_CHILDREN = {
    "line": [(0, 2), (2, 1)],
    "triangle": [(0, 3, 5), (3, 1, 4), (5, 4, 2), (3, 4, 5)],
    "quad": [(0, 4, 8, 7), (4, 1, 5, 8), (8, 5, 2, 6), (7, 8, 6, 3)],
    # Four corner tetrahedra, then the residual octahedron split along the fixed
    # interior diagonal 4-9 with the remaining ring 6->7->8->5.
    "tetra": [
        (0, 4, 6, 7),
        (4, 1, 5, 8),
        (6, 5, 2, 9),
        (7, 8, 9, 3),
        (4, 9, 6, 7),
        (4, 9, 7, 8),
        (4, 9, 8, 5),
        (4, 9, 5, 6),
    ],
    "wedge": [
        (0, 6, 8, 12, 15, 17),
        (6, 1, 7, 15, 13, 16),
        (8, 7, 2, 17, 16, 14),
        (6, 7, 8, 15, 16, 17),
        (12, 15, 17, 3, 9, 11),
        (15, 13, 16, 9, 4, 10),
        (17, 16, 14, 11, 10, 5),
        (15, 16, 17, 9, 10, 11),
    ],
    "hexahedron": [
        (0, 8, 24, 11, 16, 20, 26, 23),
        (8, 1, 9, 24, 20, 17, 21, 26),
        (11, 24, 10, 3, 23, 26, 22, 19),
        (24, 9, 2, 10, 26, 21, 18, 22),
        (16, 20, 26, 23, 4, 12, 25, 15),
        (20, 17, 21, 26, 12, 5, 13, 25),
        (23, 26, 22, 19, 15, 25, 14, 7),
        (26, 21, 18, 22, 25, 13, 6, 14),
    ],
}

SUPPORTED_TYPES = tuple(_CHILDREN)


# --------------------------------------------------------------------------- #
# pure-numpy reference (rectangular blocks only; ragged is C++-core only)      #
# --------------------------------------------------------------------------- #
def _check_type(cell_type, data):
    if data.ndim != 2:
        raise ValueError(
            f"refine: cannot refine ragged cell block '{cell_type}' "
            "(polygon/polyhedron blocks have no same-type subdivision template)"
        )
    if cell_type in _CHILDREN:
        return
    if cell_type == "pyramid":
        raise ValueError(
            "refine: cannot refine cell type 'pyramid' into same-type children "
            "(its uniform refinement is 6 pyramids + 4 tetrahedra; simplexify "
            "the mesh first)"
        )
    raise ValueError(
        f"refine: cannot refine cell type '{cell_type}' into same-type children "
        "(higher-order cells have no same-type subdivision template; linearize "
        "the mesh first)"
    )


def _refine_once_py(mesh):
    n = len(mesh.points)
    blocks = [(cb.type, np.asarray(cb.data)) for cb in mesh.cells]
    for cell_type, data in blocks:
        _check_type(cell_type, data)

    # --- entity keys, in the same (block, cell, slot) order the C++ core uses,
    # so the numbering matches it exactly.
    keys = []
    for cell_type, data in blocks:
        edges = _EDGES[cell_type]
        faces = _QUAD_FACES[cell_type]
        for row in data:
            for a, b in edges:
                x, y = int(row[a]), int(row[b])
                keys.append((-1, -1, x, y) if x < y else (-1, -1, y, x))
            for f in faces:
                keys.append(tuple(sorted(int(row[i]) for i in f)))

    # Serial first-seen dedup: the id order is a pure function of the key order.
    node_id = {}
    new_nodes = []
    slot_id = np.empty(len(keys), dtype=np.int64)
    for i, key in enumerate(keys):
        found = node_id.get(key)
        if found is None:
            found = n + len(new_nodes)
            node_id[key] = found
            new_nodes.append(key)
        slot_id[i] = found

    body_base = n + len(new_nodes)
    n_bodies = sum(len(data) for cell_type, data in blocks if cell_type in _HAS_BODY)
    n_out = body_base + n_bodies

    def _extend(values):
        """Append the deduped-entity and body means to a per-point array.

        The input dtype is preserved, matching the C++ core (which writes the
        mean back through the destination array's own dtype).
        """
        out = np.empty((n_out,) + values.shape[1:], dtype=values.dtype)
        out[:n] = values
        for i, key in enumerate(new_nodes):
            corners = key[2:] if key[0] < 0 else key
            out[n + i] = values[list(corners)].mean(axis=0)
        body = body_base
        for cell_type, data in blocks:
            if cell_type not in _HAS_BODY:
                continue
            nc = _NUM_CORNERS[cell_type]
            for row in data:
                out[body] = values[[int(x) for x in row[:nc]]].mean(axis=0)
                body += 1
        return out

    out_points = _extend(np.asarray(mesh.points))

    # --- child connectivity
    new_cells = []
    cell_maps = []
    slot = 0
    body = body_base
    for cell_type, data in blocks:
        edges = _EDGES[cell_type]
        faces = _QUAD_FACES[cell_type]
        nslots = len(edges) + len(faces)
        ncorners = _NUM_CORNERS[cell_type]
        children = _CHILDREN[cell_type]
        ncells = len(data)

        local = np.empty(
            (ncells, ncorners + nslots + (1 if cell_type in _HAS_BODY else 0)),
            dtype=np.int64,
        )
        local[:, :ncorners] = data[:, :ncorners]
        if nslots:
            local[:, ncorners : ncorners + nslots] = slot_id[
                slot : slot + ncells * nslots
            ].reshape(ncells, nslots)
        if cell_type in _HAS_BODY:
            local[:, -1] = np.arange(body, body + ncells, dtype=np.int64)
            body += ncells
        slot += ncells * nslots

        # Parent-major / child-minor: one gather of every child's node columns,
        # reshaped so parent c's children land at rows [c*k, (c+1)*k).
        conn = local[:, [i for child in children for i in child]].reshape(
            ncells * len(children), len(children[0])
        )
        new_cells.append((cell_type, conn))
        cell_maps.append(np.arange(ncells, dtype=np.int64) * len(children))

    point_data = {}
    for key, value in mesh.point_data.items():
        value = np.asarray(value)
        if value.shape[:1] != (n,):
            point_data[key] = value.copy()
            continue
        point_data[key] = _extend(value)

    cell_data = {}
    for key, blk in mesh.cell_data.items():
        out_blk = []
        for b, value in enumerate(blk):
            if value is None or b >= len(blocks):
                out_blk.append(None if value is None else np.asarray(value).copy())
                continue
            value = np.asarray(value)
            nchildren = len(_CHILDREN[blocks[b][0]])
            if value.shape[:1] != (len(blocks[b][1]),):
                out_blk.append(value.copy())
                continue
            out_blk.append(np.repeat(value, nchildren, axis=0))
        cell_data[key] = out_blk

    out = Mesh(
        out_points,
        new_cells,
        point_data=point_data,
        cell_data=cell_data,
        field_data={
            k: (v.copy() if isinstance(v, np.ndarray) else v)
            for k, v in mesh.field_data.items()
        },
    )
    return out, np.arange(n, dtype=np.int64), cell_maps


def _refine_py(mesh, levels, record_parent_ids):
    if levels <= 0:
        out = Mesh(
            np.asarray(mesh.points).copy(),
            [(cb.type, np.asarray(cb.data).copy()) for cb in mesh.cells],
            point_data={k: np.asarray(v).copy() for k, v in mesh.point_data.items()},
            cell_data={
                k: [None if v is None else np.asarray(v).copy() for v in blk]
                for k, blk in mesh.cell_data.items()
            },
            field_data=dict(mesh.field_data),
        )
        point_map = np.arange(len(mesh.points), dtype=np.int64)
        cell_maps = [np.arange(len(cb.data), dtype=np.int64) for cb in mesh.cells]
    else:
        out, point_map, cell_maps = _refine_once_py(mesh)
        for _ in range(levels - 1):
            out, next_pm, next_cm = _refine_once_py(out)
            point_map = next_pm[point_map]
            cell_maps = [next_cm[b][cm] for b, cm in enumerate(cell_maps)]

    if record_parent_ids:
        out.cell_data["refine:parent_cell"] = _parent_ids(out, cell_maps)
    return out, point_map, cell_maps


def _parent_ids(out, cell_maps):
    """Per output cell, the index of its ORIGINAL ancestor within its block."""
    blocks = []
    for b, cb in enumerate(out.cells):
        ncells = len(cb.data)
        ids = np.full(ncells, -1, dtype=np.int64)
        first = np.asarray(cell_maps[b], dtype=np.int64) if b < len(cell_maps) else None
        if first is not None and len(first):
            ends = np.append(first[1:], ncells)
            for p, (lo, hi) in enumerate(zip(first, ends)):
                ids[lo:hi] = p
        blocks.append(ids)
    return blocks


# --------------------------------------------------------------------------- #
# set remapping (never reaches the C++ core)                                   #
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
                first = np.asarray(cell_maps[b], dtype=np.int64)
                total = len(out.cells[b].data) if b < len(out.cells) else 0
                # The end of parent c's child range is the next parent's start.
                ends = np.append(first[1:], total) if len(first) else first
                expanded = []
                for c in np.asarray(idx, dtype=np.int64):
                    if c < 0 or c >= len(first):
                        continue
                    expanded.append(np.arange(first[c], ends[c], dtype=np.int64))
                remapped.append(
                    np.concatenate(expanded)
                    if expanded
                    else np.empty(0, dtype=np.int64)
                )
            new_cs[name] = remapped
        out.cell_sets = new_cs


# --------------------------------------------------------------------------- #
# public API                                                                   #
# --------------------------------------------------------------------------- #
def refine(mesh, levels: int = 1, record_parent_ids: bool = False) -> Mesh:
    """Uniformly refine a mesh, subdividing every cell into same-type children.

    :param mesh: the mesh to refine (unchanged).
    :param levels: how many times to apply the subdivision templates. ``0`` (or
        less) returns an unchanged copy.
    :param record_parent_ids: attach an int64 ``refine:parent_cell`` ``cell_data``
        array naming, per output cell, the **original** input cell it descends
        from within its own block.
    :returns: the refined mesh.
    :raises ValueError: on a higher-order, ``pyramid``, or ragged cell block,
        none of which can be subdivided into same-type children.
    """
    levels = int(levels)

    out = None
    point_map = None
    cell_maps = None
    try:
        from . import _core

        res = _core.refine(mesh, levels, bool(record_parent_ids))
        out = res["mesh"]
        point_map = np.asarray(res["point_map"])
        cell_maps = [np.asarray(a) for a in res["cell_maps"]]
    except ValueError:
        # A genuine user error (a block with no same-type subdivision) must not
        # fall through to the numpy path and silently succeed.
        raise
    except Exception:
        out = None
    if out is None:
        out, point_map, cell_maps = _refine_py(mesh, levels, record_parent_ids)

    _remap_sets(mesh, out, point_map, cell_maps)
    return out
