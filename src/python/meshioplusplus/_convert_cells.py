"""Conversion of a mesh's element representation (linearize / simplexify / elevate).

A dependency-free mesh *operation* (not a file format): it rewrites which cell
*types* a mesh is built from, leaving the mesh it describes intact. The C++ core
(``_core.convert_cells``) does the work; this module is the thin shim (try C++,
fall back to a pure-numpy reference) plus the ``point_sets``/``cell_sets``
remapping the core cannot see across the Python boundary.

Three modes:

* ``"linearize"`` — every higher-order cell becomes its linear base
  (``tetra10`` -> ``tetra``, ...), dropping the higher-order nodes and pruning
  the points that become unreferenced.
* ``"simplexify"`` — every cell is decomposed into simplices of the same
  topological dimension (quad -> 2 triangles, hexahedron -> 6 tetra, ...),
  reusing the parent's own corner nodes. Higher-order input is linearized first.
* ``"elevate"`` — every linear cell is promoted to its serendipity quadratic
  counterpart (``triangle`` -> ``triangle6``, ...), adding one new node per
  unique edge at the edge midpoint.

Public API:

* :func:`convert_cells` — convert a mesh's element representation.
"""

from __future__ import annotations

import numpy as np

from ._mesh import Mesh

# --------------------------------------------------------------------------- #
# tables (Python twins of the ones in src/cpp/src/operations/convert_cells.cpp)    #
# --------------------------------------------------------------------------- #

#: Higher-order cell type -> its linear base. Corner nodes always come first in
#: the meshio/VTK orderings, so linearizing is "keep the leading N columns".
_LINEAR_BASE = {
    "line3": "line",
    "line4": "line",
    "triangle6": "triangle",
    "triangle10": "triangle",
    "quad8": "quad",
    "quad9": "quad",
    "quad16": "quad",
    "tetra10": "tetra",
    "tetra20": "tetra",
    "hexahedron20": "hexahedron",
    "hexahedron24": "hexahedron",
    "hexahedron27": "hexahedron",
    "hexahedron64": "hexahedron",
    "wedge15": "wedge",
    "wedge18": "wedge",
    "pyramid13": "pyramid",
    "pyramid14": "pyramid",
}

#: Corner counts of the linear bases (the number of leading columns to keep).
_NUM_CORNERS = {
    "vertex": 1,
    "line": 2,
    "triangle": 3,
    "quad": 4,
    "tetra": 4,
    "pyramid": 5,
    "wedge": 6,
    "hexahedron": 8,
}

#: Decomposition templates, all positively oriented for a well-oriented parent.
#: The hexahedron entry is the canonical Freudenthal fan around the main
#: diagonal 0-6 (every tetra contains the edge (0, 6)).
_SIMPLEX_TEMPLATES = {
    "quad": ("triangle", [[0, 1, 2], [0, 2, 3]]),
    "pyramid": ("tetra", [[0, 1, 2, 4], [0, 2, 3, 4]]),
    "wedge": ("tetra", [[0, 1, 2, 3], [1, 2, 3, 4], [2, 3, 4, 5]]),
    "hexahedron": (
        "tetra",
        [
            [0, 1, 2, 6],
            [0, 2, 3, 6],
            [0, 3, 7, 6],
            [0, 7, 4, 6],
            [0, 4, 5, 6],
            [0, 5, 1, 6],
        ],
    ),
}

#: Cell types that are already simplices of their own dimension.
_SIMPLEX_TYPES = {"vertex", "line", "triangle", "tetra"}

#: Linear type -> (quadratic target, corner count, edges in mid-node order).
#: Target node ``ncorners + k`` is the midpoint of edge ``k``.
_ELEVATE = {
    "line": ("line3", 2, [(0, 1)]),
    "triangle": ("triangle6", 3, [(0, 1), (1, 2), (2, 0)]),
    "quad": ("quad8", 4, [(0, 1), (1, 2), (2, 3), (3, 0)]),
    "tetra": ("tetra10", 4, [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]),
    "hexahedron": (
        "hexahedron20",
        8,
        [
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
    ),
    "wedge": (
        "wedge15",
        6,
        [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)],
    ),
    "pyramid": (
        "pyramid13",
        5,
        [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 4), (2, 4), (3, 4)],
    ),
}

#: Full-Lagrange targets that need face/body-centre nodes — out of scope here.
_FULL_LAGRANGE = {"quad9", "hexahedron27", "wedge18", "pyramid14"}

MODES = ("linearize", "simplexify", "elevate")


# --------------------------------------------------------------------------- #
# pure-numpy reference (rectangular blocks only; ragged is C++-core only)      #
# --------------------------------------------------------------------------- #
def _block_data(cb):
    data = np.asarray(cb.data)
    if data.ndim != 2:
        raise NotImplementedError(
            "convert_cells numpy fallback: ragged/polyhedron blocks are C++-core only"
        )
    return data


def _copy_field_data(mesh):
    return {
        k: (v.copy() if isinstance(v, np.ndarray) else v)
        for k, v in mesh.field_data.items()
    }


def _linearize_py(mesh, record_parent_ids):
    n = len(mesh.points)
    new_cells = []
    used = np.zeros(n, dtype=bool)
    for cb in mesh.cells:
        data = _block_data(cb)
        base = _LINEAR_BASE.get(cb.type, cb.type)
        if base == cb.type:
            conn = data.astype(np.int64)
        else:
            conn = data[:, : _NUM_CORNERS[base]].astype(np.int64)
        new_cells.append((base, conn))
        if conn.size:
            used[conn.reshape(-1)] = True

    new_point = np.full(n, -1, dtype=np.int64)
    point_src = np.nonzero(used)[0]
    new_point[point_src] = np.arange(len(point_src), dtype=np.int64)

    pts = np.asarray(mesh.points)
    out_points = pts[point_src] if len(point_src) else pts[:0]
    new_cells = [(t, new_point[c] if c.size else c) for t, c in new_cells]

    point_data = {}
    for key, value in mesh.point_data.items():
        value = np.asarray(value)
        point_data[key] = value[point_src] if value.shape[:1] == (n,) else value.copy()

    cell_data = {
        key: [None if v is None else np.asarray(v).copy() for v in blk]
        for key, blk in mesh.cell_data.items()
    }
    if record_parent_ids:
        cell_data["convert:parent_cell"] = [
            np.arange(len(cb.data), dtype=np.int64) for cb in mesh.cells
        ]

    out = Mesh(
        out_points,
        new_cells,
        point_data=point_data,
        cell_data=cell_data,
        field_data=_copy_field_data(mesh),
    )
    cell_maps = [np.arange(len(cb.data), dtype=np.int64) for cb in mesh.cells]
    return out, new_point, cell_maps


def _simplexify_py(mesh, record_parent_ids):
    # Higher-order input is linearized first (which also prunes the mid nodes
    # the decomposition would otherwise orphan).
    if any(_LINEAR_BASE.get(cb.type, cb.type) != cb.type for cb in mesh.cells):
        prepared, point_map, _ = _linearize_py(mesh, False)
    else:
        prepared, point_map = mesh, np.arange(len(mesh.points), dtype=np.int64)

    new_cells = []
    first_child = []
    parents_per_block = []
    for cb in prepared.cells:
        data = _block_data(cb)
        if cb.type == "polygon":
            # (n-2)-triangle fan around node 0. The numpy path only ever sees
            # the rectangular (uniform n-gon) shape; jagged rows raise above.
            n = data.shape[1]
            rows = [[0, k, k + 1] for k in range(1, n - 1)]
            tpl = ("triangle", rows)
        else:
            tpl = _SIMPLEX_TEMPLATES.get(cb.type)
        if tpl is None:
            new_cells.append((cb.type, data.astype(np.int64)))
            first_child.append(np.arange(len(data), dtype=np.int64))
            parents_per_block.append(np.arange(len(data), dtype=np.int64))
            continue
        child_type, rows = tpl
        per_cell = len(rows)
        # (ncells, per_cell, nodes) -> (ncells * per_cell, nodes), parent-major.
        conn = np.stack([data[:, r] for r in rows], axis=1).reshape(
            len(data) * per_cell, len(rows[0])
        )
        new_cells.append((child_type, conn.astype(np.int64)))
        first_child.append(np.arange(len(data), dtype=np.int64) * per_cell)
        parents_per_block.append(
            np.repeat(np.arange(len(data), dtype=np.int64), per_cell)
        )

    point_data = {k: np.asarray(v).copy() for k, v in prepared.point_data.items()}

    cell_data = {}
    for key, blk in prepared.cell_data.items():
        out_blk = []
        for b, value in enumerate(blk):
            if value is None or b >= len(parents_per_block):
                out_blk.append(None if value is None else np.asarray(value).copy())
                continue
            value = np.asarray(value)
            if value.shape[:1] == (len(prepared.cells[b].data),):
                out_blk.append(value[parents_per_block[b]])
            else:
                out_blk.append(value.copy())
        cell_data[key] = out_blk
    if record_parent_ids:
        cell_data["convert:parent_cell"] = [p.copy() for p in parents_per_block]

    out = Mesh(
        np.asarray(prepared.points).copy(),
        new_cells,
        point_data=point_data,
        cell_data=cell_data,
        field_data=_copy_field_data(prepared),
    )
    return out, point_map, first_child


def _elevate_py(mesh, record_parent_ids):
    n = len(mesh.points)
    for cb in mesh.cells:
        if cb.type in _FULL_LAGRANGE:
            raise ValueError(
                f"convert_cells: cannot elevate to the full-Lagrange type '{cb.type}' "
                "(face/body-centre nodes are not supported in this version)"
            )

    # Assign one new node id per unique edge, in first-appearance order — the
    # same order the C++ core's serial dedup pass produces.
    edge_id = {}
    per_block_mid = []
    for cb in mesh.cells:
        data = _block_data(cb)
        spec = _ELEVATE.get(cb.type)
        if spec is None:
            per_block_mid.append(None)
            continue
        _, _, edges = spec
        mids = np.empty((len(data), len(edges)), dtype=np.int64)
        for c in range(len(data)):
            for e, (a, b) in enumerate(edges):
                key = (int(data[c, a]), int(data[c, b]))
                if key[0] > key[1]:
                    key = (key[1], key[0])
                found = edge_id.get(key)
                if found is None:
                    found = n + len(edge_id)
                    edge_id[key] = found
                mids[c, e] = found
        per_block_mid.append(mids)

    endpoints = np.asarray(list(edge_id.keys()), dtype=np.int64).reshape(-1, 2)

    pts = np.asarray(mesh.points)
    if len(endpoints):
        mid_pts = 0.5 * (
            pts[endpoints[:, 0]].astype(np.float64)
            + pts[endpoints[:, 1]].astype(np.float64)
        )
        out_points = np.concatenate([pts, mid_pts.astype(pts.dtype)], axis=0)
    else:
        out_points = pts.copy()

    new_cells = []
    for cb, mids in zip(mesh.cells, per_block_mid):
        data = _block_data(cb)
        if mids is None:
            new_cells.append((cb.type, data.astype(np.int64)))
            continue
        target, ncorners, _ = _ELEVATE[cb.type]
        conn = np.concatenate([data[:, :ncorners].astype(np.int64), mids], axis=1)
        new_cells.append((target, conn))

    point_data = {}
    for key, value in mesh.point_data.items():
        value = np.asarray(value)
        if value.shape[:1] != (n,) or not len(endpoints):
            point_data[key] = value.copy()
            continue
        mid = 0.5 * (
            value[endpoints[:, 0]].astype(np.float64)
            + value[endpoints[:, 1]].astype(np.float64)
        )
        point_data[key] = np.concatenate([value, mid.astype(value.dtype)], axis=0)

    cell_data = {
        key: [None if v is None else np.asarray(v).copy() for v in blk]
        for key, blk in mesh.cell_data.items()
    }
    if record_parent_ids:
        cell_data["convert:parent_cell"] = [
            np.arange(len(cb.data), dtype=np.int64) for cb in mesh.cells
        ]

    out = Mesh(
        out_points,
        new_cells,
        point_data=point_data,
        cell_data=cell_data,
        field_data=_copy_field_data(mesh),
    )
    cell_maps = [np.arange(len(cb.data), dtype=np.int64) for cb in mesh.cells]
    return out, np.arange(n, dtype=np.int64), cell_maps


def _convert_cells_py(mesh, mode, record_parent_ids):
    if mode == "linearize":
        return _linearize_py(mesh, record_parent_ids)
    if mode == "simplexify":
        return _simplexify_py(mesh, record_parent_ids)
    return _elevate_py(mesh, record_parent_ids)


# --------------------------------------------------------------------------- #
# set remapping                                                                #
# --------------------------------------------------------------------------- #
def _remap_sets(src, out, point_map, cell_maps):
    """Remap the shim-only ``point_sets``/``cell_sets`` onto the converted mesh.

    ``cell_maps[b][c]`` is the index of parent ``c``'s *first* child, and a
    parent's children are contiguous, so a cell set entry expands to the whole
    ``[first[c], first[c + 1])`` range.
    """
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
def convert_cells(
    mesh,
    mode: str = "linearize",
    record_parent_ids: bool = False,
) -> Mesh:
    """Convert the element representation of a mesh.

    Parameters
    ----------
    mesh :
        the mesh to convert (unmodified).
    mode :
        ``"linearize"`` (default), ``"simplexify"``, or ``"elevate"`` — see the
        module docstring.
    record_parent_ids :
        when true, attach an Int64 ``convert:parent_cell`` cell_data array
        recording, per output cell, the index of the input cell it came from
        within its own block. Mainly useful under ``"simplexify"``.

    Returns
    -------
    Mesh
        the converted mesh. ``point_sets``/``cell_sets`` are remapped (under
        ``"simplexify"`` a cell set entry expands to the parent's children).

    Raises
    ------
    ValueError
        for an unknown ``mode``, a polyhedron block under ``"simplexify"``, or a
        full-Lagrange target (``quad9``/``hexahedron27``) under ``"elevate"``.
    """
    if mode not in MODES:
        raise ValueError(
            f"convert_cells: unknown mode '{mode}' (expected one of {', '.join(MODES)})"
        )

    out = None
    point_map = None
    cell_maps = None
    try:
        from . import _core

        res = _core.convert_cells(mesh, mode, bool(record_parent_ids))
        out = res["mesh"]
        point_map = np.asarray(res["point_map"])
        cell_maps = [np.asarray(a) for a in res["cell_maps"]]
    except ValueError:
        # A genuine user error (unsupported construct for this mode) must not
        # fall through to the numpy path and silently succeed.
        raise
    except Exception:
        out = None
    used_cpp = out is not None
    if out is None:
        out, point_map, cell_maps = _convert_cells_py(mesh, mode, record_parent_ids)

    # The C++ core carries named regions across itself (and therefore
    # `point_sets`/`cell_sets`, which are views over them), so remapping again
    # here would apply the maps twice. Only the numpy fallback needs it.
    if not used_cpp:
        _remap_sets(mesh, out, point_map, cell_maps)
    return out
