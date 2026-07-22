"""Surface/boundary extraction: the general form of skin extraction.

A *facet* (a face of a 3D cell, or an edge of a 2D cell) is on the boundary
when its sorted corner-node key occurs exactly once across all cells of the
chosen dimension. :func:`extract_surface` picks the dimension automatically:
boundary faces for a mesh with volume cells, boundary edges for a surface mesh.

The public entry point uses the C++ core when available and falls back to the
pure-numpy implementation here. The face tables are shared with
``_skin.py``; the edge tables below are the Python twin of
``src/cpp/include/meshioplusplus/detail/cell_edges.hpp`` — KEEP THEM IN SYNC.
"""

from __future__ import annotations

import numpy as np

from ._common import warn
from ._mesh import Mesh
from ._skin import _CELL_FACES, _is_volume_type

# Surface cell type -> list of (edge_type, num_corners, local node indices).
# The Python twin of detail/cell_edges.hpp — keep the two in sync.
_CELL_EDGES = {
    "triangle": [
        ("line", 2, (0, 1)),
        ("line", 2, (1, 2)),
        ("line", 2, (2, 0)),
    ],
    "triangle6": [
        ("line3", 2, (0, 1, 3)),
        ("line3", 2, (1, 2, 4)),
        ("line3", 2, (2, 0, 5)),
    ],
    "quad": [
        ("line", 2, (0, 1)),
        ("line", 2, (1, 2)),
        ("line", 2, (2, 3)),
        ("line", 2, (3, 0)),
    ],
    "quad8": [
        ("line3", 2, (0, 1, 4)),
        ("line3", 2, (1, 2, 5)),
        ("line3", 2, (2, 3, 6)),
        ("line3", 2, (3, 0, 7)),
    ],
    "quad9": [
        ("line3", 2, (0, 1, 4)),
        ("line3", 2, (1, 2, 5)),
        ("line3", 2, (2, 3, 6)),
        ("line3", 2, (3, 0, 7)),
    ],
}

_FACE_OUT_TYPES = ["triangle", "triangle6", "quad", "quad8", "quad9"]
_FACE_OUT_NODES = {"triangle": 3, "triangle6": 6, "quad": 4, "quad8": 8, "quad9": 9}
_EDGE_OUT_TYPES = ["line", "line3"]
_EDGE_OUT_NODES = {"line": 2, "line3": 3}

_SURFACE_PREFIXES = ("triangle", "quad", "polygon")


def _is_surface_type(cell_type: str) -> bool:
    return cell_type.startswith(_SURFACE_PREFIXES) or cell_type in (
        "VTK_LAGRANGE_TRIANGLE",
        "VTK_LAGRANGE_QUADRILATERAL",
    )


def _has_surface_extractable_cells(mesh) -> bool:
    if any(block.type in _CELL_FACES for block in mesh.cells):
        return True
    if any(_is_volume_type(block.type) for block in mesh.cells):
        return False
    return any(block.type in _CELL_EDGES for block in mesh.cells)


def _extract_surface_py(mesh, record_parent_ids: bool = False) -> Mesh:
    """Pure-numpy reference implementation of :func:`extract_surface`."""
    face_mode = any(_is_volume_type(block.type) for block in mesh.cells)
    if face_mode:
        tables, out_types, out_nodes = _CELL_FACES, _FACE_OUT_TYPES, _FACE_OUT_NODES
        same_dim, noun, dim, key_w, conn_w = _is_volume_type, "volume", 3, 4, 9
    elif any(_is_surface_type(block.type) for block in mesh.cells):
        tables, out_types, out_nodes = _CELL_EDGES, _EDGE_OUT_TYPES, _EDGE_OUT_NODES
        same_dim, noun, dim, key_w, conn_w = _is_surface_type, "surface", 2, 4, 3
    else:
        raise ValueError(
            "extract_surface: mesh contains no supported 2D or 3D cell block"
        )

    out_index = {name: i for i, name in enumerate(out_types)}

    keys_parts = []
    conn_parts = []
    type_parts = []
    parent_parts = []
    any_supported = False
    global_base = 0
    for block in mesh.cells:
        nc = (
            np.asarray(block.data).shape[0] if block.type in tables else len(block.data)
        )
        base = global_base
        global_base += nc
        if not same_dim(block.type):
            continue
        if block.type not in tables:
            warn(
                f"extract_surface: {noun} cell block '{block.type}' is not "
                "supported; skipping it."
            )
            continue
        any_supported = True
        conn = np.asarray(block.data, dtype=np.int64)
        facets = tables[block.type]
        nf = len(facets)
        keys = np.full((nc, nf, key_w), -1, dtype=np.int64)
        fconn = np.full((nc, nf, conn_w), -1, dtype=np.int64)
        tidx = np.empty(nf, dtype=np.int64)
        for j, (facet_type, num_corners, ids) in enumerate(facets):
            keys[:, j, :num_corners] = np.sort(conn[:, ids[:num_corners]], axis=1)
            fconn[:, j, : len(ids)] = conn[:, ids]
            tidx[j] = out_index[facet_type]
        keys_parts.append(keys.reshape(nc * nf, key_w))
        conn_parts.append(fconn.reshape(nc * nf, conn_w))
        type_parts.append(np.tile(tidx, nc))
        parent_parts.append(np.repeat(base + np.arange(nc, dtype=np.int64), nf))

    if not any_supported:
        raise ValueError(
            f"extract_surface: mesh contains no supported {dim}D {noun} cell block"
        )

    keys = np.concatenate(keys_parts)
    fconn = np.concatenate(conn_parts)
    type_idx = np.concatenate(type_parts)
    parents = np.concatenate(parent_parts)

    _, inverse, counts = np.unique(
        keys, axis=0, return_inverse=True, return_counts=True
    )
    boundary = counts[inverse.reshape(-1)] == 1

    out_blocks = []
    out_parents = []
    for name in out_types:
        sel = boundary & (type_idx == out_index[name])
        if not np.any(sel):
            continue
        out_blocks.append((name, fconn[sel][:, : out_nodes[name]]))
        out_parents.append(parents[sel])

    used = np.unique(np.concatenate([c.reshape(-1) for _, c in out_blocks]))
    cells = [
        (name, np.searchsorted(used, c).astype(np.int64)) for name, c in out_blocks
    ]
    points = mesh.points[used]
    point_data = {
        key: value[used]
        for key, value in mesh.point_data.items()
        if np.asarray(value).shape[:1] == (len(mesh.points),)
    }
    cell_data = {}
    if record_parent_ids:
        cell_data["surface:parent_cell"] = [
            p.reshape(-1, 1).astype(np.int64) for p in out_parents
        ]
    return Mesh(points, cells, point_data=point_data, cell_data=cell_data)


def extract_surface(mesh, record_parent_ids: bool = False) -> Mesh:
    """Extract the boundary of a mesh's highest-dimension cells as a new mesh.

    Returns boundary faces (``triangle``/``quad`` and quadratic variants) for a
    mesh with 3D volume cells, or boundary edges (``line``/``line3``) for a 2D
    surface mesh. Facets whose sorted corner-node key occurs exactly once are
    kept; points are compacted to the referenced subset. With
    ``record_parent_ids=True`` a ``cell_data`` array ``"surface:parent_cell"``
    records, for each output facet, the global index of the owning input cell.
    Raises ``ValueError`` when the mesh has no supported 2D/3D cell block.
    """
    try:
        from . import _core

        return _core.extract_surface(mesh, record_parent_ids)
    except ValueError:
        raise
    except Exception:
        pass
    return _extract_surface_py(mesh, record_parent_ids)
