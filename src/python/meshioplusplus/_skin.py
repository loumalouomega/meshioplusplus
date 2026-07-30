"""Boundary-skin extraction: derive the surface mesh of a 3D volume mesh.

Implements the face-hashing algorithm of Kratos Multiphysics'
``SkinDetectionProcess``: a face whose sorted corner-node key occurs exactly
once across all volume cells is a boundary face; faces seen twice are
interior and cancel out.

The public entry point is :func:`extract_skin`, which uses the C++ core when
available and falls back to the pure-numpy implementation here (the format
shim pattern). The face tables below are the Python twin of
``src/cpp/include/meshioplusplus/detail/cell_faces.hpp`` — KEEP THEM IN SYNC:
every face row lists the corner nodes first (wound outward), then the
mid-edge nodes (mid of corner k -> corner k+1), then the face center, so each
row is itself a valid meshio/VTK ``triangle6``/``quad8``/``quad9`` ordering.
"""

from __future__ import annotations

import numpy as np

from ._common import warn
from ._mesh import Mesh

# Volume cell type -> list of (face_type, num_corners, local node indices).
# The Python twin of detail/cell_faces.hpp — keep the two in sync.
_CELL_FACES = {
    "tetra": [
        ("triangle", 3, (0, 1, 3)),
        ("triangle", 3, (1, 2, 3)),
        ("triangle", 3, (2, 0, 3)),
        ("triangle", 3, (0, 2, 1)),
    ],
    "tetra10": [
        ("triangle6", 3, (0, 1, 3, 4, 8, 7)),
        ("triangle6", 3, (1, 2, 3, 5, 9, 8)),
        ("triangle6", 3, (2, 0, 3, 6, 7, 9)),
        ("triangle6", 3, (0, 2, 1, 6, 5, 4)),
    ],
    "hexahedron": [
        ("quad", 4, (0, 4, 7, 3)),
        ("quad", 4, (1, 2, 6, 5)),
        ("quad", 4, (0, 1, 5, 4)),
        ("quad", 4, (3, 7, 6, 2)),
        ("quad", 4, (0, 3, 2, 1)),
        ("quad", 4, (4, 5, 6, 7)),
    ],
    "hexahedron20": [
        ("quad8", 4, (0, 4, 7, 3, 16, 15, 19, 11)),
        ("quad8", 4, (1, 2, 6, 5, 9, 18, 13, 17)),
        ("quad8", 4, (0, 1, 5, 4, 8, 17, 12, 16)),
        ("quad8", 4, (3, 7, 6, 2, 19, 14, 18, 10)),
        ("quad8", 4, (0, 3, 2, 1, 11, 10, 9, 8)),
        ("quad8", 4, (4, 5, 6, 7, 12, 13, 14, 15)),
    ],
    # Mid-face numbering follows vtkTriQuadraticHexahedron::Faces; corrected in
    # v9.9.0 (20/22/23 were a permuted 3-cycle). Kept in sync with the C++ twin
    # detail/cell_faces.cpp -- see its comment for the full derivation.
    "hexahedron27": [
        ("quad9", 4, (0, 4, 7, 3, 16, 15, 19, 11, 20)),
        ("quad9", 4, (1, 2, 6, 5, 9, 18, 13, 17, 21)),
        ("quad9", 4, (0, 1, 5, 4, 8, 17, 12, 16, 22)),
        ("quad9", 4, (3, 7, 6, 2, 19, 14, 18, 10, 23)),
        ("quad9", 4, (0, 3, 2, 1, 11, 10, 9, 8, 24)),
        ("quad9", 4, (4, 5, 6, 7, 12, 13, 14, 15, 25)),
    ],
    "wedge": [
        ("triangle", 3, (0, 2, 1)),
        ("triangle", 3, (3, 4, 5)),
        ("quad", 4, (0, 1, 4, 3)),
        ("quad", 4, (1, 2, 5, 4)),
        ("quad", 4, (2, 0, 3, 5)),
    ],
    "wedge15": [
        ("triangle6", 3, (0, 2, 1, 8, 7, 6)),
        ("triangle6", 3, (3, 4, 5, 9, 10, 11)),
        ("quad8", 4, (0, 1, 4, 3, 6, 13, 9, 12)),
        ("quad8", 4, (1, 2, 5, 4, 7, 14, 10, 13)),
        ("quad8", 4, (2, 0, 3, 5, 8, 12, 11, 14)),
    ],
    "wedge18": [
        ("triangle6", 3, (0, 2, 1, 8, 7, 6)),
        ("triangle6", 3, (3, 4, 5, 9, 10, 11)),
        ("quad9", 4, (0, 1, 4, 3, 6, 13, 9, 12, 15)),
        ("quad9", 4, (1, 2, 5, 4, 7, 14, 10, 13, 16)),
        ("quad9", 4, (2, 0, 3, 5, 8, 12, 11, 14, 17)),
    ],
    "pyramid": [
        ("quad", 4, (0, 3, 2, 1)),
        ("triangle", 3, (0, 1, 4)),
        ("triangle", 3, (1, 2, 4)),
        ("triangle", 3, (2, 3, 4)),
        ("triangle", 3, (3, 0, 4)),
    ],
    "pyramid13": [
        ("quad8", 4, (0, 3, 2, 1, 8, 7, 6, 5)),
        ("triangle6", 3, (0, 1, 4, 5, 10, 9)),
        ("triangle6", 3, (1, 2, 4, 6, 11, 10)),
        ("triangle6", 3, (2, 3, 4, 7, 12, 11)),
        ("triangle6", 3, (3, 0, 4, 8, 9, 12)),
    ],
    "pyramid14": [
        ("quad9", 4, (0, 3, 2, 1, 8, 7, 6, 5, 13)),
        ("triangle6", 3, (0, 1, 4, 5, 10, 9)),
        ("triangle6", 3, (1, 2, 4, 6, 11, 10)),
        ("triangle6", 3, (2, 3, 4, 7, 12, 11)),
        ("triangle6", 3, (3, 0, 4, 8, 9, 12)),
    ],
}

# Canonical output block order (matches the C++ extractor).
_OUT_TYPES = ["triangle", "triangle6", "quad", "quad8", "quad9"]
_OUT_TYPE_INDEX = {name: i for i, name in enumerate(_OUT_TYPES)}
_OUT_TYPE_NODES = {"triangle": 3, "triangle6": 6, "quad": 4, "quad8": 8, "quad9": 9}
_LINEARIZED = {
    "triangle": "triangle",
    "triangle6": "triangle",
    "quad": "quad",
    "quad8": "quad",
    "quad9": "quad",
}

# All volume cell type names (topological dimension 3), for warning about
# unsupported ones.
_VOLUME_PREFIXES = ("tetra", "hexahedron", "wedge", "pyramid", "polyhedron")


def _is_volume_type(cell_type: str) -> bool:
    return cell_type.startswith(_VOLUME_PREFIXES) or cell_type in (
        "VTK_LAGRANGE_TETRAHEDRON",
        "VTK_LAGRANGE_HEXAHEDRON",
        "VTK_LAGRANGE_WEDGE",
        "VTK_LAGRANGE_PYRAMID",
    )


def _has_skinnable_cells(mesh) -> bool:
    """Whether the mesh has at least one cell block the extractor supports."""
    return any(block.type in _CELL_FACES for block in mesh.cells)


def _extract_skin_py(
    mesh, linearize: bool = False, record_parent_ids: bool = False
) -> Mesh:
    """Pure-numpy reference implementation of :func:`extract_skin`.

    With ``record_parent_ids`` the result carries an ``Int64``
    ``"surface:parent_cell"`` cell-data array giving, per output facet, the
    global input-cell index of the cell that owns it -- the same contract as
    :func:`meshioplusplus.extract_surface`, and the twin of the C++
    ``detail::surface_extract(..., recordParentIds=True)`` the SVG/TikZ writers
    use to colour a projected volume mesh by cell data.
    """
    # `base` counts over EVERY block, including the ones skipped below, because
    # "surface:parent_cell" indexes the input mesh's cells globally. This
    # mirrors the C++ extractor advancing global_cell_base before its `continue`.
    supported = []
    base = 0
    for block in mesh.cells:
        block_base = base
        base += (
            np.asarray(block.data).shape[0]
            if block.type in _CELL_FACES
            else len(block.data)
        )
        if not _is_volume_type(block.type):
            continue
        if block.type in _CELL_FACES:
            supported.append((block, block_base))
        else:
            warn(
                f"extract_skin: volume cell block '{block.type}' is not supported; "
                "skipping it."
            )
    if not supported:
        raise ValueError(
            "extract_skin: mesh contains no supported 3D volume cell block"
        )

    # Enumerate every face of every cell in block-major / cell-major /
    # face-slot order (matching the C++ extractor), building for each face
    # its sorted-corner key, its node ids (padded to 9), and its output type.
    keys_parts = []
    conn_parts = []
    type_parts = []
    parent_parts = []
    for block, block_base in supported:
        conn = np.asarray(block.data, dtype=np.int64)
        nc = conn.shape[0]
        faces = _CELL_FACES[block.type]
        nf = len(faces)
        keys = np.full((nc, nf, 4), -1, dtype=np.int64)
        fconn = np.full((nc, nf, 9), -1, dtype=np.int64)
        tidx = np.empty(nf, dtype=np.int64)
        for j, (face_type, num_corners, ids) in enumerate(faces):
            keys[:, j, :num_corners] = np.sort(conn[:, ids[:num_corners]], axis=1)
            fconn[:, j, : len(ids)] = conn[:, ids]
            out_type = _LINEARIZED[face_type] if linearize else face_type
            tidx[j] = _OUT_TYPE_INDEX[out_type]
        keys_parts.append(keys.reshape(nc * nf, 4))
        conn_parts.append(fconn.reshape(nc * nf, 9))
        type_parts.append(np.tile(tidx, nc))
        parent_parts.append(np.repeat(block_base + np.arange(nc, dtype=np.int64), nf))
    keys = np.concatenate(keys_parts)
    fconn = np.concatenate(conn_parts)
    type_idx = np.concatenate(type_parts)
    parents = np.concatenate(parent_parts)

    # Boundary faces = sorted-corner keys occurring exactly once.
    _, inverse, counts = np.unique(
        keys, axis=0, return_inverse=True, return_counts=True
    )
    boundary = counts[inverse.reshape(-1)] == 1

    # Partition into output blocks in the fixed canonical order, keeping
    # enumeration order within each type.
    out_blocks = []
    out_parents = []
    for name in _OUT_TYPES:
        sel = boundary & (type_idx == _OUT_TYPE_INDEX[name])
        if not np.any(sel):
            continue
        n_nodes = (
            (3 if name == "triangle" else 4) if linearize else _OUT_TYPE_NODES[name]
        )
        out_blocks.append((name, fconn[sel][:, :n_nodes]))
        out_parents.append(parents[sel])

    # Compaction: keep only referenced points, in ascending original-id order
    # (np.unique is sorted — the C++ extractor relies on the same order).
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


def extract_skin(mesh, linearize: bool = False) -> Mesh:
    """Extract the boundary skin of a volume mesh as a new surface mesh.

    Keeps the faces of the supported 3D volume cells (tetra, hexahedron,
    wedge, pyramid, tetra10, hexahedron20/27, wedge15, pyramid13/14) whose
    sorted corner-node key occurs exactly once. Points are compacted to the
    referenced subset (``point_data`` rows follow); ``cell_data``,
    ``field_data``, and sets are dropped. With ``linearize=True`` only the
    corner nodes are emitted (``triangle``/``quad`` output even for
    higher-order input). Raises ``ValueError`` when the mesh has no
    supported volume cell block.
    """
    try:
        from . import _core

        return _core.extract_skin(mesh, linearize)
    except ValueError:
        raise
    except Exception:
        pass
    return _extract_skin_py(mesh, linearize)
