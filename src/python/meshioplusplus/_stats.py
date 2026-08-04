"""Geometric statistics of a mesh (read-only), complementing the topological
``info``.

A dependency-free mesh *operation* (not a file format): it reports the bounding
box + extents, centroid, per-cell-type counts, total surface area (of 2D cells
plus the boundary of 3D cells), total signed & unsigned volume of 3D cells, and
the count of inverted cells — without modifying the mesh. The C++ core
(``_core.compute_stats``) does the work; this module is the thin shim (try C++,
fall back to a pure-numpy reference).

Public API:

* :func:`compute_stats` — return a dict of geometric statistics.
"""

from __future__ import annotations

import numpy as np

# Linear-corner face node lists (matching detail/cell_faces.hpp), outward-wound.
_FACES = {
    "tetra": [[0, 1, 3], [1, 2, 3], [2, 0, 3], [0, 2, 1]],
    "hexahedron": [
        [0, 4, 7, 3],
        [1, 2, 6, 5],
        [0, 1, 5, 4],
        [3, 7, 6, 2],
        [0, 3, 2, 1],
        [4, 5, 6, 7],
    ],
    "wedge": [[0, 2, 1], [3, 4, 5], [0, 1, 4, 3], [1, 2, 5, 4], [2, 0, 3, 5]],
    "pyramid": [[0, 3, 2, 1], [0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
}

_CORNERS = {
    "vertex": 1,
    "line": 2,
    "triangle": 3,
    "quad": 4,
    "tetra": 4,
    "hexahedron": 8,
    "wedge": 6,
    "pyramid": 5,
}

_DIM = {
    "vertex": 0,
    "line": 1,
    "triangle": 2,
    "quad": 2,
    "tetra": 3,
    "hexahedron": 3,
    "wedge": 3,
    "pyramid": 3,
}


def _base_type(cell_type):
    for base in _CORNERS:
        if cell_type == base or cell_type.startswith(base):
            return base
    return None


def _coords_3d(points):
    p = np.asarray(points, dtype=np.float64)
    if p.ndim == 1:
        p = p.reshape(-1, 1)
    out = np.zeros((p.shape[0], 3), dtype=np.float64)
    out[:, : min(p.shape[1], 3)] = p[:, : min(p.shape[1], 3)]
    return out


def _poly_area(corners):
    # Newell normal magnitude / 2 for a corner polygon (n, k, 3).
    s = np.zeros((corners.shape[0], 3), dtype=np.float64)
    k = corners.shape[1]
    for i in range(k):
        s += np.cross(corners[:, i], corners[:, (i + 1) % k])
    return 0.5 * np.linalg.norm(s, axis=1)


def _signed_volume(corner_coords, faces):
    # Divergence-theorem volume over the outward faces, each fanned about its
    # OWN corner average -- the twin of detail/polyhedron.hpp's poly_measure.
    #
    # Before v9.16.0 both sides fanned about the face's first node instead. That
    # apex is a function of the cell's storage rather than of the face, so two
    # cells sharing a warped quad could split it along different diagonals; the
    # corner average is a function of the face alone and cannot. See
    # doc/polyhedra.md. The two agree exactly for planar faces and differ only
    # on warped ones, so this must move in lockstep with the C++ side or a
    # Python-fallback build (Windows CI) reports different volumes from a
    # native one.
    nc, k, _ = corner_coords.shape
    # Recentre on the cell's corner average: V = sum(x . A)/3 only telescopes
    # because sum(A) == 0, so far from the origin the raw form cancels away its
    # own significant digits.
    origin = corner_coords.sum(axis=1) / float(k)
    local = corner_coords - origin[:, None, :]
    vol6 = np.zeros(nc, dtype=np.float64)
    for f in faces:
        m = len(f)
        fc = np.zeros((nc, 3), dtype=np.float64)
        for idx in f:
            fc = fc + local[:, idx]
        fc = fc / float(m)
        for i in range(m):
            a = local[:, f[i]]
            b = local[:, f[(i + 1) % m]]
            vol6 += np.einsum("ij,ij->i", fc, np.cross(a, b))
    return vol6 / 6.0


def _compute_stats_py(mesh):
    pts = _coords_3d(mesh.points)
    n = pts.shape[0]
    dim = np.asarray(mesh.points).shape[1] if np.asarray(mesh.points).ndim == 2 else 0

    rep = {
        "num_points": int(n),
        "num_cells": 0,
        "bbox_min": (0.0, 0.0, 0.0),
        "bbox_max": (0.0, 0.0, 0.0),
        "extent": (0.0, 0.0, 0.0),
        "centroid": (0.0, 0.0, 0.0),
        "cell_type_counts": {},
        "total_area": 0.0,
        "signed_volume": 0.0,
        "unsigned_volume": 0.0,
        "num_inverted": 0,
    }
    if n > 0:
        mn = pts.min(axis=0)
        mx = pts.max(axis=0)
        ct = pts.mean(axis=0)
        keep = [i < dim for i in range(3)]
        mn = [float(mn[i]) if keep[i] else 0.0 for i in range(3)]
        mx = [float(mx[i]) if keep[i] else 0.0 for i in range(3)]
        ct = [float(ct[i]) if keep[i] else 0.0 for i in range(3)]
        rep["bbox_min"] = tuple(mn)
        rep["bbox_max"] = tuple(mx)
        rep["extent"] = tuple(mx[i] - mn[i] for i in range(3))
        rep["centroid"] = tuple(ct)

    counts = {}
    any_3d = False
    total_area = 0.0
    signed_v = 0.0
    unsigned_v = 0.0
    inverted = 0
    for cb in mesh.cells:
        rep["num_cells"] += len(cb.data)
        counts[cb.type] = counts.get(cb.type, 0) + len(cb.data)
        base = _base_type(cb.type)
        if base is None:
            continue
        data = np.asarray(cb.data)
        if data.ndim != 2 or len(data) == 0:
            continue
        cc = _CORNERS[base]
        cdim = _DIM[base]
        corner_coords = pts[data[:, :cc]]  # (nc, cc, 3)
        if cdim == 2 and base in ("triangle", "quad"):
            total_area += float(_poly_area(corner_coords).sum())
        elif cdim == 3 and base in _FACES:
            any_3d = True
            v = _signed_volume(corner_coords, _FACES[base])
            signed_v += float(v.sum())
            unsigned_v += float(np.abs(v).sum())
            inverted += int(np.count_nonzero(v < 0.0))

    if any_3d:
        try:
            from ._surface import extract_surface

            surf = extract_surface(mesh)
            sp = _coords_3d(surf.points)
            for cb in surf.cells:
                base = _base_type(cb.type)
                if base in ("triangle", "quad"):
                    data = np.asarray(cb.data)
                    if data.ndim == 2 and len(data):
                        cc = _CORNERS[base]
                        total_area += float(_poly_area(sp[data[:, :cc]]).sum())
        except Exception:
            pass

    rep["cell_type_counts"] = counts
    rep["total_area"] = total_area
    rep["signed_volume"] = signed_v
    rep["unsigned_volume"] = unsigned_v
    rep["num_inverted"] = inverted
    return rep


def compute_stats(mesh) -> dict:
    """Compute geometric statistics of ``mesh`` (read-only).

    Returns
    -------
    dict
        with keys ``num_points``, ``num_cells``, ``bbox_min``/``bbox_max``/
        ``extent``/``centroid`` (3-tuples), ``cell_type_counts`` (``{type:
        count}``), ``total_area``, ``signed_volume``, ``unsigned_volume``, and
        ``num_inverted``.
    """
    try:
        from . import _core

        res = _core.compute_stats(mesh)
        return {
            "num_points": int(res["num_points"]),
            "num_cells": int(res["num_cells"]),
            "bbox_min": tuple(res["bbox_min"]),
            "bbox_max": tuple(res["bbox_max"]),
            "extent": tuple(res["extent"]),
            "centroid": tuple(res["centroid"]),
            "cell_type_counts": {k: int(v) for k, v in res["cell_type_counts"].items()},
            "total_area": float(res["total_area"]),
            "signed_volume": float(res["signed_volume"]),
            "unsigned_volume": float(res["unsigned_volume"]),
            "num_inverted": int(res["num_inverted"]),
        }
    except Exception:
        pass
    return _compute_stats_py(mesh)
