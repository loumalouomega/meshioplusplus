"""The feature edges of a surface, as a mesh of ``line`` cells.

The numpy twin of ``operations/feature_edges.hpp`` and of the crease test in
``detail/crease_edges.hpp``, which ``decimate`` and ``smooth`` pin nodes with.
Per edge of the surface's face rings:

- used once: **boundary**;
- used three or more times: **non-manifold**;
- used by two faces: their unit normals are compared -- after negating the
  second when both faces walk the edge the same way (an **inconsistent**
  pair) -- and the edge is **sharp** when ``dot < cos(feature_angle)``.

With any volume cell present, the surface is the skin of the volume cells.
"""

from __future__ import annotations

import numpy as np

from ._fallback import core_op_declined
from ._regions import block_bases

#: ``feature:kind`` values.
FEATURE, BOUNDARY, NON_MANIFOLD, INCONSISTENT = 1, 2, 3, 4

_PI = 3.14159265358979323846


def ring_unit_normals(xyz, start, nodes):
    """Per ring, Newell's unit normal (sum of ``cross(p_k, p_k+1)``), or zero."""
    start = np.asarray(start, dtype=np.int64)
    nodes = np.asarray(nodes, dtype=np.int64)
    nf = len(start) - 1
    out = np.zeros((nf, 3))
    if nf == 0:
        return out
    sizes = np.diff(start)
    for s in np.unique(sizes):
        s = int(s)
        faces = np.nonzero(sizes == s)[0]
        ring = nodes[start[faces][:, None] + np.arange(s)[None, :]]
        acc = np.zeros((len(faces), 3))
        for k in range(s):
            a = xyz[ring[:, k]]
            b = xyz[ring[:, (k + 1) % s]]
            acc = acc + np.column_stack(
                [
                    a[:, 1] * b[:, 2] - a[:, 2] * b[:, 1],
                    a[:, 2] * b[:, 0] - a[:, 0] * b[:, 2],
                    a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0],
                ]
            )
        norm = np.sqrt(
            acc[:, 0] * acc[:, 0] + acc[:, 1] * acc[:, 1] + acc[:, 2] * acc[:, 2]
        )
        ok = norm >= 1e-300
        inv = np.zeros_like(norm)
        inv[ok] = 1.0 / norm[ok]
        out[faces] = acc * inv[:, None]
    return out


def crease_edges(start, nodes, unit_normals, feature_angle):
    """Classify the edges of a polygonal surface.

    :returns: a dict of equal-length arrays -- ``lo``, ``hi``, ``uses``,
        ``inconsistent``, ``sharp``, ``angle`` (degrees, NaN unless a pair of
        two non-degenerate faces) -- for every boundary, non-manifold,
        inconsistent or sharp edge, ascending in ``(lo, hi)``.
    """
    start = np.asarray(start, dtype=np.int64)
    nodes = np.asarray(nodes, dtype=np.int64)
    normals = np.asarray(unit_normals, dtype=np.float64).reshape(-1, 3)
    nf = len(start) - 1
    sizes = np.diff(start) if nf > 0 else np.empty(0, dtype=np.int64)
    face_of = np.repeat(np.arange(nf, dtype=np.int64), sizes)
    pos = np.arange(len(nodes), dtype=np.int64)
    nxt = pos + 1
    if nf > 0:
        last = start[1:] - 1
        nxt[last[sizes > 0]] = start[:-1][sizes > 0]
    u = nodes
    v = nodes[nxt] if len(nodes) else nodes
    keep = u != v
    if nf > 0:
        keep &= sizes[face_of] >= 2
    u, v, face = u[keep], v[keep], face_of[keep]
    lo = np.minimum(u, v)
    hi = np.maximum(u, v)
    fwd = u < v
    order = np.lexsort((fwd, face, hi, lo))
    lo, hi, face, fwd = lo[order], hi[order], face[order], fwd[order]

    empty = {
        "lo": np.empty(0, dtype=np.int64),
        "hi": np.empty(0, dtype=np.int64),
        "uses": np.empty(0, dtype=np.int32),
        "inconsistent": np.empty(0, dtype=bool),
        "sharp": np.empty(0, dtype=bool),
        "angle": np.empty(0),
    }
    if len(lo) == 0:
        return empty
    head = np.ones(len(lo), dtype=bool)
    head[1:] = (lo[1:] != lo[:-1]) | (hi[1:] != hi[:-1])
    first = np.nonzero(head)[0]
    uses = np.diff(np.append(first, len(lo)))
    elo, ehi = lo[first], hi[first]

    inconsistent = np.zeros(len(first), dtype=bool)
    sharp = np.zeros(len(first), dtype=bool)
    angle = np.full(len(first), np.nan)
    report = uses != 2
    pair = np.nonzero(uses == 2)[0]
    a = first[pair]
    b = a + 1
    same_face = face[a] == face[b]
    drop = np.zeros(len(first), dtype=bool)
    drop[pair[same_face]] = True
    pair, a, b = pair[~same_face], a[~same_face], b[~same_face]
    inconsistent[pair] = fwd[a] == fwd[b]
    nx = normals[face[a]]
    ny = normals[face[b]]
    zx = (nx[:, 0] == 0.0) & (nx[:, 1] == 0.0) & (nx[:, 2] == 0.0)
    zy = (ny[:, 0] == 0.0) & (ny[:, 1] == 0.0) & (ny[:, 2] == 0.0)
    dot = nx[:, 0] * ny[:, 0] + nx[:, 1] * ny[:, 1] + nx[:, 2] * ny[:, 2]
    dot = np.where(inconsistent[pair], -dot, dot)
    defined = ~zx & ~zy
    cos_thr = np.cos(feature_angle * _PI / 180.0)
    sharp[pair] = defined & (dot < cos_thr)
    angle[pair] = np.where(
        defined, np.arccos(np.clip(dot, -1.0, 1.0)) * (180.0 / _PI), np.nan
    )
    report[pair] = inconsistent[pair] | sharp[pair]
    report &= ~drop
    return {
        "lo": elo[report],
        "hi": ehi[report],
        "uses": uses[report].astype(np.int32),
        "inconsistent": inconsistent[report],
        "sharp": sharp[report],
        "angle": angle[report],
    }


def crease_nodes(edges, n):
    """Boolean mask of the endpoints of every sharp or non-manifold edge."""
    mask = np.zeros(n, dtype=bool)
    crease = edges["sharp"] | (edges["uses"] >= 3)
    mask[edges["lo"][crease]] = True
    mask[edges["hi"][crease]] = True
    return mask


def _coords3(mesh):
    pts = np.asarray(mesh.points, dtype=np.float64)
    if pts.ndim == 1:
        pts = pts.reshape(-1, 1)
    xyz = np.zeros((len(pts), 3))
    xyz[:, : min(3, pts.shape[1])] = pts[:, :3]
    return xyz


def _selection(mesh, region, total):
    if not region:
        return np.ones(total, dtype=bool)
    mask = np.zeros(total, dtype=bool)
    found = False
    for r in getattr(mesh, "regions", ()) or ():
        if r.name == region and r.kind == "cell":
            e = np.asarray(r.entries, dtype=np.int64).ravel()
            mask[e[(e >= 0) & (e < total)]] = True
            found = True
            break
    if not found:
        names = sorted({r.name for r in getattr(mesh, "regions", ()) or ()})
        raise ValueError(
            f"meshio++: feature_edges: no cell region named '{region}' "
            f"(available: {', '.join(names) if names else 'none'})"
        )
    return mask


def _is_polyhedron(block):
    return block.type.startswith("polyhedron")


def _rings(mesh, selected):
    """``(start, nodes)`` of the examined surface's face rings."""
    from ._skin import _CELL_FACES
    from ._surface import _CELL_EDGES

    n = len(mesh.points)
    bases = block_bases(mesh.cells)
    rings = []

    volume = [b for b in mesh.cells if b.type in _CELL_FACES or _is_polyhedron(b)]
    if volume:
        if any(_is_polyhedron(b) for b in mesh.cells):
            raise NotImplementedError(
                "meshio++: feature_edges: the numpy fallback does not handle polyhedron "
                "cells; install a build with the compiled meshioplusplus._core extension"
            )
        keys, faces = [], []
        for base, block in zip(bases, mesh.cells):
            defs = _CELL_FACES.get(block.type)
            if not defs:
                continue
            data = np.asarray(block.data, dtype=np.int64)
            sel = selected[base : base + len(data)]
            data = data[sel]
            for _name, ncorner, local in defs:
                corners = data[:, list(local[:ncorner])]
                for row in corners:
                    faces.append(row)
                    keys.append(tuple(sorted(row.tolist())))
        counts = {}
        for k in keys:
            counts[k] = counts.get(k, 0) + 1
        rings = [f for f, k in zip(faces, keys) if counts[k] == 1]
    else:
        for base, block in zip(bases, mesh.cells):
            data = block.data
            if block.type.startswith("polygon"):
                rows = list(data)
            else:
                defs = _CELL_EDGES.get(block.type)
                if not defs or len(defs) < 3:
                    continue
                corners = [nodes[0] for _name, _k, nodes in defs]
                rows = list(np.asarray(data, dtype=np.int64)[:, corners])
            for c, row in enumerate(rows):
                if selected[base + c] and len(row) >= 3:
                    rings.append(np.asarray(row, dtype=np.int64))

    rings = [r for r in rings if len(r) and r.min() >= 0 and r.max() < n]
    start = np.zeros(len(rings) + 1, dtype=np.int64)
    if rings:
        start[1:] = np.cumsum([len(r) for r in rings])
        nodes = np.concatenate(rings).astype(np.int64)
    else:
        nodes = np.empty(0, dtype=np.int64)
    return start, nodes


def _feature_edges_py(mesh, angle, want, region):
    from ._common import warn
    from ._mesh import CellBlock, Mesh
    from ._regions import Region

    bases = block_bases(mesh.cells)
    total = int(bases[-1]) + len(mesh.cells[-1].data) if len(mesh.cells) else 0
    selected = _selection(mesh, region, total)
    start, nodes = _rings(mesh, selected)
    xyz = _coords3(mesh)
    edges = crease_edges(start, nodes, ring_unit_normals(xyz, start, nodes), angle)

    boundary = edges["uses"] == 1
    non_manifold = edges["uses"] >= 3
    kind = np.zeros(len(edges["lo"]), dtype=np.int32)
    for flag, mask, value in (
        (want["feature"], edges["sharp"], FEATURE),
        (want["inconsistent"], edges["inconsistent"], INCONSISTENT),
        (want["boundary"], boundary, BOUNDARY),
        (want["non_manifold"], non_manifold, NON_MANIFOLD),
    ):
        if flag:
            kind[mask] = value  # later assignments take precedence
    keep = kind != 0
    report = {
        "num_feature": int(edges["sharp"].sum()),
        "num_boundary": int(boundary.sum()),
        "num_non_manifold": int(non_manifold.sum()),
        "num_inconsistent": int(edges["inconsistent"].sum()),
    }
    conn = np.column_stack([edges["lo"][keep], edges["hi"][keep]]).astype(np.int64)
    out = Mesh(
        np.array(mesh.points, copy=True),
        [CellBlock("line", conn.reshape(-1, 2))],
        point_data={k: np.array(v, copy=True) for k, v in mesh.point_data.items()},
        cell_data={
            "feature:kind": [kind[keep]],
            "feature:angle": [edges["angle"][keep].astype(np.float64)],
        },
        field_data={k: np.array(v, copy=True) for k, v in mesh.field_data.items()},
    )
    regions = getattr(mesh, "regions", None) or []
    kept = [r.copy() for r in regions if r.kind == "point"]
    if len(kept) != len(regions):
        warn(
            "feature_edges: Cell and Side regions name input cells, which the edge mesh "
            "does not have; they are dropped (Point regions are kept)"
        )
    out.regions = [Region(r.name, r.kind, r.entries, r.dim, r.tag) for r in kept]
    return out, report


def feature_edges(
    mesh,
    feature_angle: float = 30.0,
    feature: bool = True,
    boundary: bool = True,
    non_manifold: bool = True,
    inconsistent: bool = True,
    region: str = "",
    return_report: bool = False,
):
    """The feature edges of a surface (or of a volume mesh's skin), as ``line`` cells.

    :param mesh: a surface mesh (triangles, quads, polygons and their quadratic
        variants, which contribute their corners) or a volume mesh, whose skin is
        examined; lines and vertices are ignored.
    :param feature_angle: the largest dihedral angle, in degrees, still treated
        as smooth; must lie in ``[0, 180]``.
    :param feature: report edges sharper than ``feature_angle``.
    :param boundary: report open edges (used by one face).
    :param non_manifold: report edges used by three or more faces.
    :param inconsistent: report face pairs that walk their shared edge the same
        way (disagree about which side is out).
    :param region: restrict to the cells of this named cell region.
    :param return_report: also return the per-category edge counts.
    :returns: a mesh with the input's points, one ``line`` block (ascending in
        its endpoint ids) and cell data ``feature:kind`` (1 feature, 2 boundary,
        3 non-manifold, 4 inconsistent -- the first selected of non-manifold,
        boundary, inconsistent, feature) and ``feature:angle`` (degrees; NaN
        where undefined). With ``return_report``, ``(mesh, report)``.
    """
    angle = float(feature_angle)
    if not (0.0 <= angle <= 180.0):
        raise ValueError(
            "meshio++: feature_edges: the feature angle must lie in [0, 180] degrees"
        )
    want = {
        "feature": bool(feature),
        "boundary": bool(boundary),
        "non_manifold": bool(non_manifold),
        "inconsistent": bool(inconsistent),
    }
    out = None
    try:
        from . import _core

        res = _core.feature_edges(
            mesh,
            angle,
            want["feature"],
            want["boundary"],
            want["non_manifold"],
            want["inconsistent"],
            str(region),
        )
        out = res["mesh"]
        report = {
            k: int(res[k])
            for k in (
                "num_feature",
                "num_boundary",
                "num_non_manifold",
                "num_inconsistent",
            )
        }
    except Exception as exc:
        if not core_op_declined(exc, "feature_edges"):
            raise
        out = None
    if out is None:
        out, report = _feature_edges_py(mesh, angle, want, str(region))
    return (out, report) if return_report else out


__all__ = ["feature_edges"]
