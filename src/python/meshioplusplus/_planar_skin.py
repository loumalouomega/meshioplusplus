"""The planar skin shared by ``voxelize`` (pixels) and ``octree`` (quadtrees).

A planar polyline from ``line`` cells, or the boundary edges of ``triangle`` /
``quad`` cells, with the segment/box overlap and signed-distance kernels that go
with it. The numpy twin of ``detail/planar_skin.hpp``.
"""

from __future__ import annotations

import numpy as np


def seg_box(centre, half, a, b):
    """Segment / box overlap in the xy plane; touching counts."""
    ax, ay = float(a[0] - centre[0]), float(a[1] - centre[1])
    bx, by = float(b[0] - centre[0]), float(b[1] - centre[1])
    if min(ax, bx) > half[0] or max(ax, bx) < -half[0]:
        return False
    if min(ay, by) > half[1] or max(ay, by) < -half[1]:
        return False
    dx, dy = bx - ax, by - ay
    nx, ny = -dy, dx
    proj = nx * ax + ny * ay
    return abs(proj) <= abs(nx) * half[0] + abs(ny) * half[1]


_PLANES = {"xy": (0, 1, 2), "xz": (0, 2, 1), "yz": (1, 2, 0)}
PLANES = ("auto",) + tuple(_PLANES)


def frame(mesh, plane="auto"):
    """The axis permutation that maps the skin's plane to xy (in-plane axes first).

    ``auto`` picks the axis along which the points have zero extent (z, then y,
    then x), else z. The numpy twin of ``detail::planar_frame``.
    """
    if plane not in PLANES:
        raise ValueError(
            f"meshio++: unknown plane '{plane}' (expected one of: {', '.join(PLANES)})"
        )
    if plane == "auto":
        pts = np.asarray(mesh.points, dtype=np.float64)
        plane = "xy"
        if len(pts):
            if pts.shape[1] == 2:
                pts = np.column_stack([pts, np.zeros(len(pts))])
            lo, hi = pts.min(axis=0), pts.max(axis=0)
            for k, name in ((2, "xy"), (1, "xz"), (0, "yz")):
                if hi[k] == lo[k]:
                    plane = name
                    break
    return _PLANES[plane]


def permute(arr, axes):
    """World -> frame coordinates of an ``(..., 3)`` array."""
    return np.asarray(arr)[..., list(axes)]


def unpermute(arr, axes):
    """Frame -> world coordinates of an ``(..., 3)`` array."""
    arr = np.asarray(arr)
    out = np.empty_like(arr)
    out[..., list(axes)] = arr
    return out


def collect(mesh, prefix, axes=(0, 1, 2)):
    points = np.asarray(mesh.points, dtype=np.float64)
    if points.shape[1] == 2:
        points = np.column_stack([points, np.zeros(len(points))])
    points = permute(points, axes)
    ids = []
    boundary = {}
    order = []
    has_lines = False
    for cb in mesh.cells:
        if cb.type == "line":
            kind = "line"
        elif cb.type in ("triangle", "quad"):
            kind = "surf"
        else:
            continue
        data = np.asarray(cb.data, dtype=np.int64)
        for row in data:
            if row.min() < 0 or row.max() >= len(points):
                raise ValueError(
                    f"{prefix}a cell references a point the mesh does not have"
                )
            if kind == "line":
                has_lines = True
                ids.append((int(row[0]), int(row[1])))
                continue
            n = len(row)
            for e in range(n):
                a, b = int(row[e]), int(row[(e + 1) % n])
                key = (min(a, b), max(a, b))
                if key not in boundary:
                    boundary[key] = [1, (a, b)]
                    order.append(key)
                else:
                    boundary[key][0] += 1
    for key in order:
        if boundary[key][0] == 1:
            ids.append(boundary[key][1])
    _orient_loops(ids)
    prims = [(points[a], points[b], points[b]) for a, b in ids]
    return prims, ids, has_lines


def _orient_loops(ids):
    """Orient every closed loop like its first segment; the twin of ``ps_orient_loops``."""
    incident = {}
    for i, (a, b) in enumerate(ids):
        incident.setdefault(a, []).append(i)
        if b != a:
            incident.setdefault(b, []).append(i)
    seen = [False] * len(ids)
    for first in range(len(ids)):
        if seen[first]:
            continue
        seen[first] = True
        if ids[first][0] == ids[first][1]:
            continue
        flips = []
        closed = False
        cur, at = first, ids[first][1]
        while True:
            inc = incident[at]
            if len(inc) != 2:
                break
            nxt = inc[1] if inc[0] == cur else inc[0]
            if nxt == first:
                closed = True
                break
            if seen[nxt]:
                break
            seen[nxt] = True
            if ids[nxt][0] != at:
                flips.append(nxt)
            at = ids[nxt][1] if ids[nxt][0] == at else ids[nxt][0]
            cur = nxt
        if closed:
            for f in flips:
                ids[f] = (ids[f][1], ids[f][0])


def require_closed(ids, prefix):
    degree = {}
    for a, b in ids:
        degree[a] = degree.get(a, 0) + 1
        degree[b] = degree.get(b, 0) + 1
    odd = sorted(v for v, d in degree.items() if d % 2 != 0)
    if odd:
        shown = ", ".join(str(v) for v in odd[:8]) + (", ..." if len(odd) > 8 else "")
        raise ValueError(
            f"{prefix}the polyline does not close, so it has no inside: "
            f"open ends at points {shown}"
        )


def signed_distance(prims, p, even_odd):
    best = float("inf")
    winding = 0
    crossings = 0
    for a, b, _c in prims:
        ax, ay, bx, by = float(a[0]), float(a[1]), float(b[0]), float(b[1])
        dx, dy = bx - ax, by - ay
        len2 = dx * dx + dy * dy
        t = ((p[0] - ax) * dx + (p[1] - ay) * dy) / len2 if len2 > 0.0 else 0.0
        t = 0.0 if t < 0.0 else (1.0 if t > 1.0 else t)
        qx = ax + t * dx - p[0]
        qy = ay + t * dy - p[1]
        d = float(np.sqrt(qx * qx + qy * qy))
        if d < best:
            best = d
        left = dx * (p[1] - ay) - (p[0] - ax) * dy
        if ay <= p[1]:
            if by > p[1] and left > 0.0:
                winding += 1
                crossings += 1
        elif by <= p[1] and left < 0.0:
            winding -= 1
            crossings += 1
    inside = (crossings % 2 != 0) if even_odd else (winding != 0)
    return -best if inside else best


def polyline_only(mesh):
    """True when the mesh has ``line`` cells and nothing but ``line``/``vertex``."""
    types = [cb.type for cb in mesh.cells]
    return "line" in types and all(t in ("line", "vertex") for t in types)
