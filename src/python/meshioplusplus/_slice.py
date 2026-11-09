"""Planar cross-section of a mesh (the intersection with a plane).

A dependency-free mesh *operation* (not a file format): it returns a new mesh
one topological dimension below the cut cells -- a 3D volume mesh cut by a plane
yields a ``triangle``/``quad`` surface, a 2D surface mesh yields a ``line`` mesh.
Unlike :func:`crop` (plane mode), which keeps whole cells on one side, ``slice``
computes the actual intersection and lowers the dimension.

The C++ core (``_core.slice``) does the work; this module is the thin shim
(try C++, fall back to a pure-numpy reference) plus the plane arithmetic. The
cutter itself -- marching tetrahedra on a simplexified input, the watertight
edge dedup, the winding and the degeneracy rule -- lives in :mod:`._marching`,
shared with :mod:`._isosurface` exactly as the C++ side shares
``detail/marching.hpp``; it reproduces the C++ result **byte for byte**, pinned
by ``tests/python/test_slice.py::test_cpp_matches_python``.

``slice`` shadows the Python built-in **only as a module attribute** -- the name
``meshioplusplus.slice`` is deliberate; the built-in is never rebound here.

Public API:

* :func:`slice` -- the planar cross-section of a mesh.
* :func:`sections` -- chained sections at parallel planes (roadmap 7.6.1).
"""

from __future__ import annotations

import numpy as np

from ._fallback import core_op_declined
from ._marching import FIXED_DIRECTION, _marching_cut, _marching_prepare


def _slice_py(mesh, origin, normal, record_parent_ids):
    origin = np.asarray(origin, dtype=np.float64).reshape(3)
    normal = np.asarray(normal, dtype=np.float64).reshape(3)
    if float(normal @ normal) <= 0.0:
        raise ValueError("meshio++: slice: normal must be non-zero")

    prep = _marching_prepare(mesh)

    # Signed distance of every simplexified node to the plane (padded z=0 in 2D).
    pts = prep.points
    dim = prep.dim
    n_in = len(pts)
    pts3 = np.zeros((n_in, 3), dtype=np.float64)
    if n_in:
        pts3[:, :dim] = pts
    dist = (
        (pts3[:, 0] - origin[0]) * normal[0]
        + (pts3[:, 1] - origin[1]) * normal[1]
        + (pts3[:, 2] - origin[2]) * normal[2]
    )

    return _marching_cut(
        prep,
        dist,
        "slice:parent_cell",
        record_parent_ids,
        orientation=FIXED_DIRECTION,
        direction=normal,
    )


def _polyline_chains(segs, npts):
    """``detail::polyline_chains``, step for step: a list of
    ``(points, segments, reversed, closed)`` in the canonical order."""
    adj = [[] for _ in range(npts)]
    for s, (a, b) in enumerate(segs):
        if a == b or a < 0 or b < 0 or a >= npts or b >= npts:
            continue
        adj[a].append((b, s))
        adj[b].append((a, s))
    for lst in adj:
        lst.sort()
    used = [False] * len(segs)

    def next_edge(p):
        for e in adj[p]:
            if not used[e[1]]:
                return e
        return None

    def walk(start, loop):
        pts, sg, rev = [start], [], []
        p = start
        while True:
            e = next_edge(p)
            if e is None:
                break
            q, s = e
            used[s] = True
            sg.append(s)
            rev.append(0 if segs[s][0] == p else 1)
            pts.append(q)
            p = q
            if (p == start) if loop else (len(adj[p]) != 2):
                break
        return pts, sg, rev, loop

    out = []
    for p in range(npts):
        if len(adj[p]) != 2 and adj[p]:
            while next_edge(p) is not None:
                out.append(walk(p, False))
    for p in range(npts):
        while next_edge(p) is not None:
            pts, sg, rev, _ = walk(p, True)
            out.append((pts, sg, rev, pts[-1] == pts[0]))
    return out


def _chain_py(section, min_length):
    """``slc_chain``: the line section walked into chains."""
    from ._mesh import CellBlock, Mesh

    if not section.cells:
        return section
    if len(section.cells) != 1 or section.cells[0].type != "line":
        raise ValueError(
            "meshio++: slice: chain needs a line section; slice a surface mesh (a "
            "volume's section is a surface)"
        )
    conn = np.asarray(section.cells[0].data, dtype=np.int64)
    pts = np.asarray(section.points)
    npts = len(pts)
    # Points at exactly one position are one point of a chain, the lowest id
    # standing in (slc_chain's rule).
    first, same = {}, []
    for i in range(npts):
        key = tuple(float(v) for v in pts[i])
        same.append(first.setdefault(key, i))
    chains = _polyline_chains([(same[int(a)], same[int(b)]) for a, b in conn], npts)
    new_of = np.full(npts, -1, dtype=np.int64)
    old_points, seg_order, conn_out, chain_id, closed = [], [], [], [], []
    kept = 0
    for cp, cs, _rev, cl in chains:
        length = 0.0
        for k in range(len(cp) - 1):
            d2 = 0.0
            for a in range(pts.shape[1]):
                d = float(pts[cp[k + 1], a]) - float(pts[cp[k], a])
                d2 += d * d
            length += float(np.sqrt(d2))
        if length < min_length:
            continue
        for p in cp:
            if new_of[p] < 0:
                new_of[p] = len(old_points)
                old_points.append(p)
        for k, s in enumerate(cs):
            seg_order.append(s)
            conn_out.append((new_of[cp[k]], new_of[cp[k + 1]]))
            chain_id.append(kept)
            closed.append(1 if cl else 0)
        kept += 1
    old = np.asarray(old_points, dtype=np.int64)
    order = np.asarray(seg_order, dtype=np.int64)
    out = Mesh(
        pts[old],
        [CellBlock("line", np.asarray(conn_out, dtype=np.int64).reshape(-1, 2))],
        point_data={k: np.asarray(v)[old] for k, v in section.point_data.items()},
        cell_data={k: [np.asarray(v[0])[order]] for k, v in section.cell_data.items()},
    )
    out.cell_data["section:chain"] = [np.asarray(chain_id, dtype=np.int64)]
    out.cell_data["section:closed"] = [np.asarray(closed, dtype=np.int64)]
    return out


def slice(
    mesh,
    origin=(0.0, 0.0, 0.0),
    normal=(0.0, 0.0, 1.0),
    record_parent_ids: bool = False,
    *,
    chain: bool = False,
    min_length: float = 0.0,
):
    """Compute the planar cross-section of a mesh.

    Parameters
    ----------
    mesh :
        the mesh to cut (unmodified).
    origin :
        a point on the cutting plane (3 floats).
    normal :
        the plane normal (3 floats; non-zero, need not be unit length).
    record_parent_ids :
        when true, attach an Int64 ``slice:parent_cell`` ``cell_data`` array
        recording, per section cell, the global (block-major) index of the input
        cell it was cut from.
    chain :
        walk a surface's line section into polylines (roadmap 7.6.1): open
        chains from their lowest-id end, then loops from their lowest point,
        segments emitted in walking order with Int64 ``section:chain`` and
        ``section:closed`` cell data, points renumbered in visiting order. A
        volume's section (a surface) is refused.
    min_length :
        with ``chain``, drop chains shorter than this.

    Returns
    -------
    Mesh
        the cross-section: a ``triangle``/``quad`` surface for a volume mesh, a
        ``line`` mesh for a 2D surface mesh, or an empty mesh when the plane
        misses the geometry. Section points are all new, so
        ``point_sets``/``cell_sets`` are **not** carried; each section cell
        inherits its parent cell's ``cell_data``, and interpolated
        ``point_data`` is promoted to Float64.

    See Also
    --------
    meshioplusplus.isosurface :
        the data-driven sibling -- the level set of a scalar field rather than
        of the distance to a plane, through the same cutter.
    """
    origin = [float(v) for v in np.asarray(origin, dtype=np.float64).reshape(3)]
    normal = [float(v) for v in np.asarray(normal, dtype=np.float64).reshape(3)]
    min_length = float(min_length)
    if not (np.isfinite(min_length) and min_length >= 0.0):
        raise ValueError("meshio++: slice: min_length must be finite and >= 0")

    out = None
    try:
        from . import _core

        out = _core.slice(
            mesh,
            origin,
            normal,
            bool(record_parent_ids),
            chain=bool(chain),
            min_length=min_length,
        )
    except Exception as exc:
        if not core_op_declined(exc, "slice"):
            raise
        out = None
    if out is None:
        out = _slice_py(mesh, origin, normal, record_parent_ids)
        if chain:
            out = _chain_py(out, min_length)
    return out


def sections(
    mesh,
    direction=(0.0, 0.0, 1.0),
    offsets=(0.0,),
    *,
    min_length: float = 0.0,
    record_parent_ids: bool = False,
):
    """Chained sections of a surface at parallel planes, merged into one line
    mesh (roadmap 7.6.1): plane ``k`` is ``n . x = offsets[k]`` with ``n`` the
    unit ``direction``; each segment carries ``section:offset`` (Float64),
    ``section:chain`` (numbered across all planes) and ``section:closed``.

    :param mesh: a surface mesh (a volume's sections are surfaces: refused).
    :param direction: the planes' normal, non-zero.
    :param offsets: the planes' signed distances from the origin.
    :param min_length: drop chains shorter than this.
    :param record_parent_ids: as :func:`slice`.
    :returns: the merged line mesh (empty when no plane meets the surface).
    """
    from ._merge import merge
    from ._mesh import Mesh

    d = np.asarray(direction, dtype=np.float64).reshape(3)
    length = float(np.sqrt(d @ d))
    if not (length > 0.0 and np.isfinite(length)):
        raise ValueError("meshio++: sections: direction must be non-zero")
    offs = [float(v) for v in np.asarray(offsets, dtype=np.float64).ravel()]
    if not offs:
        raise ValueError("meshio++: sections: give at least one offset")
    if not all(np.isfinite(offs)):
        raise ValueError("meshio++: sections: offsets must be finite")
    try:
        from . import _core

        return _core.sections(
            mesh,
            [float(v) for v in d],
            offs,
            float(min_length),
            bool(record_parent_ids),
        )
    except Exception as exc:
        if not core_op_declined(exc, "sections"):
            raise
    n = d / length
    parts, base = [], 0
    for off in offs:
        part = slice(
            mesh,
            n * off,
            n,
            record_parent_ids,
            chain=True,
            min_length=min_length,
        )
        if not part.cells or len(part.cells[0].data) == 0:
            continue
        chains = np.asarray(part.cell_data["section:chain"][0], dtype=np.int64)
        part.cell_data["section:chain"] = [chains + base]
        base += int(chains.max()) + 1
        part.cell_data["section:offset"] = [np.full(len(chains), off, dtype=np.float64)]
        parts.append(part)
    if not parts:
        return Mesh(np.zeros((0, 3)), [])
    if len(parts) == 1:
        return parts[0]
    return merge(parts, source_tag=False)
