"""Periodic node pairs: which master node each slave node maps onto.

``match_periodic_nodes(mesh, slave, master, translate=..., rotate=...)`` maps
every node of the *slave* region through an affine transform and finds the
*master* node it lands on within ``atol`` -- the pairs a periodic boundary
condition ties together (a Kratos periodic condition, a Gmsh ``$Periodic``
section). The contract is ``operations/periodic.hpp``'s: nearest master wins,
ties to the lower id; a master claimed twice is an error; slave nodes that map
onto themselves and are also master nodes (on a rotation axis) are *fixed
points*, counted and not paired.

The transform uses :func:`meshioplusplus.transform`'s vocabulary --
``translate``, ``rotate=(axis, degrees)``, ``matrix`` -- plus ``origin``, the
point a rotation turns about.
"""

from __future__ import annotations

import numpy as np

from ._fallback import core_op_declined
from ._transform import _build_matrix, _translation

_PREFIX = "meshio++: match_periodic_nodes: "


def periodic_matrix(translate=None, rotate=None, matrix=None, origin=None):
    """The row-major 4x4 slave-to-master matrix the keywords describe.

    ``rotate`` turns about ``origin`` (default the coordinate origin), then
    ``translate`` applies; ``matrix`` overrides both.
    """
    if matrix is None and translate is None and rotate is None:
        raise ValueError(f"{_PREFIX}give a transform: translate, rotate or matrix")
    if matrix is not None:
        return _build_matrix(None, None, None, matrix, None)
    m = np.eye(4)
    if rotate is not None:
        o = (
            np.zeros(3)
            if origin is None
            else np.asarray(origin, dtype=np.float64).reshape(3)
        )
        m = (
            _translation(*o)
            @ _build_matrix(None, None, rotate, None, None)
            @ _translation(*(-o))
        )
    if translate is not None:
        m = _build_matrix(translate, None, None, None, None) @ m
    return m


def _selector(spec):
    if isinstance(spec, str):
        return {"name": spec}
    sel = {"name": str(spec["name"])}
    for key in ("kind", "dim", "tag"):
        if spec.get(key) is not None:
            sel[key] = spec[key]
    return sel


def _region_points(mesh, region):
    from ._facets import facet_nodes
    from ._regions import block_bases

    n = len(mesh.points)
    e = np.asarray(region.entries, dtype=np.int64)
    if region.kind == "point":
        ids = e.ravel()
    elif region.kind == "cell":
        bases = block_bases(mesh.cells)
        out = []
        for c in e.ravel():
            b = int(np.searchsorted(bases, c, side="right") - 1)
            if b < 0 or b >= len(mesh.cells) or c - bases[b] >= len(mesh.cells[b].data):
                continue
            out.extend(
                int(x) for x in np.asarray(mesh.cells[b].data[c - bases[b]]).ravel()
            )
        ids = np.asarray(out, dtype=np.int64)
    else:
        out = []
        for cell, facet in e.reshape(-1, 2):
            hit = facet_nodes(mesh, int(cell), int(facet))
            if hit is not None:
                out.extend(hit[1])
        ids = np.asarray(out, dtype=np.int64)
    ids = ids[(ids >= 0) & (ids < n)]
    return np.unique(ids)


def _match_py(mesh, slave, master, m, atol, require_complete):
    from ._region_ops import _find

    regions = getattr(mesh, "regions", None) or []
    s_ids = _region_points(mesh, regions[_find(regions, slave)])
    m_ids = _region_points(mesh, regions[_find(regions, master)])
    pts = np.asarray(mesh.points, dtype=np.float64)
    xyz = np.zeros((len(pts), 3))
    xyz[:, : min(3, pts.shape[1])] = pts[:, :3]
    p = xyz[s_ids]
    q = p @ m[:3, :3].T + m[:3, 3]
    mp = xyz[m_ids]
    match = np.full(len(s_ids), -1, dtype=np.int64)
    resid = np.zeros(len(s_ids))
    for i in range(len(s_ids)):
        d2 = np.sum((mp - q[i]) ** 2, axis=1)
        ok = np.nonzero(d2 <= atol * atol)[0]
        if len(ok):
            k = ok[np.argmin(d2[ok])]  # first minimum: the lowest master id
            match[i] = k
            resid[i] = np.sqrt(d2[k])
    pairs, unmatched, fixed, claimed, max_res = [], [], 0, {}, 0.0
    for i, k in enumerate(match):
        s = int(s_ids[i])
        if k < 0:
            unmatched.append(s)
            continue
        mm = int(m_ids[k])
        if mm == s:
            fixed += 1
            continue
        if mm in claimed:
            raise ValueError(
                f"{_PREFIX}slave nodes {claimed[mm]} and {s} both map onto master node "
                f"{mm} (the tolerance is too loose, or the regions overlap)"
            )
        claimed[mm] = s
        pairs.append((s, mm))
        max_res = max(max_res, float(resid[i]))
    if unmatched and require_complete:
        first = ", ".join(str(x) for x in unmatched[:8]) + (
            ", ..." if len(unmatched) > 8 else ""
        )
        raise ValueError(
            f"{_PREFIX}{len(unmatched)} slave node(s) have no master node within {atol:g} "
            f"({first}); check the transform, or pass require_complete=False"
        )
    pairs = np.asarray(pairs, dtype=np.int64).reshape(-1, 2)
    return {
        "slave": pairs[:, 0].copy(),
        "master": pairs[:, 1].copy(),
        "unmatched": np.asarray(unmatched, dtype=np.int64),
        "num_fixed": fixed,
        "max_residual": max_res,
    }


def match_periodic_nodes(
    mesh,
    slave,
    master,
    translate=None,
    rotate=None,
    matrix=None,
    origin=None,
    atol: float = 1e-8,
    require_complete: bool = True,
    return_report: bool = False,
):
    """The master node each node of the ``slave`` region maps onto.

    :param mesh: the mesh carrying both regions.
    :param slave: the region whose nodes are mapped: a name, or a
        ``{"name", "kind", "dim", "tag"}`` selector matching exactly one region.
        Point regions contribute their entries, Cell regions their cells' nodes,
        Side regions their facets' nodes.
    :param master: the region they map onto, likewise.
    :param translate: a ``(dx, dy, dz)`` slave-to-master translation.
    :param rotate: ``(axis, degrees)`` -- ``axis`` is ``"x"``/``"y"``/``"z"`` or a
        3-vector -- turning about ``origin``; applied before ``translate``.
    :param matrix: a row-major 4x4 (or flat 16) affine matrix, overriding both.
    :param origin: the point ``rotate`` turns about (default the origin).
    :param atol: the largest distance, after the transform, at which two nodes
        match; must be positive.
    :param require_complete: raise when a slave node has no master.
    :param return_report: also return ``{"unmatched", "num_fixed",
        "max_residual"}``.
    :returns: an ``(k, 2)`` int64 array of ``[slave, master]`` node pairs,
        ascending in the slave id; with ``return_report``, ``(pairs, report)``.
    """
    m = periodic_matrix(translate, rotate, matrix, origin)
    if not (float(atol) > 0.0) or not np.isfinite(float(atol)):
        raise ValueError(f"{_PREFIX}the tolerance must be positive")
    s_sel, m_sel = _selector(slave), _selector(master)
    res = None
    try:
        from . import _core

        res = _core.match_periodic_nodes(
            mesh,
            s_sel,
            m_sel,
            [float(x) for x in m.reshape(-1)],
            float(atol),
            bool(require_complete),
        )
    except Exception as exc:
        if not core_op_declined(exc, "match_periodic_nodes"):
            raise
        res = None
    if res is None:
        res = _match_py(mesh, s_sel, m_sel, m, float(atol), bool(require_complete))
    pairs = np.column_stack(
        [
            np.asarray(res["slave"], dtype=np.int64),
            np.asarray(res["master"], dtype=np.int64),
        ]
    ).reshape(-1, 2)
    if not return_report:
        return pairs
    report = {
        "unmatched": np.asarray(res["unmatched"], dtype=np.int64),
        "num_fixed": int(res["num_fixed"]),
        "max_residual": float(res["max_residual"]),
    }
    return pairs, report


__all__ = ["match_periodic_nodes"]
