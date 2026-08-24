"""Mesh smoothing (``meshioplusplus.smooth``).

A dependency-free mesh *operation* (not a file format): it relaxes point
coordinates toward their edge-neighbour centroids, improving element shape
without touching connectivity, topology or any data value. Only the point
coordinates change.

The C++ core (``_core.smooth``) is used when available and this module's pure
numpy reference (``_smooth_py``) is the fallback. Both implement the same two
operators -- Laplacian, which smooths strongly but shrinks, and Taubin, which
alternates a positive and a slightly larger negative pass and so does not.

``point_sets``/``cell_sets`` pass through unchanged: this operation never adds,
removes or renumbers a node or a cell, so there is nothing to remap -- unlike
``refine``/``crop``/``split``, whose shims must rebuild them. They are copied
here anyway because the core cannot see them across the Python boundary.
"""

import numpy as np

from ._mesh import Mesh
from ._skin import _CELL_FACES
from ._surface import _CELL_EDGES

# Interior edge topology, mirroring cell_refine_edges in
# src/cpp/src/detail/cell_subdivision.cpp. KEEP IN SYNC with that table -- the two
# must agree or the C++ core and this fallback will smooth differently.
#
# These are the seven linear types. Every other type (the higher-order family,
# the VTK-Lagrange types, ``custom``) has no known edge topology, so its nodes
# are pinned rather than guessed at; see the module docs of
# src/cpp/include/meshioplusplus/detail/node_adjacency.hpp for why.
_REFINE_EDGES = {
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
    "wedge": [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)],
    "pyramid": [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 4), (2, 4), (3, 4)],
}

_METHODS = ("laplacian", "taubin")


def _coords3(mesh):
    """Points as an (n, 3) float64 array, z-padded for a 2D mesh."""
    pts = np.asarray(mesh.points, dtype=np.float64)
    if pts.ndim != 2:
        raise ValueError("meshio++: smooth: points must be a 2-D array")
    n, dim = pts.shape
    xyz = np.zeros((n, 3), dtype=np.float64)
    xyz[:, : min(dim, 3)] = pts[:, : min(dim, 3)]
    return xyz, dim


def _edge_graph(mesh, n):
    """CSR edge adjacency with ascending, deduplicated rows.

    Mirrors ``detail::build_node_adjacency(..., NodeAdjacencyKind::Edge)``.
    """
    src = []
    dst = []

    def add(a, b):
        src.append(np.asarray(a, dtype=np.int64).ravel())
        dst.append(np.asarray(b, dtype=np.int64).ravel())

    for block in mesh.cells:
        data = block.data
        if (
            block.type in _REFINE_EDGES
            and isinstance(data, np.ndarray)
            and data.ndim == 2
        ):
            for a, b in _REFINE_EDGES[block.type]:
                add(data[:, a], data[:, b])
                add(data[:, b], data[:, a])
        elif not isinstance(data, np.ndarray) or data.ndim != 2:
            # A ragged polygon block: each row is a closed ring of nodes.
            for row in data:
                row = np.asarray(row, dtype=np.int64).ravel()
                if row.size >= 2:
                    add(row, np.roll(row, -1))
                    add(np.roll(row, -1), row)
        # else: unknown edge topology -- contributes nothing; pinned below.

    if not src:
        return np.zeros(n + 1, dtype=np.int64), np.zeros(0, dtype=np.int64)

    a = np.concatenate(src)
    b = np.concatenate(dst)
    keep = (a >= 0) & (a < n) & (b >= 0) & (b < n) & (a != b)
    a, b = a[keep], b[keep]
    order = np.lexsort((b, a))
    a, b = a[order], b[order]
    if a.size:
        uniq = np.ones(a.size, dtype=bool)
        uniq[1:] = (a[1:] != a[:-1]) | (b[1:] != b[:-1])
        a, b = a[uniq], b[uniq]
    xadj = np.zeros(n + 1, dtype=np.int64)
    xadj[1:] = np.cumsum(np.bincount(a, minlength=n))
    return xadj, b


def _has_volume_cells(mesh):
    return any(block.type in _CELL_FACES for block in mesh.cells)


def _boundary_and_features(mesh, xyz, n, face_mode, preserve_features, feature_angle):
    """Boundary-node mask and, optionally, the feature-node mask."""
    table = _CELL_FACES if face_mode else _CELL_EDGES
    keys = []
    facets = []  # (corner ids array,) per facet, aligned with keys
    for block in mesh.cells:
        defs = table.get(block.type)
        data = block.data
        if not defs or not isinstance(data, np.ndarray) or data.ndim != 2:
            continue
        for _name, ncorner, nodes in defs:
            corners = data[:, list(nodes[:ncorner])]
            key = np.sort(corners, axis=1)
            pad = np.full((key.shape[0], 4 - key.shape[1]), -1, dtype=key.dtype)
            keys.append(np.hstack([pad, key]) if pad.size else key)
            facets.append(corners)
    boundary = np.zeros(n, dtype=bool)
    if not keys:
        return boundary, np.zeros(n, dtype=bool)

    allkeys = np.vstack(keys)
    allfacets = (
        np.vstack([f for f in facets])
        if len({f.shape[1] for f in facets}) == 1
        else None
    )
    _uniq, inverse, counts = np.unique(
        allkeys, axis=0, return_inverse=True, return_counts=True
    )
    once = counts[inverse] == 1

    # Mark boundary nodes.
    offset = 0
    bfacets = []
    for corners in facets:
        m = once[offset : offset + corners.shape[0]]
        offset += corners.shape[0]
        sel = corners[m]
        if sel.size:
            boundary[sel.ravel()] = True
            bfacets.append(sel)
    del allfacets

    features = np.zeros(n, dtype=bool)
    if not preserve_features or not bfacets:
        return boundary, features

    # Facet normals: Newell for a face, the in-plane perpendicular for a 2D edge.
    normals = []
    nodes_per = []
    for sel in bfacets:
        k = sel.shape[1]
        if k >= 3:
            p = xyz[sel]
            nrm = np.zeros((sel.shape[0], 3), dtype=np.float64)
            for i in range(k):
                nrm += np.cross(p[:, i], p[:, (i + 1) % k])
        else:
            p0, p1 = xyz[sel[:, 0]], xyz[sel[:, 1]]
            nrm = np.column_stack(
                [p1[:, 1] - p0[:, 1], p0[:, 0] - p1[:, 0], np.zeros(sel.shape[0])]
            )
        length = np.linalg.norm(nrm, axis=1, keepdims=True)
        nrm = np.divide(nrm, length, out=np.zeros_like(nrm), where=length > 0)
        normals.append(nrm)
        nodes_per.append(sel)

    allnorm = np.vstack(normals)
    node_ids = np.concatenate([s.ravel() for s in nodes_per])
    facet_ids = np.concatenate(
        [
            np.repeat(np.arange(s.shape[0]) + off, s.shape[1])
            for off, s in zip(
                np.cumsum([0] + [s.shape[0] for s in nodes_per[:-1]]), nodes_per
            )
        ]
    )
    cos_thr = np.cos(np.deg2rad(feature_angle))
    order = np.argsort(node_ids, kind="stable")
    node_ids, facet_ids = node_ids[order], facet_ids[order]
    starts = np.searchsorted(node_ids, np.arange(n), side="left")
    ends = np.searchsorted(node_ids, np.arange(n), side="right")
    for i in range(n):
        s, e = starts[i], ends[i]
        if e - s < 2:
            continue
        nn = allnorm[facet_ids[s:e]]
        if (nn @ nn.T).min() < cos_thr:
            features[i] = True
    return boundary, features


def _pinned_unknown_topology(mesh, n):
    """Nodes of blocks whose edge topology is unknown are pinned."""
    pinned = np.zeros(n, dtype=bool)
    for block in mesh.cells:
        data = block.data
        ragged = not isinstance(data, np.ndarray) or data.ndim != 2
        if ragged or block.type in _REFINE_EDGES:
            continue
        ids = np.asarray(data, dtype=np.int64).ravel()
        ids = ids[(ids >= 0) & (ids < n)]
        pinned[ids] = True
    return pinned


def _smooth_py(
    mesh,
    method="taubin",
    iterations=10,
    lambda_=-1.0,
    mu=-0.34,
    fix_boundary=True,
    preserve_features=True,
    feature_angle=30.0,
    guard_inversion=True,
    frozen=None,
):
    """Pure-numpy reference implementation of :func:`smooth`.

    The inversion guard is **not** implemented here: it is a discrete branch on
    the sign of a cell measure, and for a cell within round-off of degenerate
    this reference and the C++ core could land on opposite sides and then
    diverge macroscopically rather than by round-off. Rather than ship a second,
    subtly different guard, the fallback raises when one is requested.
    """
    if method == "odt":
        # Unlike laplacian/taubin, ODT is C++-core only unconditionally, not
        # merely under guard_inversion=True: an unguarded ODT step can invert
        # a tet near a concave boundary, so the method is only meaningful
        # guarded, and the guard itself is the same discrete-branch-on-a-sign
        # this reference already declines to duplicate above.
        raise NotImplementedError(
            "meshio++: smooth: method='odt' is C++-core only and has no "
            "pure-numpy fallback -- install a build with the compiled "
            "meshioplusplus._core extension"
        )
    if method not in _METHODS:
        raise ValueError(
            f"meshio++: smooth: unknown method {method!r} "
            f"(expected 'laplacian', 'taubin' or 'odt')"
        )
    if guard_inversion:
        raise NotImplementedError(
            "meshio++: smooth: the inversion guard is C++-core only; "
            "pass guard_inversion=False to use the numpy fallback"
        )
    taubin = method == "taubin"
    lam = lambda_
    if lam < 0.0:
        lam = 0.33 if taubin else 0.5
    if not 0.0 < lam < 1.0:
        raise ValueError(f"meshio++: smooth: lambda must lie in (0, 1); got {lam}")
    if taubin and not mu < -lam:
        raise ValueError(
            f"meshio++: smooth: taubin requires mu < -lambda < 0; got mu={mu}, lambda={lam}"
        )

    xyz, dim = _coords3(mesh)
    n = xyz.shape[0]
    xadj, adj = _edge_graph(mesh, n)

    pinned = _pinned_unknown_topology(mesh, n)
    if frozen is not None:
        ids = np.asarray(frozen, dtype=np.int64).ravel()
        if ids.size and (ids.min() < 0 or ids.max() >= n):
            raise ValueError("meshio++: smooth: frozen node id is out of range")
        pinned[ids] = True
    if fix_boundary:
        boundary, features = _boundary_and_features(
            mesh, xyz, n, _has_volume_cells(mesh), preserve_features, feature_angle
        )
        pinned |= boundary
        pinned |= features

    deg = np.diff(xadj)
    movable = (~pinned) & (deg > 0)
    passes = ([lam, mu] if taubin else [lam]) * max(int(iterations), 0)
    original = xyz.copy()
    cur = xyz
    for factor in passes:
        if adj.size == 0:
            break  # no edges at all: every node is isolated and holds still
        # reduceat needs in-range start indices; a zero-degree node's start may
        # be adj.size, so clamp and zero those rows afterwards (they are pinned
        # by `movable` anyway).
        starts = np.minimum(xadj[:-1], adj.size - 1)
        sums = np.add.reduceat(cur[adj], starts, axis=0)
        sums[deg == 0] = 0.0
        centroid = sums / np.maximum(deg, 1)[:, None]
        nxt = cur + factor * (centroid - cur)
        cur = np.where(movable[:, None], nxt, cur)

    disp = np.linalg.norm(cur - original, axis=1)
    extent = original.max(axis=0) - original.min(axis=0) if n else np.zeros(3)
    diag = float(np.linalg.norm(extent)) or 1.0
    out = Mesh(
        np.ascontiguousarray(cur[:, :dim]).astype(
            np.asarray(mesh.points).dtype, copy=False
        ),
        [(b.type, b.data) for b in mesh.cells],
    )
    out.point_data = {k: v for k, v in mesh.point_data.items()}
    out.cell_data = {k: list(v) for k, v in mesh.cell_data.items()}
    out.field_data = {k: v for k, v in mesh.field_data.items()}
    report = {
        "num_nodes_moved": int((disp > 1e-12 * diag).sum()),
        "max_displacement": float(disp.max()) if n else 0.0,
        "num_skipped_inversion": 0,
    }
    return out, report


def smooth(
    mesh,
    method="taubin",
    iterations=10,
    lambda_=-1.0,
    mu=-0.34,
    fix_boundary=True,
    preserve_features=True,
    feature_angle=30.0,
    guard_inversion=True,
    frozen=None,
    return_report=False,
):
    """Smooth a mesh's point coordinates, leaving topology and data intact.

    :param mesh: the mesh to smooth (never modified).
    :param method: ``"taubin"`` (default, shrink-free), ``"laplacian"``
        (stronger per pass, but shrinks the mesh), or ``"odt"``
        (optimal-Delaunay-triangulation smoothing: each free interior tet
        vertex moves toward the volume-weighted average of its incident
        tets' circumcenters -- **tet-only**, raises on any other block; see
        ``doc/smooth.md#odt-smoothing``. C++-core only, no numpy fallback,
        unlike the other two methods).
    :param iterations: how many iterations to run; for Taubin one iteration is
        two passes (``+lambda_`` then ``mu``).
    :param lambda_: relaxation factor in ``(0, 1)``. **Negative means "this
        method's own default"** -- ``0.5`` for Laplacian, ``0.33`` for
        Taubin, ``0.9`` for ODT (a near-full step, damped rather than
        exactly ``1.0``). The trailing underscore is because ``lambda`` is a
        Python keyword.
    :param mu: Taubin's un-shrinking factor; must satisfy ``mu < -lambda_ < 0``.
        Silently ignored by ``"laplacian"`` and ``"odt"``.
    :param fix_boundary: pin nodes lying on a boundary facet (default on --
        without it the domain's own shape changes).
    :param preserve_features: additionally pin boundary nodes where incident
        boundary facets meet at more than ``feature_angle``, so corners survive.
    :param feature_angle: that angle, in degrees, between facet normals.
    :param guard_inversion: never let a move turn a valid cell inverted.
    :param frozen: optional extra node ids to pin, as an index array or the name
        of one of ``mesh.point_sets``.
    :param return_report: also return the run summary dict.
    :returns: the smoothed mesh, or ``(mesh, report)`` if ``return_report``.
    """
    if isinstance(frozen, str):
        try:
            frozen = np.asarray(mesh.point_sets[frozen], dtype=np.int64)
        except KeyError:
            raise ValueError(
                f"meshio++: smooth: no point_set named {frozen!r} "
                f"(available: {sorted(mesh.point_sets)})"
            ) from None

    out = None
    report = None
    try:
        from . import _core

        res = _core.smooth(
            mesh,
            method,
            int(iterations),
            float(lambda_),
            float(mu),
            bool(fix_boundary),
            bool(preserve_features),
            float(feature_angle),
            bool(guard_inversion),
            None if frozen is None else np.asarray(frozen, dtype=np.int64),
        )
        out = res["mesh"]
        report = {
            "num_nodes_moved": res["num_nodes_moved"],
            "max_displacement": res["max_displacement"],
            "num_skipped_inversion": res["num_skipped_inversion"],
        }
    except (ValueError, TypeError):
        # A genuine user error (bad method name, lambda out of range, a frozen
        # id out of bounds) must not fall through to the numpy path and there
        # either succeed silently or raise something less precise.
        raise
    except Exception:
        out = None
    if out is None:
        out, report = _smooth_py(
            mesh,
            method,
            iterations,
            lambda_,
            mu,
            fix_boundary,
            preserve_features,
            feature_angle,
            guard_inversion,
            frozen,
        )

    # Sets are carried through untouched: smooth never renumbers a node or cell.
    if getattr(mesh, "point_sets", None):
        out.point_sets = {k: np.asarray(v).copy() for k, v in mesh.point_sets.items()}
    if getattr(mesh, "cell_sets", None):
        out.cell_sets = {
            k: [np.asarray(a).copy() for a in v] for k, v in mesh.cell_sets.items()
        }

    return (out, report) if return_report else out
