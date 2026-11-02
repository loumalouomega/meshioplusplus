"""The convex hull and the planar Delaunay triangulation of a point set
(roadmap 7.5.1).

``convex_hull(points)`` is the closed hull of scattered points (quickhull in
3-D, Andrew's monotone chain in 2-D) and ``delaunay_2d(points)`` the Delaunay
triangulation of their ``x, y`` coordinates: an ``xyz``, ``pcd`` or ``ply``
scan turned into a hull, a footprint or a terrain. See ``doc/hull.md``.

**These operations are C++-core only, with no numpy fallback** (the
``geodesic_distance`` precedent): every decision is an exact geometric
predicate, which floating-point numpy cannot reproduce on degenerate input
(cocircular, coplanar or collinear points), so a twin would disagree exactly
where it matters. This module raises when ``_core`` cannot be used.
"""

from __future__ import annotations

import numpy as np

__all__ = ["convex_hull", "delaunay_2d", "HULL_POINT_ID_NAME"]

#: Point data on a hull: the input point each hull vertex is. Twins
#: ``kHullPointIdName``.
HULL_POINT_ID_NAME = "hull:point_id"


def _as_mesh(points_or_mesh):
    from ._mesh import Mesh

    if hasattr(points_or_mesh, "points") and hasattr(points_or_mesh, "cells"):
        return points_or_mesh
    pts = np.asarray(points_or_mesh, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] not in (2, 3):
        raise ValueError(
            "meshio++: expected a mesh or an (n, 2) or (n, 3) array of points, got "
            f"shape {pts.shape}"
        )
    return Mesh(pts, [])


def _core_op(name):
    try:
        from . import _core
    except ImportError:
        _core = None
    if _core is None or not hasattr(_core, name):
        raise NotImplementedError(
            f"meshio++: {name}: this operation is C++-core only and has no pure-numpy "
            "fallback -- install a build with the compiled meshioplusplus._core "
            "extension"
        )
    return getattr(_core, name)


def convex_hull(points, *, dimension: int = 0, return_report: bool = False):
    """The convex hull of a point set.

    In 3-D a closed ``triangle`` surface wound outward; in 2-D (the ``x, y``
    coordinates) a counter-clockwise closed ``line`` loop. Only the extreme
    points are vertices; the output holds them alone, in ascending input
    order, with their point data and ``hull:point_id`` (the input point each
    one is). Every decision is exact, so coplanar and cospherical points are
    handled without a tolerance.

    :param points: a mesh (its points; its cells are ignored) or an
        ``(n, 2)``/``(n, 3)`` array.
    :param dimension: 2, 3, or 0 to take 2 for 2-D points and for 3-D points
        whose ``z`` are all equal, else 3.
    :param return_report: also return ``{"dimension", "num_vertices",
        "num_facets", "area", "volume", "perimeter"}`` (``area`` is the
        enclosed area in 2-D and the surface area in 3-D).
    :returns: the hull mesh, or ``(mesh, report)``.
    :raises ValueError: on a non-finite coordinate, fewer than three points,
        or points that span no hull (collinear in 2-D, coplanar in 3-D).
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    res = _core_op("convex_hull")(_as_mesh(points), int(dimension))
    out = res.pop("mesh")
    return (out, res) if return_report else out


def delaunay_2d(points, *, return_report: bool = False):
    """The Delaunay triangulation of a point set's ``x, y`` coordinates.

    The output keeps every input point -- its coordinates (a terrain keeps
    its ``z``) and its point data -- with one counter-clockwise ``triangle``
    block covering the convex hull. Every decision is exact; among cocircular
    points the choice is deterministic. A point that coincides in ``x, y``
    with an earlier one is in no triangle.

    :param points: a mesh (its points; its cells are replaced) or an
        ``(n, 2)``/``(n, 3)`` array.
    :param return_report: also return ``{"num_triangles", "num_duplicates",
        "num_boundary"}``.
    :returns: the triangulated mesh, or ``(mesh, report)``.
    :raises ValueError: on a non-finite coordinate or one outside the exact
        range (a magnitude above 2**128 or below 2**-128), fewer than three
        distinct points, or collinear points.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    res = _core_op("delaunay_2d")(_as_mesh(points))
    out = res.pop("mesh")
    return (out, res) if return_report else out
