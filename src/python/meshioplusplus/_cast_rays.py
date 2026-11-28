"""Rays against a surface and point location in a mesh (roadmap 9.3).

``cast_rays(mesh, origins, directions)`` reports the nearest surface hit of
each ray; ``locate_points(mesh, points)`` reports the cell that contains each
point and the barycentric weights inside it. Both return dictionaries of NumPy
arrays. See ``doc/ray_queries.md``.

**C++-core only, with no numpy fallback**: the search structures are the
core's, and ``locate_points`` shares ``interpolate``'s simplex search. These
functions raise when ``_core`` cannot be used.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

import numpy as np

__all__ = ["cast_rays", "locate_points", "SurfaceQuery"]


def cast_rays(
    mesh,
    origins,
    directions,
    *,
    max_distance: float = 0.0,
    region: Optional[str] = None,
    all_hits: bool = False,
) -> Dict[str, Any]:
    """Surface hits of every ray (nearest by default, or all crossings with ``all_hits=True``).

    :param mesh: a mesh with triangle, quad or polygon cells (a volume mesh is
        refused: take its ``extract_surface`` first).
    :param origins: ``(n, 3)`` ray origins.
    :param directions: one direction ``(3,)`` for every ray, or ``(n, 3)``
        (not necessarily unit vectors).
    :param max_distance: ignore hits farther than this; ``0`` for no limit.
    :param region: restrict the surface to this named cell region.
    :param all_hits: when ``True``, return every crossing along each ray in CSR
        layout (adding ``ray_offset`` of length ``n + 1``).
    :returns: a dict of ``hit`` (``(n,)`` bool), ``distance`` (``(n,)`` or
        ``(m,)``), ``point`` and ``normal`` (``(n, 3)`` or ``(m, 3)``),
        ``cell`` (``(n,)`` or ``(m,)`` int64), ``num_hits``, and (when
        ``all_hits=True``) ``ray_offset`` (``(n + 1,)`` int64).
    :raises ValueError: on a zero or non-finite direction, a volume or
        higher-order block, an unknown region, or a mesh with no surface cell.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _as_mesh, _core_op

    return _core_op("cast_rays")(
        _as_mesh(mesh),
        np.ascontiguousarray(origins, dtype=np.float64).reshape(-1, 3),
        np.ascontiguousarray(directions, dtype=np.float64),
        float(max_distance),
        "" if region is None else str(region),
        bool(all_hits),
    )


def locate_points(mesh, points, *, tolerance: float = 1e-12) -> Dict[str, Any]:
    """The cell containing each point, and its barycentric weights.

    The mesh is reduced to simplices (hexahedra, prisms, pyramids, polyhedra and
    quadrilaterals are split; higher-order cells are linearised) and searched as
    ``interpolate`` searches it: the first containing simplex in ascending order
    wins, so a point on a shared face belongs to the lowest-numbered cell.

    :param mesh: a volume mesh, or a planar triangle mesh in the xy-plane.
    :param points: ``(n, 3)`` points.
    :param tolerance: a barycentric weight at least ``-tolerance`` counts as
        inside.
    :returns: a dict of ``cell`` (``(n,)`` int64, the global index of the input
        cell, -1 outside), ``nodes`` (``(n, 4)`` int64, the point ids of the
        containing simplex, -1 padding), ``barycentric`` (``(n, 4)`` weights with
        respect to those nodes, 0 padding) and ``num_located``.
    :raises ValueError: on a negative or non-finite tolerance, a non-finite
        coordinate, or a mesh with no tetrahedron or triangle after
        simplexification.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _as_mesh, _core_op

    return _core_op("locate_points")(
        _as_mesh(mesh),
        np.ascontiguousarray(points, dtype=np.float64).reshape(-1, 3),
        float(tolerance),
    )


class SurfaceQuery:
    """A surface's triangle soup and hierarchies, built once and queried many
    times (roadmap 9.3.2).

    ``SurfaceQuery(mesh, region=None).cast_rays(origins, directions)`` returns
    exactly what ``cast_rays(mesh, origins, directions, region=region)`` does,
    without rebuilding the hierarchy per batch; ``closest_points(points)``
    returns the nearest surface point of each query point. The handle copies
    the surface: later changes to ``mesh`` do not reach it.

    :param mesh: a mesh with triangle, quad or polygon cells.
    :param region: restrict the surface to this named cell region.
    :raises ValueError: on an unknown region or a mesh with no surface cell.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """

    def __init__(self, mesh, region: Optional[str] = None):
        from ._hull import _as_mesh, _core_op

        cls = _core_op("SurfaceQuery")
        self._query = cls(_as_mesh(mesh), "" if region is None else str(region))

    @property
    def num_triangles(self) -> int:
        """The surface's triangles (quads and polygons fanned)."""
        return int(self._query.num_triangles)

    def cast_rays(
        self, origins, directions, *, max_distance: float = 0.0, all_hits: bool = False
    ) -> Dict[str, Any]:
        """:func:`cast_rays` against this surface; same arguments and result."""
        return self._query.cast_rays(
            np.ascontiguousarray(origins, dtype=np.float64).reshape(-1, 3),
            np.ascontiguousarray(directions, dtype=np.float64),
            float(max_distance),
            bool(all_hits),
        )

    def closest_points(self, points) -> Dict[str, Any]:
        """The nearest surface point of each of ``points`` (``(n, 3)``): a dict
        of ``point`` (``(n, 3)``), ``distance`` (``(n,)``, unsigned) and
        ``cell`` (``(n,)`` int64, the global input cell)."""
        return self._query.closest_points(
            np.ascontiguousarray(points, dtype=np.float64).reshape(-1, 3)
        )


def _xyz(points) -> "np.ndarray":
    """A mesh's points as ``(n, 3)`` float64 (a planar mesh gets ``z = 0``)."""
    xyz = np.asarray(points, dtype=np.float64)
    if xyz.ndim == 2 and xyz.shape[1] < 3:
        xyz = np.pad(xyz, ((0, 0), (0, 3 - xyz.shape[1])))
    return xyz


def _vertex_mesh(points, point_data):
    from ._mesh import Mesh

    n = len(points)
    return Mesh(
        points,
        [("vertex", np.arange(n, dtype=np.int64).reshape(-1, 1))],
        point_data=point_data,
    )


def _hits_mesh(res):
    """The vertex mesh of the rays that hit (the CLI and MCP output)."""
    if "ray_offset" in res:
        offsets = res["ray_offset"]
        counts = np.diff(offsets)
        ray_index = np.repeat(np.arange(len(counts), dtype=np.int64), counts)
        return _vertex_mesh(
            res["point"],
            {
                "ray:index": ray_index,
                "ray:distance": res["distance"],
                "ray:cell": res["cell"],
                "ray:normal": res["normal"],
            },
        )
    hit = res["hit"]
    return _vertex_mesh(
        res["point"][hit],
        {
            "ray:index": np.flatnonzero(hit).astype(np.int64),
            "ray:distance": res["distance"][hit],
            "ray:cell": res["cell"][hit],
            "ray:normal": res["normal"][hit],
        },
    )


def _located_mesh(xyz, res):
    """The vertex mesh of the points that fell inside a cell."""
    inside = res["cell"] >= 0
    return _vertex_mesh(
        xyz[inside],
        {
            "locate:cell": res["cell"][inside],
            "locate:nodes": res["nodes"][inside],
            "locate:barycentric": res["barycentric"][inside],
        },
    )


def _directions_of(rays, direction):
    """The direction argument of ``cast_rays`` for a rays mesh: ``direction``
    (three numbers) when given, else the point array ``directions``."""
    if direction is not None:
        return [float(v) for v in direction]
    if "directions" in rays.point_data:
        return np.asarray(rays.point_data["directions"], dtype=np.float64)
    raise ValueError(
        "cast_rays: give a direction, or a point array 'directions' in the rays file"
    )
