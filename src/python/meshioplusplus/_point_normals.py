"""Normals of a point cloud by local plane fits (roadmap 7.5.2, first part).

``estimate_normals(points, k=16)`` fits a plane to each point and its ``k``
nearest neighbours and orients the normals consistently (Hoppe et al.'s
minimum spanning tree propagation, then outward), writing
``point_data["normals"]`` -- the array the PCD and XYZ writers emit as normal
columns. See ``doc/point_normals.md``.

**C++-core only, with no numpy fallback** (the ``geodesic_distance``
precedent): the propagation walks a minimum spanning tree whose ties a second
implementation would have to break bit for bit the same way, or flip whole
regions differently. This module raises when ``_core`` cannot be used.
"""

from __future__ import annotations

__all__ = ["estimate_normals", "VARIATION_NAME"]

#: Point data: the neighbourhood's surface variation. Twins
#: ``kNormalsVariationName``.
VARIATION_NAME = "normals:variation"


def estimate_normals(
    points,
    *,
    k: int = 16,
    orientation: str = "consistent",
    viewpoint=(0.0, 0.0, 0.0),
    record_variation: bool = False,
    return_report: bool = False,
):
    """Normals of a point cloud: a plane fitted to each point's neighbourhood.

    Returns the mesh (or a new one around an array of points) with
    ``point_data["normals"]``, unit Float64 ``(n, 3)``. Points, cells and
    the other arrays pass through.

    :param points: a mesh (its points; cells ride along) or an ``(n, 2)`` /
        ``(n, 3)`` array.
    :param k: neighbours per fit (the point itself is added); at least 2.
    :param orientation: ``"consistent"`` (Hoppe et al.'s propagation over the
        neighbour graph, then each connected piece turned to point away from
        its centroid -- outward for a closed scan; a flat patch faces +z),
        ``"viewpoint"`` (towards ``viewpoint``, the scanner) or ``"none"``
        (each normal's largest component positive).
    :param record_variation: also attach ``normals:variation``, the surface
        variation ``l0 / (l0 + l1 + l2)`` of the neighbourhood (0 on a plane).
    :param return_report: also return ``{"num_degenerate",
        "num_components"}``.
    :returns: the mesh, or ``(mesh, report)``. A point whose neighbourhood
        spans no plane gets a ``NaN`` normal (``num_degenerate``).
    :raises ValueError: on fewer than three points, ``k`` below 2, an unknown
        orientation, or a non-finite coordinate.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _as_mesh, _core_op

    res = _core_op("estimate_normals")(
        _as_mesh(points),
        int(k),
        str(orientation),
        [float(v) for v in viewpoint],
        bool(record_variation),
    )
    out = res.pop("mesh")
    return (out, res) if return_report else out
