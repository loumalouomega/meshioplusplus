"""Geodesic distance on a triangle surface by the heat method (roadmap 7.3.1).

A dependency-free mesh *operation*: Crane, Weischedel and Wardetzky's
"Geodesics in Heat" (ACM TOG 32(5), 2013) over the surface's cotangent
stiffness and lumped mass -- a short heat flow from the sources, its normalised
gradient, and a Poisson solve for the distance. See ``doc/geodesic.md``.

**This operation is C++-core only, with no numpy fallback** (the
``sobolev_deform`` precedent): both solves stop on a rounded test, so a second
implementation could stop a sweep early or late and disagree beyond the last
ulp. This module raises when ``_core`` cannot be used.
"""

from __future__ import annotations

import numpy as np

__all__ = ["geodesic_distance", "DISTANCE_NAME"]

#: Point data: the surface distance to the nearest source. Twins
#: ``kGeodesicDistanceName``.
DISTANCE_NAME = "geodesic:distance"


def _source_ids(mesh, sources):
    """``sources`` as point ids and a region name (either may be empty)."""
    if sources is None:
        return [], ""
    if isinstance(sources, str):
        sets = getattr(mesh, "point_sets", None) or {}
        if sources in sets and not any(
            getattr(r, "name", None) == sources for r in getattr(mesh, "regions", [])
        ):
            return np.asarray(sets[sources], dtype=np.int64).ravel().tolist(), ""
        return [], sources
    if np.ndim(sources) == 0:
        return [int(sources)], ""
    return np.asarray(sources, dtype=np.int64).ravel().tolist(), ""


def geodesic_distance(
    mesh,
    sources=None,
    *,
    source_region=None,
    time_factor: float = 1.0,
    tolerance: float = 1e-10,
    max_iterations: int = 0,
    return_report: bool = False,
):
    """The distance along a triangle surface from ``sources`` to every point.

    Returns the mesh with ``point_data["geodesic:distance"]``. The heat method
    is first-order accurate and smooths the distance on the scale of
    ``sqrt(time_factor) * h`` (``h`` the mean edge length): on a unit sphere
    meshed at ``h = 0.03`` the pole-to-equator distance comes out within about
    1% of pi/2. A point in no triangle, and every point of a connected piece
    with no source, gets ``NaN``.

    :param mesh: a surface: its ``triangle`` blocks (never modified). A volume
        cell, a quad, a polygon or a higher-order cell is refused naming the
        operation that makes triangles; line and vertex blocks ride along.
    :param sources: the source points: a point id, a sequence of them, or the
        name of a region (its points: a point region's entries, a cell
        region's cells' points, a side region's facets' points) or of a
        ``point_set``.
    :param source_region: a region whose points are sources too.
    :param time_factor: the heat time is ``time_factor * h**2`` (1, Crane et
        al.'s choice; larger is smoother and, on smooth surfaces, closer).
    :param tolerance: the relative tolerance of both solves.
    :param max_iterations: the iteration cap of each solve; 0 picks one from
        the point count.
    :param return_report: also return ``{"num_sources", "num_unreached",
        "num_degenerate", "max_distance", "time_step", "heat_iterations",
        "poisson_iterations", "converged"}``.
    :returns: the mesh with the distance, or ``(mesh, report)``.
    :raises ValueError: without a source, with a source outside the mesh or in
        no triangle, or with an out-of-scope cell type.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    ids, region = _source_ids(mesh, sources)
    if source_region:
        if region and region != source_region:
            raise ValueError(
                "meshio++: geodesic_distance: name one source region, not two "
                f"({region!r} and {source_region!r})"
            )
        region = str(source_region)
    try:
        from . import _core
    except ImportError:
        _core = None
    if _core is None or not hasattr(_core, "geodesic_distance"):
        raise NotImplementedError(
            "meshio++: geodesic_distance: this operation is C++-core only and has no "
            "pure-numpy fallback -- install a build with the compiled "
            "meshioplusplus._core extension"
        )
    res = _core.geodesic_distance(
        mesh,
        ids,
        region,
        float(time_factor),
        float(tolerance),
        int(max_iterations),
    )
    out = res.pop("mesh")
    if not res["converged"]:
        from ._common import warn

        warn(
            "geodesic_distance: a solve hit its iteration cap "
            f"(heat {res['heat_iterations']}, Poisson {res['poisson_iterations']}); "
            "the distance is approximate -- raise max_iterations"
        )
    return (out, res) if return_report else out
