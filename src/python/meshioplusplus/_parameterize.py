"""A UV atlas of a triangle surface (roadmap 7.3.3).

``parameterize`` cuts a surface into charts, flattens each without flipping a
triangle (least-squares conformal maps, else mean-value coordinates, the
kernels ``fit_cad`` uses), and shelf-packs them into the unit square. The
result carries ``point_data["uv"]``, which the OBJ, glTF and USD writers turn
into texture coordinates. See ``doc/parameterize.md``.

**This operation is C++-core only, with no numpy fallback** (the
``geodesic_distance`` precedent): the conformal solve is iterative, so a second
implementation would not agree to the last bit. This module raises when
``_core`` cannot be used.
"""

from __future__ import annotations

__all__ = ["parameterize", "UV_NAME", "PARENT_POINT_NAME", "CHART_NAME"]

#: Point data: the texture coordinates, Float64 ``(n, 2)`` in ``[0, 1]``.
#: Twins ``kUvName``.
UV_NAME = "uv"
#: Point data: the input point each output point copies. Twins
#: ``kParameterizeParentPointName``.
PARENT_POINT_NAME = "parameterize:parent_point"
#: Cell data: each triangle's chart (-1 on other blocks). Twins
#: ``kParameterizeChartName``.
CHART_NAME = "parameterize:chart"


def parameterize(
    mesh,
    *,
    padding: float = 0.02,
    max_depth: int = 8,
    max_stretch: float = 2.0,
    return_report: bool = False,
):
    """A UV atlas of a triangle surface.

    Returns the surface with ``point_data["uv"]`` (two coordinates per point in
    the unit square). Each edge-connected piece that is a disc is one chart; a
    closed or holed piece is cut into discs, and one that still does not
    flatten without a flipped triangle, or whose map stretches more than
    ``max_stretch``, is halved along the direction its normals spread most, up
    to ``max_depth`` times. Every chart keeps its area
    (one texel density), its triangles' winding, and is turned to its
    principal axes before packing.

    A point on a seam belongs to two charts and is split: the first copy keeps
    its index, the others are appended, and
    ``point_data["parameterize:parent_point"]`` names the input point of each,
    so point regions, lines and vertex cells stay valid and ``uv`` holds one
    value per point. Point data is copied to the copies, cell data and regions
    ride along, and ``cell_data["parameterize:chart"]`` gives each triangle's
    chart.

    :param mesh: a surface: its ``triangle`` blocks (never modified). A volume
        cell, a quad, a polygon or a higher-order cell is refused naming the
        operation that makes triangles; line and vertex blocks ride along.
    :param padding: the gap between charts, as a fraction of the square root of
        the surface's area.
    :param max_depth: how many times a piece that does not flatten may be
        halved; one still folding at the cap keeps its map and counts its
        flipped triangles.
    :param max_stretch: a chart whose L2 stretch exceeds this (1 keeps
        lengths) is halved again while ``max_depth`` allows; ``float("inf")``
        keeps every chart that does not fold.
    :param return_report: also return ``{"num_charts", "num_flipped",
        "stretch", "packing_efficiency", "num_seam_points"}``: ``stretch`` is
        Sander et al.'s area-weighted L2 stretch (1 for a map that keeps
        lengths), ``packing_efficiency`` the fraction of the unit square the
        triangles cover.
    :returns: the mesh with ``uv``, or ``(mesh, report)``.
    :raises ValueError: with an out-of-scope cell type, no triangle, or a
        negative ``padding`` or ``max_depth``, or a ``max_stretch`` below 1.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    try:
        from . import _core
    except ImportError:
        _core = None
    if _core is None or not hasattr(_core, "parameterize"):
        raise NotImplementedError(
            "meshio++: parameterize: this operation is C++-core only and has no "
            "pure-numpy fallback -- install a build with the compiled "
            "meshioplusplus._core extension"
        )
    res = _core.parameterize(mesh, float(padding), int(max_depth), float(max_stretch))
    out = res.pop("mesh")
    return (out, res) if return_report else out
