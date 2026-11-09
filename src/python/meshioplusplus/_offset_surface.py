"""Offset a surface, thicken it into a shell, hollow a solid (roadmap 7.8.4).
C++-core only (see ``doc/offset_surface.md``)."""

from __future__ import annotations

__all__ = ["offset_surface", "thicken", "hollow"]

_KEYS = ("num_not_visible", "num_clamped", "num_folded", "num_rim_quads")


def _run(mesh, kind, amount, side, smooth_iterations, return_report):
    from ._hull import _core_op

    res = _core_op("offset_surface")(
        mesh, kind, float(amount), str(side), int(smooth_iterations)
    )
    out = res.pop("mesh")
    if not return_report:
        return out
    report = {k: int(res[k]) for k in _KEYS}
    report["max_move"] = float(res["max_move"])
    report["warnings"] = [str(w) for w in res["warnings"]]
    return out, report


def offset_surface(
    mesh, distance: float, *, smooth_iterations: int = 0, return_report=False
):
    """Move each point of a triangle and quad surface along its offset
    direction (the one ``boundary_layer`` grows along) by ``distance`` over
    the least cosine to its faces, so every face's plane moves by
    ``distance``: a cube offset by ``d`` is the cube grown by ``d``.

    :param distance: along the faces' normals (their winding); negative
        against them.
    :param smooth_iterations: smoothing passes over the directions.
    :returns: the moved mesh (connectivity, data and regions kept); with
        ``return_report`` ``(mesh, report)``: ``num_not_visible``,
        ``num_clamped``, ``num_folded``, ``num_rim_quads``, ``max_move`` and
        ``warnings``.
    :raises ValueError: on a non-finite distance or a cell that is not a
        triangle or a quad.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    return _run(mesh, "offset", distance, "inward", smooth_iterations, return_report)


def thicken(
    mesh,
    thickness: float,
    *,
    side: str = "inward",
    smooth_iterations: int = 0,
    return_report=False,
):
    """Turn a surface into a closed shell ``thickness`` thick: the surface
    offset to the two sides of the wall (``"inward"``: it and its offset by
    ``-thickness``; ``"outward"``: by ``+thickness``; ``"both"``: by
    ``+-thickness / 2``), the lower layer reversed and the boundary edges
    joined across by quads; the Int64 cell data ``thicken:side`` is 0 on the
    upper layer, 1 on the lower and 2 on the rim. See :func:`offset_surface`
    for the report.
    """
    return _run(mesh, "thicken", thickness, side, smooth_iterations, return_report)


def hollow(mesh, wall: float, *, smooth_iterations: int = 0, return_report=False):
    """Give a closed surface an inner wall: the surface and its offset by
    ``-wall``, reversed (an open surface is refused). See
    :func:`offset_surface` for the report.
    """
    return _run(mesh, "hollow", wall, "inward", smooth_iterations, return_report)
