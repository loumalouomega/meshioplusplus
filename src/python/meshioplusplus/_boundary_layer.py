"""Prismatic boundary layers grown from a surface (roadmap 6.9).

From a triangle and quad surface, ``num_layers`` layers of cells grow along
the surface's normals: wedges from triangles, hexahedra from quads, the first
``first_height`` thick and each next ``growth`` times the one before. Each
point grows along its *most normal* normal, the direction that sees all its
faces best (A. Aubry and R. Löhner, *On the "most normal" normal*, 2009), its
offsets stretched so every face's layers keep their thickness at ridges and
corners; a ray against the surface stops layers before they meet across a gap
(``collision_fraction`` of the distance), the fraction smoothed over
neighbours; and the points of a cell that would fold have their layers halved
and at last collapsed, so their cells become pyramids and tetrahedra and no
cell is written folded. The design is ``doc/boundary_layer_plan.md``.

The C++ core does the work; there is no numpy twin.
"""

from __future__ import annotations

import math

__all__ = ["boundary_layer"]

_PREFIX = "meshio++: boundary_layer: "

_REPORT_INTS = (
    "num_cells",
    "num_wedges",
    "num_hexahedra",
    "num_pyramids",
    "num_tetrahedra",
    "num_shrunk",
    "num_collapsed",
    "num_not_visible",
)
_REPORT_FLOATS = ("min_scaled_jacobian", "min_thickness", "max_thickness")


def boundary_layer(
    mesh,
    first_height: float,
    *,
    num_layers: int = 3,
    growth: float = 1.2,
    reverse: bool = False,
    collision_fraction: float = 0.45,
    smooth_iterations: int = 5,
    max_shrink: int = 6,
    return_report: bool = False,
):
    """Grow boundary layers from a triangle and quad surface.

    :param mesh: the surface: ``triangle`` and ``quad`` blocks (lines and
        vertices are skipped), manifold and consistently oriented; an open
        surface is fine.
    :param first_height: the first layer's thickness (positive).
    :param num_layers: how many layers (1 to 1000).
    :param growth: each layer's thickness over the one before (1 to 10).
    :param reverse: grow against the faces' normals (right-hand rule) rather
        than along them.
    :param collision_fraction: the share of the distance to the surface a
        point's ray hits that its layers may take, in ``(0, 1]``.
    :param smooth_iterations: passes smoothing that share over neighbours
        (never raising it).
    :param max_shrink: halvings of a point's layers, where a cell would fold,
        before they collapse.
    :param return_report: also return the report.
    :returns: the mesh (``wedge``, ``hexahedron``, ``pyramid`` and ``tetra``
        cells, then the wall and cap faces, with ``boundary_layer:layer`` cell
        data, the surface's cell data, point data and regions carried, and the
        ``boundary_layer:cells`` / ``:wall`` / ``:cap`` regions), or
        ``(mesh, report)``: ``num_cells``, ``num_wedges``,
        ``num_hexahedra``, ``num_pyramids``, ``num_tetrahedra``,
        ``num_shrunk``, ``num_collapsed``, ``num_not_visible``,
        ``min_scaled_jacobian``, ``min_thickness``, ``max_thickness``,
        ``cap_closed`` and ``warnings``.
    :raises ValueError: for an option out of range, or a surface with no
        triangle or quad, another cell type, a non-manifold edge or faces that
        disagree on an edge's direction.
    """
    first_height = float(first_height)
    if not (math.isfinite(first_height) and first_height > 0.0):
        raise ValueError(f"{_PREFIX}first_height must be positive and finite")
    if isinstance(num_layers, bool) or int(num_layers) != num_layers:
        raise ValueError(f"{_PREFIX}num_layers must be an integer")
    if not 1 <= int(num_layers) <= 1000:
        raise ValueError(f"{_PREFIX}num_layers must be between 1 and 1000")
    growth = float(growth)
    if not (math.isfinite(growth) and 1.0 <= growth <= 10.0):
        raise ValueError(f"{_PREFIX}growth must be between 1 and 10")
    collision_fraction = float(collision_fraction)
    if not 0.0 < collision_fraction <= 1.0:
        raise ValueError(f"{_PREFIX}collision_fraction must be in (0, 1]")
    try:
        from . import _core
    except ImportError as exc:  # pragma: no cover - the numpy-only install
        raise NotImplementedError(
            "meshio++: boundary_layer: needs the compiled C++ core "
            "(meshioplusplus._core); it has no numpy twin"
        ) from exc
    res = _core.boundary_layer(
        mesh,
        first_height,
        growth=growth,
        num_layers=int(num_layers),
        reverse=bool(reverse),
        collision_fraction=collision_fraction,
        smooth_iterations=int(smooth_iterations),
        max_shrink=int(max_shrink),
    )
    out = res["mesh"]
    if not return_report:
        return out
    report = {k: int(res[k]) for k in _REPORT_INTS}
    report.update({k: float(res[k]) for k in _REPORT_FLOATS})
    report["cap_closed"] = bool(res["cap_closed"])
    report["warnings"] = [str(w) for w in res["warnings"]]
    return out, report
