"""Smooth subdivision surfaces: Loop (triangles) and Catmull-Clark (polygonal -> quads).

Thin Python shim over the C++ core (`operations/subdivide_surface.hpp`, roadmap 9.2).
"""

from __future__ import annotations

from ._fallback import core_op_declined
from ._mesh import Mesh


def subdivide_surface(
    mesh: Mesh,
    scheme: str = "loop",
    *,
    levels: int = 1,
    crease_angle: float = 60.0,
    preserve_boundary: bool = True,
    interpolate_point_data: bool = True,
    record_parent_ids: bool = False,
    return_report: bool = False,
):
    """Apply smooth Loop or Catmull-Clark subdivision to a surface mesh.

    Parameters
    ----------
    mesh : Mesh
        A 2-D surface mesh (triangles for ``"loop"``; triangles, quads, or
        polygons for ``"catmull_clark"``).
    scheme : {"loop", "catmull_clark"}, default "loop"
        Subdivision scheme.
    levels : int, default 1
        Number of subdivision passes (>= 0).
    crease_angle : float, default 60.0
        Dihedral angle threshold in degrees for automatic sharp creases
        (<= 0 or >= 180 disables sharp-edge detection).
    preserve_boundary : bool, default True
        Apply Hoppe crease/boundary rules to open boundary edges (corners pinned).
    interpolate_point_data : bool, default True
        Interpolate floating-point ``point_data`` using the subdivision stencil.
    record_parent_ids : bool, default False
        Attach Int64 ``subdivide:parent_cell`` cell data.
    return_report : bool, default False
        When ``True``, return ``(mesh, report)`` with ``num_crease_edges``,
        ``num_boundary_edges``, and ``warnings``.
    """
    try:
        from . import _core
    except ImportError as exc:
        raise RuntimeError(
            "subdivide_surface requires the compiled C++ core (_core)"
        ) from exc

    try:
        out, report = _core.subdivide_surface(
            mesh,
            scheme,
            int(levels),
            float(crease_angle),
            bool(preserve_boundary),
            bool(interpolate_point_data),
            bool(record_parent_ids),
        )
    except Exception as exc:
        if not core_op_declined(exc, "subdivide_surface"):
            raise
        raise RuntimeError(
            f"subdivide_surface: C++ core declined input: {exc}"
        ) from exc

    return (out, report) if return_report else out
