"""Move a mesh's points back onto the CAD model it was made from (roadmap
6.10.1).

A mesh :func:`meshioplusplus.mesh_cad_surface` wrote names, cell by cell, the
CAD face it lies on (``cad:face`` on its triangles and quads) and the CAD edge
it runs along (``cad:edge`` on its lines). An operation that knows no CAD
(smoothing, an adaptation, a remesh that keeps those arrays) moves points off
the surfaces; :func:`project_to_cad` puts them back: a point on one edge's
lines onto that edge, a point on two or more edges (a corner) onto the lowest
of them, a point whose cells all name one face onto that face's surface (the
nearest point, by Newton on the first-order condition with the surface's
second derivatives, from the nearest few of a grid of the patch's points; L.
Ma and W. T. Hewitt, *Point inversion and projection for NURBS curve and
surface*, 2003). A point on two faces and no edge is ambiguous and stays. A
projection outside its face's trimmed region is kept and reported, never
clamped. Data, regions and property sets pass through.

The C++ core does the work; there is no numpy twin.
"""

from __future__ import annotations

import math

from ._mesh_cad import _brep_of, _core_or_raise

__all__ = ["project_to_cad"]

_PREFIX = "meshio++: project_to_cad: "


def project_to_cad(
    mesh,
    model,
    *,
    face_array: str = "cad:face",
    edge_array: str = "cad:edge",
    max_distance: float = 0.0,
    trim_samples: int = 16,
    return_report: bool = False,
):
    """Move every point of ``mesh`` onto the face or edge of ``model`` its
    cells name.

    :param mesh: a surface mesh whose triangles and quads carry
        ``face_array`` (each cell's CAD face index, -1 for none) and whose
        lines may carry ``edge_array`` (each line's CAD edge index), as
        :func:`meshioplusplus.mesh_cad_surface` writes them.
    :param model: the CAD model: a mesh read from a CAD JSON, STEP or IGES
        file (carrying ``info["brep"]``), or the B-rep dict itself.
    :param face_array: the cell data naming each cell's face.
    :param edge_array: the cell data naming each line's edge.
    :param max_distance: a point farther than this from where it would go
        stays and is counted (0, the default, is no limit).
    :param trim_samples: samples per knot span when testing whether a
        projection lies in its face's trimmed region.
    :param return_report: also return the report.
    :returns: the mesh with its points moved, or ``(mesh, report)``. The
        report has ``num_projected``, ``num_on_edges``, ``num_on_curves``
        (of those on edges, onto the edge's own 3-D curve), ``num_corners``,
        ``num_ambiguous``, ``num_too_far``, ``num_failed``,
        ``num_outside_trim``, ``max_move``, ``outside_trim`` (the points
        projected outside their trims) and ``warnings``.
    :raises ValueError: when the mesh has no ``face_array`` cell data, its
        points are not three-dimensional, an index is out of range, or an
        option is invalid.
    """
    brep = _brep_of(model)
    max_distance = float(max_distance)
    if not (math.isfinite(max_distance) and max_distance >= 0.0):
        raise ValueError(f"{_PREFIX}max_distance must be finite and not negative")
    if isinstance(trim_samples, bool) or int(trim_samples) != trim_samples:
        raise ValueError(f"{_PREFIX}trim_samples must be an integer")
    if not 2 <= int(trim_samples) <= 1024:
        raise ValueError(f"{_PREFIX}trim_samples must be between 2 and 1024")
    _core = _core_or_raise("project_to_cad")
    res = _core.project_to_cad(
        mesh,
        brep,
        face_array=str(face_array),
        edge_array=str(edge_array),
        max_distance=max_distance,
        trim_samples=int(trim_samples),
    )
    out = res["mesh"]
    info = getattr(mesh, "info", None)
    if isinstance(info, dict) and info:
        out.info = dict(info)
    if not return_report:
        return out
    report = {
        k: int(res[k])
        for k in (
            "num_projected",
            "num_on_edges",
            "num_on_curves",
            "num_corners",
            "num_ambiguous",
            "num_too_far",
            "num_failed",
            "num_outside_trim",
        )
    }
    report["max_move"] = float(res["max_move"])
    report["outside_trim"] = [int(i) for i in res["outside_trim"]]
    report["warnings"] = [str(w) for w in res["warnings"]]
    return out, report
