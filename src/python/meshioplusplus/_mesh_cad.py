"""Mesh generation from CAD: conforming surface triangulation of a B-rep.

The surface mesher of ``operations/mesh_cad.hpp``. ``mesh_cad_surface``
triangulates the faces of a B-rep (``mesh.info["brep"]`` of a Kratos CAD
JSON, STEP or IGES read, or the bare dict) at one global target edge size,
with no CAD kernel: every edge is sampled once and shared by the faces that
use it, so the result conforms; each face is a constrained Delaunay
triangulation of its parameter domain, lifted onto its surface. A face that
cannot be meshed is reported by index and reason, never meshed wrong. With
``min_angle`` every face is refined by Delaunay refinement (Ruppert, Chew)
until its triangles' smallest angles on the surface reach the bound; corners
under 60 degrees and poles that keep triangles below it are reported.
``smooth_iterations`` adds guarded Laplacian or ODT smoothing.
``triangulate_polygon_2d`` exposes the planar constrained Delaunay
triangulation underneath, with the same ``min_angle``. There is no numpy twin: without the compiled core
both raise ``NotImplementedError``. See ``doc/mesh_cad.md``.
"""

from __future__ import annotations

import math

import numpy as np

_PREFIX = "meshio++: mesh_cad_surface: "

#: Why a face was left out, as the report's ``status`` names it.
FACE_STATUSES = (
    "ok",
    "empty",
    "open_loop",
    "self_intersecting",
    "degenerate",
    "projection_failed",
    "budget_exceeded",
    "unsupported",
    "other",
)

#: How smoothing moves an interior point.
SMOOTH_METHODS = ("laplacian", "odt")

_REPORT_INTS = (
    "num_points",
    "num_triangles",
    "num_edge_segments",
    "num_faces",
    "num_failed",
    "num_edges",
    "num_open_edges",
    "num_non_manifold_edges",
    "num_misoriented_edges",
    "num_degenerate_dropped",
    "num_refinement_points",
    "num_faces_below_bound",
    "num_smoothed",
    "num_curved_straightened",
)
_REPORT_FLOATS = (
    "tolerance",
    "min_angle_deg",
    "max_angle_deg",
    "min_edge",
    "max_edge",
    "max_edge_gap",
)


def _brep_of(mesh):
    """The B-rep dict a mesh carries (``mesh.info["brep"]``), or a bare one."""
    if isinstance(mesh, dict):
        if "faces" not in mesh:
            raise ValueError(f"{_PREFIX}the dict is not a B-rep (it has no 'faces')")
        return mesh
    info = getattr(mesh, "info", None) or {}
    brep = info.get("brep") if isinstance(info, dict) else None
    if brep is None:
        raise ValueError(
            f"{_PREFIX}the mesh carries no B-rep (mesh.info['brep']): read a CAD JSON, "
            "STEP or IGES file"
        )
    return brep


def _core_or_raise(what):
    try:
        from . import _core
    except ImportError as exc:
        raise NotImplementedError(
            f"meshio++: {what}: needs the compiled C++ core (meshioplusplus._core); "
            "it has no numpy twin"
        ) from exc
    return _core


def mesh_cad_surface(
    mesh,
    target_size: float,
    *,
    max_chordal_deviation: float = 0.0,
    min_points_per_edge: int = 0,
    tolerance: float = 0.0,
    record_face_ids: bool = True,
    edges: bool = True,
    match_unreferenced_trims: bool = True,
    strict: bool = False,
    max_points: int = 50_000_000,
    min_angle: float = 0.0,
    smooth_iterations: int = 0,
    smooth_method: str = "laplacian",
    return_report: bool = False,
):
    """Triangulate a B-rep's faces into one conforming surface mesh.

    :param mesh: a mesh carrying ``info["brep"]`` (a CAD JSON, STEP or IGES
        read), or the B-rep dict itself.
    :param target_size: the target edge length ``h`` (required, > 0).
    :param max_chordal_deviation: the largest distance a boundary segment may
        sit from its edge (its sagitta; OpenCASCADE's linear deflection);
        0 is off.
    :param min_points_per_edge: at least this many points inside every edge.
    :param tolerance: edges shorter than this are degenerate (a pole) and
        unreferenced trims closer than this are paired; 0 is ``1e-7`` times
        the diagonal of the control nets' box.
    :param record_face_ids: ``cad:face`` cell data and a ``face:<k>`` cell
        region per face.
    :param edges: a ``line`` block of the discretised edges, with ``cad:edge``.
    :param match_unreferenced_trims: pair trims no edge references whose 3-D
        images coincide, for models without edge topology.
    :param strict: raise ``RuntimeError`` when any face fails, instead of
        reporting it.
    :param max_points: the most points in the result.
    :param min_angle: refine until every triangle's smallest angle on the
        surface is at least this (degrees, under 60); 0 keeps the unrefined
        mesh. Up to 20.7 degrees it is guaranteed away from corners under 60
        degrees; above, a target whose shortfall is reported.
    :param smooth_iterations: passes of guarded smoothing; 0 is none.
    :param smooth_method: ``"laplacian"`` or ``"odt"``.
    :param return_report: also return the report.
    :returns: the surface mesh (``triangle`` and ``line`` blocks), or
        ``(mesh, report)``. The report has the counts (``num_failed``,
        ``num_open_edges``, ``num_non_manifold_edges``,
        ``num_misoriented_edges``, ...), the extreme angles and edge lengths,
        ``max_edge_gap`` (how far two uses of one edge disagree), ``faces``
        (one dict per face: ``status``, ``reason``, counts and extremes,
        ``below_bound``, ``bound_deg``, ``min_input_angle_deg`` and its
        ``corners``), ``warnings``, refinement's ``num_refinement_points``,
        ``num_faces_below_bound`` and ``num_smoothed``, and
        ``angle_histogram`` (every angle on the surface in 18 bins of 10
        degrees).
    :raises ValueError: on an option out of range or a mesh with no B-rep.
    :raises RuntimeError: with ``strict=True``, naming the first failed face.
    :raises NotImplementedError: without the compiled core.
    """
    brep = _brep_of(mesh)
    target_size = float(target_size)
    if not (math.isfinite(target_size) and target_size > 0.0):
        raise ValueError(f"{_PREFIX}target_size must be positive and finite")
    for name, value in (
        ("max_chordal_deviation", max_chordal_deviation),
        ("tolerance", tolerance),
    ):
        value = float(value)
        if not (math.isfinite(value) and value >= 0.0):
            raise ValueError(f"{_PREFIX}{name} must be finite and not negative")
    if int(min_points_per_edge) < 0:
        raise ValueError(f"{_PREFIX}min_points_per_edge must not be negative")
    if int(max_points) < 1:
        raise ValueError(f"{_PREFIX}max_points must be positive")
    min_angle = float(min_angle)
    if not (math.isfinite(min_angle) and 0.0 <= min_angle < 60.0):
        raise ValueError(f"{_PREFIX}min_angle must lie in [0, 60) degrees")
    if int(smooth_iterations) < 0:
        raise ValueError(f"{_PREFIX}smooth_iterations must not be negative")
    if smooth_method not in SMOOTH_METHODS:
        raise ValueError(f"{_PREFIX}smooth_method must be one of {SMOOTH_METHODS}")

    _core = _core_or_raise("mesh_cad_surface")
    res = _core.mesh_cad_surface(
        brep,
        target_size,
        max_chordal_deviation=float(max_chordal_deviation),
        min_points_per_edge=int(min_points_per_edge),
        tolerance=float(tolerance),
        record_face_ids=bool(record_face_ids),
        edges=bool(edges),
        match_unreferenced_trims=bool(match_unreferenced_trims),
        strict=bool(strict),
        max_points=int(max_points),
        min_angle=min_angle,
        smooth_iterations=int(smooth_iterations),
        smooth_method=smooth_method,
    )
    out = res["mesh"]
    if not return_report:
        return out
    report = {k: int(res[k]) for k in _REPORT_INTS}
    report.update({k: float(res[k]) for k in _REPORT_FLOATS})
    report["faces"] = [
        {**f, "corners": [dict(c) for c in f["corners"]]} for f in res["faces"]
    ]
    report["warnings"] = [str(w) for w in res["warnings"]]
    report["angle_histogram"] = [int(b) for b in res["angle_histogram"]]
    return out, report


def triangulate_polygon_2d(
    outer,
    holes=(),
    target_size: float = 0.0,
    *,
    min_angle: float = 0.0,
    max_points: int = 1_000_000,
):
    """The constrained Delaunay triangulation of a polygon with holes.

    :param outer: the boundary, ``(n, 2)``, a closed polyline whose last
        point is not repeated (either orientation).
    :param holes: the holes, each ``(m, 2)``; the even-odd rule decides the
        domain, so their orientation does not matter.
    :param target_size: 0 keeps the boundary as given and adds no interior
        point; ``h > 0`` splits boundary segments longer than ``h`` and fills
        the inside with an equilateral lattice of spacing ``h``.
    :param min_angle: refine (Ruppert's algorithm, splitting the boundary)
        until every angle is at least this (degrees, under 60); 0 keeps the
        unrefined triangulation. Corners under 60 degrees keep the triangles
        they force.
    :param max_points: the most points refinement adds.
    :returns: a ``triangle`` mesh in the plane z = 0, triangles
        counter-clockwise.
    :raises ValueError: for a negative size, a ``min_angle`` outside
        ``[0, 60)``, non-finite points, a coordinate
        outside the range the exact predicates hold in (zero, or a magnitude
        in ``[2**-128, 2**128]``), or a boundary that is degenerate, open or
        crosses itself.
    :raises NotImplementedError: without the compiled core.
    """
    loops = []
    for i, loop in enumerate([outer, *holes]):
        a = np.asarray(loop, dtype=np.float64)
        if a.ndim != 2 or a.shape[1] != 2:
            raise ValueError(
                f"meshio++: triangulate_polygon_2d: loop {i} must be an (n, 2) array"
            )
        loops.append(np.ascontiguousarray(a))
    target_size = float(target_size)
    if not (math.isfinite(target_size) and target_size >= 0.0):
        raise ValueError(
            "meshio++: triangulate_polygon_2d: target_size must be finite and not negative"
        )
    min_angle = float(min_angle)
    if not (math.isfinite(min_angle) and 0.0 <= min_angle < 60.0):
        raise ValueError(
            "meshio++: triangulate_polygon_2d: min_angle must lie in [0, 60) degrees"
        )
    if int(max_points) < 0:
        raise ValueError(
            "meshio++: triangulate_polygon_2d: max_points must not be negative"
        )
    _core = _core_or_raise("triangulate_polygon_2d")
    return _core.triangulate_polygon_2d(
        loops, target_size, min_angle=min_angle, max_points=int(max_points)
    )


__all__ = ["mesh_cad_surface", "triangulate_polygon_2d"]
