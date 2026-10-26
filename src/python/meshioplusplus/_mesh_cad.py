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

#: The most rings of points ``pole_rings`` seeds around a pole.
_MAX_POLE_RINGS = 8

#: How smoothing moves an interior point.
SMOOTH_METHODS = ("laplacian", "odt")
METRICS = ("mean", "varying")

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
    order: int = 1,
    curvature_deviation: float = 0.0,
    cells_across_gaps: float = 0.0,
    grading: float = 0.0,
    pole_rings: int = 0,
    face_sizes=None,
    edge_sizes=None,
    background=None,
    background_array: str = "size",
    virtual_edges: bool = False,
    metric: str = "mean",
    reparam_condition_number: float = 0.0,
    model_tolerances: bool = False,
    virtual_faces: bool = False,
    lop_corners: bool = False,
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
    :param order: 1 for linear triangles; 2 or 3 for curved ``triangle6`` or
        ``triangle10`` elements (and ``line3``, ``line4``) whose nodes lie on
        the surface and the edges. An element whose Jacobian is not positive
        everywhere inside is written straight-sided and counted
        (``num_curved_straightened``).
    :param curvature_deviation: size each region by its curvature so that a
        chord sits at most this far from the surface, ``2 sqrt(2 d / k - d^2)``
        for the largest principal curvature ``k``; 0 is off. Never above
        ``target_size`` nor below a 64th of it.
    :param cells_across_gaps: size each region so that at least this many
        elements span its local feature size (a thin wall's thickness, a
        slot's width); 0 is off.
    :param grading: the most the size may grow from one element to the next,
        ``|grad h| <= grading - 1``; 0 is 1.5. Acts only with a size source on.
    :param pole_rings: rings of points seeded around each pole (a degenerate
        trim), at about ``target_size`` apart, so the pole is no longer the
        apex of a fan of thin triangles; 0 keeps the fan.
    :param face_sizes: sizes given per face, ``{face index: size}`` or a
        sequence of ``(face index, size)`` pairs (indices into
        ``brep["faces"]``); every point of the face is at most that size.
    :param edge_sizes: the same per edge (indices into ``brep["edges"]``).
    :param background: a mesh whose one-component ``point_data``
        ``background_array`` is the size at its points, read at its nearest
        surface point (a volume mesh is skinned first).
    :param background_array: the name of that array.
    :param virtual_edges: sample chains of short edges meeting smoothly as
        one, so a sliver of an edge does not force tiny triangles; the chains
        are reported, the model is not changed.
    :param metric: ``"mean"`` meshes each face in its mean metric;
        ``"varying"`` follows the surface's metric as it varies over the face.
        Given sizes, like the curvature and gap sizes, join one size field
        graded by ``grading``; they only refine: nothing is coarser than
        ``target_size``.
    :param reparam_condition_number: mesh a face whose metric's condition
        number over its unrefined triangulation (the largest stretch of the
        surface against the face's plane over the smallest) exceeds this in
        a flattening of that triangulation (a conformal map, or mean-value
        coordinates on a circle) instead of its own parameters, every point
        still lifted onto the surface; 0 is off, otherwise above 1. Order 1
        only.
    :param model_tolerances: accept the gaps the model itself declares, the
        ``"tolerance"`` and ``"end_tolerances"`` in an edge's ``extra`` JSON
        (a STEP or IGES read writes OpenCASCADE's), where a closed edge's
        uses start and how far its uses and corners may lie apart.
    :param virtual_faces: mesh a face narrower than ``target_size`` together
        with the neighbour it meets tangentially along the most edge, through
        their common flattening; each triangle keeps the face it lies in, and
        the edges they share leave the ``line`` block. Order 1 only.
    :param lop_corners: protect every acute corner by corner lopping (S. E.
        Pav and N. J. Walkington, "Delaunay refinement by corner lopping",
        IMR 2005) in place of concentric shells: a ball of radius
        ``(sqrt 2 - 1)`` times the corner's distance to the nearest other
        corner of the face, its edges split where they leave it, halved while
        a vertex lies inside, and refinement leaves the triangles whose point
        falls in it. Each corner is reported with the angle bound the paper
        proves (Corollary 3) for the face's turning on a planar face, and
        whether the run meets it. Only with refinement (``min_angle`` or a
        size).
    :param return_report: also return the report.
    :returns: the surface mesh (``triangle`` and ``line`` blocks, or their
        curved forms with ``order`` 2 or 3), or
        ``(mesh, report)``. The report has the counts (``num_failed``,
        ``num_open_edges``, ``num_non_manifold_edges``,
        ``num_misoriented_edges``, ...), the extreme angles and edge lengths,
        ``max_edge_gap`` (how far two uses of one edge disagree), ``faces``
        (one dict per face: ``status``, ``reason``, counts and extremes,
        ``below_bound``, ``bound_deg``, ``min_input_angle_deg`` and its
        ``corners``), ``warnings``, refinement's ``num_refinement_points``,
        ``num_faces_below_bound`` and ``num_smoothed``,
        ``num_curved_straightened``,
        ``angle_histogram`` (every angle on the surface in 18 bins of 10
        degrees), ``virtual_edges`` (the edge indices of each chain) and, per
        face, ``limited_by_edges`` (the edges whose sampling is finer than the
        face's given size allows), ``reparametrised`` (per face over
        ``reparam_condition_number``: ``face``, ``condition_number``, the
        flattening's ``condition_number_disc``, ``taken``, ``method``
        ``"conformal"`` or ``"mean_value"``, and ``reason`` when not taken),
        ``merged_faces`` (per sliver found: ``sliver``, ``host``, ``width``,
        the shared ``edges``, ``merged`` and ``reason``) and
        ``tolerance_uses`` (per edge and place where a declared tolerance
        accepted a gap: ``edge``, ``declared``, ``gap``, ``where``) and
        ``lopped_corners`` (per corner lopped: ``face``, ``point``,
        ``angle_deg``, ``lfs0``, ``lop_radius``, ``turning_deg``,
        ``stated_bound_deg``, ``num_forced``, ``min_forced_deg``,
        ``within_two_radii``, ``guaranteed`` and ``reason`` when not).
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
    if isinstance(order, bool) or int(order) not in (1, 2, 3):
        raise ValueError(f"{_PREFIX}order must be 1, 2 or 3")
    for name, value in (
        ("curvature_deviation", curvature_deviation),
        ("cells_across_gaps", cells_across_gaps),
    ):
        value = float(value)
        if not (math.isfinite(value) and value >= 0.0):
            raise ValueError(f"{_PREFIX}{name} must be finite and not negative")
    grading = float(grading)
    if not (math.isfinite(grading) and (grading == 0.0 or grading > 1.0)):
        raise ValueError(
            f"{_PREFIX}grading must be 0 (the default, 1.5) or greater than 1"
        )
    if isinstance(pole_rings, bool) or not 0 <= int(pole_rings) <= _MAX_POLE_RINGS:
        raise ValueError(f"{_PREFIX}pole_rings must lie in [0, {_MAX_POLE_RINGS}]")
    reparam_condition_number = float(reparam_condition_number)
    if not (
        math.isfinite(reparam_condition_number)
        and (reparam_condition_number == 0.0 or reparam_condition_number > 1.0)
    ):
        raise ValueError(
            f"{_PREFIX}reparam_condition_number must be 0 (off) or greater than 1"
        )
    if reparam_condition_number > 0.0 and int(order) > 1:
        raise ValueError(f"{_PREFIX}reparam_condition_number is order 1 only")
    if virtual_faces and int(order) > 1:
        raise ValueError(f"{_PREFIX}virtual_faces is order 1 only")

    face_pairs = _size_pairs(face_sizes, "face_sizes", len(brep["faces"]))
    edge_pairs = _size_pairs(edge_sizes, "edge_sizes", len(brep.get("edges", [])))
    if metric not in METRICS:
        raise ValueError(f"{_PREFIX}metric must be one of {METRICS}")
    if background is not None:
        data = getattr(background, "point_data", {})
        if background_array not in data:
            raise ValueError(
                f"{_PREFIX}background: it has no point_data {background_array!r}"
            )
        values = np.asarray(data[background_array])
        if values.ndim != 1 and not (values.ndim == 2 and values.shape[1] == 1):
            raise ValueError(
                f"{_PREFIX}background: its point_data {background_array!r} must "
                "have one component"
            )

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
        order=int(order),
        grading=grading,
        curvature_deviation=float(curvature_deviation),
        cells_across_gaps=float(cells_across_gaps),
        pole_rings=int(pole_rings),
        face_sizes=face_pairs,
        edge_sizes=edge_pairs,
        background=background,
        background_array=str(background_array),
        virtual_edges=bool(virtual_edges),
        metric=metric,
        reparam_condition_number=reparam_condition_number,
        model_tolerances=bool(model_tolerances),
        virtual_faces=bool(virtual_faces),
        lop_corners=bool(lop_corners),
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
    report["virtual_edges"] = [[int(e) for e in c] for c in res["virtual_edges"]]
    for f in report["faces"]:
        f["limited_by_edges"] = [int(e) for e in f["limited_by_edges"]]
    report["reparametrised"] = [dict(r) for r in res["reparametrised"]]
    report["merged_faces"] = [
        {**m, "edges": [int(e) for e in m["edges"]]} for m in res["merged_faces"]
    ]
    report["tolerance_uses"] = [dict(t) for t in res["tolerance_uses"]]
    report["lopped_corners"] = [dict(c) for c in res["lopped_corners"]]
    return out, report


def _size_pairs(sizes, name, count):
    """``{index: size}`` or ``(index, size)`` pairs as sorted pairs, each
    index in ``[0, count)`` and each size finite and positive."""
    if sizes is None:
        return []
    items = sizes.items() if isinstance(sizes, dict) else sizes
    pairs = []
    for item in items:
        try:
            index, size = item
        except (TypeError, ValueError):
            raise ValueError(f"{_PREFIX}{name} must hold (index, size) pairs") from None
        if isinstance(index, bool) or int(index) != index:
            raise ValueError(f"{_PREFIX}{name}: index {index!r} is not an integer")
        index, size = int(index), float(size)
        if not 0 <= index < count:
            raise ValueError(
                f"{_PREFIX}{name}: index {index} is out of range [0, {count})"
            )
        if not (math.isfinite(size) and size > 0.0):
            raise ValueError(
                f"{_PREFIX}{name}: the size of {index} must be finite and positive"
            )
        pairs.append((index, size))
    return sorted(pairs)


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
