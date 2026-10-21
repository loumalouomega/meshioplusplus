"""Turn the fitted patches of a closed skin into a CAD solid.

The last step of the skin-to-CAD chain :func:`segment` -> :func:`fit_primitives`
-> :func:`fit_cad`, in the core (``operations/fit_cad.hpp``): through the
optional OpenCASCADE backend when the core was built with
``-DMESHIOPLUSPLUS_WITH_OCCT=ON``, else through the kernel-free native backend,
which writes the fitted patches as analytic STEP faces without trimming them.
There is no numpy twin of either: without the compiled core the function
raises ``NotImplementedError``.
"""

from __future__ import annotations

import math
import operator
import os

from ._fit_primitives import TYPES

_PREFIX = "meshio++: fit_cad: "
_FLAG = "-DMESHIOPLUSPLUS_WITH_OCCT=ON"

_REPORT_FLOATS = (
    "tolerance",
    "volume",
    "skin_volume",
    "volume_ratio",
    "max_deviation",
)
_REPORT_INTS = (
    "num_faces",
    "num_plane",
    "num_sphere",
    "num_cylinder",
    "num_cone",
    "num_torus",
    "num_faceted_faces",
    "num_patches",
    "num_demoted",
    "num_exact_edges",
    "num_approx_edges",
    "num_solids",
    "num_bspline",
    "num_split_patches",
    "num_fillet",
    "num_helical",
    "num_variable_fillet",
    "num_blend",
    "num_rounded_thread",
    "num_vertex_blend",
    "num_chamfer",
    "num_conic_blend",
    "num_tapered_thread",
    "num_variable_pitch_thread",
)

# The refinement's report keys, added only when ``refine_segmentation`` is set.
_REFINE_INTS = ("refine_iterations", "num_refine_splits", "num_refine_merges")
_REFINE_FLOATS = ("initial_explained_area_fraction", "explained_area_fraction")
# The general freeform fit's, only with ``freeform_general``; the extended
# surfaces, only with ``extend_surfaces``; the snapped junctions, only with
# ``snap_junctions``; the Hausdorff distances, only with ``face_deviation``.
_FREEFORM_INTS = ("num_freeform_general", "num_freeform_refinements")
_FREEFORM_FLOATS = ("avg_freeform_poles",)
_DEVIATION_FLOATS = ("hausdorff_max", "hausdorff_rms")

# The Phase 3 options and their defaults: passed to the core only when one
# differs, so a default call is the one the binding always took.
_PHASE3_DEFAULTS = {
    "freeform_general": False,
    "freeform_tolerance": 0.0,
    "freeform_smoothing": 1e-12,
    "max_freeform_poles": 40,
    "max_freeform_iterations": 30,
    "max_freeform_splits": 2,
    "sew_tolerance": 0.0,
    "extend_surfaces": False,
    "snap_junctions": False,
    "face_deviation": False,
}

# The freeform-mode options (v17.12.0) and their defaults, passed to the core
# only when one differs; their report keys appear only outside the legacy mode.
_FREEFORM_MODES = ("legacy", "adaptive", "network", "auto")
_PARAMETERIZATIONS = ("auto", "projection", "lscm", "harmonic")
_EXTRAORDINARY_VERTEX = ("g0", "g1", "g1x", "g1xx")
_LOCAL_REFINE = ("knots", "thb")
_MODE_DEFAULTS = {
    "freeform_mode": "legacy",
    "freeform_parameterization": "auto",
    "extraordinary_vertex": "g1x",
    "local_refine": "knots",
    "freeform_g1_tolerance_deg": 1.0,
    "freeform_max_control_points": 400,
    "layout_target_patch_size": 0.1,
    "layout_max_valence_deviation": 2,
    "resegment_rounds": 2,
    "layout_snap_to_features": True,
    "validate_by_surface_sampling": True,
}
_MODE_INTS = (
    "num_freeform_adaptive",
    "num_network_patches",
    "num_network_edges",
    "num_extraordinary_vertices",
    "num_resegmented",
    "num_escalations",
)
_MODE_FLOATS = (
    "facet_area_fraction",
    "max_g1_angle_deg",
    "max_freeform_distance",
    "rms_freeform_distance",
)


def has_occt() -> bool:
    """Whether the compiled core carries the OpenCASCADE backend."""
    try:
        from . import _core
    except ImportError:
        return False
    return bool(getattr(_core, "__has_occt__", False))


def _types_mask(types, prefix: str = _PREFIX) -> int:
    """The ``1 << type`` bit mask of primitive type names (``None``: all of them)."""
    if types is None:
        return sum(1 << t for t in range(1, 6))
    if isinstance(types, str):
        types = [types]
    mask = 0
    for t in types:
        if t not in TYPES[1:]:
            raise ValueError(
                f"{prefix}unknown primitive type {t!r} "
                f"(expected one of {', '.join(TYPES[1:])})"
            )
        mask |= 1 << TYPES.index(t)
    return mask


#: The default seed of the refinement's RANSAC split: the bytes of "meshio++".
DEFAULT_RANSAC_SEED = 0x6D657368696F2B2B


def _check_seed(seed, prefix: str) -> int:
    """``seed`` as an int in ``[0, 2**64)``, else ``ValueError``."""
    if isinstance(seed, bool):
        raise ValueError(f"{prefix}ransac_seed must be an integer")
    try:
        value = operator.index(seed)
    except TypeError:
        raise ValueError(f"{prefix}ransac_seed must be an integer") from None
    if not 0 <= value < 1 << 64:
        raise ValueError(f"{prefix}ransac_seed must lie in [0, 2**64)")
    return value


def fit_cad(
    mesh,
    path=None,
    *,
    tolerance: float | None = None,
    exact_edges: bool = True,
    types=None,
    mesh_deflection: float | None = None,
    max_demotion_rounds: int = 8,
    freeform: bool = True,
    freeform_degree: int = 3,
    freeform_max_control: int = 64,
    freeform_split: bool = True,
    fillets: bool = True,
    threads: bool = True,
    return_report: bool = False,
    merge_faces: bool = True,
    backend: str = "auto",
    refine_segmentation: bool = False,
    max_refine_iterations: int = 10,
    refine_epsilon: float = 1e-3,
    normal_weight: float = 1.0,
    allow_merge: bool = True,
    allow_split: bool = True,
    ransac_seed: int = DEFAULT_RANSAC_SEED,
    refine_tolerance: float = 1e-3,
    relative_tolerance: bool = True,
    freeform_general: bool = False,
    freeform_tolerance: float = 0.0,
    freeform_smoothing: float = 1e-12,
    max_freeform_poles: int = 40,
    max_freeform_iterations: int = 30,
    max_freeform_splits: int = 2,
    sew_tolerance: float = 0.0,
    extend_surfaces: bool = False,
    snap_junctions: bool = False,
    face_deviation: bool = False,
    freeform_mode: str = "legacy",
    freeform_parameterization: str = "auto",
    extraordinary_vertex: str = "g1x",
    local_refine: str = "knots",
    freeform_g1_tolerance_deg: float = 1.0,
    freeform_max_control_points: int = 400,
    layout_target_patch_size: float = 0.1,
    layout_max_valence_deviation: int = 2,
    resegment_rounds: int = 2,
    layout_snap_to_features: bool = True,
    validate_by_surface_sampling: bool = True,
):
    """Build a CAD solid from the fitted patches of a closed skin.

    Every patch :func:`fit_primitives` accepted becomes one plane, sphere,
    cylinder, cone or torus face trimmed by the patch boundary. A patch no
    primitive explains is tried, with ``fillets``, as a rolling-ball fillet
    when it runs between two analytic faces (of one radius, else of a radius
    varying along it), and, with ``threads``, as a helical thread face: about
    the axis of a cylinder or cone beside it, or of the screw the patch itself
    fits, with a straight flank or any other profile (a rounded root or
    crest). Failing those it becomes, with ``freeform``, a B-spline face
    fitted over the patch flattened onto its best-fit plane or unrolled around
    a cylinder, and one both fold is, with ``freeform_split``, split into
    pieces that flatten, each its own face under a new patch label. A patch
    between two fitted faces, one of them freeform, is then tried, with
    ``fillets``, as a rolling-ball blend. The rest of the skin becomes planar
    facets. The faces are sewn into a solid that
    OpenCASCADE checks, and a face that fails the check or strays further than
    ``tolerance`` from the skin is demoted to facets, so a feature the fitter
    cannot explain comes back faceted, never silently wrong.

    :param mesh: a closed, consistently wound surface mesh, normally the
        output of :func:`fit_primitives`. Without ``segment:patch``,
        :func:`segment` runs first with ``method="primitive"`` (so a fillet
        comes apart from the planes it is tangent to) and :func:`fit_primitives`
        at the same absolute tolerance; without the ``fit:patch_*`` table,
        :func:`fit_primitives` runs with its default options.
    :param path: write the solid here as AP214 STEP (millimetres, a pinned
        header timestamp); ``None`` writes nothing.
    :param tolerance: the largest distance a face may lie from the skin; by
        default three times the largest ``fit:patch_tolerance`` of an accepted
        patch (``1e-3`` of the bounding-box diagonal when none was accepted).
    :param exact_edges: try the exact intersection of two analytic surfaces
        before approximating the edge between them with a B-spline.
    :param types: which accepted types become analytic faces, from ``plane``,
        ``sphere``, ``cylinder``, ``cone``, ``torus`` (default all); an empty
        sequence facets everything.
    :param mesh_deflection: the chordal deflection of the re-meshed solid; by
        default a quarter of the tolerance.
    :param max_demotion_rounds: how many times failed faces are demoted and
        the solid rebuilt before every patch is faceted.
    :param freeform: try a patch no primitive explains as a B-spline face.
    :param freeform_degree: that B-spline's highest degree, 3 to 8.
    :param freeform_max_control: the most grid points per direction its fit
        samples, 4 to 1024.
    :param freeform_split: split a freeform patch that folds over both its
        flattenings into pieces that do not.
    :param fillets: try a patch between two fitted faces as a rolling-ball
        fillet or blend, then as a chamfer or a conic blend, and a corner among
        three or more blends as a vertex blend.
    :param threads: try a patch as a helical thread face, of constant or
        varying pitch and cut on a taper.
    :param return_report: also return the report.
    :param merge_faces: once the solid is built, merge faces on one surface
        and join collinear edges (``ShapeUpgrade_UnifySameDomain``), kept only
        when the solid stays valid with the same volume.
    :param backend: ``"occt"`` builds with OpenCASCADE (trimmed, sewn and
        checked faces, freeform, fillet and thread faces, the re-meshed solid);
        ``"native"`` needs no CAD kernel: every fitted plane, sphere, cylinder,
        cone and torus patch becomes one analytic STEP face, bounded by line,
        circle and polyline edges, and the rest coplanar polygons (the
        returned mesh is the skin itself with the ``fit_cad:*`` cell data, its
        volumes NaN, and ``brep`` empty in the legacy freeform mode and the
        solid as trimmed NURBS faces outside it); ``"auto"`` (the default) is
        ``"occt"`` when the core carries it, else ``"native"``.
    :param refine_segmentation: first refine the patches with
        :func:`refine_patches` (variational shape approximation: patch
        boundaries move, freeform patches are split by RANSAC, neighbours one
        primitive explains merge), after any segmentation and fit ``fit_cad``
        runs itself, so more of the skin becomes analytic faces. Off by
        default, so the output is what it was before the option existed.
    :param max_refine_iterations: the refinement's most iterations, ``>= 0``.
    :param refine_epsilon: stop refining when an iteration adds less than this
        to the share of the area primitives explain, ``>= 0``.
    :param normal_weight: the weight of the normal term in the refinement's
        cell cost, ``>= 0``.
    :param allow_merge: let the refinement merge neighbouring patches one
        primitive explains.
    :param allow_split: let the refinement split a primitive out of a patch
        no primitive explains.
    :param ransac_seed: the seed of the refinement's RANSAC split, in
        ``[0, 2**64)``.
    :param refine_tolerance: the distance every point of a refined patch must
        lie within from its proxy.
    :param relative_tolerance: ``refine_tolerance`` is a fraction of the
        skin's bounding-box diagonal (default) rather than a length.
    :param freeform_general: try a patch the plane and unrolled-cylinder
        freeform fit leaves unexplained as a general freeform B-spline face: a
        disc parametrisation (least-squares conformal, else mean-value) and an
        iterative least-squares fit with parameter correction and knot
        insertion, the patch halved up to ``max_freeform_splits`` times when
        that fails. Off by default.
    :param freeform_tolerance: the general fit's target at the vertices,
        ``>= 0``; 0 takes a quarter of the tolerance. The face is accepted by
        the same test as any freeform face, so a larger target never loosens
        it.
    :param freeform_smoothing: the general fit's thin-plate weight, relative
        to the sample count on a unit parameter square, ``>= 0``.
    :param max_freeform_poles: the most poles per direction of a general
        freeform face, 4 to 1024.
    :param max_freeform_iterations: the general fit's iteration budget,
        ``>= 1``.
    :param max_freeform_splits: how many times a patch the general fit cannot
        explain is halved, ``>= 0``.
    :param sew_tolerance: the tolerance the faces are sewn at, ``>= 0``; 0
        (the default) sews at the tolerance itself. OpenCASCADE only.
    :param extend_surfaces: extend every freeform B-spline surface past its
        patch before the edges are built, so an intersection with a neighbour
        just outside the fitted area is still found. OpenCASCADE only.
    :param snap_junctions: place the B-rep vertices where analytic faces meet
        by a least-squares solve on the surfaces' tangent planes before the
        averaged projection. OpenCASCADE only.
    :param face_deviation: measure the two-sided distance between each face
        and the skin, as the field data ``fit_cad:face_hausdorff`` (per face,
        an ``(n, 2)`` array: the face's re-mesh to the skin, and the skin it
        stands for to the solid) and ``fit_cad:face_rms``. OpenCASCADE only.
    :param freeform_mode: how the patches no primitive or feature fit explains
        are treated (v17.12.0): ``"legacy"`` (the default) keeps every output
        byte as before; ``"adaptive"`` runs the general freeform fit with the
        parametrisation ladder, surface-sampling validation and the escalation
        ladder before demotion; ``"network"`` fits a quad-layout B-spline patch
        network to every residual region; ``"auto"`` is adaptive first and a
        network for what still fails. Every mode but ``"legacy"`` implies
        ``freeform_general``.
    :param freeform_parameterization: ``"auto"`` (projection, then
        least-squares conformal maps, then mean-value coordinates, then a
        split), ``"projection"``, ``"lscm"`` or ``"harmonic"``.
    :param extraordinary_vertex: a patch network's continuity at an
        extraordinary vertex: ``"g0"``, ``"g1"``, ``"g1x"`` (the default) or
        ``"g1xx"``.
    :param local_refine: ``"knots"`` (full knot lines, the default) or
        ``"thb"`` (truncated hierarchical B-splines, flattened; experimental).
    :param freeform_g1_tolerance_deg: the largest normal jump, in degrees,
        across an edge between two network faces, ``> 0``.
    :param freeform_max_control_points: the most poles of one freeform face
        outside the legacy mode, 16 to 1048576.
    :param layout_target_patch_size: a quad-layout patch's target size,
        relative to the bounding-box diagonal, ``> 0``.
    :param layout_max_valence_deviation: the largest departure from valence 4
        of a layout vertex, 1 to 3.
    :param resegment_rounds: re-segmentation rounds of a residual region in
        the escalation ladder, ``>= 0``.
    :param layout_snap_to_features: snap layout edges to feature edges.
    :param validate_by_surface_sampling: outside the legacy mode, accept a
        freeform face only once samples of its surface lie within the
        tolerance of the skin too.
    :returns: the solid re-meshed as triangles with cell data ``fit_cad:face``,
        ``fit_cad:type`` (0 a facet, 6 a freeform B-spline, 7 a fillet,
        blend or conic blend, 8 a helical thread face, 9 a vertex blend, 10 a
        chamfer, else as ``fit:type``) and
        ``fit_cad:patch``, and the solid as NURBS faces,
        loops and edges in ``mesh.info["brep"]`` (face ``i`` of
        ``fit_cad:face`` has id ``i + 1``), so writing it as ``.cad.json`` gives
        Kratos CAD JSON. With ``return_report``, ``(mesh, report)``, the report
        holding ``valid``, the tolerance, the solid's and the skin's volume and
        their ratio, the Hausdorff ``max_deviation``, the faces by type
        (``num_bspline`` the freeform ones, ``num_fillet``, ``num_helical``),
        of those the fillets of varying radius (``num_variable_fillet``), the
        blends beside a freeform face (``num_blend``) and the thread faces of
        a profile that is not straight (``num_rounded_thread``), the vertex
        blends (``num_vertex_blend``), chamfers (``num_chamfer``), conic blends
        (``num_conic_blend``) and tapered and variable-pitch thread faces
        (``num_tapered_thread``, ``num_variable_pitch_thread``),
        ``num_split_patches``, the demoted patches and the exact and
        approximated edges. With ``refine_segmentation`` the report also holds
        ``refined``, ``refinement_kept`` (whether the solid was built from the
        refined patches: the OpenCASCADE backend also builds it from the
        patches as given and keeps that one when it is better),
        ``refine_iterations``, ``num_refine_splits``,
        ``num_refine_merges`` and the area share primitives explain before and
        after (``initial_explained_area_fraction``,
        ``explained_area_fraction``), and the mesh the ``refine:*`` field data
        of :func:`refine_patches`. The report always holds ``closed`` (every
        shell's edges border two faces); with ``freeform_general`` also
        ``num_freeform_general``, ``num_freeform_refinements`` and
        ``avg_freeform_poles``; with ``extend_surfaces``
        ``num_extended_edges``; with ``snap_junctions``
        ``num_snapped_junctions``; and with ``face_deviation`` the two-sided
        ``hausdorff_max`` and ``hausdorff_rms`` (NaN when not measured, as from
        the native backend). Outside the legacy ``freeform_mode`` the report
        also holds ``freeform_mode``, ``brep_checked`` (whether OpenCASCADE's
        B-rep check ran), ``facet_area_fraction`` (the share of the result's
        area left on facet faces), ``num_freeform_adaptive``,
        ``num_network_patches``, ``num_network_edges``,
        ``num_extraordinary_vertices``, ``num_resegmented``,
        ``num_escalations``, ``max_g1_angle_deg``, ``max_freeform_distance``,
        ``rms_freeform_distance`` and ``faces``: one dict per patch tried as a
        surface (``patch``, ``face``, ``kind``, ``param_method``, ``num_ctrl``,
        ``iterations``, ``refinements``, ``max_dist``, ``rms_dist``,
        ``reverse_max_dist``, ``max_g1_angle_deg``, ``escalation_steps`` and
        ``demotion_reason``), which the mesh also carries as the
        ``fit_cad:face_table`` field data.
    :raises NotImplementedError: without the compiled core, or with
        ``backend="occt"`` when the core was built without OpenCASCADE.
    """
    opts = {
        "path": "" if path is None else os.fspath(path),
        "tolerance": 0.0 if tolerance is None else float(tolerance),
        "exact_edges": bool(exact_edges),
        "types": _types_mask(types),
        "mesh_deflection": 0.0 if mesh_deflection is None else float(mesh_deflection),
        "max_demotion_rounds": int(max_demotion_rounds),
        "freeform": bool(freeform),
        "freeform_degree": int(freeform_degree),
        "freeform_max_control": int(freeform_max_control),
        "freeform_split": bool(freeform_split),
        "fillets": bool(fillets),
        "threads": bool(threads),
        "merge_faces": bool(merge_faces),
        "backend": str(backend),
    }
    if opts["backend"] not in ("auto", "occt", "native"):
        raise ValueError(f"{_PREFIX}backend must be 'auto', 'occt' or 'native'")
    if tolerance is not None and not (
        opts["tolerance"] > 0.0 and math.isfinite(opts["tolerance"])
    ):
        raise ValueError(f"{_PREFIX}the tolerance must be finite and > 0")
    if mesh_deflection is not None and not (
        opts["mesh_deflection"] > 0.0 and math.isfinite(opts["mesh_deflection"])
    ):
        raise ValueError(f"{_PREFIX}mesh_deflection must be finite and > 0")
    if opts["max_demotion_rounds"] < 0:
        raise ValueError(f"{_PREFIX}max_demotion_rounds must be >= 0")
    if not 3 <= opts["freeform_degree"] <= 8:
        raise ValueError(f"{_PREFIX}freeform_degree must be 3 to 8")
    if not 4 <= opts["freeform_max_control"] <= 1024:
        raise ValueError(f"{_PREFIX}freeform_max_control must be 4 to 1024")
    refine = {
        "refine_segmentation": bool(refine_segmentation),
        "refine_tolerance": float(refine_tolerance),
        "relative_tolerance": bool(relative_tolerance),
        "normal_weight": float(normal_weight),
        "refine_epsilon": float(refine_epsilon),
        "max_refine_iterations": int(max_refine_iterations),
        "allow_split": bool(allow_split),
        "allow_merge": bool(allow_merge),
        "ransac_seed": _check_seed(ransac_seed, _PREFIX),
    }
    if not (
        refine["refine_tolerance"] > 0.0 and math.isfinite(refine["refine_tolerance"])
    ):
        raise ValueError(f"{_PREFIX}refine_tolerance must be finite and > 0")
    if not (refine["normal_weight"] >= 0.0 and math.isfinite(refine["normal_weight"])):
        raise ValueError(f"{_PREFIX}normal_weight must be finite and >= 0")
    if not (
        refine["refine_epsilon"] >= 0.0 and math.isfinite(refine["refine_epsilon"])
    ):
        raise ValueError(f"{_PREFIX}refine_epsilon must be finite and >= 0")
    if refine["max_refine_iterations"] < 0:
        raise ValueError(f"{_PREFIX}max_refine_iterations must be >= 0")
    phase3 = {
        "freeform_general": bool(freeform_general),
        "freeform_tolerance": float(freeform_tolerance),
        "freeform_smoothing": float(freeform_smoothing),
        "max_freeform_poles": int(max_freeform_poles),
        "max_freeform_iterations": int(max_freeform_iterations),
        "max_freeform_splits": int(max_freeform_splits),
        "sew_tolerance": float(sew_tolerance),
        "extend_surfaces": bool(extend_surfaces),
        "snap_junctions": bool(snap_junctions),
        "face_deviation": bool(face_deviation),
    }
    for key in ("freeform_tolerance", "freeform_smoothing", "sew_tolerance"):
        if not (phase3[key] >= 0.0 and math.isfinite(phase3[key])):
            raise ValueError(f"{_PREFIX}{key} must be finite and >= 0")
    if not 4 <= phase3["max_freeform_poles"] <= 1024:
        raise ValueError(f"{_PREFIX}max_freeform_poles must be 4 to 1024")
    if phase3["max_freeform_iterations"] < 1:
        raise ValueError(f"{_PREFIX}max_freeform_iterations must be >= 1")
    if phase3["max_freeform_splits"] < 0:
        raise ValueError(f"{_PREFIX}max_freeform_splits must be >= 0")
    modes = {
        "freeform_mode": str(freeform_mode),
        "freeform_parameterization": str(freeform_parameterization),
        "extraordinary_vertex": str(extraordinary_vertex),
        "local_refine": str(local_refine),
        "freeform_g1_tolerance_deg": float(freeform_g1_tolerance_deg),
        "freeform_max_control_points": int(freeform_max_control_points),
        "layout_target_patch_size": float(layout_target_patch_size),
        "layout_max_valence_deviation": int(layout_max_valence_deviation),
        "resegment_rounds": int(resegment_rounds),
        "layout_snap_to_features": bool(layout_snap_to_features),
        "validate_by_surface_sampling": bool(validate_by_surface_sampling),
    }
    for key, words in (
        ("freeform_mode", _FREEFORM_MODES),
        ("freeform_parameterization", _PARAMETERIZATIONS),
        ("extraordinary_vertex", _EXTRAORDINARY_VERTEX),
        ("local_refine", _LOCAL_REFINE),
    ):
        if modes[key] not in words:
            raise ValueError(
                f"{_PREFIX}{key} must be one of {', '.join(repr(w) for w in words)}"
            )
    for key in ("freeform_g1_tolerance_deg", "layout_target_patch_size"):
        if not (modes[key] > 0.0 and math.isfinite(modes[key])):
            raise ValueError(f"{_PREFIX}{key} must be finite and > 0")
    if not 16 <= modes["freeform_max_control_points"] <= 1048576:
        raise ValueError(f"{_PREFIX}freeform_max_control_points must be 16 to 1048576")
    if not 1 <= modes["layout_max_valence_deviation"] <= 3:
        raise ValueError(f"{_PREFIX}layout_max_valence_deviation must be 1 to 3")
    if modes["resegment_rounds"] < 0:
        raise ValueError(f"{_PREFIX}resegment_rounds must be >= 0")

    try:
        from . import _core
    except ImportError as exc:
        raise NotImplementedError(
            f"{_PREFIX}needs the compiled core; there is no pure-Python CAD kernel"
        ) from exc
    if opts["backend"] == "occt" and not getattr(_core, "__has_occt__", False):
        raise NotImplementedError(
            f"{_PREFIX}this build has no OpenCASCADE backend (the PyPI wheels carry "
            f"it); rebuild with {_FLAG} (OpenCASCADE 7.6 or later), or use "
            "backend='native'"
        )
    res = _core.fit_cad(
        mesh,
        opts["path"],
        opts["tolerance"],
        opts["exact_edges"],
        opts["types"],
        opts["mesh_deflection"],
        opts["max_demotion_rounds"],
        opts["freeform"],
        opts["freeform_degree"],
        opts["freeform_max_control"],
        opts["freeform_split"],
        opts["fillets"],
        opts["threads"],
        opts["merge_faces"],
        opts["backend"],
        # The defaults leave the core's refinement off; passed only when it
        # runs, so a default call is the one the binding always took.
        **(refine if refine["refine_segmentation"] else {}),
        **{k: v for k, v in phase3.items() if v != _PHASE3_DEFAULTS[k]},
        **{k: v for k, v in modes.items() if v != _MODE_DEFAULTS[k]},
    )
    out = res["mesh"]
    info = out.info if isinstance(getattr(out, "info", None), dict) else {}
    info["brep"] = res["brep"]
    out.info = info
    if not return_report:
        return out
    report = {
        "valid": bool(res["valid"]),
        "backend": str(res["backend"]),
        "closed": bool(res["closed"]),
    }
    report.update({k: float(res[k]) for k in _REPORT_FLOATS})
    report.update({k: int(res[k]) for k in _REPORT_INTS})
    if refine["refine_segmentation"]:
        report["refined"] = bool(res["refined"])
        report["refinement_kept"] = bool(res["refinement_kept"])
        report.update({k: int(res[k]) for k in _REFINE_INTS})
        report.update({k: float(res[k]) for k in _REFINE_FLOATS})
    if phase3["freeform_general"]:
        report.update({k: int(res[k]) for k in _FREEFORM_INTS})
        report.update({k: float(res[k]) for k in _FREEFORM_FLOATS})
    if phase3["extend_surfaces"]:
        report["num_extended_edges"] = int(res["num_extended_edges"])
    if phase3["snap_junctions"]:
        report["num_snapped_junctions"] = int(res["num_snapped_junctions"])
    if phase3["face_deviation"]:
        report.update({k: float(res[k]) for k in _DEVIATION_FLOATS})
    if modes["freeform_mode"] != "legacy":
        report["freeform_mode"] = str(res["freeform_mode"])
        report["brep_checked"] = bool(res["brep_checked"])
        report.update({k: int(res[k]) for k in _MODE_INTS})
        report.update({k: float(res[k]) for k in _MODE_FLOATS})
        report["faces"] = [dict(row) for row in res["faces"]]
    return out, report


__all__ = ["fit_cad", "has_occt"]
