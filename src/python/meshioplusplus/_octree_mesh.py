"""Mesh the inside of a skin on a balanced octree.

The volume mesher of ``operations/octree_mesh.hpp``: a size field (curvature,
sharp features, local thickness, an optional array) refines a linear octree
that is 2:1-balanced across faces, edges and corners, and inside and outside
come from the fast generalised winding number, so a slightly dirty skin (holes,
self-intersections, non-manifold edges) still meshes. The output is castellated
hexahedra with surrogate-boundary tags, conforming graded tetrahedra fitted to
the skin (the default), or all-hexahedra. There is no numpy twin: without the
compiled core the function raises ``NotImplementedError``. See
``doc/octree_mesh.md``.
"""

from __future__ import annotations

import math

_PREFIX = "meshio++: octree_mesh: "

#: The outputs, as ``output=`` names them.
OUTPUTS = ("castellated", "tet", "hex")
#: How a tetrahedral output conforms to the skin.
CONFORMITIES = ("stuffing", "cleaving")
#: How the size field may vary along the skin.
SKIN_GRADINGS = ("graded", "uniform")
#: Which leaves a castellated output keeps.
FILLS = ("all", "inside", "cover")

_REPORT_INTS = (
    "max_level",
    "finest_level",
    "surface_level",
    "num_leaves",
    "num_points",
    "num_cells",
    "num_hexahedra",
    "num_tetrahedra",
    "num_inside",
    "num_boundary_cells",
    "num_surrogate_faces",
    "num_hanging",
    "num_cut_points",
    "num_warped",
    "num_warp_refused",
    "num_inverted",
    "num_boundary_edges_open",
    "num_non_manifold_boundary_edges",
    "num_corners",
    "num_corners_snapped",
    "num_feature_chains",
    "num_feature_chains_missed",
    "num_feature_moves",
    "num_optimized_moves",
    "num_flips",
    "num_winding_rechecks",
)
_REPORT_FLOATS = (
    "root_size",
    "min_size",
    "max_size",
    "min_dihedral",
    "max_dihedral",
    "min_scaled_jacobian",
    "hausdorff_to_skin",
)


def _choice(name, value, allowed):
    value = str(value)
    if value not in allowed:
        raise ValueError(
            f"{_PREFIX}unknown {name} {value!r} (expected one of {allowed})"
        )
    return value


def _nonneg(name, value):
    value = float(value)
    if not (math.isfinite(value) and value >= 0.0):
        raise ValueError(f"{_PREFIX}{name} must be finite and >= 0")
    return value


def octree_mesh(
    mesh,
    *,
    output: str = "tet",
    min_size: float = 0.0,
    max_size: float = 0.0,
    max_level: int = -1,
    bounds=None,
    padding: float = 0.0,
    padding_relative: float = 0.1,
    curvature_angle: float = 15.0,
    feature_angle: float = 30.0,
    corner_angle: float = 30.0,
    feature_size: float = 0.0,
    gap_cells: float = 2.0,
    gap_outside: bool = False,
    size_array=None,
    surface_region=None,
    grading_layers: int = 1,
    skin_grading: str = "graded",
    surface_size: float = 0.0,
    winding_threshold: float = 0.5,
    winding_beta: float = 2.0,
    conformity: str = "stuffing",
    alpha_long: float = 0.24799,
    alpha_short: float = 0.40173,
    preserve_features: bool = True,
    fit_iterations: int = 5,
    dihedral_floor_guard: bool = True,
    fill: str = "inside",
    keep_outside: bool = False,
    attach_distance: bool = True,
    attach_quality: bool = False,
    pillow: bool = True,
    untangle_iterations: int = 50,
    max_cells: int = 20000000,
    max_tets: int = 100000000,
    return_report: bool = False,
):
    """Mesh the inside of a skin on a balanced octree.

    The skin is any surface (``triangle``, ``quad``, ``polygon`` blocks; the
    boundary of a volume mesh is taken), closed or slightly dirty. Leaves are
    no larger than the size field: ``max_size`` away from the skin, and on it
    the smallest of the curvature size (a leaf edge subtends at most
    ``curvature_angle`` degrees), ``feature_size`` along edges sharper than
    ``feature_angle``, the local thickness over ``gap_cells`` and the
    ``size_array`` values, never below ``min_size``. The tree is 2:1-balanced
    across faces, edges and corners, with ``grading_layers`` rings of leaves
    around the skin at its level. A point is inside where the fast generalised
    winding number exceeds ``winding_threshold``.

    :param mesh: the skin.
    :param output: ``"castellated"`` (hexahedral leaves with surrogate-boundary
        tags), ``"tet"`` (conforming tetrahedra fitted to the skin, the
        default) or ``"hex"`` (all-hexahedra).
    :param min_size: the smallest leaf edge; 0: ``max_size / 32``.
    :param max_size: the largest leaf edge; 0: an eighth of the skin's
        bounding-box diagonal.
    :param max_level: a cap on the finest level (at most 19); -1: from the
        sizes.
    :param bounds: ``(xmin, ymin, zmin, xmax, ymax, zmax)`` to mesh instead of
        the padded bounding box.
    :param padding: padding on every side, in world units.
    :param padding_relative: padding on every side, as a fraction of the
        diagonal.
    :param curvature_angle: the angle in degrees a leaf edge may subtend on
        the skin; 0 disables the curvature size.
    :param feature_angle: edges sharper than this dihedral (degrees) are
        features; 0 disables them.
    :param corner_angle: a feature vertex where the feature turns by more than
        this is a corner.
    :param feature_size: the size along feature edges; 0: ``min_size``.
    :param gap_cells: leaves across the local thickness; 0 disables it.
    :param gap_outside: also resolve gaps outside the skin.
    :param size_array: a ``point_data`` array of target sizes on the skin.
    :param surface_region: restrict the skin to this cell region.
    :param grading_layers: rings of leaves around every cut leaf at its level.
    :param skin_grading: ``"graded"`` (the size field also varies along the
        skin) or ``"uniform"`` (one level along the skin, so the stuffing
        angle bound applies).
    :param surface_size: the skin size under uniform grading, and the cell
        size of the hexahedral output; 0: the smallest the field asks for.
    :param winding_threshold: inside where the winding number exceeds this.
    :param winding_beta: accuracy of the fast winding number; 0 sums every
        triangle.
    :param conformity: ``"stuffing"`` or ``"cleaving"`` (tetrahedral output).
    :param alpha_long: stuffing's threshold on long lattice edges.
    :param alpha_short: stuffing's threshold on short lattice edges.
    :param preserve_features: snap the boundary to feature corners and edges.
    :param fit_iterations: rounds of optimisation (ODT relocation, sliding
        along the skin, flips); 0 disables them.
    :param dihedral_floor_guard: undo optimisation moves that lower the
        minimum dihedral or raise the maximum.
    :param fill: castellated output: keep ``"all"`` leaves, the ``"inside"``
        ones, or ``"cover"`` (inside or cut).
    :param keep_outside: tetrahedral output: also mesh the outside, as the
        ``outside`` region.
    :param attach_distance: attach ``sdf:distance`` (castellated: and
        ``sdf:closest_point``).
    :param attach_quality: attach the ``quality:*`` cell data.
    :param pillow: hexahedral output: a layer of hexahedra along the boundary.
    :param untangle_iterations: hexahedral output: untangling iterations.
    :param max_cells: refuse an octree with more leaves.
    :param max_tets: refuse an output with more cells.
    :param return_report: also return the report.
    :returns: the volume mesh, or ``(mesh, report)``. The report holds the
        frame (``origin``, ``root_size``, ``root_cells``, ``max_level``,
        ``finest_level``, ``surface_level``, ``min_size``, ``max_size``), the
        tree (``num_leaves``, ``leaves_per_level``), the counts
        (``num_points``, ``num_cells``, ``num_hexahedra``, ``num_tetrahedra``,
        ``num_inside``, ``num_boundary_cells``, ``num_surrogate_faces``,
        ``num_hanging``, ``num_cut_points``, ``num_warped``,
        ``num_warp_refused``), the quality (``min_dihedral``,
        ``max_dihedral``, ``min_scaled_jacobian``, ``num_inverted``,
        ``hausdorff_to_skin``, ``num_boundary_edges_open``,
        ``num_non_manifold_boundary_edges``, ``quality`` with its histograms),
        the features (``num_corners``, ``num_corners_snapped``,
        ``num_feature_chains``, ``num_feature_chains_missed``,
        ``num_feature_moves``), the optimisation (``num_optimized_moves``,
        ``num_flips``), ``guarantee_applies``, ``num_winding_rechecks``, ``timings`` (seconds
        per stage) and ``warnings``.
    :raises ValueError: on an option out of range.
    :raises NotImplementedError: without the compiled core.
    """
    opts = {
        "output": _choice("output", output, OUTPUTS),
        "min_size": _nonneg("min_size", min_size),
        "max_size": _nonneg("max_size", max_size),
        "max_level": int(max_level),
        "bounds": None if bounds is None else tuple(float(b) for b in bounds),
        "padding": _nonneg("padding", padding),
        "padding_relative": _nonneg("padding_relative", padding_relative),
        "curvature_angle": float(curvature_angle),
        "feature_angle": float(feature_angle),
        "corner_angle": float(corner_angle),
        "feature_size": _nonneg("feature_size", feature_size),
        "gap_cells": _nonneg("gap_cells", gap_cells),
        "gap_outside": bool(gap_outside),
        "size_array": "" if size_array is None else str(size_array),
        "surface_region": "" if surface_region is None else str(surface_region),
        "grading_layers": int(grading_layers),
        "skin_grading": _choice("skin_grading", skin_grading, SKIN_GRADINGS),
        "surface_size": _nonneg("surface_size", surface_size),
        "winding_threshold": float(winding_threshold),
        "winding_beta": _nonneg("winding_beta", winding_beta),
        "conformity": _choice("conformity", conformity, CONFORMITIES),
        "alpha_long": float(alpha_long),
        "alpha_short": float(alpha_short),
        "preserve_features": bool(preserve_features),
        "fit_iterations": int(fit_iterations),
        "dihedral_floor_guard": bool(dihedral_floor_guard),
        "fill": _choice("fill", fill, FILLS),
        "keep_outside": bool(keep_outside),
        "attach_distance": bool(attach_distance),
        "attach_quality": bool(attach_quality),
        "pillow": bool(pillow),
        "untangle_iterations": int(untangle_iterations),
        "max_cells": int(max_cells),
        "max_tets": int(max_tets),
    }
    if opts["bounds"] is not None and len(opts["bounds"]) != 6:
        raise ValueError(f"{_PREFIX}bounds must have six entries")
    if (
        opts["min_size"] > 0.0
        and opts["max_size"] > 0.0
        and opts["min_size"] > opts["max_size"]
    ):
        raise ValueError(f"{_PREFIX}min_size must not exceed max_size")
    if not -1 <= opts["max_level"] <= 19:
        raise ValueError(f"{_PREFIX}max_level must be -1 (automatic) or in [0, 19]")
    if not 0.0 <= opts["curvature_angle"] < 180.0:
        raise ValueError(f"{_PREFIX}curvature_angle must be in [0, 180) degrees")
    for name in ("feature_angle", "corner_angle"):
        if not 0.0 <= opts[name] <= 180.0:
            raise ValueError(f"{_PREFIX}{name} must be in [0, 180] degrees")
    if not 0 <= opts["grading_layers"] <= 8:
        raise ValueError(f"{_PREFIX}grading_layers must be in [0, 8]")
    if not 0.0 < opts["winding_threshold"] < 1.0:
        raise ValueError(f"{_PREFIX}winding_threshold must be in (0, 1)")
    for name in ("alpha_long", "alpha_short"):
        if not 0.0 < opts[name] < 0.5:
            raise ValueError(f"{_PREFIX}{name} must be in (0, 0.5)")
    for name in ("fit_iterations", "untangle_iterations"):
        if opts[name] < 0:
            raise ValueError(f"{_PREFIX}{name} must be >= 0")
    for name in ("max_cells", "max_tets"):
        if opts[name] < 1:
            raise ValueError(f"{_PREFIX}{name} must be >= 1")

    try:
        from . import _core
    except ImportError as exc:
        raise NotImplementedError(
            f"{_PREFIX}needs the compiled C++ core (meshioplusplus._core); "
            "it has no numpy twin"
        ) from exc
    res = _core.octree_mesh(mesh, **opts)
    out = res["mesh"]
    if not return_report:
        return out
    report = {"output": str(res["output"])}
    report["origin"] = [float(v) for v in res["origin"]]
    report["root_cells"] = [int(v) for v in res["root_cells"]]
    report.update({k: int(res[k]) for k in _REPORT_INTS})
    report.update({k: float(res[k]) for k in _REPORT_FLOATS})
    report["leaves_per_level"] = [int(v) for v in res["leaves_per_level"]]
    report["guarantee_applies"] = bool(res["guarantee_applies"])
    report["quality"] = res["quality"]
    report["timings"] = {str(k): float(v) for k, v in res["timings"]}
    report["warnings"] = [str(w) for w in res["warnings"]]
    return out, report


__all__ = ["octree_mesh"]
