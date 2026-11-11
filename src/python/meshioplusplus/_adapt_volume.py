"""Size-field tetrahedral remeshing (``meshioplusplus.adapt_volume``).

A mesh *operation*: turn a valid, over-refined tetrahedral mesh (typically
:func:`octree_mesh`'s graded output) into one with fewer, better-shaped
tetrahedra, towards an isotropic size field with a gradation and a bound on how
far the boundary may move. Each iteration collapses the edges shorter than 4/5
in the metric, splits those longer than 4/3, swaps (2-3, 3-2, edge removal,
multi-face removal) and smooths (ODT or AMIPS inside; tangential moves projected
back onto the input boundary on it). Corners stay, ridge vertices move along
their chain, surface vertices within their patch; material interfaces are kept
as internal boundaries. Every edit is guarded by exact orientation tests and a
quality floor, and the boundary stays within ``hausd`` of the input's.

The algorithm is from published descriptions only (Dapogny, Dobrzynski and Frey
2014; Klingner and Shewchuk 2007; Freitag and Ollivier-Gooch 1997; Hu et al.
2020; Jiang et al. 2022; Botsch and Kobbelt 2004; Alliez et al. 2005); see
``doc/adapt_volume.md``.

**C++-core only, with no numpy fallback**: every edit is a discrete branch on
exact predicates and near-tie comparisons, so a second implementation could
diverge into a different mesh. Without ``_core`` this raises
:class:`NotImplementedError`.
"""

from __future__ import annotations

__all__ = ["adapt_volume"]

_PREFIX = "meshio++: adapt_volume: "

_REPORT_INTS = (
    "num_tets_before",
    "num_tets_after",
    "num_points_before",
    "num_points_after",
    "num_collapsed",
    "num_split",
    "num_swapped",
    "num_smoothed",
    "num_rejected_by_link",
    "num_rejected_by_inversion",
    "num_rejected_by_quality",
    "num_rejected_by_length",
    "num_rejected_by_class",
    "num_rejected_by_hausdorff",
    "num_rejected_by_normal",
    "num_corners",
    "num_ridge_vertices",
    "num_chains",
    "num_patches",
    "iterations",
)
_REPORT_FLOATS = (
    "min_dihedral_before",
    "min_dihedral_after",
    "max_dihedral_before",
    "max_dihedral_after",
    "mean_quality_before",
    "mean_quality_after",
    "max_boundary_deviation",
    "hmin",
    "hmax",
    "hausd",
    "volume_before",
    "volume_after",
)


def adapt_volume_report(res):
    """The plain-Python report of a ``_core`` adapt_volume result dict."""
    report = {k: int(res[k]) for k in _REPORT_INTS}
    report.update({k: float(res[k]) for k in _REPORT_FLOATS})
    report["region_volume_change"] = {
        str(label): float(dv) for label, dv in res["region_volume_change"]
    }
    report["timings"] = {str(k): float(v) for k, v in res["timings"]}
    report["warnings"] = [str(w) for w in res["warnings"]]
    return report


def adapt_volume_options(
    hmin=0.0,
    hmax=0.0,
    hgrad=1.3,
    hausd=0.0,
    hausd_relative=False,
    size_array=None,
    curvature_sizing=False,
    feature_angle=30.0,
    corner_angle=30.0,
    max_normal_deviation=30.0,
    preserve_boundary=False,
    preserve_features=True,
    quality_metric="min_dihedral",
    quality_min=-1.0,
    max_iterations=10,
    target_min_dihedral=0.0,
    frozen=None,
    material_arrays=None,
    collapse=True,
    split=True,
    swap=True,
    smooth=True,
    boundary_sizing="input",
):
    """Validate and normalise :func:`adapt_volume`'s options into the keyword
    dictionary ``_core`` reads (shared with ``octree_mesh(post_adapt=True)``)."""
    if not float(hgrad) >= 1.0:
        raise ValueError(f"{_PREFIX}hgrad must be >= 1")
    if int(max_iterations) < 0:
        raise ValueError(f"{_PREFIX}max_iterations must be >= 0")
    if not 0.0 < float(max_normal_deviation) <= 180.0:
        raise ValueError(f"{_PREFIX}max_normal_deviation must be in (0, 180] degrees")
    if str(quality_metric) not in ("min_dihedral", "mean_ratio", "amips"):
        raise ValueError(
            f"{_PREFIX}quality_metric must be 'min_dihedral', 'mean_ratio' or 'amips'"
        )
    if str(boundary_sizing) not in ("input", "geometry"):
        raise ValueError(f"{_PREFIX}boundary_sizing must be 'input' or 'geometry'")
    if float(hmin) > 0.0 and float(hmax) > 0.0 and float(hmin) > float(hmax):
        raise ValueError(f"{_PREFIX}hmin must not exceed hmax")
    return {
        "hmin": float(hmin),
        "hmax": float(hmax),
        "hgrad": float(hgrad),
        "hausd": float(hausd),
        "hausd_relative": bool(hausd_relative),
        "size_array": None if size_array is None else str(size_array),
        "curvature_sizing": bool(curvature_sizing),
        "feature_angle": float(feature_angle),
        "corner_angle": float(corner_angle),
        "max_normal_deviation": float(max_normal_deviation),
        "preserve_boundary": bool(preserve_boundary),
        "preserve_features": bool(preserve_features),
        "quality_metric": str(quality_metric),
        "quality_min": float(quality_min),
        "max_iterations": int(max_iterations),
        "target_min_dihedral": float(target_min_dihedral),
        "frozen": None if frozen is None else [1 if v else 0 for v in frozen],
        "material_arrays": (
            None if material_arrays is None else [str(n) for n in material_arrays]
        ),
        "collapse": bool(collapse),
        "split": bool(split),
        "swap": bool(swap),
        "smooth": bool(smooth),
        "boundary_sizing": str(boundary_sizing),
    }


def adapt_volume(
    mesh,
    *,
    hmin=0.0,
    hmax=0.0,
    hgrad=1.3,
    hausd=0.0,
    hausd_relative=False,
    size_array=None,
    curvature_sizing=False,
    feature_angle=30.0,
    corner_angle=30.0,
    max_normal_deviation=30.0,
    preserve_boundary=False,
    preserve_features=True,
    quality_metric="min_dihedral",
    quality_min=-1.0,
    max_iterations=10,
    target_min_dihedral=0.0,
    frozen=None,
    material_arrays=None,
    collapse=True,
    split=True,
    swap=True,
    smooth=True,
    boundary_sizing="input",
    return_report=False,
):
    """Remesh a tetrahedral mesh towards a graded isotropic size field.

    :param mesh: a tetra-only volume mesh (one or more ``tetra`` blocks). Never
        modified.
    :param hmin: the smallest target edge length; ``<= 0`` takes the smallest
        boundary size. It also raises the boundary sizes, which otherwise keep
        the input boundary's own resolution.
    :param hmax: the largest target edge length; ``<= 0`` takes an eighth of
        the bounding-box diagonal.
    :param hgrad: the gradation: neighbouring sizes differ by at most
        ``(hgrad - 1) |uv|``.
    :param hausd: how far the boundary may move from the input boundary;
        ``<= 0`` takes 1e-3 of the bounding-box diagonal.
    :param hausd_relative: read ``hausd`` as a fraction of the diagonal.
    :param size_array: a one-component float ``point_data`` name of target
        sizes, or ``None``.
    :param curvature_sizing: bound the boundary size by ``sqrt(8 hausd R)``.
    :param feature_angle: boundary edges whose faces turn by more than this
        many degrees are ridges (with ``preserve_features``).
    :param corner_angle: a ridge vertex where its chain turns by more than this
        is a corner.
    :param max_normal_deviation: the most a boundary face may turn in one edit.
    :param preserve_boundary: freeze the boundary bit for bit.
    :param preserve_features: keep sharp edges and corners.
    :param quality_metric: ``"min_dihedral"``, ``"mean_ratio"`` or
        ``"amips"`` (AMIPS smoothing inside, the mean ratio guarded).
    :param quality_min: an edit may lower its star's least quality only down to
        this; ``< 0`` takes 10 degrees, or 0.2.
    :param max_iterations: collapse / split / swap / smooth cycles.
    :param target_min_dihedral: stop once the least dihedral angle reaches it
        (``0`` disables).
    :param frozen: a per-point mask of vertices that never move or go.
    :param material_arrays: integer ``cell_data`` names kept by parent; ``None``
        takes the known material names.
    :param collapse: run the collapse pass (and ``split``, ``swap``,
        ``smooth`` theirs), for debugging and ablation.
    :param boundary_sizing: where the boundary sizes come from: ``"input"``
        (the default) keeps the input boundary's own resolution (the mean
        length of each boundary vertex's edges), ``"geometry"`` sizes it by
        curvature, ``sqrt(8 hausd R)``, and ``hmax`` where it is flat, as MMG
        does, so a planar or gently curved boundary coarsens too.
    :param return_report: also return the report: counts before and after,
        the edits, the refusals by guard (``num_rejected_by_*``), the dihedral
        range, the mean quality, ``max_boundary_deviation`` (the symmetric
        Hausdorff distance to the input boundary), ``hmin``/``hmax``/``hausd``
        used, the volumes, ``region_volume_change`` per material class, the
        classification counts, ``iterations``, ``timings`` and ``warnings``.
    :returns: the adapted mesh, or ``(mesh, report)``.
    :raises ValueError: on a non-``tetra`` or ragged block, no tetrahedra, a
        tangled or non-manifold volume, or a bad option.
    :raises NotImplementedError: without the compiled core.
    """
    opts = adapt_volume_options(
        hmin=hmin,
        hmax=hmax,
        hgrad=hgrad,
        hausd=hausd,
        hausd_relative=hausd_relative,
        size_array=size_array,
        curvature_sizing=curvature_sizing,
        feature_angle=feature_angle,
        corner_angle=corner_angle,
        max_normal_deviation=max_normal_deviation,
        preserve_boundary=preserve_boundary,
        preserve_features=preserve_features,
        quality_metric=quality_metric,
        quality_min=quality_min,
        max_iterations=max_iterations,
        target_min_dihedral=target_min_dihedral,
        frozen=frozen,
        material_arrays=material_arrays,
        collapse=collapse,
        split=split,
        swap=swap,
        smooth=smooth,
        boundary_sizing=boundary_sizing,
    )
    try:
        from . import _core
    except ImportError as exc:
        raise NotImplementedError(
            f"{_PREFIX}this operation is C++-core only and has no numpy "
            "fallback (its edits are discrete branches on exact predicates, so a "
            "second implementation could diverge) -- install a build with the "
            "compiled meshioplusplus._core extension"
        ) from exc
    res = _core.adapt_volume(mesh, **opts)
    out = res["mesh"]
    if not return_report:
        return out
    return out, adapt_volume_report(res)
