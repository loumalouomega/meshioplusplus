"""Refine a segmentation until each patch is explained by one surface.

Between :func:`segment` and :func:`fit_primitives` and :func:`fit_cad` in the
skin-to-CAD chain, in the core (``operations/refine_patches.hpp``): a hybrid
variational shape approximation (Cohen-Steiner, Alliez and Desbrun 2004, with
the primitive proxies of Wu and Kobbelt 2005) moves patch boundaries, splits a
primitive out of a patch no primitive explains by seeded RANSAC (Schnabel,
Wahl and Klein 2007) and merges neighbours one primitive explains. There is no
numpy twin: without the compiled core the function raises
``NotImplementedError``.
"""

from __future__ import annotations

import math

from ._fit_cad import DEFAULT_RANSAC_SEED, _check_seed, _types_mask

_PREFIX = "meshio++: refine_patches: "

#: Why refinement stopped, as the report's ``stop`` names it.
STOPS = ("unchanged", "small_gain", "iteration_cap", "energy_rose")

_REPORT_INTS = (
    "num_patches",
    "num_iterations",
    "num_splits",
    "num_merges",
    "num_absorbed",
    "num_primitive",
    "num_freeform",
    "num_facets",
)
_REPORT_FLOATS = (
    "tolerance",
    "explained_area_fraction",
    "initial_explained_area_fraction",
    "hausdorff_max",
    "hausdorff_rms",
    "initial_energy",
    "energy",
)


def refine_patches(
    mesh,
    *,
    tolerance: float = 1e-3,
    relative: bool = True,
    normal_weight: float = 1.0,
    epsilon: float = 1e-3,
    feature_angle: float = 30.0,
    max_iterations: int = 10,
    min_patch_cells: int = 3,
    ransac_seed: int = DEFAULT_RANSAC_SEED,
    types=None,
    allow_split: bool = True,
    allow_merge: bool = True,
    prefix: str = "patch",
    return_report: bool = False,
):
    """Refine the ``segment:patch`` labels of a surface by variational shape
    approximation.

    :func:`segment` draws its boundaries once: a fillet strip claimed by the
    plane beside it, a cylinder cut in two at a smooth edge, or a patch with a
    cylinder inside it stays that way. ``refine_patches`` explains each patch
    by a proxy -- a plane, sphere, cylinder, cone or torus that holds every
    point within ``tolerance``, else a freeform quadric height field -- and
    iterates: the patches grow back from their best-explained cells, cheapest
    first (a Lloyd step, undone if it raises the energy); a freeform patch is
    searched for a primitive inside it by seeded RANSAC (``allow_split``);
    neighbours one primitive explains together are merged (``allow_merge``);
    a freeform patch smaller than ``min_patch_cells`` joins a neighbour.
    Growth never crosses a feature edge, so patches stay edge-connected. See
    ``doc/refine_patches.md``.

    :param mesh: a surface mesh carrying ``segment:patch`` (from
        :func:`segment`).
    :param tolerance: the distance every point of a patch must lie within from
        its primitive; a fraction of the surface's bounding-box diagonal when
        ``relative``.
    :param relative: ``tolerance`` is relative to the diagonal (default)
        rather than a length.
    :param normal_weight: the weight of the normal term in a cell's cost,
        ``>= 0``.
    :param epsilon: stop when an iteration adds less than this to the share of
        the area primitives explain, ``>= 0``.
    :param feature_angle: the dihedral angle in degrees above which an edge is
        a feature edge no patch crosses, in ``[0, 180]`` (``segment``'s
        default).
    :param max_iterations: the most iterations, ``>= 0`` (``0`` rewrites the
        outputs of the input partition).
    :param min_patch_cells: a freeform or facets patch with fewer cells joins
        a neighbour, ``>= 1``.
    :param ransac_seed: the seed of the RANSAC split, in ``[0, 2**64)``; the
        output is byte-identical for a seed on every backend and thread count.
    :param types: the primitives a proxy may be, from ``plane``, ``sphere``,
        ``cylinder``, ``cone``, ``torus`` (default all); an empty sequence
        leaves every patch freeform.
    :param allow_split: split a primitive out of a freeform patch.
    :param allow_merge: merge neighbouring patches one primitive explains.
    :param prefix: the refined patches' cell regions are named
        ``<prefix>:<k>``.
    :param return_report: also return the report.
    :returns: the input -- points, connectivity, block order, data and
        regions untouched -- with ``segment:patch``, ``segment:adjacency`` and
        the patch regions rewritten for the refined patches, the ``fit:*`` data
        of :func:`fit_primitives` at the absolute tolerance (so
        :func:`fit_cad` reads it unchanged), and per-patch field data
        ``refine:patch_id``, ``refine:patch_proxy`` (0 facets, 1-5 as
        ``fit:type``, 6 freeform), ``refine:patch_iterations``,
        ``refine:patch_rms``, ``refine:patch_max``, ``refine:patch_cells``,
        ``refine:patch_demotion`` (0 none, 1 too small, 2 no primitive within
        the tolerance), plus ``refine:explained_area_fraction`` and
        ``refine:hausdorff`` (the largest and the RMS distance of a point from
        its proxy). With ``return_report``, ``(mesh, report)``, the report
        holding ``num_patches``, ``num_iterations`` (Lloyd steps kept),
        ``num_splits``, ``num_merges``, ``num_absorbed``, ``num_primitive``,
        ``num_freeform``, ``num_facets``, the absolute ``tolerance``,
        ``explained_area_fraction`` and ``initial_explained_area_fraction``,
        ``hausdorff_max`` and ``hausdorff_rms``, ``initial_energy`` and
        ``energy``, and ``stop`` (``"unchanged"``, ``"small_gain"``,
        ``"iteration_cap"`` or ``"energy_rose"``).
    :raises ValueError: on an option out of range, or without
        ``segment:patch``.
    :raises NotImplementedError: without the compiled core.
    """
    opts = {
        "tolerance": float(tolerance),
        "relative": bool(relative),
        "normal_weight": float(normal_weight),
        "epsilon": float(epsilon),
        "feature_angle": float(feature_angle),
        "max_iterations": int(max_iterations),
        "min_patch_cells": int(min_patch_cells),
        "ransac_seed": _check_seed(ransac_seed, _PREFIX),
        "types": _types_mask(types, _PREFIX),
        "allow_split": bool(allow_split),
        "allow_merge": bool(allow_merge),
        "prefix": str(prefix),
    }
    if not (opts["tolerance"] > 0.0 and math.isfinite(opts["tolerance"])):
        raise ValueError(f"{_PREFIX}the tolerance must be finite and > 0")
    if not (opts["normal_weight"] >= 0.0 and math.isfinite(opts["normal_weight"])):
        raise ValueError(f"{_PREFIX}normal_weight must be finite and >= 0")
    if not (opts["epsilon"] >= 0.0 and math.isfinite(opts["epsilon"])):
        raise ValueError(f"{_PREFIX}epsilon must be finite and >= 0")
    if not 0.0 <= opts["feature_angle"] <= 180.0:
        raise ValueError(f"{_PREFIX}the feature angle must lie in [0, 180] degrees")
    if opts["max_iterations"] < 0:
        raise ValueError(f"{_PREFIX}max_iterations must be >= 0")
    if opts["min_patch_cells"] < 1:
        raise ValueError(f"{_PREFIX}min_patch_cells must be >= 1")
    if not opts["prefix"]:
        raise ValueError(f"{_PREFIX}the region prefix must not be empty")
    if "segment:patch" not in getattr(mesh, "cell_data", {}):
        raise ValueError(
            f"{_PREFIX}the mesh has no 'segment:patch' cell data; run segment first"
        )

    try:
        from . import _core
    except ImportError as exc:
        raise NotImplementedError(
            f"{_PREFIX}needs the compiled C++ core (meshioplusplus._core); "
            "it has no numpy twin"
        ) from exc
    res = _core.refine_patches(
        mesh,
        opts["tolerance"],
        opts["relative"],
        opts["normal_weight"],
        opts["epsilon"],
        opts["feature_angle"],
        opts["max_iterations"],
        opts["min_patch_cells"],
        opts["ransac_seed"],
        opts["types"],
        opts["allow_split"],
        opts["allow_merge"],
        opts["prefix"],
    )
    out = res["mesh"]
    if not return_report:
        return out
    report = {k: int(res[k]) for k in _REPORT_INTS}
    report.update({k: float(res[k]) for k in _REPORT_FLOATS})
    report["stop"] = str(res["stop"])
    return out, report


__all__ = ["refine_patches"]
