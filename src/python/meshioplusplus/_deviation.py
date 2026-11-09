"""Deviation analysis of a scan against its intended geometry (roadmap 7.4).

Every scan point's signed distance to a reference surface (positive outside,
excess material; negative inside, missing material), or to a fitted primitive,
and the statistics an inspection reads from it: min, max, mean, RMS, sigma,
percentiles, a histogram, the share within a tolerance and a suggested colour
range. See ``doc/deviation.md``.

The C++ core does the work; there is no numpy twin.
"""

from __future__ import annotations

import math

import numpy as np

__all__ = ["deviation", "deviation_regions", "exaggerate"]

_PREFIX = "meshio++: deviation: "
_TYPES = ("none", "plane", "sphere", "cylinder", "cone", "torus")
_REPORT_INTS = ("num_points", "num_measured", "num_outside_range", "num_without_normal")
_REPORT_FLOATS = (
    "min",
    "max",
    "mean",
    "rms",
    "sigma",
    "share_within",
    "suggested_range",
)


def _core_module():
    try:
        from . import _core
    except ImportError as exc:  # pragma: no cover - the numpy-only install
        raise NotImplementedError(
            f"{_PREFIX}needs the compiled C++ core (meshioplusplus._core); "
            "it has no numpy twin"
        ) from exc
    return _core


def primitive_reference(reference, patch=None):
    """``(type, params)`` of a primitive reference.

    ``reference`` is a ``(type, params)`` pair (``type`` a name or 1-5, ten
    ``fit:patch_params`` entries), a dict with ``type`` and ``params``, or a
    mesh ``fit_primitives`` wrote, whose row for ``patch`` (a patch label, or
    the only row when ``None``) is taken.
    """
    if isinstance(reference, dict):
        t, p = reference["type"], reference["params"]
    elif isinstance(reference, (tuple, list)) and len(reference) == 2:
        t, p = reference
    else:
        fd = reference.field_data
        if "fit:patch_params" not in fd:
            raise ValueError(
                f"{_PREFIX}the reference mesh carries no fit:patch_params "
                "(run fit_primitives first)"
            )
        ids = np.asarray(fd["fit:patch_id"]).ravel()
        if patch is None:
            if len(ids) != 1:
                raise ValueError(
                    f"{_PREFIX}the fitted mesh has {len(ids)} patches; name one with patch="
                )
            row = 0
        else:
            hits = np.nonzero(ids == int(patch))[0]
            if len(hits) == 0:
                raise ValueError(f"{_PREFIX}no fitted patch {patch!r}")
            row = int(hits[0])
        t = int(np.asarray(fd["fit:patch_type"]).ravel()[row])
        p = np.asarray(fd["fit:patch_params"], dtype=np.float64).reshape(-1, 10)[row]
    if isinstance(t, str):
        if t not in _TYPES[1:]:
            raise ValueError(f"{_PREFIX}unknown primitive type {t!r}")
        t = _TYPES.index(t)
    t = int(t)
    if not 1 <= t <= 5:
        raise ValueError(f"{_PREFIX}the primitive has no fitted type ({t})")
    p = [float(x) for x in np.asarray(p, dtype=np.float64).ravel()]
    if len(p) != 10:
        raise ValueError(f"{_PREFIX}a primitive takes ten parameters, got {len(p)}")
    return t, p


def _is_mesh(obj):
    return hasattr(obj, "points") and hasattr(obj, "cells")


def _report(res):
    report = {k: int(res[k]) for k in _REPORT_INTS}
    report.update({k: float(res[k]) for k in _REPORT_FLOATS})
    report["percentiles"] = {str(k): float(v) for k, v in res["percentiles"].items()}
    report["histogram_edges"] = [float(x) for x in res["histogram_edges"]]
    report["histogram_counts"] = [int(x) for x in res["histogram_counts"]]
    report["warnings"] = [str(w) for w in res["warnings"]]
    return report


def deviation(
    scan,
    reference,
    *,
    patch=None,
    tolerance: float = 0.0,
    bins: int = 20,
    range=None,  # noqa: A002 - the documented keyword
    max_normal_angle: float = 90.0,
    reference_region: str | None = None,
    attach: bool = True,
    return_report: bool = False,
):
    """Signed deviation of every point of ``scan`` from ``reference``.

    :param scan: a point cloud or a surface; its points are measured.
    :param reference: a surface mesh, or a fitted primitive: a
        ``(type, params)`` pair, a ``{"type", "params"}`` dict, or, with
        ``patch``, a mesh ``fit_primitives`` wrote.
    :param patch: measure against this patch label's fitted primitive of the
        ``reference`` mesh (``fit:patch_params``) instead of its surface.
    :param tolerance: report the share within ``+-tolerance`` (when positive).
    :param bins: histogram bins.
    :param range: ``(lo, hi)`` of the histogram; ``None`` takes ``+-max|d|``.
    :param max_normal_angle: the far-side rule's angle in degrees: a reference
        triangle whose normal turns away from the scan point's by more is
        skipped (180 turns it off). Not used against a primitive.
    :param reference_region: restrict the reference surface to a ``Cell`` region.
    :param attach: attach the Float64 point data ``deviation:distance``.
    :param return_report: also return the report.
    :returns: ``scan`` with ``deviation:distance``; with ``return_report``
        ``(mesh, report)``.
    """
    _core = _core_module()
    lo, hi = (0.0, 0.0) if range is None else (float(range[0]), float(range[1]))
    if range is not None and not lo < hi:
        raise ValueError(f"{_PREFIX}range must be (lo, hi) with lo < hi")
    tol = float(tolerance)
    if not (math.isfinite(tol) and tol >= 0.0):
        raise ValueError(f"{_PREFIX}tolerance must be finite and >= 0")
    kwargs = dict(
        tolerance=tol,
        bins=int(bins),
        range_min=lo,
        range_max=hi,
        max_normal_angle=float(max_normal_angle),
        reference_region=reference_region or "",
        attach=bool(attach),
    )
    if _is_mesh(reference) and patch is None:
        res = _core.deviation(scan, reference, **kwargs)
    else:
        t, p = primitive_reference(reference, patch)
        res = _core.deviation(
            scan, None, primitive_type=t, primitive_params=p, **kwargs
        )
    out = res["mesh"]
    return (out, _report(res)) if return_report else out


def deviation_regions(
    mesh,
    tolerance: float,
    *,
    field: str = "deviation:distance",
    prefix: str = "deviation:",
    min_area: float = 0.0,
    return_report: bool = False,
):
    """Connected regions of a measured surface beyond ``+-tolerance``.

    Cells joined across shared edges whose deviation (the mean of their
    corners' ``field`` values) exceeds ``+tolerance`` are excess material, below
    ``-tolerance`` missing material (roadmap 7.4.3).

    :param mesh: a surface carrying the point data ``field`` (``deviation`` writes
        ``deviation:distance``).
    :param tolerance: the threshold, positive.
    :param field: the signed deviation's point data.
    :param prefix: region names are ``<prefix>excess_<k>`` and
        ``<prefix>missing_<k>``.
    :param min_area: drop regions smaller than this area.
    :param return_report: also return the report.
    :returns: ``mesh`` with one named cell region per record and the cell data
        ``deviation:region`` (-1 within tolerance); with ``return_report``
        ``(mesh, report)``, ``report["regions"]`` a list of dicts (``name``,
        ``sign``, ``num_cells``, ``area``, ``extreme``, ``mean``, ``volume``,
        ``centroid``, ``bounds``), largest ``|volume|`` first.
    """
    _core = _core_module()
    tol = float(tolerance)
    if not (math.isfinite(tol) and tol > 0.0):
        raise ValueError(
            f"{_PREFIX}deviation_regions: tolerance must be finite and > 0"
        )
    res = _core.deviation_regions(
        mesh,
        tolerance=tol,
        field=str(field),
        prefix=str(prefix),
        min_area=float(min_area),
    )
    out = res["mesh"]
    if not return_report:
        return out
    report = {
        "regions": [
            {
                "name": str(r["name"]),
                "sign": int(r["sign"]),
                "num_cells": int(r["num_cells"]),
                "area": float(r["area"]),
                "extreme": float(r["extreme"]),
                "mean": float(r["mean"]),
                "volume": float(r["volume"]),
                "centroid": [float(x) for x in r["centroid"]],
                "bounds": [float(x) for x in r["bounds"]],
            }
            for r in res["regions"]
        ],
        "num_excess": int(res["num_excess"]),
        "num_missing": int(res["num_missing"]),
        "num_dropped": int(res["num_dropped"]),
        "warnings": [str(w) for w in res["warnings"]],
    }
    return out, report


def exaggerate(
    mesh,
    *,
    field: str = "deviation:distance",
    scale: float = 10.0,
    clamp: float = 0.0,
    smooth_fraction: float = 0.01,
    return_report: bool = False,
):
    """Displace each point along its normal by ``scale`` times its deviation.

    A mesh for a figure or an animation (roadmap 7.4.4): the deviation is
    clamped to ``+-clamp`` (0: the colour range ``deviation`` suggests, the
    95th percentile of ``|d|`` rounded up to two significant figures), then
    averaged over the points within ``smooth_fraction`` of the bounding-box
    diagonal so scan noise does not become fur. The measured scan is left
    untouched; connectivity and data are copied.

    :returns: the exaggerated mesh; with ``return_report`` ``(mesh, report)``
        with ``clamp``, ``radius``, ``num_moved`` and ``warnings``.
    """
    _core = _core_module()
    res = _core.exaggerate(
        mesh,
        field=str(field),
        scale=float(scale),
        clamp=float(clamp),
        smooth_fraction=float(smooth_fraction),
    )
    out = res["mesh"]
    if not return_report:
        return out
    return out, {
        "clamp": float(res["clamp"]),
        "radius": float(res["radius"]),
        "num_moved": int(res["num_moved"]),
        "warnings": [str(w) for w in res["warnings"]],
    }
