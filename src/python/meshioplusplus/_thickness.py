"""Wall thickness and clearance per facet, by a ray along its normal (roadmap 3.1.4).

From each 2-D cell's area-weighted centroid a ray runs along its unit normal
through a bounding-volume hierarchy of the target surface. Inward (the
default for one surface) it stops at the first wall it leaves through, the far
side of the wall: the thickness. Outward (the default against a second
surface) it stops at the first wall that faces it: the clearance. The surface
must be consistently oriented with outward normals. See ``doc/thickness.md``.

The C++ core does the work; there is no numpy twin.
"""

from __future__ import annotations

import math

__all__ = ["thickness"]

_PREFIX = "meshio++: thickness: "
_DIRECTIONS = ("auto", "inward", "outward")
_REPORT_INTS = ("num_cells", "num_hit", "num_missed", "num_below", "min_cell")
_REPORT_FLOATS = ("min", "max", "mean", "max_distance", "p2", "p95")
_METHODS = ("ray", "sphere", "cone")


def _non_negative(value, name: str) -> float:
    v = float(value)
    if not math.isfinite(v) or v < 0.0:
        raise ValueError(
            f"{_PREFIX}{name} must be finite and non-negative, got {value!r}"
        )
    return v


def _method(value) -> str:
    m = str(value).lower()
    if m not in _METHODS:
        raise ValueError(f"{_PREFIX}method must be one of {_METHODS}, got {value!r}")
    return m


def thickness(
    mesh,
    other=None,
    *,
    direction: str | None = None,
    max_distance: float = 0.0,
    threshold: float | None = None,
    region: str | None = None,
    region_other: str | None = None,
    attach: bool = True,
    method: str = "ray",
    cone_angle: float = 15.0,
    cone_rays: int = 16,
    two_sided: bool = False,
    max_normal_deviation: float = 90.0,
    return_report: bool = False,
):
    """Measure the wall thickness or clearance of each facet of ``mesh``.

    :param mesh: a consistently oriented surface mesh (triangles, quads,
        polygons; other cells are not measured).
    :param other: the surface to reach; ``None`` measures ``mesh`` against
        itself.
    :param direction: ``"inward"`` (against the normal, to the far side of the
        wall), ``"outward"`` (along it, across a gap) or ``None``/``"auto"``:
        inward for one surface, outward for two.
    :param max_distance: the longest ray; 0 takes the diagonal of both
        surfaces' bounding box.
    :param threshold: count the cells whose distance is below it.
    :param region: measure only this named ``Cell`` region of ``mesh``.
    :param region_other: reach only this named ``Cell`` region of the target
        (``other``, else ``mesh``).
    :param attach: attach the cell data ``thickness:distance`` (Float64, NaN
        for a miss) and ``thickness:hit_cell`` (Int64, -1 for a miss).
    :param method: ``"ray"`` (one ray along the normal), ``"sphere"`` (twice
        the distance from the ray's midpoint to the nearest surface outside the
        facet's one-ring: the perpendicular distance where the walls diverge)
        or ``"cone"`` (the shortest of ``cone_rays`` rays within
        ``cone_angle`` degrees of the normal: a chamfer's narrow spot).
    :param cone_angle: the cone's half-angle in degrees, in (0, 90).
    :param cone_rays: the cone's rays, the central one included.
    :param two_sided: rays both ways along the normal, either face of every
        triangle counting (a scan's winding is unreliable); with
        ``max_distance`` 0 the ray is then a fifth of the diagonal.
    :param max_normal_deviation: step over triangles more edge-on to the ray
        than this, in degrees in (0, 90]; 90 accepts any facing the right way.
    :param return_report: also return the report.
    :returns: ``mesh``, bit-identical but for the two arrays; with
        ``return_report`` ``(mesh, report)``, where ``report`` has
        ``direction``, ``method``, ``num_cells``, ``num_hit``, ``num_missed``,
        ``num_below``, ``min``, ``max``, ``mean``, ``p2``, ``p95`` (the 2nd and
        95th percentiles of the hits), ``min_cell``, ``min_point``,
        ``max_distance`` and ``warnings``.
    :raises ValueError: for a volume, polyhedron or higher-order block, an
        unknown region, an unknown direction or a negative distance.
    """
    try:
        from . import _core
    except ImportError as exc:  # pragma: no cover - the numpy-only install
        raise NotImplementedError(
            "meshio++: thickness: needs the compiled C++ core "
            "(meshioplusplus._core); it has no numpy twin"
        ) from exc
    direction = "auto" if direction is None else str(direction).lower()
    if direction not in _DIRECTIONS:
        raise ValueError(
            f"{_PREFIX}direction must be one of {_DIRECTIONS}, got {direction!r}"
        )
    res = _core.thickness(
        mesh,
        other,
        direction=direction,
        max_distance=_non_negative(max_distance, "max_distance"),
        threshold=0.0 if threshold is None else _non_negative(threshold, "threshold"),
        region=region or "",
        region_target=region_other or "",
        attach=bool(attach),
        method=_method(method),
        cone_angle=float(cone_angle),
        cone_rays=int(cone_rays),
        two_sided=bool(two_sided),
        max_normal_deviation=float(max_normal_deviation),
    )
    out = res["mesh"]
    if not return_report:
        return out
    report = {"direction": str(res["direction"]), "method": str(res["method"])}
    report.update({k: int(res[k]) for k in _REPORT_INTS})
    report.update({k: float(res[k]) for k in _REPORT_FLOATS})
    report["min_point"] = tuple(float(c) for c in res["min_point"])
    report["warnings"] = [str(w) for w in res["warnings"]]
    return out, report
