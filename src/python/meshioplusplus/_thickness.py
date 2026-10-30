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
_REPORT_FLOATS = ("min", "max", "mean", "max_distance")


def _non_negative(value, name: str) -> float:
    v = float(value)
    if not math.isfinite(v) or v < 0.0:
        raise ValueError(
            f"{_PREFIX}{name} must be finite and non-negative, got {value!r}"
        )
    return v


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
    :param return_report: also return the report.
    :returns: ``mesh``, bit-identical but for the two arrays; with
        ``return_report`` ``(mesh, report)``, where ``report`` has
        ``direction``, ``num_cells``, ``num_hit``, ``num_missed``,
        ``num_below``, ``min``, ``max``, ``mean``, ``min_cell``,
        ``min_point``, ``max_distance`` and ``warnings``.
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
    )
    out = res["mesh"]
    if not return_report:
        return out
    report = {"direction": str(res["direction"])}
    report.update({k: int(res[k]) for k in _REPORT_INTS})
    report.update({k: float(res[k]) for k in _REPORT_FLOATS})
    report["min_point"] = tuple(float(c) for c in res["min_point"])
    report["warnings"] = [str(w) for w in res["warnings"]]
    return out, report
