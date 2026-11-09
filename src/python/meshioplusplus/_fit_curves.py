"""Lines, arcs and circles from chained sections (roadmap 7.6.2).

The sketch a drawing, a DXF file or a CAD profile is made of: each chain of a
chained line mesh (``slice(chain=True)`` or ``sections``) is split into the
runs a line or a circular arc explains within a tolerance, in its section
plane's 2-D frame. C++-core only (see ``doc/fit_curves.md``).
"""

from __future__ import annotations

import json

import numpy as np

__all__ = ["fit_curves", "sketch_json"]

_PREFIX = "meshio++: fit_curves: "


def fit_curves(
    sections,
    *,
    tolerance: float = 0.0,
    snap_angle: float = 0.0,
    normal=None,
    min_arc_points: int = 4,
    min_arc_sweep: float = 5.0,
    return_report: bool = False,
):
    """Fit lines, arcs and circles to the chains of a chained line mesh.

    :param sections: a line mesh with ``section:chain`` cell data
        (``slice(chain=True)`` or ``sections``).
    :param tolerance: the largest residual a piece may have; 0 takes 1e-3 of
        the bounding-box diagonal.
    :param snap_angle: snap lines within this many degrees of the frame's axes
        or their diagonals (0: no snapping; below 22.5).
    :param normal: the sketch plane's normal; ``None`` takes the least-spread
        direction of each section's points.
    :param min_arc_points: an arc spans at least this many chain points.
    :param min_arc_sweep: and at least this many degrees.
    :param return_report: also return the report.
    :returns: the input with ``sketch:piece`` cell data and the ``sketch:*``
        field data; with ``return_report`` ``(mesh, report)``, where
        ``report["pieces"]`` lists each piece (``type``, ``start``, ``end``,
        ``centre``, ``radius``, ``sweep`` in degrees, ``residual``,
        ``num_points``, ``chain``, ``frame``) and ``report["frames"]`` each
        frame's row-major 4x4 (columns: in-plane axes, normal, origin).
    :raises ValueError: on a bad option or an input that is not a chained
        line mesh.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _core_op

    n = (0.0, 0.0, 0.0) if normal is None else tuple(float(v) for v in normal)
    if len(n) != 3:
        raise ValueError(f"{_PREFIX}normal must have three components")
    res = _core_op("fit_curves")(
        sections,
        float(tolerance),
        float(snap_angle),
        n,
        int(min_arc_points),
        float(min_arc_sweep),
    )
    out = res.pop("mesh")
    if not return_report:
        return out
    report = {
        "pieces": [
            {
                "type": str(p["type"]),
                "start": [float(v) for v in p["start"]],
                "end": [float(v) for v in p["end"]],
                "centre": [float(v) for v in p["centre"]],
                "radius": float(p["radius"]),
                "sweep": float(p["sweep"]),
                "residual": float(p["residual"]),
                "num_points": int(p["num_points"]),
                "chain": int(p["chain"]),
                "frame": int(p["frame"]),
            }
            for p in res["pieces"]
        ],
        "frames": [
            np.asarray(f, dtype=np.float64).reshape(4, 4).tolist()
            for f in res["frames"]
        ],
        "warnings": [str(w) for w in res["warnings"]],
    }
    for k in ("num_lines", "num_arcs", "num_circles"):
        report[k] = int(res[k])
    for k in ("max_residual", "tolerance"):
        report[k] = float(res[k])
    return out, report


def sketch_json(report) -> str:
    """The sketch of a ``fit_curves`` report as a JSON string."""
    return json.dumps({"frames": report["frames"], "pieces": report["pieces"]})
