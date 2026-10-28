"""Sweep a surface into a volume, or a curve into a surface (roadmap 3.1.2, 3.1.3).

:func:`extrude` copies the profile's points along a vector to ``num_layers``
levels and joins consecutive levels: a triangle sweeps to a ``wedge``, a quad
to a ``hexahedron``, a line (when the mesh has no triangle or quad) to a
``quad``. :func:`revolve` turns them around an axis instead; points on the
axis stay put and the cells around them collapse into pyramids, tetrahedra
and wedges. The output follows the parametric primitives: the swept cells
only, layer by layer, with ``Side`` regions ``bottom``/``top`` (``start``/
``end`` for a partial turn) and ``side``, and the input's regions, data and
field data carried. See ``doc/extrude.md``.

The C++ core does the work (it shares ``boundary_layer``'s layer sweep); there
is no numpy twin.
"""

from __future__ import annotations

import math

__all__ = ["extrude", "revolve"]

_REPORT_INTS = (
    "num_cells",
    "num_hexahedra",
    "num_wedges",
    "num_pyramids",
    "num_tetrahedra",
    "num_quads",
    "num_triangles",
    "num_flipped",
    "num_on_axis",
)


def _vector(value, name: str, prefix: str) -> tuple:
    try:
        v = tuple(float(c) for c in value)
    except TypeError:
        raise ValueError(f"{prefix}{name} must be three numbers") from None
    if len(v) != 3 or not all(math.isfinite(c) for c in v):
        raise ValueError(f"{prefix}{name} must be three finite numbers")
    return v


def _count(value, name: str, prefix: str) -> int:
    if isinstance(value, bool) or int(value) != value:
        raise ValueError(f"{prefix}{name} must be an integer")
    return int(value)


def _core_module(op: str):
    try:
        from . import _core
    except ImportError as exc:  # pragma: no cover - the numpy-only install
        raise NotImplementedError(
            f"meshio++: {op}: needs the compiled C++ core "
            "(meshioplusplus._core); it has no numpy twin"
        ) from exc
    return _core


def _report(res: dict, revolve: bool) -> dict:
    report = {k: int(res[k]) for k in _REPORT_INTS}
    if not revolve:
        del report["num_on_axis"]
    report["min_scaled_jacobian"] = float(res["min_scaled_jacobian"])
    if revolve:
        report["full"] = bool(res["full"])
    report["warnings"] = [str(w) for w in res["warnings"]]
    return report


def extrude(
    mesh,
    direction=(0.0, 0.0, 1.0),
    *,
    num_layers: int = 1,
    ratio: float = 1.0,
    offsets=None,
    record_layers: bool = False,
    max_cells: int = 20_000_000,
    return_report: bool = False,
):
    """Sweep a surface (or a curve) along a vector.

    :param mesh: the profile: ``triangle`` and ``quad`` cells, or, without
        them, ``line`` cells. Lower-dimensional cells (boundary lines,
        vertices) are not swept; their ``Cell`` regions become ``Side``
        regions of what they sweep to.
    :param direction: the whole displacement of the last level (not zero).
        A face is turned to point along it first; a face or line parallel to
        it is refused.
    :param num_layers: the number of layers (at least 1); with ``offsets``,
        1 or their count.
    :param ratio: each layer's thickness over the one below's (positive; 1
        is uniform).
    :param offsets: explicit level fractions of ``direction``, increasing,
        positive and ending at 1; not with a ``ratio`` other than 1.
    :param record_layers: attach the Int32 ``extrude:layer`` cell data.
    :param max_cells: refuse more than this many profile cells times layers
        (non-positive disables the budget).
    :param return_report: also return the report.
    :returns: the swept mesh (``wedge``/``hexahedron`` cells, or ``quad``
        cells for a curve; ``bottom``, ``top`` and ``side`` regions), or
        ``(mesh, report)`` with ``num_cells``, ``num_hexahedra``,
        ``num_wedges``, ``num_pyramids``, ``num_tetrahedra``, ``num_quads``,
        ``num_triangles``, ``num_flipped``, ``min_scaled_jacobian`` and
        ``warnings``.
    :raises ValueError: for a bad option, a face or line parallel to the
        direction, or an exceeded budget.
    :raises NotImplementedError: for higher-order, polygon, polyhedron or
        volume cells.
    """
    prefix = "meshio++: extrude: "
    direction = _vector(direction, "direction", prefix)
    num_layers = _count(num_layers, "num_layers", prefix)
    ratio = float(ratio)
    offs = [] if offsets is None else [float(v) for v in offsets]
    res = _core_module("extrude").extrude(
        mesh,
        direction=direction,
        num_layers=num_layers,
        ratio=ratio,
        offsets=offs,
        record_layers=bool(record_layers),
        max_cells=_count(max_cells, "max_cells", prefix),
    )
    out = res["mesh"]
    return (out, _report(res, False)) if return_report else out


def revolve(
    mesh,
    axis=(0.0, 1.0, 0.0),
    *,
    axis_point=(0.0, 0.0, 0.0),
    angle: float = 360.0,
    num_layers: int = 16,
    tolerance: float = 1e-9,
    record_layers: bool = False,
    max_cells: int = 20_000_000,
    return_report: bool = False,
):
    """Sweep a surface (or a curve) around an axis.

    :param mesh: the profile, as for :func:`extrude`.
    :param axis: the axis direction (not zero); the default turns an
        ``(r, z)`` profile drawn in the xy plane about the y axis.
    :param axis_point: a point of the axis.
    :param angle: degrees, ``0 < |angle| <= 360``; the right-hand rule about
        ``axis``. A full turn closes on itself.
    :param num_layers: steps of the turn; at least 3 for a full turn, and no
        step may reach 180 degrees.
    :param tolerance: a point within this fraction of the profile's
        bounding-box diagonal of the axis is on it, and stays put.
    :param record_layers: attach the Int32 ``revolve:layer`` cell data.
    :param max_cells: as for :func:`extrude`.
    :param return_report: also return the report.
    :returns: the swept mesh, with ``start``, ``end`` (a partial turn) and
        ``side`` regions, or ``(mesh, report)``: the :func:`extrude` keys,
        ``num_on_axis`` and ``full``.
    :raises ValueError: for a bad option, a face that crosses the axis or
        turns within its own plane, a line on the axis, or an exceeded budget.
    :raises NotImplementedError: as for :func:`extrude`.
    """
    prefix = "meshio++: revolve: "
    res = _core_module("revolve").revolve(
        mesh,
        axis=_vector(axis, "axis", prefix),
        axis_point=_vector(axis_point, "axis_point", prefix),
        angle=float(angle),
        num_layers=_count(num_layers, "num_layers", prefix),
        tolerance=float(tolerance),
        record_layers=bool(record_layers),
        max_cells=_count(max_cells, "max_cells", prefix),
    )
    out = res["mesh"]
    return (out, _report(res, True)) if return_report else out
