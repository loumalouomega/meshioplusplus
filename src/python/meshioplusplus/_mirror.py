"""Reflect a mesh across a plane, keeping every cell right side out (roadmap 3.2.5).

A reflecting :func:`~meshioplusplus.transform` leaves every cell inverted;
:func:`mirror` reflects the points and permutes each cell's nodes so its
orientation is valid again (mid-edge, mid-face and centre nodes stay in their
slots; polygon and polyhedron rings are reversed), and moves ``Side`` region
entries to their facets' new numbers. With ``keep_original`` the result is the
input followed by its image, the points on the plane shared (``weld``), so a
half model becomes a whole one. See ``doc/mirror.md``.

The C++ core does the work; there is no numpy twin.
"""

from __future__ import annotations

import math

__all__ = ["mirror"]

_PREFIX = "meshio++: mirror: "
_REPORT_INTS = (
    "num_flipped_cells",
    "num_points_on_plane",
    "num_dropped_cells",
    "num_dropped_sides",
)


def _vector(value, name: str) -> tuple:
    try:
        v = tuple(float(c) for c in value)
    except TypeError:
        raise ValueError(f"{_PREFIX}{name} must be three numbers") from None
    if len(v) != 3 or not all(math.isfinite(c) for c in v):
        raise ValueError(f"{_PREFIX}{name} must be three finite numbers")
    return v


def mirror(
    mesh,
    origin=(0.0, 0.0, 0.0),
    normal=(1.0, 0.0, 0.0),
    *,
    keep_original: bool = False,
    weld: bool = True,
    tolerance: float = 1e-9,
    reflect_vector_data: bool = False,
    return_report: bool = False,
):
    """Reflect a mesh across the plane through ``origin`` with ``normal``.

    :param mesh: any mesh. Two-dimensional points stay two-dimensional, so the
        normal's z component must then be 0.
    :param origin: a point of the plane.
    :param normal: the plane's normal (any non-zero length).
    :param keep_original: return the input followed by its image; each block
        holds its own cells, then their reflected copies.
    :param weld: with ``keep_original``, share the points on the plane
        between both halves: a copy lying wholly on the plane is dropped, and
        ``Side`` entries on it (now interior) are dropped, their regions kept.
    :param tolerance: "on the plane": within this fraction of the points'
        bounding-box diagonal.
    :param reflect_vector_data: reflect float 3-column (vector) and 9-column
        (tensor) point and cell data on the reflected copies.
    :param return_report: also return the report.
    :returns: the mirrored mesh, or ``(mesh, report)`` with
        ``num_flipped_cells``, ``num_points_on_plane``, ``num_dropped_cells``,
        ``num_dropped_sides`` and ``warnings``.
    :raises ValueError: for a zero normal, a bad tolerance, or a plane that
        leaves the points' dimension.
    :raises NotImplementedError: for a cell type with no orientation flip
        (``triangle10`` and beyond, the VTK Lagrange and Bezier families).
    """
    try:
        from . import _core
    except ImportError as exc:  # pragma: no cover - the numpy-only install
        raise NotImplementedError(
            "meshio++: mirror: needs the compiled C++ core "
            "(meshioplusplus._core); it has no numpy twin"
        ) from exc
    res = _core.mirror(
        mesh,
        origin=_vector(origin, "origin"),
        normal=_vector(normal, "normal"),
        keep_original=bool(keep_original),
        weld=bool(weld),
        tolerance=float(tolerance),
        reflect_vector_data=bool(reflect_vector_data),
    )
    out = res["mesh"]
    if not return_report:
        return out
    report = {k: int(res[k]) for k in _REPORT_INTS}
    report["warnings"] = [str(w) for w in res["warnings"]]
    return out, report
