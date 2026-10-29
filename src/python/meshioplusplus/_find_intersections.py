"""Where a surface intersects itself, or two surfaces interpenetrate (roadmap 3.1.2).

The surfaces are reduced to triangles (quads and polygons fanned from their
first corner), candidate pairs come from a bucket grid over their bounding
boxes, and every pair is decided by an exact triangle-triangle test on
Shewchuk's ``orient3d``/``orient2d``: every reported hit is exact. On one
surface, triangles that share corners are an intersection only when they meet
beyond them; on two surfaces any contact counts, touching included. See
``doc/find_intersections.md``.

The C++ core does the work; there is no numpy twin.
"""

from __future__ import annotations

__all__ = ["find_intersections"]

_PREFIX = "meshio++: find_intersections: "
_REPORT_INTS = (
    "num_pairs",
    "num_triangle_pairs",
    "num_coplanar",
    "num_segments",
    "num_cells_a",
    "num_cells_b",
)


def find_intersections(
    mesh,
    other=None,
    *,
    region: str | None = None,
    region_other: str | None = None,
    attach_flags: bool = True,
    segments: bool = True,
    max_pairs: int = 0,
    return_report: bool = False,
):
    """Find where ``mesh`` intersects itself, or where it and ``other`` meet.

    :param mesh: a surface mesh (triangles, quads, polygons; lines and
        vertices are ignored). Run :func:`~meshioplusplus.extract_surface` on a
        volume mesh first.
    :param other: a second surface; ``None`` checks ``mesh`` against itself.
    :param region: restrict ``mesh`` to this named ``Cell`` region.
    :param region_other: restrict ``other`` to this named ``Cell`` region.
    :param attach_flags: attach the Int64 cell data ``intersections:count``,
        how many cells each cell intersects.
    :param segments: build the line mesh of the intersection segments.
    :param max_pairs: keep at most this many pairs in ``report["pairs"]``
        (0 keeps all); the counts, flags and segments still cover every pair.
    :param return_report: also return the report.
    :returns: ``mesh``, bit-identical but for the flags; with
        ``return_report`` ``(mesh, report)``, where ``report`` has ``pairs``
        (an ``(n, 2)`` Int64 array of global cell ids, the cell of ``mesh``
        first), ``segments`` (a line mesh with ``intersections:cell_a`` and
        ``intersections:cell_b`` cell data, or ``None``), ``other`` (``other``
        with its flags, or ``None``), ``num_pairs``, ``num_triangle_pairs``,
        ``num_coplanar``, ``num_segments``, ``num_cells_a``, ``num_cells_b``,
        ``truncated``, ``clean`` and ``warnings``.
    :raises ValueError: for a volume, polyhedron or higher-order block, an
        unknown region or a negative ``max_pairs``.
    """
    try:
        from . import _core
    except ImportError as exc:  # pragma: no cover - the numpy-only install
        raise NotImplementedError(
            "meshio++: find_intersections: needs the compiled C++ core "
            "(meshioplusplus._core); it has no numpy twin"
        ) from exc
    max_pairs = int(max_pairs)
    if max_pairs < 0:
        raise ValueError(f"{_PREFIX}max_pairs must be >= 0, got {max_pairs}")
    res = _core.find_intersections(
        mesh,
        other,
        region=region or "",
        region_other=region_other or "",
        attach_flags=bool(attach_flags),
        segments=bool(segments),
        max_pairs=max_pairs,
    )
    out = res["mesh"]
    if not return_report:
        return out
    report = {k: int(res[k]) for k in _REPORT_INTS}
    report["pairs"] = res["pairs"]
    report["segments"] = res["segments"]
    report["other"] = res["mesh_b"]
    report["truncated"] = bool(res["truncated"])
    report["clean"] = bool(res["clean"])
    report["warnings"] = [str(w) for w in res["warnings"]]
    return out, report
