"""Materialize side regions as lower-dimensional cell blocks and back (roadmap §8.4).

Public API:
    sides_to_cells
    cells_to_sides
"""

from __future__ import annotations

from ._mesh import Mesh


def sides_to_cells(
    mesh: Mesh,
    region: str | None = None,
    *,
    return_report: bool = False,
):
    """Materialize ``Side`` region(s) as lower-dimensional cell blocks carrying
    ``cell_data["meshio:side_of"]`` ``(global_cell, local_facet)``.

    Any ``field_data["exodus:sset:<region>:<var>"]`` values are transferred onto
    the materialized facet blocks as ``cell_data["<var>"]``.
    """
    from . import _core

    out, report = _core.sides_to_cells(mesh, "" if region is None else str(region))
    if return_report:
        return out, dict(report)
    return out


def cells_to_sides(
    mesh: Mesh,
    region_name: str | None = None,
    *,
    drop_facet_blocks: bool = False,
    return_report: bool = False,
):
    """Inverse of :func:`sides_to_cells`: rebuild ``Side`` regions (and
    ``exodus:sset:<region>:<var>`` field arrays) from blocks carrying
    ``cell_data["meshio:side_of"]``.
    """
    from . import _core

    out, report = _core.cells_to_sides(
        mesh,
        "" if region_name is None else str(region_name),
        bool(drop_facet_blocks),
    )
    if return_report:
        return out, dict(report)
    return out
