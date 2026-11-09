"""Reverse the winding of cells (roadmap 7.8.2). C++-core only (see
``doc/flip_normals.md``)."""

from __future__ import annotations

__all__ = ["flip_normals"]


def flip_normals(mesh, region: str | None = None, *, return_report: bool = False):
    """Reverse the winding of every cell of ``mesh``, or of the cell region
    ``region``, with no point moved; ``Side`` region entries move with their
    cells. Flipping twice gives the input back.

    :returns: the flipped mesh; with ``return_report`` ``(mesh, report)``,
        ``report["num_flipped"]`` the cells flipped.
    :raises ValueError: on an unknown region.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _core_op

    res = _core_op("flip_normals")(mesh, region or "")
    out = res["mesh"]
    return (out, {"num_flipped": int(res["num_flipped"])}) if return_report else out
