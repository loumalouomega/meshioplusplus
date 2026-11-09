"""Drop a scan's crumbs (roadmap 7.8.1): the connected components whose
bounding-box diagonal, area or cell count is below a threshold. C++-core only
(see ``doc/remove_islands.md``)."""

from __future__ import annotations

__all__ = ["remove_islands"]


def remove_islands(
    mesh,
    *,
    connectivity: str = "edge",
    min_diagonal: float = 0.01,
    min_area: float = 0.0,
    min_cells: int = 0,
    relative: bool = True,
    dry_run: bool = False,
    return_report: bool = False,
):
    """Remove the small connected components of ``mesh``.

    :param connectivity: ``"edge"`` (cells sharing a facet) or ``"point"``.
    :param min_diagonal: a component whose bounding-box diagonal is below
        this goes; with ``relative`` a fraction of the whole mesh's diagonal.
    :param min_area: likewise for the area of its 2-D cells (0: off).
    :param min_cells: likewise for its cell count (0: off).
    :param relative: diagonal and area as fractions of the whole mesh's.
    :param dry_run: remove nothing; label the cells with ``islands:component``
        and ``islands:removed`` instead.
    :param return_report: also return the report: ``num_components``,
        ``num_removed``, ``num_cells_removed``, ``num_points_removed`` and per
        component ``component_cells``, ``component_diagonal``,
        ``component_area`` and ``component_removed``.
    :raises ValueError: on a negative or non-finite threshold.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _core_op

    res = _core_op("remove_islands")(
        mesh,
        str(connectivity),
        float(min_diagonal),
        float(min_area),
        int(min_cells),
        bool(relative),
        bool(dry_run),
    )
    out = res.pop("mesh")
    if not return_report:
        return out
    report = {
        k: int(res[k])
        for k in (
            "num_components",
            "num_removed",
            "num_cells_removed",
            "num_points_removed",
        )
    }
    report["component_cells"] = [int(v) for v in res["component_cells"]]
    report["component_diagonal"] = [float(v) for v in res["component_diagonal"]]
    report["component_area"] = [float(v) for v in res["component_area"]]
    report["component_removed"] = [bool(v) for v in res["component_removed"]]
    return out, report
