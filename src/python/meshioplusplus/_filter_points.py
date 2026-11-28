"""Outlier removal and thinning of a point cloud (roadmap 9.9).

``filter_points(points, statistical_k=..., radius=..., voxel_size=...,
min_distance=...)`` applies up to four steps in a fixed order -- statistical
outlier removal, radius outlier removal, voxel-grid thinning and Poisson-disk
thinning -- and returns the surviving points as a point cloud. See
``doc/filter_points.md``.

**C++-core only, with no numpy fallback**: the neighbour searches are the
core's (exact, ties to the lower id), and the voxel and Poisson-disk steps are
defined by their visiting order, which a second implementation would have to
repeat. This module raises when ``_core`` cannot be used.
"""

from __future__ import annotations

__all__ = ["filter_points"]


def filter_points(
    points,
    *,
    statistical_k: int = 0,
    std_ratio: float = 1.0,
    radius: float = 0.0,
    min_neighbors: int = 4,
    voxel_size: float = 0.0,
    min_distance: float = 0.0,
    voxel_mode: str = "nearest",
    return_report: bool = False,
):
    """Remove outliers from, and thin, a point cloud.

    The steps run in this order, each seeing the points the one before kept;
    a step is off at its default of 0, and at least one must be on.

    :param points: a mesh (its points; cells are dropped) or an ``(n, 3)`` array.
    :param statistical_k: statistical outlier removal (Rusu et al. 2008): a
        point whose mean distance to its ``statistical_k`` nearest other points
        exceeds ``mean + std_ratio * sigma`` of those means over all points goes.
    :param std_ratio: standard deviations above the mean for the step above.
    :param radius: radius outlier removal: a point with fewer than
        ``min_neighbors`` other points within ``radius`` goes.
    :param min_neighbors: other points needed within ``radius``.
    :param voxel_size: voxel thinning: one point per occupied cube of this
        side; cubes are visited in ascending (x, y, z) key order.
    :param min_distance: Poisson-disk thinning (greedy dart throwing in index
        order, no random numbers): a point is kept when no kept point lies closer.
    :param voxel_mode: ``"nearest"`` (default: keep the real point nearest the
        cube's mean) or ``"average"`` (emit the mean coordinates and float
        arrays, with ``point_map`` set to ``-1``).
    :param return_report: also return ``{"point_map", "num_input",
        "num_removed_statistical", "num_removed_radius", "num_removed_voxel",
        "num_removed_poisson", "warnings"}``; ``point_map`` is each output row's
        input point (or ``-1`` for averaged voxel points).
    :returns: the point cloud (one ``vertex`` block, every per-point array
        gathered, Point regions following), or ``(mesh, report)``.
    :raises ValueError: on a negative or non-finite parameter, a non-finite
        coordinate, no step enabled, or a voxel size too small for the extent.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _as_mesh, _core_op

    res = _core_op("filter_points")(
        _as_mesh(points),
        int(statistical_k),
        float(std_ratio),
        float(radius),
        int(min_neighbors),
        float(voxel_size),
        float(min_distance),
        str(voxel_mode),
    )
    out = res.pop("mesh")
    return (out, res) if return_report else out
