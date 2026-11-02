"""A triangle surface from an oriented point cloud (roadmap 7.5.2).

``reconstruct_surface(cloud)`` contours the signed distance to the samples'
tangent planes (Hoppe et al. 1992) on a band of lattice cells about the
samples, by marching tetrahedra, and returns one ``triangle`` block. The
normals come from ``point_data["normals"]`` (as ``estimate_normals`` writes
them, and as the PLY, PCD and XYZ readers read them), or are estimated. See
``doc/reconstruct_surface.md``.

**C++-core only, with no numpy fallback** (the ``estimate_normals``
precedent): the nearest-sample ties and the cutter's serial numbering would
have to be replicated bit for bit. This module raises when ``_core`` cannot be
used.
"""

from __future__ import annotations

import numpy as np

__all__ = ["reconstruct_surface"]


def reconstruct_surface(
    points,
    *,
    normals="normals",
    k: int = 8,
    cell_size: float = 0.0,
    resolution: int = 0,
    max_distance: float = 0.0,
    estimate_normals: bool = True,
    normals_k: int = 16,
    max_cells: int = 20_000_000,
    return_report: bool = False,
):
    """A triangle surface through an oriented point cloud.

    The field at each lattice node blends the signed distances to the
    tangent planes of its ``k`` nearest samples (``k=1`` is Hoppe et al.'s
    field exactly); its zero level set, cut by marching tetrahedra, is the
    surface. A node whose projection onto its nearest sample's plane lies
    farther than ``max_distance`` from that sample is undefined, which keeps
    an open scan open instead of capping it.

    :param points: a mesh (its points; cells are ignored) or an ``(n, 3)``
        array.
    :param normals: the point-data name of the normals, or an ``(n, 3)``
        array of them. Points whose normal is not finite or is zero are
        skipped. When the named array is absent the normals are estimated
        (``estimate_normals``, ``normals_k`` neighbours, consistent
        orientation) unless ``estimate_normals=False``.
    :param k: tangent planes blended per lattice node; at least 1.
    :param cell_size: the lattice cell size; 0 takes the cloud's extent over
        ``resolution`` when that is given, else the sample spacing (the median
        distance from a sample to its fourth-nearest other sample).
    :param resolution: cells along the cloud's longest axis, with
        ``cell_size`` 0.
    :param max_distance: Hoppe's trim distance; 0 takes twice the larger of
        the spacing and the cell size.
    :param max_cells: the most lattice cells the band about the samples may
        hold.
    :param return_report: also return ``{"num_samples", "spacing",
        "cell_size", "max_distance", "num_lattice_points", "num_undefined",
        "num_triangles", "num_components", "estimated_normals"}``.
    :returns: the surface (points and one ``triangle`` block, wound the way
        the normals point, with the input's length unit), or
        ``(mesh, report)``.
    :raises ValueError: on points without three coordinates, fewer than three
        samples, coincident samples, missing normals with
        ``estimate_normals=False``, both ``cell_size`` and ``resolution``, a
        negative size or distance, or a band over ``max_cells``.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _as_mesh, _core_op
    from ._mesh import Mesh

    mesh = _as_mesh(points)
    name = normals
    if not isinstance(normals, str):
        arr = np.asarray(normals, dtype=np.float64)
        if arr.shape != (len(mesh.points), 3):
            raise ValueError(
                "meshio++: reconstruct_surface: normals must be a point-data name or an "
                f"(n, 3) array with n = {len(mesh.points)}, got shape {arr.shape}"
            )
        mesh = Mesh(
            mesh.points,
            [],
            point_data={"normals": arr},
            field_data=dict(mesh.field_data),
        )
        name = "normals"
    res = _core_op("reconstruct_surface")(
        mesh,
        str(name),
        int(k),
        float(cell_size),
        int(resolution),
        float(max_distance),
        bool(estimate_normals),
        int(normals_k),
        int(max_cells),
    )
    out = res.pop("mesh")
    res["estimated_normals"] = bool(res["estimated_normals"])
    return (out, res) if return_report else out
