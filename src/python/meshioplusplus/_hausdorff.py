"""The (sampled) Hausdorff distance between two surfaces.

``hausdorff_distance(a, b)`` samples each surface -- its vertices, plus the
centroids of the ``s * s`` sub-triangles of every triangle when
``face_samples = s > 0`` -- and measures every sample's unsigned distance to the
*other* surface. The one-sided distances are the largest sample distance each
way; the Hausdorff distance is the larger of the two. Vertex sampling alone is
exact when the farthest point is a vertex, and a lower bound otherwise.

The C++ core (``operations/hausdorff.hpp``) does the work; the numpy twin below
is a brute-force reference (every sample against every triangle) for
pure-Python installs.
"""

from __future__ import annotations

import numpy as np

from ._fallback import core_op_declined

_PREFIX = "meshio++: hausdorff_distance: "


def _has_volume(mesh):
    from ._mesh import topological_dimension

    return any(
        cb.type.startswith("polyhedron") or topological_dimension.get(cb.type, 0) == 3
        for cb in mesh.cells
    )


def _surface(mesh, region, side):
    from ._normals import _region_mask
    from ._sdf import _soup

    if _has_volume(mesh):
        if region:
            raise ValueError(
                f"{_PREFIX}mesh {side} is a volume mesh, whose skin is compared; a region "
                "selects surface cells (run extract_surface first)"
            )
        from ._skin import extract_skin

        mesh = extract_skin(mesh, linearize=True)
    points, verts, corners, source = _soup(mesh)
    if region:
        mask = _region_mask(mesh, region)
        if not mask.any() and not any(
            r.name == region and r.kind == "cell"
            for r in getattr(mesh, "regions", ()) or ()
        ):
            raise ValueError(f"{_PREFIX}no cell region named '{region}'")
        keep = mask[source]
        verts, corners = verts[keep], corners[keep]
    if len(verts) == 0:
        raise ValueError(f"{_PREFIX}mesh {side} has no surface triangles to measure")
    return points, verts, corners


def _samples(points, verts, corners, s):
    used = np.unique(verts.ravel())
    out = [points[used]]
    if s > 0:
        a = corners[:, 0]
        ab = corners[:, 1] - a
        ac = corners[:, 2] - a
        uv = []
        for i in range(s):
            for j in range(s - i):
                uv.append((i + 1.0 / 3.0, j + 1.0 / 3.0))
                if i + j + 1 < s:
                    uv.append((i + 2.0 / 3.0, j + 2.0 / 3.0))
        uv = np.asarray(uv)
        pts = (
            a[:, None, :]
            + ab[:, None, :] * (uv[None, :, 0:1] / s)
            + ac[:, None, :] * (uv[None, :, 1:2] / s)
        )
        out.append(pts.reshape(-1, 3))
    return np.concatenate(out)


def _one_sided(samples, corners):
    from ._sdf import _closest_points

    d = np.empty(len(samples))
    for q in range(len(samples)):
        _point, dist2, _feature = _closest_points(samples[q], corners)
        d[q] = np.sqrt(dist2.min())
    worst = int(np.argmax(d)) if len(d) else 0
    return {
        "max": float(d[worst]) if len(d) else 0.0,
        "mean": float(d.mean()) if len(d) else 0.0,
        "rms": float(np.sqrt(np.mean(d * d))) if len(d) else 0.0,
        "worst": tuple(float(x) for x in samples[worst]) if len(d) else (0.0, 0.0, 0.0),
        "n": len(d),
    }


def _hausdorff_py(a, b, face_samples, region_a, region_b):
    pa, va, ca = _surface(a, region_a, "A")
    pb, vb, cb = _surface(b, region_b, "B")
    ab = _one_sided(_samples(pa, va, ca, face_samples), cb)
    ba = _one_sided(_samples(pb, vb, cb, face_samples), ca)
    return {
        "distance": max(ab["max"], ba["max"]),
        "a_to_b": ab["max"],
        "b_to_a": ba["max"],
        "mean_a_to_b": ab["mean"],
        "rms_a_to_b": ab["rms"],
        "mean_b_to_a": ba["mean"],
        "rms_b_to_a": ba["rms"],
        "num_samples_a": ab["n"],
        "num_samples_b": ba["n"],
        "worst_point_a": ab["worst"],
        "worst_point_b": ba["worst"],
    }


def hausdorff_distance(
    a,
    b,
    face_samples: int = 0,
    region_a: str = "",
    region_b: str = "",
    grid_cell_size: float = 0.0,
) -> dict:
    """The (sampled) Hausdorff distance between the surfaces of two meshes.

    :param a: a surface mesh, or a volume mesh whose skin is compared.
    :param b: the same, for the other side.
    :param face_samples: 0 samples the vertices only; ``s > 0`` also samples the
        centroids of the ``s * s`` sub-triangles of every triangle, tightening
        the (lower-bound) estimate.
    :param region_a: restrict ``a`` to this named cell region of surface cells.
    :param region_b: the same, for ``b``.
    :param grid_cell_size: bucket size of the nearest-triangle search; 0 picks
        one automatically.
    :returns: ``{"distance", "a_to_b", "b_to_a", "mean_a_to_b", "rms_a_to_b",
        "mean_b_to_a", "rms_b_to_a", "num_samples_a", "num_samples_b",
        "worst_point_a", "worst_point_b"}``.
    """
    s = int(face_samples)
    if s < 0:
        raise ValueError(f"{_PREFIX}face samples must be >= 0")
    if grid_cell_size < 0:
        raise ValueError(f"{_PREFIX}the grid cell size must be >= 0")
    try:
        from . import _core

        return dict(
            _core.hausdorff_distance(
                a, b, s, str(region_a), str(region_b), float(grid_cell_size)
            )
        )
    except Exception as exc:
        if not core_op_declined(exc, "hausdorff_distance"):
            raise
    return _hausdorff_py(a, b, s, str(region_a), str(region_b))


__all__ = ["hausdorff_distance"]
