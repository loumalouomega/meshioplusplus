"""Inverted, degenerate and overlapping tetrahedra of a mesh.

A thin wrapper over ``operations/check_tetrahedra.hpp``, with a numpy twin for a
pure-Python install. The tolerance and the overlap rule are the ones the
volume fix of :func:`~meshioplusplus.octree_mesh` uses.
"""

from __future__ import annotations

import math

import numpy as np

from ._fallback import core_op_declined

_PREFIX = "meshio++: check_tetrahedra: "

#: The faces of a tetrahedron, each opposite the node of the same index.
_FACES = ((0, 1, 3), (1, 2, 3), (2, 0, 3), (0, 2, 1))


def check_tetrahedra(mesh, zero_volume_factor: float = 1e-6) -> dict:
    """Report whether the tetrahedra of ``mesh`` are usable.

    The ``tetra`` and ``tetra10`` cells (a quadratic tetrahedron by its four
    corners) are examined; the mesh is not changed. A tetrahedron is

    - **inverted** when its signed volume is below minus the tolerance,
    - **degenerate** when its volume is within the tolerance, which is
      ``zero_volume_factor`` times the cube of its longest edge plus the effect
      of six significant digits of written coordinates;

    and two tetrahedra **overlap** when they share a face (no other tetrahedron
    uses it) and their opposite nodes lie on the same side of it. A conforming
    mesh of a valid solid has none, whatever its materials.

    :param mesh: any mesh.
    :param zero_volume_factor: the volume tolerance, as a fraction of the
        longest edge cubed; finite and non-negative.
    :returns: a dict with ``num_tetrahedra``, ``num_inverted``,
        ``num_degenerate``, ``num_overlaps``, ``min_volume``, ``max_volume``,
        ``valid`` (no inverted or degenerate tetrahedron and no overlap) and
        ``warnings``.
    :raises ValueError: on a bad ``zero_volume_factor``.
    """
    factor = float(zero_volume_factor)
    if not (math.isfinite(factor) and factor >= 0.0):
        raise ValueError(f"{_PREFIX}zero_volume_factor must be finite and non-negative")
    try:
        from . import _core

        res = _core.check_tetrahedra(mesh, factor)
    except Exception as exc:
        if not core_op_declined(exc, "check_tetrahedra"):
            raise
        res = _check_tetrahedra_py(mesh, factor)
    return {
        "num_tetrahedra": int(res["num_tetrahedra"]),
        "num_inverted": int(res["num_inverted"]),
        "num_degenerate": int(res["num_degenerate"]),
        "num_overlaps": int(res["num_overlaps"]),
        "min_volume": float(res["min_volume"]),
        "max_volume": float(res["max_volume"]),
        "valid": bool(res["valid"]),
        "warnings": [str(w) for w in res["warnings"]],
    }


def _check_tetrahedra_py(mesh, factor):
    points = np.asarray(mesh.points, dtype=float)
    if points.ndim != 2 or points.shape[1] < 3:
        pad = np.zeros((len(points), 3))
        pad[:, : points.shape[1] if points.ndim == 2 else 0] = points
        points = pad
    blocks = []
    warnings = []
    for index, block in enumerate(mesh.cells):
        if block.type not in ("tetra", "tetra10") or len(block.data) == 0:
            continue
        try:
            data = np.asarray(block.data, dtype=np.int64)
        except (ValueError, TypeError):
            warnings.append(
                f"block {index} ('{block.type}') is ragged and was not examined"
            )
            continue
        blocks.append(data[:, :4])
    result = {
        "num_tetrahedra": 0,
        "num_inverted": 0,
        "num_degenerate": 0,
        "num_overlaps": 0,
        "min_volume": 0.0,
        "max_volume": 0.0,
        "valid": True,
        "warnings": warnings,
    }
    if not blocks:
        warnings.append("the mesh has no tetrahedra; nothing was examined")
        return result
    tets = np.concatenate(blocks)
    if tets.min() < 0 or tets.max() >= len(points):
        raise ValueError(f"{_PREFIX}a cell refers to a point out of range")

    corners = points[tets]  # (n, 4, 3)
    a = corners[:, 0]
    volume = (
        np.einsum(
            "ij,ij->i",
            np.cross(corners[:, 1] - a, corners[:, 2] - a),
            corners[:, 3] - a,
        )
        / 6.0
    )
    longest = np.zeros(len(tets))
    for i in range(4):
        for j in range(i + 1, 4):
            longest = np.maximum(
                longest, np.linalg.norm(corners[:, j] - corners[:, i], axis=1)
            )
    tolerance = (
        factor * longest**3 + 1e-6 * np.abs(corners).max(axis=(1, 2)) * longest**2
    )
    inverted = volume < -tolerance
    degenerate = ~inverted & (volume <= tolerance)

    # Faces as sorted node triples with their tetrahedron and opposite node.
    keys = np.concatenate([np.sort(tets[:, list(f)], axis=1) for f in _FACES])
    opposite = np.concatenate([tets[:, i] for i in range(4)])
    order = np.lexsort((keys[:, 2], keys[:, 1], keys[:, 0]))
    keys, opposite = keys[order], opposite[order]
    # The faces used by exactly two tetrahedra: runs of length two.
    change = np.any(keys[1:] != keys[:-1], axis=1)
    starts = np.concatenate(([0], np.flatnonzero(change) + 1))
    counts = np.diff(np.concatenate((starts, [len(keys)])))
    pair = starts[counts == 2]
    overlaps = 0
    if len(pair):
        origin = points[keys[pair, 0]]
        normal = np.cross(
            points[keys[pair, 1]] - origin, points[keys[pair, 2]] - origin
        )
        side0 = np.einsum("ij,ij->i", normal, points[opposite[pair]] - origin)
        side1 = np.einsum("ij,ij->i", normal, points[opposite[pair + 1]] - origin)
        overlaps = int(np.count_nonzero(side0 * side1 > 0.0))

    result.update(
        num_tetrahedra=len(tets),
        num_inverted=int(inverted.sum()),
        num_degenerate=int(degenerate.sum()),
        num_overlaps=overlaps,
        min_volume=float(volume.min()),
        max_volume=float(volume.max()),
    )
    result["valid"] = not (
        result["num_inverted"] or result["num_degenerate"] or overlaps
    )
    return result


__all__ = ["check_tetrahedra"]
