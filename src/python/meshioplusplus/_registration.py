"""Rigid registration by iterative closest points (roadmap 7.5.3).

``register_rigid(source, target)`` returns the 4x4 matrix that lays the
source onto the target -- apply it with ``transform(source, matrix=m)`` --
and a report. Point-to-point (Horn's closed form) or point-to-plane
(linearised, applied as an exact rotation); closest points on the target's
surface, or its nearest points for a cloud. See ``doc/registration.md``.

**C++-core only, with no numpy fallback** (the ``geodesic_distance``
precedent): the iteration stops on a rounded test, so a second
implementation could stop a step early or late and return a different
matrix in the last digits. This module raises when ``_core`` cannot be used.
"""

from __future__ import annotations

import numpy as np

__all__ = ["register_rigid"]

_METHODS = ("point_to_point", "point_to_plane")


def register_rigid(
    source,
    target,
    *,
    method: str = "point_to_point",
    max_iterations: int = 50,
    tolerance: float = 1e-10,
    max_distance: float = 0.0,
    initial=None,
):
    """The rigid motion laying ``source`` onto ``target``, by iterative
    closest points (Besl and McKay 1992).

    :param source: the mesh or ``(n, 3)`` points to move (only its points are
        read).
    :param target: the mesh to lay it onto: closest points on its surface when
        it has surface cells (a volume cell is refused naming
        ``extract_surface``), else its nearest points (a cloud).
    :param method: ``"point_to_point"`` or ``"point_to_plane"`` (the distance
        to the target's tangent planes; usually far fewer iterations; a cloud
        target needs ``normals`` point data, see ``estimate_normals``).
    :param max_iterations: the iteration cap.
    :param tolerance: stop when a step rotates by less than this (radians)
        and moves the source's centre by less than this times its
        bounding-box diagonal.
    :param max_distance: drop pairs farther apart than this; 0 keeps all.
    :param initial: a starting 4x4 rigid matrix (row-major), or ``None``.
    :returns: ``(matrix, report)``: the 4x4 ``numpy`` matrix mapping source
        coordinates onto the target's (``transform(source, matrix=matrix)``
        applies it), and ``{"initial_rms", "rms", "iterations", "num_pairs",
        "converged"}``. ICP finds the nearest local optimum: start roughly
        aligned.
    :raises ValueError: on an unknown method, an empty or non-finite input,
        a target with volume cells, ``point_to_plane`` against a cloud with no
        ``normals``, or no pair within ``max_distance``.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _as_mesh, _core_op

    if method not in _METHODS:
        raise ValueError(
            f"meshio++: register_rigid: unknown method '{method}' (expected "
            "'point_to_point' or 'point_to_plane')"
        )
    if initial is None:
        start = np.eye(4)
    else:
        start = np.asarray(initial, dtype=np.float64).reshape(4, 4)
    res = _core_op("register_rigid")(
        _as_mesh(source),
        _as_mesh(target),
        method,
        int(max_iterations),
        float(tolerance),
        float(max_distance),
        [float(v) for v in start.ravel()],
    )
    matrix = np.asarray(res.pop("matrix"), dtype=np.float64).reshape(4, 4)
    if not res["converged"]:
        from ._common import warn

        warn(
            f"register_rigid: no convergence in {res['iterations']} iterations "
            f"(rms {res['rms']:.3g}); raise max_iterations or start closer"
        )
    return matrix, res
