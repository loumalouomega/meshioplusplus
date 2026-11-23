"""Streamlines of a point vector field, as a mesh of line segments (roadmap 15.4).

A dependency-free mesh *operation*: each line is traced from its seed by the
field's own velocity (``dx/dt = v``) with fixed-step RK4 or adaptive
Dormand-Prince 5(4), over the piecewise-linear interpolation of the 2- or
3-component point array on the mesh's cells (a volume mesh in its volume, a
surface on its surface). See ``doc/streamlines.md``.

**This operation is C++-core only, with no numpy fallback** (the
``geodesic_distance`` precedent): the adaptive step control branches on rounded
error estimates, which a second implementation could decide differently. This
module raises when ``_core`` cannot be used.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "streamlines",
    "TIME_NAME",
    "SPEED_NAME",
    "SEED_ID_NAME",
    "STOP_REASONS",
]

#: Point data: the time since the seed (negative backward). Twins
#: ``kStreamlineTimeName``.
TIME_NAME = "IntegrationTime"
#: Point data: the field's speed at the point.
SPEED_NAME = "Speed"
#: Point and cell data: the index of the seed of the line.
SEED_ID_NAME = "SeedId"
#: Why a line stopped in one direction; the index is the ``end_*`` code.
STOP_REASONS = ("left_domain", "stagnation", "max_length", "max_steps")


def streamlines(
    mesh,
    array,
    *,
    seeds: int = 40,
    seed_points=None,
    direction: str = "both",
    integrator: str = "rk4",
    step: float = 0.0,
    max_length: float = 0.5,
    max_steps: int = 0,
    terminal_speed: float = 0.0,
    tolerance: float = 1e-6,
    return_report: bool = False,
):
    """Streamlines of the point vector field ``array``, as a new mesh of lines.

    The result carries the points of every line, one ``line`` segment between
    consecutive points, the point data ``IntegrationTime`` (the field's time
    since the seed, negative backward), ``Speed``, ``SeedId`` and the field
    itself under its own name, and the cell data ``SeedId``.

    :param mesh: a mesh with triangle, quad, tetra, hexahedron, wedge or
        pyramid cells (never modified).
    :param array: the name of a 2- or 3-component ``point_data`` array (a
        2-component field is ``(u, v, 0)``).
    :param seeds: about this many seeds (simplex centroids ranked at equal
        measure along a golden-ratio sequence), used without ``seed_points``.
    :param seed_points: explicit seeds, shape ``(n, 3)`` (or ``(n, 2)``); a seed
        outside the mesh makes no line and is counted.
    :param direction: ``"both"``, ``"forward"`` or ``"backward"``.
    :param integrator: ``"rk4"`` (fixed step) or ``"rk45"`` (adaptive).
    :param step: the step as a length at the fastest speed; 0 is 0.005 of the
        mesh diagonal.
    :param max_length: the longest path in each direction, as a fraction of the
        diagonal.
    :param max_steps: the step limit in each direction; 0 is 20000.
    :param terminal_speed: a line stops where the speed falls to this; 0 is
        ``1e-9`` of the fastest speed.
    :param tolerance: the adaptive integrator's tolerance, a fraction of the
        diagonal.
    :param return_report: also return ``{"num_lines", "num_seeds_outside",
        "seed_ids", "end_backward", "end_forward", "dim"}``; the ``end_*``
        entries are indices into :data:`STOP_REASONS`.
    :returns: the mesh of lines, or ``(mesh, report)``.
    :raises ValueError: on an unknown or non 2-/3-component array, a mesh with
        no cell the tracer takes, or invalid options.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    try:
        from . import _core
    except ImportError:
        _core = None
    if _core is None or not hasattr(_core, "streamlines"):
        raise NotImplementedError(
            "meshio++: streamlines: this operation is C++-core only and has no "
            "pure-numpy fallback -- install a build with the compiled "
            "meshioplusplus._core extension"
        )
    flat = []
    if seed_points is not None:
        pts = np.asarray(seed_points, dtype=float)
        if pts.ndim == 1 and pts.size in (2, 3):
            pts = pts.reshape(1, -1)
        if pts.ndim != 2 or pts.shape[1] not in (2, 3):
            raise ValueError(
                "meshio++: streamlines: seed_points must have shape (n, 2) or (n, 3)"
            )
        if pts.shape[1] == 2:
            pts = np.column_stack([pts, np.zeros(len(pts))])
        flat = pts.ravel().tolist()
    res = _core.streamlines(
        mesh,
        str(array),
        int(seeds),
        flat,
        str(direction),
        str(integrator),
        float(step),
        float(max_length),
        int(max_steps),
        float(terminal_speed),
        float(tolerance),
    )
    out = res.pop("mesh")
    return (out, res) if return_report else out
