"""Boundary-driven mesh morphing (roadmap 7.4.1): move the boundary (or
chosen points) by a prescribed displacement and carry the motion into the
rest of the mesh by a harmonic extension with Jacobian-based stiffening,
under an inversion guard. See ``doc/morph.md``.

**C++-core only, with no numpy fallback** (the ``geodesic_distance``
precedent): the extension is a conjugate-gradient solve that stops on a
rounded test, so a second implementation could stop a step early or late.
"""

from __future__ import annotations

import numpy as np

__all__ = ["morph"]

_TEMP = "morph:displacement"


def morph(
    mesh,
    displacement,
    *,
    fix_boundary: bool = True,
    constrained=None,
    stiffening: float = 1.0,
    steps: int = 1,
    allow_inverted: bool = False,
    tolerance: float = 1e-10,
    max_iterations: int = 0,
    return_report: bool = False,
):
    """Move a mesh by a prescribed displacement of its boundary (or of chosen
    points), carried smoothly into the rest of the mesh.

    The **constrained** points -- the boundary (``fix_boundary``), the points
    ``constrained`` names, and any point no surface or volume cell uses --
    move by ``displacement`` exactly; every other point moves by the harmonic
    extension of those values over the mesh's P1 stiffness (every linear
    cell type, decomposed into simplices for the assembly), each simplex
    stiffened by ``(mean measure / its measure) ** stiffening`` so small cells
    move nearly rigidly. Give the boundary points that should stay put a zero
    displacement.

    :param mesh: a mesh of linear volume (or surface) cells.
    :param displacement: the name of a point-data array, or an ``(n, dim)`` /
        ``(n, 3)`` array, one row per point. Only its constrained rows are
        read.
    :param fix_boundary: constrain the boundary points.
    :param constrained: more points to constrain: the name of a
        ``point_set``, a point-data array whose nonzero entries pin, or a
        boolean mask with one entry per point.
    :param stiffening: the stiffening exponent (0: plain Laplace).
    :param steps: apply the displacement in this many increments, the
        stiffness reassembled on the moved mesh each time (large motions).
    :param allow_inverted: return a mesh with inverted cells instead of
        raising.
    :param tolerance: the relative residual of each solve.
    :param max_iterations: each solve's cap; 0 picks one from the size.
    :param return_report: also return ``{"num_constrained", "num_free",
        "num_inverted", "iterations", "max_displacement", "converged"}``.
    :returns: the moved mesh (connectivity, data and regions unchanged), or
        ``(mesh, report)``.
    :raises ValueError: on a missing or misshapen displacement, a negative
        stiffening, ``steps`` below 1, a mesh with no surface or volume cell,
        a polyhedron or higher-order cell, or -- unless ``allow_inverted`` --
        a result with inverted cells.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _core_op
    from ._mesh import Mesh

    fn = _core_op("morph")
    n = len(mesh.points)
    work = mesh
    array = displacement
    if not isinstance(displacement, str):
        d = np.asarray(displacement, dtype=np.float64)
        if d.ndim != 2 or d.shape[0] != n:
            raise ValueError(
                "meshio++: morph: the displacement must be (n, dim) or (n, 3), one "
                f"row per point; got shape {d.shape}"
            )
        work = Mesh(
            mesh.points,
            mesh.cells,
            point_data={**mesh.point_data, _TEMP: d},
            cell_data=dict(mesh.cell_data),
            field_data=dict(mesh.field_data),
            point_sets=dict(mesh.point_sets),
            cell_sets=dict(mesh.cell_sets),
        )
        if getattr(mesh, "regions", None):
            work.regions = list(mesh.regions)
        array = _TEMP
    mask = []
    constrained_array = ""
    if isinstance(constrained, str):
        sets = getattr(mesh, "point_sets", None) or {}
        if constrained in sets:
            m = np.zeros(n, dtype=np.uint8)
            m[np.asarray(sets[constrained], dtype=np.int64).ravel()] = 1
            mask = m.tolist()
        else:
            constrained_array = constrained
    elif constrained is not None:
        m = np.asarray(constrained).ravel()
        if m.shape[0] != n:
            raise ValueError(
                "meshio++: morph: the constrained mask must have one entry per point"
            )
        mask = (m != 0).astype(np.uint8).tolist()
    res = fn(
        work,
        array,
        bool(fix_boundary),
        mask,
        constrained_array,
        float(stiffening),
        int(steps),
        bool(allow_inverted),
        float(tolerance),
        int(max_iterations),
    )
    out = res.pop("mesh")
    if array == _TEMP:
        out.point_data.pop(_TEMP, None)
    if not res["converged"]:
        from ._common import warn

        warn("morph: a solve hit its iteration cap; the motion is approximate")
    return (out, res) if return_report else out
