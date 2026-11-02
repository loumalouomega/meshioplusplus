"""Untangling inverted tetrahedra and hexahedra (roadmap 7.4.2): move nodes
until every corner has a positive signed volume, by Escobar et al.'s
untangling and smoothing functional. See ``doc/untangle.md``.

**C++-core only, with no numpy fallback** (the ``geodesic_distance``
precedent): a node-by-node descent whose line search accepts on a rounded
comparison, so a second implementation could take a different step.
"""

from __future__ import annotations

import numpy as np

__all__ = ["untangle"]


def untangle(
    mesh,
    *,
    fix_boundary: bool = True,
    frozen=None,
    rings: int = 2,
    smooth: bool = False,
    max_sweeps: int = 100,
    return_report: bool = False,
):
    """Move nodes until every tetrahedron and hexahedron corner has a positive
    signed volume.

    Minimises Escobar et al.'s untangling and smoothing functional node by
    node: a regularised signed volume keeps it finite across inversions and
    pushes every corner positive. By default only the nodes within ``rings``
    edge rings of an inverted cell move, and it stops once none is left;
    ``smooth=True`` moves every free node until the functional settles.

    :param mesh: a volume mesh; its tetrahedra and hexahedra are measured.
        The nodes of other volume cells, and the boundary nodes with
        ``fix_boundary``, stay where they are.
    :param fix_boundary: keep the nodes of faces used by one cell fixed.
    :param frozen: more nodes to keep fixed: a ``point_set`` name, a point-data
        array (nonzero pins) or a boolean mask, one entry per point.
    :param rings: edge rings around the inverted cells whose nodes may move.
    :param smooth: untangle and smooth the whole mesh.
    :param max_sweeps: sweeps at most.
    :param return_report: also return ``{"num_inverted_before",
        "num_inverted_after", "min_quality_before", "min_quality_after",
        "sweeps", "moves"}`` (quality: the least corner scaled Jacobian).
    :returns: the mesh with its nodes moved, or ``(mesh, report)``.
    :raises ValueError: on a mesh with no tetrahedron or hexahedron, or bad
        options.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _core_op

    fn = _core_op("untangle")
    n = len(mesh.points)
    mask = []
    if isinstance(frozen, str):
        sets = getattr(mesh, "point_sets", None) or {}
        if frozen in sets:
            m = np.zeros(n, dtype=np.uint8)
            m[np.asarray(sets[frozen], dtype=np.int64).ravel()] = 1
            mask = m.tolist()
        elif frozen in mesh.point_data:
            mask = (np.asarray(mesh.point_data[frozen]).ravel() != 0).astype(np.uint8)
            mask = mask.tolist()
        else:
            raise ValueError(
                f"meshio++: untangle: no point_set or point data named '{frozen}'"
            )
    elif frozen is not None:
        m = np.asarray(frozen).ravel()
        if m.shape[0] != n:
            raise ValueError(
                "meshio++: untangle: the frozen mask must have one entry per point"
            )
        mask = (m != 0).astype(np.uint8).tolist()
    res = fn(mesh, bool(fix_boundary), mask, int(rings), bool(smooth), int(max_sweeps))
    out = res.pop("mesh")
    if res["num_inverted_after"]:
        from ._common import warn

        warn(
            f"untangle: {res['num_inverted_after']} cell(s) still inverted after "
            f"{res['sweeps']} sweep(s); free more nodes (fix_boundary=False, more rings) "
            "or raise max_sweeps"
        )
    return (out, res) if return_report else out
