"""Hanging nodes, T-junctions and partly overlapping facets of a mesh.

A thin wrapper over ``operations/check_conformity.hpp``. The check has no numpy
twin: its geometric tests (points on facets, coplanar overlap) would be a second
implementation of the same tolerances to keep in step, and the wheels always
carry the core, so the Python entry point names the core when it is missing
rather than answering differently.
"""

from __future__ import annotations

import math

_PREFIX = "meshio++: check_conformity: "

#: ``conformity:hanging`` values.
CONFORMING, HANGING_NODE, T_JUNCTION = 0, 1, 2

_REPORT_INTS = (
    "dimension",
    "num_hanging_nodes",
    "num_t_junctions",
    "num_partial_overlaps",
    "num_unwelded_facets",
    "num_boundary_facets",
    "num_non_manifold_facets",
)


def check_conformity(
    mesh, tolerance: float = 1e-9, attach_flags: bool = True, return_report=False
):
    """Report whether ``mesh`` is conforming, and where it is not.

    The mesh's volume cells (else its 2-D cells) are examined; lines and
    vertices are ignored. Only facets one cell uses can hold a defect:

    - a **hanging node** (3-D) is a point some cell uses that lies on a face of a
      volume cell -- inside it or on one of its edges -- without being one of
      that cell's nodes;
    - a **T-junction** (2-D) is a point some cell uses on the interior of an
      edge of a 2-D cell that does not use it;
    - a **partial overlap** is a pair of unmatched facets in one plane (on one
      line) that overlap over a positive area (length) without being one facet;
    - an **unwelded** pair is two unmatched facets whose corners coincide but
      whose point ids differ.

    :param mesh: any mesh.
    :param tolerance: distances at or below this fraction of the bounding-box
        diagonal of the examined facets count as zero.
    :param attach_flags: attach the Int64 ``conformity:hanging`` point data (0,
        1 for a hanging node, 2 for a T-junction).
    :param return_report: also return the counts and the verdict.
    :returns: the mesh, bit-identical but for the flags; with ``return_report``
        ``(mesh, report)``, where ``report`` has ``dimension``,
        ``num_hanging_nodes``, ``num_t_junctions``, ``num_partial_overlaps``,
        ``num_unwelded_facets``, ``num_boundary_facets``,
        ``num_non_manifold_facets``, ``conforming``, ``interface_mode``
        (``"conforming"``, ``"weld"`` or ``"proximity"``: which
        ``find_interface`` mode suits the mesh) and ``warnings``.
    """
    tol = float(tolerance)
    if not (math.isfinite(tol) and tol >= 0.0):
        raise ValueError(f"{_PREFIX}the tolerance must be finite and non-negative")
    try:
        from . import _core
    except ImportError as exc:
        raise NotImplementedError(
            f"{_PREFIX}needs the compiled C++ core (meshioplusplus._core); "
            "it has no numpy twin"
        ) from exc
    res = _core.check_conformity(mesh, tol, bool(attach_flags))
    out = res["mesh"]
    if not return_report:
        return out
    report = {k: int(res[k]) for k in _REPORT_INTS}
    report["conforming"] = bool(res["conforming"])
    report["interface_mode"] = str(res["interface_mode"])
    report["warnings"] = [str(w) for w in res["warnings"]]
    return out, report


__all__ = ["check_conformity"]
