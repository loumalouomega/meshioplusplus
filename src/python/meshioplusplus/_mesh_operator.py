"""The mesh's own finite-element operators as sparse matrices (roadmap 15.2.1).

The cotangent Laplacian, the lumped and the consistent mass and the
per-triangle gradient of a triangle surface -- the P1 operators
``geodesic_distance``, ``laplacian_eigenmodes`` and ``sobolev_deform`` solve
with -- returned as data. See ``doc/mesh_operator.md``.

**C++-core only, with no numpy fallback** (the ``geodesic_distance``
precedent): the operators are assembled by the core's shared surface kernel,
and a second implementation would only be a second place to keep in step.
"""

from __future__ import annotations

import numpy as np

__all__ = ["mesh_operator", "OPERATORS"]

#: The operators, in the core's enumerator order (``MeshOperatorKind``).
OPERATORS = ("laplacian", "mass_lumped", "mass_consistent", "gradient")


def mesh_operator(mesh, kind="laplacian", *, sparse=None, return_report=False):
    """A P1 operator of a triangle surface as a compressed-row matrix.

    ``kind`` is one of ``"laplacian"`` (the cotangent stiffness: symmetric,
    positive semidefinite, rows summing to 0, natural boundary),
    ``"mass_lumped"`` (diagonal, a third of each triangle's area per corner),
    ``"mass_consistent"`` (``A / 12 (1 + delta_ij)`` per triangle) or
    ``"gradient"`` (``3 T x n``: row ``3 t + k`` is component ``k`` of the
    gradient on triangle ``t``, ``T`` the triangles of the ``triangle`` blocks
    in input order; the rows of a zero-area triangle are empty).

    :param sparse: ``True`` returns a ``scipy.sparse.csr_matrix`` (an
        ``ImportError`` without scipy), ``False`` the tuple ``(indptr, indices,
        values, shape)``, and ``None`` (default) the matrix when scipy imports
        and the tuple otherwise.
    :param return_report: also return ``{"num_degenerate"}``.
    :raises ValueError: on an unknown kind, a volume cell, a quad, a polygon or
        a higher-order surface cell (the message names the operation that
        fixes it), no triangle, or every triangle degenerate.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    if kind not in OPERATORS:
        raise ValueError(
            f"meshio++: mesh_operator: unknown kind '{kind}' "
            f"(expected one of {', '.join(OPERATORS)})"
        )
    try:
        from . import _core
    except ImportError:
        _core = None
    if _core is None or not hasattr(_core, "mesh_operator"):
        raise NotImplementedError(
            "meshio++: mesh_operator: this operation is C++-core only and has no "
            "pure-numpy fallback -- install a build with the compiled "
            "meshioplusplus._core extension"
        )
    res = _core.mesh_operator(mesh, OPERATORS.index(kind))
    shape = (int(res["rows"]), int(res["cols"]))
    indptr, indices, values = res["indptr"], res["indices"], res["values"]
    if sparse is None:
        try:
            import scipy.sparse  # noqa: F401

            sparse = True
        except ImportError:
            sparse = False
    if sparse:
        from scipy.sparse import csr_matrix

        matrix = csr_matrix((values, indices, indptr), shape=shape)
    else:
        matrix = (
            np.asarray(indptr),
            np.asarray(indices),
            np.asarray(values),
            shape,
        )
    report = {"num_degenerate": int(res["num_degenerate"])}
    return (matrix, report) if return_report else matrix
