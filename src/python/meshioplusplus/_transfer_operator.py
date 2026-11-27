"""Inter-mesh sparse CSR transfer operators (prolongation / restriction; roadmap 15.5).

C++-core only. Builds the coarse-to-fine prolongation matrix ``P`` (or its
transpose restriction ``R = P^T``) from exact nested ``refine:entity`` stencils
or geometric barycentric location.
"""

from __future__ import annotations

import numpy as np

__all__ = ["transfer_operator", "TRANSFER_MODES"]

TRANSFER_MODES = ("auto", "nested", "barycentric")


def transfer_operator(
    coarse,
    fine,
    *,
    mode: str = "auto",
    transpose: bool = False,
    tolerance: float = 1e-12,
    fallback_nearest: bool = False,
    sparse=None,
    return_report: bool = False,
):
    """Build sparse CSR prolongation ``P`` (coarse -> fine) or restriction ``R = P^T``.

    :param coarse: coarse source mesh.
    :param fine: fine target mesh.
    :param mode: ``"auto"``, ``"nested"``, or ``"barycentric"``.
    :param transpose: if True, return restriction ``R = P^T`` (``n_coarse x n_fine``).
    :param tolerance: geometric containment tolerance.
    :param fallback_nearest: map unlocated fine points to their nearest coarse vertex.
    :param sparse: ``True`` for ``scipy.sparse.csr_matrix``, ``False`` for
        ``(indptr, indices, values, shape)``, ``None`` auto-selects scipy if installed.
    :param return_report: also return ``{"num_unlocated"}``.
    """
    if mode not in TRANSFER_MODES:
        raise ValueError(
            f"meshio++: transfer_operator: unknown mode '{mode}' "
            f"(expected one of {', '.join(TRANSFER_MODES)})"
        )
    try:
        from . import _core
    except ImportError:
        _core = None
    if _core is None or not hasattr(_core, "transfer_operator"):
        raise NotImplementedError(
            "meshio++: transfer_operator: requires the compiled meshioplusplus._core extension"
        )
    res = _core.transfer_operator(
        coarse,
        fine,
        str(mode),
        bool(transpose),
        float(tolerance),
        bool(fallback_nearest),
    )
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
    report = {"num_unlocated": int(res["num_unlocated"])}
    return (matrix, report) if return_report else matrix
