"""Sparse inter-mesh transfer operators for hierarchical GNN pooling / multigrid (roadmap 15.5.2)."""

from __future__ import annotations

import numpy as np

from .._transfer_operator import transfer_operator

__all__ = ["transfer_matrix"]


def transfer_matrix(
    coarse,
    fine,
    *,
    mode: str = "auto",
    transpose: bool = False,
    tolerance: float = 1e-12,
    fallback_nearest: bool = False,
    format: str = "scipy",
):
    """Return a sparse inter-mesh prolongation (or restriction) operator for GNN pooling.

    :param format: ``"scipy"`` (``scipy.sparse.csr_matrix``), ``"torch"``
        (``torch.sparse_coo_tensor``), or ``"csr"`` (raw ``(indptr, indices, values, shape)``).
    """
    if format == "csr":
        return transfer_operator(
            coarse,
            fine,
            mode=mode,
            transpose=transpose,
            tolerance=tolerance,
            fallback_nearest=fallback_nearest,
            sparse=False,
        )
    if format == "scipy":
        return transfer_operator(
            coarse,
            fine,
            mode=mode,
            transpose=transpose,
            tolerance=tolerance,
            fallback_nearest=fallback_nearest,
            sparse=True,
        )
    if format == "torch":
        import torch

        indptr, indices, values, shape = transfer_operator(
            coarse,
            fine,
            mode=mode,
            transpose=transpose,
            tolerance=tolerance,
            fallback_nearest=fallback_nearest,
            sparse=False,
        )
        counts = np.diff(indptr)
        rows = np.repeat(np.arange(shape[0], dtype=np.int64), counts)
        idx = torch.as_tensor(np.vstack([rows, indices]), dtype=torch.long)
        val = torch.as_tensor(values, dtype=torch.float32)
        return torch.sparse_coo_tensor(idx, val, size=shape).coalesce()
    raise ValueError(
        f"meshio++: physicsnemo.transfer_matrix: unknown format '{format}' "
        "(expected 'scipy', 'torch', or 'csr')"
    )
