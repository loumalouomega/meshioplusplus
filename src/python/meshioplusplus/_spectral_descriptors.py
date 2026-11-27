"""Intrinsic spectral shape descriptors (HKS and WKS; roadmap 15.6).

C++-core only, built directly on ``laplacian_eigenmodes`` (or reusing
precomputed ``laplacian:modes`` and ``laplacian:eigenvalues`` on the mesh).
"""

from __future__ import annotations

import numpy as np

__all__ = ["spectral_descriptors", "HKS_NAME", "WKS_NAME", "SCALES_NAME"]

HKS_NAME = "spectral:hks"
WKS_NAME = "spectral:wks"
SCALES_NAME = "spectral:scales"


def spectral_descriptors(
    mesh,
    kind: str = "hks",
    *,
    num_modes: int = 32,
    num_scales: int = 16,
    scale_invariant: bool = False,
    reuse_eigenmodes: bool = True,
    wks_variance: float = 0.0,
    tolerance: float = 1e-8,
    max_iterations: int = 0,
    output: str = "",
    return_report: bool = False,
):
    """Compute Heat Kernel Signature (HKS) or Wave Kernel Signature (WKS).

    :param mesh: triangle surface mesh.
    :param kind: ``"hks"`` or ``"wks"``.
    :param num_modes: eigenmodes to compute when not reusing precomputed ones.
    :param num_scales: number of log-spaced scales.
    :param scale_invariant: normalize HKS per scale by its spectral trace.
    :param reuse_eigenmodes: reuse ``laplacian:modes`` / ``laplacian:eigenvalues`` if present.
    :param wks_variance: WKS log-energy Gaussian std (0 = automatic).
    :param tolerance: eigensolver tolerance when computing eigenmodes.
    :param max_iterations: eigensolver iteration cap (0 = default).
    :param output: output point_data array name (default ``spectral:hks`` / ``spectral:wks``).
    :param return_report: also return ``{"scales", "num_modes_used", "reused_eigenmodes"}``.
    """
    try:
        from . import _core
    except ImportError:
        _core = None
    if _core is None or not hasattr(_core, "spectral_descriptors"):
        raise NotImplementedError(
            "meshio++: spectral_descriptors: requires the compiled meshioplusplus._core extension"
        )
    res = _core.spectral_descriptors(
        mesh,
        kind,
        int(num_modes),
        int(num_scales),
        bool(scale_invariant),
        bool(reuse_eigenmodes),
        float(wks_variance),
        float(tolerance),
        int(max_iterations),
        str(output),
    )
    report = {
        "scales": np.asarray(res["scales"], dtype=float),
        "num_modes_used": int(res["num_modes_used"]),
        "reused_eigenmodes": bool(res["reused_eigenmodes"]),
    }
    return (res["mesh"], report) if return_report else res["mesh"]
