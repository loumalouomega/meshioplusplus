"""Feature-preserving surface denoising (bilateral and guided normal filtering; roadmap 9.12).

C++-core only.
"""

from __future__ import annotations

__all__ = ["denoise_surface", "DENOISE_METHODS"]

DENOISE_METHODS = ("bilateral", "guided")


def denoise_surface(
    mesh,
    *,
    method: str = "bilateral",
    normal_iterations: int = 10,
    vertex_iterations: int = 10,
    sigma_spatial: float = 0.0,
    sigma_normal: float = 0.35,
    fix_boundary: bool = True,
    guard_foldover: bool = True,
    region: str = "",
    return_report: bool = False,
):
    """Denoise a surface mesh while preserving sharp dihedral features."""
    if method not in DENOISE_METHODS:
        raise ValueError(
            f"meshio++: denoise_surface: unknown method '{method}' "
            f"(expected one of {', '.join(DENOISE_METHODS)})"
        )
    try:
        from . import _core
    except ImportError:
        _core = None
    if _core is None or not hasattr(_core, "denoise_surface"):
        raise NotImplementedError(
            "meshio++: denoise_surface: requires the compiled meshioplusplus._core extension"
        )
    res = _core.denoise_surface(
        mesh,
        str(method),
        int(normal_iterations),
        int(vertex_iterations),
        float(sigma_spatial),
        float(sigma_normal),
        bool(fix_boundary),
        bool(guard_foldover),
        str(region),
    )
    report = {
        "num_nodes_moved": int(res["num_nodes_moved"]),
        "max_displacement": float(res["max_displacement"]),
        "num_skipped_foldover": int(res["num_skipped_foldover"]),
        "sigma_spatial_used": float(res["sigma_spatial_used"]),
    }
    return (res["mesh"], report) if return_report else res["mesh"]
