"""Curve-skeleton and centerline extraction via mean-curvature flow contraction (roadmap 9.13).

C++-core only.
"""

from __future__ import annotations

__all__ = [
    "extract_skeleton",
    "SKELETON_RADIUS_NAME",
    "SKELETON_DEGREE_NAME",
    "SKELETON_SURFACE_MAP_NAME",
]

SKELETON_RADIUS_NAME = "skeleton:radius"
SKELETON_DEGREE_NAME = "skeleton:degree"
SKELETON_SURFACE_MAP_NAME = "skeleton:surface_map"


def extract_skeleton(
    mesh,
    *,
    max_iterations: int = 10,
    contraction_weight: float = 1.0,
    attraction_weight: float = 1.0,
    volume_tolerance: float = 1e-4,
    merge_tolerance: float = 0.0,
    recenter: bool = True,
    return_report: bool = False,
):
    """Extract a 1-D curve skeleton / centerline from a closed triangle surface."""
    try:
        from . import _core
    except ImportError:
        _core = None
    if _core is None or not hasattr(_core, "extract_skeleton"):
        raise NotImplementedError(
            "meshio++: extract_skeleton: requires the compiled meshioplusplus._core extension"
        )
    res = _core.extract_skeleton(
        mesh,
        int(max_iterations),
        float(contraction_weight),
        float(attraction_weight),
        float(volume_tolerance),
        float(merge_tolerance),
        bool(recenter),
    )
    report = {
        "contraction_iterations": int(res["contraction_iterations"]),
        "volume_ratio": float(res["volume_ratio"]),
        "num_nodes": int(res["num_nodes"]),
        "num_segments": int(res["num_segments"]),
        "num_branches": int(res["num_branches"]),
        "num_endpoints": int(res["num_endpoints"]),
    }
    return (res["mesh"], report) if return_report else res["mesh"]
