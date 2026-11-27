"""Shortest geodesic paths and multi-source geodesic Voronoi regions (roadmap 9.14).

C++-core only.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "shortest_path",
    "geodesic_voronoi",
    "PATH_DISTANCE_NAME",
    "PATH_SOURCE_CELL_NAME",
    "GEODESIC_NEAREST_SOURCE_NAME",
]

PATH_DISTANCE_NAME = "path:distance"
PATH_SOURCE_CELL_NAME = "path:source_cell"
GEODESIC_NEAREST_SOURCE_NAME = "geodesic:nearest_source"


def shortest_path(
    mesh,
    *,
    source: int = -1,
    target: int = -1,
    source_point=None,
    target_point=None,
    method: str = "dijkstra",
    time_factor: float = 1.0,
    return_report: bool = False,
):
    """Extract a shortest path between two landmarks as a 1-D ``line`` mesh."""
    try:
        from . import _core
    except ImportError:
        _core = None
    if _core is None or not hasattr(_core, "shortest_path"):
        raise NotImplementedError(
            "meshio++: shortest_path: requires the compiled meshioplusplus._core extension"
        )
    sp = (
        np.asarray(source_point, dtype=float).ravel().tolist()
        if source_point is not None
        else []
    )
    tp = (
        np.asarray(target_point, dtype=float).ravel().tolist()
        if target_point is not None
        else []
    )
    res = _core.shortest_path(
        mesh,
        int(source),
        int(target),
        sp,
        tp,
        str(method),
        float(time_factor),
    )
    report = {
        "total_length": float(res["total_length"]),
        "num_segments": int(res["num_segments"]),
        "source_vertex": int(res["source_vertex"]),
        "target_vertex": int(res["target_vertex"]),
        "reached": bool(res["reached"]),
    }
    return (res["mesh"], report) if return_report else res["mesh"]


def geodesic_voronoi(
    mesh,
    *,
    sources=None,
    seed_points=None,
    attach_cell_labels: bool = True,
    return_report: bool = False,
):
    """Partition a mesh into multi-source geodesic Voronoi regions."""
    try:
        from . import _core
    except ImportError:
        _core = None
    if _core is None or not hasattr(_core, "geodesic_voronoi"):
        raise NotImplementedError(
            "meshio++: geodesic_voronoi: requires the compiled meshioplusplus._core extension"
        )
    src_list = [int(s) for s in sources] if sources is not None else []
    sp_flat = (
        np.asarray(seed_points, dtype=float).ravel().tolist()
        if seed_points is not None
        else []
    )
    res = _core.geodesic_voronoi(
        mesh,
        src_list,
        sp_flat,
        bool(attach_cell_labels),
    )
    report = {
        "num_sources": int(res["num_sources"]),
        "num_unreached": int(res["num_unreached"]),
        "max_distance": float(res["max_distance"]),
    }
    return (res["mesh"], report) if return_report else res["mesh"]
