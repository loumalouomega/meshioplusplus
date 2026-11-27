"""Parametric FFD lattice and RBF control-point morphing (roadmap 9.15).

C++-core only.
"""

from __future__ import annotations

import numpy as np

__all__ = ["deform_ffd", "deform_rbf", "RBF_KERNELS"]

RBF_KERNELS = ("tps", "wendland_c2")


def deform_ffd(
    mesh,
    grid_dims,
    displacements,
    *,
    degrees=(-1, -1, -1),
    box_min=None,
    box_max=None,
    region: str = "",
    return_report: bool = False,
):
    """Deform ``mesh`` via a tensor-product Bezier / B-spline FFD control lattice."""
    try:
        from . import _core
    except ImportError:
        _core = None
    if _core is None or not hasattr(_core, "deform_ffd"):
        raise NotImplementedError(
            "meshio++: deform_ffd: requires the compiled meshioplusplus._core extension"
        )
    gd = [int(x) for x in grid_dims]
    deg = [int(x) for x in degrees]
    disp = np.asarray(displacements, dtype=float).ravel().tolist()
    bmin = (
        np.asarray(box_min, dtype=float).ravel().tolist() if box_min is not None else []
    )
    bmax = (
        np.asarray(box_max, dtype=float).ravel().tolist() if box_max is not None else []
    )
    res = _core.deform_ffd(mesh, gd, disp, deg, bmin, bmax, str(region))
    report = {
        "num_nodes_moved": int(res["num_nodes_moved"]),
        "num_inverted": int(res["num_inverted"]),
        "max_displacement": float(res["max_displacement"]),
    }
    return (res["mesh"], report) if return_report else res["mesh"]


def deform_rbf(
    mesh,
    control_points,
    displacements,
    *,
    kernel: str = "tps",
    radius: float = 0.0,
    polynomial: bool = True,
    regularization: float = 0.0,
    region: str = "",
    return_report: bool = False,
):
    """Deform ``mesh`` by RBF interpolation of sparse control-point displacements."""
    if kernel not in RBF_KERNELS:
        raise ValueError(
            f"meshio++: deform_rbf: unknown kernel '{kernel}' "
            f"(expected one of {', '.join(RBF_KERNELS)})"
        )
    try:
        from . import _core
    except ImportError:
        _core = None
    if _core is None or not hasattr(_core, "deform_rbf"):
        raise NotImplementedError(
            "meshio++: deform_rbf: requires the compiled meshioplusplus._core extension"
        )
    cp = np.asarray(control_points, dtype=float).ravel().tolist()
    disp = np.asarray(displacements, dtype=float).ravel().tolist()
    res = _core.deform_rbf(
        mesh,
        cp,
        disp,
        str(kernel),
        float(radius),
        bool(polynomial),
        float(regularization),
        str(region),
    )
    report = {
        "num_nodes_moved": int(res["num_nodes_moved"]),
        "num_inverted": int(res["num_inverted"]),
        "max_displacement": float(res["max_displacement"]),
    }
    return (res["mesh"], report) if return_report else res["mesh"]
