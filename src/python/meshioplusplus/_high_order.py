"""High-order cell evaluation, Lagrange <-> Bezier change of basis, and critical points (roadmap 15.3).

C++-core only.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "evaluate_cells",
    "change_basis",
    "critical_points",
    "CRITICAL_VALUE_NAME",
    "CRITICAL_TYPE_NAME",
    "CRITICAL_CELL_NAME",
]

CRITICAL_VALUE_NAME = "critical:value"
CRITICAL_TYPE_NAME = "critical:type"
CRITICAL_CELL_NAME = "critical:cell"


def evaluate_cells(
    mesh,
    cells,
    reference_points,
    *,
    fields=None,
    derivatives: bool = False,
):
    """Evaluate physical positions, fields, and optional physical gradients at reference points."""
    try:
        from . import _core
    except ImportError:
        _core = None
    if _core is None or not hasattr(_core, "evaluate_cells"):
        raise NotImplementedError(
            "meshio++: evaluate_cells: requires the compiled meshioplusplus._core extension"
        )
    cell_ids = [int(c) for c in np.asarray(cells, dtype=np.int64).ravel()]
    ref = np.asarray(reference_points, dtype=float)
    if ref.ndim == 2 and ref.shape[1] == 2:
        ref = np.column_stack([ref, np.zeros(ref.shape[0], dtype=float)])
    ref_flat = ref.ravel().tolist()
    field_list = [str(f) for f in fields] if fields is not None else []
    return _core.evaluate_cells(mesh, cell_ids, ref_flat, field_list, bool(derivatives))


def change_basis(
    mesh,
    target: str = "bezier",
    *,
    transform_fields: bool = True,
    return_report: bool = False,
):
    """Convert between Lagrange and Bezier cell representations via exact transformation matrices."""
    try:
        from . import _core
    except ImportError:
        _core = None
    if _core is None or not hasattr(_core, "change_basis"):
        raise NotImplementedError(
            "meshio++: change_basis: requires the compiled meshioplusplus._core extension"
        )
    res = _core.change_basis(mesh, str(target), bool(transform_fields))
    report = {"num_converted": int(res["num_converted"])}
    return (res["mesh"], report) if return_report else res["mesh"]


def critical_points(
    mesh,
    array: str,
    *,
    tolerance: float = 1e-10,
    max_iterations: int = 30,
    merge_tolerance: float = 1e-6,
    seed_subdivisions: int = 0,
    return_report: bool = False,
):
    """Locate and classify critical points (minima, maxima, saddles) of a scalar point_data field."""
    try:
        from . import _core
    except ImportError:
        _core = None
    if _core is None or not hasattr(_core, "critical_points"):
        raise NotImplementedError(
            "meshio++: critical_points: requires the compiled meshioplusplus._core extension"
        )
    res = _core.critical_points(
        mesh,
        str(array),
        float(tolerance),
        int(max_iterations),
        float(merge_tolerance),
        int(seed_subdivisions),
    )
    report = {
        "num_minima": int(res["num_minima"]),
        "num_maxima": int(res["num_maxima"]),
        "num_saddles": int(res["num_saddles"]),
        "num_degenerate": int(res["num_degenerate"]),
    }
    return (res["mesh"], report) if return_report else res["mesh"]
