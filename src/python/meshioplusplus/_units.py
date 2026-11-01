"""Length units as metadata (roadmap 7.1.4).

A mesh whose file states its length unit carries it as
``field_data["meshio:length_unit"]``: one float64, the number of **metres per
mesh unit** (``0.001`` for millimetres, ``0.0254`` for inches). The readers of
formats that state a unit fill it (UNV's dataset 164, Radioss's ``/BEGIN`` work
unit, USD's ``metersPerUnit``), the writers of formats that can state one honour
it (UNV, Radioss, USD), and :func:`~meshioplusplus.transform`'s ``scale_units``
updates it, so a UNV file in millimetres converted to USD gets
``metersPerUnit = 0.001`` and a mesh scaled from millimetres to metres says so.
A mesh without the key has no stated unit: nothing is assumed.

The C++ core uses the same key (``detail::kLengthUnitName``) in the UNV and
Radioss readers and writers, the pipeline's ``Transform`` step and the native
CLI's ``transform --scale-units``.
"""

from __future__ import annotations

import numpy as np

LENGTH_UNIT_KEY = "meshio:length_unit"


def length_unit(mesh):
    """The metres per unit ``mesh`` states, or ``None`` when it states none."""
    fd = getattr(mesh, "field_data", None) or {}
    value = fd.get(LENGTH_UNIT_KEY)
    if value is None:
        return None
    flat = np.asarray(value, dtype=float).ravel()
    if flat.size != 1 or not np.isfinite(flat[0]) or flat[0] <= 0.0:
        return None
    return float(flat[0])


def set_length_unit(mesh, metres_per_unit) -> None:
    """Record (or, with ``None``, remove) the length unit ``mesh`` is in."""
    if metres_per_unit is None:
        mesh.field_data.pop(LENGTH_UNIT_KEY, None)
        return
    value = float(metres_per_unit)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(
            f"meshio++: a length unit is a positive number of metres, not {value!r}"
        )
    mesh.field_data[LENGTH_UNIT_KEY] = np.array([value])
