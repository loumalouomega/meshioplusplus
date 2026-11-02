"""Length units as metadata (roadmap 7.1.4).

A mesh whose file states its length unit carries it as
``field_data["meshio:length_unit"]``: one float64, the number of **metres per
mesh unit** (``0.001`` for millimetres, ``0.0254`` for inches). The readers of
formats that state a unit fill it (UNV's dataset 164, Radioss's ``/BEGIN`` work
unit, USD's ``metersPerUnit``, CGNS's ``DimensionalUnits_t``, MED's ``UNI``), the
writers of formats that can state one honour it (the same five and STEP), and :func:`~meshioplusplus.transform`'s ``scale_units``
updates it, so a UNV file in millimetres converted to USD gets
``metersPerUnit = 0.001`` and a mesh scaled from millimetres to metres says so.
A mesh without the key has no stated unit: nothing is assumed.

The C++ core uses the same key (``detail::kLengthUnitName``) in the UNV,
Radioss, CGNS and MED readers and writers, the STEP writer, the pipeline's ``Transform`` step and the native
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


#: CGNS ``DimensionalUnits_t`` length names (SIDS), metres per unit.
CGNS_LENGTH_UNITS = {
    "Meter": 1.0,
    "Centimeter": 0.01,
    "Millimeter": 0.001,
    "Foot": 0.3048,
    "Inch": 0.0254,
}

#: MED coordinate unit names (``UNI``) recognised as lengths, metres per unit;
#: the first of a value is the one written.
MED_LENGTH_UNITS = {
    "m": 1.0,
    "cm": 0.01,
    "mm": 0.001,
    "um": 1e-6,
    "µm": 1e-6,
    "km": 1000.0,
    "in": 0.0254,
    "ft": 0.3048,
}


def unit_name(metres_per_unit, table):
    """The first name in ``table`` whose value matches ``metres_per_unit``
    (relative 1e-12), or ``None``."""
    if metres_per_unit is None:
        return None
    for name, value in table.items():
        if abs(value - metres_per_unit) <= 1e-12 * value:
            return name
    return None


def med_length_unit(unit_coords):
    """The metres per unit a MED ``UNI`` string names: the same known length
    name for every axis (written 16 characters per axis), else ``None``."""
    names = set((unit_coords or "").replace("\x00", " ").split())
    if len(names) != 1:
        return None
    return MED_LENGTH_UNITS.get(names.pop())


def med_unit_coords(mesh, dim):
    """The MED ``UNI`` string for ``mesh``'s length unit (its name padded to 16
    characters per axis), or ``""`` when it states none MED can name."""
    name = unit_name(length_unit(mesh), MED_LENGTH_UNITS)
    return "".join(f"{name:<16}" for _ in range(dim)) if name else ""
