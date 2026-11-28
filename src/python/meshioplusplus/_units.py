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


#: Prefix for per-array physical units in ``field_data``: ``meshio:unit:<array>``.
FIELD_UNIT_PREFIX = "meshio:unit:"

#: ``field_data`` key for a mesh's Coordinate Reference System (WKT2 or ``EPSG:<code>``).
CRS_KEY = "meshio:crs"


#: The unit names the core's ``field_unit_from_name`` knows, in its order
#: (``src/cpp/src/operations/units.cpp``; ``tests/python/test_units_table.py``
#: holds the two to each other): ``(name, scale_to_SI, (m, kg, s, K, A, mol, cd))``.
_UNIT_TABLE = (
    ("1", 1.0, (0, 0, 0, 0, 0, 0, 0)),
    ("dimensionless", 1.0, (0, 0, 0, 0, 0, 0, 0)),
    ("m", 1.0, (1, 0, 0, 0, 0, 0, 0)),
    ("mm", 1e-3, (1, 0, 0, 0, 0, 0, 0)),
    ("cm", 1e-2, (1, 0, 0, 0, 0, 0, 0)),
    ("km", 1e3, (1, 0, 0, 0, 0, 0, 0)),
    ("in", 0.0254, (1, 0, 0, 0, 0, 0, 0)),
    ("ft", 0.3048, (1, 0, 0, 0, 0, 0, 0)),
    ("kg", 1.0, (0, 1, 0, 0, 0, 0, 0)),
    ("g", 1e-3, (0, 1, 0, 0, 0, 0, 0)),
    ("t", 1e3, (0, 1, 0, 0, 0, 0, 0)),
    ("s", 1.0, (0, 0, 1, 0, 0, 0, 0)),
    ("ms", 1e-3, (0, 0, 1, 0, 0, 0, 0)),
    ("K", 1.0, (0, 0, 0, 1, 0, 0, 0)),
    ("A", 1.0, (0, 0, 0, 0, 1, 0, 0)),
    ("Pa", 1.0, (-1, 1, -2, 0, 0, 0, 0)),
    ("kPa", 1e3, (-1, 1, -2, 0, 0, 0, 0)),
    ("MPa", 1e6, (-1, 1, -2, 0, 0, 0, 0)),
    ("GPa", 1e9, (-1, 1, -2, 0, 0, 0, 0)),
    ("bar", 1e5, (-1, 1, -2, 0, 0, 0, 0)),
    ("psi", 6894.757293168361, (-1, 1, -2, 0, 0, 0, 0)),
    ("N", 1.0, (1, 1, -2, 0, 0, 0, 0)),
    ("kN", 1e3, (1, 1, -2, 0, 0, 0, 0)),
    ("m/s", 1.0, (1, 0, -1, 0, 0, 0, 0)),
    ("mm/s", 1e-3, (1, 0, -1, 0, 0, 0, 0)),
    ("m/s^2", 1.0, (1, 0, -2, 0, 0, 0, 0)),
    ("mm/s^2", 1e-3, (1, 0, -2, 0, 0, 0, 0)),
    ("kg/m^3", 1.0, (-3, 1, 0, 0, 0, 0, 0)),
    ("t/mm^3", 1e12, (-3, 1, 0, 0, 0, 0, 0)),
    ("J", 1.0, (2, 1, -2, 0, 0, 0, 0)),
    ("W", 1.0, (2, 1, -3, 0, 0, 0, 0)),
)


def _unit_text_key(text):
    """``text`` trimmed of blanks and NULs, a digit after a letter read as an
    exponent (``"kg/m3"`` -> ``"kg/m^3"``)."""
    text = str(text).strip(" \t\x00")
    out = []
    for c in text:
        if c.isdigit() and out and out[-1].isalpha():
            out.append("^")
        out.append(c)
    return "".join(out)


def field_unit_from_text(text):
    """The 8-float unit a file's unit string names, or ``None``: the core's
    ``field_unit_from_text`` (exact name, then ignoring case)."""
    key = _unit_text_key(text)
    if not key:
        return None
    for name, scale, exps in _UNIT_TABLE:
        if key == name:
            return np.asarray([scale, *exps], dtype=np.float64)
    low = key.lower()
    for name, scale, exps in _UNIT_TABLE:
        if low == name.lower():
            return np.asarray([scale, *exps], dtype=np.float64)
    return None


def field_unit_name(unit):
    """The table name of an 8-float unit (or a ``field_unit`` pair), or ``""``:
    the core's ``field_unit_name``."""
    if unit is None:
        return ""
    if isinstance(unit, tuple) and len(unit) == 2:
        flat = [float(unit[0]), *[float(x) for x in unit[1]]]
    else:
        flat = [float(x) for x in np.asarray(unit, dtype=np.float64).ravel()]
    if len(flat) != 8:
        return ""
    for name, scale, exps in _UNIT_TABLE:
        if abs(scale - flat[0]) <= 1e-12 * max(scale, flat[0]) and tuple(
            float(e) for e in exps
        ) == tuple(flat[1:]):
            return name
    return ""


def field_unit(mesh, name: str):
    """Return ``(scale_to_SI, (m, kg, s, K, A, mol, cd))`` for ``name``, or ``None``."""
    fd = getattr(mesh, "field_data", None) or {}
    raw = fd.get(f"{FIELD_UNIT_PREFIX}{name}")
    if raw is None:
        return None
    flat = np.asarray(raw, dtype=np.float64).ravel()
    if flat.size != 8 or not np.all(np.isfinite(flat)) or flat[0] <= 0.0:
        return None
    return float(flat[0]), tuple(float(x) for x in flat[1:])


def set_field_unit(mesh, name: str, unit) -> None:
    """Attach (or, with ``None``, remove) a physical unit for array ``name``.

    ``unit`` may be a unit string (e.g. ``"Pa"``, ``"MPa"``, ``"m/s"``, ``"kg/m^3"``),
    a ``(scale_to_SI, exponents_7)`` pair, or an 8-element float sequence.
    """
    key = f"{FIELD_UNIT_PREFIX}{name}"
    if unit is None:
        mesh.field_data.pop(key, None)
        return
    if isinstance(unit, str):
        from . import _core

        arr = np.asarray(_core.field_unit_from_name(unit), dtype=np.float64)
    elif isinstance(unit, tuple) and len(unit) == 2 and not np.isscalar(unit[1]):
        scale, exps = unit
        arr = np.asarray([float(scale), *[float(x) for x in exps]], dtype=np.float64)
    else:
        arr = np.asarray(unit, dtype=np.float64).ravel()
    if arr.size != 8 or not np.all(np.isfinite(arr)) or arr[0] <= 0.0:
        raise ValueError(
            f"meshio++: set_field_unit('{name}'): expected a unit name or 8 finite floats "
            "with positive SI scale"
        )
    mesh.field_data[key] = arr


def convert_units(
    mesh,
    system: str = "SI",
    *,
    strict: bool = True,
    return_report: bool = False,
):
    """Convert mesh coordinates and declared physical fields to a consistent unit system.

    Supported systems: ``"SI"``, ``"mm-t-s"``, ``"mm-kg-ms"``, ``"cm-g-s"``, ``"mm-g-ms"``.
    When ``strict=True`` (default), any ``point_data`` or ``cell_data`` array lacking a
    ``meshio:unit:<name>`` entry raises :class:`ValueError`; when ``strict=False``,
    undeclared arrays are left untouched and recorded in the report's ``warnings``.
    """
    from . import _core

    res = _core.convert_units(mesh, str(system), bool(strict))
    out = res.pop("mesh")
    return (out, res) if return_report else out


def mesh_crs(mesh) -> str | None:
    """Return the CRS string (WKT2 or ``EPSG:<code>``) stored in ``field_data["meshio:crs"]``."""
    fd = getattr(mesh, "field_data", None) or {}
    raw = fd.get(CRS_KEY)
    if raw is None:
        return None
    if isinstance(raw, str):
        return raw or None
    arr = np.asarray(raw, dtype=np.uint8).ravel()
    if arr.size == 0:
        return None
    text = bytes(arr.tolist()).rstrip(b"\x00").decode("utf-8", errors="strict")
    return text or None


def set_mesh_crs(mesh, crs: str | None) -> None:
    """Record (or, with ``None``/``""``, remove) a CRS string on ``mesh``."""
    if not crs:
        mesh.field_data.pop(CRS_KEY, None)
        return
    data = str(crs).encode("utf-8")
    mesh.field_data[CRS_KEY] = np.frombuffer(data, dtype=np.uint8).copy()
