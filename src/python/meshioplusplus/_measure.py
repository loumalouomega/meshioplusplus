"""Elements, dimensions and GD&T tolerances measured from a scan (roadmap
7.2), and the 3-2-1 datum alignment on them (roadmap 7.3.4).

C++-core only (see ``doc/measure.md``); ``measure_table`` reads a mesh's
element table in Python.
"""

from __future__ import annotations

import json
import os

import numpy as np

__all__ = ["measure", "measure_table", "align_datums"]

_PREFIX = "meshio++: measure: "
TABLE_KEY = "measure:table"


def _spec_text(spec) -> str:
    if isinstance(spec, dict):
        return json.dumps(spec)
    if isinstance(spec, (str, os.PathLike)):
        text = os.fspath(spec)
        if text.lstrip().startswith("{"):
            return text
        with open(text, encoding="utf-8") as f:
            return f.read()
    raise TypeError(
        f"{_PREFIX}spec is a dict, a JSON string or the path of a JSON file"
    )


def measure(mesh, spec):
    """Measure ``mesh`` by ``spec``: elements, dimensions and tolerances.

    :param mesh: the scanned surface (its cell regions, or seeds, pick the
        elements to fit).
    :param spec: a dict, its JSON text, or the path of a JSON file: ``elements``
        (fitted with ``"fit"`` on a ``"region"``, from ``"seeds"`` or a
        ``"patch"``; constructed with ``"construct"`` and ``"of"``; or a
        literal ``"point"``), ``dimensions`` (``distance``, ``angle``, with
        ``nominal`` and ``tolerance``), ``tolerances`` (flatness,
        cylindricity, sphericity, circularity, parallelism, perpendicularity,
        angularity, coaxiality, concentricity, with a ``limit``) and the form
        ``criterion``. See ``doc/measure.md``.
    :returns: ``(mesh, report)``: the input with the element table
        (``measure:table`` field data) and each fitted element's cells as the
        region ``measure:<name>``; the report has ``elements``,
        ``dimensions``, ``tolerances`` (each with ``value`` and ``pass``),
        ``num_fail`` and ``warnings``.
    :raises ValueError: on a malformed spec or an element that cannot be fitted
        or constructed.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _core_op

    res = _core_op("measure")(mesh, _spec_text(spec))
    return res["mesh"], json.loads(res["report"])


def measure_table(mesh):
    """The element table ``measure`` stored on ``mesh`` (a list of dicts:
    ``name``, ``kind``, ``point``, ``dir``, ``radius``, ``angle``,
    ``form_error``, ``measured``, ``source``, ``inputs``); empty without one."""
    raw = getattr(mesh, "field_data", {}).get(TABLE_KEY)
    if raw is None:
        return []
    text = np.asarray(raw, dtype=np.uint8).tobytes().decode("utf-8")
    return json.loads(text)["elements"]


def align_datums(mesh, primary: str, secondary: str, origin: str | None = None):
    """The 3-2-1 datum alignment of a measured mesh (roadmap 7.3.4).

    :param mesh: a mesh carrying an element table (``measure``).
    :param primary: the element (a plane, or an axis) whose normal or
        direction becomes Z.
    :param secondary: the element whose normal or direction, projected, becomes X.
    :param origin: a point-like element, or a third plane (the three planes'
        meeting point); ``None`` takes the primary's point.
    :returns: ``(mesh, matrix)``: the mesh moved, its element table with it,
        and the row-major 4x4 matrix ``transform`` applies.
    :raises ValueError: on an unknown element, a secondary parallel to the
        primary, or planes that do not meet in a point.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _core_op

    res = _core_op("align_datums")(mesh, str(primary), str(secondary), origin or "")
    return res["mesh"], np.asarray(res["matrix"], dtype=np.float64).reshape(4, 4)
