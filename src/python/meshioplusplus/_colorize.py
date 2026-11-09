"""Bake a scalar field into 8-bit RGB colours through a colormap (roadmap 7.4.4).

The colours are written as UInt8 ``red``, ``green`` and ``blue`` arrays on the
field's location, the names PLY's vertex and face colour properties take; glTF
takes them with ``color_by="rgb"``. ``steps`` bands the map into flat colours
at fixed steps. See ``doc/colorize.md``.

The C++ core does the work; there is no numpy twin.
"""

from __future__ import annotations

__all__ = ["colorize"]

_PREFIX = "meshio++: colorize: "


def colorize(
    mesh,
    field: str,
    *,
    cmap: str = "viridis",
    vmin: float | None = None,
    vmax: float | None = None,
    symmetric: bool = False,
    steps: int = 0,
    nan_color: str = "#808080",
    component: int | None = None,
    return_report: bool = False,
):
    """Bake ``field`` (point or cell data) into ``red``/``green``/``blue`` UInt8 data.

    :param field: the array to colour by; a vector field is reduced to its
        magnitude unless ``component`` picks one entry.
    :param cmap: a colormap name (viridis, coolwarm, turbo, magma, ...).
    :param vmin, vmax: the range; ``None`` takes the field's finite min and max.
    :param symmetric: with the automatic range, ``+-max|v|`` (zero mid-map).
    :param steps: 0 or 1 continuous, 2 or more bands at fixed steps.
    :param nan_color: ``#rrggbb`` for a non-finite value.
    :param return_report: also return ``{"vmin", "vmax", "location"}``.
    :returns: the mesh with the three channels; with ``return_report``
        ``(mesh, report)``.
    """
    try:
        from . import _core
    except ImportError as exc:  # pragma: no cover - the numpy-only install
        raise NotImplementedError(
            f"{_PREFIX}needs the compiled C++ core (meshioplusplus._core); "
            "it has no numpy twin"
        ) from exc
    if (vmin is None) != (vmax is None):
        raise ValueError(f"{_PREFIX}give both vmin and vmax, or neither")
    lo, hi = (0.0, 0.0) if vmin is None else (float(vmin), float(vmax))
    if vmin is not None and not lo < hi:
        raise ValueError(f"{_PREFIX}vmin must be below vmax")
    res = _core.colorize(
        mesh,
        str(field),
        cmap=str(cmap),
        vmin=lo,
        vmax=hi,
        symmetric=bool(symmetric),
        steps=int(steps),
        nan_color=str(nan_color),
        component=-1 if component is None else int(component),
    )
    out = res["mesh"]
    if not return_report:
        return out
    return out, {
        "vmin": float(res["vmin"]),
        "vmax": float(res["vmax"]),
        "location": str(res["location"]),
    }
