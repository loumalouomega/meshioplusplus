from typing import Union

from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from .._projection import ISO_AZIMUTH, ISO_ELEVATION
from ._svg import write as _py_write


def write(
    filename,
    mesh,
    float_fmt: str = ".3f",
    stroke_width: Union[str, None] = None,
    image_width: Union[int, float, None] = 100,
    fill: str = "#c8c5bd",
    stroke: str = "#000080",
    azimuth: float = ISO_AZIMUTH,
    elevation: float = ISO_ELEVATION,
    roll: float = 0.0,
    color_by: Union[str, None] = None,
    component: Union[int, None] = None,
    cmap: str = "viridis",
    vmin: Union[float, None] = None,
    vmax: Union[float, None] = None,
    nan_color: str = "#808080",
    colorbar: bool = False,
):
    """Write an SVG (C++ core for real file paths, Python fallback).

    Flat (2D or all-z~0) meshes draw as before; a genuinely 3D mesh renders
    the boundary skin of its volume cells (or its surface cells as-is)
    through an orthographic camera given by ``azimuth``/``elevation``/
    ``roll`` in degrees — default the classic CAD isometric view — painted
    back-to-front.

    ``color_by`` names a ``point_data`` or ``cell_data`` array to colour the
    faces by: point data uses the mean of a face's corner values, cell data
    its owning cell's value (found through ``"surface:parent_cell"`` for a
    projected volume mesh). Multi-component arrays reduce to ``component`` or
    to their magnitude. The range is ``vmin``..``vmax``, defaulting to the
    finite range of the *drawn* faces; non-finite values draw in
    ``nan_color``. ``colorbar`` appends a gradient bar, widening the viewBox
    (and only the viewBox). With ``color_by`` unset the output is
    byte-identical to previous releases.
    """
    if not is_buffer(filename, "w"):
        try:
            # Positional, and the order is load-bearing: it must match the
            # py::arg list in bindings/python/_core.cpp exactly.
            _core.svg_write(
                str(filename),
                mesh,
                float_fmt,
                stroke_width,
                image_width,
                fill,
                stroke,
                azimuth,
                elevation,
                roll,
                color_by or "",
                component,
                cmap,
                vmin,
                vmax,
                nan_color,
                colorbar,
            )
            return
        except Exception:
            pass
    return _py_write(
        filename,
        mesh,
        float_fmt=float_fmt,
        stroke_width=stroke_width,
        image_width=image_width,
        fill=fill,
        stroke=stroke,
        azimuth=azimuth,
        elevation=elevation,
        roll=roll,
        color_by=color_by,
        component=component,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        nan_color=nan_color,
        colorbar=colorbar,
    )


register_format("svg", [".svg"], None, {"svg": write})

__all__ = ["write"]
