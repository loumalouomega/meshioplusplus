from typing import Union

from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from .._projection import ISO_AZIMUTH, ISO_ELEVATION
from ._tikz import read as _py_read
from ._tikz import write as _py_write


def write(
    filename,
    mesh,
    float_fmt: str = ".6f",
    standalone: bool = True,
    line_width: Union[str, None] = None,
    fill: str = "gray!30",
    draw: str = "black",
    scale: Union[int, float, None] = None,
    azimuth: float = ISO_AZIMUTH,
    elevation: float = ISO_ELEVATION,
    roll: float = 0.0,
    color_by: Union[str, None] = None,
    component: Union[int, None] = None,
    cmap: str = "viridis",
    vmin: Union[float, None] = None,
    vmax: Union[float, None] = None,
    nan_color: str = "gray",
    colorbar: bool = False,
):
    """Write a TikZ figure (C++ core for real file paths, Python fallback).

    Flat (2D or all-z~0) meshes draw as before; a genuinely 3D mesh renders
    the boundary skin of its volume cells (or its surface cells as-is)
    through an orthographic camera given by ``azimuth``/``elevation``/
    ``roll`` in degrees — default the classic CAD isometric view — drawn
    back-to-front.

    ``color_by`` names a ``point_data`` or ``cell_data`` array to colour the
    faces by: point data uses the mean of a face's corner values, cell data
    its owning cell's value (found through ``"surface:parent_cell"`` for a
    projected volume mesh). Multi-component arrays reduce to ``component`` or
    to their magnitude. The range is ``vmin``..``vmax``, defaulting to the
    finite range of the *drawn* faces; non-finite values draw in
    ``nan_color``. With ``color_by`` unset the output is byte-identical to
    previous releases.
    """
    if not is_buffer(filename, "w"):
        try:
            # Positional, and the order is load-bearing: it must match the
            # py::arg list in bindings/python/_core.cpp exactly.
            _core.tikz_write(
                str(filename),
                mesh,
                float_fmt,
                standalone,
                line_width,
                fill,
                draw,
                scale,
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
        except Exception as exc:
            if not core_declined(exc, "tikz", "write", filename):
                raise
    return _py_write(
        filename,
        mesh,
        float_fmt=float_fmt,
        standalone=standalone,
        line_width=line_width,
        fill=fill,
        draw=draw,
        scale=scale,
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


def read(filename, points_only=False, arrays=None, strict=True):
    """Read a TikZ drawing written by :func:`write` back as a planar mesh.

    A drawing is not a mesh, so only the structure meshio++'s own writer emits is
    read (C++ core for real file paths, Python twin otherwise); anything else is a
    ``ReadError`` naming the construct. A file without meshio++'s marker comment is
    refused unless ``strict=False``. ``z`` is always 0, and the read warns once that
    depth is lost. See ``doc/formats/tikz.md``.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.tikz_read(str(filename), points_only, arrays, strict)
        except Exception as exc:
            if not core_declined(exc, "tikz", "read", filename):
                raise
    return _py_read(filename, points_only=points_only, arrays=arrays, strict=strict)


register_format("tikz", [".tikz"], read, {"tikz": write})

__all__ = ["read", "write"]
