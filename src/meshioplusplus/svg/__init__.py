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
):
    """Write an SVG (C++ core for real file paths, Python fallback).

    Flat (2D or all-z~0) meshes draw as before; a genuinely 3D mesh renders
    the boundary skin of its volume cells (or its surface cells as-is)
    through an orthographic camera given by ``azimuth``/``elevation``/
    ``roll`` in degrees — default the classic CAD isometric view — painted
    back-to-front.
    """
    if not is_buffer(filename, "w"):
        try:
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
    )


register_format("svg", [".svg"], None, {"svg": write})

__all__ = ["write"]
