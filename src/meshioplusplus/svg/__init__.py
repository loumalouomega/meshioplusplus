from typing import Union

from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._svg import write as _py_write


def write(
    filename,
    mesh,
    float_fmt: str = ".3f",
    stroke_width: Union[str, None] = None,
    image_width: Union[int, float, None] = 100,
    fill: str = "#c8c5bd",
    stroke: str = "#000080",
):
    """Write an SVG (C++ core for real file paths, Python fallback)."""
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
    )


register_format("svg", [".svg"], None, {"svg": write})

__all__ = ["write"]
