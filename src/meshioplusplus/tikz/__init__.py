from typing import Union

from .. import _core
from .._files import is_buffer
from .._helpers import register_format
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
):
    """Write a TikZ figure (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "w"):
        try:
            _core.tikz_write(
                str(filename),
                mesh,
                float_fmt,
                standalone,
                line_width,
                fill,
                draw,
                scale,
            )
            return
        except Exception:
            pass
    return _py_write(
        filename,
        mesh,
        float_fmt=float_fmt,
        standalone=standalone,
        line_width=line_width,
        fill=fill,
        draw=draw,
        scale=scale,
    )


register_format("tikz", [".tikz"], None, {"tikz": write})

__all__ = ["write"]
