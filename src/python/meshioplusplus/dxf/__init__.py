from .. import _core  # noqa: F401 (the shims look it up at call time)
from .._helpers import register_format
from .._shims import plain_shims
from ._dxf import read as _py_read  # noqa: F401
from ._dxf import write as _py_write  # noqa: F401

read, write = plain_shims(
    __name__,
    "dxf",
    """Read a DXF drawing's lines, arcs, circles, polylines and splines (C++ core for
real file paths, Python fallback).""",
    """Write a line mesh or a ``fit_curves`` sketch as ASCII DXF R2000 (C++ core for real
file paths, Python fallback).""",
)


register_format("dxf", [".dxf"], read, {"dxf": write})

__all__ = ["read", "write"]
