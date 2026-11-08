from .. import _core  # noqa: F401 (the shims look it up at call time)
from .._helpers import register_format
from .._shims import plain_shims
from ._wkt import read as _py_read  # noqa: F401
from ._wkt import write as _py_write  # noqa: F401

read, write = plain_shims(
    __name__,
    "wkt",
    """Read a WKT TIN file (C++ core for real file paths, Python fallback).""",
    """Write a WKT TIN file (C++ core for real file paths, Python fallback).""",
)


register_format("wkt", [".wkt"], read, {"wkt": write})

__all__ = ["read", "write"]
