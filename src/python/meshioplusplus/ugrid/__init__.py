from .. import _core  # noqa: F401 (the shims look it up at call time)
from .._helpers import register_format
from .._shims import plain_shims
from ._ugrid import read as _py_read  # noqa: F401
from ._ugrid import write as _py_write  # noqa: F401

read, write = plain_shims(
    __name__,
    "ugrid",
    """Read an AFLR UGRID file (C++ core for real file paths, Python fallback).""",
    """Write an AFLR UGRID file (C++ core for real file paths, Python fallback).""",
)


register_format("ugrid", [".ugrid"], read, {"ugrid": write})

__all__ = ["read", "write"]
