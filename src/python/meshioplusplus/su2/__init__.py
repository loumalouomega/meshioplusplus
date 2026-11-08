from .. import _core  # noqa: F401 (the shims look it up at call time)
from .._helpers import register_format
from .._shims import plain_shims
from ._su2 import read as _py_read  # noqa: F401
from ._su2 import write as _py_write  # noqa: F401

read, write = plain_shims(
    __name__,
    "su2",
    """Read an SU2 mesh file (C++ core for real file paths, Python fallback).""",
    """Write an SU2 mesh file (C++ core for real file paths, Python fallback).""",
)


register_format("su2", [".su2"], read, {"su2": write})

__all__ = ["read", "write"]
