from .. import _core  # noqa: F401 (the shims look it up at call time)
from .._helpers import register_format
from .._shims import plain_shims
from ._avsucd import read as _py_read  # noqa: F401
from ._avsucd import write as _py_write  # noqa: F401

read, write = plain_shims(
    __name__,
    "avsucd",
    """Read an AVS-UCD file (C++ core for real file paths, Python fallback).""",
    """Write an AVS-UCD file (C++ core for real file paths, Python fallback).""",
)


register_format("avsucd", [".avs"], read, {"avsucd": write})

__all__ = ["read", "write"]
