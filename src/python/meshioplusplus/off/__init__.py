from .. import _core  # noqa: F401 (the shims look it up at call time)
from .._helpers import register_format
from .._shims import plain_shims
from ._off import read as _py_read  # noqa: F401
from ._off import write as _py_write  # noqa: F401

read, write = plain_shims(
    __name__,
    "off",
    """Read an OFF file (C++ core for real file paths, Python fallback).""",
    """Write an OFF file (C++ core for real file paths, Python fallback).""",
)


register_format("off", [".off"], read, {"off": write})

__all__ = ["read", "write"]
