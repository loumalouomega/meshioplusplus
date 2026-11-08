from .. import _core  # noqa: F401 (the shims look it up at call time)
from .._helpers import register_format
from .._shims import plain_shims
from ._ip import read as _py_read  # noqa: F401
from ._ip import write as _py_write  # noqa: F401

read, write = plain_shims(
    __name__,
    "ip",
    """Read an ANSYS Fluent interpolation file (C++ core, Python fallback).""",
    """Write an ANSYS Fluent (version 3) interpolation file (C++ core, Python fallback).""",
)


register_format("ip", [".ip"], read, {"ip": write})

__all__ = ["read", "write"]
