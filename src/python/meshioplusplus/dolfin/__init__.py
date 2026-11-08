from .. import _core  # noqa: F401 (the shims look it up at call time)
from .._helpers import register_format
from .._shims import plain_shims
from ._dolfin import read as _py_read  # noqa: F401
from ._dolfin import write as _py_write  # noqa: F401

read, write = plain_shims(
    __name__,
    "dolfin",
    """Read a DOLFIN XML file (C++ core, Python fallback).""",
    """Write a DOLFIN XML file (C++ core, Python fallback).""",
)


register_format("dolfin-xml", [".xml"], read, {"dolfin-xml": write})

__all__ = ["read", "write"]
