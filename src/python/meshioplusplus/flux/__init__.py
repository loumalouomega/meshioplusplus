from .. import _core  # noqa: F401 (the shims look it up at call time)
from .._helpers import register_format
from .._shims import plain_shims
from ._flux import read as _py_read  # noqa: F401
from ._flux import write as _py_write  # noqa: F401

read, write = plain_shims(
    __name__,
    "flux",
    """Read a FLUX .pf3 file (C++ core for real file paths, Python fallback).""",
    """Write a FLUX .pf3 file (C++ core for real file paths, Python fallback).""",
)


register_format("flux", [".pf3"], read, {"flux": write})

__all__ = ["read", "write"]
