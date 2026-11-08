from .. import _core  # noqa: F401 (the shims look it up at call time)
from .._helpers import register_format
from .._shims import plain_shims
from ._freefem import read as _py_read  # noqa: F401
from ._freefem import write as _py_write  # noqa: F401

read, write = plain_shims(
    __name__,
    "freefem",
    """Read a FreeFem++ .msh file (C++ core for real file paths, Python fallback).""",
    """Write a FreeFem++ .msh file (C++ core for real file paths, Python fallback).""",
)


# `.msh` is shared with ansys and gmsh; on auto-detection those are tried first
# and freefem is attempted last. Pass file_format="freefem" to be explicit.
register_format("freefem", [".msh"], read, {"freefem": write})

__all__ = ["read", "write"]
