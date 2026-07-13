from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._freefem import read as _py_read
from ._freefem import write as _py_write


def read(filename):
    """Read a FreeFem++ .msh file (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.freefem_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh):
    """Write a FreeFem++ .msh file (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "w"):
        try:
            _core.freefem_write(str(filename), mesh)
            return
        except Exception:
            pass
    return _py_write(filename, mesh)


# `.msh` is shared with ansys and gmsh; on auto-detection those are tried first
# and freefem is attempted last. Pass file_format="freefem" to be explicit.
register_format("freefem", [".msh"], read, {"freefem": write})

__all__ = ["read", "write"]
