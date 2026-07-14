from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._ugrid import read as _py_read
from ._ugrid import write as _py_write


def read(filename):
    """Read an AFLR UGRID file (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.ugrid_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh):
    """Write an AFLR UGRID file (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "w"):
        try:
            _core.ugrid_write(str(filename), mesh)
            return
        except Exception:
            pass
    return _py_write(filename, mesh)


register_format("ugrid", [".ugrid"], read, {"ugrid": write})

__all__ = ["read", "write"]
