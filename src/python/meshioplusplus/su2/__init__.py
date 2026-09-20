from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._su2 import read as _py_read
from ._su2 import write as _py_write


def read(filename):
    """Read an SU2 mesh file (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.su2_read(str(filename))
        except Exception as exc:
            if not core_declined(exc, "su2", "read", filename):
                raise
    return _py_read(filename)


def write(filename, mesh):
    """Write an SU2 mesh file (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "w"):
        try:
            _core.su2_write(str(filename), mesh)
            return
        except Exception as exc:
            if not core_declined(exc, "su2", "write", filename):
                raise
    return _py_write(filename, mesh)


register_format("su2", [".su2"], read, {"su2": write})

__all__ = ["read", "write"]
