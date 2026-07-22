from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._ip import read as _py_read
from ._ip import write as _py_write


def read(filename):
    """Read an ANSYS Fluent interpolation file (C++ core, Python fallback)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.ip_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh):
    """Write an ANSYS Fluent (version 3) interpolation file (C++ core, Python fallback)."""
    if not is_buffer(filename, "w"):
        try:
            _core.ip_write(str(filename), mesh)
            return
        except Exception:
            pass
    return _py_write(filename, mesh)


register_format("ip", [".ip"], read, {"ip": write})

__all__ = ["read", "write"]
