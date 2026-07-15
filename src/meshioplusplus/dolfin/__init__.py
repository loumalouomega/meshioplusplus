from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._dolfin import read as _py_read
from ._dolfin import write as _py_write


def read(filename):
    """Read a DOLFIN XML file (C++ core, Python fallback)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.dolfin_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh):
    """Write a DOLFIN XML file (C++ core, Python fallback)."""
    if not is_buffer(filename, "w"):
        try:
            _core.dolfin_write(str(filename), mesh)
            return
        except Exception:
            pass
    return _py_write(filename, mesh)


register_format("dolfin-xml", [".xml"], read, {"dolfin-xml": write})

__all__ = ["read", "write"]
