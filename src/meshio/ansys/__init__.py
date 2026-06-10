from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._ansys import read as _py_read
from ._ansys import write as _py_write


def read(filename):
    """Read an Ansys/Fluent .msh file (C++ core, Python fallback)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.ansys_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh, binary=True):
    """Write an Ansys/Fluent .msh file (C++ core, Python fallback)."""
    if not is_buffer(filename, "w"):
        try:
            _core.ansys_write(str(filename), mesh, binary)
            return
        except Exception:
            pass
    return _py_write(filename, mesh, binary)


register_format("ansys", [".msh"], read, {"ansys": write})

__all__ = ["read", "write"]
