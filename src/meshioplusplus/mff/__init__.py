from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._mff import read as _py_read
from ._mff import write as _py_write


def read(filename):
    """Read a Modulef Formatted Field (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.mff_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh, float_fmt=".16e"):
    """Write a Modulef Formatted Field (C++ core for real file paths, Python fallback)."""
    if float_fmt == ".16e" and not is_buffer(filename, "w"):
        try:
            _core.mff_write(str(filename), mesh)
            return
        except Exception:
            pass
    return _py_write(filename, mesh, float_fmt)


register_format("mff", [".mff"], read, {"mff": write})

__all__ = ["read", "write"]
