from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._mfm import read as _py_read
from ._mfm import write as _py_write


def read(filename):
    """Read a Modulef Formatted Mesh (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.mfm_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh, float_fmt=".16e"):
    """Write a Modulef Formatted Mesh (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "w"):
        try:
            _core.mfm_write(str(filename), mesh, float_fmt)
            return
        except Exception:
            pass
    return _py_write(filename, mesh, float_fmt)


register_format("mfm", [".mfm"], read, {"mfm": write})

__all__ = ["read", "write"]
