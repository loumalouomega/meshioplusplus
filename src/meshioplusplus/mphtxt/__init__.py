from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._mphtxt import read as _py_read
from ._mphtxt import write as _py_write


def read(filename):
    """Read a COMSOL .mphtxt file (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.mphtxt_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh):
    """Write a COMSOL .mphtxt file (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "w"):
        try:
            _core.mphtxt_write(str(filename), mesh)
            return
        except Exception:
            pass
    return _py_write(filename, mesh)


register_format("mphtxt", [".mphtxt"], read, {"mphtxt": write})

__all__ = ["read", "write"]
