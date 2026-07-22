from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._wkt import read as _py_read
from ._wkt import write as _py_write


def read(filename):
    """Read a WKT TIN file (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.wkt_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh):
    """Write a WKT TIN file (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "w"):
        try:
            _core.wkt_write(str(filename), mesh)
            return
        except Exception:
            pass
    return _py_write(filename, mesh)


register_format("wkt", [".wkt"], read, {"wkt": write})

__all__ = ["read", "write"]
