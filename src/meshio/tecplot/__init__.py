from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._tecplot import read as _py_read
from ._tecplot import write as _py_write


def read(filename):
    """Read a Tecplot ASCII file (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.tecplot_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh):
    """Write a Tecplot ASCII file (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "w"):
        try:
            _core.tecplot_write(str(filename), mesh)
            return
        except Exception:
            pass
    return _py_write(filename, mesh)


register_format("tecplot", [".dat", ".tec"], read, {"tecplot": write})

__all__ = ["read", "write"]
