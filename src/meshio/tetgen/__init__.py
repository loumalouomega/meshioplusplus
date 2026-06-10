from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._tetgen import read as _py_read
from ._tetgen import write as _py_write


def read(filename):
    """Read a TetGen .node/.ele pair (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.tetgen_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh, float_fmt=".16e"):
    """Write a TetGen .node/.ele pair (C++ core for real file paths, Python fallback)."""
    if float_fmt == ".16e" and not is_buffer(filename, "w"):
        try:
            _core.tetgen_write(str(filename), mesh)
            return
        except Exception:
            pass
    return _py_write(filename, mesh, float_fmt)


register_format("tetgen", [".ele", ".node"], read, {"tetgen": write})

__all__ = ["read", "write"]
