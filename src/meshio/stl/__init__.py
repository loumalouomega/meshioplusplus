from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._stl import read as _py_read
from ._stl import write as _py_write


def _cpp_writable(mesh):
    # STL is triangle-only; defer anything else to the Python writer (which warns
    # and discards as before).
    return all(c.type == "triangle" for c in mesh.cells)


def read(filename):
    """Read an STL file (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.stl_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh, binary=False):
    """Write an STL file (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "w") and _cpp_writable(mesh):
        try:
            _core.stl_write(str(filename), mesh, binary)
            return
        except Exception:
            pass
    return _py_write(filename, mesh, binary=binary)


register_format("stl", [".stl"], read, {"stl": write})

__all__ = ["read", "write"]
