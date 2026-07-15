from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._obj import read as _py_read
from ._obj import write as _py_write


def _cpp_writable(mesh):
    return all(c.type in ("triangle", "quad", "polygon") for c in mesh.cells)


def read(filename):
    """Read an OBJ file (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.obj_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh):
    """Write an OBJ file (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "w") and _cpp_writable(mesh):
        try:
            _core.obj_write(str(filename), mesh)
            return
        except Exception:
            pass
    return _py_write(filename, mesh)


register_format("obj", [".obj"], read, {"obj": write})

__all__ = ["read", "write"]
