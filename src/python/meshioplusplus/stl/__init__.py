from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from .._skin import _has_skinnable_cells
from ._stl import read as _py_read
from ._stl import write as _py_write


def _cpp_writable(mesh, skin):
    # STL is triangle-only, but with skin=True the C++ core also accepts a
    # volume mesh (it extracts and triangulates the boundary skin itself);
    # defer anything else to the Python writer (which warns and discards as
    # before, or performs the same skin extraction in numpy).
    if all(c.type == "triangle" for c in mesh.cells):
        return True
    return skin and _has_skinnable_cells(mesh)


def read(filename):
    """Read an STL file (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.stl_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh, binary=False, skin=True):
    """Write an STL file (C++ core for real file paths, Python fallback).

    With ``skin=True`` (the default) a mesh containing supported 3D volume
    cells writes its extracted boundary skin (quads triangulated);
    ``skin=False`` keeps the legacy behavior (only existing ``triangle``
    blocks are written, everything else is discarded with a warning).
    """
    if not is_buffer(filename, "w") and _cpp_writable(mesh, skin):
        try:
            _core.stl_write(str(filename), mesh, binary, skin)
            return
        except Exception:
            pass
    return _py_write(filename, mesh, binary=binary, skin=skin)


register_format("stl", [".stl"], read, {"stl": write})

__all__ = ["read", "write"]
