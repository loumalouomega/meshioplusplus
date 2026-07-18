from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._ply import read as _py_read
from ._ply import write as _py_write


def read(filename):
    """Read a PLY file (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.ply_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh, binary=True, skin=True):
    """Write a PLY file (C++ core for real file paths, Python fallback).

    With ``skin=True`` (the default) a mesh containing supported 3D volume
    cells writes its extracted boundary skin (triangles/quads, points
    compacted); ``skin=False`` keeps the legacy behavior (volume cells are
    skipped with a warning).
    """
    if not is_buffer(filename, "w"):
        try:
            _core.ply_write(str(filename), mesh, binary, skin)
            return
        except Exception:
            pass
    return _py_write(filename, mesh, binary=binary, skin=skin)


register_format("ply", [".ply"], read, {"ply": write})

__all__ = ["read", "write"]
