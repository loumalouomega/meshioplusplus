from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._flac3d import read as _py_read
from ._flac3d import write as _py_write


def read(filename):
    """Read a FLAC3D .f3grid file (C++ core for the common path, Python fallback)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.flac3d_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh, float_fmt: str = ".16e", binary: bool = False):
    """Write a FLAC3D .f3grid file (C++ core for the common path, Python fallback)."""
    # The C++ writer covers points + zone/face cells; cell groups (cell_sets)
    # are left to the reference Python writer.
    if not mesh.cell_sets and not is_buffer(filename, "w"):
        try:
            _core.flac3d_write(str(filename), mesh, float_fmt, binary)
            return
        except Exception:
            pass
    return _py_write(filename, mesh, float_fmt, binary)


register_format("flac3d", [".f3grid"], read, {"flac3d": write})

__all__ = ["read", "write"]
