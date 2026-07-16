from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._dex import read as _py_read
from ._dex import write as _py_write


def read(filename):
    """Read a FLUX field file (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.dex_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh, piece="PIECE"):
    """Write a FLUX field file (C++ core for real file paths, Python fallback)."""
    if piece == "PIECE" and not is_buffer(filename, "w"):
        try:
            _core.dex_write(str(filename), mesh)
            return
        except Exception:
            pass
    return _py_write(filename, mesh, piece=piece)


register_format("dex", [".dex"], read, {"dex": write})

__all__ = ["read", "write"]
