from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._avsucd import read as _py_read
from ._avsucd import write as _py_write


def read(filename):
    """Read an AVS-UCD file (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.avsucd_read(str(filename))
        except Exception as exc:
            if not core_declined(exc, "avsucd", "read", filename):
                raise
    return _py_read(filename)


def write(filename, mesh):
    """Write an AVS-UCD file (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "w"):
        try:
            _core.avsucd_write(str(filename), mesh)
            return
        except Exception as exc:
            if not core_declined(exc, "avsucd", "write", filename):
                raise
    return _py_write(filename, mesh)


register_format("avsucd", [".avs"], read, {"avsucd": write})

__all__ = ["read", "write"]
