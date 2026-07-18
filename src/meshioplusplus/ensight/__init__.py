from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._ensight import read as _py_read
from ._ensight import write as _py_write


def read(filename):
    """Read an EnSight Gold .case/.geo pair (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.ensight_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh, binary=True):
    """Write an EnSight Gold .case/.geo pair (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "w"):
        try:
            _core.ensight_write(str(filename), mesh, binary)
            return
        except Exception:
            pass
    return _py_write(filename, mesh, binary=binary)


register_format("ensight", [".case", ".geo"], read, {"ensight": write})

__all__ = ["read", "write"]
