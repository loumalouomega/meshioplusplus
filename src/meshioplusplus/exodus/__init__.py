from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._exodus import read as _py_read
from ._exodus import write as _py_write

_HAS_NETCDF = getattr(_core, "__has_netcdf__", False)


def read(filename):
    """Read an Exodus II file (C++ core when built with netCDF, Python fallback)."""
    if _HAS_NETCDF and not is_buffer(filename, "r"):
        try:
            return _core.exodus_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh):
    """Write an Exodus II file (C++ core when built with netCDF, Python fallback)."""
    # Node sets (point_sets) live outside the conversion layer -> Python.
    if _HAS_NETCDF and not mesh.point_sets and not is_buffer(filename, "w"):
        try:
            _core.exodus_write(str(filename), mesh)
            return
        except Exception:
            pass
    return _py_write(filename, mesh)


register_format("exodus", [".e", ".exo", ".ex2"], read, {"exodus": write})

__all__ = ["read", "write"]
