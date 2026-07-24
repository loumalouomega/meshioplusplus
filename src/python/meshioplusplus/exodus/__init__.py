from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._exodus import read as _py_read
from ._exodus import write as _py_write

_HAS_NETCDF = getattr(_core, "__has_netcdf__", False)


def read(filename, time_step=0):
    """Read an Exodus II file (C++ core when built with netCDF, Python fallback).

    :param time_step: which step of a multi-step file to materialize; 0 (the
        default) is the first, negative counts from the end. Out of range is an
        error naming the available count.
    """
    if _HAS_NETCDF and not is_buffer(filename, "r"):
        try:
            return _core.exodus_read(str(filename), time_step=time_step)
        except Exception:
            # The standard shim contract: any refusal from the core retries on
            # the Python twin. A genuine user error (an out-of-range time step)
            # is not lost -- the twin implements the same rule and raises the
            # same way, so declining costs a slower answer, not a wrong one.
            pass
    return _py_read(filename, time_step=time_step)


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
