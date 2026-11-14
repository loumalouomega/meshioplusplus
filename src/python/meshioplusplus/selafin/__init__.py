"""TELEMAC Selafin results (roadmap 2.4.3.4): C++ core for real file paths, the
Python reference otherwise."""

from .. import _core
from .._fallback import core_declined
from .._helpers import register_format
from ._selafin import read as _py_read
from ._selafin import time_values as _py_time_values
from ._selafin import write as _py_write


def read(filename, time_step=0):
    """Read a Selafin file's mesh and one time step (0 the first, negatives from the
    end): each variable is point data, its unit in ``field_data["selafin:unit:*"]``;
    ``IPOBO`` is ``selafin:ipobo``; a 3-D file's ``ELEVATION Z`` is the points' z."""
    try:
        return _core.selafin_read(str(filename), int(time_step))
    except Exception as exc:
        if not core_declined(exc, "selafin", "read", filename):
            raise
    return _py_read(filename, time_step=time_step)


def time_values(filename):
    """Every complete step's time."""
    try:
        return list(_core.selafin_time_values(str(filename)))
    except Exception as exc:
        if not core_declined(exc, "selafin", "read", filename):
            raise
    return _py_time_values(filename)


def write(filename, mesh):
    """Write a 2-D triangle or quad mesh as a one-step double-precision Selafin file."""
    try:
        _core.selafin_write(str(filename), mesh)
        return
    except Exception as exc:
        if not core_declined(exc, "selafin", "write", filename):
            raise
    _py_write(filename, mesh)


register_format("selafin", [".slf", ".srf"], read, {"selafin": write})

__all__ = ["read", "write", "time_values"]
