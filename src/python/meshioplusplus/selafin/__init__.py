"""TELEMAC Selafin results (roadmap 2.4.3.4): C++ core for real file paths, the
Python reference otherwise."""

from .. import _core  # noqa: F401
from .._helpers import register_format
from .._shims import series_shims
from ._selafin import read as _py_read  # noqa: F401
from ._selafin import time_values as _py_time_values  # noqa: F401
from ._selafin import write as _py_write  # noqa: F401

read, write, time_values = series_shims(
    __name__,
    "selafin",
    """Read a Selafin file's mesh and one time step (0 the first, negatives from the
    end): each variable is point data, its unit in ``field_data["selafin:unit:*"]``;
    ``IPOBO`` is ``selafin:ipobo``; a 3-D file's ``ELEVATION Z`` is the points' z.""",
    """Write a 2-D triangle or quad mesh as a one-step double-precision Selafin file.""",
    """Every complete step's time.""",
    selective=False,
)

register_format("selafin", [".slf", ".srf"], read, {"selafin": write})

__all__ = ["read", "write", "time_values"]
