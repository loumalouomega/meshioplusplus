from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._th import read as _py_read
from ._th import time_values as _py_time_values


def read(filename, points_only=False, arrays=None, time_step=0):
    """Read one output of an OpenRadioss / Radioss time-history file
    (``<run>T01``...).

    The steps are the file's outputs. Each is a mesh with no points whose
    field data holds ``meshio:time``, the global variables
    (``radioss_th:global:kinetic_energy``...), the part and subset variables
    (``radioss_th:part:<id>:IE``...) and, per TH group,
    ``radioss_th:<kind>:<id>`` (entities x variables) with its ``:ids`` and
    ``:variables`` (the variable codes).
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.radioss_th_read(
                str(filename),
                points_only=points_only,
                arrays=arrays,
                time_step=time_step,
            )
        except Exception as exc:
            if not core_declined(exc, "radioss_th", "read", filename):
                raise
    return _py_read(filename, points_only, arrays, time_step)


def time_values(filename):
    """The time of every complete output."""
    try:
        return list(_core.radioss_th_time_values(str(filename)))
    except Exception as exc:
        if not core_declined(exc, "radioss_th", "read", filename):
            raise
    return _py_time_values(filename)


# Found by its file name (a stem, then T and two digits) or its title record,
# never an extension.
register_format("radioss_th", [], read, {})

__all__ = ["read", "time_values"]
