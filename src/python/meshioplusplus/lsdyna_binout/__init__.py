from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._binout import read as _py_read
from ._binout import time_values as _py_time_values


def read(filename, points_only=False, arrays=None, time_step=0):
    """Read one output of an LS-DYNA binout (``binout``, ``binout%001``...).

    The steps are the ``nodout`` outputs: each is a point cloud of the
    ``nodout`` nodes at their coordinates, with displacement, rotation,
    velocity and acceleration as point data and ``lsdyna:nid`` their ids.
    Every other database's variables at the same time (``glstat``,
    ``matsum``, ``elout/beam``...) are field data
    ``binout:<database>:<variable>``. Without ``nodout`` the steps are the
    first database's outputs and the mesh has no points.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.lsdyna_binout_read(
                str(filename),
                points_only=points_only,
                arrays=arrays,
                time_step=time_step,
            )
        except Exception as exc:
            if not core_declined(exc, "lsdyna_binout", "read", filename):
                raise
    return _py_read(filename, points_only, arrays, time_step)


def time_values(filename):
    """The time of every output (``nodout``'s, else the first database's)."""
    try:
        return list(_core.lsdyna_binout_time_values(str(filename)))
    except Exception as exc:
        if not core_declined(exc, "lsdyna_binout", "read", filename):
            raise
    return _py_time_values(filename)


# Found by its file name (``binout``...) or its LSDA header, never an extension.
register_format("lsdyna_binout", [], read, {})

__all__ = ["read", "time_values"]
