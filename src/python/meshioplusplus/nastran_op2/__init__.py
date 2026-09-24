from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._op2 import read as _py_read
from ._op2 import time_values as _py_time_values


def read(filename, points_only=False, arrays=None, time_step=0):
    """Read one step of a Nastran OP2 result file (MSC or NX, 32- or 64-bit).

    GEOM1 GRIDs are the points, GEOM2 elements the cells (``nastran:eid``,
    ``nastran:pid``, one region per property id); without GEOM tables the
    sibling input deck (``<stem>.bdf``/``.dat``/``.nas``/``.blk``) gives the mesh.
    Each (subcase, mode or time) of the displacement, eigenvector, velocity,
    acceleration, SPC/MPC force, applied load and temperature tables and of the
    real SORT1 stress and strain of rods, bars, shear panels, shells and solids
    is a step; ``time_step`` picks one (negative from the end) and its
    frequency, time or eigenvalue is ``field_data["meshio:time"]``. Other
    tables are skipped with a warning naming them.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.nastran_op2_read(
                str(filename),
                points_only=points_only,
                arrays=arrays,
                time_step=time_step,
            )
        except Exception as exc:
            if not core_declined(exc, "nastran_op2", "read", filename):
                raise
    return _py_read(
        filename, points_only=points_only, arrays=arrays, time_step=time_step
    )


def time_values(filename):
    """The time, frequency or eigenvalue of every step."""
    try:
        return list(_core.nastran_op2_time_values(str(filename)))
    except Exception as exc:
        if not core_declined(exc, "nastran_op2", "read", filename):
            raise
    return _py_time_values(filename)


register_format("nastran_op2", [".op2"], read, {})

__all__ = ["read", "time_values"]
