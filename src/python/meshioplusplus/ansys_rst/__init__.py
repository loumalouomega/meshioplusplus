from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._ansys_rst import read as _py_read
from ._ansys_rst import time_values as _py_time_values


def read(filename, points_only=False, arrays=None, time_step=0, lenient=False):
    """Read one result set of an Ansys MAPDL results file (``.rst``, ``.rth``).

    The geometry records become cells exactly as a ``.cdb`` deck's elements do
    (degenerate shapes resolved, ``ansys:*`` cell data, components as point and
    cell regions). ``time_step`` picks the result set (negative counts from the
    end); its time, or frequency in a modal analysis, is
    ``field_data["meshio:time"]``, with ``ansys:load_step``, ``ansys:substep``
    and ``ansys:cumulative``. The set's nodal DOF solution is point data: ``U``,
    ``ROT``, ``A`` and ``V`` as vectors rotated to the global axes, other DOFs
    (``TEMP``, ``PRES`` ...) as scalars; NaN where a node has no value.
    ``lenient`` skips elements whose type has no meshio++ cell.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.ansys_rst_read(
                str(filename),
                points_only=points_only,
                arrays=arrays,
                time_step=time_step,
                lenient=lenient,
            )
        except Exception as exc:
            if not core_declined(exc, "ansys_rst", "read", filename):
                raise
    return _py_read(
        filename,
        points_only=points_only,
        arrays=arrays,
        time_step=time_step,
        lenient=lenient,
    )


def time_values(filename):
    """The time (or frequency) of every result set."""
    try:
        return list(_core.ansys_rst_time_values(str(filename)))
    except Exception as exc:
        if not core_declined(exc, "ansys_rst", "read", filename):
            raise
    return _py_time_values(filename)


register_format("ansys_rst", [".rst", ".rth"], read, {})

__all__ = ["read", "time_values"]
