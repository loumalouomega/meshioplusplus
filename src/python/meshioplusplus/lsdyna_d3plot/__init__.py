from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._d3plot import read as _py_read
from ._d3plot import time_values as _py_time_values


def read(filename, points_only=False, arrays=None, time_step=0):
    """Read one state of an LS-DYNA d3plot family.

    Open the base file (``d3plot``): its numbered members (``d3plot01``...)
    are found beside it. ``time_step`` picks the state (negative counts from
    the end); its time is ``field_data["meshio:time"]``. Nodal variables are
    point data (``displacement`` is relative to the initial coordinates, which
    stay the points), element variables cell data (NaN on blocks without them;
    per integration point or layer as ``(cells, points * components)`` with
    ``field_data["lsdyna_d3plot:layout:<name>"]``), the deletion flags the int8
    mask ``lsdyna:alive`` and the globals field data. Parts are cell regions.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.lsdyna_d3plot_read(
                str(filename),
                points_only=points_only,
                arrays=arrays,
                time_step=time_step,
            )
        except Exception as exc:
            if not core_declined(exc, "lsdyna_d3plot", "read", filename):
                raise
    return _py_read(
        filename, points_only=points_only, arrays=arrays, time_step=time_step
    )


def time_values(filename):
    """The times of the family's states."""
    try:
        return list(_core.lsdyna_d3plot_time_values(str(filename)))
    except Exception as exc:
        if not core_declined(exc, "lsdyna_d3plot", "read", filename):
            raise
    return _py_time_values(filename)


# Found by its file name (``d3plot``) or its control block, never an extension.
register_format("lsdyna_d3plot", [], read, {})

__all__ = ["read", "time_values"]
