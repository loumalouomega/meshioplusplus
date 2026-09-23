from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._xplt import read as _py_read
from ._xplt import time_values as _py_time_values


def read(filename, points_only=False, arrays=None, time_step=0, lenient=False):
    """Read one state of an FEBio plot file (``.xplt``).

    Domains become cell blocks and cell regions; node sets, element sets and
    surfaces become point, cell and side regions. ``time_step`` picks the state
    (negative counts from the end); its time is ``field_data["meshio:time"]``.
    Nodal variables are point data, per-element ones cell data (NaN on domains
    without them), per-element-node ones point data averaged over the elements
    sharing a node, global ones field data. Surface and edge variables are not
    read. A compressed file needs zlib, which the Python fallback always has.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.xplt_read(
                str(filename),
                points_only=points_only,
                arrays=arrays,
                time_step=time_step,
                lenient=lenient,
            )
        except Exception as exc:
            if not core_declined(exc, "xplt", "read", filename):
                raise
    return _py_read(
        filename,
        points_only=points_only,
        arrays=arrays,
        time_step=time_step,
        lenient=lenient,
    )


def time_values(filename):
    """The times of the plot file's states."""
    try:
        return list(_core.xplt_time_values(str(filename)))
    except Exception as exc:
        if not core_declined(exc, "xplt", "read", filename):
            raise
    return _py_time_values(filename)


register_format("xplt", [".xplt"], read, {})

__all__ = ["read", "time_values"]
