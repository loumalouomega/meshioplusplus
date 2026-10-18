from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._write import SeriesWriter as _PySeriesWriter
from ._write import write as _py_write
from ._xplt import read as _py_read
from ._xplt import time_values as _py_time_values


def read(filename, points_only=False, arrays=None, time_step=0, lenient=False):
    """Read one state of an FEBio plot file (``.xplt``).

    Domains become cell blocks and cell regions; node sets, element sets and
    surfaces become point, cell and side regions. ``time_step`` picks the state
    (negative counts from the end); its time is ``field_data["meshio:time"]``.
    Nodal variables are point data, per-element ones cell data (NaN on domains
    without them), per-element-node ones point data averaged over the elements
    sharing a node, global ones field data. Surface and edge variables are cell
    data on the facet and segment cells (NaN elsewhere). A compressed file needs zlib, which the Python fallback always has.
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


def write(filename, mesh, compress=False):
    """Write an FEBio plot file holding one state: the mesh and its data, at
    ``field_data["meshio:time"]``. Cell blocks become domains (one element type
    each), point regions node sets and cell regions element sets; point data
    are nodal variables, cell data per-element variables and numeric field data
    global ones. Surfaces, edges and facet sets are written from their markers
    and side regions, and cell data with finite values on a surface's or edge's
    cells is a surface or edge variable. Vertices are not written. The file is
    little-endian, plot version 0x35; ``compress=True`` stores each state as a
    zlib stream (as FEBio's ``compression`` option does), which needs zlib in the
    C++ core or falls back to the Python writer."""
    if not is_buffer(filename, "w"):
        try:
            return _core.xplt_write(str(filename), mesh, bool(compress))
        except Exception as exc:
            if not core_declined(exc, "xplt", "write", filename):
                raise
    return _py_write(filename, mesh, compress)


class SeriesWriter:
    """Write a time series into one plot file: the mesh once and a state per
    step; ``write_sequence(path, steps, file_format="xplt")`` uses it.

    >>> with meshioplusplus.xplt.SeriesWriter("run.xplt") as w:
    ...     for time, mesh in steps:
    ...         w.write(time, mesh)
    """

    def __init__(self, filename, compress=False):
        self._core = self._py = None
        if not is_buffer(filename, "w"):
            try:
                self._core = _core.XpltSeriesWriter(str(filename), bool(compress))
            except Exception as exc:
                if not core_declined(exc, "xplt", "write", filename):
                    raise
        if self._core is None:
            self._py = _PySeriesWriter(filename, compress)

    def write(self, time, mesh):
        if self._core is not None:
            self._core.write(float(time), mesh)
        else:
            self._py.write(time, mesh)

    def close(self):
        if self._core is not None:
            if self._core.num_steps():
                self._core.finalize()
            self._core = None
        if self._py is not None:
            self._py.close()
            self._py = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False


register_format("xplt", [".xplt"], read, {"xplt": write})

__all__ = ["read", "write", "time_values", "SeriesWriter"]
