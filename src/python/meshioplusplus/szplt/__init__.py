from .. import _core
from .._helpers import register_format
from ._szplt import read as _py_read
from ._szplt import time_values as _py_time_values
from ._write import SeriesWriter as _PySeriesWriter
from ._write import write as _py_write

_HAS_TECIO = getattr(_core, "__has_tecio__", False)


def read(filename, time_step=0):
    """Read a Tecplot SZL (``.szplt``) file through TecIO.

    It reads exactly as the same data saved as ``.plt`` would: one cell block
    and Cell region per zone of the selected step, ``time_step`` picking a
    solution time. Needs a core built with ``MESHIOPLUSPLUS_WITH_TECIO=ON``, or a
    shared TecIO named by ``MESHIOPLUSPLUS_TECIO_LIBRARY`` (Tecplot 360 ships
    one); without either, save the file as ``.plt`` in Tecplot.
    """
    # No fallback from a core built with TecIO: the native reader is the whole
    # of this format, and the ctypes twin needs a second, shared TecIO.
    if _HAS_TECIO:
        return _core.szplt_read(str(filename), time_step=time_step)
    return _py_read(filename, time_step)


def time_values(filename):
    """The distinct solution times, ascending (empty for a static file)."""
    if _HAS_TECIO:
        return list(_core.szplt_time_values(str(filename)))
    return _py_time_values(filename)


def write(filename, mesh):
    """Write a Tecplot SZL (``.szplt``) file through TecIO.

    The zones are the ASCII and ``.plt`` writers': one per cell block, the first
    owning the coordinates and the nodal variables and the later ones sharing
    them, doubles throughout, wedges and pyramids as degenerate bricks. Needs a
    core built with ``MESHIOPLUSPLUS_WITH_TECIO=ON``, or a shared TecIO named by
    ``MESHIOPLUSPLUS_TECIO_LIBRARY``.
    """
    if _HAS_TECIO:
        return _core.szplt_write(str(filename), mesh)
    return _py_write(filename, mesh)


class SeriesWriter:
    """Write a time series into one ``.szplt`` file, a strand of zones per step;
    ``write_sequence(path, steps, file_format="szplt")`` uses it.

    >>> with meshioplusplus.szplt.SeriesWriter("run.szplt") as w:
    ...     for time, mesh in steps:
    ...         w.write(time, mesh)
    """

    def __init__(self, filename):
        if _HAS_TECIO:
            self._core, self._py = _core.SzpltSeriesWriter(str(filename)), None
        else:
            self._core, self._py = None, _PySeriesWriter(filename)

    def write(self, time, mesh):
        if self._core is not None:
            self._core.write(float(time), mesh)
        else:
            self._py.write(time, mesh)

    def close(self):
        if self._core is not None:
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


register_format("szplt", [".szplt"], read, {"szplt": write})

__all__ = ["read", "write", "time_values", "SeriesWriter"]
