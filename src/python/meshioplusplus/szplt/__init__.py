from .. import _core
from .._helpers import register_format
from ._szplt import read as _py_read
from ._szplt import time_values as _py_time_values

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


register_format("szplt", [".szplt"], read, {})

__all__ = ["read", "time_values"]
