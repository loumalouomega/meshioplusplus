from .. import _core
from .._fallback import core_declined
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
    if _HAS_TECIO:
        try:
            return _core.szplt_read(str(filename), time_step=time_step)
        except Exception as exc:
            if not core_declined(exc, "szplt", "read", filename):
                raise
    return _py_read(filename, time_step)


def time_values(filename):
    """The distinct solution times, ascending (empty for a static file)."""
    if _HAS_TECIO:
        try:
            return list(_core.szplt_time_values(str(filename)))
        except Exception as exc:
            if not core_declined(exc, "szplt", "read", filename):
                raise
    return _py_time_values(filename)


register_format("szplt", [".szplt"], read, {})

__all__ = ["read", "time_values"]
