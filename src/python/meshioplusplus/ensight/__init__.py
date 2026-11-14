from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from .._series_writer import CoreSeriesWriter
from ._ensight import SeriesWriter as _PySeriesWriter
from ._ensight import read as _py_read
from ._ensight import write as _py_write


def read(filename, time_step: int = 0):
    """Read an EnSight Gold .case/.geo pair (C++ core for real file paths, Python fallback).

    ``time_step`` selects one step of a transient .case file (0 = first,
    negative counts from the end), resolved the same way the C API/Fortran/
    Julia/R/WASM surfaces do -- see :func:`meshioplusplus.ensight.read`'s C++
    counterpart, ``read_ensight``; the Python twin resolves it identically.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.ensight_read(str(filename), time_step)
        except Exception as exc:
            if not core_declined(exc, "ensight", "read", filename):
                raise
    return _py_read(filename, time_step)


def write(filename, mesh, binary=True, fortran=False):
    """Write an EnSight Gold .case/.geo pair (C++ core for real file paths, Python fallback).

    ``fortran=True`` writes Fortran binary: the C-binary records, each framed as
    a Fortran sequential unformatted record, under a ``Fortran Binary`` header.
    """
    if not is_buffer(filename, "w"):
        try:
            _core.ensight_write(str(filename), mesh, binary, fortran)
            return
        except Exception as exc:
            if not core_declined(exc, "ensight", "write", filename):
                raise
    return _py_write(filename, mesh, binary=binary, fortran=fortran)


class SeriesWriter(CoreSeriesWriter):
    """Write a time series as one EnSight Gold case: the geometry once, one numbered
    variable file per step (``name.0000.scl``, ...) and a ``TIME`` section in the
    ``.case`` file, which is written when the writer closes. Every step must have the
    first step's points, cell blocks and variables.

    >>> with meshioplusplus.ensight.SeriesWriter("run.case") as w:
    ...     for t, mesh in steps:
    ...         w.write(t, mesh)
    """

    def __init__(self, filename, binary=True, fortran=False):
        self._core = self._py = None
        if not is_buffer(filename, "w"):
            try:
                self._core = _core.EnsightSeriesWriter(str(filename), binary, fortran)
                return
            except Exception as exc:
                if not core_declined(exc, "ensight", "write", filename):
                    raise
        self._py = _PySeriesWriter(filename, binary, fortran)


register_format("ensight", [".case", ".geo"], read, {"ensight": write})

__all__ = ["read", "write", "SeriesWriter"]
