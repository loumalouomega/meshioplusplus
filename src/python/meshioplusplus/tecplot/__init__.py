from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from .._series_writer import CoreSeriesWriter
from ._tecplot import SeriesWriter as _PySeriesWriter
from ._tecplot import read as _py_read
from ._tecplot import write as _py_write


def read(filename, time_step: int = 0):
    """Read a Tecplot ASCII (``.dat``/``.tec``) or binary (``.plt``) file (C++
    core for real file paths, Python fallback).

    Every zone of the selected step is a cell block and a Cell region;
    ``time_step`` selects one step of a transient (``SOLUTIONTIME``/
    ``STRANDID``) file's timeline (0 = first, negative counts from the end),
    resolved the same way the C API/Fortran/Julia/R/WASM surfaces do.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.tecplot_read(str(filename), time_step)
        except Exception as exc:
            if not core_declined(exc, "tecplot", "read", filename):
                raise
    return _py_read(filename, time_step)


def write(filename, mesh):
    """Write a Tecplot ASCII file (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "w"):
        try:
            _core.tecplot_write(str(filename), mesh)
            return
        except Exception as exc:
            if not core_declined(exc, "tecplot", "write", filename):
                raise
    return _py_write(filename, mesh)


class SeriesWriter(CoreSeriesWriter):
    """Write a time series into one Tecplot ASCII file: every step's zones, each with
    ``STRANDID = 1`` and the step's ``SOLUTIONTIME`` (what ``read(..., time_step=k)``
    selects by). A step is written as a file of its own and its zones appended, so each
    step repeats the coordinates; every step must have the same variables. A ``.plt``
    path is refused.

    >>> with meshioplusplus.tecplot.SeriesWriter("run.dat") as w:
    ...     for t, mesh in steps:
    ...         w.write(t, mesh)
    """

    def __init__(self, filename):
        self._core = self._py = None
        if not is_buffer(filename, "w"):
            try:
                self._core = _core.TecplotSeriesWriter(str(filename))
                return
            except Exception as exc:
                if not core_declined(exc, "tecplot", "write", filename):
                    raise
        self._py = _PySeriesWriter(filename)


register_format("tecplot", [".dat", ".tec", ".plt"], read, {"tecplot": write})

__all__ = ["read", "write", "SeriesWriter"]
