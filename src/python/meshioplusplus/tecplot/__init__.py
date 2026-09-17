from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._tecplot import read as _py_read
from ._tecplot import write as _py_write


def read(filename, time_step: int = 0):
    """Read a Tecplot ASCII file (C++ core for real file paths, Python fallback).

    ``time_step`` selects one zone of a transient (``SOLUTIONTIME``/
    ``STRANDID``) file's timeline (0 = first, negative counts from the end),
    resolved the same way the C API/Fortran/Julia/R/WASM surfaces do -- see
    :func:`meshioplusplus.tecplot.read`'s C++ counterpart, ``read_tecplot``.
    A non-default value forces the C++ path (the Python reference has no
    transient-zone support) and re-raises rather than silently falling back.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.tecplot_read(str(filename), time_step)
        except Exception:
            if time_step:
                raise
    return _py_read(filename)


def write(filename, mesh):
    """Write a Tecplot ASCII file (C++ core for real file paths, Python fallback)."""
    if not is_buffer(filename, "w"):
        try:
            _core.tecplot_write(str(filename), mesh)
            return
        except Exception:
            pass
    return _py_write(filename, mesh)


register_format("tecplot", [".dat", ".tec"], read, {"tecplot": write})

__all__ = ["read", "write"]
