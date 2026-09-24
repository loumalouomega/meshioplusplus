from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
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


register_format("tecplot", [".dat", ".tec", ".plt"], read, {"tecplot": write})

__all__ = ["read", "write"]
