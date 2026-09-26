from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._ensight import read as _py_read
from ._ensight import write as _py_write


def read(filename, time_step: int = 0):
    """Read an EnSight Gold .case/.geo pair (C++ core for real file paths, Python fallback).

    ``time_step`` selects one step of a transient .case file (0 = first,
    negative counts from the end), resolved the same way the C API/Fortran/
    Julia/R/WASM surfaces do -- see :func:`meshioplusplus.ensight.read`'s C++
    counterpart, ``read_ensight``. The Python reference reader is geometry-
    only (no VARIABLE support), so a non-default value forces the C++ path
    and re-raises rather than silently falling back to a mesh with no data.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.ensight_read(str(filename), time_step)
        except Exception as exc:
            if time_step:
                raise
            if not core_declined(exc, "ensight", "read", filename):
                raise
    return _py_read(filename)


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


register_format("ensight", [".case", ".geo"], read, {"ensight": write})

__all__ = ["read", "write"]
