from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._th import SeriesWriter as _PySeriesWriter
from ._th import read as _py_read
from ._th import time_values as _py_time_values
from ._th import write as _py_write


def read(filename, points_only=False, arrays=None, time_step=0):
    """Read one output of an OpenRadioss / Radioss time-history file
    (``<run>T01``...).

    The steps are the file's outputs. Each is a mesh with no points whose
    field data holds ``meshio:time``, the global variables
    (``radioss_th:global:kinetic_energy``...), the part and subset variables
    (``radioss_th:part:<id>:IE``...) and, per TH group,
    ``radioss_th:<kind>:<id>`` (entities x variables) with its ``:ids`` and
    ``:variables`` (the variable codes).
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.radioss_th_read(
                str(filename),
                points_only=points_only,
                arrays=arrays,
                time_step=time_step,
            )
        except Exception as exc:
            if not core_declined(exc, "radioss_th", "read", filename):
                raise
    return _py_read(filename, points_only, arrays, time_step)


def time_values(filename):
    """The time of every complete output."""
    try:
        return list(_core.radioss_th_time_values(str(filename)))
    except Exception as exc:
        if not core_declined(exc, "radioss_th", "read", filename):
            raise
    return _py_time_values(filename)


def write(filename, mesh):
    """Write a mesh's ``radioss_th:*`` field data as a one-output
    OpenRadioss / Radioss time-history file (``<run>T01``).

    The inverse of :func:`read`: ``meshio:time`` is the output's time,
    ``radioss_th:global:<name>``, ``radioss_th:part:<id>:<name>`` and
    ``radioss_th:subset:<id>:<name>`` are the scalar variables (``var<code>``
    for a code with no name), and per TH group ``radioss_th:<kind>:<id>``
    (entities x variables) with its optional ``:ids`` and ``:variables``
    (else 1..n) is a record. Everything is stored as big-endian ``float32``.
    With ``radioss_th:unit_factors`` the file is version 4021, else 3040
    (the version OpenRadioss writes).
    Points, cells and other field data are not stored (the latter is dropped
    with a warning).
    """
    if not is_buffer(filename, "w"):
        try:
            _core.radioss_th_write(str(filename), mesh)
            return
        except Exception as exc:
            if not core_declined(exc, "radioss_th", "write", filename):
                raise
    return _py_write(filename, mesh)


class SeriesWriter:
    """Write a time series into one time-history file: the descriptions with
    the first step, then one output per step, which ``read`` gives back by
    ``time_step``. Every step must have the first step's variables (the same
    names, ids and codes); a step's time is the argument, not its
    ``meshio:time``. ``write_sequence(path, steps, file_format="radioss_th")``
    uses it.

    >>> with meshioplusplus.radioss_th.SeriesWriter("runT01") as w:
    ...     for time, mesh in steps:
    ...         w.write(time, mesh)

    The C++ core writes the file; the Python twin, which writes the same bytes,
    when the core declines to open it.
    """

    def __init__(self, filename):
        self._core = self._py = None
        if not is_buffer(filename, "w"):
            try:
                self._core = _core.RadiossThSeriesWriter(str(filename))
            except Exception as exc:
                if not core_declined(exc, "radioss_th", "write", filename):
                    raise
        if self._core is None:
            self._py = _PySeriesWriter(filename)

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


# Found by its file name (a stem, then T and two digits) or its title record,
# never an extension.
register_format("radioss_th", [], read, {"radioss_th": write})

__all__ = ["read", "write", "time_values", "SeriesWriter"]
