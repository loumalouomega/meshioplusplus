from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._binout import read as _py_read
from ._binout import time_values as _py_time_values
from ._write import SeriesWriter as _PySeriesWriter
from ._write import write as _py_write


def read(filename, points_only=False, arrays=None, time_step=0):
    """Read one output of an LS-DYNA binout (``binout``, ``binout%001``...).

    The steps are the ``nodout`` outputs: each is a point cloud of the
    ``nodout`` nodes at their coordinates, with displacement, rotation,
    velocity and acceleration as point data and ``lsdyna:nid`` their ids.
    Every other database's variables at the same time (``glstat``,
    ``matsum``, ``elout/beam``...) are field data
    ``binout:<database>:<variable>``. Without ``nodout`` the steps are the
    first database's outputs and the mesh has no points.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.lsdyna_binout_read(
                str(filename),
                points_only=points_only,
                arrays=arrays,
                time_step=time_step,
            )
        except Exception as exc:
            if not core_declined(exc, "lsdyna_binout", "read", filename):
                raise
    return _py_read(filename, points_only, arrays, time_step)


def time_values(filename):
    """The time of every output (``nodout``'s, else the first database's)."""
    try:
        return list(_core.lsdyna_binout_time_values(str(filename)))
    except Exception as exc:
        if not core_declined(exc, "lsdyna_binout", "read", filename):
            raise
    return _py_time_values(filename)


def write(filename, mesh, split_size=0, word_size=8, big_endian=False):
    """Write a mesh as a one-output LS-DYNA binout (``binout``).

    The inverse of :func:`read`: a single LSDA file. ``nodout`` is the points
    (``x_coordinate``...), their ids (``lsdyna:nid``, else 1 to n), the vectors
    ``displacement``, ``rotation``, ``velocity``, ``rotational_velocity``,
    ``acceleration`` and ``rotational_acceleration`` (``(n, 3)``), any other
    per-point array under its own name and ``binout:nodout:<variable>`` field
    data. Every other database is ``binout:<database>:<variable>`` field data
    (``binout:glstat:kinetic_energy``, ``binout:elout/beam:axial``) with its
    ``binout:<database>:time`` and ``binout:<database>:ids``. Reals are stored
    as ``float32`` and integers as ``int32``. The output's time is
    ``meshio:time`` (0 when absent). Cells, regions and other arrays are
    dropped with a warning. ``split_size`` (bytes, 0 for none) starts a
    continuation file (``binout%001``, ``%002``...) once the data written to the
    current one would pass it; an output stays whole in one file, and ``read``
    merges them. ``word_size`` (8, or 4 for a file under 4 GiB) sets the bytes of every
    length and offset in the LSDA records and ``big_endian`` their byte order and the
    data's.
    """
    if not is_buffer(filename, "w"):
        try:
            _core.lsdyna_binout_write(
                str(filename), mesh, split_size, word_size, big_endian
            )
            return
        except Exception as exc:
            if not core_declined(exc, "lsdyna_binout", "write", filename):
                raise
    return _py_write(filename, mesh, split_size, word_size, big_endian)


class SeriesWriter:
    """Write a time series into one binout: ``nodout``'s ids once and an output
    per step, and any other database's output whenever its own
    ``binout:<database>:time`` changes, so a database written less often than
    ``nodout`` keeps its own cadence (the inverse of what ``read`` attaches to a
    step). A step's time is the argument, not its ``meshio:time``. Every step
    must have the first step's points and each database's variables.
    ``write_sequence(path, steps, file_format="lsdyna_binout")`` uses it.

    >>> with meshioplusplus.lsdyna_binout.SeriesWriter("binout") as w:
    ...     for time, mesh in steps:
    ...         w.write(time, mesh)

    ``split_size``, ``word_size`` and ``big_endian`` are as in
    :func:`write`. The C++ core writes the file; the Python twin, which writes
    the same bytes, when the core declines to open it.
    """

    def __init__(self, filename, split_size=0, word_size=8, big_endian=False):
        self._core = self._py = None
        if not is_buffer(filename, "w"):
            try:
                self._core = _core.LsdynaBinoutSeriesWriter(
                    str(filename), split_size, word_size, big_endian
                )
            except Exception as exc:
                if not core_declined(exc, "lsdyna_binout", "write", filename):
                    raise
        if self._core is None:
            self._py = _PySeriesWriter(filename, split_size, word_size, big_endian)

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


# Found by its file name (``binout``...) or its LSDA header, never an extension.
register_format("lsdyna_binout", [], read, {"lsdyna_binout": write})

__all__ = ["read", "write", "time_values", "SeriesWriter"]
