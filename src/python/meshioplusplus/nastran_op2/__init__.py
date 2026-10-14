from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._op2 import read as _py_read
from ._op2 import time_values as _py_time_values
from ._write import SeriesWriter as _PySeriesWriter
from ._write import write as _py_write


def read(filename, points_only=False, arrays=None, time_step=0):
    """Read one step of a Nastran OP2 result file (MSC or NX, 32- or 64-bit).

    GEOM1 GRIDs are the points, GEOM2 elements the cells (``nastran:eid``,
    ``nastran:pid``, one region per property id); without GEOM tables the
    sibling input deck (``<stem>.bdf``/``.dat``/``.nas``/``.blk``) gives the mesh.
    Each (subcase, mode or time) of the displacement, eigenvector, velocity,
    acceleration, SPC/MPC force, applied load and temperature tables and of the
    real SORT1 stress and strain of rods, bars, shear panels, shells and solids
    is a step; ``time_step`` picks one (negative from the end) and its
    frequency, time or eigenvalue is ``field_data["meshio:time"]``. Other
    tables are skipped with a warning naming them.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.nastran_op2_read(
                str(filename),
                points_only=points_only,
                arrays=arrays,
                time_step=time_step,
            )
        except Exception as exc:
            if not core_declined(exc, "nastran_op2", "read", filename):
                raise
    return _py_read(
        filename, points_only=points_only, arrays=arrays, time_step=time_step
    )


def time_values(filename):
    """The time, frequency or eigenvalue of every step."""
    try:
        return list(_core.nastran_op2_time_values(str(filename)))
    except Exception as exc:
        if not core_declined(exc, "nastran_op2", "read", filename):
            raise
    return _py_time_values(filename)


def write(filename, mesh):
    """Write a mesh as a one-step Nastran OP2 result file.

    The inverse of :func:`read` for the core of the format: 32-bit
    little-endian words, the MSC dialect, SORT1 real tables. The model is the
    points (``GEOM1`` GRIDs, numbered from 1), the cells (``GEOM2`` element
    cards, with ``nastran:eid`` and ``nastran:pid`` when the mesh has them) and
    one property card per property id. The results are the point data
    ``DISPLACEMENT``, ``EIGENVECTOR``, ``VELOCITY``, ``ACCELERATION``,
    ``SPC_FORCE``, ``MPC_FORCE``, ``APPLIED_LOAD`` (each ``(n, 3)`` with a
    ``_ROT`` twin) and ``TEMPERATURE``, and the centre ``STRESS:<M>`` and
    ``STRAIN:<M>`` cell data of shells, solids and rods. The step is
    ``nastran:subcase``, ``nastran:analysis``, ``nastran:mode`` and
    ``meshio:time`` of the field data. Other arrays are dropped with a warning.
    """
    if not is_buffer(filename, "w"):
        try:
            _core.nastran_op2_write(str(filename), mesh)
            return
        except Exception as exc:
            if not core_declined(exc, "nastran_op2", "write", filename):
                raise
    return _py_write(filename, mesh)


class SeriesWriter:
    """Write a time series into one OP2 file: the model once, and every
    table's block per step, which ``read`` gives back by ``time_step``. The
    file is assembled when the writer closes. A step's time is the argument,
    not its ``meshio:time``; its subcase, analysis and mode come from its
    ``nastran:*`` field data. Every step must have the first step's cells and
    arrays, and its own subcase, analysis and mode or time.
    ``write_sequence(path, steps, file_format="nastran_op2")`` uses it.

    >>> with meshioplusplus.nastran_op2.SeriesWriter("run.op2") as w:
    ...     for time, mesh in steps:
    ...         w.write(time, mesh)

    The C++ core writes the file; the Python twin, which writes the same bytes,
    when the core declines to open it.
    """

    def __init__(self, filename):
        self._core = self._py = None
        if not is_buffer(filename, "w"):
            try:
                self._core = _core.NastranOp2SeriesWriter(str(filename))
            except Exception as exc:
                if not core_declined(exc, "nastran_op2", "write", filename):
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


register_format("nastran_op2", [".op2"], read, {"nastran_op2": write})

__all__ = ["read", "write", "time_values", "SeriesWriter"]
