from .. import _core  # noqa: F401
from .._helpers import register_format
from .._shims import series_shims
from ._op2 import read as _py_read  # noqa: F401
from ._op2 import time_values as _py_time_values  # noqa: F401
from ._write import SeriesWriter as _PySeriesWriter  # noqa: F401
from ._write import write as _py_write  # noqa: F401

read, write, time_values, SeriesWriter = series_shims(
    __name__,
    "nastran_op2",
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
    """,
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
    """,
    """The time, frequency or eigenvalue of every step.""",
    core_writer="NastranOp2SeriesWriter",
    writer_doc="""Write a time series into one OP2 file: the model once, and every
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
    """,
)

register_format("nastran_op2", [".op2"], read, {"nastran_op2": write})

__all__ = ["read", "write", "time_values", "SeriesWriter"]
