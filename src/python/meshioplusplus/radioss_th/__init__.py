from .. import _core  # noqa: F401
from .._helpers import register_format
from .._shims import series_shims
from ._th import SeriesWriter as _PySeriesWriter  # noqa: F401
from ._th import read as _py_read  # noqa: F401
from ._th import time_values as _py_time_values  # noqa: F401
from ._th import write as _py_write  # noqa: F401

read, write, time_values, SeriesWriter = series_shims(
    __name__,
    "radioss_th",
    """Read one output of an OpenRadioss / Radioss time-history file
    (``<run>T01``...).

    The steps are the file's outputs. Each is a mesh with no points whose
    field data holds ``meshio:time``, the global variables
    (``radioss_th:global:kinetic_energy``...), the part and subset variables
    (``radioss_th:part:<id>:IE``...) and, per TH group,
    ``radioss_th:<kind>:<id>`` (entities x variables) with its ``:ids`` and
    ``:variables`` (the variable codes).
    """,
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
    """,
    """The time of every complete output.""",
    core_writer="RadiossThSeriesWriter",
    writer_doc="""Write a time series into one time-history file: the descriptions with
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
    """,
)

# Found by its file name (a stem, then T and two digits) or its title record,
# never an extension.
register_format("radioss_th", [], read, {"radioss_th": write})

__all__ = ["read", "write", "time_values", "SeriesWriter"]
