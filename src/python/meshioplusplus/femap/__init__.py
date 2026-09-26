from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._femap import SeriesWriter as _PySeriesWriter
from ._femap import read as _py_read
from ._femap import time_values as _py_time_values
from ._femap import write as _py_write


def read(filename, points_only=False, arrays=None, time_step=0):
    """Read a Femap neutral file (``.neu``).

    Nodes (403) and elements (404) become the mesh, with the element property and
    type as the ``femap:property``/``femap:type`` cell data; properties (402) name
    the ``property_<id>`` cell regions and groups (408) become point and cell
    regions. Every output set (450) is a step: ``time_step`` selects one (0 =
    first, negative counts from the end), its value is
    ``field_data["meshio:time"]``, and its output vectors (451/1051) become point
    or cell data named by their titles, NaN where a vector has no value.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.femap_read(str(filename), points_only, arrays, time_step)
        except Exception as exc:
            if not core_declined(exc, "femap", "read", filename):
                raise
    return _py_read(
        filename, points_only=points_only, arrays=arrays, time_step=time_step
    )


def time_values(filename):
    """The value (time, frequency, load factor...) of every output set, in file
    order."""
    if not is_buffer(filename, "r"):
        try:
            return list(_core.femap_time_values(str(filename)))
        except Exception as exc:
            if not core_declined(exc, "femap", "read", filename):
                raise
    return _py_time_values(filename)


def write(filename, mesh):
    """Write a Femap 8.2 neutral file (``.neu``): header, properties, nodes,
    elements, groups and one output set of results.

    Element properties and types come from ``femap:property``/``femap:type`` (else
    1 and a type derived from the cell); other point and cell regions become
    groups. Numeric point and cell data become the output set's nodal and
    elemental vectors, one per component (``<name>_0``...); its id is
    ``femap:set`` and its value ``meshio:time``. Side regions, other data and
    cell types with no Femap topology are dropped with a warning.
    """
    if not is_buffer(filename, "w"):
        try:
            _core.femap_write(str(filename), mesh)
            return
        except Exception as exc:
            if not core_declined(exc, "femap", "write", filename):
                raise
    return _py_write(filename, mesh)


class SeriesWriter:
    """Write a time series into one neutral file (v16.17.0): the mesh once,
    then an output set (450) with its vectors (451) per step, which ``read``
    gives back by ``time_step``. A step's set id is its ``femap:set`` when
    positive and unused, else the next free one; its value is the step's time.
    Every step must have the first step's cells; a step whose points moved is
    written with the first step's, with a warning. ``write_sequence(path.neu,
    steps)`` uses it.

    >>> with meshioplusplus.femap.SeriesWriter("run.neu") as w:
    ...     for time, mesh in steps:
    ...         w.write(time, mesh)

    The C++ core writes the file; the Python twin, which writes the same bytes,
    when the core declines to open it.
    """

    def __init__(self, filename):
        self._core = self._py = None
        if not is_buffer(filename, "w"):
            try:
                self._core = _core.FemapSeriesWriter(str(filename))
            except Exception as exc:
                if not core_declined(exc, "femap", "write", filename):
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


register_format("femap", [".neu"], read, {"femap": write})

__all__ = ["read", "write", "time_values", "SeriesWriter"]
