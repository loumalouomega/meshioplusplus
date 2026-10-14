from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._abaqus_fil import read as _py_read
from ._abaqus_fil import time_values as _py_time_values
from ._write import SeriesWriter as _PySeriesWriter
from ._write import write as _py_write


def read(filename, points_only=False, arrays=None, time_step=0):
    """Read an Abaqus results file (``.fil``), binary or ASCII.

    Nodes (1901) and elements (1900/1990) become the mesh, with their labels as
    ``abaqus:id`` point and cell data; node and element sets (1931-1934) become
    point and cell regions, long names resolved through the 1940 labels. Every
    increment is a step: ``time_step`` selects one (0 = first, negative counts
    from the end), its total time is ``field_data["meshio:time"]``, and its
    nodal and element records become point and cell data named by their Abaqus
    identifier (``U``, ``RF``, ``S``, ``E``, ...).
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.abaqus_fil_read(str(filename), points_only, arrays, time_step)
        except Exception as exc:
            if not core_declined(exc, "abaqus_fil", "read", filename):
                raise
    return _py_read(
        filename, points_only=points_only, arrays=arrays, time_step=time_step
    )


def time_values(filename):
    """The total time of every increment, in file order."""
    if not is_buffer(filename, "r"):
        try:
            return list(_core.abaqus_fil_time_values(str(filename)))
        except Exception as exc:
            if not core_declined(exc, "abaqus_fil", "read", filename):
                raise
    return _py_time_values(filename)


def write(filename, mesh, binary=False):
    """Write a mesh as a one-increment Abaqus results file (``.fil``).

    The inverse of :func:`read`. The model is the points (``abaqus:id`` labels,
    else 1 to n), the cells (``abaqus:id`` labels, else 1 to n; each cell type
    is the Abaqus element named in the docs, ``C3D8R`` for a hexahedron...) and
    the point and cell regions as node and element sets. The increment holds the
    point data named like an Abaqus nodal result (``U``, ``V``, ``A``, ``RF``,
    ``key_<n>``...) and the cell data named like an element result (``S``, ``E``,
    ``MISES``, ``NFORC``, ``key_<n>``...): an array with an
    ``abaqus:layout:<name>`` of ``(points, components)`` is per integration
    point (per element node for ``NFORC``), any other one value set per element.
    Its time is ``meshio:time`` (0 when absent). ``binary=True`` writes the
    little-endian binary layout instead of ASCII. Cells with no Abaqus element,
    side regions and other arrays are dropped with a warning.
    """
    if not is_buffer(filename, "w"):
        try:
            _core.abaqus_fil_write(str(filename), mesh, bool(binary))
            return
        except Exception as exc:
            if not core_declined(exc, "abaqus_fil", "write", filename):
                raise
    return _py_write(filename, mesh, binary=binary)


class SeriesWriter:
    """Write a time series into one results file: the model with the first step
    and an increment per step, which ``read`` gives back by ``time_step``. A
    step's time is the argument, not its ``meshio:time``; its step, increment,
    step time and procedure are its ``abaqus:*`` field data. Every step must
    have the first step's points and cell blocks.
    ``write_sequence(path, steps, file_format="abaqus_fil")`` uses it.

    >>> with meshioplusplus.abaqus_fil.SeriesWriter("job.fil") as w:
    ...     for time, mesh in steps:
    ...         w.write(time, mesh)

    The C++ core writes the file; the Python twin, which writes the same bytes,
    when the core declines to open it.
    """

    def __init__(self, filename, binary=False):
        self._core = self._py = None
        if not is_buffer(filename, "w"):
            try:
                self._core = _core.AbaqusFilSeriesWriter(str(filename), bool(binary))
            except Exception as exc:
                if not core_declined(exc, "abaqus_fil", "write", filename):
                    raise
        if self._core is None:
            self._py = _PySeriesWriter(filename, binary=binary)

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


register_format("abaqus_fil", [".fil"], read, {"abaqus_fil": write})

__all__ = ["read", "write", "time_values", "SeriesWriter"]
