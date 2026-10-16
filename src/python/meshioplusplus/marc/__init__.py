from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._marc import read as _py_read
from ._marc import read_t19 as _py_read_t19
from ._marc import time_values as _py_time_values
from ._marc import write as _py_write
from ._t19 import SeriesWriter as _PySeriesWriter
from ._t19 import write as _py_write_t19


def read(filename):
    """Read an MSC Marc input deck (``.dat``).

    ``COORDINATES`` and ``CONNECTIVITY`` of the model definition make the mesh,
    in fixed (``EXTENDED`` or not) or free format; ``marc:element`` and
    ``marc:type`` cell data keep each cell's element number and Marc element
    type. ``DEFINE ELEMENT SET`` and ``DEFINE NODE SET`` become cell and point
    regions (ranges, ``AND``/``EXCEPT``/``INTERSECT`` and set names resolved).
    A ``.dat`` file that is not a Marc deck raises ``ReadError``, so Tecplot's
    ``.dat`` reader is tried next. The C++ core reads the whole format; the
    Python reader is the fallback for buffers and for anything the C++ path
    raises on.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.marc_read(str(filename))
        except Exception as exc:
            if not core_declined(exc, "marc", "read", filename):
                raise
    return _py_read(filename)


def read_t19(filename, points_only=False, arrays=None, time_step=0):
    """Read one increment of an MSC Marc formatted post file (``.t19``).

    The mesh and sets are as in the input deck. ``time_step`` picks the
    increment (negative counts from the end); its time (a frequency or buckling
    factor in a modal, harmonic or buckling analysis) is
    ``field_data["meshio:time"]``, with ``marc:increment`` and
    ``marc:subincrement``. Nodal vectors are point data named as the file names
    them (``Displacement``, ``Reaction Force`` ...); element post codes are cell
    data per integration point, ``(cells, points * components)`` flattened
    point-major with ``field_data["marc:layout:<name>"]``, a stress or strain
    tensor's six codes combined into one ``xx yy zz xy yz zx`` array.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.marc_t19_read(
                str(filename),
                points_only=points_only,
                arrays=arrays,
                time_step=time_step,
            )
        except Exception as exc:
            if not core_declined(exc, "marc_t19", "read", filename):
                raise
    return _py_read_t19(
        filename, points_only=points_only, arrays=arrays, time_step=time_step
    )


def time_values(filename):
    """The time of every increment of a ``.t19`` post file."""
    try:
        return list(_core.marc_t19_time_values(str(filename)))
    except Exception as exc:
        if not core_declined(exc, "marc_t19", "read", filename):
            raise
    return _py_time_values(filename)


def write(filename, mesh):
    """Write an MSC Marc input deck (``.dat``, ``EXTENDED`` format).

    The parameter section (``TITLE``, ``SIZING``, one ``ELEMENTS`` line per
    type, ``END``), ``CONNECTIVITY``, ``COORDINATES``, a ``DEFINE NODE SET`` per
    point region and a ``DEFINE ELEMENT SET`` per cell region, and ``END
    OPTION``. A cell's Marc type is its ``marc:type`` when that fits it, else a
    default (solids and shells in 3-D, plane-strain elements in a planar mesh; a
    pyramid as a degenerate brick); element numbers are ``marc:element`` when
    positive and unique. ``.dat`` is Tecplot's for a write by extension, so name
    the format: ``meshioplusplus.write("job.dat", mesh, file_format="marc")``.
    Both engines write the same bytes.
    """
    if not is_buffer(filename, "w"):
        try:
            _core.marc_write(str(filename), mesh)
            return
        except Exception as exc:
            if not core_declined(exc, "marc", "write", filename):
                raise
    return _py_write(filename, mesh)


def write_t19(filename, mesh):
    """Write an MSC Marc formatted post file (``.t19``) with one increment: the
    model (nodes, elements with a Marc type, node and element sets), the mesh's
    point arrays of 1, 2, 3, 6 or 9 components as nodal results and its cell
    arrays as element post codes (per integration point where
    ``field_data["marc:layout:<name>"]`` says so, a tensor's six components as
    six consecutive codes), at ``field_data["meshio:time"]``. Both engines write
    the same bytes; the file is read back by meshio++'s own reader and has not
    been opened in Mentat.
    """
    if not is_buffer(filename, "w"):
        try:
            _core.marc_t19_write(str(filename), mesh)
            return
        except Exception as exc:
            if not core_declined(exc, "marc_t19", "write", filename):
                raise
    return _py_write_t19(filename, mesh)


class SeriesWriter:
    """Write a time series into one post file, an increment per mesh;
    ``write_sequence(path, steps, file_format="marc_t19")`` uses it. Every mesh
    must have the first one's points, cells and arrays.

    >>> with meshioplusplus.marc.SeriesWriter("job.t19") as w:
    ...     for time, mesh in steps:
    ...         w.write(time, mesh)
    """

    def __init__(self, filename):
        self._core = self._py = None
        if not is_buffer(filename, "w"):
            try:
                self._core = _core.MarcT19SeriesWriter(str(filename))
            except Exception as exc:
                if not core_declined(exc, "marc_t19", "write", filename):
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


register_format("marc", [".dat"], read, {"marc": write})
register_format("marc_t19", [".t19"], read_t19, {"marc_t19": write_t19})

__all__ = ["read", "read_t19", "time_values", "write", "write_t19", "SeriesWriter"]
