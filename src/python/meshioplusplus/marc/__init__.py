from .. import _core  # noqa: F401
from .._helpers import register_format
from .._shims import (
    plain_shims,
    series_read_shim,
    series_writer_class,
    time_values_shim,
    write_shim,
)
from ._marc import read as _py_read  # noqa: F401
from ._marc import read_t19 as _py_read_t19  # noqa: F401
from ._marc import time_values as _py_time_values  # noqa: F401
from ._marc import write as _py_write  # noqa: F401
from ._t19 import SeriesWriter as _PySeriesWriter  # noqa: F401
from ._t19 import write as _py_write_t19  # noqa: F401

read, write = plain_shims(
    __name__,
    "marc",
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
    """,
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
    """,
)

read_t19 = series_read_shim(
    __name__,
    "marc_t19",
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
    """,
    name="read_t19",
    py_attr="_py_read_t19",
)

time_values = time_values_shim(
    __name__, "marc_t19", """The time of every increment of a ``.t19`` post file."""
)

write_t19 = write_shim(
    __name__,
    "marc_t19",
    """Write an MSC Marc formatted post file (``.t19``) with one increment: the
    model (nodes, elements with a Marc type, node and element sets), the mesh's
    point arrays of 1, 2, 3, 6 or 9 components as nodal results and its cell
    arrays as element post codes (per integration point where
    ``field_data["marc:layout:<name>"]`` says so, a tensor's six components as
    six consecutive codes), at ``field_data["meshio:time"]``. Both engines write
    the same bytes; the file is read back by meshio++'s own reader and has not
    been opened in Mentat.
    """,
    name="write_t19",
    py_attr="_py_write_t19",
)

SeriesWriter = series_writer_class(
    __name__,
    "marc_t19",
    "MarcT19SeriesWriter",
    """Write a time series into one post file, an increment per mesh;
    ``write_sequence(path, steps, file_format="marc_t19")`` uses it. Every mesh
    must have the first one's points, cells and arrays.

    >>> with meshioplusplus.marc.SeriesWriter("job.t19") as w:
    ...     for time, mesh in steps:
    ...         w.write(time, mesh)
    """,
)

register_format("marc", [".dat"], read, {"marc": write})
register_format("marc_t19", [".t19"], read_t19, {"marc_t19": write_t19})

__all__ = ["read", "read_t19", "time_values", "write", "write_t19", "SeriesWriter"]
