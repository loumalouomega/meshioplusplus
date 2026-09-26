from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._marc import read as _py_read
from ._marc import read_t19 as _py_read_t19
from ._marc import time_values as _py_time_values
from ._marc import write as _py_write


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


register_format("marc", [".dat"], read, {"marc": write})
register_format("marc_t19", [".t19"], read_t19, {})

__all__ = ["read", "read_t19", "time_values", "write"]
