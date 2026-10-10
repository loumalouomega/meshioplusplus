from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._frd import read as _py_read
from ._frd import read_dat
from ._frd import write as _py_write


def read(filename, points_only=False, arrays=None, time_step=0, derived=False):
    """Read a CalculiX result file (``.frd``, ASCII short or long format).

    The mesh is the one ``ccx`` wrote, which is not the ``.inp`` mesh: shells and beams
    are expanded into solids. Every increment is a step: ``time_step`` selects one (0 =
    first, negative counts from the end) and the step's value is
    ``field_data["meshio:time"]``, with ``frd:step`` and ``frd:analysis`` (0 static, 1
    time, 2 frequency, 3 buckling, 4 user) beside it. Each ``-4`` result becomes a
    ``point_data`` array named as the file names it (``DISP``, ``STRESS``, ``TOSTRAIN``,
    ``NDTEMP``, ``FORC``...), NaN at nodes the block has no value for; six-component
    tensors keep the file's order ``xx yy zz xy yz zx``.

    ``derived=True`` adds ``<NAME>_mises`` and ``<NAME>_principal`` (ascending min, mid,
    max) beside each ``STRESS``/``TOSTRAIN``/``MESTRAIN`` tensor. The C++ core reads the
    whole format; the Python reader is the fallback for buffers and for anything the
    C++ path raises on.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.frd_read(
                str(filename), points_only, arrays, time_step, derived
            )
        except Exception as exc:
            if not core_declined(exc, "frd", "read", filename):
                raise
    return _py_read(
        filename,
        points_only=points_only,
        arrays=arrays,
        time_step=time_step,
        derived=derived,
    )


def write(filename, mesh, long_ids=True):
    """Write a CalculiX result file (``.frd``, ASCII).

    The inverse of :func:`read`: the points and the twelve cell types ``ccx`` writes
    (others are dropped with a warning), ``frd:group`` / ``frd:material`` as each
    element's group and material, and the point data as one ``-4`` result block per
    array under its own name, a step's ``meshio:time``, ``frd:step`` and
    ``frd:analysis`` in the ``100C`` header. ``.frd`` is nodal, holds no sets and
    prints six digits, so cell data, regions and other field data are dropped with a
    warning; an array whose name is longer than eight characters, or that is not one
    value or vector per point, is dropped too. ``long_ids=True`` (what ``ccx``
    writes) uses ``I10`` ids, ``False`` the short ``I5`` form (99999 nodes or
    elements at most).
    """
    if not is_buffer(filename, "w"):
        try:
            return _core.frd_write(str(filename), mesh, long_ids)
        except Exception as exc:
            if not core_declined(exc, "frd", "write", filename):
                raise
    return _py_write(filename, mesh, long_ids=long_ids)


register_format("frd", [".frd"], read, {"frd": write})

__all__ = ["read", "read_dat", "write"]
