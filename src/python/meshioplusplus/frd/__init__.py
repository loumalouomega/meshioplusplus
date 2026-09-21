from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._frd import read as _py_read


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


register_format("frd", [".frd"], read, {})

__all__ = ["read"]
