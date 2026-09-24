from .. import _core
from .._exceptions import ReadError
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._anim import read as _py_read


def read(filename, time_step=0):
    """Read one OpenRadioss / Radioss animation file (``<run>A001``...).

    Each file is one state of the run: its time is ``field_data["meshio:time"]``,
    so a glob over the run's files (``runA*``) is a transient sequence. The
    points are the state's coordinates; nodal scalars and vectors (velocity,
    displacement, contact forces...) are point data; element scalars, tensors
    (xx yy zz xy yz zx) and 1-D force/moment sets are cell data, NaN where an
    element family has none. Parts are cell regions; ``radioss:part``,
    ``radioss:alive`` (0 once an element is deleted) and, when the file has them, ids, masses, materials and
    properties are cell data. The file holds one step: ``time_step`` may be 0
    or -1.
    """
    if time_step not in (0, -1):
        raise ReadError(
            f"Radioss animation: time_step {time_step} is out of range; "
            f"'{filename}' holds 1 step"
        )
    if not is_buffer(filename, "r"):
        try:
            return _core.radioss_anim_read(str(filename))
        except Exception as exc:
            if not core_declined(exc, "radioss_anim", "read", filename):
                raise
    return _py_read(filename)


# Found by its file name (a stem, then A and three or more digits) or its magic
# number, never an extension.
register_format("radioss_anim", [], read, {})

__all__ = ["read"]
