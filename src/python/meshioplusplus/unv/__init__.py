from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._unv import read as _py_read
from ._unv import time_values
from ._unv import write as _py_write


def read(filename, points_only=False, arrays=None, time_step=0):
    """Read an I-DEAS Universal File (``.unv`` / ``.uff``).

    Nodes (2411/781/15), elements (2412/780), permanent groups (as ``mesh.regions``,
    so also ``point_sets``/``cell_sets``), units (164) and results: 2414, 55, 56 and
    the functions of 58/58b. Every result block and every 58 abscissa sample belongs
    to a step (a time, mode or frequency); ``time_step`` selects one (0 = first,
    negative counts from the end) and ``read_metadata(...)["time_values"]`` lists
    them. The C++ core reads the file; the Python reader is the fallback for
    buffers and for anything the C++ path declines.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.unv_read(str(filename), points_only, arrays, time_step)
        except Exception as exc:
            if not core_declined(exc, "unv", "read", filename):
                raise
    return _py_read(
        filename, points_only=points_only, arrays=arrays, time_step=time_step
    )


def write(filename, mesh, code_aster=False, node_dataset=2411):
    """Write an I-DEAS Universal File.

    Nodes, elements, the mesh's point and cell regions as 2467 groups, and
    ``point_data`` / ``cell_data`` as results (dataset 2414, or the legacy 55/56
    with ``code_aster=True``). ``node_dataset`` selects ``2411`` (default) or
    ``781``. The C++ core writes the file; the Python writer is the fallback for
    buffer targets.
    """
    if not is_buffer(filename, "w"):
        try:
            _core.unv_write(
                str(filename), mesh, code_aster=code_aster, node_dataset=node_dataset
            )
            return
        except Exception as exc:
            if not core_declined(exc, "unv", "write", filename):
                raise
    return _py_write(filename, mesh, code_aster=code_aster, node_dataset=node_dataset)


register_format("unv", [".unv", ".uff"], read, {"unv": write})

__all__ = ["read", "write", "time_values"]
